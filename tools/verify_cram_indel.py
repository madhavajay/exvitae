#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


@dataclass
class VariantSite:
    yaml_path: Path
    rsid: str
    kind: str
    grch37_chrom: str | None
    grch37_start: int | None
    grch37_end: int | None
    grch38_chrom: str | None
    grch38_start: int | None
    grch38_end: int | None
    ref: str
    alts: list[str]


@dataclass
class NormalizedCall:
    chrom: str
    pos: int
    ref: str
    alt: str
    qual: str
    filt: str
    info: str
    fmt: str
    sample: str


@dataclass
class PileupSummary:
    total_reads: int
    spanning_reads: int
    reference_like_reads: int
    deletion_1bp_reads: int
    insertion_1bp_reads: int
    insertion_2bp_reads: int
    other_indel_reads: int


@dataclass
class ReadDebug:
    read_name: str
    start: int
    cigar: str
    spans_locus: bool
    classification: str
    relevant_indels: list[str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Verify a YAML-defined indel in a CRAM using samtools/bcftools."
    )
    parser.add_argument("--variant", required=True, type=Path, help="Variant YAML path")
    parser.add_argument("--cram", required=True, type=Path, help="Indexed CRAM path")
    parser.add_argument("--reference", required=True, type=Path, help="Reference FASTA path")
    parser.add_argument(
        "--region",
        help="Explicit region override such as 19:45679760-45679800",
    )
    parser.add_argument(
        "--window",
        type=int,
        default=14,
        help="Half-window padding around the variant locus when deriving region",
    )
    parser.add_argument(
        "--dump-reads",
        action="store_true",
        help="Include a samtools view dump for the chosen region",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Emit machine-readable JSON summary",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Show detailed per-read and mpileup debug information",
    )
    return parser.parse_args()


def require_binary(name: str) -> str:
    path = shutil.which(name)
    if not path:
        raise SystemExit(f"missing required binary: {name}")
    return path


def load_variant(path: Path) -> VariantSite:
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"variant file is not a mapping: {path}")

    identifiers = data.get("identifiers") or {}
    rsids = identifiers.get("rsids") or []
    alleles = data.get("alleles") or {}
    coords = data.get("coordinates") or {}
    grch37 = coords.get("grch37") or {}
    grch38 = coords.get("grch38") or {}

    rsid = rsids[0] if rsids else path.stem
    kind = str(alleles.get("kind") or "")
    ref = str(alleles.get("ref") or "")
    alts = [str(value) for value in (alleles.get("alts") or [])]
    if not ref or not alts:
        raise SystemExit(f"variant file is missing ref/alts: {path}")

    return VariantSite(
        yaml_path=path,
        rsid=rsid,
        kind=kind,
        grch37_chrom=str(grch37.get("chrom")) if grch37.get("chrom") is not None else None,
        grch37_start=as_int(grch37.get("start") or grch37.get("pos")),
        grch37_end=as_int(grch37.get("end") or grch37.get("pos")),
        grch38_chrom=str(grch38.get("chrom")) if grch38.get("chrom") is not None else None,
        grch38_start=as_int(grch38.get("start") or grch38.get("pos")),
        grch38_end=as_int(grch38.get("end") or grch38.get("pos")),
        ref=ref,
        alts=alts,
    )


def as_int(value: Any) -> int | None:
    if value is None:
        return None
    return int(value)


def detect_assembly(reference: Path) -> str:
    name = reference.name.lower()
    full = str(reference).lower()
    if any(token in name or token in full for token in ("grch38", "hg38", "assembly38")):
        return "grch38"
    if any(token in name or token in full for token in ("grch37", "hg19", "assembly37")):
        return "grch37"
    raise SystemExit(
        f"could not infer assembly from reference path {reference}; use a GRCh37/GRCh38-named file"
    )


def choose_region(variant: VariantSite, assembly: str, window: int, override: str | None) -> str:
    if override:
        return override
    chrom, start, end = variant_locus(variant, assembly)
    return f"{chrom}:{max(1, start - window)}-{end + window}"


def variant_locus(variant: VariantSite, assembly: str) -> tuple[str, int, int]:
    if assembly == "grch38":
        if variant.grch38_chrom and variant.grch38_start and variant.grch38_end:
            return variant.grch38_chrom, variant.grch38_start, variant.grch38_end
    if assembly == "grch37":
        if variant.grch37_chrom and variant.grch37_start and variant.grch37_end:
            return variant.grch37_chrom, variant.grch37_start, variant.grch37_end
    raise SystemExit(f"variant {variant.rsid} does not have {assembly} coordinates")


def run_command(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, capture_output=True, text=True, check=False)


def run_pipeline(cram: Path, reference: Path, region: str) -> tuple[str, str, int]:
    command = (
        f"bcftools mpileup -f {shell_quote(reference)} -r {shell_quote(region)} -Ou {shell_quote(cram)} | "
        f"bcftools call -mv -Ou | "
        f"bcftools norm -f {shell_quote(reference)} -Ov"
    )
    completed = subprocess.run(
        ["bash", "-lc", command],
        capture_output=True,
        text=True,
        check=False,
    )
    return completed.stdout, completed.stderr, completed.returncode


def shell_quote(value: str | Path) -> str:
    text = str(value)
    return "'" + text.replace("'", "'\"'\"'") + "'"


def parse_vcf_calls(vcf_text: str) -> list[NormalizedCall]:
    calls: list[NormalizedCall] = []
    for line in vcf_text.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 10:
            continue
        chrom, pos, _id, ref, alts, qual, filt, info, fmt, sample = fields[:10]
        for alt in alts.split(","):
            calls.append(
                NormalizedCall(
                    chrom=chrom,
                    pos=int(pos),
                    ref=ref,
                    alt=alt,
                    qual=qual,
                    filt=filt,
                    info=info,
                    fmt=fmt,
                    sample=sample,
                )
            )
    return calls


def summarize_pileup_from_sam(
    sam_text: str, locus_start: int, locus_end: int
) -> PileupSummary:
    total_reads = 0
    spanning_reads = 0
    reference_like_reads = 0
    deletion_1bp_reads = 0
    insertion_1bp_reads = 0
    insertion_2bp_reads = 0
    other_indel_reads = 0

    for line in sam_text.splitlines():
        if not line or line.startswith("@"):
            continue
        fields = line.split("\t")
        if len(fields) < 6:
            continue
        total_reads += 1
        start_pos = int(fields[3])
        cigar = fields[5]
        spans, indels = analyze_cigar_over_locus(start_pos, cigar, locus_start, locus_end)
        if spans:
            spanning_reads += 1

        relevant_indel = False
        for op, length, anchor in indels:
            if anchor < (locus_start - 1) or anchor > locus_end:
                continue
            relevant_indel = True
            if op == "D" and length == 1:
                deletion_1bp_reads += 1
            elif op == "I" and length == 1:
                insertion_1bp_reads += 1
            elif op == "I" and length == 2:
                insertion_2bp_reads += 1
            else:
                other_indel_reads += 1

        if spans and not relevant_indel:
            reference_like_reads += 1

    return PileupSummary(
        total_reads=total_reads,
        spanning_reads=spanning_reads,
        reference_like_reads=reference_like_reads,
        deletion_1bp_reads=deletion_1bp_reads,
        insertion_1bp_reads=insertion_1bp_reads,
        insertion_2bp_reads=insertion_2bp_reads,
        other_indel_reads=other_indel_reads,
    )


def collect_read_debug(
    sam_text: str, locus_start: int, locus_end: int
) -> list[ReadDebug]:
    rows: list[ReadDebug] = []
    for line in sam_text.splitlines():
        if not line or line.startswith("@"):
            continue
        fields = line.split("\t")
        if len(fields) < 6:
            continue
        read_name = fields[0]
        start_pos = int(fields[3])
        cigar = fields[5]
        spans, indels = analyze_cigar_over_locus(start_pos, cigar, locus_start, locus_end)
        relevant = []
        classification = "off_locus"
        for op, length, anchor in indels:
            if anchor < (locus_start - 1) or anchor > locus_end:
                continue
            relevant.append(f"{op}{length}@{anchor}")
        if spans and not relevant:
            classification = "reference_like"
        elif relevant:
            if any(token.startswith("D1@") for token in relevant):
                classification = "deletion_1bp"
            elif any(token.startswith("I1@") for token in relevant):
                classification = "insertion_1bp"
            elif any(token.startswith("I2@") for token in relevant):
                classification = "insertion_2bp"
            else:
                classification = "other_indel"
        elif spans:
            classification = "spans_locus"
        rows.append(
            ReadDebug(
                read_name=read_name,
                start=start_pos,
                cigar=cigar,
                spans_locus=spans,
                classification=classification,
                relevant_indels=relevant,
            )
        )
    return rows


def analyze_cigar_over_locus(
    start_pos: int, cigar: str, locus_start: int, locus_end: int
) -> tuple[bool, list[tuple[str, int, int]]]:
    ref_pos = start_pos
    covers_start = False
    covers_end = False
    indels: list[tuple[str, int, int]] = []

    for length_text, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        length = int(length_text)
        if op in {"M", "=", "X"}:
            block_start = ref_pos
            block_end = ref_pos + length - 1
            if block_start <= locus_start <= block_end:
                covers_start = True
            if block_start <= locus_end <= block_end:
                covers_end = True
            ref_pos += length
        elif op == "D":
            anchor = ref_pos - 1
            indels.append((op, length, anchor))
            deleted_start = ref_pos
            deleted_end = ref_pos + length - 1
            if deleted_start <= locus_start <= deleted_end:
                covers_start = True
            if deleted_start <= locus_end <= deleted_end:
                covers_end = True
            ref_pos += length
        elif op == "I":
            anchor = ref_pos - 1
            indels.append((op, length, anchor))
        elif op == "N":
            ref_pos += length

    return covers_start and covers_end, indels


def classify_call(variant: VariantSite, calls: list[NormalizedCall]) -> tuple[str, NormalizedCall | None]:
    expected = {alt: f"{len(alt)}T" for alt in variant.alts}
    for call in calls:
        if call.ref == variant.ref and call.alt in expected:
            return expected[call.alt], call
    for call in calls:
        if same_length_delta(variant.ref, call.ref, call.alt):
            label = label_from_delta(variant.ref, call.alt)
            if label:
                return label, call
    return "none", calls[0] if calls else None


def same_length_delta(variant_ref: str, call_ref: str, call_alt: str) -> bool:
    return len(call_ref) - len(call_alt) == len(variant_ref) - len(call_alt)


def label_from_delta(variant_ref: str, alt: str) -> str | None:
    if alt and set(alt) == {"T"}:
        return f"{len(alt)}T"
    delta = len(alt) - len(variant_ref)
    if delta == -1:
        return f"{len(alt)}T"
    if delta in {1, 2}:
        return f"{len(alt)}T"
    return None


def run_faidx(reference: Path, region: str) -> subprocess.CompletedProcess[str]:
    return run_command(["samtools", "faidx", str(reference), region])


def run_view(cram: Path, region: str) -> subprocess.CompletedProcess[str]:
    raise AssertionError("use run_view_with_reference")


def run_view_with_reference(
    cram: Path, reference: Path, region: str
) -> subprocess.CompletedProcess[str]:
    return run_command(["samtools", "view", "-T", str(reference), str(cram), region])


def run_mpileup_with_reference(
    cram: Path, reference: Path, region: str
) -> subprocess.CompletedProcess[str]:
    return run_command(
        ["samtools", "mpileup", "-f", str(reference), "-r", region, str(cram)]
    )


def normalize_region_to_reference(reference: Path, region: str) -> str:
    raw_chrom, span = region.split(":", 1)
    fai = Path(str(reference) + ".fai")
    contigs = []
    if fai.exists():
        for line in fai.read_text(encoding="utf-8").splitlines():
            fields = line.split("\t")
            if fields:
                contigs.append(fields[0])
    candidates = [raw_chrom]
    if raw_chrom.startswith("chr"):
        candidates.append(raw_chrom.removeprefix("chr"))
    else:
        candidates.append(f"chr{raw_chrom}")
    for candidate in candidates:
        if candidate in contigs:
            return f"{candidate}:{span}"
    return region


def emit_text(
    *,
    variant: VariantSite,
    assembly: str,
    locus: tuple[str, int, int],
    region: str,
    faidx: subprocess.CompletedProcess[str],
    reads: subprocess.CompletedProcess[str] | None,
    pileup_summary: PileupSummary | None,
    read_debug: list[ReadDebug] | None,
    mpileup_debug: subprocess.CompletedProcess[str] | None,
    pipeline_stdout: str,
    pipeline_stderr: str,
    pipeline_code: int,
    calls: list[NormalizedCall],
    match_label: str,
    matched_call: NormalizedCall | None,
) -> None:
    print(f"variant={variant.rsid}")
    print(f"kind={variant.kind}")
    print(f"assembly={assembly}")
    print(f"locus={locus[0]}:{locus[1]}-{locus[2]}")
    print(f"region={region}")
    print(f"ref={variant.ref}")
    print(f"alts={','.join(variant.alts)}")
    print()
    print("== reference slice ==")
    print((faidx.stdout or "").rstrip() or "<empty>")
    if faidx.stderr.strip():
        print(faidx.stderr.rstrip())
    print()
    if reads is not None:
        print("== reads ==")
        body = (reads.stdout or "").splitlines()
        if body:
            print("\n".join(body[:40]))
            if len(body) > 40:
                print(f"... ({len(body) - 40} more reads)")
        else:
            print("<empty>")
        if reads.stderr.strip():
            print(reads.stderr.rstrip())
        print()
    if pileup_summary is not None:
        print("== pileup summary ==")
        print(f"total_reads={pileup_summary.total_reads}")
        print(f"spanning_reads={pileup_summary.spanning_reads}")
        print(f"reference_like_reads={pileup_summary.reference_like_reads}")
        print(f"deletion_1bp_reads={pileup_summary.deletion_1bp_reads}")
        print(f"insertion_1bp_reads={pileup_summary.insertion_1bp_reads}")
        print(f"insertion_2bp_reads={pileup_summary.insertion_2bp_reads}")
        print(f"other_indel_reads={pileup_summary.other_indel_reads}")
        print()
    if read_debug is not None:
        print("== read classifications ==")
        interesting = [
            row for row in read_debug if row.spans_locus or row.relevant_indels
        ]
        for row in interesting[:80]:
            suffix = f" indels={','.join(row.relevant_indels)}" if row.relevant_indels else ""
            print(
                f"{row.read_name} start={row.start} cigar={row.cigar} "
                f"spans_locus={str(row.spans_locus).lower()} class={row.classification}{suffix}"
            )
        if len(interesting) > 80:
            print(f"... ({len(interesting) - 80} more classified reads)")
        print()
    if mpileup_debug is not None:
        print("== samtools mpileup ==")
        print((mpileup_debug.stdout or "").rstrip() or "<empty>")
        if mpileup_debug.stderr.strip():
            print(mpileup_debug.stderr.rstrip())
        print()
    print("== normalized calls ==")
    if pipeline_stdout.strip():
        print(pipeline_stdout.rstrip())
    else:
        print("<empty>")
    if pipeline_stderr.strip():
        print()
        print("== caller stderr ==")
        print(pipeline_stderr.rstrip())
    print()
    print(f"pipeline_exit_code={pipeline_code}")
    print(f"present_alt={match_label}")
    if matched_call is not None:
        print(
            "matched_call="
            f"{matched_call.chrom}:{matched_call.pos} {matched_call.ref}>{matched_call.alt} "
            f"sample={matched_call.sample}"
        )
    else:
        print("matched_call=<none>")
    print(f"call_count={len(calls)}")


def emit_json(
    *,
    variant: VariantSite,
    assembly: str,
    locus: tuple[str, int, int],
    region: str,
    faidx: subprocess.CompletedProcess[str],
    reads: subprocess.CompletedProcess[str] | None,
    pileup_summary: PileupSummary | None,
    read_debug: list[ReadDebug] | None,
    mpileup_debug: subprocess.CompletedProcess[str] | None,
    pipeline_stdout: str,
    pipeline_stderr: str,
    pipeline_code: int,
    calls: list[NormalizedCall],
    match_label: str,
    matched_call: NormalizedCall | None,
) -> None:
    payload = {
        "variant": variant.rsid,
        "kind": variant.kind,
        "assembly": assembly,
        "locus": f"{locus[0]}:{locus[1]}-{locus[2]}",
        "region": region,
        "ref": variant.ref,
        "alts": variant.alts,
        "reference_slice": faidx.stdout,
        "reads": reads.stdout if reads is not None else None,
        "pileup_summary": pileup_summary.__dict__ if pileup_summary else None,
        "read_debug": [row.__dict__ for row in read_debug] if read_debug is not None else None,
        "samtools_mpileup": mpileup_debug.stdout if mpileup_debug is not None else None,
        "pipeline_exit_code": pipeline_code,
        "pipeline_stdout": pipeline_stdout,
        "pipeline_stderr": pipeline_stderr,
        "calls": [call.__dict__ for call in calls],
        "present_alt": match_label,
        "matched_call": matched_call.__dict__ if matched_call else None,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


def main() -> int:
    args = parse_args()
    require_binary("samtools")
    require_binary("bcftools")

    variant = load_variant(args.variant)
    assembly = detect_assembly(args.reference)
    locus = variant_locus(variant, assembly)
    region = normalize_region_to_reference(
        args.reference, choose_region(variant, assembly, args.window, args.region)
    )

    faidx = run_faidx(args.reference, region)
    reads = (
        run_view_with_reference(args.cram, args.reference, region)
        if args.dump_reads
        else None
    )
    pileup_summary = (
        summarize_pileup_from_sam(reads.stdout, locus[1], locus[2]) if reads is not None else None
    )
    read_debug = (
        collect_read_debug(reads.stdout, locus[1], locus[2]) if args.debug and reads is not None else None
    )
    mpileup_debug = (
        run_mpileup_with_reference(args.cram, args.reference, region) if args.debug else None
    )
    pipeline_stdout, pipeline_stderr, pipeline_code = run_pipeline(
        args.cram, args.reference, region
    )
    calls = parse_vcf_calls(pipeline_stdout)
    match_label, matched_call = classify_call(variant, calls)

    if args.json:
        emit_json(
            variant=variant,
            assembly=assembly,
            locus=locus,
            region=region,
            faidx=faidx,
            reads=reads,
            pileup_summary=pileup_summary,
            read_debug=read_debug,
            mpileup_debug=mpileup_debug,
            pipeline_stdout=pipeline_stdout,
            pipeline_stderr=pipeline_stderr,
            pipeline_code=pipeline_code,
            calls=calls,
            match_label=match_label,
            matched_call=matched_call,
        )
    else:
        emit_text(
            variant=variant,
            assembly=assembly,
            locus=locus,
            region=region,
            faidx=faidx,
            reads=reads,
            pileup_summary=pileup_summary,
            read_debug=read_debug,
            mpileup_debug=mpileup_debug,
            pipeline_stdout=pipeline_stdout,
            pipeline_stderr=pipeline_stderr,
            pipeline_code=pipeline_code,
            calls=calls,
            match_label=match_label,
            matched_call=matched_call,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
