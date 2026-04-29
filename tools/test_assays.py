#!/usr/bin/env python3
from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import html
import io
import json
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from fnmatch import fnmatch
import shlex
import shutil
import subprocess
import sys
import textwrap
import time
import yaml
import zipfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parent.parent
WORKSPACE_ROOT = REPO_ROOT.parent
BS = REPO_ROOT / "bioscript" / "bs"
ASSAY_DIR = REPO_ROOT / "assays"
DATA_DIR = REPO_ROOT / "test-data"
OUTPUT_DIR = REPO_ROOT / "test-output"
RESULTS_CSV = OUTPUT_DIR / "results.csv"
REPORT_DIR = OUTPUT_DIR / "report"
REPORT_RUNS_DIR = REPORT_DIR / "runs"
REPORT_DATA_DIR = REPORT_DIR / "data"
AGGREGATE_OBSERVATIONS_TSV = REPORT_DATA_DIR / "observations.tsv"
AGGREGATE_REPORTS_JSONL = REPORT_DATA_DIR / "reports.jsonl"
AGGREGATE_REPORTS_JSON = REPORT_DATA_DIR / "reports.json"
GENERATED_SCRIPT_DIR = REPORT_DIR / ".generated"
DEFAULT_REF_FA = DATA_DIR / "1k-genomes" / "ref" / "GRCh38_full_analysis_set_plus_decoy_hla.fa"
DEFAULT_REF_FAI = DATA_DIR / "1k-genomes" / "ref" / "GRCh38_full_analysis_set_plus_decoy_hla.fa.fai"
HARD_TIMEOUT_SECONDS = 60
BIOSCRIPT_BINARY = REPO_ROOT / "bioscript" / "rust" / "target" / "debug" / "bioscript"


@dataclass
class VariantDefinition:
    name: str
    fields: dict[str, Any]
    source: str
    path: Path | None = None
    manifest_name: str | None = None


OBSERVATION_HEADERS = [
    "participant_id",
    "assay_id",
    "assay_version",
    "variant_key",
    "rsid",
    "assembly",
    "chrom",
    "pos_start",
    "pos_end",
    "ref",
    "alt",
    "kind",
    "match_status",
    "coverage_status",
    "call_status",
    "genotype",
    "genotype_display",
    "zygosity",
    "ref_count",
    "alt_count",
    "depth",
    "genotype_quality",
    "allele_balance",
    "outcome",
    "evidence_type",
    "evidence_raw",
    "facets",
]


def ensure_bioscript_binary() -> Path:
    env = os.environ.copy()
    xcrun = shutil.which("xcrun")
    if xcrun:
        sdkroot = env.get("SDKROOT")
        if not sdkroot:
            probe = subprocess.run([xcrun, "--show-sdk-path"], cwd=REPO_ROOT, capture_output=True, text=True, check=False)
            if probe.returncode == 0:
                sdkroot = probe.stdout.strip()
                env["SDKROOT"] = sdkroot
        if sdkroot and "BINDGEN_EXTRA_CLANG_ARGS" not in env:
            env["BINDGEN_EXTRA_CLANG_ARGS"] = f"-isysroot {sdkroot}"
    result = subprocess.run(
        [
            "cargo",
            "build",
            "--manifest-path",
            str(REPO_ROOT / "bioscript" / "rust" / "Cargo.toml"),
            "-p",
            "bioscript-cli",
        ],
        cwd=REPO_ROOT,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(result.stderr.strip() or result.stdout.strip() or "failed to build bioscript")
    return BIOSCRIPT_BINARY


def find_assays() -> list[Path]:
    assays = sorted(ASSAY_DIR.rglob("*.py"))
    seen_dirs: set[Path] = {p.parent for p in assays}
    for yaml_dir in sorted(ASSAY_DIR.rglob("*.yaml")):
        parent = yaml_dir.parent
        if parent not in seen_dirs and is_yaml_assay_dir(parent):
            seen_dirs.add(parent)
            assays.append(parent)
    return sorted(assays)


def find_inputs(base_dir: Path) -> list[Path]:
    results: list[Path] = []
    for path in sorted(base_dir.rglob("*")):
        if not path.is_file():
            continue
        if path.name == ".DS_Store":
            continue
        suffixes = "".join(path.suffixes[-2:])
        allowed = (
            path.suffix == ".zip"
            or path.suffix == ".txt"
            or path.suffix == ".vcf"
            or path.suffix == ".cram"
            or suffixes == ".txt.zip"
            or suffixes == ".vcf.gz"
        )
        if not allowed:
            continue
        if path.suffix in {".fa", ".fasta", ".fai", ".crai", ".tbi", ".blake3"}:
            continue
        if path.name.startswith(("Homo_sapiens", "GRCh38")):
            continue
        results.append(path)
    return results


def detect_format(input_path: Path) -> str:
    name = input_path.name
    if name.endswith(".cram"):
        return "cram"
    if name.endswith(".vcf.gz") or name.endswith(".vcf"):
        return "vcf"
    if name.endswith(".txt.zip") or name.endswith(".zip"):
        return "zip"
    if name.endswith(".txt"):
        return "text"
    return "unknown"


def detect_source(input_path: Path) -> str:
    try:
        rel = input_path.relative_to(DATA_DIR)
        parts = rel.parts[:2]
        return "/".join(parts) if parts else rel.parent.as_posix()
    except ValueError:
        return input_path.parent.name


def participant_from_path(input_path: Path) -> str:
    name = input_path.name
    for suffix in (".txt.zip", ".vcf.gz", ".final.cram", ".zip", ".txt", ".vcf", ".cram"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return input_path.stem


def assay_name_from_path(assay_path: Path) -> str:
    return assay_path.stem


def output_path_for(assay_path: Path, input_path: Path) -> Path:
    return OUTPUT_DIR / assay_name_from_path(assay_path) / detect_source(input_path) / f"{participant_from_path(input_path)}.tsv"


def trace_path_for(assay_path: Path, input_path: Path) -> Path:
    return REPORT_DATA_DIR / assay_name_from_path(assay_path) / detect_source(input_path) / f"{participant_from_path(input_path)}.trace.tsv"


def genome_reads_path_for(assay_path: Path, input_path: Path) -> Path:
    return REPORT_DATA_DIR / assay_name_from_path(assay_path) / detect_source(input_path) / f"{participant_from_path(input_path)}.genome_reads.tsv"


def timing_path_for(assay_path: Path, input_path: Path) -> Path:
    return REPORT_DATA_DIR / assay_name_from_path(assay_path) / detect_source(input_path) / f"{participant_from_path(input_path)}.timings.tsv"


def observations_path_for(assay_path: Path, input_path: Path) -> Path:
    return REPORT_DATA_DIR / assay_name_from_path(assay_path) / detect_source(input_path) / f"{participant_from_path(input_path)}.observations.tsv"


def report_json_path_for(assay_path: Path, input_path: Path) -> Path:
    return REPORT_DATA_DIR / assay_name_from_path(assay_path) / detect_source(input_path) / f"{participant_from_path(input_path)}.report.json"


def details_json_path_for(assay_path: Path, input_path: Path) -> Path:
    return REPORT_DATA_DIR / assay_name_from_path(assay_path) / detect_source(input_path) / f"{participant_from_path(input_path)}.json"


def detail_page_path_for(assay_path: Path, input_path: Path) -> Path:
    return REPORT_RUNS_DIR / assay_name_from_path(assay_path) / detect_source(input_path) / f"{participant_from_path(input_path)}.html"


def detail_markdown_path_for(assay_path: Path, input_path: Path) -> Path:
    return REPORT_RUNS_DIR / assay_name_from_path(assay_path) / detect_source(input_path) / f"{participant_from_path(input_path)}.md"


def find_reference_for(input_path: Path) -> Path | None:
    for candidate in sorted(input_path.parent.glob("*.fa")) + sorted(input_path.parent.glob("*.fasta")):
        if candidate.is_file():
            return candidate
    return None


def needs_reference(input_path: Path) -> bool:
    return detect_format(input_path) == "cram"


def read_tsv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def read_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8")


def read_tsv_text_preview(path: Path, max_lines: int = 80) -> str:
    if not path.exists():
        return ""
    lines = path.read_text(encoding="utf-8").splitlines()
    return "\n".join(lines[:max_lines])


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    if path.is_dir():
        for child in sorted(path.rglob("*.yaml")):
            digest.update(child.read_bytes())
        return digest.hexdigest()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_text(*args: str) -> str:
    result = subprocess.run(
        ["git", *args],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        return ""
    return result.stdout.strip()


def assay_version_metadata(assay_path: Path) -> dict[str, str]:
    rel_path = os.path.relpath(assay_path, REPO_ROOT)
    commit = git_text("rev-parse", "HEAD")
    short_commit = git_text("rev-parse", "--short", "HEAD")
    status = git_text("status", "--short", "--", rel_path)
    return {
        "git_commit": commit,
        "git_commit_short": short_commit,
        "git_status": status,
        "file_sha256": sha256_file(assay_path),
    }


def render_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=True, sort_keys=True)


def literal_value(node: ast.AST) -> Any:
    return ast.literal_eval(node)


def extract_variant_definitions(assay_path: Path) -> list[VariantDefinition]:
    source = assay_path.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(assay_path))
    variants: list[VariantDefinition] = []
    for node in tree.body:
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        target = node.targets[0]
        if not isinstance(target, ast.Name):
            continue
        call = node.value
        if not isinstance(call, ast.Call):
            continue
        func = call.func
        if not (
            isinstance(func, ast.Attribute)
            and isinstance(func.value, ast.Name)
            and func.value.id == "bioscript"
            and func.attr == "variant"
        ):
            continue
        fields: dict[str, Any] = {}
        for keyword in call.keywords:
            if keyword.arg is None:
                continue
            fields[keyword.arg] = literal_value(keyword.value)
        snippet = ast.get_source_segment(source, node) or target.id
        variants.append(
            VariantDefinition(
                name=target.id,
                fields=fields,
                source=snippet,
                path=assay_path,
                manifest_name=target.id,
            )
        )
    return variants


def safe_identifier(value: str) -> str:
    chars = []
    for char in value:
        if char.isascii() and (char.isalnum() or char == "_"):
            chars.append(char)
        else:
            chars.append("_")
    identifier = "".join(chars).strip("_")
    while "__" in identifier:
        identifier = identifier.replace("__", "_")
    if not identifier:
        identifier = "variant"
    if identifier[0].isdigit():
        identifier = f"v_{identifier}"
    return identifier


def yaml_to_variant_definition(yaml_path: Path) -> VariantDefinition | None:
    data = yaml.safe_load(yaml_path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        return None
    schema = data.get("schema", "")
    if not isinstance(schema, str) or "variant" not in schema:
        return None

    coords = data.get("coordinates", {}) or {}
    alleles_block = data.get("alleles", {}) or {}
    identifiers = data.get("identifiers", {}) or {}
    rsids = identifiers.get("rsids", []) or []
    rsid = rsids[0] if rsids else None

    def fmt_coord(c: dict) -> str | None:
        if not c:
            return None
        chrom = c.get("chrom")
        if not chrom:
            return None
        pos = c.get("pos")
        start = c.get("start")
        end = c.get("end")
        if pos is not None:
            return f"{chrom}:{pos}-{pos}"
        if start is not None and end is not None:
            return f"{chrom}:{start}-{end}"
        if start is not None:
            return f"{chrom}:{start}-{start}"
        return None

    fields: dict[str, Any] = {}
    if rsid:
        fields["rsid"] = rsid
    grch37 = fmt_coord(coords.get("grch37", {}))
    grch38 = fmt_coord(coords.get("grch38", {}))
    if grch37:
        fields["grch37"] = grch37
    if grch38:
        fields["grch38"] = grch38
    ref = alleles_block.get("ref")
    if ref:
        fields["ref"] = ref
    alts = alleles_block.get("alts", [])
    if alts:
        fields["alt"] = alts[0]
    kind = alleles_block.get("kind")
    if kind:
        kind_map = {"snv": "snp", "deletion": "deletion", "insertion": "insertion", "indel": "indel"}
        fields["kind"] = kind_map.get(kind, kind)

    name_str = data.get("name", "") or yaml_path.stem
    safe_name = safe_identifier(name_str)
    source_text = yaml_path.read_text(encoding="utf-8")
    return VariantDefinition(
        name=safe_name,
        fields=fields,
        source=source_text,
        path=yaml_path,
        manifest_name=name_str,
    )


def yaml_to_panel_variant_definitions(yaml_path: Path) -> list[VariantDefinition]:
    data = yaml.safe_load(yaml_path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        return []
    schema = data.get("schema", "")
    if schema != "bioscript:panel:1.0":
        return []
    members = data.get("members", []) or []
    if not isinstance(members, list):
        return []
    variants: list[VariantDefinition] = []
    for member in members:
        if not isinstance(member, dict):
            continue
        if member.get("kind") != "variant":
            continue
        path_value = member.get("path")
        if not isinstance(path_value, str) or not path_value.strip():
            continue
        target = (yaml_path.parent / path_value).resolve()
        vd = yaml_to_variant_definition(target)
        if vd:
            variants.append(vd)
    return variants


def extract_yaml_variant_definitions(directory: Path) -> list[VariantDefinition]:
    variants = []
    for yaml_file in sorted(directory.glob("*.yaml")):
        vd = yaml_to_variant_definition(yaml_file)
        if vd:
            variants.append(vd)
    return variants


def is_yaml_variant_file(path: Path) -> bool:
    if not path.is_file():
        return False
    if path.suffix not in {".yaml", ".yml"}:
        return False
    return yaml_to_variant_definition(path) is not None


def is_yaml_panel_file(path: Path) -> bool:
    if not path.is_file():
        return False
    if path.suffix not in {".yaml", ".yml"}:
        return False
    return bool(yaml_to_panel_variant_definitions(path))


def is_yaml_assay_dir(path: Path) -> bool:
    if not path.is_dir():
        return False
    return any(path.glob("*.yaml"))


def bioscript_literal(value: Any) -> str:
    if value is None:
        return "None"
    if isinstance(value, bool):
        return "True" if value else "False"
    if isinstance(value, (int, float)):
        return repr(value)
    if isinstance(value, str):
        return repr(value)
    if isinstance(value, list):
        return "[" + ", ".join(bioscript_literal(item) for item in value) + "]"
    raise TypeError(f"unsupported bioscript literal: {value!r}")


def build_probe_script(variants: list[VariantDefinition]) -> str:
    blocks: list[str] = []
    for variant in variants:
        args = ",\n".join(
            f"    {key}={bioscript_literal(value)}" for key, value in variant.fields.items()
        )
        blocks.append(f"{variant.name} = bioscript.variant(\n{args}\n)\n")

    row_lines = []
    for variant in variants:
        row_lines.append(
            textwrap.dedent(
                f"""\
                {{
                    "variant_name": {variant.name!r},
                    "rsid": {bioscript_literal(variant.fields.get("rsid"))},
                    "kind": {bioscript_literal(variant.fields.get("kind"))},
                    "grch37": {bioscript_literal(variant.fields.get("grch37"))},
                    "grch38": {bioscript_literal(variant.fields.get("grch38"))},
                    "ref": {bioscript_literal(variant.fields.get("ref"))},
                    "alt": {bioscript_literal(variant.fields.get("alt"))},
                    "genotype": genotypes.lookup_variant({variant.name}),
                }},"""
            )
        )

    row_block = textwrap.indent("\n".join(row_lines), " " * 8)

    return textwrap.dedent(
        """\
        """
    ) + "\n".join(blocks) + "\n" + textwrap.dedent(
        f"""\
        def main():
            genotypes = bioscript.load_genotypes(input_file)
            rows = [
        {row_block}
            ]
            bioscript.write_tsv(output_file, rows)


        if __name__ == "__main__":
            main()
        """
    )


def ensure_probe_script(assay_path: Path, variants: list[VariantDefinition]) -> Path:
    GENERATED_SCRIPT_DIR.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256((assay_path.as_posix() + "\n" + render_json([v.fields for v in variants])).encode("utf-8")).hexdigest()[:16]
    target = GENERATED_SCRIPT_DIR / f"{assay_path.stem}.{digest}.probe.py"
    target.write_text(build_probe_script(variants), encoding="utf-8")
    return target


def relative_to(path: Path, root: Path) -> str:
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return path.as_posix()


def format_duration(ms: int) -> str:
    if ms >= 1000:
        return f"{ms / 1000:.2f}s"
    return f"{ms}ms"


def summarize_output(rows: list[dict[str, str]]) -> str:
    if not rows:
        return ""
    if len(rows) > 1:
        return summarize_multirow_output(rows)
    row = rows[0]
    visible = {k: v for k, v in row.items() if k != "participant_id" and v not in ("", None)}
    if len(visible) == 1:
        return next(iter(visible.values()))
    if visible:
        return "; ".join(f"{key}={value}" for key, value in visible.items())
    return "; ".join(f"{key}={value}" for key, value in row.items() if value)


def summarize_multirow_output(rows: list[dict[str, str]]) -> str:
    matched = 0
    unresolved = 0
    sample_parts: list[str] = []
    for row in rows:
        state = classify_output_row(row)
        if state == "matched":
            matched += 1
        else:
            unresolved += 1
        if len(sample_parts) < 3:
            name = row.get("name") or row.get("variant_name") or row.get("rsid") or "variant"
            genotype = normalize_evidence_value(row.get("genotype"))
            sample_parts.append(f"{name}={genotype}")
    summary = f"{len(rows)} rows; {matched} matched; {unresolved} unresolved"
    if sample_parts:
        summary += f"; sample: {', '.join(sample_parts)}"
    return summary


def classify_output_row(row: dict[str, str]) -> str:
    genotype = normalize_evidence_value(row.get("genotype"))
    if genotype != "MISSING":
        return "matched"
    matched_rsid = (row.get("matched_rsid") or "").strip()
    if matched_rsid:
        depth_text = (row.get("depth") or "").strip()
        if depth_text.isdigit() and int(depth_text) > 0:
            return "matched"
    return "unresolved"


def normalize_evidence_value(value: str | None) -> str:
    if value is None:
        return "MISSING"
    cleaned = str(value).strip()
    if cleaned in {"", "--", "not_found", "None"}:
        return "MISSING"
    return cleaned


def classify_variant_call(genotype: str | None, variant: VariantDefinition) -> str:
    call = normalize_evidence_value(genotype)
    if call == "MISSING":
        return "MISSING"
    ref = str(variant.fields.get("ref", "") or "").upper()
    alt = str(variant.fields.get("alt", "") or "").upper()
    kind = str(variant.fields.get("kind", "") or "").lower()
    normalized = call.upper().replace("/", "")
    if kind == "deletion":
        if "D" in normalized:
            return "VAR"
        if normalized and set(normalized) <= {"I", "-"}:
            return "REF"
        return normalized
    if alt and alt in normalized:
        return "VAR"
    if ref and normalized and set(normalized) <= set(ref):
        return "REF"
    return normalized


def format_rsid(rsid_value: Any) -> str:
    if isinstance(rsid_value, list):
        return "/".join(str(item) for item in rsid_value if item)
    return str(rsid_value or "")


def format_variant_target(variant: VariantDefinition) -> str:
    rsid = format_rsid(variant.fields.get("rsid"))
    ref = str(variant.fields.get("ref", "") or "")
    alt = str(variant.fields.get("alt", "") or "")
    kind = str(variant.fields.get("kind", "") or "")
    if ref or alt:
        change = f"{ref}>{alt}"
        if kind and kind != "snp":
            change = f"{change} ({kind})"
        return f"{rsid} {change}".strip()
    return rsid or variant.name


def format_info_display(variants: list[VariantDefinition]) -> str:
    if not variants:
        return ""
    return "; ".join(format_variant_target(variant) for variant in variants)


def format_evidence_display(genome_rows: list[dict[str, str]]) -> str:
    if not genome_rows:
        return ""
    parts = []
    for row in genome_rows:
        name = row.get("variant_name") or row.get("rsid") or "variant"
        evidence = normalize_evidence_value(row.get("genotype"))
        parts.append(f"{name}={evidence}")
    return "; ".join(parts)


def format_result_display(
    record_status: str,
    output_rows: list[dict[str, str]],
    genome_rows: list[dict[str, str]],
    variants: list[VariantDefinition],
) -> str:
    if output_rows:
        if len(output_rows) > 1:
            return summarize_multirow_output(output_rows)
        row = output_rows[0]
        visible = {k: v for k, v in row.items() if k != "participant_id" and v not in ("", None)}
        if len(visible) == 1:
            return str(next(iter(visible.values())))
    if genome_rows and variants:
        parts = []
        variant_by_name = {variant.name: variant for variant in variants}
        for row in genome_rows:
            name = row.get("variant_name") or row.get("rsid") or "variant"
            variant = variant_by_name.get(name)
            if variant is None:
                continue
            parts.append(f"{name}={classify_variant_call(row.get('genotype'), variant)}")
        if parts:
            return "; ".join(parts)
    if output_rows:
        return summarize_output(output_rows)
    if record_status == "skip":
        return "SKIPPED"
    return ""


def observation_variant_key(variant: VariantDefinition) -> str:
    return variant.manifest_name or variant.name


def variant_lookup_keys(variant: VariantDefinition) -> set[str]:
    keys = {variant.name, observation_variant_key(variant)}
    rsid = variant.fields.get("rsid")
    if isinstance(rsid, list):
        keys.update(str(item) for item in rsid if item)
    elif rsid:
        keys.add(str(rsid))
    if variant.path is not None:
        keys.add(variant.path.name)
        keys.add(variant.path.stem)
        keys.add(variant.path.as_posix())
        keys.add(relative_to(variant.path, REPO_ROOT))
    return {key for key in keys if key}


def build_variant_lookup(variants: list[VariantDefinition]) -> dict[str, VariantDefinition]:
    lookup: dict[str, VariantDefinition] = {}
    for variant in variants:
        for key in variant_lookup_keys(variant):
            lookup.setdefault(key, variant)
    return lookup


def variant_for_result_row(
    row: dict[str, str],
    variants: list[VariantDefinition],
    lookup: dict[str, VariantDefinition],
) -> VariantDefinition | None:
    for key in [
        row.get("name"),
        row.get("variant_name"),
        row.get("path"),
        Path(row.get("path", "")).name if row.get("path") else "",
        row.get("rsid"),
        row.get("matched_rsid"),
    ]:
        if key and key in lookup:
            return lookup[key]
    return variants[0] if len(variants) == 1 else None


def parse_locus(text: str | None) -> tuple[str, int, int] | None:
    if not text or ":" not in text:
        return None
    chrom, rest = text.split(":", 1)
    if "-" in rest:
        start_text, end_text = rest.split("-", 1)
    else:
        start_text = end_text = rest
    try:
        return chrom, int(start_text), int(end_text)
    except ValueError:
        return None


def preferred_locus(variant: VariantDefinition, row: dict[str, str]) -> tuple[str, str, int, int]:
    assembly = (row.get("assembly") or "").strip()
    if assembly.lower() == "grch37":
        preferred = [("GRCh37", variant.fields.get("grch37")), ("GRCh38", variant.fields.get("grch38"))]
    else:
        preferred = [("GRCh38", variant.fields.get("grch38")), ("GRCh37", variant.fields.get("grch37"))]
    for assembly_label, locus_text in preferred:
        locus = parse_locus(str(locus_text or ""))
        if locus:
            chrom, start, end = locus
            return assembly_label, chrom, start, end
    return assembly.upper() if assembly else "unknown", "", 0, 0


def observation_kind(variant: VariantDefinition, row: dict[str, str]) -> str:
    raw = (row.get("kind") or variant.fields.get("kind") or "").strip().lower()
    if raw == "variant":
        raw = str(variant.fields.get("kind") or "").strip().lower()
    return {
        "snv": "snp",
        "snp": "snp",
        "insertion": "ins",
        "ins": "ins",
        "deletion": "del",
        "del": "del",
        "indel": "hybrid",
    }.get(raw, raw or "snp")


def parse_int_or_none(value: Any) -> int | None:
    text = "" if value is None else str(value).strip()
    if not text:
        return None
    try:
        return int(text)
    except ValueError:
        return None


def normalize_observation_genotype(raw: str, ref: str, alt: str, kind: str) -> tuple[str, str, str]:
    display = normalize_evidence_value(raw)
    if display == "MISSING":
        return "./.", ".", "unknown"
    if display in {"0/0", "0|0"}:
        return "0/0", display, "hom_ref"
    if display in {"0/1", "1/0", "0|1", "1|0"}:
        return "0/1", display, "het"
    if display in {"1/1", "1|1"}:
        return "1/1", display, "hom_alt"

    cleaned = display.upper().replace("/", "").replace("|", "")
    ref_clean = ref.upper()
    alt_clean = alt.upper()
    if kind == "del" and cleaned == "ID":
        return "0/1", display, "het"
    if kind == "del" and cleaned == "DD":
        return "1/1", display, "hom_alt"
    if kind == "del" and cleaned in {"II", "--"}:
        return "0/0", display, "hom_ref"
    if ref_clean and cleaned == ref_clean * 2:
        return "0/0", display, "hom_ref"
    if alt_clean and cleaned == alt_clean * 2:
        return "1/1", display, "hom_alt"
    if ref_clean and alt_clean and ref_clean in cleaned and alt_clean in cleaned:
        return "0/1", display, "het"
    return display, display, "unknown"


def evidence_type_for(input_format: str, row: dict[str, str]) -> str:
    if not normalize_evidence_value(row.get("genotype")) == "MISSING":
        if input_format == "vcf":
            return "vcf_record"
        if input_format == "cram":
            return "mpileup"
    return ""


def facets_for_record(record: dict[str, Any]) -> str:
    return ""


def observation_from_row(
    row: dict[str, str],
    variant: VariantDefinition,
    record: dict[str, Any],
) -> dict[str, Any]:
    assembly, chrom, pos_start, pos_end = preferred_locus(variant, row)
    ref = str(row.get("ref") or variant.fields.get("ref") or "")
    alt = str(row.get("alt") or variant.fields.get("alt") or "")
    kind = observation_kind(variant, row)
    genotype, genotype_display, zygosity = normalize_observation_genotype(
        row.get("genotype", ""),
        ref,
        alt,
        kind,
    )
    ref_count = parse_int_or_none(row.get("ref_count"))
    alt_count = parse_int_or_none(row.get("alt_count"))
    depth = parse_int_or_none(row.get("depth"))
    genotype_missing = genotype == "./."
    matched = bool((row.get("matched_rsid") or row.get("rsid") or "").strip()) or not genotype_missing
    match_status = "found" if matched else "not_found"
    coverage_status = "unknown"
    if depth is not None:
        coverage_status = "covered" if depth > 0 else "not_covered"
    elif not genotype_missing:
        coverage_status = "covered"
    call_status = "called" if not genotype_missing else "no_call"
    if coverage_status == "not_covered":
        outcome = "not_covered"
    elif genotype_missing:
        outcome = "no_call"
    elif zygosity == "hom_ref":
        outcome = "reference"
    elif zygosity in {"het", "hom_alt"}:
        outcome = "variant"
    else:
        outcome = "unknown"
    evidence_raw = row.get("evidence") or ""
    if not evidence_raw and depth is not None:
        evidence_raw = f"{chrom}:{pos_start}-{pos_end} depth={depth}"
        if ref_count is not None:
            evidence_raw += f" ref_count={ref_count}"
        if alt_count is not None:
            evidence_raw += f" alt_count={alt_count}"
    rsid = row.get("matched_rsid") or row.get("rsid") or variant.fields.get("rsid")
    if isinstance(rsid, list):
        rsid = rsid[0] if rsid else ""
    return {
        "participant_id": record["participant"],
        "assay_id": record["assay"],
        "assay_version": "1.0",
        "variant_key": observation_variant_key(variant),
        "rsid": rsid or "",
        "assembly": assembly,
        "chrom": chrom,
        "pos_start": pos_start or "",
        "pos_end": pos_end or "",
        "ref": ref,
        "alt": alt,
        "kind": kind,
        "match_status": match_status,
        "coverage_status": coverage_status,
        "call_status": call_status,
        "genotype": genotype,
        "genotype_display": genotype_display,
        "zygosity": zygosity,
        "ref_count": ref_count if ref_count is not None else "",
        "alt_count": alt_count if alt_count is not None else "",
        "depth": depth if depth is not None else "",
        "genotype_quality": "",
        "allele_balance": (alt_count / depth) if alt_count is not None and depth else "",
        "outcome": outcome,
        "evidence_type": evidence_type_for(record["input_format"], row),
        "evidence_raw": evidence_raw,
        "facets": facets_for_record(record),
    }


def build_observations(
    record: dict[str, Any],
    output_rows: list[dict[str, str]],
    genome_rows: list[dict[str, str]],
    variants: list[VariantDefinition],
) -> list[dict[str, Any]]:
    lookup = build_variant_lookup(variants)
    source_rows = output_rows if output_rows and len(output_rows) > 1 else genome_rows
    observations = []
    seen_keys: set[str] = set()
    for row in source_rows:
        variant = variant_for_result_row(row, variants, lookup)
        if variant is None:
            continue
        observation = observation_from_row(row, variant, record)
        observations.append(observation)
        seen_keys.add(observation["variant_key"])
    if not observations and variants:
        for variant in variants:
            observation = observation_from_row({}, variant, record)
            observations.append(observation)
            seen_keys.add(observation["variant_key"])
    return observations


def write_observations(path: Path, observations: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OBSERVATION_HEADERS, delimiter="\t")
        writer.writeheader()
        for observation in observations:
            writer.writerow({key: observation.get(key, "") for key in OBSERVATION_HEADERS})


def read_observations(path: Path) -> list[dict[str, str]]:
    return read_tsv_rows(path)


def report_status_from(record: dict[str, Any], observations: list[dict[str, Any]]) -> str:
    if record["status"] == "fail":
        return "failed"
    if record["status"] == "skip":
        return "partial"
    if any(observation.get("outcome") in {"no_call", "not_covered", "unknown", "partial"} for observation in observations):
        return "partial"
    return "complete"


def input_profile(record: dict[str, Any]) -> dict[str, Any]:
    input_source = str(record.get("input_source", "") or "")
    input_format = str(record.get("input_format", "") or "")
    parts = [part for part in input_source.split("/") if part]
    provider = parts[0] if parts else ""
    source_version = parts[1] if len(parts) > 1 else ""
    source_label = provider
    version_label = source_version
    data_type = input_format

    if input_format in {"zip", "text"}:
        data_type = "snp_genotype"
    elif input_format in {"vcf", "cram"}:
        data_type = "sequencing"

    if provider == "23andme" and source_version.startswith("v") and source_version[1:].isdigit():
        source_label = "23andMe"
        version_label = f"chip version {source_version[1:]}"
    elif provider == "1k-genomes":
        source_label = "1000 Genomes"
        version_label = source_version
    elif provider:
        source_label = provider.replace("-", " ").title()

    display_parts = []
    if data_type:
        display_parts.append(data_type.replace("_", " "))
    if input_format:
        display_parts.append(input_format)
    if source_label:
        source_text = source_label
        if version_label:
            source_text += f" {version_label}"
        display_parts.append(source_text)

    return {
        "file_name": record.get("input_file", ""),
        "file_path": record.get("input_file_path", ""),
        "source": input_source,
        "format": input_format,
        "data_type": data_type,
        "provider": provider,
        "source_version": source_version,
        "source_label": source_label,
        "version_label": version_label,
        "display": " · ".join(part for part in display_parts if part),
    }


def build_report_object(
    record: dict[str, Any],
    observations: list[dict[str, Any]],
    started_at: str,
    finished_at: str,
    probe_duration_ms: int,
) -> dict[str, Any]:
    warnings = []
    if record.get("error"):
        warnings.append(record["error"])
    if record.get("genome_reads_error"):
        warnings.append(f"genome-reads probe: {record['genome_reads_error']}")
    n_variant = sum(1 for observation in observations if observation.get("outcome") == "variant")
    n_called = sum(1 for observation in observations if observation.get("call_status") == "called")
    fields = [
        {
            "key": "result",
            "label": "Result",
            "value": record.get("result_display") or record.get("result") or "",
            "value_type": "string",
            "format": "plain_text",
        },
        {
            "key": "summary",
            "label": "Summary",
            "value": summarize_output(read_tsv_rows(REPO_ROOT / record["output_file"])) if record.get("output_file") else "",
            "value_type": "string",
            "format": "plain_text",
        },
    ]
    return {
        "schema": "bioscript:report:1.0",
        "version": "1.0",
        "participant_id": record["participant"],
        "assay_id": record["assay"],
        "assay_version": "1.0",
        "input": input_profile(record),
        "report_status": report_status_from(record, observations),
        "derived_from": [observation["variant_key"] for observation in observations],
        "facets": facets_for_record(record) or None,
        "timing": {
            "started_at": started_at,
            "finished_at": finished_at,
            "duration_ms": record["duration_ms"],
            "breakdown": {
                "load_ms": None,
                "query_ms": record["duration_ms"],
                "compute_ms": probe_duration_ms if probe_duration_ms else None,
                "write_ms": None,
            },
        },
        "fields": fields,
        "metrics": {
            "n_variants": n_variant,
            "n_sites_tested": len(observations),
            "n_sites_called": n_called,
            "n_sites_missing": max(0, len(observations) - n_called),
        },
        "warnings": warnings or None,
        "notes": None,
    }


def write_aggregate_report_data(records: list[dict[str, Any]]) -> None:
    REPORT_DATA_DIR.mkdir(parents=True, exist_ok=True)
    all_observations: list[dict[str, Any]] = []
    reports: list[dict[str, Any]] = []
    for record in records:
        observations_file = record.get("observations_file")
        if observations_file:
            all_observations.extend(read_observations(REPO_ROOT / observations_file))
        report_json_file = record.get("report_json_file")
        if report_json_file:
            report_path = REPO_ROOT / report_json_file
            if report_path.exists():
                reports.append(json.loads(report_path.read_text(encoding="utf-8")))

    write_observations(AGGREGATE_OBSERVATIONS_TSV, all_observations)
    with AGGREGATE_REPORTS_JSONL.open("w", encoding="utf-8") as handle:
        for report in reports:
            handle.write(json.dumps(report, ensure_ascii=True, sort_keys=True))
            handle.write("\n")
    AGGREGATE_REPORTS_JSON.write_text(
        json.dumps(
            {
                "schema": "bioscript:report-set:1.0",
                "version": "1.0",
                "reports": reports,
            },
            indent=2,
        ),
        encoding="utf-8",
    )


def html_table(rows: list[dict[str, Any]], empty_message: str) -> str:
    if not rows:
        return f"<p>{html.escape(empty_message)}</p>"
    headers = list(rows[0].keys())
    out = ["<table>", "<thead><tr>"]
    out.extend(f"<th>{html.escape(str(header))}</th>" for header in headers)
    out.append("</tr></thead><tbody>")
    for row in rows:
        out.append("<tr>")
        for header in headers:
            value = row.get(header, "")
            text = "" if value is None else str(value)
            out.append(f"<td>{html.escape(text)}</td>")
        out.append("</tr>")
    out.append("</tbody></table>")
    return "".join(out)


def markdown_table(rows: list[dict[str, Any]], empty_message: str) -> str:
    if not rows:
        return empty_message
    headers = list(rows[0].keys())
    head = "| " + " | ".join(headers) + " |"
    sep = "| " + " | ".join("---" for _ in headers) + " |"
    body = []
    for row in rows:
        values = []
        for header in headers:
            value = row.get(header, "")
            text = "" if value is None else str(value)
            values.append(text.replace("|", "\\|").replace("\n", " "))
        body.append("| " + " | ".join(values) + " |")
    return "\n".join([head, sep, *body])


def read_json_file(path: Path) -> Any:
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def report_value_display(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=True, sort_keys=True)
    return str(value)


def render_report_html(report: dict[str, Any] | None) -> str:
    if not report:
        return "<p>No report JSON emitted.</p>"
    fields = report.get("fields") if isinstance(report.get("fields"), list) else []
    field_rows = [
        {
            "key": item.get("key", ""),
            "label": item.get("label", ""),
            "value": report_value_display(item.get("value")),
            "value_type": item.get("value_type", ""),
            "format": item.get("format", ""),
        }
        for item in fields
        if isinstance(item, dict)
    ]
    metrics = report.get("metrics") if isinstance(report.get("metrics"), dict) else {}
    metric_rows = [{"metric": key, "value": report_value_display(value)} for key, value in metrics.items()]
    timing = report.get("timing") if isinstance(report.get("timing"), dict) else {}
    breakdown = timing.get("breakdown") if isinstance(timing.get("breakdown"), dict) else {}
    timing_rows = [
        {"field": "started_at", "value": timing.get("started_at", "")},
        {"field": "finished_at", "value": timing.get("finished_at", "")},
        {"field": "duration_ms", "value": timing.get("duration_ms", "")},
        {"field": "load_ms", "value": breakdown.get("load_ms", "")},
        {"field": "query_ms", "value": breakdown.get("query_ms", "")},
        {"field": "compute_ms", "value": breakdown.get("compute_ms", "")},
        {"field": "write_ms", "value": breakdown.get("write_ms", "")},
    ]
    summary_rows = [
        {"field": "schema", "value": report.get("schema", "")},
        {"field": "version", "value": report.get("version", "")},
        {"field": "participant_id", "value": report.get("participant_id", "")},
        {"field": "assay_id", "value": report.get("assay_id", "")},
        {"field": "assay_version", "value": report.get("assay_version", "")},
        {"field": "report_status", "value": report.get("report_status", "")},
        {"field": "derived_from", "value": ", ".join(str(item) for item in report.get("derived_from", []) or [])},
        {"field": "facets", "value": report.get("facets", "") or ""},
    ]
    warning_rows = [{"warning": item} for item in report.get("warnings") or []]
    note_rows = [{"note": item} for item in report.get("notes") or []]
    raw_json = html.escape(json.dumps(report, indent=2, ensure_ascii=True, sort_keys=True))
    return "\n".join(
        [
            "<h2>Report</h2>",
            html_table(summary_rows, "No report summary."),
            "<h2>Fields</h2>",
            html_table(field_rows, "No report fields."),
            "<h2>Metrics</h2>",
            html_table(metric_rows, "No report metrics."),
            "<h2>Timing</h2>",
            html_table(timing_rows, "No report timing."),
            "<h2>Warnings</h2>",
            html_table(warning_rows, "No warnings."),
            "<h2>Notes</h2>",
            html_table(note_rows, "No notes."),
            "<h2>Raw JSON</h2>",
            f"<pre>{raw_json}</pre>",
        ]
    )


def write_results_csv(records: list[dict[str, Any]]) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "timestamp",
        "assay",
        "assay_path",
        "input_file",
        "input_source",
        "input_format",
        "participant",
        "status",
        "result",
        "evidence",
        "info",
        "duration_ms",
        "exit_code",
        "error",
        "output_file",
        "trace_file",
        "timing_file",
        "observations_file",
        "report_json_file",
        "details_file",
    ]
    with RESULTS_CSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for record in records:
            row = {name: record.get(name, "") for name in fieldnames}
            row["result"] = record.get("result_display", record.get("result", ""))
            row["evidence"] = record.get("evidence_display", "")
            row["info"] = record.get("info_display", "")
            writer.writerow(row)


def list_mode(assays: list[Path], inputs: list[Path]) -> int:
    print("\nAssays:")
    for assay in assays:
        print(f"  {relative_to(assay, REPO_ROOT)}")
    print("\nTest data:")
    for input_path in inputs:
        print(f"  {relative_to(input_path, REPO_ROOT)}")
    return 0


def maybe_reference_args(input_path: Path, runtime_root: Path) -> tuple[list[str], str | None]:
    if not needs_reference(input_path):
        return [], None
    reference = find_reference_for(input_path) or (DEFAULT_REF_FA if DEFAULT_REF_FA.exists() else None)
    if reference is None:
        return [], "reference genome not available"
    args = ["--reference-file", relative_to(reference, runtime_root)]
    reference_index = reference.with_suffix(reference.suffix + ".fai")
    if reference_index.exists():
        args.extend(["--reference-index", relative_to(reference_index, runtime_root)])
    cram_index = Path(f"{input_path}.crai")
    if cram_index.exists():
        args.extend(["--input-index", relative_to(cram_index, runtime_root)])
    tbi_index = Path(f"{input_path}.tbi")
    if tbi_index.exists():
        args.extend(["--input-index", relative_to(tbi_index, runtime_root)])
    return args, None


def runtime_root_for(assay_path: Path, input_path: Path, output_path: Path) -> Path:
    roots = [str(input_path.parent), str(output_path.parent), str(WORKSPACE_ROOT)]
    if str(assay_path).startswith(str(WORKSPACE_ROOT)):
        roots.append(str(assay_path.parent))
    return Path(os.path.commonpath(roots))


def run_subprocess(cmd: list[str], timeout_seconds: int = HARD_TIMEOUT_SECONDS) -> tuple[int, str, str, int, bool]:
    start = time.time_ns()
    try:
        completed = subprocess.run(
            cmd,
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=False,
            timeout=timeout_seconds,
        )
        duration_ms = int((time.time_ns() - start) / 1_000_000)
        return completed.returncode, completed.stdout, completed.stderr, duration_ms, False
    except subprocess.TimeoutExpired as exc:
        duration_ms = int((time.time_ns() - start) / 1_000_000)
        stdout = exc.stdout or ""
        stderr = exc.stderr or ""
        stderr = (stderr + f"\nprocess exceeded {timeout_seconds}s hard timeout").strip()
        return 124, stdout, stderr, duration_ms, True


def open_genotype_text(input_path: Path):
    if detect_format(input_path) == "zip":
        archive = zipfile.ZipFile(input_path)
        entries = [info for info in archive.infolist() if not info.is_dir()]
        entries.sort(key=lambda info: (not info.filename.lower().endswith((".txt", ".csv", ".tsv")), info.filename))
        if not entries:
            archive.close()
            raise RuntimeError(f"zip archive {input_path} does not contain a genotype entry")
        handle = archive.open(entries[0], "r")
        text_handle = open_text_wrapper(handle)
        return archive, text_handle
    return None, input_path.open("r", encoding="utf-8-sig", errors="replace", newline="")


def open_text_wrapper(handle):
    return io.TextIOWrapper(handle, encoding="utf-8-sig", errors="replace", newline="")


def detect_text_delimiter(lines: list[str]) -> str:
    for line in lines:
        trimmed = line.strip()
        if not trimmed:
            continue
        if "\t" in line:
            return "\t"
        if "," in line:
            return ","
        if len(trimmed.split()) > 1:
            return " "
    return "\t"


def split_fields(line: str, delimiter: str) -> list[str]:
    if delimiter == "\t":
        return [field.strip() for field in line.split("\t")]
    if delimiter == ",":
        return next(csv.reader([line]))
    return line.split()


def normalize_name(name: str) -> str:
    return "".join(ch for ch in name.strip().lower() if ch not in " _-")


def looks_like_header(fields: list[str]) -> bool:
    if not fields:
        return False
    return normalize_name(fields[0]) in {"rsid", "name", "snp", "marker", "id", "snpid"}


def genotype_from_row(fields: list[str], header: list[str]) -> tuple[str | None, str | None]:
    row = {normalize_name(key): value.strip() for key, value in zip(header, fields)}
    rsid = next((row.get(key) for key in ("rsid", "name", "snp", "marker", "id", "snpid") if row.get(key)), None)
    genotype = next(
        (row.get(key) for key in ("genotype", "gt", "result", "results", "call", "calls", "yourcode", "code", "genotypevalue", "variation") if row.get(key)),
        None,
    )
    if genotype is None:
        allele1 = next((row.get(key) for key in ("allele1", "allelea", "allele_a", "allele1top") if row.get(key)), "")
        allele2 = next((row.get(key) for key in ("allele2", "alleleb", "allele_b", "allele2top") if row.get(key)), "")
        genotype = f"{allele1}{allele2}" if (allele1 or allele2) else None
    return rsid, genotype


def normalize_genotype(value: str | None) -> str:
    if not value:
        return ""
    cleaned = value.strip().replace(" ", "").upper()
    if cleaned in {"", "NA", "N/A", "#N/A", "NONE"}:
        return "--"
    if "/" in cleaned:
        parts = cleaned.split("/")
        if any(part in {"", "-"} for part in parts):
            return "ID"
        return "".join(parts)
    return cleaned


def extract_genome_reads_direct(input_path: Path, variants: list[VariantDefinition], output_path: Path) -> tuple[bool, str]:
    if detect_format(input_path) not in {"zip", "text"}:
        return False, "direct streaming genome-read extraction only supports zip/text genotype inputs"

    targets: dict[str, list[int]] = {}
    resolved: dict[int, str] = {}
    for idx, variant in enumerate(variants):
        rsid_value = variant.fields.get("rsid")
        rsids = rsid_value if isinstance(rsid_value, list) else [rsid_value]
        for rsid in [item for item in rsids if isinstance(item, str) and item]:
            targets.setdefault(rsid, []).append(idx)

    archive = None
    handle = None
    try:
        archive, handle = open_genotype_text(input_path)
        buffered: list[str] = []
        header: list[str] | None = None
        delimiter = "\t"
        for _ in range(8):
            line = handle.readline()
            if not line:
                break
            buffered.append(line.rstrip("\n\r"))
        delimiter = detect_text_delimiter(buffered)
        for line in buffered + [raw.rstrip("\n\r") for raw in handle]:
            trimmed = line.strip()
            if not trimmed:
                continue
            if trimmed.startswith(("#", "//")):
                candidate = trimmed.lstrip("#/ ").strip()
                fields = split_fields(candidate, delimiter)
                if looks_like_header(fields):
                    header = fields
                continue
            fields = split_fields(line, delimiter)
            if not fields:
                continue
            if header is None:
                if looks_like_header(fields):
                    header = fields
                    continue
                header = ["rsid", "chromosome", "position", "genotype"][: len(fields)]
            rsid, genotype = genotype_from_row(fields, header)
            if not rsid or rsid not in targets:
                continue
            normalized = normalize_genotype(genotype)
            for variant_idx in targets[rsid]:
                resolved.setdefault(variant_idx, normalized)
            if len(resolved) == len(variants):
                break
    except Exception as exc:  # noqa: BLE001
        return False, str(exc)
    finally:
        if handle is not None:
            handle.close()
        if archive is not None:
            archive.close()

    rows = []
    for idx, variant in enumerate(variants):
        rows.append(
            {
                "variant_name": variant.name,
                "rsid": json.dumps(variant.fields.get("rsid")) if isinstance(variant.fields.get("rsid"), list) else str(variant.fields.get("rsid", "")),
                "kind": str(variant.fields.get("kind", "")),
                "grch37": str(variant.fields.get("grch37", "")),
                "grch38": str(variant.fields.get("grch38", "")),
                "ref": str(variant.fields.get("ref", "")),
                "alt": str(variant.fields.get("alt", "")),
                "genotype": resolved.get(idx, ""),
            }
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle_out:
        writer = csv.DictWriter(handle_out, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return True, ""


def run_probe(
    assay_path: Path,
    input_path: Path,
    participant: str,
    variants: list[VariantDefinition],
    output_path: Path,
) -> tuple[bool, str]:
    if not variants:
        return False, "assay does not declare bioscript.variant definitions"
    direct_ok, direct_error = extract_genome_reads_direct(input_path, variants, output_path)
    if direct_ok:
        return True, ""
    if detect_format(input_path) in {"zip", "text"}:
        return False, direct_error
    probe_script = ensure_probe_script(assay_path, variants)
    runtime_root = runtime_root_for(assay_path, input_path, output_path)
    cmd = [
        str(BIOSCRIPT_BINARY),
        str(probe_script),
        "--root",
        str(runtime_root),
        "--input-file",
        relative_to(input_path, runtime_root),
        "--output-file",
        relative_to(output_path, runtime_root),
        "--participant-id",
        participant,
        "--auto-index",
        "--max-duration-ms",
        "30000",
        "--max-memory-bytes",
        str(64 * 1024 * 1024),
    ]
    reference_args, reference_error = maybe_reference_args(input_path, runtime_root)
    if reference_error:
        return False, reference_error
    cmd.extend(reference_args)
    code, stdout, stderr, _duration, _timed_out = run_subprocess(cmd, timeout_seconds=45)
    if code != 0:
        summary = stderr.strip().splitlines()[-1] if stderr.strip() else stdout.strip().splitlines()[-1] if stdout.strip() else "probe failed"
        return False, summary
    return True, ""


def render_detail_html(record: dict[str, Any]) -> str:
    summary_rows = [
        {"field": "status", "value": record["status"]},
        {"field": "result", "value": record["result_display"]},
        {"field": "evidence", "value": record.get("evidence_display", "")},
        {"field": "info", "value": record.get("info_display", "")},
        {"field": "assay", "value": record["assay_path"]},
        {"field": "assay_sha256", "value": record["assay_sha256"]},
        {"field": "git_commit", "value": record["git_commit_short"] or record["git_commit"]},
        {"field": "git_status", "value": record["git_status"] or "clean"},
        {"field": "ran_at_utc", "value": record["timestamp"]},
        {"field": "input_file", "value": record["input_file_path"]},
        {"field": "input_format", "value": record["input_format"]},
        {"field": "output_file", "value": record["output_file"]},
        {"field": "trace_file", "value": record["trace_file"]},
        {"field": "timing_file", "value": record["timing_file"]},
        {"field": "observations_file", "value": record.get("observations_file", "")},
        {"field": "report_json_file", "value": record.get("report_json_file", "")},
    ]
    trace_rows = read_tsv_rows(REPO_ROOT / record["trace_file"]) if record["trace_file"] else []
    timing_rows = read_tsv_rows(REPO_ROOT / record["timing_file"]) if record["timing_file"] else []
    genome_rows = read_tsv_rows(REPO_ROOT / record["genome_reads_file"]) if record["genome_reads_file"] else []
    observation_rows = read_tsv_rows(REPO_ROOT / record["observations_file"]) if record.get("observations_file") else []
    output_rows = read_tsv_rows(REPO_ROOT / record["output_file"]) if record["output_file"] else []
    report_data = read_json_file(REPO_ROOT / record["report_json_file"]) if record.get("report_json_file") else None
    summary_table = html_table(summary_rows, "No summary available.")
    report_table = render_report_html(report_data)
    genome_table = html_table(genome_rows, "No genome reads captured.")
    observation_table = html_table(observation_rows, "No observations emitted.")
    output_table = html_table(output_rows, "No assay output.")
    trace_table = html_table(trace_rows, "No trace available.")
    timing_table = html_table(timing_rows, "No timing data available.")
    stdout_text = html.escape(record.get("stdout", ""))
    stderr_text = html.escape(record.get("stderr", ""))
    command = html.escape(" ".join(shlex.quote(part) for part in record["command"]))
    title = html.escape(f"{record['assay']} / {record['participant']}")
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{title}</title>
  <style>
    :root {{
      --bg: #f5f1e8;
      --panel: #fffdf8;
      --ink: #1c1b19;
      --muted: #6d685f;
      --line: #d8cfbf;
      --accent: #0b6e4f;
      --warn: #9f2a2a;
      --shadow: 0 18px 40px rgba(40, 28, 8, 0.08);
    }}
    body {{
      margin: 0;
      font-family: "Iowan Old Style", "Palatino Linotype", "Book Antiqua", serif;
      color: var(--ink);
      background:
        radial-gradient(circle at top left, rgba(11, 110, 79, 0.08), transparent 28rem),
        linear-gradient(180deg, #fbf7ef 0%, var(--bg) 100%);
    }}
    main {{ max-width: 1200px; margin: 0 auto; padding: 2rem; }}
    .card {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 18px;
      box-shadow: var(--shadow);
      padding: 1.25rem;
      margin-bottom: 1.25rem;
    }}
    h1, h2 {{ margin: 0 0 0.8rem; }}
    p.meta {{ color: var(--muted); margin-top: 0.35rem; }}
    .tabs {{ display: flex; flex-wrap: wrap; gap: 0.5rem; margin-bottom: 1rem; }}
    .tabs button {{
      border: 1px solid var(--line);
      background: #efe7d9;
      color: var(--ink);
      padding: 0.55rem 0.85rem;
      border-radius: 999px;
      cursor: pointer;
    }}
    .tabs button.active {{ background: var(--accent); color: white; border-color: var(--accent); }}
    .panel {{ display: none; }}
    .panel.active {{ display: block; }}
    table {{
      width: 100%;
      border-collapse: collapse;
      font-family: ui-monospace, "SFMono-Regular", Menlo, monospace;
      font-size: 0.9rem;
    }}
    th, td {{ border-bottom: 1px solid var(--line); text-align: left; vertical-align: top; padding: 0.5rem; }}
    th {{ background: #f2ebde; }}
    pre {{
      margin: 0;
      padding: 1rem;
      background: #171512;
      color: #f5f0e6;
      border-radius: 14px;
      overflow: auto;
      font-family: ui-monospace, "SFMono-Regular", Menlo, monospace;
      font-size: 0.85rem;
    }}
    a {{ color: var(--accent); }}
  </style>
</head>
<body>
  <main>
    <div class="card">
      <h1>{title}</h1>
      <p class="meta">{html.escape(record['status'].upper())} in {html.escape(record['duration_display'])} on {html.escape(record['timestamp'])}</p>
      <p><a href="{html.escape(record['index_href'])}">Back to report index</a></p>
    </div>
    <div class="card">
      <div class="tabs">
        <button class="tab active" data-target="summary">Summary</button>
        <button class="tab" data-target="report">Report</button>
        <button class="tab" data-target="genome">Genome Reads</button>
        <button class="tab" data-target="observations">Observations</button>
        <button class="tab" data-target="output">Assay Output</button>
        <button class="tab" data-target="trace">Trace</button>
        <button class="tab" data-target="timings">Timings</button>
        <button class="tab" data-target="logs">Logs</button>
      </div>
      <section id="summary" class="panel active">{summary_table}</section>
      <section id="report" class="panel">{report_table}</section>
      <section id="genome" class="panel">{genome_table}</section>
      <section id="observations" class="panel">{observation_table}</section>
      <section id="output" class="panel">{output_table}</section>
      <section id="trace" class="panel">{trace_table}</section>
      <section id="timings" class="panel">{timing_table}</section>
      <section id="logs" class="panel">
        <h2>Command</h2>
        <pre>{command}</pre>
        <h2>stdout</h2>
        <pre>{stdout_text}</pre>
        <h2>stderr</h2>
        <pre>{stderr_text}</pre>
      </section>
    </div>
  </main>
  <script>
    const tabs = document.querySelectorAll('.tab');
    const panels = document.querySelectorAll('.panel');
    tabs.forEach((tab) => {{
      tab.addEventListener('click', () => {{
        tabs.forEach((item) => item.classList.remove('active'));
        panels.forEach((item) => item.classList.remove('active'));
        tab.classList.add('active');
        document.getElementById(tab.dataset.target).classList.add('active');
      }});
    }});
  </script>
</body>
</html>
"""


def render_detail_markdown(record: dict[str, Any]) -> str:
    trace_rows = read_tsv_rows(REPO_ROOT / record["trace_file"]) if record["trace_file"] else []
    timing_rows = read_tsv_rows(REPO_ROOT / record["timing_file"]) if record["timing_file"] else []
    genome_rows = read_tsv_rows(REPO_ROOT / record["genome_reads_file"]) if record["genome_reads_file"] else []
    observation_rows = read_tsv_rows(REPO_ROOT / record["observations_file"]) if record.get("observations_file") else []
    output_rows = read_tsv_rows(REPO_ROOT / record["output_file"]) if record["output_file"] else []
    lines = [
        f"# {record['assay']} / {record['participant']}",
        "",
        f"- Status: `{record['status']}`",
        f"- Result: `{record['result_display']}`",
        f"- Evidence: `{record.get('evidence_display', '')}`",
        f"- Info: `{record.get('info_display', '')}`",
        f"- Duration: `{record['duration_display']}`",
        f"- Ran at: `{record['timestamp']}`",
        f"- Input file: `{record['input_file_path']}`",
        f"- Assay file: `{record['assay_path']}`",
        f"- Assay SHA-256: `{record['assay_sha256']}`",
        f"- Git commit: `{record['git_commit_short'] or record['git_commit']}`",
        f"- Git status: `{record['git_status'] or 'clean'}`",
        f"- Trace file: `{record['trace_file']}`",
        f"- Timing file: `{record['timing_file']}`",
        f"- Observations file: `{record.get('observations_file', '')}`",
        f"- Report JSON file: `{record.get('report_json_file', '')}`",
        "",
        "## Genome Reads",
        "",
        markdown_table(genome_rows, "No genome reads captured."),
        "",
        "## Observations",
        "",
        markdown_table(observation_rows, "No observations emitted."),
        "",
        "## Assay Output",
        "",
        markdown_table(output_rows, "No assay output."),
        "",
        "## Trace",
        "",
        markdown_table(trace_rows[:50], "No trace available."),
        "",
        "## Timings",
        "",
        markdown_table(timing_rows, "No timing data available."),
        "",
        "## stdout",
        "",
        "```text",
        record.get("stdout", ""),
        "```",
        "",
        "## stderr",
        "",
        "```text",
        record.get("stderr", ""),
        "```",
    ]
    return "\n".join(lines)


def render_call_cell(text: str) -> str:
    """Render a Result/Evidence cell, colorizing REF / VAR / MISSING tokens."""
    if not text:
        return ""
    parts = [part.strip() for part in text.split(";")]
    rendered = []
    for part in parts:
        if not part:
            continue
        key, sep, value = part.partition("=")
        if not sep:
            rendered.append(html.escape(part))
            continue
        value_stripped = value.strip()
        token_class = None
        upper = value_stripped.upper()
        if upper == "REF":
            token_class = "token-ref"
        elif upper == "VAR":
            token_class = "token-var"
        elif upper == "MISSING":
            token_class = "token-missing"
        value_html = html.escape(value_stripped)
        if token_class:
            value_html = f'<span class="{token_class}">{value_html}</span>'
        rendered.append(f"{html.escape(key.strip())}={value_html}")
    return "; ".join(rendered)


def read_jsonl_file(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        value = json.loads(line)
        if isinstance(value, dict):
            rows.append(value)
    return rows


def render_report_field_item(field: dict[str, Any]) -> str:
    key = html.escape(str(field.get("key", "")))
    label = html.escape(str(field.get("label", "") or field.get("key", "")))
    value = html.escape(report_value_display(field.get("value")))
    value_type = html.escape(str(field.get("value_type", "")))
    field_format = str(field.get("format", "plain_text"))
    css = "report-field report-field-badge" if field_format == "badge" else "report-field"
    return (
        f'<div class="{css}">'
        f'<div class="report-field-label">{label}</div>'
        f'<div class="report-field-value">{value}</div>'
        f'<div class="report-field-meta">{key} · {value_type} · {html.escape(field_format)}</div>'
        "</div>"
    )


def render_report_card(report: dict[str, Any]) -> str:
    fields = report.get("fields") if isinstance(report.get("fields"), list) else []
    metrics = report.get("metrics") if isinstance(report.get("metrics"), dict) else {}
    warnings = report.get("warnings") or []
    notes = report.get("notes") or []
    timing = report.get("timing") if isinstance(report.get("timing"), dict) else {}
    input_info = report.get("input") if isinstance(report.get("input"), dict) else {}
    field_items = "".join(
        render_report_field_item(field)
        for field in fields
        if isinstance(field, dict)
    ) or '<p class="meta">No report fields.</p>'
    metric_items = "".join(
        f'<span class="metric"><strong>{html.escape(str(key))}</strong> {html.escape(report_value_display(value))}</span>'
        for key, value in metrics.items()
    )
    warning_items = "".join(f"<li>{html.escape(str(item))}</li>" for item in warnings)
    note_items = "".join(f"<li>{html.escape(str(item))}</li>" for item in notes)
    derived = ", ".join(str(item) for item in report.get("derived_from", []) or [])
    status = str(report.get("report_status", ""))
    status_class = f"report-status-{status.replace('_', '-')}" if status else ""
    input_display = str(input_info.get("display") or "")
    input_file = str(input_info.get("file_name") or "")
    input_source = str(input_info.get("source") or "")
    raw_json = html.escape(json.dumps(report, indent=2, ensure_ascii=True, sort_keys=True))
    return f"""
      <article class="report-card">
        <div class="report-card-head">
          <div>
            <h3>{html.escape(str(report.get('participant_id', '')))}</h3>
            <p class="meta">{html.escape(str(report.get('assay_id', '')))} · version {html.escape(str(report.get('assay_version', '')))}</p>
          </div>
          <span class="report-status {status_class}">{html.escape(status)}</span>
        </div>
        <div class="input-profile">
          <span>{html.escape(input_display or 'unknown input')}</span>
          <small>{html.escape(input_file)}{html.escape(' · ' + input_source if input_source else '')}</small>
        </div>
        <div class="report-fields">{field_items}</div>
        <div class="metrics">{metric_items}</div>
        <dl class="report-details">
          <dt>Derived from</dt><dd>{html.escape(derived)}</dd>
          <dt>Facets</dt><dd>{html.escape(str(report.get('facets') or ''))}</dd>
          <dt>Timing</dt><dd>{html.escape(str(timing.get('duration_ms', '')))} ms</dd>
        </dl>
        {'<h4>Warnings</h4><ul>' + warning_items + '</ul>' if warning_items else ''}
        {'<h4>Notes</h4><ul>' + note_items + '</ul>' if note_items else ''}
        <details>
          <summary>Raw report JSON</summary>
          <pre>{raw_json}</pre>
        </details>
      </article>
    """


def rsid_link_html(rsid: str) -> str:
    rsid = rsid.strip()
    if not rsid.startswith("rs"):
        return html.escape(rsid)
    href = f"https://www.ncbi.nlm.nih.gov/snp/{html.escape(rsid)}"
    return f'<a href="{href}" target="_blank" rel="noopener noreferrer">{html.escape(rsid)}</a>'


def highlighted_genotype_html(genotype_display: str, alt: str) -> str:
    genotype_display = str(genotype_display or "")
    alt = str(alt or "")
    if not genotype_display:
        return ""
    if not alt or alt.startswith("<") or len(alt) > len(genotype_display):
        return html.escape(genotype_display)
    rendered: list[str] = []
    index = 0
    while index < len(genotype_display):
        if genotype_display[index : index + len(alt)].upper() == alt.upper():
            rendered.append(f'<span class="alt-hit">{html.escape(genotype_display[index:index + len(alt)])}</span>')
            index += len(alt)
        else:
            rendered.append(html.escape(genotype_display[index]))
            index += 1
    return "".join(rendered)


def observation_cell_html(observation: dict[str, str], header: str, outcome_class: str) -> str:
    value = str(observation.get(header, ""))
    classes = [f"col-{header.replace('_', '-')}"]
    if header == "outcome" and outcome_class:
        classes.append(outcome_class)
    class_attr = f' class="{" ".join(classes)}"'
    if header == "rsid":
        content = rsid_link_html(value)
    elif header == "genotype_display":
        content = highlighted_genotype_html(value, str(observation.get("alt", "")))
    else:
        content = html.escape(value)
    return f"<td{class_attr}>{content}</td>"


def render_index_html(records: list[dict[str, Any]]) -> str:
    observation_rows = read_tsv_rows(AGGREGATE_OBSERVATIONS_TSV)
    report_rows = read_jsonl_file(AGGREGATE_REPORTS_JSONL)
    observation_html_rows = []
    observation_headers = [
        "participant_id",
        "rsid",
        "ref",
        "alt",
        "genotype_display",
        "outcome",
        "assay_id",
        "variant_key",
        "kind",
        "genotype",
        "zygosity",
        "assembly",
        "chrom",
        "pos_start",
        "pos_end",
        "depth",
        "allele_balance",
        "facets",
    ]
    for observation in observation_rows:
        outcome = observation.get("outcome", "")
        outcome_class = f"outcome-{outcome.replace('_', '-')}" if outcome else ""
        observation_html_rows.append(
            "<tr>"
            + "".join(
                observation_cell_html(observation, header, outcome_class)
                for header in observation_headers
            )
            + "</tr>"
        )
    observation_header_cells = "".join(
        f'<th class="col-{header.replace("_", "-")}"><button class="sort-button" data-column="{idx}" type="button">{html.escape(header)}</button></th>'
        for idx, header in enumerate(observation_headers)
    )
    report_cards = "".join(render_report_card(report) for report in report_rows)
    total = len(records)
    passed = sum(record["status"] == "pass" for record in records)
    failed = sum(record["status"] == "fail" for record in records)
    skipped = sum(record["status"] == "skip" for record in records)
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Assay Test Report</title>
  <style>
    :root {{
      --bg: #efe7da;
      --paper: #fffdf8;
      --ink: #181614;
      --muted: #625c53;
      --line: #d4c6b0;
      --accent: #144f8a;
      --pass: #1b7f5d;
      --fail: #9e2f2f;
      --skip: #8b6b21;
      --shadow: 0 24px 60px rgba(20, 18, 14, 0.12);
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      color: var(--ink);
      background:
        linear-gradient(130deg, rgba(20, 79, 138, 0.10), transparent 42%),
        radial-gradient(circle at top right, rgba(27, 127, 93, 0.12), transparent 26rem),
        var(--bg);
      font-family: Georgia, "Times New Roman", serif;
    }}
    main {{ max-width: 1320px; margin: 0 auto; padding: 2rem; }}
    .hero, .panel {{
      background: var(--paper);
      border: 1px solid var(--line);
      border-radius: 24px;
      box-shadow: var(--shadow);
    }}
    .hero {{ padding: 1.5rem; margin-bottom: 1.25rem; }}
    .stats {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
      gap: 0.75rem;
      margin-top: 1rem;
    }}
    .stat {{
      background: linear-gradient(180deg, rgba(255,255,255,0.9), rgba(243,236,224,0.95));
      border: 1px solid var(--line);
      border-radius: 18px;
      padding: 0.9rem 1rem;
    }}
    .stat strong {{ display: block; font-size: 1.6rem; }}
    .links {{ display: flex; flex-wrap: wrap; gap: 0.75rem; margin-top: 1rem; }}
    .links a {{
      display: inline-block;
      border: 1px solid var(--line);
      border-radius: 999px;
      background: #fffdfa;
      padding: 0.55rem 0.85rem;
      text-decoration: none;
    }}
    .view-tabs {{
      display: flex;
      flex-wrap: wrap;
      gap: 0.5rem;
      margin-bottom: 1rem;
    }}
    .view-tab {{
      border: 1px solid var(--line);
      border-radius: 999px;
      background: #fffdfa;
      color: var(--ink);
      cursor: pointer;
      font: inherit;
      padding: 0.65rem 0.95rem;
    }}
    .view-tab.active {{
      background: var(--accent);
      border-color: var(--accent);
      color: white;
    }}
    .view-panel {{ display: none; }}
    .view-panel.active {{ display: block; }}
    .panel {{ padding: 1rem; }}
    input {{
      width: 100%;
      border: 1px solid var(--line);
      border-radius: 999px;
      padding: 0.9rem 1rem;
      font-size: 1rem;
      margin-bottom: 1rem;
      background: #fffdfa;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      font-family: ui-monospace, "SFMono-Regular", Menlo, monospace;
      font-size: 0.9rem;
    }}
    th, td {{
      text-align: left;
      padding: 0.55rem 0.45rem;
      border-bottom: 1px solid var(--line);
      vertical-align: top;
    }}
    th {{ background: #f3ebdd; position: sticky; top: 0; }}
    .sort-button {{
      appearance: none;
      border: 0;
      background: transparent;
      color: inherit;
      cursor: pointer;
      display: inline-flex;
      align-items: center;
      gap: 0.2rem;
      font: inherit;
      font-weight: 700;
      padding: 0;
      text-align: left;
      width: 100%;
    }}
    .sort-button::after {{
      color: var(--muted);
      content: "";
      font-size: 0.85rem;
      line-height: 1;
    }}
    .sort-button.sorted-asc::after {{ content: "^"; }}
    .sort-button.sorted-desc::after {{ content: "v"; }}
    a {{ color: var(--accent); }}
    .table-wrap {{ max-height: 72vh; overflow: auto; border-radius: 16px; }}
    .col-participant-id {{
      font-size: 0.74rem;
      line-height: 1.2;
      max-width: 8.5rem;
      overflow-wrap: anywhere;
    }}
    .col-rsid {{ white-space: nowrap; }}
    .col-ref, .col-alt {{
      max-width: 4.5rem;
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
      width: 1%;
    }}
    .col-genotype-display {{
      font-weight: 700;
      white-space: nowrap;
    }}
    .alt-hit {{
      background: rgba(209, 120, 35, 0.28);
      border-radius: 4px;
      color: #8f4306;
      display: inline-block;
      padding: 0 0.22em;
    }}
    .status-pass {{ color: var(--pass); font-weight: 700; }}
    .status-fail {{ color: var(--fail); font-weight: 700; }}
    .status-skip {{ color: var(--skip); font-weight: 700; }}
    .token-ref {{
      background: rgba(27, 127, 93, 0.16);
      color: #0f6347;
      border-radius: 4px;
      padding: 0 0.35em;
      font-weight: 700;
    }}
    .token-var {{
      background: rgba(209, 120, 35, 0.18);
      color: #a65412;
      border-radius: 4px;
      padding: 0 0.35em;
      font-weight: 700;
    }}
    .token-missing {{
      background: rgba(120, 113, 101, 0.12);
      color: #8a8277;
      border-radius: 4px;
      padding: 0 0.35em;
    }}
    .outcome-variant {{ color: #a65412; font-weight: 700; }}
    .outcome-reference {{ color: #0f6347; font-weight: 700; }}
    .outcome-no-call, .outcome-not-covered, .outcome-unknown {{
      color: #8a8277;
      font-weight: 700;
    }}
    .report-grid {{
      display: grid;
      gap: 1rem;
    }}
    .report-card {{
      background: #fffdfa;
      border: 1px solid var(--line);
      border-radius: 16px;
      padding: 1rem;
    }}
    .report-card h3 {{ margin: 0; }}
    .report-card h4 {{ margin: 1rem 0 0.4rem; }}
    .report-card-head {{
      align-items: flex-start;
      display: flex;
      gap: 1rem;
      justify-content: space-between;
    }}
    .report-status {{
      border: 1px solid var(--line);
      border-radius: 999px;
      font-family: ui-monospace, "SFMono-Regular", Menlo, monospace;
      font-weight: 700;
      padding: 0.35rem 0.65rem;
    }}
    .report-status-complete {{ color: var(--pass); }}
    .report-status-partial {{ color: var(--skip); }}
    .report-status-failed {{ color: var(--fail); }}
    .input-profile {{
      align-items: baseline;
      display: flex;
      flex-wrap: wrap;
      gap: 0.4rem 0.75rem;
      margin-top: 0.65rem;
    }}
    .input-profile span {{
      background: #f3ebdd;
      border-radius: 999px;
      font-family: ui-monospace, "SFMono-Regular", Menlo, monospace;
      font-size: 0.86rem;
      padding: 0.3rem 0.6rem;
    }}
    .input-profile small {{
      color: var(--muted);
      font-family: ui-monospace, "SFMono-Regular", Menlo, monospace;
      overflow-wrap: anywhere;
    }}
    .report-fields {{
      display: grid;
      gap: 0.75rem;
      grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
      margin: 1rem 0;
    }}
    .report-field {{
      border: 1px solid var(--line);
      border-radius: 12px;
      padding: 0.75rem;
    }}
    .report-field-badge .report-field-value {{
      background: rgba(20, 79, 138, 0.12);
      border-radius: 999px;
      display: inline-block;
      font-weight: 700;
      margin-top: 0.35rem;
      padding: 0.3rem 0.65rem;
    }}
    .report-field-label {{ color: var(--muted); font-size: 0.9rem; }}
    .report-field-value {{
      font-family: ui-monospace, "SFMono-Regular", Menlo, monospace;
      margin-top: 0.25rem;
      overflow-wrap: anywhere;
    }}
    .report-field-meta {{
      color: var(--muted);
      font-family: ui-monospace, "SFMono-Regular", Menlo, monospace;
      font-size: 0.75rem;
      margin-top: 0.4rem;
    }}
    .metrics {{
      display: flex;
      flex-wrap: wrap;
      gap: 0.5rem;
      margin-bottom: 1rem;
    }}
    .metric {{
      background: #f3ebdd;
      border-radius: 999px;
      padding: 0.35rem 0.65rem;
    }}
    .report-details {{
      display: grid;
      gap: 0.35rem 1rem;
      grid-template-columns: max-content 1fr;
      margin: 0;
    }}
    .report-details dt {{ color: var(--muted); }}
    .report-details dd {{ margin: 0; overflow-wrap: anywhere; }}
    p.meta {{ color: var(--muted); }}
  </style>
</head>
<body>
  <main>
    <section class="hero">
      <h1>Local Assay Test Report</h1>
      <p class="meta">Generated {html.escape(datetime.now(timezone.utc).isoformat())} from {html.escape(str(REPO_ROOT))}</p>
      <div class="stats">
        <div class="stat"><span>Total runs</span><strong>{total}</strong></div>
        <div class="stat"><span>Passed</span><strong>{passed}</strong></div>
        <div class="stat"><span>Failed</span><strong>{failed}</strong></div>
        <div class="stat"><span>Skipped</span><strong>{skipped}</strong></div>
      </div>
      <div class="links">
        <a href="{html.escape(relative_to(AGGREGATE_OBSERVATIONS_TSV, REPORT_DIR))}">Aggregate observations TSV</a>
        <a href="{html.escape(relative_to(AGGREGATE_REPORTS_JSONL, REPORT_DIR))}">Aggregate reports JSONL</a>
        <a href="{html.escape(relative_to(AGGREGATE_REPORTS_JSON, REPORT_DIR))}">Aggregate reports JSON</a>
      </div>
    </section>
    <section class="panel">
      <div class="view-tabs">
        <button class="view-tab active" data-view="observations-view" type="button">Observations</button>
        <button class="view-tab" data-view="reports-view" type="button">Reports</button>
      </div>
      <div id="observations-view" class="view-panel active">
        <h2>Observations</h2>
        <p class="meta">One row per participant, assay, and variant observation.</p>
        <input id="observation-filter" type="search" placeholder="Filter observations by participant, assay, rsid, variant, genotype, outcome, facets">
        <div class="table-wrap">
          <table id="observations">
            <thead>
              <tr>
                {observation_header_cells}
              </tr>
            </thead>
            <tbody>
              {''.join(observation_html_rows)}
            </tbody>
          </table>
        </div>
      </div>
      <div id="reports-view" class="view-panel">
        <h2>Reports</h2>
        <p class="meta">One report per participant and assay, rendered from bioscript:report:1.0 fields.</p>
        <input id="report-filter" type="search" placeholder="Filter reports by participant, assay, status, field, value, metric">
        <div id="reports" class="report-grid">
          {report_cards or '<p class="meta">No reports emitted.</p>'}
        </div>
      </div>
    </section>
  </main>
  <script>
    const viewTabs = Array.from(document.querySelectorAll('.view-tab'));
    const viewPanels = Array.from(document.querySelectorAll('.view-panel'));
    viewTabs.forEach((tab) => {{
      tab.addEventListener('click', () => {{
        viewTabs.forEach((item) => item.classList.remove('active'));
        viewPanels.forEach((item) => item.classList.remove('active'));
        tab.classList.add('active');
        document.getElementById(tab.dataset.view).classList.add('active');
      }});
    }});

    const observationFilter = document.getElementById('observation-filter');
    const observationTable = document.getElementById('observations');
    const observationBody = observationTable.querySelector('tbody');
    const sortButtons = Array.from(observationTable.querySelectorAll('.sort-button'));
    let observationRows = Array.from(observationBody.querySelectorAll('tr'));

    function visibleObservationRows() {{
      return observationRows.filter((row) => row.style.display !== 'none');
    }}

    observationFilter.addEventListener('input', () => {{
      const query = observationFilter.value.toLowerCase();
      observationRows.forEach((row) => {{
        row.style.display = row.textContent.toLowerCase().includes(query) ? '' : 'none';
      }});
    }});

    function cellValue(row, column) {{
      return (row.children[column]?.textContent || '').trim();
    }}

    function compareValues(left, right, direction) {{
      const leftNumber = Number(left);
      const rightNumber = Number(right);
      const bothNumeric = left !== '' && right !== '' && Number.isFinite(leftNumber) && Number.isFinite(rightNumber);
      const result = bothNumeric
        ? leftNumber - rightNumber
        : left.localeCompare(right, undefined, {{ numeric: true, sensitivity: 'base' }});
      return direction === 'asc' ? result : -result;
    }}

    sortButtons.forEach((button) => {{
      button.addEventListener('click', () => {{
        const column = Number(button.dataset.column);
        const current = button.dataset.direction || 'desc';
        const direction = current === 'asc' ? 'desc' : 'asc';
        sortButtons.forEach((item) => {{
          item.classList.remove('sorted-asc', 'sorted-desc');
          delete item.dataset.direction;
        }});
        button.dataset.direction = direction;
        button.classList.add(direction === 'asc' ? 'sorted-asc' : 'sorted-desc');
        observationRows = observationRows.sort((left, right) =>
          compareValues(cellValue(left, column), cellValue(right, column), direction)
        );
        observationRows.forEach((row) => observationBody.appendChild(row));
      }});
    }});

    const reportFilter = document.getElementById('report-filter');
    const reportCards = Array.from(document.querySelectorAll('.report-card'));
    reportFilter.addEventListener('input', () => {{
      const query = reportFilter.value.toLowerCase();
      reportCards.forEach((card) => {{
        card.style.display = card.textContent.toLowerCase().includes(query) ? '' : 'none';
      }});
    }});
  </script>
</body>
</html>
"""


def render_index_markdown(records: list[dict[str, Any]]) -> str:
    rows = [
        {
            "status": record["status"],
            "assay": record["assay"],
            "participant": record["participant"],
            "result": record["result_display"],
            "evidence": record["evidence_display"],
            "info": record["info_display"],
            "duration": record["duration_display"],
            "detail": record["detail_markdown_href"],
        }
        for record in records
    ]
    passed = sum(record["status"] == "pass" for record in records)
    failed = sum(record["status"] == "fail" for record in records)
    skipped = sum(record["status"] == "skip" for record in records)
    return "\n".join(
        [
            "# Local Assay Test Report",
            "",
            f"- Total runs: `{len(records)}`",
            f"- Passed: `{passed}`",
            f"- Failed: `{failed}`",
            f"- Skipped: `{skipped}`",
            "",
            markdown_table(rows, "No runs."),
            "",
            f"Aggregate observations TSV: `{relative_to(AGGREGATE_OBSERVATIONS_TSV, REPO_ROOT)}`",
            f"Aggregate reports JSONL: `{relative_to(AGGREGATE_REPORTS_JSONL, REPO_ROOT)}`",
            f"Aggregate reports JSON: `{relative_to(AGGREGATE_REPORTS_JSON, REPO_ROOT)}`",
            "",
            f"HTML report: `{relative_to(REPORT_DIR / 'index.html', REPO_ROOT)}`",
        ]
    )


def generate_report(records: list[dict[str, Any]]) -> None:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_RUNS_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DATA_DIR.mkdir(parents=True, exist_ok=True)
    write_aggregate_report_data(records)

    index_html = REPORT_DIR / "index.html"
    index_md = REPORT_DIR / "README.md"

    for record in records:
        detail_path = REPO_ROOT / record["detail_file"]
        detail_path.parent.mkdir(parents=True, exist_ok=True)
        md_path = REPO_ROOT / record["detail_markdown_file"]
        md_path.parent.mkdir(parents=True, exist_ok=True)
        record["index_href"] = relative_to(index_html, detail_path.parent)
        detail_path.write_text(render_detail_html(record), encoding="utf-8")
        md_path.write_text(render_detail_markdown(record), encoding="utf-8")

    index_html.write_text(render_index_html(records), encoding="utf-8")
    index_md.write_text(render_index_markdown(records), encoding="utf-8")


def run_one(assay_path: Path, input_path: Path, variants: list[VariantDefinition], debug: bool) -> dict[str, Any]:
    run_started = time.monotonic()
    assay = assay_name_from_path(assay_path)
    participant = participant_from_path(input_path)
    source = detect_source(input_path)
    output_path = output_path_for(assay_path, input_path)
    trace_path = trace_path_for(assay_path, input_path)
    genome_reads_path = genome_reads_path_for(assay_path, input_path)
    timing_path = timing_path_for(assay_path, input_path)
    observations_path = observations_path_for(assay_path, input_path)
    report_json_path = report_json_path_for(assay_path, input_path)
    detail_json_path = details_json_path_for(assay_path, input_path)
    detail_file = detail_page_path_for(assay_path, input_path)
    detail_markdown_file = detail_markdown_path_for(assay_path, input_path)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    trace_path.parent.mkdir(parents=True, exist_ok=True)
    genome_reads_path.parent.mkdir(parents=True, exist_ok=True)
    observations_path.parent.mkdir(parents=True, exist_ok=True)
    report_json_path.parent.mkdir(parents=True, exist_ok=True)
    detail_json_path.parent.mkdir(parents=True, exist_ok=True)
    for stale_path in (output_path, trace_path, genome_reads_path, timing_path, observations_path, report_json_path, detail_json_path, detail_file, detail_markdown_file):
        if stale_path.exists():
            stale_path.unlink()

    run_started_at = datetime.now(timezone.utc)
    timestamp = run_started_at.replace(microsecond=0).isoformat().replace("+00:00", "Z")
    runtime_root = runtime_root_for(assay_path, input_path, output_path)
    version = assay_version_metadata(assay_path)

    # For YAML variant assays, generate a probe script; panel YAML runs directly in bioscript.
    if assay_path.is_dir() or is_yaml_variant_file(assay_path):
        script_path = ensure_probe_script(assay_path, variants)
    else:
        script_path = assay_path

    command = [
        str(BIOSCRIPT_BINARY),
        str(script_path),
        "--root",
        str(runtime_root),
        "--input-file",
        relative_to(input_path, runtime_root),
        "--output-file",
        relative_to(output_path, runtime_root),
        "--participant-id",
        participant,
        "--trace-report",
        relative_to(trace_path, runtime_root),
        "--max-duration-ms",
        "30000",
        "--max-memory-bytes",
        str(64 * 1024 * 1024),
    ]
    if detect_format(input_path) == "cram":
        command.append("--auto-index")
    if debug:
        command.extend(["--timing-report", relative_to(timing_path, runtime_root)])
    reference_args, reference_error = maybe_reference_args(input_path, runtime_root)
    if reference_error:
        return {
            "timestamp": timestamp,
            "assay": assay,
            "assay_path": relative_to(assay_path, REPO_ROOT),
            "input_file": input_path.name,
            "input_file_path": relative_to(input_path, REPO_ROOT),
            "input_source": source,
            "input_format": detect_format(input_path),
            "participant": participant,
            "status": "skip",
            "result": "",
            "result_display": "SKIPPED",
            "evidence": "",
            "evidence_display": "",
            "info": format_info_display(variants),
            "info_display": format_info_display(variants),
            "duration_ms": 0,
            "duration_display": "0ms",
            "exit_code": "",
            "error": reference_error,
            "output_file": relative_to(output_path, REPO_ROOT),
            "trace_file": "",
            "timing_file": "",
            "genome_reads_file": "",
            "observations_file": "",
            "report_json_file": "",
            "details_file": relative_to(detail_json_path, REPO_ROOT),
            "detail_file": relative_to(detail_file, REPO_ROOT),
            "detail_markdown_file": relative_to(detail_markdown_file, REPO_ROOT),
            "detail_href": relative_to(detail_file, REPORT_DIR),
            "detail_markdown_href": relative_to(detail_markdown_file, REPORT_DIR),
            "command": command,
            "stdout": "",
            "stderr": "",
            "assay_sha256": version["file_sha256"],
            "git_commit": version["git_commit"],
            "git_commit_short": version["git_commit_short"],
            "git_status": version["git_status"],
        }

    command.extend(reference_args)
    exit_code, stdout, stderr, duration_ms, timed_out = run_subprocess(command)
    output_rows = read_tsv_rows(output_path)
    result = summarize_output(output_rows)
    error = ""
    status = "pass" if exit_code == 0 else "fail"
    if exit_code != 0:
        lines = [line for line in (stderr or stdout).splitlines() if line.strip()]
        error = lines[-1] if lines else "bioscript command failed"
        if timed_out:
            error = f"timed out after {HARD_TIMEOUT_SECONDS}s"

    probe_ok = False
    probe_error = ""
    probe_duration_ms = 0
    if exit_code == 0:
        probe_started = time.monotonic()
        probe_ok, probe_error = run_probe(assay_path, input_path, participant, variants, genome_reads_path)
        probe_duration_ms = int((time.monotonic() - probe_started) * 1000)
    if debug:
        timing_path.parent.mkdir(parents=True, exist_ok=True)
        existing = timing_path.read_text(encoding="utf-8") if timing_path.exists() else "stage\tduration_ms\tdetail\n"
        if not existing.endswith("\n"):
            existing += "\n"
        existing += f"runner_assay_subprocess\t{duration_ms}\tcommand execution wall clock\n"
        if exit_code == 0:
            existing += f"runner_genome_reads\t{probe_duration_ms}\tgenome reads extraction\n"
        existing += f"runner_total\t{int((time.monotonic() - run_started) * 1000)}\tfull run_one wall clock\n"
        timing_path.write_text(existing, encoding="utf-8")

    genome_rows = read_tsv_rows(genome_reads_path) if probe_ok and genome_reads_path.exists() else []
    info_display = format_info_display(variants)
    evidence_display = format_evidence_display(genome_rows)
    result_display = format_result_display(status, output_rows, genome_rows, variants)

    record = {
        "timestamp": timestamp,
        "assay": assay,
        "assay_path": relative_to(assay_path, REPO_ROOT),
        "input_file": input_path.name,
        "input_file_path": relative_to(input_path, REPO_ROOT),
        "input_source": source,
        "input_format": detect_format(input_path),
        "participant": participant,
        "status": status,
        "result": result,
        "result_display": result_display,
        "evidence": evidence_display,
        "evidence_display": evidence_display,
        "info": info_display,
        "info_display": info_display,
        "duration_ms": duration_ms,
        "duration_display": format_duration(duration_ms),
        "exit_code": exit_code,
        "error": error,
        "output_file": relative_to(output_path, REPO_ROOT),
        "trace_file": relative_to(trace_path, REPO_ROOT) if trace_path.exists() else "",
        "timing_file": relative_to(timing_path, REPO_ROOT) if debug and timing_path.exists() else "",
        "genome_reads_file": relative_to(genome_reads_path, REPO_ROOT) if probe_ok and genome_reads_path.exists() else "",
        "genome_reads_error": probe_error,
        "observations_file": relative_to(observations_path, REPO_ROOT),
        "report_json_file": relative_to(report_json_path, REPO_ROOT),
        "details_file": relative_to(detail_json_path, REPO_ROOT),
        "detail_file": relative_to(detail_file, REPO_ROOT),
        "detail_markdown_file": relative_to(detail_markdown_file, REPO_ROOT),
        "detail_href": relative_to(detail_file, REPORT_DIR),
        "detail_markdown_href": relative_to(detail_markdown_file, REPORT_DIR),
        "command": command,
        "stdout": stdout,
        "stderr": stderr,
        "assay_sha256": version["file_sha256"],
        "git_commit": version["git_commit"],
        "git_commit_short": version["git_commit_short"],
        "git_status": version["git_status"],
    }
    observations = build_observations(record, output_rows, genome_rows, variants)
    write_observations(observations_path, observations)
    finished_at = datetime.now(timezone.utc)
    report_object = build_report_object(
        record,
        observations,
        run_started_at.isoformat().replace("+00:00", "Z"),
        finished_at.isoformat().replace("+00:00", "Z"),
        probe_duration_ms,
    )
    report_json_path.write_text(json.dumps(report_object, indent=2), encoding="utf-8")
    detail_json_path.write_text(json.dumps(record, indent=2), encoding="utf-8")
    return record


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run bioscript assays against local test data and generate static reports.")
    parser.add_argument("assay_path", nargs="?", help="Single assay path to run (shorthand for --assay)")
    parser.add_argument("--assay", help="Single assay path to run")
    parser.add_argument("--jobs", default="auto", help="Parallel jobs to run across independent assay x input tasks, or 'auto' for CPU count")
    parser.add_argument("-k", "--filter", help="Pytest-style keyword filter for discovered assays, e.g. 'apol1 and not glp1'")
    parser.add_argument("--ignore-assay", help="Comma-separated glob patterns to exclude discovered assay paths")
    parser.add_argument("--input", help="Input file or directory to run against")
    parser.add_argument("--only", help="Comma-separated glob patterns to include matching input paths")
    parser.add_argument("--exclude", help="Comma-separated glob patterns to exclude matching input paths")
    parser.add_argument("--debug", action="store_true", help="Emit per-stage timing reports and include them in the output report")
    parser.add_argument("--list", action="store_true", help="List assays and test data")
    parser.add_argument("--report", action="store_true", help="Open the generated HTML report in the default browser after running")
    args = parser.parse_args(argv)
    if args.assay and args.assay_path and args.assay != args.assay_path:
        parser.error("use either positional assay path or --assay, not both")
    if not args.assay and args.assay_path:
        args.assay = args.assay_path
    return args


def resolve_jobs(value: str) -> int:
    if value.strip().lower() == "auto":
        return max(1, os.cpu_count() or 1)
    try:
        jobs = int(value)
    except ValueError as exc:
        raise ValueError(f"invalid --jobs value {value!r}; expected integer or 'auto'") from exc
    if jobs < 1:
        raise ValueError(f"invalid --jobs value {value!r}; expected integer >= 1")
    return jobs


def resolve_assays(filter_assay: str | None, keyword_filter: str | None, ignore_assay: str | None) -> list[Path]:
    if filter_assay:
        path = Path(filter_assay)
        if not path.is_absolute():
            cwd_candidate = (Path.cwd() / path).resolve()
            repo_candidate = (REPO_ROOT / path).resolve()
            path = cwd_candidate if cwd_candidate.exists() else repo_candidate
        return [path.resolve()]
    assays = find_assays()
    if keyword_filter:
        assays = [path for path in assays if keyword_expr_matches(assay_keywords(path), keyword_filter)]
    if ignore_assay:
        assays = [path for path in assays if not matches_any_pattern(path, ignore_assay)]
    return assays


def matches_any_pattern(path: Path, patterns: str | None) -> bool:
    if not patterns:
        return False
    rel = relative_to(path, REPO_ROOT)
    for pattern in [item.strip() for item in patterns.split(",") if item.strip()]:
        if fnmatch(rel, pattern) or fnmatch(path.as_posix(), pattern) or fnmatch(path.name, pattern):
            return True
    return False


def assay_keywords(assay_path: Path) -> str:
    rel = relative_to(assay_path, REPO_ROOT)
    return " ".join(
        part.lower()
        for part in [
            assay_path.stem,
            assay_path.name,
            rel,
            str(assay_path.parent.relative_to(REPO_ROOT)) if assay_path.is_relative_to(REPO_ROOT) else assay_path.parent.as_posix(),
        ]
        if part
    )


def tokenize_keyword_expr(expr: str) -> list[str]:
    tokens: list[str] = []
    current: list[str] = []
    for ch in expr:
        if ch in "()":
            if current:
                tokens.append("".join(current))
                current = []
            tokens.append(ch)
            continue
        if ch.isspace():
            if current:
                tokens.append("".join(current))
                current = []
            continue
        current.append(ch)
    if current:
        tokens.append("".join(current))
    return tokens


def keyword_expr_matches(text: str, expr: str) -> bool:
    tokens = tokenize_keyword_expr(expr.lower())
    if not tokens:
        return True
    position = 0

    def parse_expr() -> bool:
        nonlocal position
        value = parse_term()
        while position < len(tokens) and tokens[position] == "or":
            position += 1
            value = value or parse_term()
        return value

    def parse_term() -> bool:
        nonlocal position
        value = parse_factor()
        while position < len(tokens) and tokens[position] == "and":
            position += 1
            value = value and parse_factor()
        return value

    def parse_factor() -> bool:
        nonlocal position
        if position >= len(tokens):
            return True
        token = tokens[position]
        if token == "not":
            position += 1
            return not parse_factor()
        if token == "(":
            position += 1
            value = parse_expr()
            if position < len(tokens) and tokens[position] == ")":
                position += 1
            return value
        if token == ")":
            return True
        position += 1
        return token in text

    value = parse_expr()
    return value and position >= len(tokens)


def resolve_inputs(filter_input: str | None, only: str | None, exclude: str | None) -> list[Path]:
    target = DATA_DIR
    if filter_input:
        path = Path(filter_input)
        target = path if path.is_absolute() else (REPO_ROOT / path)
    target = Path(os.path.abspath(target))
    inputs = [target] if target.is_file() else find_inputs(target)
    if only:
        inputs = [path for path in inputs if matches_any_pattern(path, only)]
    if exclude:
        inputs = [path for path in inputs if not matches_any_pattern(path, exclude)]
    return inputs


def print_run_status(record: dict[str, Any]) -> None:
    label = f"{record['assay']} <- {record['input_source']}/{record['input_file']}"
    if record["status"] == "pass":
        suffix = f" -> {record['result']}" if record["result"] else ""
        print(f"[PASS] {label} ({record['duration_display']}){suffix}")
        if record.get("genome_reads_error"):
            print(f"       genome-reads probe: {record['genome_reads_error']}")
    elif record["status"] == "skip":
        print(f"[SKIP] {label} ({record['error']})")
    else:
        print(f"[FAIL] {label} (exit {record['exit_code']}, {record['duration_display']})")
        lines = [line for line in record.get("stderr", "").splitlines() if line.strip()]
        for line in lines[-5:]:
            print(f"       {line}")


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    try:
        jobs = resolve_jobs(args.jobs)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    if not BS.exists():
        print(f"bioscript/bs not found at {BS}", file=sys.stderr)
        return 1
    try:
        ensure_bioscript_binary()
    except RuntimeError as exc:
        print(str(exc), file=sys.stderr)
        return 1

    assays = resolve_assays(args.assay, args.filter, args.ignore_assay)
    inputs = resolve_inputs(args.input, args.only, args.exclude)
    if not assays:
        print(f"No assays found in {ASSAY_DIR}", file=sys.stderr)
        return 1
    if not inputs:
        target = args.input or str(DATA_DIR)
        print(f"No test data files found in {target}", file=sys.stderr)
        return 1
    if args.list:
        return list_mode(assays, inputs)

    print("\nBioScript Assay Test Runner")
    print("===========================")
    print(f"Bioscript: {BS}")
    print(f"Assays:    {ASSAY_DIR}")
    print(f"Data:      {DATA_DIR}")
    print(f"Output:    {OUTPUT_DIR}")
    print(f"Report:    {REPORT_DIR}")
    print(f"Jobs:      {jobs}")
    print("")
    print(f"[INFO] Running {len(assays)} assay(s) x {len(inputs)} input(s)")
    print("")

    if REPORT_DIR.exists():
        shutil.rmtree(REPORT_DIR)
    if RESULTS_CSV.exists():
        RESULTS_CSV.unlink()

    variant_cache: dict[Path, list[VariantDefinition]] = {}
    for assay in assays:
        if assay.is_dir():
            variant_cache[assay] = extract_yaml_variant_definitions(assay)
        elif is_yaml_variant_file(assay):
            vd = yaml_to_variant_definition(assay)
            variant_cache[assay] = [vd] if vd else []
        elif is_yaml_panel_file(assay):
            variant_cache[assay] = yaml_to_panel_variant_definitions(assay)
        else:
            variant_cache[assay] = extract_variant_definitions(assay)
    task_specs = [
        (assay_path, input_path, variant_cache[assay_path], args.debug)
        for assay_path in assays
        for input_path in inputs
    ]
    records: list[dict[str, Any]] = []
    passed = failed = skipped = 0

    if jobs == 1:
        current_assay: Path | None = None
        for assay_path, input_path, variants, debug in task_specs:
            if current_assay != assay_path:
                if current_assay is not None:
                    print("")
                assay_name = assay_name_from_path(assay_path)
                print(f"--- {assay_name} ({relative_to(assay_path, REPO_ROOT)}) ---")
                print("")
                current_assay = assay_path
            record = run_one(assay_path, input_path, variants, debug)
            records.append(record)
            print_run_status(record)
            if record["status"] == "pass":
                passed += 1
            elif record["status"] == "fail":
                failed += 1
            else:
                skipped += 1
        if current_assay is not None:
            print("")
    else:
        print(f"[INFO] Parallel mode enabled with {jobs} workers")
        print("")
        future_to_index = {}
        with ThreadPoolExecutor(max_workers=jobs) as executor:
            for index, (assay_path, input_path, variants, debug) in enumerate(task_specs):
                future = executor.submit(run_one, assay_path, input_path, variants, debug)
                future_to_index[future] = index
            ordered_records: list[dict[str, Any] | None] = [None] * len(task_specs)
            for future in as_completed(future_to_index):
                index = future_to_index[future]
                record = future.result()
                ordered_records[index] = record
                print_run_status(record)
        records = [record for record in ordered_records if record is not None]
        passed = sum(1 for record in records if record["status"] == "pass")
        failed = sum(1 for record in records if record["status"] == "fail")
        skipped = sum(1 for record in records if record["status"] == "skip")

    write_results_csv(records)
    generate_report(records)

    print("===========================")
    print(f"Results: {passed} passed, {failed} failed, {skipped} skipped")
    print(f"CSV:     {relative_to(RESULTS_CSV, REPO_ROOT)}")
    print(f"OBS:     {relative_to(AGGREGATE_OBSERVATIONS_TSV, REPO_ROOT)}")
    print(f"REPORTS: {relative_to(AGGREGATE_REPORTS_JSONL, REPO_ROOT)}")
    print(f"HTML:    {relative_to(REPORT_DIR / 'index.html', REPO_ROOT)}")
    print(f"MD:      {relative_to(REPORT_DIR / 'README.md', REPO_ROOT)}")

    if args.report:
        report_url = (REPORT_DIR / "index.html").resolve().as_uri()
        print(f"Opening: {report_url}")
        import webbrowser
        webbrowser.open(report_url)

    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
