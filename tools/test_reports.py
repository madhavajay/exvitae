#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


REPO_ROOT = Path(__file__).resolve().parent.parent
# Default: the bioscript shell wrapper nested next to this tools/ dir. Can be
# overridden via `--bioscript /path/to/bs` or `BIOSCRIPT_BIN=/path/to/bs`,
# which lets the same sample suite run against any bioscript checkout (e.g.
# a working branch in another workspace) for parity testing.
DEFAULT_BS = REPO_ROOT / "bioscript" / "bs"
DEFAULT_SAMPLES = REPO_ROOT / "samples.yaml"
DEFAULT_PRIVATE_SAMPLES = REPO_ROOT / "samples.private.yaml"


@dataclass
class TestCase:
    source: Path
    assay: Path
    sample: dict[str, Any]
    defaults: dict[str, Any]

    @property
    def id(self) -> str:
        return str(self.sample["id"])

    @property
    def output_dir(self) -> Path:
        base = resolve_path(str(self.sample.get("output_dir") or self.defaults.get("output_dir") or "test-output/report-tests"))
        return base / safe_name(assay_run_id(self.assay)) / safe_name(self.id)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run YAML-defined bioscript report tests.")
    parser.add_argument("--samples", default=str(DEFAULT_SAMPLES), help="Primary samples YAML file.")
    parser.add_argument("--private-samples", default=str(DEFAULT_PRIVATE_SAMPLES), help="Optional private samples YAML file.")
    parser.add_argument("--sample", action="append", help="Only run sample id matching this value. Can be repeated.")
    parser.add_argument("--assay", action="append", help="Only run assay path/name matching this value. Can be repeated.")
    parser.add_argument("--list", action="store_true", help="List resolved test cases without running them.")
    parser.add_argument("--no-private", action="store_true", help="Do not load samples.private.yaml.")
    parser.add_argument("--keep-going", action="store_true", help="Continue after failures.")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without running them.")
    parser.add_argument("--generate-only", action="store_true", help="Run reports but skip YAML assertions.")
    parser.add_argument(
        "--bioscript",
        default=os.environ.get("BIOSCRIPT_BIN", str(DEFAULT_BS)),
        help="Path to the bioscript binary/wrapper to run (env: BIOSCRIPT_BIN).",
    )
    args = parser.parse_args()
    # Stash on env so run_case can read it without threading the path through
    # every helper. This keeps the diff to test_reports.py minimal.
    os.environ["BIOSCRIPT_BIN_RESOLVED"] = args.bioscript

    config_files = [Path(args.samples)]
    private_path = Path(args.private_samples)
    if not args.no_private and private_path.exists():
        config_files.append(private_path)

    cases = load_cases(config_files)
    cases = filter_cases(cases, sample_filters=args.sample or [], assay_filters=args.assay or [])

    if args.list:
        for case in cases:
            print(f"{case.source.name}: {case.id} x {display_path(case.assay)} -> {display_path(case.output_dir)}")
        return 0

    if not cases:
        print("no report test cases selected", file=sys.stderr)
        return 2

    failures: list[str] = []
    for idx, case in enumerate(cases, start=1):
        label = f"[{idx}/{len(cases)}] {case.id} x {display_path(case.assay)}"
        print(label)
        try:
            run_case(case, dry_run=args.dry_run, assert_results=not args.generate_only)
        except Exception as err:
            failures.append(f"{label}: {err}")
            print(f"FAIL {label}: {err}", file=sys.stderr)
            if not args.keep_going:
                break
        else:
            print(f"PASS {label}")

    if failures:
        print("\nfailures:", file=sys.stderr)
        for failure in failures:
            print(f"  - {failure}", file=sys.stderr)
        return 1
    return 0


def load_cases(paths: list[Path]) -> list[TestCase]:
    cases: list[TestCase] = []
    for path in paths:
        config = load_yaml(path)
        defaults = dict(config.get("defaults") or {})
        assays = resolve_assays(config.get("assays"))
        samples = config.get("samples") or []
        if not isinstance(samples, list):
            raise ValueError(f"{path}: samples must be a list")
        for sample in samples:
            if not isinstance(sample, dict):
                raise ValueError(f"{path}: each sample must be a map")
            if not sample.get("id"):
                raise ValueError(f"{path}: sample is missing id")
            if not sample.get("input_file"):
                raise ValueError(f"{path}: sample {sample.get('id')} is missing input_file")
            for assay in assays:
                cases.append(TestCase(source=path, assay=assay, sample=sample, defaults=defaults))
    return cases


def load_yaml(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(path)
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ValueError(f"{path}: top-level YAML must be a map")
    return data


def resolve_assays(value: Any) -> list[Path]:
    if value is None:
        candidates = sorted((REPO_ROOT / "assays" / "pgx").glob("pgx-1*/manifest.yaml"))
        candidates.extend(sorted((REPO_ROOT / "assays" / "pgx").glob("pgx-1*/pgx-1*.zip")))
        return candidates
    if not isinstance(value, list):
        raise ValueError("assays must be a list")
    assays: list[Path] = []
    for item in value:
        text = str(item)
        matches = sorted(REPO_ROOT.glob(text)) if any(ch in text for ch in "*?[") else []
        if matches:
            assays.extend(matches)
        else:
            assays.append(resolve_path(text))
    return assays


def filter_cases(cases: list[TestCase], sample_filters: list[str], assay_filters: list[str]) -> list[TestCase]:
    def matches(value: str, filters: list[str]) -> bool:
        return not filters or any(f == value or f in value for f in filters)

    return [
        case
        for case in cases
        if matches(case.id, sample_filters)
        and (not assay_filters or matches(str(case.assay), assay_filters) or matches(case.assay.name, assay_filters))
    ]


def run_case(case: TestCase, dry_run: bool = False, assert_results: bool = True) -> None:
    input_file = require_file(case, "input_file")
    output_dir = case.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        os.environ.get("BIOSCRIPT_BIN_RESOLVED", str(DEFAULT_BS)),
        "report",
        str(case.assay),
        "--root",
        str(report_root(case)),
        "--input-file",
        str(input_file),
        "--output-dir",
        str(output_dir),
    ]
    add_optional_path_arg(cmd, case, "input_index", "--input-index")
    add_optional_path_arg(cmd, case, "reference_file", "--reference-file")
    add_optional_path_arg(cmd, case, "reference_index", "--reference-index")
    add_optional_value_arg(cmd, case, "input_format", "--input-format")
    add_optional_value_arg(cmd, case, "sample_sex", "--sample-sex")

    if bool_value(case, "detect_sex", default=True):
        cmd.append("--detect-sex")
    if bool_value(case, "html", default=False):
        cmd.append("--html")
    if bool_value(case, "allow_md5_mismatch", default=False):
        cmd.append("--allow-md5-mismatch")

    max_ms = case.sample.get("analysis_max_duration_ms", case.defaults.get("analysis_max_duration_ms"))
    if max_ms:
        cmd.extend(["--analysis-max-duration-ms", str(max_ms)])

    for arg in case.sample.get("extra_args") or []:
        cmd.append(str(arg))

    print("  " + " ".join(shell_quote(part) for part in cmd))
    if dry_run:
        return

    result = subprocess.run(cmd, cwd=REPO_ROOT, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        message = result.stderr.strip() or result.stdout.strip() or f"exit {result.returncode}"
        raise RuntimeError(message)

    if assert_results:
        assert_outputs(case, output_dir)


def assert_outputs(case: TestCase, output_dir: Path) -> None:
    assertions = merged_assertions(case)
    observations_path = output_dir / "observations.tsv"
    reports_path = output_dir / "reports.jsonl"
    if not observations_path.is_file():
        raise AssertionError(f"missing {display_path(observations_path)}")
    if not reports_path.is_file():
        raise AssertionError(f"missing {display_path(reports_path)}")

    observations = read_tsv(observations_path)
    reports = read_jsonl(reports_path)

    observation_assertions = assertions.get("observations") or {}
    report_assertions = assertions.get("reports") or {}

    min_observation_rows = int(observation_assertions.get("min_rows", 0))
    if len(observations) < min_observation_rows:
        raise AssertionError(f"expected at least {min_observation_rows} observations, got {len(observations)}")
    for expected in observation_assertions.get("rows") or []:
        if not any(row_matches(row, expected) for row in observations):
            raise AssertionError(f"no observation row matched {expected}")

    min_report_rows = int(report_assertions.get("min_rows", 0))
    if len(reports) < min_report_rows:
        raise AssertionError(f"expected at least {min_report_rows} report rows, got {len(reports)}")
    required_status = report_assertions.get("require_status")
    if required_status:
        allowed_statuses = {str(value) for value in required_status} if isinstance(required_status, list) else {str(required_status)}
        for report in reports:
            if report.get("report_status") not in allowed_statuses:
                expected = ", ".join(sorted(allowed_statuses))
                raise AssertionError(f"expected report_status in [{expected}], got {report.get('report_status')}")
    for expected in report_assertions.get("contains") or []:
        if not any(json_contains(report, expected) for report in reports):
            raise AssertionError(f"no report JSON matched {expected}")


def merged_assertions(case: TestCase) -> dict[str, Any]:
    merged = dict(case.defaults.get("assertions") or {})
    sample_assertions = case.sample.get("assertions") or {}
    for key, value in sample_assertions.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            child = dict(merged[key])
            child.update(value)
            merged[key] = child
        else:
            merged[key] = value
    return merged


def row_matches(row: dict[str, str], expected: dict[str, Any]) -> bool:
    for key, expected_value in expected.items():
        actual = row.get(str(key))
        if isinstance(expected_value, list):
            if actual not in {str(v) for v in expected_value}:
                return False
        elif actual != str(expected_value):
            return False
    return True


def json_contains(report: dict[str, Any], expected: dict[str, Any]) -> bool:
    values = values_at_path(report, str(expected["path"]))
    if "equals" in expected:
        return any(value == expected["equals"] for value in values)
    if "contains" in expected:
        needle = str(expected["contains"])
        return any(needle in str(value) for value in values)
    if "exists" in expected:
        return bool(values) is bool(expected["exists"])
    raise ValueError(f"unsupported report assertion: {expected}")


def values_at_path(value: Any, path: str) -> list[Any]:
    current = [value]
    for part in path.split("."):
        next_values: list[Any] = []
        is_list = part.endswith("[]")
        key = part[:-2] if is_list else part
        for item in current:
            if isinstance(item, dict) and key in item:
                child = item[key]
                if is_list and isinstance(child, list):
                    next_values.extend(child)
                else:
                    next_values.append(child)
        current = next_values
    return current


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def require_file(case: TestCase, key: str) -> Path:
    path = resolve_path(str(case.sample[key]))
    if not path.is_file():
        raise FileNotFoundError(f"{case.id}: {key} does not exist: {display_path(path)}")
    return path


def add_optional_path_arg(cmd: list[str], case: TestCase, key: str, flag: str) -> None:
    value = case.sample.get(key)
    if value:
        path = resolve_path(str(value))
        if not path.is_file():
            raise FileNotFoundError(f"{case.id}: {key} does not exist: {display_path(path)}")
        cmd.extend([flag, str(path)])


def add_optional_value_arg(cmd: list[str], case: TestCase, key: str, flag: str) -> None:
    value = case.sample.get(key)
    if value:
        cmd.extend([flag, str(value)])


def bool_value(case: TestCase, key: str, default: bool) -> bool:
    value = case.sample.get(key, case.defaults.get(key, default))
    return bool(value)


def report_root(case: TestCase) -> Path:
    configured = case.sample.get("root") or case.defaults.get("root") or os.environ.get("BIOSCRIPT_TEST_ROOT")
    if configured:
        return resolve_path(str(configured))
    home = Path.home()
    try:
        REPO_ROOT.relative_to(home)
    except ValueError:
        return REPO_ROOT.parent
    return home


def resolve_path(value: str) -> Path:
    expanded = os.path.expandvars(os.path.expanduser(value))
    path = Path(expanded)
    if path.is_absolute():
        return path
    return REPO_ROOT / path


def display_path(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "-", value).strip("-")


def assay_run_id(path: Path) -> str:
    if path.name in {"manifest.yaml", "package.yaml"}:
        return f"{path.parent.name}-{path.stem}"
    return path.stem if path.is_file() else path.name


def shell_quote(value: str) -> str:
    if re.fullmatch(r"[A-Za-z0-9_./:=,+@%-]+", value):
        return value
    return "'" + value.replace("'", "'\"'\"'") + "'"


if __name__ == "__main__":
    raise SystemExit(main())
