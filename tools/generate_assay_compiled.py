#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "bioscript"))

from assay_compiled import load_assay_package_compiled_yaml, write_assay_package_compiled


def iter_assay_roots(input_paths: list[Path]) -> list[Path]:
    assay_roots: set[Path] = set()

    for raw_path in input_paths:
        path = raw_path.resolve()
        if path.is_file():
            if path.name != "assay.yaml":
                raise RuntimeError(f"{path} is a file, but not an assay.yaml file")
            assay_roots.add(path.parent)
            continue

        if not path.is_dir():
            raise RuntimeError(f"{path} does not exist")

        if (path / "assay.yaml").exists():
            assay_roots.add(path)
            continue

        for manifest_path in path.rglob("assay.yaml"):
            assay_roots.add(manifest_path.parent)

    if not assay_roots:
        raise RuntimeError("No assay packages found")

    return sorted(assay_roots)


def check_assay_compiled(assay_root: Path) -> tuple[bool, str]:
    target_path = assay_root / "assay.compiled.yaml"
    expected = load_assay_package_compiled_yaml(assay_root)
    if not target_path.exists():
        return False, f"missing {target_path}"
    current = target_path.read_text(encoding="utf-8")
    if current != expected:
        return False, f"out of date {target_path}"
    return True, f"ok {target_path}"


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate or check assay.compiled.yaml files for bioscript assay packages.")
    parser.add_argument(
        "paths",
        nargs="+",
        type=Path,
        help="Assay package directories, assay.yaml files, or directories containing assay packages.",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Check whether assay.compiled.yaml is present and up to date instead of writing files.",
    )
    args = parser.parse_args()

    assay_roots = iter_assay_roots(args.paths)

    if args.check:
        all_ok = True
        for assay_root in assay_roots:
            ok, message = check_assay_compiled(assay_root)
            print(message)
            all_ok = all_ok and ok
        return 0 if all_ok else 1

    for assay_root in assay_roots:
        written_path = write_assay_package_compiled(assay_root)
        print(written_path)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
