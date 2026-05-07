#!/usr/bin/env python3
"""Format YAML files with the repository's canonical PyYAML style."""

from __future__ import annotations

import argparse
import difflib
import sys
from pathlib import Path
from typing import Any

import yaml


EXCLUDED_DIRS = {
    ".git",
    ".mypy_cache",
    ".pytest_cache",
    "__pycache__",
    "target",
    "test-output",
}

DEFAULT_PATHS = (Path("projects"),)


class CleanDumper(yaml.SafeDumper):
    def increase_indent(self, flow: bool = False, indentless: bool = False):
        return super().increase_indent(flow, False)


def needs_quotes(data: str) -> bool:
    if data == "":
        return True
    try:
        parsed = yaml.safe_load(data)
    except yaml.YAMLError:
        return False
    return not isinstance(parsed, str) or parsed != data


def str_representer(dumper: yaml.SafeDumper, data: str) -> yaml.nodes.ScalarNode:
    if "\n" in data:
        return dumper.represent_scalar("tag:yaml.org,2002:str", data, style="|")
    if needs_quotes(data):
        return dumper.represent_scalar("tag:yaml.org,2002:str", data, style='"')
    return dumper.represent_scalar("tag:yaml.org,2002:str", data)


CleanDumper.add_representer(str, str_representer)


def yaml_paths(paths: list[Path]) -> list[Path]:
    files: list[Path] = []
    for path in paths:
        if path.is_dir():
            for child in path.rglob("*"):
                if any(part in EXCLUDED_DIRS for part in child.parts):
                    continue
                if child.suffix.lower() in {".yaml", ".yml"}:
                    files.append(child)
        elif path.suffix.lower() in {".yaml", ".yml"}:
            files.append(path)
    return sorted(set(files))


def formatted_yaml(value: Any) -> str:
    text = yaml.dump(
        value,
        Dumper=CleanDumper,
        default_flow_style=False,
        sort_keys=False,
        allow_unicode=True,
        width=120,
    )
    return text if text.endswith("\n") else f"{text}\n"


def format_file(path: Path, check: bool, diff: bool) -> bool:
    original = path.read_text(encoding="utf-8")
    data = yaml.safe_load(original)
    formatted = formatted_yaml(data)
    if original == formatted:
        return False

    if check:
        print(f"would format {path}")
        if diff:
            sys.stdout.writelines(
                difflib.unified_diff(
                    original.splitlines(keepends=True),
                    formatted.splitlines(keepends=True),
                    fromfile=str(path),
                    tofile=str(path),
                )
            )
        return True

    path.write_text(formatted, encoding="utf-8")
    print(f"formatted {path}")
    return True


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "paths",
        nargs="*",
        type=Path,
        help="YAML files or directories to format (defaults to projects/)",
    )
    parser.add_argument("--check", action="store_true", help="fail if any YAML file would be reformatted")
    parser.add_argument("--diff", action="store_true", help="print unified diffs in --check mode")
    args = parser.parse_args(argv)

    changed = False
    paths = args.paths or list(DEFAULT_PATHS)

    for path in yaml_paths(paths):
        try:
            changed = format_file(path, args.check, args.diff) or changed
        except (OSError, yaml.YAMLError) as exc:
            print(f"{path}: {exc}", file=sys.stderr)
            return 2

    if args.check and changed:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
