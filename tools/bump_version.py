#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
VERSION_RE = re.compile(r"^\d+\.\d+\.\d+(?:-[0-9A-Za-z.-]+)?$")
IMAGE_RE = re.compile(r"ghcr\.io/(?:openmined|madhavajay)/exvitae:[0-9A-Za-z][0-9A-Za-z.-]*")
IMAGE_REPOSITORY = "ghcr.io/madhavajay/exvitae"


def current_version() -> str:
    return (ROOT / "VERSION").read_text(encoding="utf-8").strip()


def bump_patch(version: str) -> str:
    major, minor, patch = version.split("-", 1)[0].split(".")
    return f"{major}.{minor}.{int(patch) + 1}"


def replace_yaml_metadata_version(text: str, version: str) -> str:
    return re.sub(r"(?m)^  version: .*$", f"  version: {version}", text, count=1)


def update_file(path: Path, version: str) -> None:
    text = path.read_text(encoding="utf-8")
    text = replace_yaml_metadata_version(text, version)
    text = IMAGE_RE.sub(f"{IMAGE_REPOSITORY}:{version}", text)
    path.write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Bump ExVitae image and flow version metadata.")
    parser.add_argument("version", nargs="?", help="New semver version. Defaults to patch bump.")
    args = parser.parse_args()

    version = args.version or bump_patch(current_version())
    if not VERSION_RE.fullmatch(version):
        parser.error(f"invalid semver version: {version}")

    (ROOT / "VERSION").write_text(f"{version}\n", encoding="utf-8")
    for path in sorted((ROOT / "flows").glob("exvitae-*/*.yaml")):
        update_file(path, version)
    for path in sorted((ROOT / "flows").glob("exvitae-*/main.nf")):
        update_file(path, version)
    print(version)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
