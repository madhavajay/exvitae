#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import zipfile

import yaml


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_yaml(path: Path) -> dict:
    data = yaml.safe_load(path.read_text())
    if not isinstance(data, dict):
        raise RuntimeError(f"{path} did not contain a mapping")
    return data


def zip_dir(source_dir: Path, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(output_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(source_dir.rglob("*")):
            if path.is_dir():
                continue
            archive.write(path, path.relative_to(source_dir))


def build_index_entry(repo_root: Path, assay_dir: Path, artifact_base_url: str, artifact_dir: Path) -> dict:
    compiled_path = assay_dir / "assay.compiled.yaml"
    assay_path = assay_dir / "assay.yaml"
    compiled = load_yaml(compiled_path)

    assay = compiled.get("assay") or {}
    compatibility = compiled.get("compatibility") or {}
    ui = compiled.get("ui") or {}

    assay_id = str(assay["id"])
    package_version = str(assay.get("package_version") or compiled.get("version") or "1.0")
    artifact_name = f"{assay_id}-{package_version}.zip"
    artifact_path = artifact_dir / artifact_name
    zip_dir(assay_dir, artifact_path)

    return {
        "id": assay_id,
        "title": assay.get("label") or assay_id,
        "summary": assay.get("summary") or "",
        "category": assay.get("category") or "traits",
        "tags": assay.get("tags") or [],
        "disclaimer": assay.get("disclaimer"),
        "package_version": package_version,
        "source_of_truth": assay.get("source_of_truth") or "package",
        "path": assay_dir.relative_to(repo_root).as_posix(),
        "assay_path": assay_path.relative_to(repo_root).as_posix(),
        "compiled_path": compiled_path.relative_to(repo_root).as_posix(),
        "artifact_url": f"{artifact_base_url.rstrip('/')}/{artifact_name}",
        "artifact_sha256": sha256_file(artifact_path),
        "artifact_size": artifact_path.stat().st_size,
        "artifact_format": "zip",
        "works_with": compatibility.get("works_with") or [],
        "assemblies": compatibility.get("assemblies") or [],
        "notes": compatibility.get("notes") or [],
        "template": ui.get("template") or "variant-panel",
        "runnable_variant_count": len(compiled.get("runnable_variants") or []),
        "unsupported_variant_count": len(compiled.get("unsupported_variants") or []),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build exvitae assay zip artifacts and refresh assays/index.json.")
    parser.add_argument("--repo-root", default=".", help="Repository root")
    parser.add_argument("--assays-root", default="assays", help="Assays root directory relative to repo root")
    parser.add_argument("--artifact-dir", default="dist/assays", help="Artifact output directory relative to repo root")
    parser.add_argument("--index-path", default="assays/index.json", help="Index output path relative to repo root")
    parser.add_argument("--artifact-base-url", required=True, help="Public base URL for published artifacts")
    args = parser.parse_args()

    repo_root = Path(args.repo_root).resolve()
    assays_root = (repo_root / args.assays_root).resolve()
    artifact_dir = (repo_root / args.artifact_dir).resolve()
    index_path = (repo_root / args.index_path).resolve()

    if artifact_dir.exists():
        shutil.rmtree(artifact_dir)
    artifact_dir.mkdir(parents=True, exist_ok=True)

    entries = []
    for compiled_path in sorted(assays_root.rglob("assay.compiled.yaml")):
        assay_dir = compiled_path.parent
        if not (assay_dir / "assay.yaml").exists():
            continue
        entries.append(build_index_entry(repo_root, assay_dir, args.artifact_base_url, artifact_dir))

    payload = {
        "schema": "exvitae:assay-index",
        "version": "1.0",
        "assays": entries,
    }
    index_path.parent.mkdir(parents=True, exist_ok=True)
    index_path.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
