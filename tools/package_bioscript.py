#!/usr/bin/env python3
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import os
import subprocess
import sys
import zipfile
from pathlib import Path
from typing import Any

import yaml


RUNNABLE_SCHEMAS = {
    "bioscript:panel:1.0",
    "bioscript:assay:1.0",
    "bioscript:variant:1.0",
    "bioscript:variant",
    "bioscript:variant-catalogue:1.0",
}

EXBUILD_FILENAME = ".exbuild"

EXCLUDE_DIRS = {
    ".bioscript-cache",
    ".git",
    ".sessions",
    "__pycache__",
    "logs",
    "tmp",
}


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError(f"{path} must contain a YAML mapping")
    return data


def resolve_package_path(package_root: Path, base_file: Path, value: str) -> Path:
    candidate = (base_file.parent / value).resolve()
    root = package_root.resolve()
    if not candidate.is_relative_to(root):
        raise ValueError(f"{base_file}: path escapes package root: {value}")
    return candidate


def collect_includes(package_root: Path, yaml_file: Path, node: Any, out: set[Path]) -> None:
    if isinstance(node, dict):
        include = node.get("include")
        if isinstance(include, str):
            include_path = resolve_package_path(package_root, yaml_file, include)
            add_yaml_closure(package_root, include_path, out)
        for value in node.values():
            collect_includes(package_root, yaml_file, value, out)
    elif isinstance(node, list):
        for value in node:
            collect_includes(package_root, yaml_file, value, out)


def add_file(package_root: Path, path: Path, out: set[Path]) -> Path:
    path = path.resolve()
    if not path.is_relative_to(package_root.resolve()):
        raise ValueError(f"path escapes package root: {path}")
    if not path.exists():
        raise FileNotFoundError(path)
    if path.is_dir():
        raise ValueError(f"package dependency is a directory, expected file: {path}")
    out.add(path)
    return path


def add_yaml_closure(package_root: Path, yaml_file: Path, out: set[Path]) -> None:
    yaml_file = add_file(package_root, yaml_file, out)
    data = load_yaml(yaml_file)
    collect_includes(package_root, yaml_file, data, out)

    schema = str(data.get("schema", ""))
    if schema == "bioscript:variant-catalogue:1.0":
        for table_key in ("variants", "findings"):
            table = data.get(table_key)
            if isinstance(table, dict) and isinstance(table.get("source"), str):
                add_file(package_root, resolve_package_path(package_root, yaml_file, table["source"]), out)
        return

    if schema not in {"bioscript:panel:1.0", "bioscript:assay:1.0"}:
        return

    for member in data.get("members", []) or []:
        if not isinstance(member, dict) or not isinstance(member.get("path"), str):
            continue
        member_path = resolve_package_path(package_root, yaml_file, member["path"])
        add_yaml_closure(package_root, member_path, out)

    for analysis in data.get("analyses", data.get("interpretations", [])) or []:
        if not isinstance(analysis, dict):
            continue
        if isinstance(analysis.get("path"), str):
            add_file(package_root, resolve_package_path(package_root, yaml_file, analysis["path"]), out)
        for derived in analysis.get("derived_from", []) or []:
            if isinstance(derived, str):
                add_yaml_closure(package_root, resolve_package_path(package_root, yaml_file, derived), out)
        for asset in analysis.get("assets", []) or []:
            if isinstance(asset, dict) and isinstance(asset.get("path"), str):
                add_file(package_root, resolve_package_path(package_root, yaml_file, asset["path"]), out)


def yaml_files(files: set[Path]) -> list[Path]:
    return sorted(path for path in files if path.suffix.lower() in {".yaml", ".yml"})


def collect_package_files(package_root: Path, manifest_path: Path, entrypoint: Path) -> set[Path]:
    files: set[Path] = {manifest_path.resolve()}
    add_yaml_closure(package_root, entrypoint, files)
    return files


def format_package_yamls(files: set[Path], formatter: Path, check_format: bool) -> None:
    yamls = yaml_files(files)
    if not yamls:
        return
    formatter = formatter.resolve()
    cmd = [str(formatter)]
    if check_format:
        cmd.append("--check")
    cmd.extend(str(path) for path in yamls)
    subprocess.run(cmd, check=True)


def validate_entrypoint(entrypoint: Path, entrypoint_schema: str, bioscript: Path) -> None:
    if entrypoint_schema == "bioscript:panel:1.0":
        subprocess.run([str(bioscript), "validate-panels", str(entrypoint)], check=True)
    elif entrypoint_schema == "bioscript:assay:1.0":
        subprocess.run([str(bioscript), "validate-assays", str(entrypoint)], check=True)
    elif entrypoint_schema in {"bioscript:variant:1.0", "bioscript:variant"}:
        subprocess.run([str(bioscript), "validate-variants", str(entrypoint)], check=True)
    elif entrypoint_schema == "bioscript:variant-catalogue:1.0":
        subprocess.run([str(bioscript), "validate-variants", str(entrypoint)], check=True)
    else:
        raise ValueError(f"{entrypoint} has unsupported schema: {entrypoint_schema}")


def package_entrypoint(manifest_path: Path, manifest: dict[str, Any], package_root: Path) -> Path:
    schema = manifest.get("schema")
    if schema == "bioscript:package:1.0":
        entrypoint_name = manifest.get("entrypoint")
        if not isinstance(entrypoint_name, str) or not entrypoint_name:
            raise ValueError(f"{manifest_path} is missing entrypoint")
        return resolve_package_path(package_root, manifest_path, entrypoint_name)
    if schema in RUNNABLE_SCHEMAS:
        return manifest_path
    raise ValueError(f"{manifest_path} is not a BioScript package or runnable manifest")


def package_enabled(manifest: dict[str, Any]) -> bool:
    schema = manifest.get("schema")
    if schema == "bioscript:package:1.0":
        return True
    package = manifest.get("package")
    return schema in RUNNABLE_SCHEMAS and isinstance(package, dict) and package.get("publish") is True


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def release_summary(manifest: dict[str, Any], files: list[Path], zip_path: Path, package_root: Path) -> dict[str, Any]:
    return {
        "schema": "bioscript:package-release:1.0",
        "version": "1.0",
        "name": manifest.get("name"),
        "label": manifest.get("label"),
        "package_version": manifest.get("version"),
        "kind": str(manifest.get("schema", "")).removeprefix("bioscript:").split(":")[0],
        "created_at": dt.datetime.now(dt.UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        "entrypoint": "manifest.yaml",
        "artifact": {
            "path": zip_path.name,
            "sha256": sha256_file(zip_path),
            "size_bytes": zip_path.stat().st_size,
        },
        "contents": {
            "files": [path.relative_to(package_root).as_posix() for path in files],
        },
    }


def write_release_yaml(manifest: dict[str, Any], files: list[Path], zip_path: Path, package_root: Path) -> Path:
    release_path = package_root / f"{package_root.name}.yaml"
    with release_path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(
            release_summary(manifest, files, zip_path, package_root),
            handle,
            sort_keys=False,
            allow_unicode=True,
            width=120,
        )
    return release_path


def find_exr_dir(start: Path) -> Path | None:
    for directory in (start.resolve(), *start.resolve().parents):
        candidate = directory / "exr"
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return directory
    return None


def run_exbuild(package_root: Path) -> None:
    exbuild_path = package_root / EXBUILD_FILENAME
    if not exbuild_path.exists():
        return
    config = load_yaml(exbuild_path)
    steps = config.get("steps") or []
    if not isinstance(steps, list):
        raise ValueError(f"{exbuild_path}: 'steps' must be a list")

    env = os.environ.copy()
    exr_dir = find_exr_dir(package_root)
    if exr_dir is not None:
        env["PATH"] = f"{exr_dir}{os.pathsep}{env.get('PATH', '')}"

    for index, step in enumerate(steps):
        if not isinstance(step, dict) or not isinstance(step.get("run"), str):
            raise ValueError(f"{exbuild_path}: step {index} must have a string 'run'")
        name = step.get("name", step["run"])
        print(f"build [{package_root.name}]: {name}")
        subprocess.run(step["run"], shell=True, cwd=package_root, env=env, check=True)


def package_manifest(
    manifest_path: Path,
    bioscript: Path,
    formatter: Path,
    dry_run: bool,
    check_format: bool,
) -> Path:
    manifest_path = manifest_path.resolve()
    package_root = manifest_path.parent
    run_exbuild(package_root)
    manifest = load_yaml(manifest_path)
    entrypoint = package_entrypoint(manifest_path, manifest, package_root)

    files = collect_package_files(package_root, manifest_path, entrypoint)
    format_package_yamls(files, formatter, check_format)

    manifest = load_yaml(manifest_path)
    entrypoint = package_entrypoint(manifest_path, manifest, package_root)
    files = collect_package_files(package_root, manifest_path, entrypoint)
    entrypoint_data = load_yaml(entrypoint)
    validate_entrypoint(entrypoint, str(entrypoint_data.get("schema")), bioscript)

    zip_path = package_root / f"{package_root.name}.zip"
    rel_files = sorted(path.relative_to(package_root) for path in files)
    print(f"{manifest_path}: {len(rel_files)} files -> {zip_path}")
    for rel in rel_files:
        print(f"  {rel}")
    if dry_run:
        return zip_path

    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for rel in rel_files:
            archive.write(package_root / rel, rel.as_posix())
    release_path = write_release_yaml(entrypoint_data, sorted(files), zip_path, package_root)
    format_package_yamls({release_path}, formatter, check_format=False)
    print(f"release: {release_path}")
    return zip_path


def iter_package_manifests(root: Path) -> list[Path]:
    manifests: list[Path] = []
    for path in root.rglob("manifest.yaml"):
        if any(part in EXCLUDE_DIRS for part in path.parts):
            continue
        try:
            data = load_yaml(path)
        except Exception:
            continue
        if package_enabled(data):
            manifests.append(path)
    return sorted(manifests)


def main() -> int:
    parser = argparse.ArgumentParser(description="Build BioScript package zips from manifest.yaml files.")
    parser.add_argument("paths", nargs="*", type=Path, help="Package manifest files or directories to scan.")
    parser.add_argument("--bioscript", type=Path, default=Path("bioscript/bs"))
    default_formatter = Path(__file__).resolve().parents[1] / "format-yaml.sh"
    parser.add_argument("--formatter", type=Path, default=default_formatter)
    parser.add_argument("--check-format", action="store_true", help="check YAML formatting instead of writing changes")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    roots = args.paths or [Path("assays")]
    manifests: list[Path] = []
    for path in roots:
        if path.is_file():
            manifests.append(path)
        else:
            manifests.extend(iter_package_manifests(path))
    if not manifests:
        print("No bioscript package manifest.yaml files found.", file=sys.stderr)
        return 1

    for manifest in sorted(set(manifests)):
        package_manifest(manifest, args.bioscript, args.formatter, args.dry_run, args.check_format)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
