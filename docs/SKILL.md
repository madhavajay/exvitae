# ExVitae Assay And Panel Skill

Use this skill when creating, editing, validating, testing, or packaging BioScript assays, panels, variant catalogues, and report tests in the ExVitae repository.

The full human guide is [How To Make A Panel And An Assay](how-to-make-a-panel-and-assay.md). This file is the concise agent workflow.

## Toolchain

ExVitae uses three layers:

- Rust: BioScript runtime, CLI, schema validation, file-format inspection, reporting, and tests live under `bioscript/rust/`.
- Python: analysis scripts such as `apol1.py`, repository tools under `tools/`, and report-test orchestration under `tools/test_reports.py`.
- Bash scripts: stable entrypoints such as `bioscript/bs`, `bioscript/test.sh`, `package.sh`, `test-report.sh`, and `bioscript/tools/fetch_test_data.sh`.

BioScript analysis code is Python-like but runs through Monty, a Rust-hosted Pythonic runtime with unrestricted Python capabilities removed. Treat it as a controlled analysis sandbox: inputs, outputs, filesystem access, and networking should only happen through BioScript-approved APIs. Assay code must not assume it can make network calls, read arbitrary host files, exfiltrate data, or perform side effects outside the sandbox.

BioScript/ExVitae also includes or wraps common bioinformatics tools around the runtime, such as `samtools`, `bcftools`, and file-preparation utilities. Exact availability is version/package dependent; assume the newer bundled toolchain is evolving and verify commands locally before documenting them as required.

Check the runtime before changing assay logic:

```bash
git submodule update --init --recursive bioscript
```

```bash
bioscript/test.sh
```

Run a fast smoke test:

```bash
cd bioscript
./bs bioscripts/hello-world.py
```

## Where To Look First

Schema docs:

- `bioscript/docs/variant-schema.md`
- `bioscript/docs/assay-schema.md`
- `bioscript/docs/panel-schema.md`
- `bioscript/docs/variant-catalogue.md`
- `bioscript/docs/file-formats.md`

Example assays and panels:

- `assays/risk/APOL1/`
- `assays/risk/longevity/`
- `assays/risk/pcsk9-ldl/`
- `assays/pgx/glp1/`
- `assays/pgx/pgx-1/`

Copy local patterns before inventing new structure. Small assays should use hand-written variant YAML. Large assays can use TSV-backed catalogues.

## Build Order

For a small assay:

1. Create `assays/<category>/<name>/`.
2. Add one `*.yaml` per variant with `schema: bioscript:variant:1.0`.
3. Validate variants with `../exr validate-variants <dir>`.
4. Add `assay.yaml` or package-entrypoint `manifest.yaml` with `schema: bioscript:assay:1.0`.
5. Add `members` pointing to variant YAML files.
6. Add Python analysis only when a derived result is needed.
7. Add `analyses[].derived_from` and `analyses[].emits`.
8. Validate the assay with `../exr validate-assays <manifest>`.
9. Run it against at least one local test input with `bioscript/bs` or `./test-report.sh`.
10. Add `package.publish: true` or a package wrapper `manifest.yaml`.
11. Package with `./package.sh <manifest.yaml>`.
12. Run the built zip with `bioscript/bs report <zip> ...`.

For a panel:

1. Create `panel.yaml` or `manifest.yaml` with `schema: bioscript:panel:1.0`.
2. Add `members` entries for `variant`, `assay`, or `variant-catalogue`.
3. Keep packageable member paths inside the package root, commonly under `assets/`.
4. Validate with `../exr validate-panels <panel-or-manifest>`.
5. Package with `./package.sh <manifest.yaml>` when `package.publish: true` is present.

For a catalogue:

1. Create `catalogue/variants.yaml`.
2. Create `catalogue/variants.tsv`.
3. Create `catalogue/findings.tsv` if findings are needed.
4. Map every executable TSV column in `catalogue/variants.yaml`.
5. Validate with `./exr variant-catalogue validate` from the workspace root.
6. Build generated variant YAML into `/private/tmp/...`.
7. Validate the generated YAML with `./exr validate-variants`.
8. Attach TSV files as analysis `assets` if Python needs them.

## Validation Commands

Run from `exvitae/`:

```bash
../exr validate-variants assays/risk/APOL1
../exr validate-assays assays/risk/APOL1/manifest.yaml
../exr validate-panels assays/risk/longevity/manifest.yaml
./package.sh --dry-run assays/risk/APOL1/manifest.yaml
```

For catalogues, run from the workspace root:

```bash
./exr variant-catalogue validate projects/apol1-catalogue-sandbox/catalogue/variants.yaml
./exr variant-catalogue build projects/apol1-catalogue-sandbox/catalogue/variants.yaml /private/tmp/apol1-catalogue-generated --clean
./exr validate-variants /private/tmp/apol1-catalogue-generated
```

## Test Data

Public/local-cache test data is fetched by:

```bash
bioscript/tools/fetch_test_data.sh --check
bioscript/tools/fetch_test_data.sh --list
bioscript/tools/fetch_test_data.sh --dataset 23andme
bioscript/tools/fetch_test_data.sh --dataset 1k-genomes --only "*.tbi,*.crai,*.fai"
```

The fetcher reads `bioscript/tools/sources.yaml`, stores files in `~/.bioscript/cache/test-data` by default, and materializes symlinks under `test-data/`.

## Report Testing

Two report-test paths exist:

- `./test-report.sh <data-alias> <assay-or-panel>` is the quick alias-based manual runner.
- `tools/test_reports.py` reads `samples.yaml` and optional `samples.private.yaml` for repeatable sample suites.

Use `test-report.sh` for a focused manual run:

```bash
./test-report.sh --list
./test-report.sh 23andme_v5 apol1 --no-open -- --html
./test-report.sh NA06985-vcf longevity --no-open
./test-report.sh 23andme_v5 assays/risk/APOL1/manifest.yaml --no-open -- --html
```

Use the YAML-driven report suite when adding samples:

```bash
./test.sh --list
./test.sh --sample 23andme-v5-hu50B3F5 --assay pgx-1
./test.sh --sample my-new-sample --assay APOL1
```

If `./test.sh` is not the report wrapper in the current branch, call the Python tool directly:

```bash
uv run --with pyyaml python tools/test_reports.py --list
uv run --with pyyaml python tools/test_reports.py --sample my-new-sample --assay APOL1
```

## Private Samples

Use `samples.private.yaml` for local/private genomes. It is intentionally ignored by git.

Minimal private sample:

```yaml
defaults:
  output_dir: test-output/report-tests-private
  analysis_max_duration_ms: 30000
  detect_sex: true
  html: true
  assertions:
    observations:
      min_rows: 1
    reports:
      min_rows: 1
      require_status:
        - complete
        - partial

assays:
  - assays/risk/APOL1/manifest.yaml

samples:
  - id: my-23andme
    input_file: /absolute/path/to/genome.zip

  - id: my-vcf
    input_file: /absolute/path/to/sample.vcf.gz
    input_index: /absolute/path/to/sample.vcf.gz.tbi
    input_format: vcf

  - id: my-cram
    input_file: /absolute/path/to/sample.cram
    input_index: /absolute/path/to/sample.cram.crai
    reference_file: /absolute/path/to/reference.fa
    reference_index: /absolute/path/to/reference.fa.fai
    sample_sex: male
```

Run only the private sample you added:

```bash
uv run --with pyyaml python tools/test_reports.py --sample my-23andme --assay APOL1
```

## Safety Rules

- Do not treat outputs as clinical results.
- Do not commit private genomes, private sample paths, or generated private reports.
- Keep variant identity, coordinates, alleles, findings, and provenance in YAML.
- Keep cross-variant derivation logic in Python.
- Always validate and run at least one report before packaging.
