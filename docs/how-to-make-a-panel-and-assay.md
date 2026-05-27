# How To Make A Panel And An Assay

This guide shows how to create a small BioScript assay, wrap it in a panel, verify it, and package it as a zip.

The tutorial uses APOL1 because it is small enough to understand by hand but still has real assay shape: two SNVs form the G1 marker, one deletion forms the G2 marker, and a Python BioScript analysis turns the observed genotypes into one derived status.

Use this repository for assay development, validation, research, education, and informational workflows only. Do not treat example outputs as medical advice or as a clinical result.

AI agents: use [SKILL.md](SKILL.md) as the agent-skill entrypoint for this workflow. This human guide explains the concepts; the skill file should contain the concise agent recipe, checks, and repo-specific commands.

## Quick Mental Model

BioScript packages are plain files:

- `variant.yaml`: one genomic variant definition.
- `assay.yaml`: one named test that observes variants and may run Python analysis logic.
- `panel.yaml`: a collection of variants and assays.
- `*.py`: BioScript-compatible Python analysis code.
- `manifest.yaml`: package entrypoint metadata used by the zip packager.
- `catalogue/variants.yaml` plus TSV files: compact TSV-backed variant catalogue for larger assays.
- `*.zip`: distributable package built from the manifest and all referenced files.

Typical shape:

```text
APOL1/
  manifest.yaml          # assay manifest and package entrypoint
  rs73885319.yaml        # variant: APOL1 G1 site 1
  rs60910145.yaml        # variant: APOL1 G1 site 2
  rs71785313.yaml        # variant: APOL1 G2 deletion
  apol1.py               # derived APOL1 status logic
  APOL1.zip              # packaged artifact
```

APOL1 assay flow:

```text
manifest.yaml
  |
  |-- members -> rs73885319.yaml, rs60910145.yaml, rs71785313.yaml
  |
  `-- analyses -> apol1.py -> output TSV/report rows
```

The manifest is the entrypoint. It points to the files the runtime needs: variant YAML files under `members`, optional assay or catalogue members in larger panels, and Python analysis scripts under `analyses`.

Execution order is:

1. Read the manifest.
2. Resolve and observe all variant members first.
3. Pass the collected observations to assay or panel analyses such as `apol1.py`.
4. Merge variant observations, analysis output, findings, and provenance into report artifacts.

Small assays should start with hand-written `variant.yaml` files. Large catalogues can use TSV-backed catalogues and optionally generate per-variant YAML during packaging or compatibility work.

## Toolchain And Reference Files

ExVitae uses three layers:

- Rust builds the BioScript runtime, CLI, schema validators, file-format inspection, reporting, and most low-level tests. The Rust workspace is under `bioscript/rust/`.
- Python is used for assay analysis scripts and repository tooling. Assay scripts such as `apol1.py` run inside BioScript, while helper tools live under `tools/`.
- Bash scripts provide stable local commands, including `bioscript/bs`, `bioscript/test.sh`, `package.sh`, `test-report.sh`, and `bioscript/tools/fetch_test_data.sh`.

BioScript analysis code is Python-like, but it is not normal unrestricted Python. It runs through Monty, a Rust-hosted Pythonic runtime, with dangerous or uncontrolled capabilities removed. The point is to let people run shared genomic analysis code while keeping control over input, output, filesystem access, and networking.

In practice that means an assay script should only see the files and helper APIs BioScript deliberately exposes. It should not be able to make arbitrary network calls, read unrelated host files, exfiltrate genetic data, or perform side effects outside the sandbox without explicit host support or permission.

BioScript/ExVitae also includes or wraps common bioinformatics tooling used around assays, such as `samtools`, `bcftools`, and related file-preparation utilities. Exact tool availability depends on the BioScript version and packaging target; a newer bundled toolchain is expected to make more of these standard genomics tools available through controlled runtime paths.

Schema docs live in:

```text
bioscript/docs/variant-schema.md
bioscript/docs/assay-schema.md
bioscript/docs/panel-schema.md
bioscript/docs/variant-catalogue.md
bioscript/docs/file-formats.md
```

Good examples to copy from:

```text
assays/risk/APOL1/
assays/risk/longevity/
assays/risk/pcsk9-ldl/
assays/pgx/glp1/
assays/pgx/pgx-1/
```

## Setup And Hello World

Before adding genetics, prove that the toolchain can build BioScript and that the runtime can run a simple script.

Install Rust if needed:

```bash
rustup --version
```

If `rustup` is missing, install the Rust toolchain from `https://rustup.rs/`, then open a new shell and check:

```bash
rustc --version
cargo --version
```

From the ExVitae repository, run the BioScript test suite:

```bash
git submodule update --init --recursive bioscript
```

```bash
bioscript/test.sh
```

The submodule command makes sure the nested BioScript checkout exists. The test command compiles the Rust workspace and runs the core BioScript tests for file formats, CLI behavior, schema validation, source-size checks, runtime security, and resource coverage.

For a faster smoke test, run the hello-world script from inside the `bioscript` directory:

```bash
cd bioscript
./bs bioscripts/hello-world.py
```

Expected output:

```text
    Finished `dev` profile [optimized + debuginfo] target(s) in 0.25s
     Running `rust/target/debug/bioscript bioscripts/hello-world.py`
hello from bioscript
2 + 3 = 5
loaded: sample input for bioscript
```

The first two Cargo lines can vary slightly by machine and build cache. The important part is that the script prints `hello from bioscript`, computes `2 + 3 = 5`, and loads `sample input for bioscript`.

That example lives at:

```text
bioscript/bioscripts/hello-world.py
bioscript/bioscripts/input.txt
```

Run it from `bioscript/` because the example intentionally reads `bioscripts/input.txt` as a relative path. Running `bioscript/bs bioscript/bioscripts/hello-world.py` from the repository root builds the runtime, but the sample then looks for `bioscripts/input.txt` in the wrong place.

It demonstrates the minimum runtime idea:

- Python syntax is used for assay logic.
- The runtime exposes safe file helpers such as `read_text` and `write_text`.
- Output is written into a controlled path rather than arbitrary host locations.

For real assays, prefer the structured helpers shown later:

```python
genotypes = bioscript.load_genotypes(input_file)
bioscript.write_tsv(output_file, rows)
```

## Recommended File Layout

For a small APOL1 assay:

```text
assays/risk/APOL1/
  manifest.yaml
  apol1.py
  rs73885319.yaml
  rs60910145.yaml
  rs71785313.yaml
```

For a panel that includes APOL1 plus other assays:

```text
assays/risk/kidney-risk/
  manifest.yaml
  panel.yaml
  assets/APOL1/manifest.yaml
  assets/APOL1/apol1.py
  assets/APOL1/rs73885319.yaml
  assets/APOL1/rs60910145.yaml
  assets/APOL1/rs71785313.yaml
```

The existing APOL1 example is already present at:

```text
assays/risk/APOL1/
```

Use that as the concrete reference while reading this guide.

Note: many assay directories in this repository use `assay.yaml` for the assay manifest and a separate `manifest.yaml` for package metadata. The current APOL1 example uses `manifest.yaml` as both the assay manifest and package entrypoint.

## Step 1: Write Variant YAML

A variant file is the atomic unit BioScript can look up in a genotype, VCF, CRAM, or BAM-derived input.

It tries to describe the same biological variant across the formats you may see in the wild:

- `coordinates` can include multiple assemblies, usually `grch37` and `grch38`.
- `alleles.ref` is the curated reference allele for this variant record.
- `alleles.observed_alts` is the broader set of alternate alleles that may appear at the site in upstream sources such as dbSNP.
- `alleles.alts` is the curated subset this assay cares about and can bind findings to.
- `findings[].alt` binds interpretation to one curated alt, or to `"*"` when the finding applies to any alternate allele at the site.

Minimal APOL1 G1 site 1 example:

```yaml
schema: "bioscript:variant:1.0"
version: "1.0"
name: "APOL1-G1-site-1-rs73885319"
gene: "APOL1"
summary: "APOL1 G1 haplotype site 1 SNV used by the APOL1 classifier."

identifiers:
  rsids:
    - "rs73885319"

coordinates:
  grch37:
    chrom: "22"
    pos: 36661906
  grch38:
    chrom: "22"
    pos: 36265860

alleles:
  kind: "snv"
  ref: "A"
  alts:
    - "G"
  observed_alts:
    - "G"

findings:
  - schema: "bioscript:risk:1.0"
    alt: "G"
    notes: "APOL1 G1 site 1; interpreted together with APOL1 G1 site 2 and G2 deletion."
```

Important format rules:

- `schema` must be `bioscript:variant:1.0`.
- Use `pos` for single-base variants.
- Use `start` and `end` for spans such as deletions.
- Store biological alleles, not display-only symbols, in variant YAML.
- Include both `grch37` and `grch38` when known so the same assay can run on multiple input assemblies.
- Use `observed_alts` for the full known alternate set at the site, and `alts` for the curated alternates that should trigger assay/report bindings.
- Keep evidence, notes, source links, and allele-specific findings in YAML rather than hiding them in Python.

APOL1 G2 is a 6 bp deletion, so its YAML uses a span:

```yaml
alleles:
  kind: "deletion"
  ref: "TTATAA"
  alts:
    - "<DEL:6>"
  deletion_length: 6
  motifs:
    - "TTATAA"
    - "ATAATT"
```

Deletion and repeat-context variants may have equivalent shifted representations. Use `equivalent_coordinates` and notes when representation is ambiguous.

Verify variant files:

```bash
../exr validate-variants assays/risk/APOL1
```

or directly with the BioScript binary:

```bash
bioscript/bs validate-variants assays/risk/APOL1
```

## Step 2: Write The Assay Manifest

An assay groups variants and declares any derived analysis logic.

Minimal APOL1 assay:

```yaml
schema: "bioscript:assay:1.0"
version: "1.0"
name: "APOL1"
label: "APOL1 Risk Assay"
tags:
  - "type:risk"
  - "gene:APOL1"

members:
  - kind: "variant"
    path: "rs73885319.yaml"
    version: "1.0"
  - kind: "variant"
    path: "rs60910145.yaml"
    version: "1.0"
  - kind: "variant"
    path: "rs71785313.yaml"
    version: "1.0"

analyses:
  - id: "apol1_status"
    kind: "bioscript"
    path: "apol1.py"
    output_format: "tsv"
    label: "APOL1 risk genotype"
    derived_from:
      - "rs73885319.yaml"
      - "rs60910145.yaml"
      - "rs71785313.yaml"
    emits:
      - key: "apol1_status"
        label: "APOL1 status"
        value_type: "string"
        format: "badge"
```

Format notes:

- `members` tells BioScript which variants to observe.
- `analyses` tells BioScript which Python files produce derived report columns.
- `derived_from` should name the variant or catalogue files the script depends on.
- `emits` documents the output columns so reports can label them.
- Keep variant identity in YAML and cross-variant classification in Python.

Verify the assay:

```bash
../exr validate-assays assays/risk/APOL1/manifest.yaml
```

## Step 3: Write The Python Analysis

The APOL1 script in `assays/risk/APOL1/apol1.py` defines query variants, loads the participant's genotype file, and writes one TSV row.

You will currently see variant details repeated in `bioscript.variant(...)` calls inside the Python script. That gives the script enough information to query a genotype input directly, and it can also fetch or resolve variant details dynamically when the YAML observation pipeline did not already pass the needed observation through. This duplication is a transitional shape and should be refactored so analysis code can rely more directly on manifest/YAML-derived observations.

The essential shape is:

```python
G1_SITE_1 = bioscript.variant(
    rsid="rs73885319",
    grch37="22:36661906-36661906",
    grch38="22:36265860-36265860",
    ref="A",
    alt="G",
    kind="snp",
)


def main():
    genotypes = bioscript.load_genotypes(input_file)
    site1 = genotypes.lookup_variant(G1_SITE_1)
    rows = [{
        "participant_id": participant_id,
        "apol1_status": "example",
        "rs73885319": site1 if site1 is not None else "missing",
    }]
    bioscript.write_tsv(output_file, rows)
```

Runtime variables and helpers:

- `input_file`: genotype, VCF, CRAM, BAM, or zip input path.
- `output_file`: primary output path.
- `participant_id`: participant or sample identifier.
- `bioscript.load_genotypes(...)`: loads supported genotype formats.
- `genotypes.lookup_variant(...)`: returns the observed genotype for one variant.
- `bioscript.write_tsv(...)`: writes list-of-dict rows as TSV.

For multi-variant assays, build a query plan:

```python
APOL1_QUERY_PLAN = bioscript.query_plan([
    G1_SITE_1,
    G1_SITE_2,
    G2_SITE,
])

site1, site2, g2 = genotypes.lookup_variants(APOL1_QUERY_PLAN)
```

APOL1-specific caution:

- G1 is a haplotype represented by two SNVs.
- G2 is a deletion that can be represented differently across file formats.
- Phase, missingness, ancestry context, and clinical interpretation require careful review.
- The tutorial logic is suitable for testing the assay pathway, not for unsupervised clinical use.

## Step 4: Make A Panel

A panel is a runnable collection. It can include variants directly and can include assays.

Minimal panel that includes APOL1:

```yaml
schema: "bioscript:panel:1.0"
version: "1.0"
name: "kidney-risk"
label: "Kidney Risk Research Panel"
tags:
  - "type:risk"
  - "topic:kidney"

members:
  - kind: "assay"
    path: "assets/APOL1/manifest.yaml"
    version: "1.0"
```

Format notes:

- Use a panel when you want one manifest to run or package several items.
- Use `kind: assay` when the member has its own analysis logic.
- Use `kind: variant` for simple independent observations.
- Use `kind: variant-catalogue` for TSV-backed catalogues when supported by the target runtime.
- For a packageable panel, keep member paths inside the panel package root, for example under `assets/`.

Verify the panel:

```bash
../exr validate-panels assays/risk/kidney-risk/panel.yaml
```

## Step 5: Run The Assay Or Panel

First fetch or check the public fixture data:

```bash
bioscript/tools/fetch_test_data.sh --check
bioscript/tools/fetch_test_data.sh --dataset 23andme
bioscript/tools/fetch_test_data.sh --dataset 1k-genomes --only "*.vcf.gz,*.tbi,*.cram,*.crai,*.fa,*.fai"
```

These fixtures are open data pulled from the [OpenMined `biovault-data` repository](https://github.com/OpenMined/biovault-data) and related public sources, including Kaggle, the Personal Genome Project (PGP), and the 1000 Genomes Project (1KGP). The fetch manifest in this repo records the exact URLs used.

Then inspect the real fixture files so you know what BioScript thinks they are:

```bash
bioscript/bs inspect test-data/23andme/v5/hu50B3F5/genome_hu50B3F5_v5_Full.zip

bioscript/bs inspect test-data/1k-genomes/vcf/NA06985.clean.vcf.gz \
  --input-index test-data/1k-genomes/vcf/NA06985.clean.vcf.gz.tbi

bioscript/bs inspect test-data/1k-genomes/aligned/NA06985.final.cram \
  --input-index test-data/1k-genomes/aligned/NA06985.final.cram.crai \
  --reference-file test-data/1k-genomes/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  --reference-index test-data/1k-genomes/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
```

`inspect` reports file family, container type, inferred assembly when available, phasing when available, selected zip entry, index status, source/vendor hints, and warnings.

Run a single manifest against one input:

```bash
bioscript/bs assays/risk/APOL1/manifest.yaml \
  --input-file test-data/23andme/v5/hu50B3F5/genome_hu50B3F5_v5_Full.zip \
  --output-file /tmp/apol1.tsv
```

Run a report against 23andMe v5:

```bash
bioscript/bs report assays/risk/APOL1/manifest.yaml \
  --input-file test-data/23andme/v5/hu50B3F5/genome_hu50B3F5_v5_Full.zip \
  --output-dir /tmp/apol1-23andme-v5-report \
  --html
```

Run a report against the NA06985 VCF:

```bash
bioscript/bs report assays/risk/APOL1/manifest.yaml \
  --input-file test-data/1k-genomes/vcf/NA06985.clean.vcf.gz \
  --input-index test-data/1k-genomes/vcf/NA06985.clean.vcf.gz.tbi \
  --input-format vcf \
  --output-dir /tmp/apol1-na06985-vcf-report \
  --html
```

Run a report against the NA06985 CRAM:

```bash
bioscript/bs report assays/risk/APOL1/manifest.yaml \
  --input-file test-data/1k-genomes/aligned/NA06985.final.cram \
  --input-index test-data/1k-genomes/aligned/NA06985.final.cram.crai \
  --reference-file test-data/1k-genomes/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  --reference-index test-data/1k-genomes/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai \
  --allow-md5-mismatch \
  --output-dir /tmp/apol1-na06985-cram-report \
  --html
```

Run the repository assay test harness:

```bash
./test-assays.sh assays/risk/APOL1
```

The test harness builds the local BioScript CLI when needed, runs assays against available local test data, and writes report artifacts under:

```text
test-output/
```

## Test Data And Report Tests

Public/local-cache test data is managed by:

```bash
bioscript/tools/fetch_test_data.sh --check
bioscript/tools/fetch_test_data.sh --list
bioscript/tools/fetch_test_data.sh --dataset 23andme
bioscript/tools/fetch_test_data.sh --dataset 1k-genomes --only "*.tbi,*.crai,*.fai"
```

The fetcher reads:

```text
bioscript/tools/sources.yaml
```

The fixture URLs point at open datasets such as the [OpenMined `biovault-data` repository](https://github.com/OpenMined/biovault-data), Kaggle-hosted data, Personal Genome Project data, and 1000 Genomes Project data. Treat the fetch manifest as the source of truth for exact filenames and upstream URLs.

By default it downloads into:

```text
~/.bioscript/cache/test-data
```

and materializes symlinks under:

```text
test-data/
```

There are two useful report-test workflows.

Use `test-report.sh` for a focused manual run against a named data alias:

```bash
./test-report.sh --list
./test-report.sh 23andme_v5 apol1 --no-open -- --html
./test-report.sh NA06985-vcf longevity --no-open
./test-report.sh 23andme_v5 assays/risk/APOL1/manifest.yaml --no-open -- --html
```

`test-report.sh` knows aliases for common inputs such as `23andme_v5`, `NA06985`, and `NA06985-vcf`, and aliases for assays or panels such as `apol1`, `apol1-zip`, `longevity`, and `pcsk9-ldl`.

Use the YAML-driven report suite for repeatable sample tests:

```bash
./test.sh --list
./test.sh --sample 23andme-v5-hu50B3F5 --assay pgx-1
```

That wrapper runs `tools/test_reports.py`, which reads:

```text
samples.yaml
samples.private.yaml
```

`samples.yaml` is for public or local-cache fixtures. `samples.private.yaml` is for local/private genomes and is intentionally ignored by git.

Minimal `samples.private.yaml`:

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
./test.sh --sample my-23andme --assay APOL1
```

If you need to call the Python runner directly:

```bash
uv run --with pyyaml python tools/test_reports.py --sample my-23andme --assay APOL1
```

## Step 6: Add A Package Manifest

The package builder looks for `manifest.yaml`.

For an assay package, `manifest.yaml` can be the assay itself if the assay includes `package.publish: true`, or it can be a package wrapper. The current APOL1 package uses the assay manifest shape directly:

```yaml
schema: "bioscript:assay:1.0"
version: "1.0"
name: "APOL1"
label: "APOL1 Risk Assay"
package:
  schema: "bioscript:package:1.0"
  publish: true
  artifact: "APOL1.zip"
members:
  - kind: "variant"
    path: "rs73885319.yaml"
    version: "1.0"
```

Packaging rules:

- The package root is the directory containing `manifest.yaml`.
- Every `members[].path` is included.
- Every `analyses[].path` is included.
- Every `analyses[].derived_from[]` is included and recursively validated.
- Every `analyses[].assets[].path` is included.
- YAML `include` fields are followed recursively.
- Paths may not escape the package root.
- `.git`, cache, logs, and temp folders are ignored.

Build the zip:

```bash
./package.sh assays/risk/APOL1/manifest.yaml
```

Package output:

```text
assays/risk/APOL1/APOL1.zip
assays/risk/APOL1/APOL1.yaml
```

The generated release YAML records:

- package name and version
- package kind
- created timestamp
- zip path
- zip SHA-256
- zip size
- exact file list

Dry run without writing a zip:

```bash
./package.sh --dry-run assays/risk/APOL1/manifest.yaml
```

Check YAML formatting during package build:

```bash
./package.sh --check-format assays/risk/APOL1/manifest.yaml
```

## TSV-Backed Catalogues

Use a variant catalogue when a panel or assay has many variants and hand-written YAML becomes noisy.

Catalogue layout:

```text
catalogue/
  variants.yaml
  variants.tsv
  findings.tsv
  interpretation_rules.tsv
```

`catalogue/variants.yaml` declares the TSV files and maps columns to BioScript roles:

```yaml
schema: "bioscript:variant-catalogue:1.0"
version: "1.0"
name: "apol1-catalogue"

variants:
  source: "variants.tsv"
  format: "tsv"
  key: "variant_id"
  columns:
    variant_id:
      role: "id"
      required: true
    rsid:
      role: "identifier.rsid"
    aliases:
      role: "identifier.aliases"
      list_separator: "|"
    kind:
      role: "alleles.kind"
    ref:
      role: "alleles.ref"
    alts:
      role: "alleles.alts"
      list_separator: "|"
    grch38_chrom:
      role: "coordinates.grch38.chrom"
    grch38_pos:
      role: "coordinates.grch38.pos"
      type: "integer"

findings:
  source: "findings.tsv"
  format: "tsv"
  key: "variant_id"
  columns:
    variant_id:
      role: "variant.id"
      required: true
    schema:
      role: "finding.schema"
    alt:
      role: "finding.alt"
    summary:
      role: "finding.summary"
```

TSV rules:

- Use a header row.
- Use tab delimiters.
- Keep cell values as strings unless the catalogue manifest declares a typed column.
- Use `|` for list-valued cells when declared by `list_separator`.
- Keep source-only curation columns if useful, but executable tools should depend only on mapped columns.

Validate and expand a catalogue:

```bash
cd /Users/madhavajay/dev/exvitae-data
./exr variant-catalogue validate projects/apol1-catalogue-sandbox/catalogue/variants.yaml
./exr variant-catalogue build projects/apol1-catalogue-sandbox/catalogue/variants.yaml /private/tmp/apol1-catalogue-generated --clean
./exr validate-variants /private/tmp/apol1-catalogue-generated
```

The APOL1 catalogue sandbox shows the intended asset-based analysis pattern:

```python
observations = bioscript.read_tsv(bioscript.context["observations_file"])
variants = bioscript.read_tsv(bioscript.context["assets"]["variants"])
findings = bioscript.read_tsv(bioscript.context["assets"]["findings"])
rules = bioscript.read_tsv(bioscript.context["assets"]["interpretation_rules"])
```

Use catalogue assets for curation tables. Use the catalogue member to tell BioScript what variants to observe. Use Python to join observations to the curation tables and emit derived rows.

## Verification Checklist

Run these from `exvitae/` unless noted.

Validate variant YAML:

```bash
../exr validate-variants assays/risk/APOL1
```

Validate assay YAML:

```bash
../exr validate-assays assays/risk/APOL1/manifest.yaml
```

Validate panel YAML:

```bash
../exr validate-panels assays/risk/longevity/manifest.yaml
```

Inspect input files:

```bash
bioscript/bs inspect test-data/23andme/v5/hu50B3F5/genome_hu50B3F5_v5_Full.zip
bioscript/bs inspect test-data/1k-genomes/vcf/NA06985.clean.vcf.gz \
  --input-index test-data/1k-genomes/vcf/NA06985.clean.vcf.gz.tbi
```

Check or fetch public test data:

```bash
bioscript/tools/fetch_test_data.sh --check
bioscript/tools/fetch_test_data.sh --dataset 23andme
bioscript/tools/fetch_test_data.sh --dataset 1k-genomes --only "*.vcf.gz,*.tbi,*.cram,*.crai,*.fa,*.fai"
```

Run one assay:

```bash
bioscript/bs assays/risk/APOL1/manifest.yaml \
  --input-file test-data/23andme/v5/hu50B3F5/genome_hu50B3F5_v5_Full.zip \
  --output-file /tmp/apol1.tsv
```

Generate a report:

```bash
bioscript/bs report assays/risk/APOL1/manifest.yaml \
  --input-file test-data/1k-genomes/vcf/NA06985.clean.vcf.gz \
  --input-index test-data/1k-genomes/vcf/NA06985.clean.vcf.gz.tbi \
  --input-format vcf \
  --output-dir /tmp/apol1-na06985-vcf-report \
  --html
```

Run the local assay test harness:

```bash
./test-assays.sh assays/risk/APOL1
```

Run a focused report test:

```bash
./test-report.sh 23andme_v5 apol1 --no-open -- --html
```

Run a YAML-defined sample from `samples.yaml` or `samples.private.yaml`:

```bash
./test.sh --sample 23andme-v5-hu50B3F5 --assay pgx-1
```

Package:

```bash
./package.sh assays/risk/APOL1/manifest.yaml
```

Verify the built zip by importing or reporting against it:

```bash
bioscript/bs report assays/risk/APOL1/APOL1.zip \
  --input-file test-data/23andme/v5/hu50B3F5/genome_hu50B3F5_v5_Full.zip \
  --output-dir /tmp/apol1-zip-23andme-v5-report \
  --html
```

## Common Mistakes

- Putting all interpretation in Python and leaving YAML without evidence or allele identity.
- Using `I` and `D` as stored deletion alleles in `variant.yaml` instead of biological sequence or symbolic deletion alleles.
- Omitting `grch37` or `grch38` coordinates when the target input may use that assembly.
- Forgetting `derived_from`, which makes package dependency collection incomplete.
- Forgetting `emits`, which makes reports harder to label.
- Pointing package paths outside the package root.
- Treating `variant-catalogue` TSV assets as the same thing as analysis assets. The catalogue says what to observe; assets are files the analysis script reads.
- Testing only against one input format. At minimum, inspect and run representative text, VCF, and CRAM inputs when the assay claims to support them.

## Useful References In This Repo

- `bioscript/docs/variant-schema.md`
- `bioscript/docs/assay-schema.md`
- `bioscript/docs/panel-schema.md`
- `bioscript/docs/variant-catalogue.md`
- `bioscript/docs/file-formats.md`
- `bioscript/bioscripts/hello-world.py`
- `bioscript/bioscripts/apol1.py`
- `assays/risk/APOL1/`
