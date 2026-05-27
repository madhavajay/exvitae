# ExVitae

Making Genomics Research into Secure Code that runs Everywhere.

## Guides

- [How to make a panel and an assay](docs/how-to-make-a-panel-and-assay.md): human-friendly tutorial for creating BioScript variants, assays, panels, tests, and package zips.
- [Agent skill](docs/SKILL.md): concise agent-facing workflow for building and validating assays in this repository.

## Project Notice

ExVitae is intended for assay development, validation, research, educational, and informational use. The materials and outputs in this repository are not medical advice, are not intended to diagnose, treat, cure, prevent, or monitor any disease or health condition, and should not be used for those purposes. This repository should only be used by researchers or by qualified, licensed physicians and other appropriately credentialed healthcare professionals acting within the scope of their training, licensure, and applicable law.

If you use this repository with human genetic data:

- You must only use data that you have the legal right and appropriate authorization to access, process, and analyze.
- You are responsible for complying with applicable privacy, consent, data-protection, and genetic-information laws in your jurisdiction.
- You should not rely on repository outputs as a substitute for review by a qualified clinician, genetic counselor, or other licensed healthcare professional.

## Disclaimer

To the fullest extent permitted by law, ExVitae and repository contributors provide this repository, its assay definitions, and its generated outputs on an "as is" and "as available" basis, without warranties of accuracy, completeness, merchantability, fitness for a particular purpose, or non-infringement.

Software can contain bugs, assay logic can be incomplete or incorrect, input files and data formats can be malformed or mislabeled, and sequencing, genotyping, imputation, reference, or annotation data can contain errors or limitations. Genetic interpretation is context-dependent and can change over time as evidence evolves.

Any output generated from this repository must be independently reviewed and double-checked before it is relied on for any research, clinical, operational, or other consequential purpose. Any potentially meaningful finding, signal, or interpretation should be validated or verified through appropriate follow-up review, screening, confirmatory testing, or other independent methods before any action is taken.

Where possible, the source of input or reference data should be identified in the relevant file, dataset, or accompanying metadata. However, much of the upstream material used in development or testing may originate from public data available on the internet or other third-party sources, and ExVitae and repository contributors do not control, warrant, or accept responsibility for the accuracy, completeness, legality, licensing, provenance, or continued availability of that upstream data.


## Data Handling

This repository includes local tooling for processing genetic data files during development and testing. Anyone operating ExVitae is responsible for:

- handling uploaded or local genetic data securely;
- limiting access to authorized users only;
- defining retention and deletion practices appropriate for the deployment or workflow in use; and
- reviewing any third-party services used in the surrounding stack before processing real user data.

This README does not itself create a hosted-service privacy policy or terms of service. If code from ExVitae is deployed as a user-facing product, those documents should be published separately and aligned with the actual infrastructure, retention windows, billing flows, and support contacts in use.

## Repository Layout

- `assays/`: assay definitions for risk, pharmacogenomics, and trait-oriented panels.
- `bioscript/`: BioScript runtime, docs, and supporting tooling used by this repo.
- `test-data/`: local fixture data in formats such as VCF, CRAM, ZIP, and direct-to-consumer text exports sourced from open datasets on the internet.
- `tools/`: helper scripts for fetching test data, running validations, and benchmarking.
