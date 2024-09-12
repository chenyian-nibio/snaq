# Changelog

## [Unreleased]

## [1.1.0] - 2024-09-11

### Changed

- To make the rule rarefy optional in the pipeline, the rrf option in the target must come before the cls option.
- Instead fetching pre-formatted taxonpath.json and names.json, these files will be created based on the download contents from NCBI taxonomy database.
- Improve the taxon_id mapping for processing the taxonomy in the result data.

## [1.0.1] - 2023-02-07

### Added

 - Yaml files for QIIME2 version 2023.2 for both linux and OSX.

### Changed

- Order of conda channels in phyloseq.yaml to make sure they work smoothly with new updates.
- Use the new QIIME2 version 2023.2 as the default environment in Snakefile.
