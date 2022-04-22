# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Several environmental variables for ResFinder settings (see README.md).
- Flags "--ignore_indels" and "--ignore_stop_codons" that will make the point mutation algorithm ignore indels and premature stop codons, respectively.
- Feature to search for genes that provide resistance to disinfectants (--disinfectant).
- Nanopore flag (--nanopore) that will use different settings for KMA optimized for Nanopore data.
- New method for calling ResFinder if installed via pip "python -m resfinder -h"

### Changed
- Recommended installation method for ResFinder. See README.md.
- It is no longer necessary to have cloned ResFinder via git in order to obtain the version number.
- json output file to contain all output results and enable the user to specify a path for the file.

### Deprecated
- It is no longer recommended to clone the repository of ResFinder, unless you are a developer. Instead install via pip is recommended.

### Removed
- Flag "--databases". ResFinder will now always be run against all databases.

## 2.4.0 - 2022-04-21 [YANKED]
