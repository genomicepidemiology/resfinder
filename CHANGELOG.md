# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [4.2.2] - 2022-09-19

### Added
- ResFinder will now complain and exit if the ResFinder database is not found, as it is necessary, even if only looking for point mutations or disinfectant genes.

### Fixed
- Issue where the application failed when run using only the pointfinder option (--point)
- Changelog version format from d.d.d to [d.d.d]


## [4.2.1] - 2022-09-12

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
- Flag "--databases". ResFinder will now always be run against all databases. Option will be removed in the next major update.
- ResFinder will in the next major update not default to database paths within the application directory. Instead the use of the appropriate environment variables is recommended or the appropriate flags.
- The flag "-db_res" will in the next major update not be supported. Instead use "--db_path_res".
- The flag "-db_res_kma" will in the next major update not be supported.Instead use "--db_path_res_kma".
- The flag "-acq" will in the next major update not be supported.Instead use "--acquired".
- The flag "-ao" will in the next major update not be supported.Instead use "--acq_overlap".
- The flag "-db_disinf" will in the next major update not be supported.Instead use "--db_path_disinf".
- The flag "-db_disinf_kma" will in the next major update not be supported.Instead use "--db_path_disinf_kma".
- The flag "-db_point" will in the next major update not be supported.Instead use "--db_path_point".
- The flag "-db_point_kma" will in the next major update not be supported.Instead use "--db_path_point_kma".
- The flag "-l_p" will in the next major update not be supported.Instead use "--min_cov_point".
- The flag "-t_p" will in the next major update not be supported.Instead use "--threshold_point".

### Fixed
- Issue in PointFinder where a phenotype depending on several mutations would not be written in the phenotypes results files.
- Output in PointFinder, where some antibiotics would be listed twice.
- Issue in json ouptput where unknown mutations were listed even if option wasn't enabled.


## [4.2.0] - 2022-04-21 [YANKED]
