# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [4.3.3] - 2023-08-22

### Fixed

- Issue where ResFinder crashed with the error: "UnboundLocalError: local variable \'hit_class\' referenced before assignment'" if the flag '--output_aln' was enabled along with FASTQ input

## [4.3.2] - 2023-05-31

### Fixed

- Issue where PointFinder part crashed with the error: "NameError: name 'combined_id' is not defined"

## [4.3.1] - 2023-03-03

### Removed

- Overlap warning from pointfinder (message was only intended for testing)

### Fixed

- dockerfile, so that it works with NextFlow
- Issue with alignment output in json
- Issue with reporting concatenated gene hits in PointFinder

## [4.3.0] - 2023-02-09

### Added

- JSON output now contains the 'software_exec' class. It contains the command used to invoke ResFinder and the value of all parameters.
- JSON output now contains 'grades' for 'seq_regions' and 'phenotypes' corresponding to the colors known from the web tool: dark green, green, grey is 3, 2, 1, respectively. Phenotypes can also be 0 if no resistance features was found. Phenotype grades will as for the webtool be the highest grade of the resistance features found.

### Changed

- It is no longer necessary to set the flag --ignore_missing_species if '--species other' is set

### Fixed

- Case where species "other" is chosen and pointfinder is run causing error phenodb.mut_type_is_defined not found
- Alignments in json output to include missing gene parts in query string

## [4.2.5] - 2023-01-23

### Fixed

- Issue with unknown mutations not beeing showed and ignored correctly.
- Issue with Identical seq_region hits in Pointfinder

## [4.2.4] - 2023-01-17

### Added

- Option (--output_aln) that will output alignment and match sequences in the json output ("seq_regions": "query_string" / "alignment_string" / "ref_string" )

### Deprecated

- The Pointfinder_prediction.txt will no longer be a part of the PointFinder output.

### Fixed

- Issue with RNA mutations so the genes in RNA_genes.txt is read correctly and mutations found.
- Issue with indels such that insertions provided both as nucleotide and amino acids will be found.

## [4.2.3] - 2022-10-13

### Added

- The ResFinder databases now has a VERSION file that will be used to determine database versions. Implemented the use of these if the python module cgelib has been updated to at least version 0.7.3

### Fixed

- Issue where the application failed when species was specified but did not have a phenotype panel.
- DisinFinder overwriting ResFinder results.
- Issue with phenotype results showing integers or nothing instead of gene names in the 'Genetic background' column.
- Dockerfile

## [4.2.2] - 2022-09-19

### Added

- ResFinder will now complain and exit if the ResFinder database is not found, as it is necessary, even if only looking for point mutations or disinfectant genes.

### Fixed

- Issue where the application failed when run using only the pointfinder option (--point)
- Issue where application would crash if using the --disinfectant and --nanopore flags.
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
