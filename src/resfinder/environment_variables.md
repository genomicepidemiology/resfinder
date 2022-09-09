# Environment Variables for ResFinder

Environment variables recognized by ResFinder, the flag they replace and the default value for the flag. Provided command line flags will always take precedence. Set environment variables takes precedence over default flag values.

Additional Environment variables can be added by appending entries to the table below. The 'Flag' entry in the table must be the double dash flag recognised by ResFinder. The 'Default Value' entry is just for information.

## Environment Variables Table

| Environment Variabel       | Flag                | Default Value  |
|----------------------------|---------------------|----------------|
| CGE_KMA                    | kmaPath             | kma            |
| CGE_BLASTN                 | blastPath           | blastn         |
| CGE_RESFINDER_RESGENE_DB   | db_path_res         | None           |
| CGE_RESFINDER_RESPOINT_DB  | db_path_point       | None           |
| CGE_RESFINDER_GENE_COV     | min_cov             | 0.60           |
| CGE_RESFINDER_GENE_ID      | threshold           | 0.80           |
| CGE_RESFINDER_POINT_COV    | min_cov_point       | 0.60           |
| CGE_RESFINDER_POINT_ID     | threshold_point     | 0.80           |
| CGE_DISINFINDER_DB         | db_path_disinf      | None           |
| CGE_DISINFINDER_DB_KMA     | db_path_disinf_kma  | kma            |
