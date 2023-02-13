# PointFinder class tests

This is only testing a very small fraction of the run_resfinder.

## Setup

```python

>>> from src.resfinder.run_resfinder import *
>>> from src.resfinder.cge.config import Config

>>> class DummyArgs():
...     def __init__(self):
...         self.inputfasta = None
...         self.inputfastq = None
...         self.outputPath = "./tests/tmp_out/"
...         self.blastPath = None
...         self.kmaPath = None
...         self.species = None
...         self.ignore_missing_species = None
...         self.db_path_res = None
...         self.db_path_res_kma = None
...         self.databases = None
...         self.acquired = True
...         self.acq_overlap = None
...         self.min_cov = None
...         self.threshold = None
...         self.point = True
...         self.db_path_point = None
...         self.db_path_point_kma = None
...         self.specific_gene = None
...         self.unknown_mut = None
...         self.min_cov_point = None
...         self.threshold_point = None
...         self.ignore_indels = None
...         self.ignore_stop_codons = None
...         self.pickle = False
...         self.nanopore = False
...         self.out_json = None
...         self.disinfectant = False
...         self.db_path_disinf = None
...         self.db_path_disinf_kma = None
...         self.output_aln = False
...         self.species = "ecoli"

>>> args = DummyArgs()
>>> conf = Config(args)

```

## get_call_parameters(conf: Config) -> dict

```python

>>> parameters = get_call_parameters(conf=conf)
>>> parameters["species"]
'escherichia coli'
>>> parameters["acquired"]
True

```

## get_software_exec_res(conf: Config) -> dict

```python

>>> soft_exec = get_software_exec_res(conf=conf)
>>> soft_exec["type"]
'software_exec'
>>> len(soft_exec["key"])
40

```
