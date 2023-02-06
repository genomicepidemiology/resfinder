# Test Isolate class

**Missing tests**

## Setup

### phenoDB object

```python

>>> import os
>>> from src.resfinder.cge.phenotype2genotype.res_profile import PhenoDB

>>> resfinder_db_path = os.environ["CGE_RESFINDER_RESGENE_DB"]
>>> assert(len(resfinder_db_path) > 0)
>>> pointfinder_db_path = os.environ["CGE_RESFINDER_RESPOINT_DB"]
>>> assert(len(pointfinder_db_path) > 0)

>>> abclassdef_file = "{}/antibiotic_classes.txt".format(resfinder_db_path)
>>> acquired_file = "{}/phenotypes.txt".format(resfinder_db_path)
>>> point_file = ("{}/escherichia_coli/phenotypes.txt"
...               .format(pointfinder_db_path))
>>> res_pheno_db = PhenoDB(abclassdef_file=abclassdef_file,
...                        acquired_file=acquired_file,
...                        point_file=point_file)

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

### Result object

More details about the setup of the result object can be found in the
std_results test documentation.

```python

# >>> rf_dat_kma = {}
# >>> rf_dat_kma["sbjct_header"] = "blaOXA-384_1_KF986263"
# >>> rf_dat_kma["perc_ident"] = 97
# >>> rf_dat_kma["HSP_length"] = 100
# >>> rf_dat_kma["sbjct_length"] = 90
# >>> rf_dat_kma["sbjct_start"] = 1
# >>> rf_dat_kma["sbjct_end"] = 90
# >>> rf_dat_kma["contig_name"] = "NA"
# >>> rf_dat_kma["query_start"] = "NA"
# >>> rf_dat_kma["query_end"] = "NA"
# >>> rf_dat_kma["perc_coverage"] = 100
# >>> rf_dat_kma["depth"] = 21
# 
# >>> rf_custom_kma = {}
# >>> rf_custom_kma["excluded"] = {}
# >>> rf_custom_kma["aminoglycoside"] = "No hit found"
# >>> rf_custom_kma["beta-lactam"] = {"unique_hit_key": rf_dat_kma}
>>> from cgecore.blaster.blaster import Blaster
>>> import os
>>> import tempfile

>>> resfinder_db_path = os.environ["CGE_RESFINDER_RESGENE_DB"]
>>> rf_custom_kma = Blaster(inputfile="tests/data/test_isolate_01.fa", databases=["beta-lactam"], db_path=resfinder_db_path, out_path="tests/tmp_out")
... #doctest: +ELLIPSIS
Found...

>>> gyrA_kma_hit = {}
>>> gyrA_kma_hit["sbjct_header"] = "gyrA_1_CP073768.1"
>>> gyrA_kma_hit["perc_ident"] = 99.92
>>> gyrA_kma_hit["HSP_length"] = 2628
>>> gyrA_kma_hit["sbjct_length"] = 2628
>>> gyrA_kma_hit["sbjct_start"] = 1
>>> gyrA_kma_hit["sbjct_end"] = 2628
>>> gyrA_kma_hit["contig_name"] = "NA"
>>> gyrA_kma_hit["query_start"] = "NA"
>>> gyrA_kma_hit["query_end"] = "NA"
>>> gyrA_kma_hit["perc_coverage"] = 100.0
>>> gyrA_kma_hit["depth"] = 21
>>> gyrA_kma_hit["mis_matches"] = [
...   [ 'sub', 81, 81, 'D', 'p.G81D', 'GGT', 'GAT', 'G', 'D' ],
...   [ 'sub', 82, 82, 'G', 'p.D82G', 'GAC', 'GGC', 'D', 'G' ] ]

>>> pf_custom_kma = {}
>>> pf_custom_kma["excluded"] = {}
>>> pf_custom_kma["gyrA"] = gyrA_kma_hit
>>> pf_custom_kma["gyrB"] = "No hit found"

>>> from cgelib.output.result import Result
>>> from src.resfinder.cge.output.gene_result import GeneResult

>>> res = Result.init_software_result(name="ResFinder", gitdir=".")
>>> res.init_database("ResFinder", ".")
>>> res.init_database("PointFinder", ".")

>>> from src.resfinder.cge.output.std_results import ResFinderResultHandler
>>> ResFinderResultHandler.standardize_results(res,
...                                            rf_custom_kma,
...                                            "ResFinder",
...                                            conf)

>>> from src.resfinder.cge.output.std_results import PointFinderResultHandler
>>> PointFinderResultHandler.standardize_results(res,
...                                              pf_custom_kma,
...                                              "PointFinder")

```

## Init Isolate

```python

>>> amr_panel_file = f"{resfinder_db_path}/phenotype_panels.txt"
>>> species = "Escherichia coli"

>>> from src.resfinder.cge.phenotype2genotype.isolate import Isolate
>>> isolate = Isolate(name="Test sample", species=species,
...                   amr_panel_file=amr_panel_file)

```

## Methods

### load_amr_panel(species, panel_file)

```python

>>> amr_set = Isolate.load_amr_panel("Campylobacter jejuni", amr_panel_file)
>>> assert("ciprofloxacin" in amr_set)
>>> assert("ampicillin" in amr_set)

```

### check_panel_name(name, panels)

```python

>>> panels = {"campylobacter", "campylobacter jejuni", "escherichia coli"}
>>> panel_match = Isolate.check_panel_name("campylobacter jejuni", None)
>>> f"{panel_match}"
'None'
>>> Isolate.check_panel_name("campylobacter jejuni", panels)
'campylobacter jejuni'
>>> Isolate.check_panel_name("campylobacter", panels)
'campylobacter'
>>> Isolate.check_panel_name("no match", panels)
False

```

### load_finder_results(std_table, phenodb, type)

```python

>>> isolate.load_finder_results(std_table=res, phenodb=res_pheno_db,
...                             type="seq_regions")
>>> isolate.load_finder_results(std_table=res, phenodb=res_pheno_db,
...                             type="seq_variations")

>>> isolate["blaB-2_AF189300"][0].ab_class.pop()
'beta-lactam'
>>> isolate["gyrA;;1;;CP073768.1_81_d"][0].ref_codon
'ggt'
>>> isolate["gyrA;;1;;CP073768.1_81_d"][0].mut_codon
'gat'
>>> isolate["gyrA;;1;;CP073768.1_81_d"][0].ab_class.pop()
'quinolone'
>>> isolate["gyrA;;1;;CP073768.1_82_g"][0].ref_codon
'gac'
>>> isolate["gyrA;;1;;CP073768.1_82_g"][0].mut_codon
'ggc'
>>> isolate["gyrA;;1;;CP073768.1_82_g"][0].ab_class.pop()
'quinolone'

```

### calc_res_profile(res_pheno_db)

```python

>>> isolate.calc_res_profile(res_pheno_db)

```

### get_phenodb_id(feat_res_dict, type)

```python

>>> feat_res_dict = res["seq_regions"]["blaB-2;;1;;AF189300"]
>>> Isolate.get_phenodb_id(feat_res_dict, "seq_regions")
('blaB-2_AF189300', '')

>>> feat_res_dict = res["seq_variations"]["gyrA;;1;;CP073768.1;;81;;d"]
>>> Isolate.get_phenodb_id(feat_res_dict, "seq_variations")
('', 'gyrA;;1;;CP073768.1_81_d')

```


## Private methods

### _get_antibiotics(line, panel_name, panels)

```python

>>> panels_dummy1 = {"key1": ["val1", "val2"]}
>>> Isolate._get_antibiotics(line="val3", panel_name="key1",
...                          panels=panels_dummy1)
>>> panels_dummy1["key1"][2]
'val3'

```

### _get_panel_name(panels, match_panel)

```python

>>> Isolate._get_panel_name(panels_dummy1, ":Panel: Campylobacter jejuni")
'campylobacter jejuni'
>>> panels_dummy1["campylobacter jejuni"]
[]

```
### _get_inclusions(panel_name, match_inclusion, inclusions)

```python

>>> inclusions = {}
>>> Isolate._get_inclusions(panel_name="campylobacter jejuni" ,
...                         line=":Include: Campylobacter",
...                         inclusions=inclusions)
True
>>> inclusions["campylobacter jejuni"]
['campylobacter']

>>> Isolate._get_inclusions(panel_name="campylobacter jejuni" ,
...                         line="colistin",
...                         inclusions=inclusions)
False

```
### _merge_inclusions(panels, inclusions)

```python
>>> panels_dummy2 = {"campylobacter": ["val2", "val3"],
...                  "campylobacter jejuni": ["val1"]}
>>> inclusions = {"campylobacter jejuni": ["campylobacter"]}
>>> Isolate._merge_inclusions(panels=panels_dummy2, inclusions=inclusions)
>>> panels_dummy2["campylobacter jejuni"]
['val1', 'val2', 'val3']

```
