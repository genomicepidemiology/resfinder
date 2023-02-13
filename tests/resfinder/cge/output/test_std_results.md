# std_results tests

## setup

```python

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

## initialize

First part just creates some dummy objects needed for testing the class. A
Result object containing the databases are created, a list of GeneResult objects
(actually only one object), and a list containing results from a PointFinder
hit.

```python

>>> from cgelib.output.result import Result
>>> from src.resfinder.cge.output.gene_result import GeneResult
>>> from cgecore.blaster.blaster import Blaster
>>> import os

>>> res = Result.init_software_result(name="ResFinder", gitdir=".")
>>> res.init_database("ResFinder", ".")
>>> res.init_database("PointFinder", ".")

```

ResFinder stores results in a dict of dicts. The internal dicts represents hits
to the sub-database the dict represent. The sub-database is for legacy reasons
named after amr classes and to some extend represents hits to genes in that amr
class. This is not a naming convention that will be guaranteed in future
releases as the structure is not optimal.

Here a dummy hit (obtained using KMA) is created, named rf\_dat\_kma. It is stored in the beta-lactam sub-database/dict in the ResFinder result dict named
rf\_custom\_kma.

```python
import tempfile

>>> resfinder_db_path = os.environ["CGE_RESFINDER_RESGENE_DB"]
>>> blast_path = os.environ["CGE_BLASTN"]
>>> rf_custom_kma = Blaster(inputfile="tests/data/test_isolate_01.fa", databases=["beta-lactam"], db_path=resfinder_db_path, out_path="tests/tmp_out", blast=blast_path)
... #doctest: +ELLIPSIS
Found...

```

The PointFinder part of ResFinder unfortunately stores its results differently than ResFinder. The results are stored in a dict (here stored in the variable 'pf\_custom\_blast'), which contains a dict entry for each gene searched for mutations in the given species (corresponding dummy below is named 'gyrA'), but only if the gene is found. If a gene is not found its not a dict but the string "No hit found". A gene dict then has an entry for each hit to the reference gene in the database. The hit entry itself is also a dict (dummy hit below is named gyrA_hit).

The actual point mutations are stored in the gene dict (dummy: gyrA) in a list
of lists, where each list in the list describes a mutation.

```python

>>> gyrA_hit = {}
>>> gyrA_hit["evalue"] = 0.0
>>> gyrA_hit["sbjct_header"] = "gyrA_1_CP073768.1"
>>> gyrA_hit["bit"] = 4843.04
>>> gyrA_hit["perc_ident"] = 99.92
>>> gyrA_hit["sbjct_length"] = 2628
>>> gyrA_hit["sbjct_start"] = 1
>>> gyrA_hit["sbjct_end"] = 2628
>>> gyrA_hit["gaps"] = 0
>>> gyrA_hit["contig_name"] = "gyrA_G81D_GAT_D82G_GGC"
>>> gyrA_hit["query_start"] = 1
>>> gyrA_hit["query_end"] = 2628
>>> gyrA_hit["HSP_length"] = 2628
>>> gyrA_hit["coverage"] = 1.0
>>> gyrA_hit["cal_score"] = 99.92
>>> gyrA_hit["hit_id"] = "gyrA_G81D_GAT_D82G_GGC:1..2628:gyrA:99.923896"
>>> gyrA_hit["strand"] = 0
>>> gyrA_hit["perc_coverage"] = 100.0

>>> gyrA = {}
>>> gyrA["mis_matches"] = [
...   [ 'sub', 81, 81, 'D', 'p.G81D', 'GGT', 'GAT', 'G', 'D' ],
...   [ 'sub', 82, 82, 'G', 'p.D82G', 'GAC', 'GGC', 'D', 'G' ] ]
>>> gyrA["hits"] = {}
>>> gyrA["hits"]["gyrA_G81D_GAT_D82G_GGC:1..2628:gyrA:99.923896"] = gyrA_hit

>>> pf_custom_blast = {}
>>> pf_custom_blast["excluded"] = {}
>>> pf_custom_blast["gyrA_1_CP073768.1"] = gyrA
>>> pf_custom_blast["gyrB"] = "No hit found"

```

Create dummy data for KMA results, which are handled slightly differently for
PointFinder.

```python

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
>>> pf_custom_kma["gyrA_1_CP073768.1"] = gyrA_kma_hit
>>> pf_custom_kma["gyrB"] = "No hit found"

>>> import copy
>>> res_kma_test = copy.deepcopy(res)

```

Create the phenoDB object.

```python

>>> import os
>>> from src.resfinder.cge.phenotype2genotype.res_profile import PhenoDB

>>> resfinder_db_path = os.environ["CGE_RESFINDER_RESGENE_DB"]
>>> assert(len(resfinder_db_path) > 0)
>>> pointfinder_db_path = os.environ["CGE_RESFINDER_RESPOINT_DB"]
>>> assert(len(pointfinder_db_path) > 0)

>>> abclassdef_file= "{}/antibiotic_classes.txt".format(resfinder_db_path)
>>> acquired_file= "{}/phenotypes.txt".format(resfinder_db_path)
>>> point_file = ("{}/escherichia_coli/phenotypes.txt"
...               .format(pointfinder_db_path))
>>> res_pheno_db = PhenoDB(abclassdef_file=abclassdef_file,
...                        acquired_file=acquired_file,
...                        point_file=point_file)

```

## Test ResFinderResultHandler

### standardize\_results(res\_collection, res, ref\_db\_name)

```python

>>> from src.resfinder.cge.output.std_results import ResFinderResultHandler
>>> ResFinderResultHandler.standardize_results(res,
...                                            rf_custom_kma,
...                                            "ResFinder",
...                                            conf)

>>> for k in res["databases"]:
...   print(k)
... #doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
ResFinder-...
PointFinder-...

>>> for k in res["seq_regions"]:
...   print(k)
blaB-2;;1;;AF189300

```

### load\_res\_profile(res\_collection, isolate)

#### Setup

Create Isolate object

```python

>>> from src.resfinder.cge.phenotype2genotype.isolate import Isolate
>>> isolate = Isolate(name="Test sample")

>>> isolate.load_finder_results(std_table=res, phenodb=res_pheno_db,
...                             type="seq_regions")
>>> isolate.load_finder_results(std_table=res, phenodb=res_pheno_db,
...                             type="seq_variations")
>>> isolate.calc_res_profile(res_pheno_db)

>>> import os.path
>>> import inspect
>>> from cgelib.utils.loaders_mixin import LoadersMixin
>>> std_results_file = inspect.getfile(ResFinderResultHandler)
>>> std_results_dir = os.path.dirname(os.path.realpath(std_results_file))
>>> amr_abbreviations_file = ("{}/../../amr_abbreviations.md"
...                           .format(std_results_dir))
>>> amr_abbreviations = LoadersMixin.load_md_table_after_keyword(
...     amr_abbreviations_file, "## Abbreviations")

```

#### Test

```python

>>> from src.resfinder.cge.output.std_results import PointFinderResultHandler
>>> ResFinderResultHandler.load_res_profile(res, isolate, amr_abbreviations)
>>> res["phenotypes"]["ampicillin"]["seq_regions"]
['blaB-2;;1;;AF189300']
>>> res["result_summary"]
'AMO_AMP_TIC'

```

## Test PointFinderResultHandler

### standardize\_results(res\_collection, res, ref\_db\_name)

```python

>>> PointFinderResultHandler.standardize_results(res,
...                                              pf_custom_blast,
...                                              "PointFinder")

>>> for k in res["databases"]:
...   print(k)
... #doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
ResFinder-...
PointFinder-...

>>> for k in res["seq_regions"]:
...   print(k)
blaB-2;;1;;AF189300
gyrA;;1;;CP073768.1

>>> for k in res["seq_variations"]:
...   print(k)
gyrA;;1;;CP073768.1;;81;;d
gyrA;;1;;CP073768.1;;82;;g

```

```python

>>> from src.resfinder.cge.output.std_results import PointFinderResultHandler
>>> PointFinderResultHandler.standardize_results(res_kma_test,
...                                              pf_custom_kma,
...                                              "PointFinder")

>>> for k in res_kma_test["databases"]:
...   print(k)
... #doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
ResFinder-...
PointFinder-...

>>> for k in res_kma_test["seq_regions"]:
...   print(k)
gyrA;;1;;CP073768.1

>>> for k in res_kma_test["seq_variations"]:
...   print(k)
gyrA;;1;;CP073768.1;;81;;d
gyrA;;1;;CP073768.1;;82;;g

```

### create\_amr\_summary\_str(res_collection, amr_abbreviations)

#### Setup

```python

>>> phenotype1 = {
...     "type": "phenotype",
...     "key": "ampicillin+clavulanic acid",
...     "category": "amr",
...     "amr_classes": ['beta-lactam'],
...     "amr_resistance": "ampicillin+clavulanic acid",
...     "amr_resistant": True,
...     "amr_species_relevant": True}
>>> phenotype2 = {
...     "type": "phenotype",
...     "key": "imipenem",
...     "category": "amr",
...     "amr_classes": ['beta-lactam'],
...     "amr_resistance": "imipenem",
...     "amr_resistant": True,
...     "amr_species_relevant": True}
>>> phenotype3 = {
...     "type": "phenotype",
...     "key": "doxycycline",
...     "category": "amr",
...     "amr_classes": ['tetracycline'],
...     "amr_resistance": "doxycycline",
...     "amr_resistant": True,
...     "amr_species_relevant": False}
>>> res.add_class(cl="phenotypes", clobber_warn=False, **phenotype1)
>>> res.add_class(cl="phenotypes", clobber_warn=False, **phenotype2)
>>> res.add_class(cl="phenotypes", clobber_warn=False, **phenotype3)

```

```python

>>> sum_str = ResFinderResultHandler.create_amr_summary_str(res, amr_abbreviations)
>>> sum_str
'AMO_AMP_AML_IMI_TIC'

```
