# result_handler tests

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
...         self.acq_overlap = 1
...         self.min_cov = None
...         self.threshold = None
...         self.min_depth = None
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

First part just creates some dummy objects needed for testing the class. A phenoDB object, a hit from a blast output, ...

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

ResFinder stores results in a BlastNAligner object from cgecore, an iterator which will return each hit as a blastN hit object

Here a dummy hit (obtained using Blast) is created, named rf\_dat\_blast.

```python

>>> rf_dat_blast = {}
>>> rf_dat_blast["templateID"] = "gyrA_1_CP073768.1"
>>> rf_dat_blast["queryID"] = "contig_name"
>>> rf_dat_blast["gene_length_XMLFile_undescribed"] = 2628
>>> rf_dat_blast["score"] = 787.0
>>> rf_dat_blast["bitscore"] = 1454.43
>>> rf_dat_blast["evalue"] = 0.0
>>> rf_dat_blast["contig_name"] = "NA"
>>> rf_dat_blast["query_start"] = "NA"
>>> rf_dat_blast["query_end"] = "NA"
>>> rf_dat_blast["perc_coverage"] = 100
>>> rf_dat_blast["depth"] = 21
>>> rf_dat_blast["type"] = "aln_hit"
>>> rf_dat_blast["template_file"] = "blastn-alignment_subject-escherichia_coli"
>>> rf_dat_blast["templateID"] = "gyrA_1_CP073768.1"
>>> rf_dat_blast["queryID"] = "gyrA_test_lastPart"
>>> rf_dat_blast["gene_length_XMLFile_undescribed"] = 2628
>>> rf_dat_blast["score"] = 787.0
>>> rf_dat_blast["bitscore"] = 1454.43
>>> rf_dat_blast["evalue"] = 0.0
>>> rf_dat_blast["n_alignments"] = None
>>> rf_dat_blast["n_identity"] = 787
>>> rf_dat_blast["n_positives"] = 787
>>> rf_dat_blast["gaps"] = 0
>>> rf_dat_blast["aln_length"] = 787
>>> rf_dat_blast["strand"] = ('Plus', 'Plus')
>>> rf_dat_blast["frame"] = (1, 1)
>>> rf_dat_blast["query_aln"] = "GCGCGGTCGTCCGATCGTCAACCTGCTGCCGCTGGAGCAGGACGAACGTATCACTGCGATCCTGCCAGTGACCGAGTTTGAAGAAGGCGTGAAAGTCTTCATGGCGACCGCTAACGGTACCGTGAAGAAAACTGTCCTCACCGAGTTCAACCGTCTGCGTACCGCCGGTAAAGTGGCGATCAAACTGGTTGACGGCGATGAGCTGATCGGCGTTGACCTGACCAGCGGCGAAGACGAAGTAATGCTGTTCTCCGCTGAAGGTAAAGTGGTGCGCTTTAAAGAGTCTTCTGTCCGTGCGATGGGCTGCAACACCACCGGTGTTCGCGGTATTCGCTTAGGTGAAGGCGATAAAGTCGTCTCTCTGATCGTGCCTCGTGGCGATGGCGCAATCCTCACCGCAACGCAAAACGGTTACGGTAAACGTACCGCAGTGGCGGAATACCCAACCAAGTCGCGTGCGACGAAAGGGGTTATCTCCATCAAGGTTACCGAACGTAACGGTTTAGTTGTTGGCGCGGTACAGGTAGATGACTGCGACCAGATCATGATGATCACCGATGCCGGTACGCTGGTACGTACTCGCGTTTCGGAAATCAGCATCGTGGGCCGTAACACCCAGGGCGTGATCCTCATCCGTACTGCGGAAGATGAAAACGTAGTGGGTCTGCAACGTGTTGCTGAACCGGTTGACGAGGAAGATCTGGATACCATCGACGGCAGTGCCGCGGAAGGGGACGATGAAATCGCTCCGGAAGTGGACGTTGACGACGAGCCAGAAGAAGAATAA"
>>> rf_dat_blast["query_start_aln"] = 1
>>> rf_dat_blast["query_end_aln"] = 787
>>> rf_dat_blast["aln_scheme"] = "|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
>>> rf_dat_blast["template_aln"] = "GCGCGGTCGTCCGATCGTCAACCTGCTGCCGCTGGAGCAGGACGAACGTATCACTGCGATCCTGCCAGTGACCGAGTTTGAAGAAGGCGTGAAAGTCTTCATGGCGACCGCTAACGGTACCGTGAAGAAAACTGTCCTCACCGAGTTCAACCGTCTGCGTACCGCCGGTAAAGTGGCGATCAAACTGGTTGACGGCGATGAGCTGATCGGCGTTGACCTGACCAGCGGCGAAGACGAAGTAATGCTGTTCTCCGCTGAAGGTAAAGTGGTGCGCTTTAAAGAGTCTTCTGTCCGTGCGATGGGCTGCAACACCACCGGTGTTCGCGGTATTCGCTTAGGTGAAGGCGATAAAGTCGTCTCTCTGATCGTGCCTCGTGGCGATGGCGCAATCCTCACCGCAACGCAAAACGGTTACGGTAAACGTACCGCAGTGGCGGAATACCCAACCAAGTCGCGTGCGACGAAAGGGGTTATCTCCATCAAGGTTACCGAACGTAACGGTTTAGTTGTTGGCGCGGTACAGGTAGATGACTGCGACCAGATCATGATGATCACCGATGCCGGTACGCTGGTACGTACTCGCGTTTCGGAAATCAGCATCGTGGGCCGTAACACCCAGGGCGTGATCCTCATCCGTACTGCGGAAGATGAAAACGTAGTGGGTCTGCAACGTGTTGCTGAACCGGTTGACGAGGAAGATCTGGATACCATCGACGGCAGTGCCGCGGAAGGGGACGATGAAATCGCTCCGGAAGTGGACGTTGACGACGAGCCAGAAGAAGAATAA"
>>> rf_dat_blast["template_start_aln"] = 1842
>>> rf_dat_blast["template_end_aln"] = 2628

```

## Test ResultHandler

Make a ResultHandler object
```python
>>> from src.resfinder.cge.output.result_handler import ResultHandler
>>> res_handler = ResultHandler(conf, res_pheno_db)

```

### keep_hit(hit_dict, current_hit, current_key, key_list)
input
    hit_dict: dict containing all the previous selcted hits and their
        information
    current_hit: a dict containing the hit to be compared.
    current_key: key of the hit to be compared
    key_list: a list of all the keys in the hit_dict representing all
        found hits to be compared against.
output
    keys_to_keep: list of keys to keep. pairwise comparison keeping the
    best hit.
    keys_to_drop: list of keys to drop. another hit was found to be better
This function will check the current hit against all the other hits in
the hit_dict looking for overlaps. it will keep the overlapping hit with
the best identity or combined coverage/identity score.
Previously a part of Blaster.compare_results()

#### setup of variables

```python
>>> import collections

>>> hitA = {'tmpl_start': 16,
...			'tmpl_end': 23,
...			'query_start': 16,
...			'query_end': 23,
...			'query_string': 'GACGGCGATGAGCTGATCGGCGTTGACCTGAC',
...			'aln_string': '               ||||||||         ',
...			'tmpl_string': "ACTGCGATCCTGCCACGACGGCAGGTTATC",
...			'aln_length': 8,
...			'contig_name': 'test_other_contig',
...			'identity': 1.0,
...			'coverage': 0.2994672754946728,
...			'gene_length': 30,
...			'hit_id': 'test_other_contig:8..15:gyrA_1_CP073768.1'}

>>> hitB = {'tmpl_start': 1,
...			'tmpl_end': 10,
...			'query_start': 1,
...			'query_end': 10,
...			'query_string': 'GATTTCAACCGTCTGCGTACCGCCGGTAAA',
...			'aln_string': '|| |||||||                    ',
...			'tmpl_string': "GAGTTCAACCTGACGACGAGCCAGAAGAAA",
...			'aln_length': 10,
...			'contig_name': 'test_overlap',
...			'identity': 0.9091,
...			'coverage': 0.2994672754946728,
...			'gene_length': 30,
...			'hit_id': 'test_overlap:1..10:gyrA_1_CP073768.1'}

>>> hitC = {'tmpl_start': 15,
...			'tmpl_end': 25,
...			'query_start': 15,
...			'query_end': 25,
...			'query_string': 'GAGTTCAACCGTCTGCGTACCTCCGGTAAA',
...			'aln_string': '              ||||||| |||     ',
...			'tmpl_string': "TCCTCACCGCAACGGCGTACCGCCGACTCC",
...			'aln_length': 11,
...			'contig_name': 'test_overlap',
...			'identity': 0.9091,
...			'coverage': 0.2994672754946728,
...			'gene_length': 30,
...			'hit_id': 'test_overlap:15..25:gyrA_1_CP073768.1'}

>>> hitD = {'tmpl_start': 1,
...			'tmpl_end': 7,
...			'query_start': 1,
...			'query_end': 7,
...			'query_string': 'GAGTTCAACCGTCTGCGTACCGCCGGTAAA',
...			'aln_string': '|||||||                       ',
...			'tmpl_string': "AGATCTGGATGGCGCAATCCTCACCGCAAC",
...			'aln_length': 8,
...			'contig_name': 'test_overlap',
...			'identity': 1.0,
...			'coverage': 0.2994672754946728,
...			'gene_length': 30,
...			'hit_id': 'test_overlap:1..7:gyrA_1_CP073768.1'}

>>> hitE = {'tmpl_start': 7,
...			'tmpl_end': 17,
...			'query_start': 7,
...			'query_end': 17,
...			'query_string': 'GAGTTCAACCGTCTGCGTACCGCCGGTAAA',
...			'aln_string': '      |||||||||||             ',
...			'tmpl_string': "AGATCTAACCGTCTGCGTCCTCACCGCAAC",
...			'aln_length': 11,
...			'contig_name': 'test_overlap',
...			'identity': 1.0,
...			'coverage': 0.2994672754946728,
...			'gene_length': 30,
...			'hit_id': 'test_overlap:1..7:gyrA_1_CP073768.1'}

>>> hitF = {'tmpl_start': 7,
...			'tmpl_end': 17,
...			'query_start': 7,
...			'query_end': 17,
...			'query_string': 'GAGTTCAACCGTCAGCGTACCGCCGGTAAA',
...			'aln_string': '      ||||||| |||             ',
...			'tmpl_string': "AGATCTAACCGTCTGCGTCCTCACCGCAAC",
...			'aln_length': 11,
...			'contig_name': 'test_overlap',
...			'identity': 0.9091,
...			'coverage': 0.2994672754946728,
...			'gene_length': 30,
...			'hit_id': 'test_overlap:1..7:gyrA_1_CP073768.1'}

>>> hitG = {'tmpl_start': 7,
...			'tmpl_end': 17,
...			'query_start': 7,
...			'query_end': 17,
...			'query_string': 'GAGTTCAATCGTCAGCGTACCGCCGGTAAA',
...			'aln_string': '      || |||| |||             ',
...			'tmpl_string': "AGATCTAACCGTCTGCGTCCTCACCGCAAC",
...			'aln_length': 11,
...			'contig_name': 'test_overlap',
...			'identity': 0.8181818182,
...			'coverage': 0.2994672754946728,
...			'gene_length': 30,
...			'hit_id': 'test_overlap:1..7:gyrA_1_CP073768.1'}

>>> gene_dict_no_overlap = collections.defaultdict(list)
>>> gene_dict_overlap = collections.defaultdict(list)
>>> gene_exact_overlap = collections.defaultdict(list)

>>> gene_dict_no_overlap["gyrA_1_CP073768.1"].append(hitD)
>>> gene_dict_no_overlap["gyrB_1_CP047010.1"].append(hitA)

>>> gene_dict_overlap["gyrA_1_CP073768.1"].append(hitB)
>>> gene_dict_overlap["gyrB_1_CP047010.1"].append(hitC)

>>> gene_exact_overlap["gyrB_1_CP047010.1"].append(hitE)
>>> gene_exact_overlap["gyrA_1_CP073768.1"].append(hitF)
>>> gene_exact_overlap["parC_1_CP084529.1"].append(hitG)


>>> keys_nooverlap = ["gyrA_1_CP073768.1"]
>>> keys_overlap = ["gyrA_1_CP073768.1", "gyrB_1_CP047010.1"]
>>> keys_exact = ["parC_1_CP084529.1", "gyrA_1_CP073768.1", "gyrB_1_CP047010.1"]

>>> current_hit = {'tmpl_start': 7,
...			'tmpl_end': 17,
...			'query_start': 7,
...			'query_end': 17,
...			'query_string': 'AACTGTCTGCG',
...			'aln_string': '||| |||||||',
...			'tmpl_string': "AACCGTCTGCG",
...			'aln_length': 11,
...			'contig_name': 'test_overlap',
...			'identity': 0.9091,
...			'coverage': 0.2994672754946728,
...			'gene_length': 30,
...			'hit_id': 'test_overlap:7..17:rpoB_1_CP073768.1'}


```

#### Test
```python
#no overlap keeping both keys
>>> no_overlap = res_handler.keep_hit(gene_dict_no_overlap, current_hit, "rpoB_1_CP073768.1",
...                      keys_nooverlap)
>>> sorted(no_overlap[0])
['gyrA_1_CP073768.1', 'rpoB_1_CP073768.1']

>>> sorted(no_overlap[1])
[]

# overlap with same cal score keeping the longest alignment or both if no difference in aln_length
>>> equal_cal = res_handler.keep_hit(gene_dict_overlap, current_hit, "rpoB_1_CP073768.1",
...                      keys_overlap)
contig test_overlap was found to hit both 
<BLANKLINE>
rpoB_1_CP073768.1 and gyrA_1_CP073768.1
hit ['rpoB_1_CP073768.1'] was kept
contig test_overlap was found to hit both 
<BLANKLINE>
rpoB_1_CP073768.1 and gyrB_1_CP047010.1
hit ['gyrB_1_CP047010.1', 'rpoB_1_CP073768.1'] was kept

>>> sorted(equal_cal[0])
['gyrB_1_CP047010.1', 'rpoB_1_CP073768.1']

>>> sorted(equal_cal[1])
['gyrA_1_CP073768.1']

#exact overlap in contig - keeping the one with best identity or both if same id. 
>>> exact_overlap = res_handler.keep_hit(gene_exact_overlap, current_hit, "rpoB_1_CP073768.1",
...                      keys_exact)
contig test_overlap was found to hit both 
<BLANKLINE>
rpoB_1_CP073768.1 and parC_1_CP084529.1
hit ['rpoB_1_CP073768.1'] was kept
contig test_overlap was found to hit both 
<BLANKLINE>
rpoB_1_CP073768.1 and gyrA_1_CP073768.1
hit ['gyrA_1_CP073768.1', 'rpoB_1_CP073768.1'] was kept
contig test_overlap was found to hit both 
<BLANKLINE>
rpoB_1_CP073768.1 and gyrB_1_CP047010.1
hit ['gyrB_1_CP047010.1'] was kept

>>> sorted(exact_overlap[0])
['gyrA_1_CP073768.1', 'gyrB_1_CP047010.1']

>>> sorted(exact_overlap[1])
['parC_1_CP084529.1', 'rpoB_1_CP073768.1']

```

### get\_query\_align(hit, contig)
Input:
    hit - an instance of a gene_dict[gene][0] dict.
            contig - the input contig corresponding to the hit
output:
    updated query string
    updated alignment string
The function will complete the query and alingment in the case that the
template gene is longer than the hit sequence. If there is corresponding
sequences in the contig these will be added otherwise '-' will be
appended accordingly. The alignment string will include spaces to the
corresponding added sequence.

### setup

```python
>>> contig1 = "TTTTTTTTTTTTTTTCGTCAACCGGGGGGGGGGGGGGGGGGGGGG"

>>> hit1 = {'tmpl_start': 8,
...			'tmpl_end': 15,
...			'query_start': 16,
...			'query_end': 23,
...			'query_string': 'CGTCAACC',
...			'aln_string': '||||||||',
...			'tmpl_string': "CGTCAACC",
...			'aln_length': 8,
...			'contig_name': 'gyrA_test_lastPart',
...			'identity': 1.0,
...			'coverage': 0.2994672754946728,
...			'gene_length': 30,
...			'hit_id': 'gyrA_test_lastPart:10..20:gyrA_1_CP073768.1'}

>>> contig2 = "TTTTCGTCAACCGGGGGGG"

>>> hit2 = {'tmpl_start': 8,
...			'tmpl_end': 15,
...			'query_start': 5,
...			'query_end': 12,
...			'query_string': 'CGTCAACC',
...			'aln_string': '||||||||',
...			'tmpl_string': "CGTCAACC",
...			'aln_length': 8,
...			'contig_name': 'gyrA_test_lastPart',
...			'identity': 1.0,
...			'coverage': 0.2994672754946728,
...			'gene_length': 30,
...			'hit_id': 'gyrA_test_lastPart:10..20:gyrA_1_CP073768.1'}

>>> contig3 = "CGTCAACC"

>>> hit3 = {'tmpl_start': 8,
...			'tmpl_end': 15,
...			'query_start': 1,
...			'query_end': 8,
...			'query_string': 'CGTCAACC',
...			'aln_string': '||||||||',
...			'tmpl_string': "CGTCAACC",
...			'aln_length': 8,
...			'contig_name': 'gyrA_test_lastPart',
...			'identity': 1.0,
...			'coverage': 0.2994672754946728,
...			'gene_length': 30,
...			'hit_id': 'gyrA_test_lastPart:10..20:gyrA_1_CP073768.1'}

```

### Test
```python
>>> res_handler.get_query_align(hit=hit1, contig=contig1)
'TTTTTTTCGTCAACCGGGGGGGGGGGGGGG'

>>> res_handler.get_query_align(hit=hit2, contig=contig2)
'---TTTTCGTCAACCGGGGGGG--------'

>>> res_handler.get_query_align(hit=hit3, contig=contig3)
'-------CGTCAACC---------------'

```

### calculate_alignment(query, template)
input:
    query sequence - string
    template sequence - string
output:
    alignment sequence - string
This function will calculate the alignment. | for match and space for
mismatch. This might also be done somewhere else in the code - could not find it

```python
>>> query = '-------CGTCAACC---------------'
>>> template = 'GTCCGATCGTCAACCTGCTGCCGCTGGAGC'
>>> res_handler.calculate_alignment(query, template)
'       ||||||||               '


```