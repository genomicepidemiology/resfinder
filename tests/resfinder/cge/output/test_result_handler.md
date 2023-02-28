# result_handler tests

## setup

```python

>>> from src.resfinder.cge.config import Config

>>> class DummyArgs():
...     def __init__(self):
...         self.inputfasta = "./tests/data/test_isolate_12.fa"
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
...         self.disinfectant = True
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
>>> disinfinder_db_path = os.environ["CGE_DISINFINDER_DB"]
>>> assert(len(disinfinder_db_path) > 0)

>>> abclassdef_file= "{}/antibiotic_classes.txt".format(resfinder_db_path)
>>> acquired_file= "{}/phenotypes.txt".format(resfinder_db_path)
>>> point_file = ("{}/escherichia_coli/phenotypes.txt"
...               .format(pointfinder_db_path))
>>> res_pheno_db = PhenoDB(abclassdef_file=abclassdef_file,
...                        acquired_file=acquired_file,
...                        point_file=point_file)

```

ResFinder stores results in a BlastNAligner object from cgecore, an iterator which will return each hit as a blastN hit object

Here a dummy hit (obtained using Blast) is created, named rf_dat_blast.

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

# different hits as found in gene_dict
>>> hitA = {'tmpl_start': 1,
...			'tmpl_end': 153,
...			'query_start': 1,
...			'query_end': 153,
...			'query_string': 'GTGTCCACACCACATCACGGCCGGCACGAGCTCGGCCAGAACTTCCTGTCCGATCGGCGCGTCATCGCCGATATCGTCGAAATCGTCTCGCGCACAAACGGTCCGATCATCGAGATCGGGGCGGGCGACGGCGCGCTGACCATACCCTTGCAA',
...			'aln_string': '|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||',
...			'tmpl_string': 'GTGTCCACACCACATCACGGCCGGCACGAGCTCGGCCAGAACTTCCTGTCCGATCGGCGCGTCATCGCCGATATCGTCGAAATCGTCTCGCGCACAAACGGTCCGATCATCGAGATCGGGGCGGGCGACGGCGCGCTGACCATACCCTTGCAA',
...			'aln_length': 153,
...			'contig_name': 'erm(38)_first1-321',
...			'identity': 1.0,
...			'coverage': 0.13178294573643412,
...			'gene_length': 1161,
...			'hit_id': 'erm(38)_first1-321:8..15:erm(38)_1_AY154657.1'}

>>> hitA2 = {'tmpl_start': 1,
...			'tmpl_end': 239,
...			'query_start': 1,
...			'query_end': 239,
...			'query_string': 'AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGAAGCTTGCTTCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCGGATGTGCCCAGAT',
...			'aln_string': '|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||',
...			'tmpl_string': 'AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGAAGCTTGCTTCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCGGATGTGCCCAGAT',
...			'aln_length': 239,
...			'contig_name': '16S-rrsB_1_CP067250.1',
...			'identity': 1.0,
...			'coverage': 0.1549935149,
...			'gene_length': 1542,
...			'hit_id': 'test_other_contig:8..15:16S-rrsB_1_CP067250.1'}

>>> hitA3 = {'tmpl_start': 61,
...			'tmpl_end': 321,
...			'query_start': 61,
...			'query_end': 321,
...			'query_string': 'GTCATCGCCGATATCGTCGAAATCGTCTCGCGCACAAACGGTCCGATCATCGAGATCGGGGCGGGCGACGGCGCGCTGACCATACCCTTGCAACGACTCGCCCGCCCGCTCACCGCCGTCGAGGTCGACGCGCGGCGCGCGCGGCGGTTGGCGCAGCGCACCGCGAGATCCGCCCCGGGGCCTGCCTCGCGGCCCACCGAGGTCGTCGCCGCCGACTTCCTGCGCTACCCACTGCCCCGCTCACCCCACGTGGTCGTGGGC',
...			'aln_string': '|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||',
...			'tmpl_string': 'GTCATCGCCGATATCGTCGAAATCGTCTCGCGCACAAACGGTCCGATCATCGAGATCGGGGCGGGCGACGGCGCGCTGACCATACCCTTGCAACGACTCGCCCGCCCGCTCACCGCCGTCGAGGTCGACGCGCGGCGCGCGCGGCGGTTGGCGCAGCGCACCGCGAGATCCGCCCCGGGGCCTGCCTCGCGGCCCACCGAGGTCGTCGCCGCCGACTTCCTGCGCTACCCACTGCCCCGCTCACCCCACGTGGTCGTGGGC',
...			'aln_length': 261,
...			'contig_name': 'erm(38)_first1-321',
...			'identity': 1.0,
...			'coverage': 0.2248062015503876,
...			'gene_length': 1161,
...			'hit_id': 'erm(38)_first1-321:8..15:erm(38)_1_AY154657.1'}

>>> hitA4 = {'tmpl_start': 202,
...			'tmpl_end': 501,
...			'query_start': 1,
...			'query_end': 300,
...			'query_string': '--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------CGGCGGTTGGCGCAGCGCACCGCGAGATCCGCCCCGGGGCCTGCCTCGCGGCCCACCGAGGTCGTCGCCGCCGACTTCCTGCGCTACCCACTGCCCCGCTCACCCCACGTGGTCGTGGGCAACCTGCCGTTCCACCTCACCACCGCGATCCTGCGGCGACTGCTGCACGGTCCGGGCTGGACCACGGCCGTGCTGCTCATGCAGTGGGAGGTGGCCCGCCGACGCGCCGCGGTGGGCGGCGCCACCATGATGACCGCCCAGTGGTGGCCGTGGTTCGAATTCGGCCTTGCCCGAAAGGTTTCCGCGGCGAGCTTCACGCCGCGGCCCGCGGTCGACGCCGGACTGCTCACCATCACGCGCCGCAGCCGGCCGCTGGTCGACGTCGCGGACCGGGCGCGT----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------',
...			'aln_string': '                                                                                                                                                                                                        ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     ',
...			'tmpl_string': 'GTGTCCACACCACATCACGGCCGGCACGAGCTCGGCCAGAACTTCCTGTCCGATCGGCGCGTCATCGCCGATATCGTCGAAATCGTCTCGCGCACAAACGGTCCGATCATCGAGATCGGGGCGGGCGACGGCGCGCTGACCATACCCTTGCAACGACTCGCCCGCCCGCTCACCGCCGTCGAGGTCGACGCGCGGCGCGCGCGGCGGTTGGCGCAGCGCACCGCGAGATCCGCCCCGGGGCCTGCCTCGCGGCCCACCGAGGTCGTCGCCGCCGACTTCCTGCGCTACCCACTGCCCCGCTCACCCCACGTGGTCGTGGGCAACCTGCCGTTCCACCTCACCACCGCGATCCTGCGGCGACTGCTGCACGGTCCGGGCTGGACCACGGCCGTGCTGCTCATGCAGTGGGAGGTGGCCCGCCGACGCGCCGCGGTGGGCGGCGCCACCATGATGACCGCCCAGTGGTGGCCGTGGTTCGAATTCGGCCTTGCCCGAAAGGTTTCCGCGGCGAGCTTCACGCCGCGGCCCGCGGTCGACGCCGGACTGCTCACCATCACGCGCCGCAGCCGGCCGCTGGTCGACGTCGCGGACCGGGCGCGTTACCAGGCGCTGGTGCACCGCGTGTTCACCGGACGCGGACACGGCATGGCGCAGATCCTGCAACGGTTGCCCACGCCGGTGCCCCGCACTTGGTTGCGGGCCAACGGGATAGCACCGAACTCCCTGCCCCGCCAGTTGTCCGCGGCGCAGTGGGCGGCGCTGTTCGAGCAGACGCGTCTAACTGGTGCCCAACGGGTCGATCGTCCACGCGATGTACAGCACGGCCGCGCTCACCGTCGCCGTGGTGGCGAAGTCGATCGCCCGGCTACGCACCACAAGCAGACCGGCCCGGTCGTCGGTCAGCGCCAACCGCAGCGCGGCCGCGACGCCGACGCCGATCCCGATGACCAGCGCACCGCGCCGCCAGTAACCCGCCACCACCAGGGCGAACGCCGCGATGAAGATCAGGCCGACCACCAGGATCGGCCATTGACCGGCGAACACCTTGCGGGCGAATTCCTTTGGCGTCACGCCAGTTTCGACTCTTCGGCTTCGACGACGTTGGTCAGCAGGAAGGCGCGGGTCAACGGGCCCACGCCACCGGGGTTGGGCGACACGTGA',
...			'aln_length': 300,
...			'contig_name': 'erm(38)_middel201-600',
...			'identity': 1.0,
...			'coverage': 0.2248062015503876,
...			'gene_length': 1161,
...			'hit_id': 'erm(38)_middel201-600:8..15:erm(38)_1_AY154657.1'}



>>> hitB = {'tmpl_start': 1,
...			'tmpl_end': 10,
...			'query_start': 1,
...			'query_end': 10,
...			'query_string': 'GATTTCAACCGTCTGCGTACCGCCGGTAAA',
...			'aln_string': '|| |||||||                    ',
...			'tmpl_string': 'GAGTTCAACCTGACGACGAGCCAGAAGAAA',
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
...			'tmpl_string': 'TCCTCACCGCAACGGCGTACCGCCGACTCC',
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
...			'tmpl_string': 'AGATCTGGATGGCGCAATCCTCACCGCAAC',
...			'aln_length': 8,
...			'contig_name': 'test_overlap',
...			'identity': 1.0,
...			'coverage': 0.2994672754946728,
...			'gene_length': 30,
...			'hit_id': 'test_overlap:1..7:gyrA_1_CP073768.1'}

>>> hitE = {'tmpl_start': 1,
...			'tmpl_end': 30,
...			'query_start': 1,
...			'query_end': 30,
...			'query_string': 'GAGTTCAACCGTCTGCGTACCGCCGGTAAA',
...			'aln_string': '      |||||||||||             ',
...			'tmpl_string': 'AGATCTAACCGTCTGCGTCCTCACCGCAAC',
...			'aln_length': 11,
...			'contig_name': 'test_overlap',
...			'identity': 1.0,
...			'coverage': 0.2994672754946728,
...			'gene_length': 30,
...			'hit_id': 'test_overlap:1..7:gyrA_1_CP073768.1'}

>>> hitF = {'tmpl_start': 1,
...			'tmpl_end': 30,
...			'query_start': 1,
...			'query_end': 30,
...			'query_string': 'GAGTTCAACCGTCAGCGTACCGCCGGTAAA',
...			'aln_string': '      ||||||| |||             ',
...			'tmpl_string': 'AGATCTAACCGTCTGCGTCCTCACCGCAAC',
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
...			'query_string': 'GAGCTCAATCGTCAGCGGACCGCCGGTAAA',
...			'aln_string': '      || |||| |||             ',
...			'tmpl_string': 'AGATCTAACCGTCTGCGTCCTCACCGCAAC',
...			'aln_length': 11,
...			'contig_name': 'test_overlap',
...			'identity': 0.8181818182,
...			'coverage': 0.2994672754946728,
...			'gene_length': 30,
...			'hit_id': 'test_overlap:1..7:gyrA_1_CP073768.1'}

>>> hitH = {'tmpl_start': 16,
...			'tmpl_end': 25,
...			'query_start': 16,
...			'query_end': 25,
...			'query_string': 'CGTCCTCACC',
...			'aln_string': '||||||||||',
...			'tmpl_string': "CGTCCTCACC",
...			'aln_length': 10,
...			'contig_name': 'test_overlapH',
...			'identity': 1.0,
...			'coverage': 0.3333333333,
...			'gene_length': 30,
...			'hit_id': 'test_overlap:1..7:gyrA_1_CP073768.1'}

>>> hitI = {'tmpl_start': 9,
...			'tmpl_end': 30,
...			'query_start': 9,
...			'query_end': 30,
...			'query_string': 'GAGCTCGTCCGTCTGCGTCCTCACCTTTTT',
...			'aln_string': '        ||||||||||||||||||||||',
...			'tmpl_string': 'AGATCTAACCGTCTGCGTCCTCACCTTTTT',
...			'aln_length': 22,
...			'contig_name': 'test_overlap',
...			'identity': 1.0,
...			'coverage': 0.7333333333333333,
...			'gene_length': 30,
...			'hit_id': 'test_overlap:1..7:gyrA_1_CP073768.1'}

>>> hitJ = {'tmpl_start': 20,
...			'tmpl_end': 30,
...			'query_start': 10,
...			'query_end': 20,
...			'query_string': 'CTCACCGCAAC',
...			'aln_string': '|||||||||||',
...			'tmpl_string': 'CTCACCGCAAC',
...			'aln_length': 11,
...			'contig_name': 'test_overlapJ',
...			'identity': 1,
...			'coverage': 0.36666666666666664,
...			'gene_length': 30,
...			'hit_id': 'test_overlap:1..7:gyrA_1_CP073768.1'}

>>> hitI_N = {'tmpl_start': 9,
...			'tmpl_end': 30,
...			'query_start': 9,
...			'query_end': 30,
...			'query_string': 'GAGCTCGTCCGTCTGCGTCCNNNCCTTTTT',
...			'aln_string': '        ||||||||||||   |||||||',
...			'tmpl_string': 'AGATCTAACCGTCTGCGTCCTCACCTTTTT',
...			'aln_length': 22,
...			'contig_name': 'test_overlapI',
...			'identity': 0.8636363636363636,
...			'coverage': 0.6333333333333333,
...			'gene_length': 30,
...			'hit_id': 'test_overlap:1..7:gyrA_1_CP073768.1'}

>>> hitK = {'tmpl_start': 141,
...         'tmpl_end': 321,
...         'query_start': 141,
...         'query_end': 321,
...         'query_string': 'ATACCCTTGCAACGACTCGCCCGCCCGCTCACCGCCGTCGAGGTCGACGCGCGGCGCGCGCGGCGGTTGGCGCAGCGCACCGCGAGATCCGCCCCGGGGCCTGCCTCGCGGCCCACCGAGGTCGTCGCCGCCGACTTCCTGCGCTACCCACTGCCCCGCTCACCCCACGTGGTCGTGGGC',
...         'aln_string': '||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||',
...         'tmpl_string': 'ATACCCTTGCAACGACTCGCCCGCCCGCTCACCGCCGTCGAGGTCGACGCGCGGCGCGCGCGGCGGTTGGCGCAGCGCACCGCGAGATCCGCCCCGGGGCCTGCCTCGCGGCCCACCGAGGTCGTCGCCGCCGACTTCCTGCGCTACCCACTGCCCCGCTCACCCCACGTGGTCGTGGGC',
...         'aln_length': 180,
...         'contig_name': 'erm(38)_first1-321',
...         'identity': 1.0,
...         'coverage': 0.15503875968992248,
...         'gene_length': 1161,
...         'hit_id': 'erm(38)_first1-321:1..7:erm(38)_1_AY154657'}

>>> hitL = {'tmpl_start': 201,
...         'tmpl_end': 341,
...         'query_start': 1,
...         'query_end': 140,
...         'query_string': '---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------CGGCGGTTGGCGCAGCGCACCGCGAGATCCGCCCCGGGGCCTGCCTCGCGGCCCACCGAGGTCGTCGCCGCCGACTTCCTGCGCTACCCACTGCCCCGCTCACCCCACGTGGTCGTGGGCAACCTGCCGTTCCACCTCACCACCGCGATCCTGCGGCGACTGCTGCACGGTCCGGGCTGGACCACGGCCGTGCTGCTCATGCAGTGGGAGGTGGCCCGCCGACGCGCCGCGGTGGGCGGCGCCACCATGATGACCGCCCAGTGGTGGCCGTGGTTCGAATTCGGCCTTGCCCGAAAGGTTTCCGCGGCGAGCTTCACGCCGCGGCCCGCGGTCGACGCCGGACTGCTCACCATCACGCGCCGCAGCCGGCCGCTGGTCGACGTCGCGGACCGGGCGCGT---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------',
...         'aln_string': '                                                                                                                                                                                                         ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    ',
...         'tmpl_string': 'GTGTCCACACCACATCACGGCCGGCACGAGCTCGGCCAGAACTTCCTGTCCGATCGGCGCGTCATCGCCGATATCGTCGAAATCGTCTCGCGCACAAACGGTCCGATCATCGAGATCGGGGCGGGCGACGGCGCGCTGACCATACCCTTGCAACGACTCGCCCGCCCGCTCACCGCCGTCGAGGTCGACGCGCGGCGCGCGCGGCGGTTGGCGCAGCGCACCGCGAGATCCGCCCCGGGGCCTGCCTCGCGGCCCACCGAGGTCGTCGCCGCCGACTTCCTGCGCTACCCACTGCCCCGCTCACCCCACGTGGTCGTGGGCAACCTGCCGTTCCACCTCACCACCGCGATCCTGCGGCGACTGCTGCACGGTCCGGGCTGGACCACGGCCGTGCTGCTCATGCAGTGGGAGGTGGCCCGCCGACGCGCCGCGGTGGGCGGCGCCACCATGATGACCGCCCAGTGGTGGCCGTGGTTCGAATTCGGCCTTGCCCGAAAGGTTTCCGCGGCGAGCTTCACGCCGCGGCCCGCGGTCGACGCCGGACTGCTCACCATCACGCGCCGCAGCCGGCCGCTGGTCGACGTCGCGGACCGGGCGCGTTACCAGGCGCTGGTGCACCGCGTGTTCACCGGACGCGGACACGGCATGGCGCAGATCCTGCAACGGTTGCCCACGCCGGTGCCCCGCACTTGGTTGCGGGCCAACGGGATAGCACCGAACTCCCTGCCCCGCCAGTTGTCCGCGGCGCAGTGGGCGGCGCTGTTCGAGCAGACGCGTCTAACTGGTGCCCAACGGGTCGATCGTCCACGCGATGTACAGCACGGCCGCGCTCACCGTCGCCGTGGTGGCGAAGTCGATCGCCCGGCTACGCACCACAAGCAGACCGGCCCGGTCGTCGGTCAGCGCCAACCGCAGCGCGGCCGCGACGCCGACGCCGATCCCGATGACCAGCGCACCGCGCCGCCAGTAACCCGCCACCACCAGGGCGAACGCCGCGATGAAGATCAGGCCGACCACCAGGATCGGCCATTGACCGGCGAACACCTTGCGGGCGAATTCCTTTGGCGTCACGCCAGTTTCGACTCTTCGGCTTCGACGACGTTGGTCAGCAGGAAGGCGCGGGTCAACGGGCCCACGCCACCGGGGTTGGGCGACACGTGA',
...         'aln_length': 140,
...         'contig_name': 'erm(38)_middel201-578',
...         'identity': 1,
...         'coverage': 0.12058570198105081,
...         'gene_length': 1161,
...         'hit_id': 'erm(38)_middel201-578:1..7:erm(38)_1_AY154657'}




```

Make a ResultHandler object
```python
>>> from src.resfinder.cge.output.result_handler import ResultHandler
>>> res_handler = ResultHandler(conf, res_pheno_db)

```

## Static methods

### calculate_alignment(query, template)
input: <br>
    query sequence - string <br>
template sequence - string <br>
output: <br>
    alignment sequence - string <br>
This function will calculate the alignment. | for match and space for
mismatch. This might also be done somewhere else in the code - could not find it

```python
>>> query = '-------CGTCAACC---------------'
>>> template = 'GTCCGATCGTCAACCTGCTGCCGCTGGAGC'
>>> res_handler.calculate_alignment(query, template)
'       ||||||||               '


```
### calculate_identity(query, subject) 
Input: <br>
    query - query sequence string <br>
    subject - subject sequence string <br>
output: <br>
the identity as int.<br> 
The function will compare the two sequences on each position and return
the identity as the relationship between shared base pairs over all bases 
The two strings must be of same length.

```python
>>> str_1 = "AAATTTTGGGGTT"
>>> str_2 = "AAATTTTGGGGTT"
>>> str_3 = "AAA--TTGGGGTT"
>>> res_handler.calculate_identity(str_1, str_2)
1.0
>>> res_handler.calculate_identity(str_2, str_3)
0.8461538461538461

```

## Private methods

### _find_overlap_start(pre_start, pre_tmpl, next_start)
Input: <br>
pre_start - int of the start position of the first hit in regards 
            to the template <br>
pre_tmpl - template sequence string <br>
next_start - int of start position of the next hit in regards to the
             template <br>
Output: <br> 
start position of the overlap - int <br>
This function will iterate through the template until the overlap_start 
        position is reached. 

```python
>>> templ = 'AAAAAAAAATTTTTTTTTTTTTTT'
>>> templ2 = 'AA---AA--TTTTTTTTTTT'
>>> templ3 = 'AGATCTAACCGTCTGCGTCCTCACCGCAAC'

>>> res_handler._find_overlap_start(1, templ, 6)
6

>>> res_handler._find_overlap_start(1, templ2, 6)
11

>>> res_handler._find_overlap_start(7, templ3, 16)
16

```

### _get_overlap_seq(overlap_start, qry_overlap_start, pre_query_start, overlap_len, pre_qry, next_qry, template)
 Input:<br> 
 overlap start - int <br>
 pre_query_start - int <br>
 overlap_len - int <br>
 pre_qry - string of query sequence <br>
 next_qry - string of next query sequence <br>
 template - string of database template <br> 
 output: <br> 
 pre_qry_overlap - string of overlap found in the query sequence<br> 
 next_qry_overlap - string of overlap foind in next sequence <br>
 tmpl_overlap   - string of overlap matching the template. <br> 
 The function will find the overlap sequences in the previous query, 
 next query and the template. 

```python
>>> pre_qry = "AAAATTCTTGGGGGT"
>>> pre_qry1 = "AATTCTTGGGGG"
>>> next_qry1 = "TTTTTGGTGGCCCCC"
>>> next_qry2 = "GGTGG"
>>> template = "AAAATTTTTGGGGGCCCCC"

>>> test_qry = "GAGCTCAATCGTCAGCGTACCGCCGGTAAA"
>>> test_next = "CGTCCTCACC"
>>> test_tmpl = "AGATCTAACCGTCTGCGTCCTCACCGCAAC"

>>> res_handler._get_overlap_seq(5, 5, 1, 10, pre_qry, next_qry1, template)
('TTCTTGGGGG', 'TTTTTGGTGG', 'TTTTTGGGGG')

>>> res_handler._get_overlap_seq(10, 10, 1, 5, pre_qry, next_qry2, template)
('GGGGG', 'GGTGG', 'GGGGG')

>>> res_handler._get_overlap_seq(10, 10, 1, 1, pre_qry, next_qry2, template)
('G', 'G', 'G')

>>> res_handler._get_overlap_seq(5, 3, 1, 10, pre_qry1, next_qry1, template)
('TTCTTGGGGG', 'TTTTTGGTGG', 'TTTTTGGGGG')

>>> res_handler._get_overlap_seq(16, 16, 7, 2, test_qry, test_next, test_tmpl)
('CG', 'CG', 'CG')

```


### _get_best_overlap(pre_seq, next_seq, tmpl)
 Input:<br>
 pre_seq - sequence string <br>
 next_seq - second sequence string <br>
 tmpl    - template sequence string <br>
 Output: <br>
 sequence string <br>
aln_seq - corresponding alignemnt of the best overlap sequence <br>
 The function returns the string with the highest identity compared to
 the template sequence. All input sequences must be of same length

```python
>>> seq_1 = 'AAATTT'
>>> seq_2 = 'AATTTT'
>>> seq_3 = 'AAAAAT'
>>> tmpl = 'AATTTT'
>>> res_handler._get_best_overlap(seq_1, seq_1, tmpl)
('AAATTT', '|| |||')

# >>> res_handler._get_best_overlap(seq_1, seq_2, tmpl)
# ('AATTTT', '||||||')

# >>> res_handler._get_best_overlap(seq_1, seq_3, tmpl)
# ('AAATTT', '|| |||') 

```

### _get_tmpl_path(db, gene_id)
Input:
    db - database name
    gene_id - header of the gene corresponding to the hit
output:     
    db_file - path to the fasta file containing the database
    template of the hit
This function will output the path to the fasta file in the database that corresponds to the hit found by Blast

```python
>>> res_handler._get_tmpl_path("PointFinder", "geneA_1_CP0449499")
... #doctest: +ELLIPSIS
'.../geneA.fsa'

>>> res_handler._get_tmpl_path("DisinFinder", "geneA_1_CP0449499")
... #doctest: +ELLIPSIS
'.../disinfectants.fsa'

>>> res_handler._get_tmpl_path("ResFinder", "geneA_1_CP0449499")
... #doctest: +ELLIPSIS
'.../all.fsa'

```

### _keep_hit(hit_dict, current_hit, current_key, key_list)
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
>>> no_overlap = res_handler._keep_hit(gene_dict_no_overlap, current_hit, "rpoB_1_CP073768.1",
...                      keys_nooverlap)
>>> sorted(no_overlap[0])
['gyrA_1_CP073768.1', 'rpoB_1_CP073768.1']

>>> sorted(no_overlap[1])
[]

# overlap with same cal score keeping the longest alignment or both if no difference in aln_length
>>> equal_cal = res_handler._keep_hit(gene_dict_overlap, current_hit, "rpoB_1_CP073768.1",
...                      keys_overlap)

>>> sorted(equal_cal[0])
['gyrB_1_CP047010.1', 'rpoB_1_CP073768.1']

>>> sorted(equal_cal[1])
['gyrA_1_CP073768.1']

#exact overlap in contig - keeping the one with best identity or both if same id. 
>>> exact_overlap = res_handler._keep_hit(gene_exact_overlap, current_hit, "rpoB_1_CP073768.1",
...                      keys_exact)

>>> sorted(exact_overlap[0])
['gyrA_1_CP073768.1', 'gyrB_1_CP047010.1']

>>> sorted(exact_overlap[1])
['parC_1_CP084529.1', 'rpoB_1_CP073768.1']

```

### _get\_query\_align(hit, contig)
Input: <br>
    hit - an instance of a gene_dict[gene][0] dict. <br>
    contig - the input contig corresponding to the hit <br>
output: <br>
    updated query string <br>
The function will complete the query in the case that the
template gene is longer than the hit sequence. If there is corresponding
sequences in the contig these will be added otherwise '-' will be
appended accordingly.

#### setup

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

#### Test
```python
>>> res_handler._get_query_align(hit=hit1, contig=contig1)
('TTTTTTTCGTCAACCGGGGGGGGGGGGGGG', '       ||||||||               ')

>>> res_handler._get_query_align(hit=hit2, contig=contig2)
('---TTTTCGTCAACCGGGGGGG--------', '       ||||||||               ')

>>> res_handler._get_query_align(hit=hit3, contig=contig3)
('-------CGTCAACC---------------', '       ||||||||               ')

```

### _complete_template(hit, db)
Input: <br>
    hit - gene_dict[gene][0] object corresponding to a dict with hit 
    info <br>
    db - database name <br>
output: <br>
    hit - a dict containing updated values of sequence, template,
    and alignment <br>
This function will add the remaining sequence to the template if
the alignment do not cover the entire gene. It will update the hit
instance with the template covering the full gene, the query string,
and the alignment string with spaces in the remaining length of the
template length.

```python
# no contig is matching as the input fasta correspondst to hitA and not hitA2. only adds '-' to query 
>>> res_handler._complete_template(hitA2, 'PointFinder')
{'tmpl_start': 1, 'tmpl_end': 239, 'query_start': 1, 'query_end': 239, 'query_string': 'AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGAAGCTTGCTTCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCGGATGTGCCCAGAT-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------', 'aln_string': '|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       ', 'tmpl_string': 'AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGAAGCTTGCTTCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCTTA', 'aln_length': 239, 'contig_name': '16S-rrsB_1_CP067250.1', 'identity': 1.0, 'coverage': 0.1549935149, 'gene_length': 1542, 'hit_id': 'test_other_contig:8..15:16S-rrsB_1_CP067250.1'}

>>> res_handler._complete_template(hitA, 'ResFinder')
{'tmpl_start': 1, 'tmpl_end': 153, 'query_start': 1, 'query_end': 153, 'query_string': 'GTGTCCACACCACATCACGGCCGGCACGAGCTCGGCCAGAACTTCCTGTCCGATCGGCGCGTCATCGCCGATATCGTCGAAATCGTCTCGCGCACAAACGGTCCGATCATCGAGATCGGGGCGGGCGACGGCGCGCTGACCATACCCTTGCAACGACTCGCCCGCCCGCTCACCGCCGTCGAGGTCGACGCGCGGCGCGCGCGGCGGTTGGCGCAGCGCACCGCGAGATCCGCCCCGGGGCCTGCCTCGCGGCCCACCGAGGTCGTCGCCGCCGACTTCCTGCGCTACCCACTGCCCCGCTCACCCCACGTGGTCGTGGGC------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------', 'aln_string': '|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                ', 'tmpl_string': 'GTGTCCACACCACATCACGGCCGGCACGAGCTCGGCCAGAACTTCCTGTCCGATCGGCGCGTCATCGCCGATATCGTCGAAATCGTCTCGCGCACAAACGGTCCGATCATCGAGATCGGGGCGGGCGACGGCGCGCTGACCATACCCTTGCAACGACTCGCCCGCCCGCTCACCGCCGTCGAGGTCGACGCGCGGCGCGCGCGGCGGTTGGCGCAGCGCACCGCGAGATCCGCCCCGGGGCCTGCCTCGCGGCCCACCGAGGTCGTCGCCGCCGACTTCCTGCGCTACCCACTGCCCCGCTCACCCCACGTGGTCGTGGGCAACCTGCCGTTCCACCTCACCACCGCGATCCTGCGGCGACTGCTGCACGGTCCGGGCTGGACCACGGCCGTGCTGCTCATGCAGTGGGAGGTGGCCCGCCGACGCGCCGCGGTGGGCGGCGCCACCATGATGACCGCCCAGTGGTGGCCGTGGTTCGAATTCGGCCTTGCCCGAAAGGTTTCCGCGGCGAGCTTCACGCCGCGGCCCGCGGTCGACGCCGGACTGCTCACCATCACGCGCCGCAGCCGGCCGCTGGTCGACGTCGCGGACCGGGCGCGTTACCAGGCGCTGGTGCACCGCGTGTTCACCGGACGCGGACACGGCATGGCGCAGATCCTGCAACGGTTGCCCACGCCGGTGCCCCGCACTTGGTTGCGGGCCAACGGGATAGCACCGAACTCCCTGCCCCGCCAGTTGTCCGCGGCGCAGTGGGCGGCGCTGTTCGAGCAGACGCGTCTAACTGGTGCCCAACGGGTCGATCGTCCACGCGATGTACAGCACGGCCGCGCTCACCGTCGCCGTGGTGGCGAAGTCGATCGCCCGGCTACGCACCACAAGCAGACCGGCCCGGTCGTCGGTCAGCGCCAACCGCAGCGCGGCCGCGACGCCGACGCCGATCCCGATGACCAGCGCACCGCGCCGCCAGTAACCCGCCACCACCAGGGCGAACGCCGCGATGAAGATCAGGCCGACCACCAGGATCGGCCATTGACCGGCGAACACCTTGCGGGCGAATTCCTTTGGCGTCACGCCAGTTTCGACTCTTCGGCTTCGACGACGTTGGTCAGCAGGAAGGCGCGGGTCAACGGGCCCACGCCACCGGGGTTGGGCGACACGTGA', 'aln_length': 153, 'contig_name': 'erm(38)_first1-321', 'identity': 1.0, 'coverage': 0.13178294573643412, 'gene_length': 1161, 'hit_id': 'erm(38)_first1-321:8..15:erm(38)_1_AY154657.1'}


```


### _gene_overlap_comparison(pre_hit, next_hit)
 Input: <br>
    pre_hit - new gene_hit to compare <br>
    next_hit - gene_hit from gene_dict to compare with the new gene_hit <br>
 Output: <br>
 combined_gene - hit information in gene_hit format. <br>
 This function will compare two hits matching the same gene and combine
 them accordingly. 

```python
>>> res_handler._gene_overlap_comparison(hitA3, hitA4)
{'tmpl_start': 61, 'tmpl_end': 501, 'query_start': 61, 'query_end': 501, 'query_string': '????????????????????????????????????????????????????????????GTCATCGCCGATATCGTCGAAATCGTCTCGCGCACAAACGGTCCGATCATCGAGATCGGGGCGGGCGACGGCGCGCTGACCATACCCTTGCAACGACTCGCCCGCCCGCTCACCGCCGTCGAGGTCGACGCGCGGCGCGCGCGGCGGTTGGCGCAGCGCACCGCGAGATCCGCCCCGGGGCCTGCCTCGCGGCCCACCGAGGTCGTCGCCGCCGACTTCCTGCGCTACCCACTGCCCCGCTCACCCCACGTGGTCGTGGGCAACCTGCCGTTCCACCTCACCACCGCGATCCTGCGGCGACTGCTGCACGGTCCGGGCTGGACCACGGCCGTGCTGCTCATGCAGTGGGAGGTGGCCCGCCGACGCGCCGCGGTGGGCGGCGCCACCATGATGACCGCCCAGTGGTGGCCGTGGTTCGAATTCGGCCTTGCCCGAAAGGTT???????????????????????????????????????????????????????????????????????????????????????????????????---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------', 'aln_string': '                                                            |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    ', 'tmpl_string': 'GTGTCCACACCACATCACGGCCGGCACGAGCTCGGCCAGAACTTCCTGTCCGATCGGCGCGTCATCGCCGATATCGTCGAAATCGTCTCGCGCACAAACGGTCCGATCATCGAGATCGGGGCGGGCGACGGCGCGCTGACCATACCCTTGCAACGACTCGCCCGCCCGCTCACCGCCGTCGAGGTCGACGCGCGGCGCGCGCGGCGGTTGGCGCAGCGCACCGCGAGATCCGCCCCGGGGCCTGCCTCGCGGCCCACCGAGGTCGTCGCCGCCGACTTCCTGCGCTACCCACTGCCCCGCTCACCCCACGTGGTCGTGGGCAACCTGCCGTTCCACCTCACCACCGCGATCCTGCGGCGACTGCTGCACGGTCCGGGCTGGACCACGGCCGTGCTGCTCATGCAGTGGGAGGTGGCCCGCCGACGCGCCGCGGTGGGCGGCGCCACCATGATGACCGCCCAGTGGTGGCCGTGGTTCGAATTCGGCCTTGCCCGAAAGGTTTCCGCGGCGAGCTTCACGCCGCGGCCCGCGGTCGACGCCGGACTGCTCACCATCACGCGCCGCAGCCGGCCGCTGGTCGACGTCGCGGACCGGGCGCGTTACCAGGCGCTGGTGCACCGCGTGTTCACCGGACGCGGACACGGCATGGCGCAGATCCTGCAACGGTTGCCCACGCCGGTGCCCCGCACTTGGTTGCGGGCCAACGGGATAGCACCGAACTCCCTGCCCCGCCAGTTGTCCGCGGCGCAGTGGGCGGCGCTGTTCGAGCAGACGCGTCTAACTGGTGCCCAACGGGTCGATCGTCCACGCGATGTACAGCACGGCCGCGCTCACCGTCGCCGTGGTGGCGAAGTCGATCGCCCGGCTACGCACCACAAGCAGACCGGCCCGGTCGTCGGTCAGCGCCAACCGCAGCGCGGCCGCGACGCCGACGCCGATCCCGATGACCAGCGCACCGCGCCGCCAGTAACCCGCCACCACCAGGGCGAACGCCGCGATGAAGATCAGGCCGACCACCAGGATCGGCCATTGACCGGCGAACACCTTGCGGGCGAATTCCTTTGGCGTCACGCCAGTTTCGACTCTTCGGCTTCGACGACGTTGGTCAGCAGGAAGGCGCGGGTCAACGGGCCCACGCCACCGGGGTTGGGCGACACGTGA', 'aln_length': 441, 'contig_name': 'erm(38)_first1-321, erm(38)_middel201-600', 'coverage': 0.7333333333333333, 'gene_length': 1161, 'identity': 1.0, 'hit_id': 'erm(38)_first1-321:8..15:erm(38)_1_AY154657.1erm(38)_middel201-600:8..15:erm(38)_1_AY154657.1'}

# since the hit is not in the database there end will not be '?' but only '-'
# overlap hit 2 extends hit 1 : 111111111122222222
>>> res_handler._gene_overlap_comparison(hitG, hitH)
{'tmpl_start': 7, 'tmpl_end': 25, 'query_start': 7, 'query_end': 25, 'query_string': '??????AATCGTCAGCGTCCTCACC-----', 'aln_string': '      || |||| |||||||||||     ', 'tmpl_string': 'AGATCTAACCGTCTGCGTCCTCACCGCAAC', 'aln_length': 19, 'contig_name': 'test_overlap, test_overlapH', 'coverage': 0.6333333333333333, 'gene_length': 30, 'identity': 0.8947368421052632, 'hit_id': 'test_overlap:1..7:gyrA_1_CP073768.1test_overlap:1..7:gyrA_1_CP073768.1'}

# hit 1 overlap within hit 2: 222222111111222222 - not combined - using hit 2. 
>>> res_handler._gene_overlap_comparison(hitH, hitI)
{'tmpl_start': 9, 'tmpl_end': 30, 'query_start': 9, 'query_end': 30, 'query_string': 'CCGTCTGCGTCCTCACCTTTTT', 'aln_string': '||||||||||||||||||||||', 'tmpl_string': 'AGATCTAACCGTCTGCGTCCTCACCTTTTT', 'aln_length': 22, 'contig_name': 'test_overlap', 'coverage': 0.7333333333333333, 'gene_length': 30, 'identity': 1.0, 'hit_id': 'test_overlap:1..7:gyrA_1_CP073768.1'}

# hit 1 overlap within hit 2 and hit 2 contains Ns as part of the overlap. 
>>> res_handler._gene_overlap_comparison(hitH, hitI_N)
{'tmpl_start': 9, 'tmpl_end': 30, 'query_start': 9, 'query_end': 30, 'query_string': '????????CCGTCTGCGTCCTCACCTTTTT', 'aln_string': '        ||||||||||||||||||||||', 'tmpl_string': 'AGATCTAACCGTCTGCGTCCTCACCTTTTT', 'aln_length': 22, 'contig_name': 'test_overlapI, test_overlapH', 'coverage': 0.7333333333333333, 'gene_length': 30, 'identity': 1.0, 'hit_id': 'test_overlap:1..7:gyrA_1_CP073768.1test_overlap:1..7:gyrA_1_CP073768.1'}

# no overlap
>>> res_handler._gene_overlap_comparison(hitJ, hitG)
{'tmpl_start': 7, 'tmpl_end': 30, 'query_start': 7, 'query_end': 30, 'query_string': '??????AATCGTCAGCGNNCTCACCGCAAC', 'aln_string': '      || |||| |||  |||||||||||', 'tmpl_string': 'AGATCTAACCGTCTGCGTCCTCACCGCAAC', 'aln_length': 24, 'contig_name': 'test_overlap, test_overlapJ', 'coverage': 0.7333333333333333, 'gene_length': 30, 'identity': 0.8333333333333334, 'hit_id': 'test_overlap:1..7:gyrA_1_CP073768.1test_overlap:1..7:gyrA_1_CP073768.1'}

```

### _find_best_blast_hit(self, aligner, db_name):
 input: <br>
 aligner: blast result in BlastNAlignment object <br>
 db_name: string of the database name ('ResFinder', 'DisinFinder', 'PointFinder')<br>
 output: <br>
 gene_dict: a dict containing all hits<br>  
 The function finds the best sequence hits combined if multiple contigs
 hit the same genes, or the best hit if the contig hits multiple genes.
 Substitute to the find_best_seq in pointfinder.py and Blaster.init from
 cgecore from previous versions. 

#### setup 

```python


```

### _filter_and_standardize_result(self, combined_hits, std_result, db_name, finder, min_coverage, min_identity):
