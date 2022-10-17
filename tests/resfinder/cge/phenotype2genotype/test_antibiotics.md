# res_profile class tests

**NOT COMPLETE** Missing lots of tests.

## Setup

### Create feature (ResGene, ResMutation) and Antibiotics objects

```python

>>> from src.resfinder.cge.phenotype2genotype.feature import ResGene
>>> from src.resfinder.cge.phenotype2genotype.feature import ResMutation
>>> rg = ResGene(unique_id="blaOXA-384_KF986263", start=1, end=90,
...              isolate=None, ab_class=["beta-lactam"],
...              ref_db="ResFinder")
>>> rm = ResMutation(unique_id="gyrA_81_d", seq_region="gyrA", pos=81,
...                  ref_codon="ggt", mut_codon="gat", ref_aa="g",
...                  mut_aa="d", isolate=None, nuc=False,
...                  ab_class=["fluoroquinolone"], ref_db="PointFinder")

```

## Initialise

```python

>>> from src.resfinder.cge.phenotype2genotype.res_profile import Antibiotics
>>> ab_m1 = Antibiotics(name="ciprofloxacin", classes=["fluoroquinolone"])
>>> ab_g1 = Antibiotics(name="carbapenem", classes=["beta-lactam"])
>>> ab_m2 = Antibiotics(name="nalidixic acid", classes=["fluoroquinolone"],
...                     feature=rm)

```

## Methods

### add_feature(self, feature)

```python

>>> ab_m1.add_feature(rg)
>>> ab_m1.features[rg.unique_id][0].unique_id
'blaOXA-384_KF986263'

>>> ab_g1.add_feature(rm)
>>> ab_g1.features[rm.unique_id][0].unique_id
'gyrA_81_d'

```

### get_gene_names(self, _list)

```python

>>> ab_m1.get_gene_names()
{'blaOXA-384_KF986263': 'Not Available'}

>>> ab_g1.get_gene_names()
{}

```
