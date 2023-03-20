ResFinder documentation
=============

ResFinder identifies acquired antimicrobial resistance genes in total or partial
sequenced isolates of bacteria.

## Important if you are updating from a previous ResFinder version

It is no longer recommended to clone the ResFinder bitbucket repository unless you plan to do development work on ResFinder.

Instead we recommend installing ResFinder using pip as described below.

There are several good reasons why the recommended installation procedure has changed, among those are the increasing size of the repository that has risen to several hundreds of megabytes, due to the long history of ResFinder. Its easier for users. And it makes sure your installation will be a tested release of the application.

## Installation
ResFinder consists of an application and 1-3 databases. The databases can be used without the application, but not the other way around. Below ResFinder, the application, will be installed first and then the databases will be installed and configured to work with ResFinder the application.

### Dependencies

ResFinder uses two external alignment tools that must be installed.

* BLAST 
* KMA

#### BLAST
If you don't want to specify the path of BLAST every time you run ResFinder, make sure that "blastn" is in you PATH or set the environment variable specified in the "Environment Variables Table" in this README.

Blastn can be obtained from:
```url
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
```

```bash
# Example of how to set the environment variable in the bash shell. Remember this is only temporary, if you want it set every time you log in you need to add this line to your .bashrc, .zshrc file.
export CGE_BLASTN="/path/to/some/dir/blastn"
```

#### KMA
If you don't want to specify the path of KMA every time you run ResFinder, make sure that KMA is in you PATH or set the environment variable specified in the "Environment Variables Table" in this README.

KMA can be obtained from:
```url
https://bitbucket.org/genomicepidemiology/kma.git
```

```bash
# Example of how to set the environment variable in the bash shell. Remember this is only temporary, if you want it set every time you log in you need to add this line to  your .bashrc, .zshrc file.
export CGE_KMA="/path/to/some/dir/kma/kma"
```

### Install ResFinder the application using pip

**Important**: This will install ResFinder in the environment where you run pip and potenitally update the python modules ResFinder depends on. It is recommended to run ResFinder in its own environment, in order to avoid breaking existing installations and prevent ResFinder from getting broken by future unrelated pip installations. This is described in the optional step below.

#### Optional: Create virtual environment ####

Go to the location where you want to store your environment.

```bash

# Create environment
python3 -m venv resfinder_env

# Activate environment
source resfinder_env/bin/activate

# When you are finished using ResFinder deactivate the environment
deactivate

```

#### Install ResFinder ####

```bash

pip install resfinder

```

#### Databases
If you don't want to specify the path to the databases every time you run ResFinder, you need to set the environment variable specified in the "Environment Variables Table" in this README.

Go to the location where you want to store the databases. Clone the datbases you need.

**Note**: We are currently working on hosting tarballed versions of the databases that can be downloaded, so that cloning can be avoided.

```bash

git clone https://bitbucket.org/genomicepidemiology/resfinder_db/
git clone https://bitbucket.org/genomicepidemiology/pointfinder_db/
git clone https://bitbucket.org/genomicepidemiology/disinfinder_db/

```

Set approximate environment variables.

```bash

# Example of how to set the environment variable in the bash shell. Remember this is only temporary, if you want it set every time you log in you need to add this line to for example your .bashrc file.
export CGE_RESFINDER_RESGENE_DB="/path/to/some/dir/resfinder_db"
export CGE_RESFINDER_RESPOINT_DB="/path/to/some/dir/pointfinder_db"
export CGE_DISINFINDER_DB="/path/to/some/dir/disinfinder_db"

```

### Install ResFinder with Docker

The ResFinder application and the 3 databases has been build into a single image on docker hub named "genomicepidemiology/resfinder". Below is an example run, where the current working directory is bound to the container "/app" path which is the container working directory.

```bash

docker run -v "$(pwd):/app" genomicepidemiology/resfinder -ifa data/test_isolate_01.fa -o test1 -s ecoli --acquired --point

```

### Test data
Test data can be found in the sub-directory tests/data

## Usage

You can run resfinder command line using python.

**NOTE**: Species should be entered with their full scientific names (e.g. "escherichia coli"), using quotation marks, not case sensitive.
          An attempt has been made to capture some deviations like "ecoli" and "e.coli", but it is far from all deviations that will be captured.


```bash

# Example of running resfinder
python -m resfinder -o path/to/outdir -s "Escherichia coli" -l 0.6 -t 0.8 --acquired --point -ifq test_isolate_01_*

# The program can be invoked with the -h option
usage: __main__.py [-h] [-ifa INPUTFASTA] [-ifq INPUTFASTQ [INPUTFASTQ ...]] [--nanopore] -o OUTPUTPATH [-j OUT_JSON] [-b BLASTPATH] [-k KMAPATH] [-s SPECIES] [--ignore_missing_species] [-db_res DB_PATH_RES]
                   [-db_res_kma DB_PATH_RES_KMA] [-acq] [-ao ACQ_OVERLAP] [-l MIN_COV] [-t THRESHOLD] [-d] [-db_disinf DB_PATH_DISINF] [-db_disinf_kma DB_PATH_DISINF_KMA] [-c] [-db_point DB_PATH_POINT]
                   [-db_point_kma DB_PATH_POINT_KMA] [-g SPECIFIC_GENE [SPECIFIC_GENE ...]] [-u] [-l_p MIN_COV_POINT] [-t_p THRESHOLD_POINT] [--ignore_indels] [--ignore_stop_codons] [-v] [--pickle]

options:
  -h, --help            show this help message and exit
  -ifa INPUTFASTA, --inputfasta INPUTFASTA
                        Input fasta file.
  -ifq INPUTFASTQ [INPUTFASTQ ...], --inputfastq INPUTFASTQ [INPUTFASTQ ...]
                        Input fastq file(s). Assumed to be single-end fastq if only one file is provided, and assumed to be paired-end data if two files are provided.
  --nanopore            If nanopore data is used
  -o OUTPUTPATH, --outputPath OUTPUTPATH
                        Output directory. If it doesnt exist, it will be created.
  -j OUT_JSON, --out_json OUT_JSON
                        Specify JSON filename and output directory. If the directory doesnt exist, it will be created.
  -b BLASTPATH, --blastPath BLASTPATH
                        Path to blastn
  -k KMAPATH, --kmaPath KMAPATH
                        Path to KMA
  -s SPECIES, --species SPECIES
                        Species in the sample
  --ignore_missing_species
                        If set, species is provided and --point flag is set, will not throw an error if no database is found for the provided species. If species is not found. Point mutations will silently be ignored.
  -db_res DB_PATH_RES, --db_path_res DB_PATH_RES
                        Path to the databases for ResFinder.
  -db_res_kma DB_PATH_RES_KMA, --db_path_res_kma DB_PATH_RES_KMA
                        Path to the ResFinder databases indexed with KMA. Defaults to the value of the --db_res flag.
  -acq, --acquired      Run resfinder for acquired resistance genes
  -ao ACQ_OVERLAP, --acq_overlap ACQ_OVERLAP
                        Genes are allowed to overlap this number of nucleotides. Default: 30.
  -l MIN_COV, --min_cov MIN_COV
                        Minimum (breadth-of) coverage of ResFinder within the range 0-1.
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold for identity of ResFinder within the range 0-1.
  -d, --disinfectant    Run resfinder for disinfectant resistance genes
  -db_disinf DB_PATH_DISINF, --db_path_disinf DB_PATH_DISINF
                        Path to the databases for DisinFinder.
  -db_disinf_kma DB_PATH_DISINF_KMA, --db_path_disinf_kma DB_PATH_DISINF_KMA
                        Path to the DisinFinder databases indexed with KMA. Defaults to the value of the --db_res flag.
  -c, --point           Run pointfinder for chromosomal mutations
  -db_point DB_PATH_POINT, --db_path_point DB_PATH_POINT
                        Path to the databases for PointFinder
  -db_point_kma DB_PATH_POINT_KMA, --db_path_point_kma DB_PATH_POINT_KMA
                        Path to the PointFinder databases indexed with KMA. Defaults to the value of the --db_path_point flag.
  -g SPECIFIC_GENE [SPECIFIC_GENE ...], --specific_gene SPECIFIC_GENE [SPECIFIC_GENE ...]
                        Specify genes existing in the database to search for - if none is specified all genes are included in the search.
  -u, --unknown_mut     Show all mutations found even if in unknown to the resistance database
  -l_p MIN_COV_POINT, --min_cov_point MIN_COV_POINT
                        Minimum (breadth-of) coverage of Pointfinder within the range 0-1. If None is selected, the minimum coverage of ResFinder will be used.
  -t_p THRESHOLD_POINT, --threshold_point THRESHOLD_POINT
                        Threshold for identity of Pointfinder within the range 0-1. If None is selected, the minimum coverage of ResFinder will be used.
  --ignore_indels       Ignore frameshift-causing indels in Pointfinder.
  --ignore_stop_codons  Ignore premature stop codons in Pointfinder.
  -v, --version         Show programs version number and exit
  --pickle              Create a pickle dump of the Isolate object. Currently needed in the CGE webserver. Dependency and this option is being removed.

```

### Environment Variables

Environment variables recognized by ResFinder, the flag they replace and the default value for the flag. Provided commandline flags will always take precedence. Set environment variables takes precedence over default flag values.

Additional Environment variables can be added by appending entries to the file named "environment_variables.md".

#### Environment Variables Table

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

### Species Abbreviations

ResFinder understands the species abbreviations listed in the Species Abbreviations Table. Additional species abbreviations can be added by appending entries to the file "species_abbreviations.md".

#### Species Abbreviations Table

| Species                       | Abbreviation            |
|-------------------------------|-------------------------|
| campylobacter jejuni          | c. jejuni               |
| campylobacter jejuni          | c.jejuni                |
| campylobacter jejuni          | c jejuni                |
| campylobacter jejuni          | cjejuni                 |
| campylobacter coli            | c. coli                 |
| campylobacter coli            | c.coli                  |
| campylobacter coli            | c coli                  |
| campylobacter coli            | ccoli                   |
| escherichia coli              | e. coli                 |
| escherichia coli              | e.coli                  |
| escherichia coli              | e coli                  |
| escherichia coli              | ecoli                   |
| salmonella enterica           | s. enterica             |
| salmonella enterica           | s.enterica              |
| salmonella enterica           | s enterica              |
| salmonella enterica           | senterica               |

### Web-server

A webserver implementing the methods is available at the [CGE
website](http://www.genomicepidemiology.org/) and can be found here:
https://cge.food.dtu.dk/services/ResFinder/

### ResFinder result files

ResFinder outputs several files. A brief description of these is given below.

* pheno_table_species.txt: table with species specific AMR phenotypes.
* pheno_table.txt: table with all AMR phenotypes.
* PointFinder_prediction.txt: tab seperated table. 1 is given to a predicted resistance against an antibiotic class, 0 is given to not resistance detected.
* PointFinder_results.txt: tab seperated table with predicted point mutations leading to antibiotic resistance.
* PointFinder_table.txt: predicted point mutations grouped into genes to which they belong.
* ResFinder_Hit_in_genome_seq.fsa: fasta sequence of resistance gene hits found in the input data (query).
* ResFinder_Resistance_gene_seq.fsa: fasta sequence of resistance gene hits found in the database (reference).
* ResFinder_results_table.txt: predicted resistance genes grouped by antibiotic class.
* ResFinder_results_tab.txt: tab seperated table with predicted resistance genes.
* ResFinder_results.txt: predicted resistance genes grouped by antibiotic class and hit alignments to reference resistance genes.
* <input_filename>.json: Output written to a CGE standardized json file. All results can be derived from this file. The format is defined here: https://bitbucket.org/genomicepidemiology/cgelib/src/master/src/cgelib/output/templates_json/beone/


Citation
=======

When using the method please cite:

ResFinder 4.0 for predictions of phenotypes from genotypes.  
Bortolaia V, Kaas RS, Ruppe E, Roberts MC, Schwarz S, Cattoir V, Philippon A, Allesoe RL, Rebelo AR, Florensa AR, Fagelhauer L,
Chakraborty T, Neumann B, Werner G, Bender JK, Stingl K, Nguyen M, Coppens J, Xavier BB, Malhotra-Kumar S, Westh H, Pinholt M,
Anjum MF, Duggett NA, Kempf I, Nykasenoja S, Olkkola S, Wieczorek K, Amaro A, Clemente L, Mossong J, Losch S, Ragimbeau C, Lund O, Aarestrup FM.
Journal of Antimicrobial Chemotherapy. 2020 Aug 11.  
PMID: 32780112			doi: 10.1093/jac/dkaa345  
[Epub ahead of print]  

References
=======

1. Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. BLAST+: architecture and applications. BMC Bioinformatics 2009; 10:421.
2. Clausen PTLC, Aarestrup FM, Lund O. Rapid and precise alignment of raw reads against redundant databases with KMA. BMC Bioinformatics 2018; 19:307.

License
=======

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
