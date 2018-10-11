===================
ResFinder
===================

This project documents ResFinder service


Documentation
=============

## What is it?

The ResFinder service contains one perl script *resfinder.pl* which is the
script of the latest version of the ResFinder service. ResFinder identifies
acquired antimicrobial resistance genes in total or partial sequenced isolates
of bacteria.

This repository also contains a python script resfinder.py which is  a new version 
of ResFinder, but not yet running on the CGE server. This program was added because
it uses a newer version of blastn,  which, in contrary from the blastall version 
that the perl script uses, is avail to download.

## Content of the repository
1. resfinder.pl - the program
2. resfinder.py - (same program with using a available blastn version - blastn-2.2.26+)
3. test.fsa     - test fasta file

## Installation

Setting up ResFinder
```bash
# Go to wanted location for resfinder
cd /path/to/some/dir
# Clone and enter the resfinder directory
git clone https://git@bitbucket.org/genomicepidemiology/resfinder.git
cd resfinder
```

```bash

# Installing up the ResFinder database
cd /path/to/some/dir

# Clone and enter the resfinder directory
git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git
cd resfinder_db

```

##Installing dependencies (for python script):

The BlastAll and FormatDB that the perl script uses are no longer available 
for downloading through ncbi. Therefor we have provided the resfinder.py 
scriot that uses Blastn instead. Note, this is not not script that is running 
on the CGE server. The CGE server is running the perl script using BlastAll


###Download Blastn and BioPython
```url
http://biopython.org/DIST/docs/install/Installation.html
```
```url
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
```
### Install the cgecore module to python3
```bash
pip3 install cgecore
```


## Test the scripts and database
```bash
cd /path/to/test_dir
cp /path/to/resfinder/test.fsa .
python3 /path/to/resfinder/resfinder.py -o . -p /path/to/resfinder_db \
-b /path/to/blastn -i test.fsa -d aminoglycoside -k 90.00 -l 0.60
```

## Usage

The program can be invoked with the -h option to get help and more information 
of the service.

```bash
                    age: resfinder.py [-h] [-i INPUTFILE] [-1 FASTQ1] [-2 FASTQ2] [-o OUT_PATH]
                    [-b BLAST_PATH] [-p DB_PATH] [-k KMA_PATH]
                    [-q DB_PATH_KMA] [-d DATABASES] [-l MIN_COV]
                    [-t THRESHOLD]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTFILE, --inputfile INPUTFILE
                        Input file
  -1 FASTQ1, --fastq1 FASTQ1
                        Raw read data file 1.
  -2 FASTQ2, --fastq2 FASTQ2
                        Raw read data file 2 (only required if data is paired-
                        end).
  -o OUT_PATH, --outputPath OUT_PATH
                        Path to blast output
  -b BLAST_PATH, --blastPath BLAST_PATH
                        Path to blast
  -p DB_PATH, --databasePath DB_PATH
                        Path to the databases
  -k KMA_PATH, --kmaPath KMA_PATH
                        Path to KMA
  -q DB_PATH_KMA, --databasePathKMA DB_PATH_KMA
                        Path to the directories containing the KMA indexed
                        databases. Defaults to the directory 'kma_indexing'
                        inside the databasePath directory.
  -d DATABASES, --databases DATABASES
                        Databases chosen to search in - if none are specified
                        all are used
  -l MIN_COV, --min_cov MIN_COV
                        Minimum coverage default 0.6
  -t THRESHOLD, --threshold THRESHOLD
                        Blast threshold for identity
                        default minimum 0.9 
```


## Web-server

A webserver implementing the methods is available at the [CGE 
website](http://www.genomicepidemiology.org/) and can be found here:
https://cge.cbs.dtu.dk/services/ResFinder/

## Documentation

The documentation available as of the date of this release can be found at
https://bitbucket.org/genomicepidemiology/resfinder/overview.


Citation
=======

When using the method please cite:

Identification of acquired antimicrobial resistance genes.
Zankari E, Hasman H, Cosentino S, Vestergaard M, Rasmussen S, Lund O, Aarestrup 
FM, Larsen MV.
J Antimicrob Chemother. 2012 Jul 10.
PMID: 22782487         doi: 10.1093/jac/dks261
[Epub ahead of print]


License
=======

Copyright (c) 2014, Ole Lund, Technical University of Denmark
All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
