#!/usr/bin/env python3
from __future__ import division
import sys
import os

import pandas as pd


class ResFinder():
    # Variables used by methods to distinguish results created by different
    # methods.
    TYPE_BLAST = "blast"
    TYPE_KMA = "kma"

    def __init__(self, db_conf_file, notes, db_path, db_path_kma,
                 databases=None, pheno_file=None):
        """
        """
        self.db_path = db_path
        self.db_path_kma = db_path_kma

        self.configured_dbs = dict()
        self.kma_db_files = None
        self.load_db_config(db_conf_file=db_conf_file)

        self.databases = []
        self.load_databases(databases=databases)

        self.phenos = dict()
        if pheno_file is None:
            self.load_notes(notes=notes)
        else:
            self.load_phenos(pheno_file=pheno_file)

        self.blast_results = None

    def load_phenos(self, pheno_file):
        df_phenos = pd.read_csv(pheno_file, sep="\t", usecols=[0,2])
        df_phenos["Gene_accession no."] = df_phenos["Gene_accession no."].str.split('_').str[0]
        df_phenos = df_phenos["Phenotype"].groupby([df_phenos["Gene_accession no."]]).apply(set).apply(list).apply(', '.join).reset_index()
        dict_phenos = df_phenos.set_index("Gene_accession no.").to_dict()['Phenotype']
        self.phenos = dict_phenos

    def load_notes(self, notes):
        with open(notes, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith("#"):
                    continue
                else:
                    tmp = line.split(":")

                    self.phenos[tmp[0]] = "%s %s" % (tmp[1], tmp[2])

                    if(tmp[2].startswith("Alternate name; ")):
                        self.phenos[tmp[2][16:]] = "%s %s" % (tmp[1], tmp[2])
        print(self.phenos)

    def load_databases(self, databases):
        """
        """
        # Check if databases and config file are correct/correponds
        if databases == '':
            sys.exit("Input Error: No database was specified!\n")
        elif databases is None:
            # Choose all available databases from the config file
            self.databases = self.configured_dbs.keys()
        else:
            # Handle multiple databases
            databases = databases.split(',')
            # Check that the ResFinder DBs are valid
            for db_prefix in databases:
                if db_prefix in self.configured_dbs:
                    self.databases.append(db_prefix)
                else:
                    sys.exit("Input Error: Provided database was not "
                             "recognised! (%s)\n" % db_prefix)

    def load_db_config(self, db_conf_file):
        """
        """
        extensions = []
        with open(db_conf_file) as f:
            for line in f:
                line = line.strip()

                if not line:
                    continue

                if line[0] == '#':
                    if 'extensions:' in line:
                        extensions = [s.strip()
                                      for s in line.split('extensions:')[-1]
                                      .split(',')]
                    continue

                tmp = line.split('\t')
                if len(tmp) != 3:
                    sys.exit(("Input Error: Invalid line in the database"
                              " config file!\nA proper entry requires 3 tab "
                              "separated columns!\n%s") % (line))

                db_prefix = tmp[0].strip()
                name = tmp[1].split('#')[0].strip()
                description = tmp[2]

                # Check if all db files are present
                for ext in extensions:
                    db = "%s/%s.%s" % (self.db_path, db_prefix, ext)
                    if not os.path.exists(db):
                        sys.exit(("Input Error: The database file (%s) "
                                  "could not be found!") % (db))

                if db_prefix not in self.configured_dbs:
                    self.configured_dbs[db_prefix] = []
                self.configured_dbs[db_prefix].append(name)

        if len(self.configured_dbs) == 0:
            sys.exit("Input Error: No databases were found in the "
                     "database config file!")

        # Loading paths for KMA databases.
        for drug in self.configured_dbs:
            kma_db = self.db_path_kma + drug
            self.kma_db_files = [kma_db + ".b", kma_db + ".length.b",
                                 kma_db + ".name.b", kma_db + ".align.b"]
