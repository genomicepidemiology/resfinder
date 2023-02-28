#!/usr/bin/env python3
from cgelib.alignment.aligner import BlastNAligner

class BlastManager():

    def __init__(self, databases, db_path_blast):
        self.databases = databases
        self.db_path = db_path_blast

    def run_BlastNAligner(self, result_files, params):
        """
        Input:
            result_files: blast output format (outfmt flag) option
            params: parameters for running blast.
        output:
            blast_aligner: a aligner object created by cgelib when calling
            blast_aligner from a BlastNAligner object
        this function will call blast with the disired parameters and output
        format.
        """
        blast_aligner = BlastNAligner(
            result_file=result_files)

        blast_aligner.set_aligner_params(**params)

        db_path = self.db_path
        dbs = []
        for db in self.databases:
            dbs.append(f"{db_path}/{db}.fsa")

        combined_stdout_stderr = blast_aligner(variable_iter="subject",
                                             values_iter=dbs)
        blast_aligner.fit_alignment()

        return blast_aligner

