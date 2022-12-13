#!/usr/bin/env python3
from cgelib.alignment.aligner import BlastNAligner

class BlastManager():

    def __init__(self, databases, db_path_blast):
        self.databases = databases
        self.db_path = db_path_blast

    # 1 generate input for BlastNALigner
    #   1.2 generate filenames - corresponding to each antimicrobial in format
#           blasttest2betalactam5.xml for xml - think we went for tsv output
    #create BlastNAligner object
    #

    def run_BlastNAligner(self, result_files, params):
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

