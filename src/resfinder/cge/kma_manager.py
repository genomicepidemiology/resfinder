#!/usr/bin/env python3
import re

from cgelib.alignment.aligner import KMAAligner
from cgelib.alignment.KMA.read_files import KMA_Result
from cgelib.alignment.read_alignment import KMAAlignment

class KMAManager():
    def __init__(self, databases, db_path_kma,
                 sample_name=''):

        self.kma_results = None
        self.databases = databases
        self.db_path = db_path_kma

        #todo: find out how to add filename to pointfinder output path.
        if sample_name != "":
            self.sample_name = f"_{sample_name}"

    def run_KMAAligner(self, result_files, params):
        """
        Input:
            result_files: a list of wanted output types
            params: parameters for running KMA.
        output:
            kma_aligner: a aligner object created by cgelib when calling
            kma_aligner from a KMAAligner object
        """
        kma_aligner = KMAAligner(
            result_file=result_files)

        if params["input_int"] == None:
            params.pop('input_int')
        if params['input_ipe'] == None:
            params.pop('input_ipe')

        kma_aligner.set_aligner_params(**params)

        db_path = self.db_path
        dbs = []
        for db in self.databases:
            dbs.append(f"{db_path}/{db}")

        combined_stdout_stderr = kma_aligner(variable_iter="template_db",
                                             values_iter=dbs)
        kma_aligner.fit_alignment()

        return kma_aligner

