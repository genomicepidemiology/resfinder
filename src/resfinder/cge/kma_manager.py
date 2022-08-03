#!/usr/bin/env python3
import re

from cgelib.alignment.aligner import KMAAligner
from cgelib.alignment.KMA.read_files import KMA_Result
from cgelib.alignment.read_alignment import KMAAlignment

class KMAManager():

    def __init__(self, min_cov, threshold, databases, db_path_kma,
                 sample_name=''):

        self.kma_results = None
        self.thresholds = dict(min_cov=min_cov*100,
                               min_id=threshold*100)
        self.databases = databases
        self.db_path = db_path_kma

        #todo: find out how to add filename to pointfinder output path.
        if sample_name != "":
            self.sample_name = f"_{sample_name}"


    def calibrate_params(self, params, resultfiles):
        params['1t1'] = params['kma_1t1']

        for file in resultfiles:
            if file == 'Result' or file == 'Fragments':
                continue
            params[file.lower()] = True

        #setting input files and -ipe or -int flag
        if params['inputfile_2'] is not None:
            params['input_ipe'] = [params['inputfile_1'], params['inputfile_2']]
        else:
            params['input_int'] = params['inputfile_1']

        #remove inputs that do not fit into kma
        params.pop('inputfile_1')
        params.pop('inputfile_2')
        params.pop('kma_1t1')

        return params


    def run_KMAAligner(self, conf, result_files, params):

        kma_aligner = KMAAligner(
            result_file=result_files)

        kma_params = self.calibrate_params(params, result_files)

        kma_aligner.set_aligner_params(**kma_params)

        db_path = self.db_path
        dbs = []
        for db in self.databases:
            dbs.append(f"{db_path}/{db}")

        combined_stdout_stderr = kma_aligner(variable_iter="template_db",
                                             values_iter=dbs)
        kma_aligner.fit_alignment()

        # self.conform_result(kma_aligner)

        return kma_aligner


    def filter_hits(self, cov, template_id, hit):
        exclude_reasons = []

        if (cov < self.thresholds['min_cov']
                or template_id < self.thresholds['min_id']):
            exclude_reasons.append(cov)
            exclude_reasons.append(template_id)

        if exclude_reasons:
            self.kma_results['excluded'][hit] = exclude_reasons


    def conform_result(self, kma_aligner):

        filenames = [x.split('/')[-1] for x in kma_aligner.alignment_files]
        out_path = [x.rsplit('/', 1)[0] for x in kma_aligner.alignment_files][0]

        # get alignments + output similar to iterator_hits
        kma_alignment = KMAAlignment(
            output_path=out_path,
            filenames=filenames,
            result_file=["Result", "Alignment"])

        self.kma_results = dict()
        self.kma_results["excluded"] = dict()

        #Needs to be after assignment. otherwise the dict will not work 
        for database in self.databases:
            self.kma_results[database] = dict()


        for kmahit in kma_alignment.parse_hits():
            hit = kmahit['templateID']
            template_id = float(kmahit['template_identity'])
            template_cov = float(kmahit['template_coverage'])
            query_aln = kmahit['query_aln']
            template_aln = kmahit['template_aln']
            aln_scheme = kmahit['aln_scheme']

            db = kmahit['template_file'].split('-', 2)[-1]

            self.filter_hits(template_cov, template_id, hit)

            self.kma_results[db][hit] = dict()
            self.kma_results[db][hit]['sbjct_length'] = float(kmahit['template_length'])
            self.kma_results[db][hit]["perc_coverage"] = template_cov
            self.kma_results[db][hit]["perc_ident"] = template_id
            self.kma_results[db][hit]["sbjct_header"] = hit
            self.kma_results[db][hit]["cal_score"] = float(kmahit['q_value'])
            self.kma_results[db][hit]["depth"] = float(kmahit['depth'])
            self.kma_results[db][hit]["p_value"] = float(kmahit['p_value'])
            self.kma_results[db][hit]["contig_name"] = "NA"
            self.kma_results[db][hit]["query_start"] = "NA"
            self.kma_results[db][hit]["query_end"] = "NA"

            start = re.search("^-*(\w+)", query_aln).start(1)
            end = re.search("\w+(-*)$", query_aln).start(1)

            self.kma_results[db][hit]['sbjct_string'] = template_aln[start:end]
            self.kma_results[db][hit]['query_string'] = query_aln[start:end]
            self.kma_results[db][hit]['homo_string'] = aln_scheme[start:end]

            # Save align start and stop positions relative to
            # subject sequence
            self.kma_results[db][hit]['sbjct_start'] = start + 1
            self.kma_results[db][hit]["sbjct_end"] = end + 1
            self.kma_results[db][hit]["HSP_length"] = end - start

            # Count gaps in the alignment
            self.kma_results[db][hit]["gaps"] = (
                    self.kma_results[db][hit]['sbjct_string'].count("-")
                    + self.kma_results[db][hit]['query_string'].count("-"))

        return self.kma_results
