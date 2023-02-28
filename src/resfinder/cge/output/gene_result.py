#!/usr/bin/env python3
import random
import string
import re

from ..phenotype2genotype.res_profile import PhenoDB


class GeneResult(dict):
    def __init__(self, res_collection, res, db_name, conf=None):
        """
            Input:
                res_collection: Result object created by the cgelib package.
                res: Custom dictionary containing information about a single hit
                     from ResFinder.
                db_name: 'ResFinder' or 'PointFinder' or 'DisinFinder'
            Method creates a seq_region dict as defined in the BeOne template:
            https://bitbucket.org/genomicepidemiology/cgelib/src/master/src/
            cgelib/output/templates_json/beone/
            from res.
        """
        self.db_name = db_name
        self["type"] = "seq_region"
        self["gene"] = True

        self["ref_id"] = res["templateID"]
        self["ref_id"] = PhenoDB.if_promoter_rename(self["ref_id"])
        self["name"], self.variant, self["ref_acc"] = (
            GeneResult._split_sbjct_header(self["ref_id"]))
        self["ref_database"] = [res_collection.get_db_key(db_name)[0]]

        alignment = res['aln_scheme']
        start = re.search(r"^[\s_-]*(\|+)", alignment).start(1)
        end = re.search(r"\|+([\s_-]*)$", alignment).start(1)
        self["alignment_length"] = end - start
        self["identity"] = float(res["template_identity"])
        self["ref_seq_lenght"] = int(res['template_length'])
        depth = res.get("depth", None)
        if depth != None:
            self["depth"] = float(depth)
        self["ref_start_pos"] = start + 1
        self["ref_end_pos"] = end + 1
        self["query_id"] = res.get("contig_name", None)  # Positional essential
        self["query_start_pos"] = res.get("query_start", None)  # Positional essential
        self["query_end_pos"] = res.get("query_end", None)  # Positional essential
        self["pmids"] = []
        self["notes"] = []

        # adding alignment patterns for Resfinder and Disinfinder results
        if conf and conf.output_aln:
            self["query_string"] = res["query_aln"]
            self["alignment_string"] = res["aln_scheme"]
            self["ref_string"] = res["template_aln"]

        coverage = res["template_coverage"]
        self["coverage"] = float(coverage)

        self["grade"] = GeneResult.calc_gene_grade(coverage=self["coverage"],
                                                   identity=self["identity"])

        self.remove_NAs()
        uniqueness = self._get_unique_gene_key(res_collection)
        self["key"] = uniqueness

    def remove_NAs(self):
        """
            Remove all entries containing NA og None as values.

            Removing None is not necessary as the Result object will ignore all
            entries with None values.
        """
        na_keys = []
        for key, val in self.items():
            if(val == "NA" or val is None):
                na_keys.append(key)
        for key in na_keys:
            del self[key]

    @staticmethod
    def calc_gene_grade(coverage: float, identity: float) -> int:
        if coverage == 100.0 and identity == 100.0:
            return 3
        elif coverage == 100.0:
            return 2
        else:
            return 1

    @staticmethod
    def get_rnd_unique_gene_key(gene_key, res_collection,
                                minimum_gene_key, delimiter):
        """
            Input:
                gene_key: None-unique key
                res_collection: Result object created by the cgelib package.
                minimum_key: Key prefix
                delimiter: String used as delimiter inside the returned key.
            Output:
                gene_key: Unique key (string)

            If gene_key is found in res_collection. Creates a unique key by
            appending a random string ton minimum_gene_key.
        """
        while(gene_key in res_collection["seq_regions"]):
            rnd_str = GeneResult.random_string(str_len=4)
            gene_key = ("{key}{deli}{rnd}"
                        .format(key=minimum_gene_key, deli=delimiter,
                                rnd=rnd_str))
        return gene_key

    @staticmethod
    def random_string(str_len=4):
        """
            Output:
                random string of length 'str_len'

            Return a random string of the provided length. The string will only
            consist of lowercase ascii letters.
        """
        letters = string.ascii_lowercase
        return ''.join(random.choice(letters) for i in range(str_len))

    @staticmethod
    def _split_sbjct_header(header):
        """
            Input:
                header: database entry header (ref_header/subject_header)
            Output:
                template: name of entry (string)
                variant: Variant interger (string) or None
                acc: Accession number given by sequnce database (string) or None

            Splits the input header by underscores and returns first list item
            as template. If list is > 1 then variant and acc will also return
            strings. If not they return None.
        """
        sbjct = header.split("_")
        template = sbjct[0]

        if(len(sbjct) > 1):
            variant = sbjct[1]
            acc = "_".join(sbjct[2:])
        else:
            variant = None
            acc = None

        return (template, variant, acc)

    def _get_unique_gene_key(self, res_collection, delimiter=";;"):
        """
            Input:
                res_collection: Result object created by the cgelib package.
                delimiter: String used as delimiter inside the returned key.
            Output:
                key: Unique key for hit

            Creates a unique key for GeneResult instance. Key format depends on
            database. If gene result is considered indentical to an existing
            gene result in the provided res_collection, it will not create a new
            key. Two restults are considered identical if they have the same
            query_id, query_start_pos and query_end_pos. If it is mapping
            results the query_* doesn't exist, and results will never be
            considered identical.
        """
        if "ref_acc" not in self:
            gene_key = ("{name}".format(**self))
        else:
            gene_key = ("{name}{deli}{var}{deli}{ref_acc}".format(
                        deli=delimiter, var=self.variant, **self))
        # Attach random string if key already exists
        minimum_gene_key = gene_key
        if gene_key in res_collection["seq_regions"]:

            query_id = self.get("query_id", "NA")

            # Query id == "NA" when FASTQ
            if(query_id == "NA"):
                res_collection["seq_regions"][gene_key]["ref_database"].extend(
                    self["ref_database"])
                gene_key = None
            # Query id != "NA" when FASTA
            elif (self["query_id"]
                    != res_collection["seq_regions"][gene_key]["query_id"]
                  or self["query_start_pos"]
                    != res_collection["seq_regions"][gene_key]["query_start_pos"]
                  or self["query_end_pos"]
                    != res_collection["seq_regions"][gene_key]["query_end_pos"]):
                gene_key = GeneResult.get_rnd_unique_gene_key(
                    gene_key, res_collection, minimum_gene_key, delimiter)
            else:
                res_collection["seq_regions"][gene_key]["ref_database"].extend(
                    self["ref_database"])
                gene_key = None

        return gene_key
