#!/usr/bin/env python3
import random
import string
import json

from ..phenotype2genotype.res_profile import PhenoDB
from ..phenotype2genotype.feature import ResGene, ResMutation

from .exceptions import DuplicateKeyError
from .gene_result import GeneResult
from .seq_variation_result import SeqVariationResult
from .phenotype_result import PhenotypeResult


def add_gene_result_if_key_not_None(gene_result, res_collection):
    '''
        Input:
            gene_result: GeneResult object (seq_region dict)
            res_collection: Result object created by the cge core module.
        adds gene_result to res_collection if gene_result['key'] is different
        from None.
            '''
    if gene_result["key"] is None:
        return
    elif gene_result["key"] not in res_collection["seq_regions"]:
        res_collection.add_class(cl="seq_regions", **gene_result)
    else:
        raise DuplicateKeyError(
            "About to overwrite dict entry. This should not be "
            "happening as all keys are supposed to be unique."
            "Non-unique key was: {}".format(gene_result["key"]))



class ResFinderResultHandler():

    @staticmethod
    def standardize_results(res_collection, res, ref_db_name, conf):
        """
            Input:
                res_collection: Result object created by the cge core module.
                res: Custom dictionary of results from ResFinder
                ref_db_name: 'ResFinder' or 'PointFinder' or 'DisinFinder'

            Method loads the given res_collection with results from res.
        """
        for db_name, db in res.results.items():
            if(db_name == "excluded"):
                continue

            if(db == "No hit found"):
                continue

            for unique_id, hit_db in db.items():
                if(unique_id in res.results["excluded"]):
                    continue
                gene_result = GeneResult(res_collection, hit_db, ref_db_name,
                                         conf, res)

                add_gene_result_if_key_not_None(gene_result, res_collection)

    @staticmethod
    def load_res_profile(res_collection, isolate, amr_abbreviations):
        """
            Input:
                res_collection: Result object created by the cge core module.
                isolate: Isolate object

            Method loads the given res_collection with results from res.
        """
        # For each antibiotic class
        for ab_class in isolate.resprofile.phenodb.antibiotics.keys():
            # For each antibiotic in current class
            for phenodb_ab in isolate.resprofile.phenodb.antibiotics[ab_class]:
                phenotype = PhenotypeResult(phenodb_ab, isolate)

                # Isolate is resistant towards the antibiotic
                if(phenodb_ab in isolate.resprofile.resistance):
                    ResFinderResultHandler._load_resistant_phenotype(
                        phenotype, isolate, phenodb_ab, res_collection)

                res_collection.add_class(cl="phenotypes", **phenotype)

        amr_sum = ResFinderResultHandler.create_amr_summary_str(
            res_collection, amr_abbreviations)
        res_collection.add(**{"result_summary": amr_sum})

    @staticmethod
    def _load_resistant_phenotype(phenotype, isolate, antibiotic,
                                  res_collection):
        phenotype.set_resistant(True)
        isolate_ab = isolate.resprofile.resistance[antibiotic]

        for feature_lst in isolate_ab.features.values():
            for feature_entry in feature_lst:
                # feature_entry is either a Feature object or a dict of Feature
                # objects.
                try:
                    for unique_id, feature in feature_entry.items():
                        ResFinderResultHandler._add_feat_to_phenotype_and_res(
                            res_collection, isolate, feature, phenotype)
                # Feature object has no items attribute.
                except AttributeError:
                    ResFinderResultHandler._add_feat_to_phenotype_and_res(
                        res_collection, isolate, feature_entry, phenotype)

    @staticmethod
    def _add_feat_to_phenotype_and_res(res_collection, isolate, feature,
                                       phenotype):
        if isinstance(feature, ResGene) or isinstance(feature, ResMutation):
            if isinstance(feature, ResMutation):
                if feature.unique_id in phenotype.get("seq_variations", ()):
                    return
            phenotype.add_feature(res_collection, isolate, feature)

    @staticmethod
    def create_amr_summary_str(res_collection, amr_abbreviations):
        amr_list = []
        for key, phenotype in res_collection["phenotypes"].items():
            if(phenotype["amr_resistant"] is True
               and phenotype["amr_species_relevant"] is True):
                amr_name = phenotype["amr_resistance"].capitalize()
                amr = amr_abbreviations.get(amr_name,
                                            [phenotype["amr_resistance"]])
                amr_list.append(amr[0])
        if(amr_list):
            out_str = "_".join(amr_list)
        else:
            out_str = ""
        return out_str


class PointFinderResultHandler():

    @staticmethod
    def standardize_results(res_collection, res, ref_db_name):
        """
            Input:
                res_collection: Result object created by the cge core module.
                res: Custom dictionary of results from PointFinder
                ref_db_name: 'ResFinder' or 'PointFinder' or 'DisinFinder'

            Method loads the given res_collection with results from res.
        """
        for gene_name, db in res.items():
            # Ignore information in excluded dict
            if(gene_name == "excluded"):
                continue

            # Ignore genes found in excluded dict
            if gene_name in res["excluded"]:
                continue
            if(isinstance(db, str)):
                if(db == "No hit found"):
                    continue
                if db.startswith("Gene found with coverage"):
                    continue

            gene_results = []

            # For BLAST results
            db_hits = db.get("hits", {})

            # For KMA results
            if(not db_hits):
                id = db["sbjct_header"]
                db_hits[id] = db

            for unique_id, hit_db in db_hits.items():

                # Ignore genes found in excluded dict
                if(unique_id in res["excluded"]):
                    continue
                gene_result = GeneResult(res_collection, hit_db, ref_db_name)

                add_gene_result_if_key_not_None(gene_result, res_collection)
                if gene_result['key'] is None:
                    continue

                gene_results.append(gene_result)

                # KMA hits
                mismatches = hit_db.get("mis_matches", None)
                # BLAST hits
                if(mismatches is None):
                    mismatches = db["mis_matches"]

                for mismatch in mismatches:
                    seq_var_result = SeqVariationResult(
                        res_collection, mismatch, gene_results, ref_db_name)
                    res_collection.add_class(cl="seq_variations",
                                                **seq_var_result)
