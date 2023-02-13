#!/usr/bin/env python3


from ..phenotype2genotype.feature import ResGene, ResMutation


class PhenotypeResult(dict):
    def __init__(self, antibiotic, isolate):
        self["type"] = "phenotype"
        self["category"] = "amr"
        self["key"] = antibiotic.name
        self["amr_classes"] = antibiotic.classes
        self["amr_resistance"] = antibiotic.name
        self["amr_resistant"] = False
        self["amr_species_relevant"] = PhenotypeResult.get_amr_relevance(
            antibiotic.name, isolate)
        self["grade"] = 0

    @staticmethod
    def get_amr_relevance(ab, isolate):
        if(isolate.amr_panel is None):
            return True
        elif(ab in isolate.amr_panel):
            return True
        else:
            return False

    def set_resistant(self, res):
        self["amr_resistant"] = res

    def add_feature(self, res_collection, isolate, feature):
        """
            Input:
                res_collection: Result object created by the cgelib package.
                isolate: Isolate object. Must have been loaded with at
                         resistance profile using Isolate.load_finder_results
                         and then Isolate.calc_res_profile.
                feature: Either ResGene or ResMutation object (inherit feature)
        """
        # Get all keys in the result that matches the feature in question.
        # Most of the time this will be a one to one relationship.
        # However if several identical features has been found in a sample,
        # they will all have different keys, but identical ref ids.
        ref_id, type = PhenotypeResult.get_ref_id_and_type(feature, isolate)
        feature_keys = PhenotypeResult.get_keys_matching_ref_id(
            ref_id, res_collection[type])

        # Add keys to phenotype results
        pheno_feat_keys = self.get(type, [])
        for feature_key in feature_keys:
            if feature_key not in pheno_feat_keys:
                pheno_feat_keys.append(feature_key)
        self[type] = pheno_feat_keys

        # Add phenotype keys to feature results
        features = res_collection[type]
        for feat_key in feature_keys:
            feat_result = features[feat_key]
            pheno_keys = feat_result.get("phenotypes", [])
            if self["key"] not in pheno_keys:
                pheno_keys.append(self["key"])
            feat_result["phenotypes"] = pheno_keys

            # Compare the grade of the feature with what is already recorded,
            # keep the greater value.
            # seq_variations has not grade value as it will always be 3
            feat_grade = feat_result.get("grade", 3)
            if feat_grade > self["grade"]:
                self["grade"] = feat_grade


        # Add unique PMIDs to feature results
        if(feature.pmids is not None):
            for pmid in feature.pmids:
                if pmid not in feat_result["pmids"]:
                    feat_result["pmids"].append(pmid)

        # Add unique Notes to feature results
        if(feature.notes is not None):
            if feature.notes not in feat_result["notes"]:
                feat_result["notes"].append(feature.notes)
        db_name = feature.ref_db
        if(type == "seq_regions"):
            db_key = res_collection.get_db_key(db_name)[0]
        elif(type == "seq_variations"):
            db_key = res_collection.get_db_key(db_name)[0]
        self["ref_database"] = [db_key]

    @staticmethod
    def get_ref_id_and_type(feature, isolate):
        """
            Input:
                feature: Either ResGene or ResMutation object (inherit feature)
                isolate: Isolate object. Must have been loaded with at
                         resistance profile using Isolate.load_finder_results
                         and then Isolate.calc_res_profile.
            Output (tuple):
                ref_id: id to identity the feature in relevant database
                type: 'seq_regions' or 'seq_variations'
        """
        type = None
        ref_id = None
        if(isinstance(feature, ResGene)):
            type = "seq_regions"
            ref_id = isolate.resprofile.phenodb.id_to_idwithvar[
                feature.unique_id]
        elif(isinstance(feature, ResMutation)):
            type = "seq_variations"
            if feature.nuc_format == '':
                ref_id = feature.aa_format
            else:
                ref_id = feature.nuc_format
        return (ref_id, type)

    @staticmethod
    def get_keys_matching_ref_id(ref_id, res_collection):
        """
            Input:
                ref_id: Key/ID to search for
                res_collection: Result object created by the cgelib package.
                                The Result object used with this method is
                                found within the main Result object and is
                                accessed either by main_result["seq_regions"]
                                or by main_result["seq_variations"].
            Output (list):
                out_keys: Returns a list of keys where the
                          input ref_id == res_collection[key]["ref_id"].
                          No matches will return an empty list.
        """
        out_keys = []
        for key, results in res_collection.items():
            if(ref_id == results["ref_id"]):
                out_keys.append(key)

        return out_keys

    def get_pmid_key(self, isolate):
        phenodb = isolate.resprofile.phenodb
        isolate_ab = isolate.resprofile.resistance[self["key"]]
        pubmed_ids = isolate_ab.get_pubmed_ids(phenodb)
        return pubmed_ids
