#!/usr/bin/env python3

import argparse
import os.path
import re
import shutil
from signal import *
import tempfile
import sys
import subprocess
from itertools import chain

from .feature import Feature, ResGene, Mutation, ResMutation
from .phenodbpoint import PhenoDBPoint
from .res_profile import PhenoDB, ResProfile, FeatureGroup
from .dbhit import DBHit


class Isolate(dict):
    """ An isolate class is a dict of Features.
    """
    NO_AB_CLASS = "No class defined"

    def __init__(self, name, species=None, amr_panel_file=None):
        self.name = name
        self.resprofile = None
        self.feature_classes = {}
        self.species = None
        if(species is not None):
            self.amr_panel = Isolate.load_amr_panel(species, amr_panel_file)
        else:
            self.amr_panel = None

    @staticmethod
    def load_amr_panel(species, panel_file):
        """
        """
        panels = {}
        inclusions = {}

        with open(panel_file, "r") as fh:
            panel_name = None
            for line in fh:

                line = line.rstrip()

                # Skip empty lines and comments
                if not line or line.startswith("#"):
                    continue

                # Get panel name
                new_panel_name = Isolate._get_panel_name(panels, line)
                if new_panel_name is not None:
                    panel_name = new_panel_name
                    continue

                if Isolate._get_inclusions(panel_name, line, inclusions):
                    continue

                # Get Antibiotics
                if(panel_name):
                    Isolate._get_antibiotics(line, panel_name, panels)

        Isolate._merge_inclusions(panels, inclusions)

        species_key = Isolate.check_panel_name(species, panels)

        return set(panels[species_key]) if species_key else None

    @staticmethod
    def _get_antibiotics(line, panel_name, panels):
        """ Stores list of antimicrobials in dict of panel names. """
        tmp_list = panels.get(panel_name, [])
        tmp_list.append(line.lower())
        panels[panel_name] = tmp_list

    @staticmethod
    def _get_panel_name(panels, line):
        match_panel = re.search(r':Panel:\s{0,1}(.+)$', line)
        if(match_panel):
            panel_name = match_panel.group(1).lower()
            panels[panel_name] = []
            return panel_name
        else:
            return None

    @staticmethod
    def _get_inclusions(panel_name, line, inclusions):
        match_inclusion = re.search(r':Include:\s{0,1}(.+)$', line)
        if(match_inclusion):
            include_panel = match_inclusion.group(1).lower()
            tmp_list = inclusions.get(panel_name, [])
            tmp_list.append(include_panel)
            inclusions[panel_name] = tmp_list
            return True
        else:
            return False

    @staticmethod
    def _merge_inclusions(panels, inclusions):
        for panel, include_list in inclusions.items():
            panel_list = panels[panel]
            include_abs = []
            for incl_panel in include_list:
                include_abs = include_abs + panels[incl_panel]
            panels[panel] = panel_list + include_abs

    @staticmethod
    def check_panel_name(name, panels):
        """ Panel names are expected to consist of "Genus species" or
            only "Genus". The method checks the name against the loaded
            panel names and returns the panel name that matches the
            given name. If no panel name matches it checks if the first
            word of the given name matches any of the genus panel
            names, and considers that a match if found.

            Returns False if no match is found
            Returns None if no panels has been loaded
            Returns the panel name that match
        """
        if(not panels):
            return None

        name = name.lower()
        if(name in panels):
            return name

        genus_name = " ".split(name)[0]
        if(genus_name in panels):
            return genus_name

        return False

    def load_resfinder_tab(self, tabbed_output, phenodb):
        with open(tabbed_output, "r") as fh:
            while(True):
                line = fh.readline()
                if(not line):
                    break

                line = line.rstrip()
                if(not line):
                    continue

                db_name = line
                second_line = fh.readline().rstrip()

                if(second_line == "No hit found"):
                    continue

                # At this point second line must be headers, and are skipped.

                res_hit = fh.readline().rstrip()

                while(res_hit):
                    hit_list = res_hit.split("\t")
                    match_length, ref_length = hit_list[2].split("/")
                    start_ref, end_ref = hit_list[4].split("..")
                    accno = hit_list[8]
                    gene_name = hit_list[0]
                    unique_id = "{0}_{1}".format(gene_name, accno)
                    hit = DBHit(name=gene_name, identity=hit_list[1],
                                match_length=match_length,
                                ref_length=ref_length, start_ref=start_ref,
                                end_ref=end_ref, acc=accno,
                                db="resfinder")

                    start_feat, end_feat = hit_list[6].split("..")

                    if(start_feat == "NA"):
                        start_feat = None
                        end_feat = None

                    phenotypes, unique_id = self.get_phenotypes(phenodb,
                                                                unique_id)
                    ab_class = set()
                    if(phenotypes):
                        for p in phenotypes:
                            for ab in p.antibiotics:
                                ab_class.update(ab.classes)
                    else:
                        ab_class.add(db_name)

                    gene_feat = ResGene(unique_id=unique_id,
                                        seq_region=hit_list[5],
                                        start=start_feat, end=end_feat,
                                        hit=hit, ab_class=ab_class)

                    if(unique_id in self):
                        temp_list = self[unique_id]
                        temp_list.append(gene_feat)
                        self[unique_id] = temp_list
                    else:
                        self[unique_id] = [gene_feat]

                    res_hit = fh.readline().rstrip()

    @staticmethod
    def get_phenodb_id(feat_res_dict, type):
        """
            Input:
                feat_res_dict: entry in the dict 'genes' or 'seq_variations'
                               from a cgecore Result object.
                type: Which dict in the Result object, either 'seq_regions' or
                      'seq_variations'
            Output: Key formatted as the ones in the PhenoDB object.

            Method is used to reformat an id from a Result object, so that it
            can be used to query a PhenoDB object.

            TODO: Make this method obsolete by harmonizing the ids in PhenoDB
                  with the ones from the Result object.
        """
        if(type == "seq_variations"):
            var_aa = feat_res_dict.get("var_aa", None)
            var_codon = feat_res_dict.get("var_codon", None)
            codon_change = feat_res_dict.get("codon_change", None)
            indels = feat_res_dict['deletion'] or feat_res_dict['insertion']

            # Not point mutation
            if(var_aa is None and var_codon is None):
                return feat_res_dict["seq_regions"][0]
            # RNA mutation(single nucleotide)
            elif(len(var_codon)==1 and codon_change is None):
                nuc_format = (f"{feat_res_dict['seq_regions'][0]}"
                              f"_{feat_res_dict['ref_start_pos']}"
                              f"_{feat_res_dict['var_codon']}")
                aa_format = ""
                return nuc_format, aa_format
            # indels
            elif (indels):
                nuc_format = (f"{feat_res_dict['seq_regions'][0]}"
                              f"_{feat_res_dict['ref_end_pos']}"
                              f"_{feat_res_dict['nuc_change']}")
                aa_format = (f"{feat_res_dict['seq_regions'][0]}"
                             f"_{feat_res_dict['ref_start_pos']}"
                             f"_{feat_res_dict['var_aa']}")
                return nuc_format, aa_format
            # Amino acid mutation
            else:
                nuc_format = ""
                aa_format = (f"{feat_res_dict['seq_regions'][0]}"
                             f"_{feat_res_dict['ref_start_pos']}"
                             f"_{feat_res_dict['var_aa']}")
                return nuc_format, aa_format

        elif(type == "seq_regions"):
            return ("{}_{}".format(feat_res_dict["name"],
                                  feat_res_dict["ref_acc"]), "")

    def load_finder_results(self, std_table, phenodb, type):
        """
            Input:
                std_table: Result object from cgecore module loaded with
                           ResFinder and/or Pointfinder results.
                phenodb: PhenoDB object
                type: 'seq_regions' or 'seq_variations', for ResFinder and
                      PointFinder results, respectively.

            Method loads Isolate object
        """
        for key, feat_info in std_table[type].items():
            # Skip genes from PointFinder database (not resistance genes).
            if(type == "seq_regions"
               and any(re.search("PointFinder", entry)
                       for entry in feat_info["ref_database"])):
                continue

            unique_id_nuc, unique_id_aa = Isolate.get_phenodb_id(feat_info,
                                                                 type)
            phenotypes, unique_id = self.get_phenotypes(phenodb, unique_id_nuc,
                                                        unique_id_aa, type)
            feat_list = self.get(unique_id, [])
            if(phenotypes):
                for p in phenotypes:
                    res_feature = self.new_res_feature(type, feat_info,
                                                       unique_id,
                                                       unique_id_aa,
                                                       unique_id_nuc, p)
                    if(res_feature not in feat_list):
                        feat_list.append(res_feature)
            else:
                res_feature = self.new_res_feature(type, feat_info,
                                                   unique_id,
                                                   unique_id_aa,
                                                   unique_id_nuc)
                feat_list.append(res_feature)
            if unique_id_nuc == "":
                self[unique_id_aa] = feat_list
            else:
                self[unique_id_nuc] = feat_list

    def get_phenotypes(self, phenodb, unique_id_nuc, unique_id_aa, type):
        """
            Input:
                phenodb: PhenoDB object
                unique_id: string in key format fitting the phenodb objects
            Output:
                returns a phenodb object based on the unique id and the unique
                id used.
            Method modifies unique id to account for mutation type (AA/NUC)
            using the pointfinder database.
        """
        if (phenodb.mut_type_is_defined
                and type == "seq_variations"):
            unique_id_aa = unique_id_aa + '_AA'
            unique_id_nuc = unique_id_nuc + '_NUC'

        phenotype_aa = phenodb.get(unique_id_aa, None)
        phenotype_nuc = phenodb.get(unique_id_nuc, None)

        if phenotype_nuc:
            phenotypes = phenotype_nuc
            unique_id_found = unique_id_nuc
        else:
            phenotypes = phenotype_aa
            unique_id_found = unique_id_aa

        return phenotypes, unique_id_found

    def new_res_feature(self, type, feat_info, unique_id, aa_format, nuc_format,
                        phenotype=None):
        ab_class = set()

        if phenotype is None:
            ab_class.add(Isolate.NO_AB_CLASS)
            db_pheno = None
        else:
            for ab in phenotype.antibiotics:
                ab_class.update(ab.classes)
            db_pheno = phenotype.res_database

        if(type == "seq_regions"):
            res_feature = self.new_res_gene(feat_info, unique_id, ab_class,
                                            phenotype)
        elif(type == "seq_variations"):
            res_feature = self.new_res_mut(feat_info, unique_id, ab_class,
                                           phenotype, nuc_format, aa_format)

        ResProfile.update_classes_dict_of_feature_sets(
            self.feature_classes, res_feature)
        return res_feature

    def new_res_mut(self, feat_info, unique_id, ab_class, phenotype, nuc_format,
                    aa_format):
        ref_aa = feat_info.get("ref_aa", None)

        if(ref_aa is None or ref_aa.upper() == "NA"):
            nucleotide_mut = True
        else:
            nucleotide_mut = False

        if phenotype:
            feat_res = ResMutation(
                unique_id=unique_id,
                nuc_format=nuc_format,
                aa_format=aa_format,
                seq_region=";;".join(feat_info["seq_regions"]),
                pos=feat_info["ref_start_pos"],
                ref_codon=feat_info["ref_codon"],
                mut_codon=feat_info["var_codon"],
                ref_aa=feat_info.get("ref_aa", None),
                mut_aa=feat_info.get("var_aa", None),
                isolate=self,
                insertion=feat_info["insertion"],
                deletion=feat_info["deletion"],
                end=feat_info["ref_end_pos"],
                nuc=nucleotide_mut,
                ab_class=ab_class,
                pmids=phenotype.pmid,
                notes=phenotype.notes,
                ref_db=phenotype.res_database
            )
        else:
            feat_res = ResMutation(
                unique_id=unique_id,
                nuc_format= nuc_format,
                seq_region=";;".join(feat_info["seq_regions"]),
                pos=feat_info["ref_start_pos"],
                ref_codon=feat_info["ref_codon"],
                mut_codon=feat_info["var_codon"],
                ref_aa=feat_info.get("ref_aa", None),
                mut_aa=feat_info.get("var_aa", None),
                isolate=self,
                insertion=feat_info["insertion"],
                deletion=feat_info["deletion"],
                end=feat_info["ref_end_pos"],
                nuc=nucleotide_mut,
                ab_class=ab_class
            )

        return feat_res

    def new_res_gene(self, gene_info, unique_id, ab_class, phenotype):
        hit = DBHit(name=gene_info["name"],
                    identity=gene_info["identity"],
                    match_length=gene_info["alignment_length"],
                    ref_length=gene_info["ref_seq_lenght"],
                    start_ref=gene_info["ref_start_pos"],
                    end_ref=gene_info["ref_end_pos"],
                    acc=gene_info["ref_acc"],
                    depth=gene_info.get("depth", None),
                    db=phenotype.res_database)

        query_start = gene_info.get("query_start_pos", None)
        query_end = gene_info.get("query_end_pos", None)
        query_id = gene_info.get("query_id", None)

        feat_res = ResGene(unique_id=unique_id,
                           seq_region=query_id,
                           start=query_start,
                           end=query_end,
                           hit=hit,
                           ab_class=ab_class,
                           pmids=phenotype.pmid,
                           notes=phenotype.notes,
                           ref_db=phenotype.res_database)
        return feat_res

    def calc_res_profile(self, phenodb):
        """
        """
        features = self.values()
        # Flatten features.
        features = list(chain.from_iterable(features))
        self.resprofile = ResProfile(features, phenodb)

    def profile_to_str_table(self, header=False):
        """
        """
        output_str = ""

        if(header):
            output_str = (
                "# ResFinder phenotype results.\n"
                "# Sample: " + self.name + "\n"
                "# \n"
                "# The phenotype 'No resistance' should be interpreted with\n"
                "# caution, as it only means that nothing in the used\n"
                "# database indicate resistance, but resistance could exist\n"
                "# from 'unknown' or not yet implemented sources.\n"
                "# \n"
                "# The 'Match' column stores one of the integers 0, 1, 2, 3.\n"
                "#      0: No match found\n"
                "#      1: Match < 100% ID AND match length < ref length\n"
                "#      2: Match = 100% ID AND match length < ref length\n"
                "#      3: Match = 100% ID AND match length = ref length\n"
                "# If several hits causing the same resistance are found,\n"
                "# the highest number will be stored in the 'Match' column.\n"
                "\n"
            )
            output_str += ("# Antimicrobial\t"
                           "Class\t"
                           "WGS-predicted phenotype\t"
                           "Match\t"
                           "Genetic background\n")

        # For each antibiotic class
        for ab_class in self.resprofile.phenodb.antibiotics.keys():
            # For each antibiotic in current class
            for ab_db in self.resprofile.phenodb.antibiotics[ab_class]:
                output_str += ("{ab:s}\t{cl:s}"
                               .format(ab=ab_db.name, cl=ab_class))

                # Isolate is resistant towards the antibiotic
                if(ab_db in self.resprofile.resistance):
                    ab = self.resprofile.resistance[ab_db]
                    output_str += "\tResistant"

                    # Find the resistance causing gene with the best match
                    # Mutations will always have best match as they are only
                    # either present or absent.
                    best_match = 0
                    for unique_id in ab.features:
                        feature_list = ab.features[unique_id]
                        for feature in feature_list:
                            if(isinstance(feature, Mutation)
                               or isinstance(feature, FeatureGroup)):
                                best_match = 3
                            # Note: Mutations do not have "hits"
                            elif(feature.hit.match_category > best_match):
                                best_match = feature.hit.match_category

                                output_str += "\t" + str(best_match)

                    gene_list = ab.get_gene_namewacc(tostring=True)
                    mut_list = ab.get_mut_namewannot(tostring=True)

                    if(gene_list and mut_list):
                        gene_mut_str = gene_list + " " + mut_list
                    elif(gene_list):
                        gene_mut_str = gene_list
                    elif(mut_list):
                        gene_mut_str = mut_list
                    else:
                        gene_mut_str = ""

                    output_str += "\t" + gene_mut_str + "\n"

                # TODO: delete elif clause.
                # Isolate is susceptibile towards the antibiotic
                elif(ab_db.name in self.resprofile.susceptibile):
                    ab = self.resprofile.susceptibile[ab_db.name]
                    # Genetic background is not written if susceptibile
                    gene_list = ""
                    # Uncomment next line to write genetic background for susc.
                    # gene_list = ab.get_gene_namewacc(tostring=True)
                    output_str += "\tNo resistance\t0\t" + gene_list + "\n"
                else:
                    output_str += "\tNo resistance\t0\t\n"

        if(self.resprofile.missing_db_features):
            output_str += ("\n# WARNING: Missing features from phenotype "
                           "database:\n")
            output_str += "# Feature_ID\tRegion\tDatabase\tHit\n"

            for feature in self.resprofile.missing_db_features:
                output_str += ("{}\t{}\t"
                               .format(feature.unique_id, feature.seq_region))

                if(feature.hit is None):
                    output_str += "\t"
                else:
                    output_str += str(feature.hit.db) + "\t" + feature.hit.name

                output_str += "\n"

        if(self.resprofile.unknown_db_features):
            output_str += ("\n# WARNING: Features with unknown phenotype \n")
            output_str += "# Feature_ID\tRegion\tDatabase\tHit\tClass\n"

            for feature in self.resprofile.unknown_db_features:
                output_str += ("{}\t{}\t"
                               .format(feature.unique_id, feature.seq_region))

                if(feature.hit is None):
                    output_str += "\t"
                else:
                    output_str += str(feature.hit.db) + "\t" + feature.hit.name

                output_str += "\t" + ",".join(feature.ab_class)

                output_str += "\n"

        return output_str
