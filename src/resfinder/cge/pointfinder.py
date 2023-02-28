#!/usr/bin/env python3
#
# Program: 	PointFinder-3.0
# Author: 	Camilla Hundahl Johnsen
#
# Dependencies: KMA or NCBI-blast together with BioPython.

import os
import re
import sys
import math

from resfinder.cge.phenotype2genotype.res_sumtable import PanelNameError


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


class GeneListError(Exception):
    """ Raise when a specified gene is not found within the gene list.
    """
    def __init__(self, message, *args):
        self.message = message
        # allow users initialize misc. arguments as any other builtin Error
        super(PanelNameError, self).__init__(message, *args)


class PointFinder():
    # Variables used by methods to distinguish results created by different
    # methods.
    TYPE_BLAST = "blast"
    TYPE_KMA = "kma"

    GENE_FILE = "genes.txt"
    RNA_GENE_FILE = "RNA_genes.txt"
    PHENOTYPE_FILE = "phenotypes.txt"
    # Deprecated variable, use phenotype (kept for legacy dbs)
    RES_OVERVIEW_FILE = "resistens-overview.txt"

    def __init__(self, db_path, species, gene_list=None, ignore_indels=False,
                 ignore_stop_codons=False):
        """
        """
        self.species = species
        self.specie_path = db_path
        self.RNA_gene_list = []
        self.ignore_indels = ignore_indels
        self.ignore_stop_codons = ignore_stop_codons

        self.gene_list = PointFinder.get_file_content(
            f"{self.specie_path}/{self.GENE_FILE}")

        if os.path.isfile(f"{self.specie_path}/{self.RNA_GENE_FILE}"):
            self.RNA_gene_list = PointFinder.get_file_content(
                f"{self.specie_path}/{self.RNA_GENE_FILE}")

        # Creat user defined gene_list if given
        if(gene_list is not None):
            self.gene_list = self.get_user_defined_gene_list(gene_list)

        # Depends on database format, current or legacy
        if os.path.exists(f"{self.specie_path}/{self.PHENOTYPE_FILE}"):

            self.known_mutations, self.drug_genes, self.known_stop_codon = (
                self.get_db_mutations(
                    f"{self.specie_path}/{self.PHENOTYPE_FILE}",
                    self.gene_list))
        elif os.path.exists(f"{self.specie_path}/{self.RES_OVERVIEW_FILE}"):
            self.known_mutations, self.drug_genes, self.known_stop_codon = (
                self.get_db_old_mutations(
                    f"{self.specie_path}/{self.RES_OVERVIEW_FILE}",
                    self.gene_list))
        else:
            raise OSError("The pointfinder db does not have neither the"
                          f"'{self.PHENOTYPE_FILE}' file (current database) or"
                          "the '{self.RES_OVERVIEW_FILE}' file "
                          "(legacy database)")

    def get_user_defined_gene_list(self, gene_list):
        genes_specified = []
        for gene in gene_list:
            # Check that the genes are valid
            if gene not in self.gene_list:
                raise(GeneListError(
                    "Input Error: Specified gene not recognised "
                    "(%s)\nChoose one or more of the following genes:"
                    "\n%s" % (gene, "\n".join(self.gene_list))))
            genes_specified.append(gene)
        # Change the gene_list to the user defined gene_list
        return genes_specified

    @staticmethod
    def get_db_names(db_root_path):
        out_lst = []
        for entry in os.scandir(db_root_path):
            if not entry.name.startswith('.') and entry.is_dir():
                out_lst.append(entry.name)
        return tuple(out_lst)



    @staticmethod
    def discard_unknown_muts(key, mismatches, phenodb):

        known_muts = PointFinder._get_known_mis_matches(
            key, mismatches, phenodb)
        return known_muts

    @staticmethod
    def _get_known_mis_matches(entry_key, mis_matches, phenodb):

        gene_db_id = entry_key.replace("_", ";;")
        known = []
        for mis_match in mis_matches:
            if mis_match[0] != 'sub':
                mis_match_key_nuc = (f"{gene_db_id}_{mis_match[2]}"
                                     f"_{mis_match[3].lower()}")
                mis_match_key_aa = (f"{gene_db_id}_{mis_match[1]}"
                                    f"_{mis_match[-1].lower()}")
            elif mis_match[4].startswith('r.'):#RNA substitution
                mis_match_key_aa = ""
                mis_match_key_nuc = (f"{gene_db_id}_{mis_match[2]}"
                                     f"_{mis_match[3].lower()}")
            else:
                mis_match_key_aa = (f"{gene_db_id}_{mis_match[2]}"
                                    f"_{mis_match[-1].lower()}")
                mis_match_key_nuc = ""

            if phenodb.mut_type_is_defined:
                mis_match_key_nuc = mis_match_key_nuc + '_NUC'
                mis_match_key_aa = mis_match_key_aa + '_AA'

            if (mis_match_key_nuc in phenodb
                    or mis_match_key_aa in phenodb):
                known.append(mis_match)
        return known

    @staticmethod
    def get_db_old_mutations(mut_db_path, gene_list):
        """
        This function opens the file resistens-overview.txt, and reads
        the content into a dict of dicts. The dict will contain
        information about all known mutations given in the database.
        This dict is returned.
        This function reads the old pointfinder database
        """

        # Initiate variables
        known_mutations = dict()
        known_stop_codon = dict()
        drug_genes = dict()
        indelflag = False
        stopcodonflag = False

        # Go throug mutation file line by line

        with open(mut_db_path, "r") as fh:
            mutdb_file = fh.readlines()
        mutdb_file = [line.strip() for line in mutdb_file]

        for line in mutdb_file:
            # Ignore headers and check where the indel section starts
            if line.startswith("#"):
                if "indel" in line.lower():
                    indelflag = True
                elif "stop codon" in line.lower():
                    stopcodonflag = True
                else:
                    stopcodonflag = False
                continue

            mutation = [data.strip() for data in line.split("\t")]

            gene_ID = mutation[0]

            # Only consider mutations in genes found in the gene list
            if gene_ID in gene_list:
                gene_name = mutation[1]
                mut_pos = int(mutation[2])
                ref_codon = mutation[3]                     # Ref_nuc (1 or 3?)
                ref_aa = mutation[4]                        # Ref_codon
                alt_aa = mutation[5].split(",")             # Res_codon
                res_drug = mutation[6].replace("\t", " ")
                pmid = mutation[7].split(",")

                # Check if stop codons are predictive for resistance
                if stopcodonflag is True:
                    if gene_ID not in known_stop_codon:
                        known_stop_codon[gene_ID] = {"pos": [],
                                                     "drug": res_drug}
                    known_stop_codon[gene_ID]["pos"].append(mut_pos)

                # Add genes associated with drug resistance to drug_genes dict
                drug_lst = res_drug.split(",")
                drug_lst = [d.strip().lower() for d in drug_lst]
                for drug in drug_lst:
                    if drug not in drug_genes:
                        drug_genes[drug] = []
                    if gene_ID not in drug_genes[drug]:
                        drug_genes[drug].append(gene_ID)

                # Initiate empty dict to store relevant mutation information
                mut_info = dict()

                # Save need mutation info with pmid cooresponding to the amino
                # acid change
                for i in range(len(alt_aa)):
                    try:
                        mut_info[alt_aa[i]] = {"gene_name": gene_name,
                                               "drug": res_drug,
                                               "pmid": pmid[i]}
                    except IndexError:
                        mut_info[alt_aa[i]] = {"gene_name": gene_name,
                                               "drug": res_drug,
                                               "pmid": "-"}

                # Add all possible types of mutations to the dict
                if gene_ID not in known_mutations:
                    known_mutations[gene_ID] = {"sub": dict(), "ins": dict(),
                                                "del": dict()}
                # Check for the type of mutation
                if indelflag is False:
                    mutation_type = "sub"
                else:
                    mutation_type = ref_aa

                # Save mutations positions with required information given in
                # mut_info
                if mut_pos not in known_mutations[gene_ID][mutation_type]:
                    known_mutations[gene_ID][mutation_type][mut_pos] = dict()
                for amino in alt_aa:
                    if (amino in known_mutations[gene_ID][mutation_type]
                                                [mut_pos]):
                        stored_mut_info = (known_mutations[gene_ID]
                                                          [mutation_type]
                                                          [mut_pos][amino])
                        if stored_mut_info["drug"] != mut_info[amino]["drug"]:
                            stored_mut_info["drug"] += "," + (mut_info[amino]
                                                                      ["drug"])
                        if stored_mut_info["pmid"] != mut_info[amino]["pmid"]:
                            stored_mut_info["pmid"] += ";" + (mut_info[amino]
                                                                      ["pmid"])

                        (known_mutations[gene_ID][mutation_type]
                                        [mut_pos][amino]) = stored_mut_info
                    else:
                        (known_mutations[gene_ID][mutation_type]
                                        [mut_pos][amino]) = mut_info[amino]

        # Check that all genes in the gene list has known mutations
        for gene in gene_list:
            if gene not in known_mutations:
                known_mutations[gene] = {"sub": dict(), "ins": dict(),
                                         "del": dict()}

        return known_mutations, drug_genes, known_stop_codon

    @staticmethod
    def get_db_mutations(mut_db_path, gene_list):
        """
        This function opens the file self.PHENOTYPE_FILE, and reads
        the content into a dict of dicts. The dict will contain
        information about all known mutations given in the database.
        This dict is returned.
        """

        # Initiate variables
        known_mutations = dict()
        known_stop_codon = dict()
        drug_genes = dict()
        indelflag = False
        stopcodonflag = False

        # Go throug mutation file line by line

        with open(mut_db_path, "r") as fh:
            mutdb_file = fh.readlines()
        mutdb_file = [line.strip() for line in mutdb_file]

        for line in mutdb_file:
            # Ignore headers and check where the indel section starts
            if line.startswith("#"):
                if "indel" in line.lower():
                    indelflag = True
                elif "stop codon" in line.lower():
                    stopcodonflag = True
                else:
                    stopcodonflag = False
                continue

            mutation = [data.strip() for data in line.split("\t")]
            gene_ID = mutation[0]
            gene_name = str(mutation[0].split("_")[0])
            # Only consider mutations in genes found in the gene list
            if gene_name in gene_list:
                mut_pos = int(mutation[3])
                ref_codon = mutation[4]                     # Ref_nuc (1 or 3?)
                ref_aa = mutation[5]                        # Ref_codon
                alt_aa = mutation[6].split(",")             # Res_codon
                res_drug = mutation[8].replace("\t", " ")
                pmid = mutation[9].split(",")

                # Add genes associated with drug resistance to drug_genes dict
                drug_lst = res_drug.split(",")
                drug_lst = [d.strip().lower() for d in drug_lst]
                for drug in drug_lst:
                    if drug not in drug_genes:
                        drug_genes[drug] = []
                    if gene_name not in drug_genes[drug]:
                        drug_genes[drug].append(gene_name)

                # Check if stop codons are predictive for resistance
                if stopcodonflag is True:
                    if gene_name not in known_stop_codon:
                        known_stop_codon[gene_name] = {"pos": [],
                                                       "drug": res_drug}
                    known_stop_codon[gene_name]["pos"].append(mut_pos)

                # Initiate empty dict to store relevant mutation information
                mut_info = dict()

                # Save needed mutation info with pmid cooresponding to the
                # amino acid change
                for i in range(len(alt_aa)):
                    try:
                        mut_info[alt_aa[i]] = {"gene_name": gene_name,
                                               "drug": drug_lst,
                                               "pmid": [pmid[i]]}
                    except IndexError:
                        mut_info[alt_aa[i]] = {"gene_name": gene_name,
                                               "drug": drug_lst,
                                               "pmid": "-"}

                # Add all possible types of mutations to the dict
                if gene_name not in known_mutations:
                    known_mutations[gene_name] = {"sub": dict(), "ins": dict(),
                                                  "del": dict()}
                # Check for the type of mutation
                if indelflag is False:
                    mutation_type = "sub"
                else:
                    mutation_type = ref_aa

                # Save mutations positions with required information given in
                # mut_info
                if mut_pos not in known_mutations[gene_name][mutation_type]:
                    known_mutations[gene_name][mutation_type][mut_pos] = dict()

                for amino in alt_aa:
                    stored_info = (known_mutations[gene_name][mutation_type]
                                   [mut_pos].get(amino, False))
                    if stored_info:
                        merged_info = PointFinder.merge_mut_info(stored_info,
                                                                 amino,
                                                                 mut_info)
                        (known_mutations[gene_name][mutation_type]
                                        [mut_pos][amino]) = merged_info
                    else:
                        (known_mutations[gene_name][mutation_type]
                                        [mut_pos][amino]) = mut_info[amino]

        # Check that all genes in the gene list has known mutations
        for gene in gene_list:
            if gene not in known_mutations:
                known_mutations[gene] = {"sub": dict(), "ins": dict(),
                                         "del": dict()}
        return known_mutations, drug_genes, known_stop_codon

    @staticmethod
    def merge_mut_info(stored_info, amino, mut_info):
        for drug in mut_info[amino]["drug"]:
            if drug not in stored_info["drug"]:
                stored_info["drug"].append(drug)

        for pmid in mut_info[amino]["pmid"]:
            if pmid not in stored_info["pmid"]:
                stored_info["pmid"].append(pmid)

        return stored_info

    def find_mismatches(self, gene, sbjct_start, sbjct_seq, qry_seq):
        """
        This function finds mis matches between two sequeces. Depending
        on the the sequence type either the function
        find_codon_mismatches or find_nucleotid_mismatches are called,
        if the sequences contains both a promoter and a coding region
        both functions are called. The function can also call itself if
        alternative overlaps is give. All found mismatches are returned
        """

        # Initiate the mis_matches list that will store all found mis matcehs
        mis_matches = []
        gene_name = gene.split('_')[0]

        # Find mis matches in RNA genes
        if gene_name in self.RNA_gene_list:
            mis_matches += PointFinder.find_nucleotid_mismatches(sbjct_start,
                                                                 sbjct_seq,
                                                                 qry_seq)
        else:
            # Check if the gene sequence is with a promoter
            regex = r"promoter-size-(\d+)(?:bp)"
            promtr_gene_objt = re.search(regex, gene)

            # Check for promoter sequences
            if promtr_gene_objt:
                # Get promoter length
                promtr_len = int(promtr_gene_objt.group(1))

                # Extract promoter sequence, while considering gaps
                # --------agt-->----
                #    ---->?
                if sbjct_start <= promtr_len:
                    # Find position in sbjct sequence where promoter ends
                    promtr_end = 0
                    nuc_count = sbjct_start - 1

                    for i in range(len(sbjct_seq)):
                        promtr_end += 1

                        if sbjct_seq[i] != "-":
                            nuc_count += 1

                        if nuc_count == promtr_len:
                            break

                    # Check if only a part of the promoter is found
                    # --------agt-->----
                    # ----
                    promtr_sbjct_start = -1
                    if nuc_count < promtr_len:
                        promtr_sbjct_start = nuc_count - promtr_len

                    # Get promoter part of subject and query
                    sbjct_promtr_seq = sbjct_seq[:promtr_end]
                    qry_promtr_seq = qry_seq[:promtr_end]

                    # For promoter part find nucleotide mis matches
                    mis_matches += PointFinder.find_nucleotid_mismatches(
                        promtr_sbjct_start, sbjct_promtr_seq, qry_promtr_seq,
                        promoter=True)

                    # Check if gene is also found
                    # --------agt-->----
                    #     -----------
                    if((sbjct_start + len(sbjct_seq.replace("-", "")))
                       > promtr_len):
                        sbjct_gene_seq = sbjct_seq[promtr_end:]
                        qry_gene_seq = qry_seq[promtr_end:]
                        sbjct_gene_start = 1

                        # Find mismatches in gene part
                        mis_matches += self.find_codon_mismatches(
                            sbjct_gene_start, sbjct_gene_seq, qry_gene_seq)

                # No promoter, only gene is found
                # --------agt-->----
                #            -----
                else:
                    sbjct_gene_start = sbjct_start - promtr_len

                    # Find mismatches in gene part
                    mis_matches += self.find_codon_mismatches(
                        sbjct_gene_start, sbjct_seq, qry_seq)

            else:
                # Find mismatches in gene
                mis_matches += self.find_codon_mismatches(
                    sbjct_start, sbjct_seq, qry_seq)

        return mis_matches

    @staticmethod
    def find_nucleotid_mismatches(sbjct_start, sbjct_seq, qry_seq,
                                  promoter=False):
        """
        This function takes two alligned sequence (subject and query),
        and the position on the subject where the alignment starts. The
        sequences are compared one nucleotide at a time. If mis matches
        are found they are saved. If a gap is found the function
        find_nuc_indel is called to find the entire indel and it is
        also saved into the list mis_matches. If promoter sequences are
        given as arguments, these are reversed the and the absolut
        value of the sequence position  used, but when mutations are
        saved the negative value and det reverse sequences are saved in
        mis_mathces.
        """

        # Initiate the mis_matches list that will store all found
        # mismatcehs
        mis_matches = []

        sbjct_start = abs(sbjct_start)
        seq_pos = sbjct_start

        # Set variables depending on promoter status
        factor = 1
        mut_prefix = "r."

        if promoter is True:
            factor = (-1)
            mut_prefix = "n."
            # Reverse promoter sequences
            sbjct_seq = sbjct_seq[::-1]
            qry_seq = qry_seq[::-1]

        # Go through sequences one nucleotide at a time
        shift = 0
        for index in range(sbjct_start - 1, len(sbjct_seq)):
            mut_name = mut_prefix
            mut = ""
            # Shift index according to gaps
            i = index + shift

            # If the end of the sequence is reached, stop
            if i == len(sbjct_seq):
                break

            sbjct_nuc = sbjct_seq[i]
            qry_nuc = qry_seq[i]

            # Check for mis matches
            if sbjct_nuc.upper() != qry_nuc.upper():

                # check for insertions and deletions
                if sbjct_nuc == "-" or qry_nuc == "-":
                    if sbjct_nuc == "-":
                        mut = "ins"
                        indel_start_pos = (seq_pos - 1) * factor + sbjct_start
                        indel_end_pos = seq_pos * factor + sbjct_start
                        indel = PointFinder.find_nuc_indel(sbjct_seq[i:],
                                                           qry_seq[i:])
                    else:
                        mut = "del"
                        indel_start_pos = seq_pos * factor + (sbjct_start - 1)
                        indel = PointFinder.find_nuc_indel(qry_seq[i:],
                                                           sbjct_seq[i:])
                        indel_end_pos = (seq_pos + len(indel) - 1) * factor \
                                        + (sbjct_start - 1)
                        seq_pos += len(indel) - 1

                    # Shift the index to the end of the indel
                    shift += len(indel) - 1

                    # Write mutation name, depending on sequnce
                    if len(indel) == 1 and mut == "del":
                        mut_name += str(indel_start_pos) + mut + indel
                    else:
                        if promoter is True:
                            # Reverse the sequence and the start and
                            # end positions
                            indel = indel[::-1]
                            temp = indel_start_pos
                            indel_start_pos = indel_end_pos
                            indel_end_pos = temp

                        mut_name += (str(indel_start_pos) + "_"
                                     + str(indel_end_pos) + mut + indel)

                    mis_matches += [[mut, indel_start_pos, indel_end_pos,
                                    indel, mut_name, mut, indel]]

                # Check for substitutions mutations
                else:
                    mut = "sub"
                    pos = seq_pos * factor + (sbjct_start - 1)
                    mut_name += (str(pos)
                                 + sbjct_nuc + ">" + qry_nuc)

                    mis_matches += [[mut, pos, pos,
                                    qry_nuc, mut_name, sbjct_nuc, qry_nuc]]

            # Increment sequence position
            if mut != "ins":
                seq_pos += 1

        return mis_matches

    @staticmethod
    def find_nuc_indel(gapped_seq, indel_seq):
        """
        This function finds the entire indel missing in from a gapped
        sequence compared to the indel_seqeunce. It is assumes that the
        sequences start with the first position of the gap.
        """
        ref_indel = indel_seq[0]
        for j in range(1, len(gapped_seq)):
            if gapped_seq[j] == "-":
                ref_indel += indel_seq[j]
            else:
                break
        return ref_indel

    @staticmethod
    def aa(codon):
        """
        This function converts a codon to an amino acid. If the codon
        is not valid an error message is given, or else, the amino acid
        is returned.

        Potential future issue: If species are added that utilizes
                                alternative translation tables.
        """
        codon = codon.upper()
        aa = {"ATT": "I", "ATC": "I", "ATA": "I",
              "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "TTA": "L",
              "TTG": "L",
              "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
              "TTT": "F", "TTC": "F",
              "ATG": "M",
              "TGT": "C", "TGC": "C",
              "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
              "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
              "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
              "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
              "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S",
              "AGC": "S",
              "TAT": "Y", "TAC": "Y",
              "TGG": "W",
              "CAA": "Q", "CAG": "Q",
              "AAT": "N", "AAC": "N",
              "CAT": "H", "CAC": "H",
              "GAA": "E", "GAG": "E",
              "GAT": "D", "GAC": "D",
              "AAA": "K", "AAG": "K",
              "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R",
              "AGG": "R",
              "TAA": "*", "TAG": "*", "TGA": "*"}
        # Translate valid codon
        try:
            amino_a = aa[codon]
        except KeyError:
            amino_a = "?"
        return amino_a

    @staticmethod
    def get_codon(seq, codon_no, start_offset):
        """
        This function takes a sequece and a codon number and returns
        the codon found in the sequence at that position
        """
        seq = seq.replace("-", "")
        codon_start_pos = int(codon_no - 1) * 3 - start_offset
        codon = seq[codon_start_pos:codon_start_pos + 3]
        return codon

    @staticmethod
    def name_insertion(sbjct_seq, codon_no, sbjct_nucs, aa_alt, start_offset):
        """
        This function is used to name a insertion mutation based on the
        HGVS recommendation.
        """
        start_codon_no = codon_no - 1

        if len(sbjct_nucs) == 3:
            start_codon_no = codon_no

        start_codon = PointFinder.get_codon(sbjct_seq, start_codon_no,
                                            start_offset)
        end_codon = PointFinder.get_codon(sbjct_seq, codon_no, start_offset)
        pos_name = "p.%s%d_%s%dins%s" % (PointFinder.aa(start_codon),
                                         start_codon_no,
                                         PointFinder.aa(end_codon),
                                         codon_no, aa_alt)
        return pos_name

    @staticmethod
    def name_deletion(sbjct_seq, sbjct_rf_indel, sbjct_nucs, codon_no, aa_alt,
                      start_offset, mutation="del"):
        """
        This function is used to name a deletion mutation based on the
        HGVS recommendation. If combination of a deletion and an
        insertion is identified the argument 'mutation' is set to
        'delins' and the mutation name will indicate that the mutation
        is a delins mutation.
        """
        del_codon = PointFinder.get_codon(sbjct_seq, codon_no, start_offset)
        pos_name = "p.%s%d" % (PointFinder.aa(del_codon), codon_no)

        # This has been changed
        if len(sbjct_rf_indel) == 3 and mutation == "del":
            return pos_name + mutation

        end_codon_no = codon_no + math.ceil(len(sbjct_nucs) / 3) - 1
        end_codon = PointFinder.get_codon(sbjct_seq, end_codon_no,
                                          start_offset)
        pos_name += "_%s%d%s" % (PointFinder.aa(end_codon), end_codon_no,
                                 mutation)
        if mutation == "delins":
            pos_name += aa_alt
        return pos_name

    @staticmethod
    def name_indel_mutation(sbjct_seq, indel, sbjct_rf_indel, qry_rf_indel,
                            codon_no, mut, start_offset):
        """
        This function serves to name the individual mutations
        dependently on the type of the mutation.
        """
        # Get the subject and query sequences without gaps
        sbjct_nucs = sbjct_rf_indel.replace("-", "")
        qry_nucs = qry_rf_indel.replace("-", "")

        # Translate nucleotides to amino acids
        aa_ref = ""
        aa_alt = ""

        for i in range(0, len(sbjct_nucs), 3):
            aa_ref += PointFinder.aa(sbjct_nucs[i:i + 3])

        for i in range(0, len(qry_nucs), 3):
            aa_alt += PointFinder.aa(qry_nucs[i:i + 3])

        # Identify the gapped sequence
        if mut == "ins":
            gapped_seq = sbjct_rf_indel
        else:
            gapped_seq = qry_rf_indel

        gap_size = gapped_seq.count("-")

        # Write mutation names
        if gap_size < 3 and len(sbjct_nucs) == 3 and len(qry_nucs) == 3:
            # Write mutation name for substitution mutation
            mut_name = "p.%s%d%s" % (PointFinder.aa(sbjct_nucs), codon_no,
                                     PointFinder.aa(qry_nucs))
        elif len(gapped_seq) == gap_size:

            if mut == "ins":
                # Write mutation name for insertion mutation
                mut_name = PointFinder.name_insertion(sbjct_seq, codon_no,
                                                      sbjct_nucs, aa_alt,
                                                      start_offset)
                aa_ref = mut
            else:
                # Write mutation name for deletion mutation
                mut_name = PointFinder.name_deletion(sbjct_seq, sbjct_rf_indel,
                                                     sbjct_nucs, codon_no,
                                                     aa_alt, start_offset,
                                                     mutation="del")
                aa_alt = mut
        # Check for delins - mix of insertion and deletion
        else:
            # Write mutation name for a mixed insertion and deletion
            # mutation
            mut_name = PointFinder.name_deletion(sbjct_seq,
                                                 sbjct_rf_indel, sbjct_nucs,
                                                 codon_no, aa_alt,
                                                 start_offset,
                                                 mutation="delins")

        # Check for frameshift
        if gapped_seq.count("-") % 3 != 0:
            # Add the frameshift tag to mutation name
            mut_name += " - Frameshift"

        return mut_name, aa_ref, aa_alt

    @staticmethod
    def get_inframe_gap(seq, nucs_needed=3):
        """
        This funtion takes a sequnece starting with a gap or the
        complementary seqeuence to the gap, and the number of
        nucleotides that the seqeunce should contain in order to
        maintain the correct reading frame. The sequence is gone
        through and the number of non-gap characters are counted. When
        the number has reach the number of needed nucleotides the indel
        is returned. If the indel is a 'clean' insert or deletion that
        starts in the start of a codon and can be divided by 3, then
        only the gap is returned.
        """
        nuc_count = 0
        gap_indel = ""
        nucs = ""

        for i in range(len(seq)):
            # Check if the character is not a gap
            if seq[i] != "-":
                # Check if the indel is a 'clean'
                # i.e. if the insert or deletion starts at the first
                # nucleotide in the codon and can be divided by 3

                if(gap_indel.count("-") == len(gap_indel)
                   and gap_indel.count("-") >= 3 and len(gap_indel) != 0):
                    return gap_indel

                nuc_count += 1
            gap_indel += seq[i]

            # if the number of nucleotides in the indel equals the amount
            # needed for the indel, the indel is returned.
            if nuc_count == nucs_needed:
                return gap_indel

        # This will only happen if the gap is in the very end of a sequence
        return gap_indel

    @staticmethod
    def get_indels(sbjct_seq, qry_seq, start_pos):
        """
        This function uses regex to find inserts and deletions in
        sequences given as arguments. A list of these indels are
        returned. The list includes, type of mutations(ins/del),
        subject codon no of found mutation, subject sequence position,
        insert/deletions nucleotide sequence, and the affected qry
        codon no.
        """

        seqs = [sbjct_seq, qry_seq]
        indels = []
        gap_obj = re.compile(r"-+")
        for i in range(len(seqs)):
            for match in gap_obj.finditer(seqs[i]):
                pos = int(match.start())
                gap = match.group()

                # Find position of the mutation corresponding to the
                # subject sequence
                sbj_pos = len(sbjct_seq[:pos].replace("-", "")) + start_pos

                # Get indel sequence and the affected sequences in
                # sbjct and qry in the reading frame
                indel = seqs[abs(i - 1)][pos:pos + len(gap)]

                # Find codon number for mutation
                codon_no = int(math.ceil((sbj_pos) / 3))
                qry_pos = len(qry_seq[:pos].replace("-", "")) + start_pos
                qry_codon = int(math.ceil((qry_pos) / 3))

                if i == 0:
                    mut = "ins"
                else:
                    mut = "del"

                indels.append([mut, codon_no, sbj_pos, indel, qry_codon])

        # Sort indels based on codon position and sequence position
        indels = sorted(indels, key=lambda x: (x[1], x[2]))

        return indels

    @staticmethod
    def _remove_insertions_in_seq(sbjct_seq, qry_seq):
        """
        Remove insertions in query sequence (ex.: assumed to be sequencing
        errors).
        """
        sbjct_seq_list = []
        qry_seq_list = []

        for i, char in enumerate(sbjct_seq):
            if char != '-':
                sbjct_seq_list.append(char)
                qry_seq_list.append(qry_seq[i])

        return ("".join(sbjct_seq_list), "".join(qry_seq_list))

    def find_codon_mismatches(self, sbjct_start, sbjct_seq, qry_seq):
        """
        This function takes two alligned sequence (subject and query),
        and the position on the subject where the alignment starts. The
        sequences are compared codon by codon. If a mismatch is
        found it is saved in 'mis_matches'. If a gap is found the
        function get_inframe_gap is used to find the indel sequence and
        keep the sequence in the correct reading frame. The function
        translate_indel is used to name indel mutations and translate
        the indels to amino acids. The function returns a list of tuples
        containing all needed information about the mutation in order
        to look it up in the database dict known mutation and the with
        the output files the the user.
        """
        mis_matches = []

        # Find start pos of first codon in frame, i_start
        codon_offset = (sbjct_start - 1) % 3
        i_start = 0

        # if codon_offset != 0:
        #     i_start = 3 - codon_offset
        if codon_offset != 0:
            i_start = (int(sbjct_start / 3) + 1) * 3

        # sbjct_start = sbjct_start + i_start

        # Set sequences in frame
        sbjct_seq = sbjct_seq[i_start:]
        qry_seq = qry_seq[i_start:]

        if self.ignore_indels:
            sbjct_seq, qry_seq = PointFinder._remove_insertions_in_seq(
                sbjct_seq, qry_seq)

        # Find codon number of the first codon in the sequence, start
        # at 0
        # codon_no = int((sbjct_start - 1) / 3)  # 1,2,3 start on 0
        if sbjct_start == 1:
            codon_no = 0
        else:
            codon_no = (int(sbjct_start / 3) + 1)

        # s_shift and q_shift are used when gaps appears
        q_shift = 0
        s_shift = 0
        mut_no = 0

        if not self.ignore_indels:
            # Find inserts and deletions in sequence
            indel_no = 0
            indels = PointFinder.get_indels(sbjct_seq, qry_seq, sbjct_start)

        # Go through sequence and save mutations when found
        for index in range(0, len(sbjct_seq), 3):
            # Count codon number
            codon_no += 1

            # Shift index according to gaps
            s_i = index + s_shift
            q_i = index + q_shift

            # Get codons
            sbjct_codon = sbjct_seq[s_i:s_i + 3]
            qry_codon = qry_seq[q_i:q_i + 3]

            if(len(sbjct_seq[s_i:].replace("-", ""))
               + len(qry_codon[q_i:].replace("-", "")) < 6):
                break

            # Check for mutations
            if sbjct_codon.upper() != qry_codon.upper():

                if(("-" in sbjct_codon or "-" in qry_codon)
                   and not self.ignore_indels):
                    # Check for codon insertions and deletions and
                    # frameshift mutations

                    # Get indel info
                    try:
                        indel_data = indels[indel_no]
                    except IndexError:
                        sys.exit("indel_data list is out of range, bug!")

                    mut = indel_data[0]
                    codon_no_indel = indel_data[1]
                    seq_pos = indel_data[2]
                    indel = indel_data[3]
                    indel_no += 1

                    # Get the affected sequence in frame for both for
                    # sbjct and qry
                    if mut == "ins":
                        sbjct_rf_indel = PointFinder.get_inframe_gap(
                            sbjct_seq[s_i:], 3)
                        qry_rf_indel = PointFinder.get_inframe_gap(
                            qry_seq[q_i:],
                            int(math.floor(len(sbjct_rf_indel) / 3) * 3))
                    else:
                        qry_rf_indel = PointFinder.get_inframe_gap(
                            qry_seq[q_i:], 3)
                        sbjct_rf_indel = PointFinder.get_inframe_gap(
                            sbjct_seq[s_i:],
                            int(math.floor(len(qry_rf_indel) / 3) * 3))

                    mut_name, aa_ref, aa_alt = PointFinder.name_indel_mutation(
                        sbjct_seq, indel, sbjct_rf_indel, qry_rf_indel,
                        codon_no, mut, sbjct_start - 1)

                    # Set index to the correct reading frame after the
                    # indel gap
                    shift_diff_before = abs(s_shift - q_shift)
                    s_shift += len(sbjct_rf_indel) - 3
                    q_shift += len(qry_rf_indel) - 3
                    shift_diff = abs(s_shift - q_shift)

                    if shift_diff_before != 0 and shift_diff % 3 == 0:

                        if s_shift > q_shift:
                            nucs_needed = (int((len(sbjct_rf_indel) / 3) * 3)
                                           + shift_diff)
                            pre_qry_indel = qry_rf_indel
                            qry_rf_indel = PointFinder.get_inframe_gap(
                                qry_seq[q_i:], nucs_needed)
                            q_shift += len(qry_rf_indel) - len(pre_qry_indel)

                        elif q_shift > s_shift:
                            nucs_needed = (int((len(qry_rf_indel) / 3) * 3)
                                           + shift_diff)
                            pre_sbjct_indel = sbjct_rf_indel
                            sbjct_rf_indel = PointFinder.get_inframe_gap(
                                sbjct_seq[s_i:], nucs_needed)
                            s_shift += (len(sbjct_rf_indel)
                                        - len(pre_sbjct_indel))

                        mut_name, aa_ref, aa_alt = (
                            PointFinder.name_indel_mutation(
                                sbjct_seq, indel, sbjct_rf_indel, qry_rf_indel,
                                codon_no, mut, sbjct_start - 1)
                        )
                        if "Frameshift" in mut_name:
                            mut_name = (mut_name.split("-")[0]
                                        + "- Frame restored")
                        if mut_name == "p.V940delins - Frame restored":
                            sys.exit()
                    mis_matches += [[mut, codon_no_indel, seq_pos, indel,
                                     mut_name, sbjct_rf_indel, qry_rf_indel,
                                     aa_ref, aa_alt]]

                    # Check if the next mutation in the indels list is
                    # in the current codon.
                    # Find the number of individul gaps in the
                    # evaluated sequence.
                    no_of_indels = (len(re.findall(r"\-\w", sbjct_rf_indel))
                                    + len(re.findall(r"\-\w", qry_rf_indel)))

                    if no_of_indels > 1:

                        for j in range(indel_no, indel_no + no_of_indels - 1):
                            try:
                                indel_data = indels[j]
                            except IndexError:
                                sys.exit("indel_data list is out of range, "
                                         "bug!")
                            mut = indel_data[0]
                            codon_no_indel = indel_data[1]
                            seq_pos = indel_data[2] + sbjct_start - 1
                            indel = indel_data[3]
                            indel_no += 1
                            mis_matches += [[mut, codon_no_indel, seq_pos,
                                             indel, mut_name, sbjct_rf_indel,
                                             qry_rf_indel, aa_ref, aa_alt]]

                    # Set codon number, and save nucleotides from out
                    # of frame mutations
                    if mut == "del":
                        codon_no += int((len(sbjct_rf_indel) - 3) / 3)
                    # If evaluated insert is only gaps codon_no should
                    # not increment
                    elif sbjct_rf_indel.count("-") == len(sbjct_rf_indel):
                        codon_no -= 1

                # Check of point mutations
                else:
                    mut = "sub"
                    aa_ref = PointFinder.aa(sbjct_codon)
                    aa_alt = PointFinder.aa(qry_codon)

                    if aa_ref != aa_alt:
                        # End search for mutation if a premature stop
                        # codon is found
                        mut_name = "p." + aa_ref + str(codon_no) + aa_alt

                        mis_matches += [[mut, codon_no, codon_no, aa_alt,
                                         mut_name, sbjct_codon, qry_codon,
                                         aa_ref, aa_alt]]
                # If a Premature stop codon occur report it an stop the
                # loop
                if not self.ignore_stop_codons:
                    try:
                        if mis_matches[-1][-1] == "*":
                            mut_name += " - Premature stop codon"
                            mis_matches[-1][4] = (
                                mis_matches[-1][4].split("-")[0]
                                + " - Premature stop codon")
                            break
                    except IndexError:
                        pass

        # Sort mutations on position
        mis_matches = sorted(mis_matches, key=lambda x: x[1])

        return mis_matches

    @staticmethod
    def get_file_content(file_path, fst_char_only=False):
        """
        This function opens a file, given as the first argument
        file_path and returns the lines of the file in a list or the
        first character of the file if fst_char_only is set to True.
        """
        with open(file_path, "r") as infile:
            line_lst = []
            for line in infile:
                line = line.strip()
                if line != "":
                    line_lst.append(line)
                if fst_char_only:
                    return line_lst[0][0]
        return line_lst
