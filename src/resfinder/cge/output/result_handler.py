#!/usr/bin/env python3

import sys
from collections import defaultdict

from resfinder.cge.output.std_results import ResFinderResultHandler,\
    PointFinderResultHandler
from resfinder.cge.resfinder import ResFinder
from resfinder.cge.pointfinder import PointFinder
from cgelib.alignment.read_alignment import KMAAlignment, BlastNAlignment
from cgecore.blaster import Blaster

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

class ResultHandler():

    def __init__(self, conf):
        self.conf = conf

    def handle_results(self, finder, res, std_result, db_name, method):

        filenames = [x.split('/')[-1] for x in res.alignment_files]
        out_path = \
            [x.rsplit('/', 1)[0] for x in res.alignment_files][0]

        if (method == ResFinder.TYPE_BLAST
            and (db_name == 'ResFinder'
                 or db_name == 'DisinFinder')):

            blast_alignment = BlastNAlignment(output_path=out_path,
                                              filenames=filenames,
                                              result_file="XML")
            best_blast_hits = self.find_best_Blast_hit(blast_alignment,)
            self.filter_and_standardize_result(best_blast_hits,
                                               std_result,
                                               db_name, finder,
                                               self.conf.rf_gene_cov,
                                               self.conf.rf_gene_id)

        elif (method == PointFinder.TYPE_BLAST
                and db_name == 'PointFinder'):

            blast_alignment = BlastNAlignment(output_path=out_path,
                                            filenames=filenames,
                                            result_file="XML")

            best_blast_hits = self.find_best_Blast_hit(blast_alignment)

            self.filter_and_standardize_result(best_blast_hits, std_result,
                                               db_name, finder,
                                               self.conf.pf_gene_cov,
                                               self.conf.pf_gene_id)

        else:
            kma_alignment = KMAAlignment(
                output_path=out_path,
                filenames=filenames,
                result_file=["Result", "Alignment"])

            if (method == PointFinder.TYPE_KMA
                and db_name == 'PointFinder'):
                for kmahit in kma_alignment.parse_hits():

                    if (float(kmahit['template_coverage'])
                            < self.conf.pf_gene_cov*100
                            or float(
                                kmahit['template_identity']) <
                            self.conf.pf_gene_id*100
                            or float(kmahit['depth']) < self.conf.min_depth):
                        continue

                    PointFinderResultHandler.standardize_results_new(std_result,
                                                                     kmahit,
                                                                     db_name,
                                                                     finder)
            else: #method == Resfinder.TYPE_KMA and (db_name 'ResFinder or db_name 'DisinFinder):
                for kmahit in kma_alignment.parse_hits():
                    if (float(kmahit['template_coverage'])
                            < self.conf.rf_gene_cov*100
                            or float(kmahit['template_identity'])
                            < self.conf.rf_gene_id*100
                            or float(kmahit['depth']) < self.conf.min_depth):
                        continue

                    ResFinderResultHandler.standardize_results_new(std_result,
                                                                   kmahit,
                                                                   db_name,
                                                                   method)

    #todo: ask alfred if this should be moved to cgelib as a blast filtering step.
    def find_best_Blast_hit(self, aligner):
        """
        input:
            aligner: blast result in BlastNAlignment object
        output:
            Finds the best sequence hits and returns adjudstet BlastNAlignment
            object, substitute to the find_best_seq in pointfinder.py
        """

        ## overlap detection from pointfinder.find_best_seq()
        # grouping hits within same gene to combine to one hit.
        gene_dict = defaultdict(list)
        for hit in aligner.parse_hits():
            gene_hit = {}
            sbjct_start = ''
            sbjct_end = ''
            query_string = ''
            sbjct_string = ''
            aln_string = ''

            strand = 0
            # If the hit is on the other strand
            if hit['template_start_aln'] > hit['template_end_aln']:
                tmp = hit['template_start_aln']
                sbjct_start = hit['template_end_aln']
                sbjct_end = tmp

                # Todo: move reversecomplement function from cgecore.Blaster
                #  to cgelib.
                query_string = Blaster.reversecomplement(
                    hit['query_aln'])
                aln_string = Blaster.reversecomplement(
                    hit['aln_scheme'])
                sbjct_string = Blaster.reversecomplement(
                    hit['template_aln'])
                strand = 1

            # Check coverage for each hit, patch together partial genes hits
            gene = hit['templateID']
            hit_id = "%s:%s..%s:%s" % (
                                    hit['queryID'], hit['template_start_aln'],
                                    hit['template_end_aln'], gene)



            query = (query_string if query_string
                                    else hit['query_aln'])
            template = (sbjct_string if sbjct_string
                                     else hit['template_aln'])

            identity = self.calculate_identity(query=query,
                                               subject=template)

            gene_hit = {"tmpl_start": sbjct_start if sbjct_start
                                    else hit['template_start_aln'],
                        "tmpl_end": sbjct_end if sbjct_end
                                    else hit['template_end_aln'],
                        "query_start": hit['query_start_aln'],
                        "query_end": hit['query_end_aln'],
                        "query_string": query,
                        "aln_string": aln_string if aln_string
                                else hit['aln_scheme'],
                        "tmpl_string": template,
                        "aln_length": hit['aln_length'],
                        "contig_name": hit['queryID'],
                        "identity": identity,
                        "coverage": hit['n_identity']/hit["gene_length_XMLFile_undescribed"],
                        # "strand": strand,
                        "gene_length": hit["gene_length_XMLFile_undescribed"],
                        "hit_id": hit_id}

            if gene in gene_dict.keys():
                combined_gene = self.gene_overlap_comparison(
                    gene_hit, gene_dict[gene][0])

                gene_dict[gene].append(combined_gene)
                gene_dict[gene].pop(0)

            else:
                keys = [k for k, v in gene_dict.items() if v[0]['contig_name'] == gene_hit['contig_name']]
                if keys:
                    keys_to_keep, old_keys_to_drop = self.keep_hit(gene_dict,
                                                                   gene_hit,
                                                                   gene, keys)
                    for drop_key in old_keys_to_drop:
                        gene_dict.pop(drop_key)

                    if gene in keys_to_keep:
                        gene_dict[gene].append(gene_hit)
                        continue
                    else:
                        continue
                gene_dict[gene].append(gene_hit)
        return gene_dict

    def gene_overlap_comparison(self, pre_hit, next_hit):

        if pre_hit['tmpl_start'] > next_hit['tmpl_start']:
            tmp = pre_hit
            pre_hit = next_hit
            next_hit = tmp

        # get info on first hit (which will be appended with additional hits:
        all_start = int(pre_hit['tmpl_start'])
        current_end = int(pre_hit['tmpl_end'])
        query_start = pre_hit['query_start']
        query_end = pre_hit['query_end']
        final_aln = pre_hit['aln_string']
        final_qry = pre_hit['query_string']
        final_tmpl = pre_hit['tmpl_string']
        first_hit_id = pre_hit['hit_id']
        gene_length = pre_hit['gene_length']

        alternative_overlaps = []
        contigs = first_hit_id
        contig_name = pre_hit['contig_name']

        pre_start = int(pre_hit['tmpl_start'])
        pre_end = int(pre_hit['tmpl_end'])
        pre_query_start = int(pre_hit['query_start'])
        pre__query_end = int(pre_hit['query_end'])
        pre_qry = pre_hit['query_string']
        pre_aln = pre_hit['aln_string']
        pre_tmpl = pre_hit['tmpl_string']
        pre_id = pre_hit['hit_id']

        next_start = int(next_hit['tmpl_start'])
        next_end = int(next_hit['tmpl_end'])
        next_query_start = int(next_hit['query_start'])
        next_query_end = int(next_hit['query_end'])
        next_aln = next_hit['aln_string']
        next_qry = next_hit['query_string']
        next_tmpl = next_hit['tmpl_string']
        next_id = next_hit['hit_id']
        next_name = next_hit['contig_name']

        contigs += next_id
        contig_name += (', ' + next_name)
        if next_start <= current_end:
            eprint("Info: {} and {} are aligning to the same gene and have "
                   "been combined to one hit".format(pre_id, next_id))
            #  ---->
            #   <--
            overlap_start = self.find_overlap_start(pre_start,
                                                    pre_tmpl,
                                                    next_start)

            # Find overlap len and add next sequence to final sequence
            if len(pre_tmpl[overlap_start:]) > len(next_tmpl):
                #  <--------->
                #     <--->
                overlap_len = len(next_tmpl)
                overlap_end_pos = next_end
            else:
                #  <--------->
                #        <--------->
                overlap_len = len(pre_tmpl[overlap_start:])
                overlap_end_pos = pre_end

                # Update current end
                current_end = int(next_end)

                # Use the entire previous sequence and add the last
                # part of the next sequence
                final_tmpl += next_tmpl[overlap_len:]
                final_qry += next_qry[overlap_len:]
                final_aln += next_aln[overlap_len:]

            # Find query overlap sequences
            pre_qry_overlap = pre_qry[overlap_start: (overlap_start
                                                      + overlap_len)]
            next_qry_overlap = next_qry[:overlap_len]
            tmpl_overlap = next_tmpl[:overlap_len]

            # If alternative query overlap exists save it
            #todo this is saved but never used - remove?
            # or figure out which is best and use that? also relates to
            # hsps of identical length - does that happen?
            if pre_qry_overlap != next_qry_overlap:
                eprint("OVERLAP WARNING:")
                eprint("{}\n{}"
                       .format(pre_qry_overlap, next_qry_overlap))

                # Save alternative overlaps
                alternative_overlaps += [(next_start,
                                          overlap_end_pos,
                                          tmpl_overlap,
                                          next_qry_overlap)]

        elif next_start > current_end:
            #  <------->
            #              <------->
            gap_size = next_start - current_end - 1
            final_qry += "N" * gap_size
            final_tmpl += "N" * gap_size
            final_aln += "-" * gap_size
            current_end = int(next_end)
            final_tmpl += next_tmpl
            final_qry += next_qry
            final_aln += next_aln

            eprint("Info: {} and {} are aligning to the same gene and have "
                   "been combined to one hit".format(pre_id, next_id))
            # Calculate new coverage
        no_call = final_qry.upper().count("N")
        # coverage of possible combined alignment length against full gene
        coverage = ((current_end - all_start + 1 - no_call)
                    / float(gene_length))

        identity = self.calculate_identity(query=final_qry,
                                           subject=final_tmpl)

        combined_gene = {"tmpl_start": all_start,
                         "tmpl_end": current_end,
                         "query_start": query_start,
                         "query_end": len(final_qry),
                         "query_string": final_qry,
                         "aln_string": final_aln,
                         "tmpl_string": final_tmpl,
                         "aln_length": (current_end - all_start + 1),
                         "contigs": contigs,
                         "contig_name": contig_name,
                         "coverage": coverage,
                         "gene_length": gene_length,
                         "identity": identity}

        return combined_gene

    def find_overlap_start(self, pre_start, pre_tmpl, next_start):

        pos_count = 0
        overlap_pos = pre_start
        for i in range(len(pre_tmpl)):
            # Stop loop if overlap_start position is reached
            if overlap_pos == next_start:
                overlap_start = pos_count
                break
            if pre_tmpl[i] != "-":
                overlap_pos += 1
            pos_count += 1
        return overlap_start

    def calculate_identity(self, query, subject):
        equal = 0
        not_equal = 0
        for i in range(len(query)):
            if query[i].upper() == subject[i].upper():
                equal += 1
            else:
                not_equal += 1
        return equal / float(equal + not_equal)

    def filter_and_standardize_result(self, combined_hits, std_result,
                                      db_name, finder, min_coverage,
                                      min_identity):
        # todo figure out what to do with alternative overlaps - will we end up
        #  keeping both if ID is the same in two hits? currently assuming only
        #  one item per gene.
        hit_list = list(combined_hits)

        for key, hit in combined_hits.items():
            if (hit[0]['coverage'] < min_coverage
                    or hit[0]['identity'] < min_identity):
                continue

            # adding the combined hit information to the first hit within
            # the correct gene and standardizing this as the result for
            # this gene.
            hit[0]['query_aln'] = hit[0]['query_string']
            hit[0]['aln_scheme'] = hit[0]['aln_string']
            hit[0]['template_aln'] = hit[0]['tmpl_string']
            hit[0]['template_coverage'] = hit[0]['coverage'] * 100
            hit[0]['template_identity'] = hit[0]['identity']
            hit[0]['template_length'] = hit[0]['gene_length']
            hit[0]['ref_start_pos'] = hit[0]['tmpl_start']
            hit[0]['templateID'] = key

            if db_name == 'PointFinder':
                PointFinderResultHandler.standardize_results_new(std_result,
                                                                 hit[0],
                                                                 db_name,
                                                                 finder)
            else:
                ResFinderResultHandler.standardize_results_new(std_result,
                                                               hit[0],
                                                               db_name)

    def keep_hit(self, hit_dict, current_hit, current_key, key_list):
        '''
        input
            hit_dict : dict containing all the different hits and their
                        information
            current_key : key in hit_dict that corresponds to the first hit
                            to be compared
            next_key : key in hit_dict that corresponds to the first hit to
                        be compared
        output
            True or False : defining whether to keep or discard the current hit
                            - if second hit is better.
        this function will check two hits against each other. If thte second hit
        is better the function will return False indicating that the hit should
        not be kept. Previously a part of Blaster.compare_results()
        '''

        current_contig = current_hit['contig_name']
        current_query_start = current_hit['query_start']
        current_query_end = current_hit['query_end']
        current_cal = (current_hit['identity'] *
                       current_hit['coverage'])
        current_length = len(current_hit['query_string'])

        keys_to_keep = []
        keys_to_drop = []
        for next_key in key_list:

            next_contig = hit_dict[next_key][0]['contig_name']
            next_query_start = hit_dict[next_key][0]['query_start']
            next_query_end = hit_dict[next_key][0]['query_end']
            next_cal = (hit_dict[next_key][0]['identity'] *
                        hit_dict[next_key][0]['coverage'])
            next_length = len(hit_dict[next_key][0]['query_string'])

            eprint("contig {} was found to hit both ".format(current_contig))
            eprint("\t{} and {}".format(current_key, next_key))

            hit_union_length = (max(current_query_end, next_query_end)
                                - min(current_query_start, next_query_start))
            hit_lengths_sum = ((current_query_end - current_query_start)
                               + (next_query_end - next_query_start))
            overlap_len = (hit_lengths_sum - hit_union_length)

            if overlap_len < self.conf.rf_overlap:
                # print("\tignore overlap ({}): {}".format(overlap_len, next_key))
                keys_to_keep.append[current_key, next_key]
                continue
            # print("\toverlap found ({}): {}".format(overlap_len, next_key))
            # TODO why do different checks for precise overlap and partial.
            #  wouldn't you assume same coverage for precise overlap.
            #  making cal score diffs be the same as identity diffs.

            # if the two hits overlap precisely, the identity is compared
            if (current_query_end == next_query_end
                and current_query_start == next_query_start):

                if (current_hit['identity']
                        > hit_dict[next_key][0]['identity']):
                    keys_to_keep.append(current_key)
                    keys_to_drop.append(next_key)
                elif (current_hit['identity']
                        < hit_dict[next_key][0]['identity']):
                    keys_to_keep.append(next_key)
                elif (current_hit['identity']
                        == hit_dict[next_key][0]['identity']):
                    keys_to_keep.append(current_key, next_key)

            # if the part of the contig do not fully overlap
            elif (hit_union_length <= hit_lengths_sum):
                if current_cal > next_cal:
                    keys_to_keep.append(current_key)
                    keys_to_drop.append(next_key)
                elif current_cal < next_cal:
                    keys_to_keep.append(next_key)
                elif current_cal == next_cal:

                    if next_length > current_length:
                        keys_to_keep.append(next_key)
                    elif next_length < current_length:
                        keys_to_keep.append(current_key)
                        keys_to_drop.append(next_key)
                    else:
                        #if both coverage and identity (cal_score) is the same,
                        # and the length is the same both hits are kept.
                        keys_to_keep.append(current_key, next_key)
        print("hit {} was kept".format(keys_to_keep))
        # TODO  If new_score == old_score but identity and coverages are not the same.
        #  which gene should be chosen?? Now they are both kept.
        return keys_to_keep, keys_to_drop
