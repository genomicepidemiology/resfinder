#!/usr/bin/env python3

import sys
from collections import defaultdict
from Bio import SeqIO
from resfinder.cge.output.std_results import ResFinderResultHandler, \
    PointFinderResultHandler
from resfinder.cge.resfinder import ResFinder
from resfinder.cge.pointfinder import PointFinder
from cgelib.alignment.read_alignment import KMAAlignment, BlastNAlignment
from cgecore.blaster import Blaster


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


class ResultHandler():

    def __init__(self, conf, pheno_db):
        self.conf = conf
        self.pheno_db = pheno_db

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
            best_blast_hits = self.find_best_Blast_hit(blast_alignment, db_name)
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

            best_blast_hits = self.find_best_Blast_hit(blast_alignment, db_name)

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
                            < self.conf.pf_gene_cov * 100
                            or float(
                                kmahit['template_identity']) <
                            self.conf.pf_gene_id * 100
                            or float(kmahit['depth']) < self.conf.min_depth):
                        continue

                    PointFinderResultHandler.standardize_results_new(std_result,
                                                                     kmahit,
                                                                     db_name,
                                                                     finder,
                                                                     self.pheno_db,
                                                                     self.conf)
            else:  # method == Resfinder.TYPE_KMA and (db_name 'ResFinder or db_name 'DisinFinder):
                for kmahit in kma_alignment.parse_hits():
                    if (float(kmahit['template_coverage'])
                            < self.conf.rf_gene_cov * 100
                            or float(kmahit['template_identity'])
                            < self.conf.rf_gene_id * 100
                            or float(kmahit['depth']) < self.conf.min_depth):
                        continue

                    ResFinderResultHandler.standardize_results_new(std_result,
                                                                   kmahit,
                                                                   db_name,
                                                                   self.conf)

    # todo: ask alfred if this should be moved to cgelib as a blast filtering step.
    def find_best_Blast_hit(self, aligner, db_name):
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
                        "coverage": hit['n_identity'] / hit[
                            "gene_length_XMLFile_undescribed"],
                        # "strand": strand,
                        "gene_length": hit["gene_length_XMLFile_undescribed"],
                        "hit_id": hit_id}

            # if gene is already added check for overlap and keep the
            # new combined hit
            if gene in gene_dict.keys():
                combined_gene = self.gene_overlap_comparison(
                    gene_hit, gene_dict[gene][0])

                gene_dict[gene].append(combined_gene)
                gene_dict[gene].pop(0)

            # check if the contig is matching another gene
            keys = [k for k, v in gene_dict.items() if v[0]['contig_name']
                    == gene_hit['contig_name']]
            if keys:
                keys_to_keep, old_keys_to_drop = self.keep_hit(gene_dict,
                                                               gene_hit,
                                                               gene, keys)
                # if gene is added to old_keys_to_drop it can be removed as it
                # is not in the gene_dict where the keys will be removed from.
                if gene in old_keys_to_drop:
                    old_keys_to_drop.remove(gene)

                for drop_key in old_keys_to_drop:
                    gene_dict.pop(drop_key)

                if gene in keys_to_keep:
                    gene_dict[gene].append(gene_hit)
                else:
                    continue
            elif gene not in gene_dict.keys():
                gene_dict[gene].append(gene_hit)

            if (len(gene_dict[gene][0]["query_string"])
                    < gene_dict[gene][0]['gene_length']):
                gene_dict[gene][0] = self.complete_template(
                    hit=gene_dict[gene][0],
                    db=db_name)
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
        final_aln = pre_hit['aln_string'].strip()
        final_qry = pre_hit['query_string'].replace('-','')
        first_hit_id = pre_hit['hit_id']
        gene_length = pre_hit['gene_length']

        # alternative_overlaps = []
        contigs = first_hit_id
        contig_name = pre_hit['contig_name']

        pre_start = int(pre_hit['tmpl_start'])
        pre_end = int(pre_hit['tmpl_end'])
        pre_query_start = int(pre_hit['query_start'])
        pre_query_end = int(pre_hit['query_end'])
        pre_qry = pre_hit['query_string'].replace('-','')
        pre_aln = pre_hit['aln_string'].strip()
        pre_tmpl = pre_hit['tmpl_string']
        pre_id = pre_hit['hit_id']

        next_start = int(next_hit['tmpl_start'])
        next_end = int(next_hit['tmpl_end'])
        next_query_start = int(next_hit['query_start'])
        next_query_end = int(next_hit['query_end'])
        next_aln = next_hit['aln_string'].strip()
        next_qry = next_hit['query_string'].replace('-','')
        next_tmpl = next_hit['tmpl_string']
        next_id = next_hit['hit_id']
        next_name = next_hit['contig_name']
        # the full template lenght is given - here next_length is the hit length
        next_length = next_end - next_start + 1

        #getting the full template.
        template = max(pre_tmpl, next_tmpl, key=len)
        final_tmpl = template

        contigs += next_id
        contig_name += (', ' + next_name)
        if next_start <= current_end:
            eprint("Info: {} and {} are aligning to the same gene and have "
                   "been combined to one hit".format(pre_id, next_id))
            #  ---->
            #   <--
            overlap_start = self.find_overlap_start(pre_start,
                                                    template,
                                                    next_start)

            # Find overlap len and add next sequence to final sequence
            if len(template[overlap_start:pre_end]) > next_length:
                #  <--------->
                #     <--->
                pass
            else:
                #  <--------->
                #        <--------->
                overlap_len = len(template[overlap_start:pre_end])
                overlap_end_pos = pre_end

                # Update current end
                current_end = int(next_end)

                # Find query overlap sequences
                pre_qry_overlap = pre_qry[(overlap_start - pre_query_start):
                                          (overlap_start + overlap_len)]
                next_qry_overlap = next_qry[:overlap_len + 1]
                tmpl_overlap = template[overlap_start - 1:(overlap_start + overlap_len)]

                # If alternative query overlap exists use the best
                if pre_qry_overlap != next_qry_overlap:
                    best_overlap = self.get_best_overlap(pre_qry_overlap,
                                          next_qry_overlap,
                                          tmpl_overlap)
                    # Todo: - do we want to print this warning - should we
                    #  identify which hit the used overlap comes from?
                    eprint("OVERLAP WARNING: The following two hits had an "
                           "overlap containing differences - the overlap sequence"
                           " with the highest identity was used.")
                    eprint("{}\n{}"
                           .format(pre_qry_overlap, next_qry_overlap))

                    # add the best overlap to the first sequence followed by
                    # the  next sequence.
                    final_qry = final_qry[:overlap_start] + best_overlap + \
                                next_qry[overlap_end_pos:]
                    final_aln = final_aln[:overlap_start] + best_overlap + \
                                next_aln[overlap_end_pos:]
                else:
                    # Use the entire previous sequence and add the last
                    # part of the next sequence
                    final_qry += next_qry[overlap_len + 1:]
                    final_aln += next_aln[overlap_len + 1:]

        elif next_start > current_end:
            #  <------->
            #              <------->
            gap_size = next_start - current_end - 1
            final_qry += "N" * gap_size
            final_aln += " " * gap_size
            current_end = int(next_end)
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
                         "identity": identity,
                         "hit_id": contigs}

        return combined_gene

    def find_overlap_start(self, pre_start, pre_tmpl, next_start):

        pos_count = 1
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

    def get_best_overlap(self, pre_seq, next_seq, tmpl):

        pre_id = self.calculate_identity(pre_seq, tmpl)
        next_id = self.calculate_identity(next_seq, tmpl)

        if pre_id > next_id:
            return pre_seq
        elif pre_id == next_id:
            return pre_seq
        else:
            return next_seq

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
                                                                 finder,
                                                                 self.pheno_db,
                                                                 self.conf)
            else:
                ResFinderResultHandler.standardize_results_new(std_result,
                                                               hit[0],
                                                               db_name,
                                                               self.conf)

    def keep_hit(self, hit_dict, current_hit, current_key, key_list):
        '''
        input
            hit_dict: dict containing all the previous selcted hits and their
                        information
            current_hit: a dict containing the hit to be compared.
            current_key: key of the hit to be compared
            key_list: a list of all the keys in the hit_dict representing all
            found hits to be compared against.
        output
            keys_to_keep: list of keys to keep. pairwise comparison keeping the
            best hit.
            keys_to_drop: list of keys to drop. another hit was found to be better

        this function will check the current hit against all the other hits in
        the hit_dict looking for overlaps of the contig. it will keep the
        overlapping hit with the best identity or combined coverage/identity
        score.
        Previously a part of Blaster.compare_results()
        '''

        current_contig = current_hit['contig_name']
        current_query_start = current_hit['query_start']
        current_query_end = current_hit['query_end']
        current_cal = (current_hit['identity'] *
                       current_hit['coverage'])
        current_length = len(current_hit['query_string'])
        current_aln_len = current_hit['aln_length']

        keys_to_keep = []
        keys_to_drop = []
        for next_key in key_list:

            next_contig = hit_dict[next_key][0]['contig_name']
            next_query_start = hit_dict[next_key][0]['query_start']
            next_query_end = hit_dict[next_key][0]['query_end']
            next_cal = (hit_dict[next_key][0]['identity'] *
                        hit_dict[next_key][0]['coverage'])
            next_length = len(hit_dict[next_key][0]['query_string'])
            next_aln_len = hit_dict[next_key][0]['aln_length']

            hit_union_length = (max(current_query_end, next_query_end)
                                - min(current_query_start, next_query_start))
            hit_lengths_sum = ((current_query_end - current_query_start)
                               + (next_query_end - next_query_start))
            overlap_len = (hit_lengths_sum - hit_union_length)

            # check if the hits overlap
            if overlap_len < self.conf.rf_overlap:
                keys_to_keep.extend([current_key, next_key])
                continue

            print("contig {} was found to hit both ".format(current_contig))
            print("\n{} and {}".format(current_key, next_key))

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
                    keys_to_drop.append(current_key)
                elif (current_hit['identity']
                      == hit_dict[next_key][0]['identity']):
                    keys_to_keep.extend([current_key, next_key])

            # if the part of the contig do not fully overlap
            elif (hit_union_length <= hit_lengths_sum):
                if current_cal > next_cal:
                    keys_to_keep.append(current_key)
                    keys_to_drop.append(next_key)
                elif current_cal < next_cal:
                    keys_to_keep.append(next_key)
                    keys_to_drop.append(current_key)

                elif current_cal == next_cal:
                    if next_aln_len > current_aln_len:
                        keys_to_keep.append(next_key)
                        keys_to_drop.append(current_key)
                    elif next_aln_len < current_aln_len:
                        keys_to_keep.append(current_key)
                        keys_to_drop.append(next_key)
                    else:
                        # if both coverage and identity (cal_score) is the same,
                        # and the length is the same both hits are kept.
                        keys_to_keep.extend([current_key, next_key])

            keys_to_print = [x for x in [current_key, next_key] if x in keys_to_keep]
            keys_to_print = set(keys_to_print) - set(keys_to_drop)
            print("hit {} was kept".format(sorted(keys_to_print)))
        # TODO  If new_score == old_score but identity and coverages are not the
        #  same. which gene should be chosen?? Keep both for now - possible fix
        #  when introducing a protein level comparison

        # keys_to_keep should not include a key that has been dropped.
        keys_to_keep = set(keys_to_keep) - set(keys_to_drop)

        return set(keys_to_keep), set(keys_to_drop)

    def complete_template(self, hit, db):
        '''
            Input:
                hit - gene_dict[gene][0] object corresponding to a dict with hit
                info
                db - database name
            output:
                hit - a dict containing updated values of sequence, template,
                and alingment
            This function will add the remaining sequence to the template if
            the alignment do not cover the entire gene. It will update the hit
            instance with the template covering the full gene, the query string,
            and the alignment string with spaces in the remaining length of the
            template length.
        '''

        gene_id = hit['hit_id'].split(":")[-1]

        # getting the db path to the fasta file containing the tmpl of the hit.
        db_file = self.get_tmpl_path(db, gene_id)

        # extending the template string to include the entire database entry
        for file in db_file:
            for seq_record in SeqIO.parse(file, "fasta"):
                if seq_record.description.replace(" ", "") == gene_id.replace(" ",""):
                    # start_seq = str(seq_record.seq)[:int(hit["tmpl_start"]) - 1]
                    # end_seq = str(seq_record.seq)[int(hit["tmpl_end"]):]
                    full_tmpl = str(seq_record.seq)
                    break
        # getting the full input contig to extend it to the length of subject
        contig = ''
        for seq_record in SeqIO.parse(self.conf.inputfasta, "fasta"):
            if (seq_record.description.replace(" ", "") == hit[
                "contig_name"].replace(" ", "")):
                contig = str(seq_record.seq)
                break

        # Extract extra sequence from query and generate the respective alignment
        query_seq = self.get_query_align(hit, contig)

        # complete the template string to be the full gene from database
        hit["tmpl_string"] = full_tmpl
        hit["query_string"] = query_seq

        test1 = len(full_tmpl)
        test2 = len(query_seq)
        aln_seq = self.calculate_alignment(query_seq, hit['tmpl_string'])

        hit["aln_string"] = aln_seq

        return hit

    def get_tmpl_path(self, db, gene_id):
        '''
            Input:
                db - database name
                gene_id - header of the gene corresponding to the hit
            output:
                db_file - path to the fasta file containing the database
                template of the hit
            This function will output the path to the fasta file in the
            database that corresponds to the hit found by Blast
        '''
        gene_name, variant, acc = gene_id.split("_")
        db_file = []
        if db == "PointFinder":
            db_path = self.conf.db_path_point
            db_file.append(db_path + "/" + gene_name + ".fsa")
        elif db == "DisinFinder":
            db_path = self.conf.db_path_disinf
            db_file.append(db_path + "/disinfectants.fsa")
        else: # db == ResFinder
            db_path = self.conf.db_path_res
            db_file.append(db_path + "/all.fsa")
        return db_file

    def get_query_align(self, hit, contig):
        '''
        Input:
            hit - an instance of a gene_dict[gene][0] dict.
            contig - the input contig corresponding to the hit
        output:
            updated query string
            updated alignment string
        The function will complete the query and alingment in the case that the
        template gene is longer than the hit sequence. If there is corresponding
        sequences in the contig these will be added otherwise '-' will be
        appended accordingly. The alignment string will include spaces to the
        corresponding added sequence.
        '''
        # getting variables from hit:
        query_seq = hit['query_string']
        aln_seq = hit['aln_string']
        sbjct_start = int(hit['tmpl_start'])
        sbjct_end = int(hit['tmpl_end'])
        query_start = int(hit['query_start'])
        query_end = int(hit['query_end'])
        length = int(hit['gene_length'])

        #if query contains insertions ('-') compared to contig
        # todo there might be more edge-cases where '-' is inserted. what about
        #  deletions?
        nr_insertions = query_seq.count('-')

        # if the subject string do not start form the beginning of the template gene
        if sbjct_start != 1:
            missing = sbjct_start - 1

            # if query is missing the start of the input contig it will be added
            if query_start >= missing:
                start_pos = query_start - missing - 1
                end_pos = query_start - 1
                chars = contig[start_pos:end_pos]
                query_seq = chars + query_seq

            else:  # the query will be extended with "-" in front to reach the gene length
                end_pos = query_start - 1
                chars = contig[0:end_pos]

                query_seq = ("-" * (missing - len(chars)) + chars + query_seq)

        # if the subject end do not reach the end of the template
        if sbjct_end < length:
            missing = length - sbjct_end - nr_insertions
            # if the query has additional sequence on the contig compared to gene length
            if missing <= (len(contig) - query_end):
                start_pos = query_end
                end_pos = query_end + missing
                chars = contig[start_pos:end_pos]
                query_seq = query_seq + chars
            else:  # the missing is longer than what is left in the contig
                start_pos = query_end
                chars = contig[start_pos:]
                query_seq = query_seq + chars + "-" * (missing - len(chars))

        return query_seq

    def calculate_alignment(self, query, temp):
        '''
        input:
            query sequence - string
            template sequence - string
        output:
            alignemnt sequence - string
        This function will calculate the alignment. | for match and space for
        mismatch. This might also be done somewhere else in the code - could
        not find it
        '''
        aln_seq = ''
        for i in range(len(query)):
            if query[i].upper() == temp[i].upper():
                aln_seq += '|'
            else:
                aln_seq += ' '
        return aln_seq
