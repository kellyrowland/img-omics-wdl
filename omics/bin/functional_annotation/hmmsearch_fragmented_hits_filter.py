#!/usr/bin/env python3

import argparse
import fileinput
import sys


parser = argparse.ArgumentParser(description="""This script expects to get
columns 1, 3, 4, 5, 6, 7, 8, 13, 14, 16, 17, 20, 21 (gene_id, gene_length,
model, accession, hmm_length, full_seq_evalue, full_seq_bitscore, evalue,
bitscore, hmm_start, hmm_end, env_start, env_end) from a domtblout output
sorted by the first column, full_seq_bitscore and full_seq_evalue
('sort -k1,1 -k7,7nr -k6,6n').
For each gene it then groups the hits to the same model (fragments) and
calculates the cumulative alignment length. If the cumulative alignment length
is not covering at least 70% of the shorter of the gene or model length the
fragmented hit gets ignored. If two fragmented hits overlap for more than 10%
the fragmented hit with the lower full sequence bitscore gets removed.
The remaining hits/fragments get printed to stdout in GFF format.""",
formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("hmmer_version", metavar="hmmer_version",
                    help="hmmer version used to create the output")
parser.add_argument("-a", "--aln_len_ratio", nargs="?",
                    metavar="min_aln_len_ratio", dest="min_aln_len_ratio",
                    type=float, default=0.7,
                    help="""the minimum required alignment length ratio on the
                    shorter gene/model""")
parser.add_argument("-o", "--overlap_ratio_cutoff", nargs="?",
                    metavar="max_overlap_ratio_cutoff",
                    dest="max_overlap_ratio_cutoff", type=float,
                    default=0.10, help="""the maximum allowed overlap ratio on
                    the shorter hit""")
parser.add_argument("-v", "--version", action="version",
                    version="%(prog)s 1.0.0")
args, input = parser.parse_known_args()



class Hit:
    def __init__(self, hit_fields):
        self.gene_name = hit_fields[0]
        self.gene_length = int(hit_fields[1])
        self.model = hit_fields[2]
        self.accession = hit_fields[3]
        if self.accession == "-":
            if self.model.startswith("SM"): # SMART
                self.accession = self.model
            else:                           # SuperFamily
                self.accession = hit_fields[2][:hit_fields[2].find("/")]
        self.model_length = int(hit_fields[4])
        self.full_seq_evalue = hit_fields[5]
        self.full_seq_bitscore = float(hit_fields[6])
        self.domain_evalue = hit_fields[7]
        self.domain_bitscore = hit_fields[8]
        self.model_start = int(hit_fields[9])
        self.model_end = int(hit_fields[10])
        self.fake_percent_id = ((self.model_end - self.model_start + 1) /
                                self.model_length) * 100
        self.gene_start = int(hit_fields[11])
        self.gene_end = int(hit_fields[12])
        self.alignment_length = self.gene_end - self.gene_start + 1


    def to_gff(self):
        print(self.gene_name + "\t" +
              args.hmmer_version + "\t" +
              self.accession + "\t" +
              str(self.gene_start) + "\t" +
              str(self.gene_end) + "\t" +
              self.domain_bitscore + "\t" +
              ".\t.\t" +                       # Strand and Phase
              "ID=" + self.gene_name + "_" +   # Start of attributes field
              str(self.gene_start) + "_" +
              str(self.gene_end) +
              ";fake_percent_id=" +
              str("%.2f" % (((self.model_end - self.model_start + 1) /
                            self.model_length) * 100)) +
              ";alignment_length=" + str(self.alignment_length) +
              ";independent_domain_e-value=" + self.domain_evalue +
              ";full_sequence_e-value=" + self.full_seq_evalue +
              ";full_sequence_bitscore=" + str(self.full_seq_bitscore) +
              ";model_start=" + str(self.model_start) +
              ";model_end=" + str(self.model_end))



class Fragmented_Hit:

    def __init__(self, hit):
        self.fragments = [hit]


    def add_hit_fragment(self, hit):
        self.fragments.append(hit)


    def calculate_cumulative_length(self):
        self.cum_aln_length = 0
        for fragment in self.fragments:
            self.cum_aln_length += fragment.alignment_length


    def get_cumulative_length(self):
        if not hasattr(self, 'cum_aln_length'):
            self.calculate_cumulative_length()
        return self.cum_aln_length


    def get_full_seq_bitscore(self):
        if not hasattr(self, 'full_Seq_bitscore'):
            self.full_seq_bitscore = self.fragments[0].full_seq_bitscore
        return self.full_seq_bitscore


    def sort_and_verify_fragments_have_same_full_seq_bitscore(self):
        self.fragments.sort(key=lambda x: x.gene_start) # Sort by gene start
        full_seq_bitscore = self.get_full_seq_bitscore()
        for i in range(1, len(self.fragments)):
            diff = full_seq_bitscore - self.fragments[i].full_seq_bitscore
            if abs(diff) > 0:
                print("Full seq bitscores for " + self.fragments[i].model +
                      " fragments hit by gene " + self.fragments[i].gene_name +
                      " are not identical (" + str(full_seq_bitscore) + " - " +
                      str(self.fragments[i].full_seq_bitscore) + ")",
                      file = sys.stderr)
                sys.exit(1)


    def output_fragments_as_gff(self):
        for fragment in self.fragments:
            fragment.to_gff()



class Fragmented_Hit_Filter:

    def __init__(self):
        self.fragmented_hits = {}


    def add_hit(self, hit):
        if hit.model not in self.fragmented_hits:
            fragmented_hit = Fragmented_Hit(hit)
            self.fragmented_hits[hit.model] = fragmented_hit
        else:
            self.fragmented_hits[hit.model].add_hit_fragment(hit)


    def get_overlap_length(self, fragment_start, fragment_end,
                         other_fragment_start, other_fragment_end):
        # First check if the overlaps overlap at all
        if fragment_start > other_fragment_end or fragment_end < other_fragment_start:
            return 0

        overlap_start = max(fragment_start, other_fragment_start)
        overlap_end = min(fragment_end, other_fragment_end)
        return overlap_end - overlap_start + 1


    def is_significantly_overlapping(self, fragment, other_fragment,
                                     overlap_size):
        overlap_ratio = overlap_size / min(fragment.get_cumulative_length(),
                                           other_fragment.get_cumulative_length())

        if overlap_ratio > args.max_overlap_ratio_cutoff:
            return True

        return False


    def get_conflicting_fragment(self, fragmented_hit, other_fragmented_hit):
        overlap_size = 0
        for fragment in fragmented_hit.fragments:
            for other_fragment in other_fragmented_hit.fragments:
                overlap_size += self.get_overlap_length(fragment.gene_start, fragment.gene_end,
                                                        other_fragment.gene_start, other_fragment.gene_end)

        if self.is_significantly_overlapping(fragmented_hit,
                                             other_fragmented_hit,
                                             overlap_size):
            if (fragmented_hit.get_full_seq_bitscore() >=
                other_fragmented_hit.get_full_seq_bitscore()):
                # Maybe make same bitscore a separate clause
                return -1
            else:
                return 1

        return 0


    def eliminate_overlaps(self):
        ignore = set()
        models = list(self.fragmented_hits.keys())
        for i in range(len(models)-1):
            if models[i] in ignore:
                continue
            fragmented_hit = self.fragmented_hits[models[i]]
            for j in range(i+1, len(models)):
                if models[j] in ignore:
                    continue
                other_fragmented_hit = self.fragmented_hits[models[j]]
                fragment_to_ignore = self.get_conflicting_fragment(fragmented_hit,
                                                              other_fragmented_hit)
                if fragment_to_ignore == 1:
                    ignore.add(models[i])
                    break
                elif fragment_to_ignore == -1:
                    ignore.add(models[j])
                # else means no conflicts
        for model in ignore:
            del self.fragmented_hits[model]


    def remove_too_short_aligned_fragments(self):
        remove = []
        for model in self.fragmented_hits.keys():
            fragmented_hit = self.fragmented_hits[model]
            gene_length = fragmented_hit.fragments[0].gene_length
            model_length = fragmented_hit.fragments[0].model_length
            if ((fragmented_hit.get_cumulative_length() /
                 min(gene_length, model_length)) < args.min_aln_len_ratio):
                remove.append(model)
        for model in remove:
            del self.fragmented_hits[model]


    def sort_each_fragment_and_check_its_validity(self):
        for model in self.fragmented_hits.keys():
            self.fragmented_hits[model].sort_and_verify_fragments_have_same_full_seq_bitscore()


    def output_hits_as_gff(self):
        for model in self.fragmented_hits:
            self.fragmented_hits[model].output_fragments_as_gff()


    def process_and_print_out_final_hits(self):
        self.sort_each_fragment_and_check_its_validity()
        self.remove_too_short_aligned_fragments()
        self.eliminate_overlaps()
        self.output_hits_as_gff()



previous_gene = ""
fragmented_hit_filter = Fragmented_Hit_Filter()
# Run over stdin
for line in fileinput.input(input):
    line = line.rstrip()
    fields = line.split()
    hit = Hit(fields)

    if hit.gene_name != previous_gene:
        # Now we look at the hits for a new gene and thus
        # the stored (fragmented) hits for the previous
        # genes need to be processed and printed out.
        fragmented_hit_filter.process_and_print_out_final_hits()
        previous_gene = hit.gene_name
        fragmented_hit_filter = Fragmented_Hit_Filter()

    fragmented_hit_filter.add_hit(hit)
fragmented_hit_filter.process_and_print_out_final_hits()

