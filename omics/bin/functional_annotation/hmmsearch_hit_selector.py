#!/usr/bin/env python3

import argparse
import fileinput


parser = argparse.ArgumentParser(description="""This script expects to get
columns 1, 3, 4 (or 5 if the accession is preferred over the name), 6, 13, 14,
16, 17, 20, 21 (gene_id, gene_length, model, hmm_length, evalue, bitscore,
hmm_start, hmm_stop, env_start, env_stop) from a domtblout output created by
hmmsearch sorted via 'sort -k1,1 -k6,6nr -k5,5n'.
First the script filters out all hits that do not cover at least 70% (default
value, can be changed via the aln_len_ratio argument) of the shorter one of the
gene or model. It'll then eliminate significant overlaps of hits to the same
gene and print to stdout the remaining hits in GFF format.
When an overlap is considered significant can be influenced by changing
the value of the over_ratio_cutoff argument.
The default is set to 0.10, which would ignore the lower scoring hit of
two hits that overlap by  more than 10% of the length of the shorter of
the two hits.""",
formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("hmmer_version", metavar="hmmer_version",
                    help="hmmer version used to create the output")
parser.add_argument("-a", "--aln_len_ratio", nargs="?",
                    metavar="min_aln_len_ratio", dest="min_aln_len_ratio",
                    type=float, default=0.7, help="""the minimum required
                    alignment length ratio on the shorter gene/model""")
parser.add_argument("-o", "--overlap_ratio", nargs="?",
                    metavar="max_overlap_ratio", dest="max_overlap_ratio",
                    type=float, default=0.10, help="""the maximum allowed
                    overlap ratio on the shorter hit""")
parser.add_argument("-v", "--version", action="version",
                    version="%(prog)s 1.0.0")
args, input = parser.parse_known_args()


def alignment_ratio_is_met(gene_length, model_length, alignment_length):
    if alignment_length >= (min(gene_length, model_length) *
                            args.min_aln_len_ratio):
        return True
    else:
        return False


def get_overlap_length(hit_a_coordinates, hit_b_coordinates):
    # First check if the overlaps overlap at all
    if hit_a_coordinates[0] > hit_b_coordinates[1] or hit_a_coordinates[1] < hit_b_coordinates[0]:
        return 0

    overlap_start = max(hit_a_coordinates[0], hit_b_coordinates[0])
    overlap_end = min(hit_a_coordinates[1], hit_b_coordinates[1])
    return overlap_end - overlap_start + 1


def is_not_significantly_overlapping(stored_hits_list, hit_coordinates):
    if not stored_hits_list:
        return True

    for stored_hit_coordinates in stored_hits_list:
        overlap = get_overlap_length(stored_hit_coordinates, hit_coordinates)
        overlap_ratio = overlap / min([
                        (stored_hit_coordinates[1] -
                         stored_hit_coordinates[0] + 1),
                        (hit_coordinates[1] - hit_coordinates[0] + 1)])
        if overlap_ratio > args.max_overlap_ratio:
            return False

    return True


previous_gene = ""
non_overlapping_hits = []
# Run over stdin
for line in fileinput.input(input):
    line = line.rstrip()
    # fields: gene, gene_len, cog, cog_len, evalue, bitscore,
    #         hmm_start, hmm_stop, env_start, env_stop
    fields = line.split()
    fields[1] = int(fields[1])
    fields[3] = int(fields[3])
    fields[6] = int(fields[6])
    fields[7] = int(fields[7])
    fields[8] = int(fields[8])
    fields[9] = int(fields[9])
    alignment_length = fields[9] - fields[8] + 1

    if fields[0] != previous_gene:
        # Now we look at the hits for a new gene
        non_overlapping_hits = []
        previous_gene = fields[0]

    if (alignment_ratio_is_met(fields[1], fields[3], alignment_length) and
        is_not_significantly_overlapping(non_overlapping_hits, [fields[8],
                                                                fields[9]])):
        non_overlapping_hits.append([fields[8], fields[9]])
        fake_percent_id = ((fields[7] - fields[6] + 1) / fields[3]) * 100
        # Print out valid hits in gff format.
        print(fields[0] + "\t" +
              args.hmmer_version + "\t" +
              fields[2] + "\t" +
              str(fields[8]) + "\t" +
              str(fields[9]) + "\t" +
              fields[5] + "\t.\t.\t" +
              "ID=" + fields[0] + "_" + str(fields[8]) + "_" + str(fields[9]) +
              ";fake_percent_id=" + str("%.2f" % fake_percent_id) +
              ";alignment_length=" + str(alignment_length) +
              ";e-value=" + fields[4] +
              ";model_start=" + str(fields[6]) + ";model_end=" + str(fields[7]))
