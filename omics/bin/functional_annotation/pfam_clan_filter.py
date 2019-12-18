#!/usr/bin/env python3

import argparse
import fileinput


parser = argparse.ArgumentParser(description="""This script expects to get
columns 1, 4, 5, 6, 13, 14, 16, 17, 20, 21 (gene_id, model, accession,
hmm_length, evalue, bitscore, hmm_start, hmm_stop, env_start, env_stop) from a
domtblout output created by hmmsearch sorted via 'sort -k1,1 -k6,6nr -k5,5n'.
As arguments it needs the hmmsearch version used (for the GFF output) and the
Pfam claninfo file that can get downloaded together with the database.
For all hits of a gene the script then checks if two hits overlap and if so
additionally checks if the two hits belong to the same class. If that's the
case the lower scoring hit gets removed.""")
parser.add_argument("hmmsearch_version", help="""hmmsearch version used to
                    create the input""")
parser.add_argument("clan_info_file", type=argparse.FileType('r'),
                    help="""the full path to the clan info file""")
parser.add_argument("-v", "--version", action="version",
                    version="%(prog)s 1.0.0")
args, input = parser.parse_known_args()


clan_lookup = {}


def create_clan_lookup():
    for line in args.clan_info_file:
        pfam, clan, rest = line.rstrip().split("\t", 2)
        if not clan:
            continue
        clan_lookup[pfam] = clan


def get_overlap_length(hit_a, hit_b):
    # First check if the overlaps overlap at all
    if hit_a[8] > hit_b[9] or hit_a[9] < hit_b[8]:
        return 0

    overlap_start = max(hit_a[8], hit_b[8])
    overlap_end = min(hit_a[9], hit_b[9])
    return overlap_end - overlap_start + 1


def remove_clan_overlaps(gene_hits):
    if not gene_hits:
        return gene_hits

    ignore = set()
    for i in range(0, len(gene_hits)-1):
        if i not in ignore:
            for j in range(i+1, len(gene_hits)):
                if j not in ignore:
                    if (get_overlap_length(gene_hits[i], gene_hits[j]) > 0
                            and gene_hits[i][2] in clan_lookup
                            and gene_hits[j][2] in clan_lookup
                            and clan_lookup[gene_hits[i][2]] == clan_lookup[gene_hits[j][2]]):
                        ignore.add(j)
    cleaned_gene_hits = [gene_hits[0]]
    for i in range(1, len(gene_hits)):
        if i not in ignore:
            cleaned_gene_hits.append(gene_hits[i])

    return cleaned_gene_hits


def output_in_gff(gene_hits):
    for hit in gene_hits:
        fake_percent_id = ((fields[7] - fields[6] + 1) / fields[3]) * 100
        alignment_length = hit[9] - hit[8] + 1
        print(hit[0] + "\t" +
              args.hmmsearch_version + "\t" +
              hit[2] + "\t" +
              str(hit[8]) + "\t" +
              str(hit[9]) + "\t" +
              hit[5] + "\t" +
              ".\t.\t" +
              "ID=" + hit[0] + "_" + str(hit[8]) + "_" + str(hit[9]) +
              ";Name=" + hit[1] +
              ";fake_percent_id=" + str("%.2f" % fake_percent_id) +
              ";alignment_length=" + str(alignment_length) +
              ";e-value=" + hit[4] +
              ";model_start=" + str(hit[6]) +
              ";model_end=" + str(hit[7]))



create_clan_lookup()
current_gene = ""
gene_hits = []
# Run over stdin
for line in fileinput.input(input):
    # fields: gene, model, accession, model_len, evalue,
    #         bitscore, hmm_start, hmm_stop, env_start, env_stop
    fields = line.rstrip().split()
    fields[2] = fields[2][:fields[2].find(".")]
    fields[3] = int(fields[3])
    fields[6] = int(fields[6])
    fields[7] = int(fields[7])
    fields[8] = int(fields[8])
    fields[9] = int(fields[9])

    if fields[0] != current_gene:
        # Now we look at the hits for a new gene
        gene_hits = remove_clan_overlaps(gene_hits)
        output_in_gff(gene_hits)
        gene_hits = []
        current_gene = fields[0]

    gene_hits.append(fields)

gene_hits = remove_clan_overlaps(gene_hits)
output_in_gff(gene_hits)
