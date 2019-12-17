#!/usr/bin/env python3

import argparse
import fileinput
from os.path import dirname


parser = argparse.ArgumentParser(description="""This script calculates the
                                 coding density and checks the length of the
                                 intergenic regions.
                                 If there are any intergenic regions above the
                                 given threshold or the coding density is
                                 below or above the given min and max coding
                                 density thresholds, the script will create a
                                 QC_NEEDED.txt file in the folder of the input
                                 fasta file and list the reason(s).
                                 Too long intergenic regions will also get
                                 reported (if they don't contain too many Ns)
                                 in a tsv file listing the sequence name, the
                                 start and stop position of the intergenig
                                 region as well as the length of the gap.""",
formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("contigs_fasta",
                    help="""path to the input fasta file containing the
                    nucleotide sequences""")
parser.add_argument("structural_annotation_gff",
                    type=argparse.FileType('r'),
                    help="""path to the gff file containing the structural
                    annotations for the contigs in contigs_fasta""")
parser.add_argument("-g", "--gap_len_cutoff",
                    type=int, default=3000,
                    help="""determines maximal allowed length of intergenic
                    regions""")
parser.add_argument("-n", "--gap_n_stretch_cutoff",
                    type=int, default=100,
                    help="""if an intergenic region contains a stretch of this
                    many Ns or more it will not get reported""")
parser.add_argument("-l", "--lower_coding_density_cutoff",
                    type=float, default=0.7,
                    help="""determines the minimum allowed coding density
                    percentage""")
parser.add_argument("-u", "--upper_coding_density_cutoff",
                    type=float, default=1.0,
                    help="""determines the maximum allowed coding density
                    percentage""")
parser.add_argument("-v", "--version", action="version",
                    version="%(prog)s 1.0.0")
args = parser.parse_args()


"""This script will only be executed for genomes. Since their sizes will only
be a few MB, we'll keep their sequences in memory for easy lookup of the
intergenic regions (to check how many Ns they contain)."""

def create_seq_lookup():
    lookup = {}
    contigs_fasta_fh = open(args.contigs_fasta, 'r')
    for line in contigs_fasta_fh:
        line = line.rstrip()
        if not line:
            continue
        if line.startswith(">"):
            seq_name = line[1:]
            lookup[seq_name] = ""
        else:
            lookup[seq_name] += line
    contigs_fasta_fh.close()
    return lookup



def initialize_coding_bp_tracker_and_return_total_bp(coding_bp_tracker,
                                                     seq_lookup):
    total_bp = 0
    for seq_name, seq in seq_lookup.items():
        seq_len = len(seq)
        total_bp += seq_len
        coding_bp_tracker[seq_name] = [0] * seq_len
    return total_bp



def gap_has_not_too_many_ns(seq_lookup, seq_name,
                            gap_start, gap_stop, n_stretch):
    intergenic_region = seq_lookup[seq_name][gap_start-1:gap_stop]
    if intergenic_region.find(n_stretch) >= 0:
        return False
    return True



def write_to_qc_needed_file(reason):
    global qc_needed_fh
    if not qc_needed_fh:
        dir_name = dirname(args.contigs_fasta)
        if dir_name == "":
            dir_name = "."
        qc_needed_file = dir_name + "/QC_NEEDED.txt"
        qc_needed_fh = open(qc_needed_file, 'w')
    qc_needed_fh.write(reason + "\n")


gap_tsv_fh = None
qc_needed_fh = None
seq_lookup = create_seq_lookup()
coding_bp_tracker = {}
total_bp = initialize_coding_bp_tracker_and_return_total_bp(coding_bp_tracker,
                                                            seq_lookup)
n_stretch = "N" * args.gap_n_stretch_cutoff
total_coding_bp = 0
previous_gene_stop = 0
big_gaps_count = 0
seq_name = ""
for line in args.structural_annotation_gff:
    fields = line.split("\t")
    if "Parent=" in fields[8]:
        continue
    if fields[0] != seq_name:
        previous_gene_stop = 0
        seq_name = fields[0]
    gene_start = int(fields[3])
    gene_stop = int(fields[4])
    for index in range(gene_start-1, gene_stop):
        coding_bp_tracker[seq_name][index] = 1
    if gene_stop < previous_gene_stop:
        continue
    gap_length = gene_start - previous_gene_stop - 1
    if (gap_length > args.gap_len_cutoff and
            gap_has_not_too_many_ns(seq_lookup, seq_name,
                                    previous_gene_stop,
                                    gene_start, n_stretch)):
        if not gap_tsv_fh:
            gap_tsv = args.contigs_fasta[0:args.contigs_fasta.rfind("_")]
            gap_tsv += "_too_long_intergenig_regions.tsv"
            gap_tsv_fh = open(gap_tsv, 'w')
        gap_tsv_fh.write(fields[0] + "\t" + str(previous_gene_stop+1) + "\t" +
                         str(gene_start-1) + "\t" + str(gap_length) + "\n")
        big_gaps_count += 1
    previous_gene_stop = gene_stop


total_coding_bp = 0
for bp_tracker in coding_bp_tracker.values():
    total_coding_bp += sum(bp_tracker)
coding_density = total_coding_bp / total_bp
if coding_density < args.lower_coding_density_cutoff:
    reason = "The coding density is only "
    reason += str("%.2f" % (coding_density * 100)) + "% and with that below "
    reason += "the minimum cutoff of "
    reason += str("%.2f" % (args.lower_coding_density_cutoff * 100)) + "%."
    write_to_qc_needed_file(reason)
elif coding_density > args.upper_coding_density_cutoff:
    reason = "The coding density is "
    reason += str("%.2f" % (coding_density * 100)) + "% and with that above "
    reason += "the maximum cutoff of "
    reason += str("%.2f" % (args.upper_coding_density_cutoff * 100)) + "%."
    write_to_qc_needed_file(reason)
if big_gaps_count > 0:
    reason = "There "
    if big_gaps_count == 1:
        reason += "is 1 too long intergenig region.\n"
    else:
        reason += "are "+str(big_gaps_count)+" too long intergenic regions "
    reason += "(see " + gap_tsv + " for details)."
    write_to_qc_needed_file(reason)
# Close files if necessary.
if gap_tsv_fh:
    gap_tsv_fh.close()
if qc_needed_fh:
    qc_needed_fh.close()

