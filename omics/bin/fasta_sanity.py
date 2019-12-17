#!/usr/bin/env python3

import argparse
import fileinput
import subprocess


parser = argparse.ArgumentParser(description="""This script will execute basic
                                 fasta sanity checks.
                                 It'll ignore all contigs shorter than a given
                                 threshold (default is 150 bp) or sequences that
                                 only consist of Ns.
                                 All non-ACGT characters will get replaced with
                                 Ns.
                                 If a sequence name prefix is provided the
                                 sequences will get renamed in the format of
                                 <sequence_name_prefix>_<0-padded_sequence_number>.
                                 All alterations, be it ignored or renamed
                                 sequences or characters replaced with N will
                                 be reported in respective tsv files.""",
formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("input_fasta",
                    help="""path to the input fasta file containing the
                    nucleotide sequences""")
parser.add_argument("output_fasta",
                    help="""output file name to write the qced nucleotide
                    sequences to (in fasta format)""")
parser.add_argument("-l", "--minimum_seq_length",
                    type=int, default=150,
                    help="""every sequence shorter than this cutoff will get
                    ignored""")
parser.add_argument("-p", "--seq_name_prefix",
                    help="""if provided, sequence names will get replaced with
                    <seq_name_prefix>_<0-padded_seq_number>""")
parser.add_argument("-v", "--version", action="version",
                    version="%(prog)s 1.0.0")
args = parser.parse_args()


if args.seq_name_prefix:
    # Determine padding width
    bash_cmd = "grep -c '^>' " + args.input_fasta + "; exit 0"
    seq_count = subprocess.check_output(bash_cmd, shell=True).decode("utf-8")
    seq_count = seq_count.rstrip()
    padding_width = len(seq_count)
    if padding_width == 1:
        padding_width = 2
seq_counter = 0
allowed = frozenset(["A", "C", "G", "T", "N"])
output_fasta_fh = None
kicked_out_seqs_fh = None
names_map_fh = None
replaced_n_count_fh = None


def write_sequence(old_seq_name, sequence, replaced_n_count):
    seq_name = old_seq_name
    if args.seq_name_prefix:
        global seq_counter
        seq_counter += 1
        seq_name = (args.seq_name_prefix + "_" +
                    str(seq_counter).zfill(padding_width))
        global names_map_fh
        if not names_map_fh:
            names_map_tsv = args.output_fasta[0:args.output_fasta.rfind("_")]
            names_map_tsv += "_contig_names_mapping.tsv"
            names_map_fh = open(names_map_tsv, 'w')
        names_map_fh.write(old_seq_name + "\t" + seq_name + "\n")
    if replaced_n_count > 0:
        global replaced_n_count_fh
        if not replaced_n_count_fh:
            replaced_n_count_tsv = (
                            args.output_fasta[0:args.output_fasta.rfind("_")])
            replaced_n_count_tsv += "_replaced_n_count.tsv"
            replaced_n_count_fh = open(replaced_n_count_tsv, 'w')
        replaced_n_count = str(replaced_n_count)
        replaced_n_count_fh.write(seq_name + "\t" + replaced_n_count + "\n")
    global output_fasta_fh
    if not output_fasta_fh:
        output_fasta_fh = open(args.output_fasta, 'w')
    output_fasta_fh.write(">" + seq_name + "\n")
    for i in range(0, len(sequence), 60):
        output_fasta_fh.write(sequence[i:i+60]  + "\n")


def write_ignore_reason(seq_name, reason):
    global kicked_out_seqs_fh
    if not kicked_out_seqs_fh:
        kicked_out_seqs_tsv = args.output_fasta[0:args.output_fasta.rfind("_")]
        kicked_out_seqs_tsv += "_ignored_contigs.tsv"
        kicked_out_seqs_fh = open(kicked_out_seqs_tsv, 'w')
    kicked_out_seqs_fh.write(seq_name + "\t" + reason + "\n")


def check_and_clean_and_output(seq_name, sequence):
    # Check length
    if len(sequence) < args.minimum_seq_length:
        write_ignore_reason(seq_name, "sequence below length cutoff of " +
                            str(args.minimum_seq_length) + " nucleotides")
        return
    # Capitalize and replace all non-ACGT letters with N
    sequence = sequence.upper()
    pre_n_count = sequence.count("N")
    sequence = "".join([c if c in allowed else "N" for c in sequence])
    # Check if we only have Ns left
    if sequence.count("N") == len(sequence):
        write_ignore_reason(seq_name, "sequence only contains Ns")
        return
    post_n_count = sequence.count("N")
    write_sequence(seq_name, sequence, (post_n_count - pre_n_count))


seq_name = ""
sequence = ""
fasta_fh = open(args.input_fasta, "r")
for line in fasta_fh:
    line = line.rstrip()
    if line.startswith(">"):
        if seq_name:
            check_and_clean_and_output(seq_name, sequence)
        name_line_split = line[1:].split()
        seq_name = name_line_split[0]
        sequence = ""
    else:
        sequence += line
check_and_clean_and_output(seq_name, sequence)

