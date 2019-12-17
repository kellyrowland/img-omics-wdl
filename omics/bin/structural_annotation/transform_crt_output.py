#!/usr/bin/env python3

import argparse
import os
import sys

parser = argparse.ArgumentParser(description="""This script takes in the CRT
                                 (1.8 or above) output and creates one gff
                                 file and one file needed by IMG that lists the
                                 repeat and spacer sequences explicitly.""")
parser.add_argument("crt_output", metavar="crt_output",
                    help="file path to the CRT output file")
parser.add_argument("crt_version", metavar="crt_version",
                    help="the CRT version used to create the output file")
args = parser.parse_args()


def get_gff_lines(seq_name, lines):
    gff_line_start = seq_name + "\t" + args.crt_version + "\t"
    start = lines[0][0]
    end = str(int(lines[-2][0]) + len(lines[-2][1]) - 1)
    gene_id = seq_name + "_" + start + "_" + end
    numbers = [s for s in lines[-1] if s.isdigit()]
    gff_print_str = gff_line_start + "CRISPR\t"
    gff_print_str += start + "\t" + end + "\t.\t.\t.\tID=" + gene_id
    gff_print_str += ";number_of_repeats=" + numbers[0]
    gff_print_str += ";average_repeat_length=" + numbers[1]
    gff_print_str += ";median_repeat_length=" + numbers[2]
    gff_print_str += ";average_spacer_length=" + numbers[3]
    gff_print_str += ";median_spacer_length=" + numbers[4] + "\n"
    gff_line_start += "repeat_unit\t"
    repeat_unit_counter = 0
    for line_split in lines[0:-1]:
        repeat_unit_counter += 1
        gff_print_str += gff_line_start
        gff_print_str += line_split[0] + "\t"
        gff_print_str += str(int(line_split[0]) + len(line_split[1]) - 1)
        gff_print_str += "\t.\t.\t.\t"
        gff_print_str += "ID=" + gene_id + "_DR" + str(repeat_unit_counter)
        gff_print_str += ";Parent=" + gene_id + "\n"
    return gff_print_str

def get_crispr_lines(seq_name, lines, element_counter):
    crispr_print_str = ""
    for fields in lines[0:-1]:
        crispr_print_str += (seq_name + "\t" + element_counter + "\t" +
                             fields[0] + "\t" + fields[1] + "\t" + fields[2] +
                             "\tc\n")
    """The last line of an element has no spacer."""
    crispr_print_str += (seq_name + "\t" + element_counter + "\t" +
                         lines[-1][0] + "\t" + lines[-1][1] + "\t\tc\n")
    return crispr_print_str

gff_out, ignore = args.crt_output.rsplit("_", 1)
crispr_out = gff_out
gff_out += "_crt.gff"
crispr_out += "_crt.crisprs"
fw_gff = open(gff_out, "w")
fw_crisprs = open(crispr_out, "w")

seq_name = ""
element_counter = 0
fr = open(args.crt_output, "r")
for line in fr:
    if line.startswith("SEQUENCE"):
        fields = line.strip().split()
        seq_name = fields[1]
        element_counter = 0
    elif line.startswith("POSITION"):
        element_counter += 1
        line = fr.readline() # Ignore next line
        lines = []
        line = fr.readline() # Line of first repeat
        while not line.startswith("---"):
            lines.append(line.split()) # Position, Repeat Seq, Spacer Seq
            line = fr.readline()
        fw_crisprs.write(get_crispr_lines(seq_name, lines,
                                         str(element_counter)))
        line = fr.readline() # Contains # of repeats, their average length, etc
        lines.append(line.split())
        fw_gff.write(get_gff_lines(seq_name, lines))
