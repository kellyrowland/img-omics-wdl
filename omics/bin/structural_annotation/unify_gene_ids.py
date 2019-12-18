#!/usr/bin/env python3

import argparse
import os

parser = argparse.ArgumentParser(description="""This script changes the gene
                                 IDs in a gff and one ore more fasta files to
                                 be of the format
                                 <sequence_id>_<gene_start>_<gene_stop>
                                 uniformly across all files.

                                 It was written to unify the gff and fasta
                                 output files created by Prodigal and
                                 GeneMark.""")
parser.add_argument("gff", metavar="gff_file",
                    help="file path to the gff file")
parser.add_argument("fasta_files", metavar="fasta_files", nargs="+",
                    help="file path(s) to the gene/protein fasta file(s)")
args = parser.parse_args()


partial_gene_info = {}
seq_to_translation_table = {}

"""Rewrite the gff file"""
tmp_file = os.path.dirname(os.path.abspath(args.gff)) + "/tmp_" + str(os.getpid()) + ".gff"
fw = open(tmp_file, "w")
fr = open(args.gff, "r")
tool_and_version = None
translation_table = None
for line in fr:
    if not line or line[0] == "#" or line == "\n":
        fw.write(line)
        if line and line.startswith("# Model Data"):          # Prodigal
            tool_and_version = line.split("version=")[1].split(";")[0]
            tool_and_version = tool_and_version.replace(".", " ", 1)
            translation_table = line.split("_table=")[1].split(";")[0]
        elif line and line.startswith("# translation table"): # GeneMark
            translation_table = line.rstrip().split(": ")[1]
        elif line and line.startswith("# GeneMark.hmm-2"):    # GeneMark
            fields = line.rstrip().split(" ")
            tool_and_version = fields[1] + " v" + fields[-1]
        continue
    fields = line.rstrip().split("\t")
    if fields[0] not in seq_to_translation_table:
        seq_to_translation_table[fields[0]] = translation_table
    fields[1] = tool_and_version
    gene_id = fields[0] + "_" + fields[3] + "_" + fields[4]
    attributes = []
    attributes.append("ID=" + gene_id)
    attributes.append("translation_table=" + translation_table)
    partial_index = fields[8].find("partial")
    if partial_index >= 0:
        partial_code_index = partial_index + 8
        partial_code = fields[8][partial_code_index:partial_code_index+2]
        if partial_code == "10":
            if fields[6] == "+":
                attributes.append("partial=5'")
                partial_gene_info[gene_id] = "partial=5'"
            else:
                attributes.append("partial=3'")
                partial_gene_info[gene_id] = "partial=3'"
        elif partial_code == "01":
            if fields[6] == "+":
                attributes.append("partial=3'")
                partial_gene_info[gene_id] = "partial=3'"
            else:
                attributes.append("partial=5'")
                partial_gene_info[gene_id] = "partial=5'"
        elif partial_code == "11":
            attributes.append("partial=5',3'")
            partial_gene_info[gene_id] = "partial=5',3'"
    fields[8] = ";".join(attributes)
    fw.write("\t".join(fields) + "\n")
fr.close()
fw.close()
os.rename(tmp_file, args.gff)


def rewrite_prodigal_line(line):
    fields = line.split(" # ")
    # Rename the sequence
    seq_id = fields[0][1:].rsplit("_", 1)[0]
    gene_id = seq_id + "_" + fields[1] + "_" + fields[2]
    fields[0] = ">" + gene_id
    # Add translation table
    fields[4] = "tt=" + seq_to_translation_table[seq_id]
    # Add partial info if necessary
    if gene_id in partial_gene_info:
        fields.append(partial_gene_info[gene_id])
    if fields[3] == "1":
        fields[3] = "+"
    else:
        fields[3] = "-"
    fields.insert(1, tool_and_version)
    """The seq name line should look like this:
        >gene_id # tool_and_version # start # stop # strand # tt # partial"""
    return (" # ".join(fields) + "\n")


def rewrite_genemark_line(line):
    new_fields = []
    fields = line.rstrip().split(" ")
    gene_id = fields[1] + "_" + fields[2] + "_" + fields[3]
    new_fields.append(">" + gene_id)
    new_fields.append(tool_and_version)
    new_fields.append(fields[2])
    new_fields.append(fields[3])
    new_fields.append(fields[4])
    new_fields.append("tt=" + seq_to_translation_table[fields[1]])
    if gene_id in partial_gene_info:
        new_fields.append(partial_gene_info[gene_id])
    new_fields.append(fields[5])
    """The seq name line should look like this:
        >gene_id # tool_and_version # start # stop # strand # tt # partial #
        gene_type"""
    return (" # ".join(new_fields) + "\n")


"""Now loop over the fasta files."""
for fasta_file in args.fasta_files:
    """Creating a tmp file for every fasta file in the same folder to avoid the
    OSError: [Errno 18] Invalid cross-device link
    problem that can happen with os.rename."""
    tmp_file = os.path.dirname(os.path.abspath(fasta_file)) + "/tmp_" + str(os.getpid()) + ".fasta"
    fw = open(tmp_file, "w")
    fr = open(fasta_file, "r")
    for line in fr:
        if line[0] == ">":
            if "Prodigal" in tool_and_version:
                line = rewrite_prodigal_line(line)
            else:
                line = rewrite_genemark_line(line)
        fw.write(line)
    fr.close()
    fw.close()
    os.rename(tmp_file, fasta_file)

