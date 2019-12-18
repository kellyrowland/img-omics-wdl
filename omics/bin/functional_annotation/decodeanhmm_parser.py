#!/usr/bin/env python3

import argparse
import fileinput
import sys

parser = argparse.ArgumentParser(description="""This script parses the
                                 decodeanhmm output and transforms it into gff
                                 format, which gets printed to stdout.
                                 It expects that decodeanhmm was run with -N 1
                                 and with the -PrintNumbers flag.
                                 The decodeanhmm version needs to be provided
                                 as an argument and the decodeanhmm output can
                                 get read from stdin or a file.""")
parser.add_argument("decodeanhmm_version", metavar="decodeanhmm_version",
                   help="decodeanhmm version")
args, input = parser.parse_known_args()

gene_id = ""
for line in fileinput.input(input):
    if line.startswith(">"):
        # We got a new sequence
        gene_id = line[1:].rstrip().split(" ")[0]
        continue
    if line.startswith("%pred"):
        info, prediciton = line.rstrip().split(": ")
        types_and_positions = prediciton.split(", ")
        if len(types_and_positions) > 1:
            # This gene contains a TMH
            for type_and_position in types_and_positions:
                feature_type, start, end = type_and_position.split(" ")
                if feature_type.lower() == "i":
                    feature_type = "Inside"
                elif feature_type.lower() == "m":
                    feature_type = "TMhelix"
                elif feature_type.lower() == "o":
                    feature_type = "Outside"
                else:
                    print("Unknown feature type: " + feature_type, sys.stderr)
                    sys.exit(1)
                gff_line = gene_id + "\t" + args.decodeanhmm_version + "\t"
                gff_line += feature_type + "\t" + start + "\t" + end
                gff_line += "\t.\t.\t.\tID=" + gene_id + "_" + start + "_" + end
                print(gff_line)
