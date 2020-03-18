#!/usr/bin/env python3

import argparse
#import sys


parser = argparse.ArgumentParser(description="""This script parses one or more
                                 tRNAscan-SE output files. It expects that the
                                 output files got created with '&>' instead of
                                 the tRNAscan-SE's -o option, so that it also
                                 contains the tRNAscan version and the used
                                 search mode.
                                 If mutliple output files are given it picks
                                 the one that contains the highest number of
                                 tRNAs with know isotype. In case two models
                                 have the same high count the one created with
                                 the higher ranking search mode gets pciked.
                                 The ranking is as follows:
                                 Bacterial > Archaeal > Eukaryotic""")
parser.add_argument("trnascan_output_files", type=argparse.FileType('r'),
                    nargs="+", metavar="trnascan_output_file",
                    help="""full path(s) to tRNAscan-SE output file(s)""")
parser.add_argument("-v", "--version", action="version",
                    version="%(prog)s 1.0.0")
args = parser.parse_args()


SEARCH_MODE_HIERACHY = {"Bacterial": 1, "Archaeal": 2, "Eukaryotic": 3}

COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


class tRNA:
    def __init__(self, fields):
        self.start = fields[0]
        self.end = fields[1]
        if int(self.start) > int(self.end):
            self.start, self.end = self.end, self.start
            self.strand = "-"
        else:
            self.strand = "+"
        self.type = fields[2]
        self.anti_codon = fields[3]
        self.has_known_isotype = True
        self.product = None
        self.codon = None
        if self.type == "Undet" or self.anti_codon == "NNN":
            self.product = "tRNA_Xxx"
            self.codon = "Pseudo"
            self.has_known_isotype = False
        self.intron_start = fields[4]
        self.intron_end = fields[5]
        self.score = fields[6]
        if len(fields) > 7:
            self.note = fields[7]
        else:
            self.note = None


    def to_gff(self, contig_name, source, search_mode):
        """These operations could be executed during object initialization, but
        they are actually only needed for the ones that get printed out."""
        intron_str = None
        if self.intron_start != "0" or self.intron_end != "0":
            if self.strand == "-":
                self.intron_start, self.intron_end = (self.intron_end,
                                                      self.intron_start)
                self.intron_start = ",".join(coord for coord in
                                             self.intron_start.split(",")[::-1])
                self.intron_end = ",".join(coord for coord in
                                             self.intron_end.split(",")[::-1])
            intron_str = ";intron_start=" + self.intron_start
            intron_str += ";intron_end=" + self.intron_end

        if not self.product:
            self.product = "tRNA_" + self.type + "_" + self.anti_codon

        if not self.codon:
            self.codon = "".join([COMPLEMENT[base]
                                  for base in self.anti_codon[::-1]])

        gff_line = (contig_name + "\t" +
                    source + "\t" +
                    "tRNA\t" +
                    self.start + "\t" +
                    self.end + "\t" +
                    self.score + "\t" +
                    self.strand + "\t" +
                    ".\t" +                   # Phase
                    "ID=" + contig_name + "_" + self.start + "_" + self.end +
                    ";product=" + self.product +
                    ";codon=" + self.codon + ";used_search_mode=" + search_mode)
        if intron_str:
            gff_line += intron_str
        if self.note:
            gff_line += ";note=" + self.note

        print(gff_line)



class Contig_Results:
    def __init__(self, contig_name, source, search_mode):
        self.tRNAs = []
        self.contig_name = contig_name
        self.source = source
        self.search_mode = search_mode


    def add_trna(self, trna):
        self.tRNAs.append(trna)


    def get_known_isotopes_count(self):
        if not hasattr(self, "known_isotopes_count"):
            self.known_isotypes_count = sum([t.has_known_isotype for t in
                                             self.tRNAs])
        return self.known_isotypes_count


    def is_better(self, other_contig_results):
        if (self.get_known_isotopes_count() >
                other_contig_results.get_known_isotopes_count()):
            return True
        elif ((self.get_known_isotopes_count() ==
               other_contig_results.get_known_isotopes_count()) and
              (SEARCH_MODE_HIERACHY[self.search_mode] <
               SEARCH_MODE_HIERACHY[other_contig_results.search_mode])):
            return True
        return False


    def to_gff(self):
        # Sort tRNA list by start coordinates first
        self.tRNAs.sort(key=lambda t: int(t.start))
        for trna in self.tRNAs:
            trna.to_gff(self.contig_name, self.source, self.search_mode)



best_contigs_results = {}
for results_file in args.trnascan_output_files:
    source = None
    search_mode = None
    at_trna_lines = False
    previous_contig_name = ""
    contig_results = None
    trnas_to_ignore = set()
    for line in results_file:
        if not line or line == "\n" or line.startswith("-"):
            at_trna_lines = False
            """Above line is to accommodate for concatenated result files
            created from a split up input fasta."""
            continue
        if not source and line.startswith("tRNAscan-SE"):
            source = line.split(" - ")[0]
            continue
        if not search_mode and line.startswith("Search Mode"):
            search_mode = line.rstrip().split()[-1]
            continue
        if not at_trna_lines and line.startswith("Name"):
            line = results_file.readline()
            at_trna_lines = True
            continue
        if at_trna_lines:
            if line.startswith("Warning"):
#                print(line, file=sys.stderr)
                """This would appear before the results of a new contig and indicate which tRNA
                on that contig has a problem. We are going to ignore that tRNA."""
                if ".trna0" in line:
                    trnas_to_ignore.add(str(int(line.split(".trna0")[1].split()[0])))
                else:
                    trnas_to_ignore.add(str(int(line.split(".tRNA")[1].split("-")[0])))
                continue
            fields = line.rstrip().split()
            if trnas_to_ignore and fields[1] in trnas_to_ignore:
                # This is (one of) the tRNA(s) to ignore
                trnas_to_ignore.remove(fields[1])
                continue
            contig_name = fields[0]
            if contig_name != previous_contig_name:
                if contig_results:
                    if (previous_contig_name not in best_contigs_results or
                          contig_results.is_better(best_contigs_results[previous_contig_name])):
                        best_contigs_results[previous_contig_name] = contig_results
                contig_results = Contig_Results(contig_name, source,
                                                search_mode)
                previous_contig_name = contig_name
#            print(f'Building tRNA from line: {line}', file=sys.stderr)
            trna = tRNA(fields[2:])
            contig_results.add_trna(trna)
    if (contig_results and (contig_results.contig_name not in best_contigs_results or
            contig_results.is_better(best_contigs_results[contig_results.contig_name]))):
        best_contigs_results[contig_name] = contig_results

# Print out final results
for contig_name in sorted(best_contigs_results.keys()):
    best_contigs_results[contig_name].to_gff()

