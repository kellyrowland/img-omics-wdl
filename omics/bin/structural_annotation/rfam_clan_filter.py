#!/usr/bin/env python3

import sys
import argparse
import fileinput


parser = argparse.ArgumentParser(description="""This script expects to get via
stdin columns 1, 3, 4, 6, 7, 8, 9, 10, 11, 15, 16 (gene_id, accession, model,
model_start, model_stop, gene_start, gene_stop, strand, trunc, bitscore,
e-value) from a tblout output created by cmsearch that got sorted via
'sort -k1,1 -k10,10nr -k11,11n' and only contains hits that achieved the
inclusion threshold (column 17 (inc) was set to '!').
As arguments it needs the cmsearch version used (for the GFF output), the Rfam
claninfo file that can get downloaded together with the database and the
Rfam-IMG-NCBI feature lookup file.
For all hits to a seqeunce the script then checks if two hits overlap and if so
additionally checks if the two hits belong to the same class. If that's the
case the lower scoring hit gets removed.""")
parser.add_argument("cmsearch_version", help="""cmsearch version used to create
                    the input""")
parser.add_argument("clan_info_file", type=argparse.FileType('r'),
                    help="""the full path to the clan info file""")
parser.add_argument("feature_lookup_file", type=argparse.FileType('r'),
                    help="""the full path to the feature lookup file""")
parser.add_argument("-v", "--version", action="version",
                    version="%(prog)s 1.0.0")
args, input = parser.parse_known_args()



class Clan_Lookup:
    def __init__(self, clan_info_file):
        self.create_clan_lookup(clan_info_file)


    def create_clan_lookup(self, clan_info_file):
        self.lookup = {}
        for line in clan_info_file:
            clan, rest = line.rstrip().split("\t", 1)
            product_names = rest.split()
            for product_name in product_names:
                self.lookup[product_name] = clan


    def belong_to_the_same_clan(self, hit_a, hit_b):
#        print("Clan check for " + hit_a.id + " and" + hit_b.id)
#        print("Acessions are " + hit_a.accession + " and" + hit_b.accession)
        if (hit_a.accession not in self.lookup or
                hit_b.accession not in self.lookup):
            return False
#        print("Clans are " + self.lookup[hit_a.accession] + " and " +
#              self.lookup[hit_b.accession])
        if self.lookup[hit_a.accession] == self.lookup[hit_b.accession]:
            return True
        return False




class Feature_Lookup:
    def __init__(self, feature_lookup_file):
        self.create_feature_lookup(feature_lookup_file)


    def create_feature_lookup(self, feature_lookup_file):
        self.lookup = {}
        # Ignore header line
        feature_lookup_file.readline()
        for line in feature_lookup_file:
            fields = line.rstrip().split("\t")
            model = fields[0]
            self.lookup[model] = {} # Create lookup entry for each model
            self.lookup[model]["description"] = fields[1]
            self.lookup[model]["entry_type"] = fields[2]
            self.lookup[model]["ncbi_locus_type"] = fields[3]
            self.lookup[model]["ncRNA_class"] = fields[4]
            self.lookup[model]["obsolete_flag"] = fields[5]
            if len(fields) == 7:
                self.lookup[model]["additional_qualifiers"] = fields[6]
            else:
                self.lookup[model]["additional_qualifiers"] = ""


    def get_ncbi_locus_type(self, model):
        if model not in self.lookup:
            print("Model '" + model + "' not in feature lookup table.",
                  file=sys.stderr)
            sys.exit(1)
        return self.lookup[model]["ncbi_locus_type"]


    def get_additional_qualifiers(self, model):
        if model not in self.lookup:
            print("Model '" + model + "' not in feature lookup table.",
                  file=sys.stderr)
            sys.exit(1)
        return self.lookup[model]["additional_qualifiers"].replace("\"", "")




class Hit:
    def __init__(self, fields):
        self.seq_name = fields[0]
        self.accession = fields[1]
        self.model = fields[2]
        self.model_start = fields[3]
        self.model_end = fields[4]
        self.seq_start = int(fields[5])
        self.seq_end = int(fields[6])
        self.strand = fields[7]
        if self.strand == "-":
            self.seq_start, self.seq_end = self.seq_end, self.seq_start
        self.trunc = fields[8]
        self.bitscore = fields[9]
        self.evalue = fields[10]
#        self.id = (self.seq_name + "_" + str(self.seq_start) + "_" +
#                   str(self.seq_end))


    def to_gff(self, feature_lookup):
        gff_line = (self.seq_name + "\t" +
                    args.cmsearch_version + "\t" +
                    feature_lookup.get_ncbi_locus_type(self.model) + "\t" +
                    str(self.seq_start) + "\t" +
                    str(self.seq_end) + "\t" +
                    self.bitscore + "\t" +
                    self.strand + "\t" +
                    ".\t" + # Phase
                    "ID=" + self.seq_name + "_" +
                    str(self.seq_start) + "_" + str(self.seq_end) +
                    ";e-value=" + self.evalue +
                    ";model=" + self.model +
                    ";accession=" + self.accession +
                    ";model_start=" + self.model_start +
                    ";model_end=" + self.model_end +
                    ";" + feature_lookup.get_additional_qualifiers(self.model))
        if self.trunc != "no":
            gff_line += ";partial=" + self.trunc.replace("&", ",")
        print(gff_line)


    def get_overlap_length(self, other_hit):
#        print("Overlap between " + self.id + " and " + other_hit.id)
        # First check if the overlaps overlap at all
        if self.seq_start > other_hit.seq_end or self.seq_end < other_hit.seq_start:
#            print("\t is 0!")
            return 0

        overlap_start = max(self.seq_start, other_hit.seq_start)
        overlap_end = min(self.seq_end, other_hit.seq_end)
#        print("\t is " + str(overlap_end - overlap_start + 1))
        return overlap_end - overlap_start + 1




class Sequence_Hits:
    def __init__(self):
        self.hits = []

    def add_hit(self, hit):
        self.hits.append(hit)

    def remove_clan_overlaps(self, clan_lookup):
        if len(self.hits) > 0:
            ignore = set()
            for i in range(0, len(self.hits)-1):
                if i not in ignore:
                    for j in range(i+1, len(self.hits)):
                        if j not in ignore:
                            if (self.hits[i].get_overlap_length(self.hits[j]) > 0
                                 and
                                 clan_lookup.belong_to_the_same_clan(self.hits[i],
                                                                     self.hits[j])):
                                ignore.add(j)
            cleaned_hits = [self.hits[0]]
            for i in range(1, len(self.hits)):
                if i not in ignore:
                    cleaned_hits.append(self.hits[i])
            self.hits = cleaned_hits

    def to_gff(self, feature_lookup):
        # Sort by start coordinate first
        self.hits.sort(key=lambda h: h.seq_start)
        for hit in self.hits:
            hit.to_gff(feature_lookup)



# Start script by building lookups
clan_lookup = Clan_Lookup(args.clan_info_file)
feature_lookup = Feature_Lookup(args.feature_lookup_file)

current_seq = ""
sequence_hits = Sequence_Hits()
# Run over stdin
for line in fileinput.input(input):
    hit = Hit(line.rstrip().split())

    if hit.seq_name != current_seq:
        # Now we look at the hits for a new sequence
        sequence_hits.remove_clan_overlaps(clan_lookup)
        sequence_hits.to_gff(feature_lookup)
        sequence_hits = Sequence_Hits()
        current_seq = hit.seq_name

    sequence_hits.add_hit(hit)

sequence_hits.remove_clan_overlaps(clan_lookup)
sequence_hits.to_gff(feature_lookup)

