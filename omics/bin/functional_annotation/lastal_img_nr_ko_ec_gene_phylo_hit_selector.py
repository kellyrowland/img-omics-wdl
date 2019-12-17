#!/usr/bin/env python3

import os
import sys
import argparse
import fileinput
from datetime import datetime

parser = argparse.ArgumentParser(description="""This script expects to get to
                                 get the blasttab+ output from a lastal run the
                                 IMG NR, as well as MD5 hash mapping file (for
                                 the sequence names in the IMG NR db) and the
                                 phylogeny mapping file for each taxon OID
                                 represented withing the the IMG NR db.
                                 It outputs three tab-delimited files listing
                                 the information IMG needs. One for the gene
                                 with KOs, another one for gene with EC number
                                 and a third one containing the phylogeny for
                                 each gene. Additionally it also outputs one
                                 gff file containing KO and EC information
                                 combined.
                                 The blasttab+ output can be provided either as
                                 a file or via STDIN.
                                 For each gene we transfer the phylogenetic
                                 assignment of the best hit to it. For a gene
                                 to get a KO assignment the top 5 hits need to
                                 contain at least two hits to genes with KO
                                 assignment. All hits to genes with KO
                                 assignment need to list the same combination
                                 of KO terms. If that's the case it reports the
                                 first hit to a gene with KO assignment that
                                 also meets the alignment length ratio, which
                                 means at least 70% of the query and subject's
                                 gene length are covered (calculated based on
                                 reported coordinates for each) for isolate
                                 projects and 70% of the shorter gene's length
                                 for metagenomes.""",
formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("project_type", choices=["isolate", "metagenome"],
                    help="either 'isolate' or 'metagenome'")
parser.add_argument("md5_mapping_file", type=argparse.FileType('r'),
                    help="filepath to the MD5 mapping file")
parser.add_argument("phylo_mapping_file", type=argparse.FileType('r'),
                    help="filepath to the phylogeny mapping")
parser.add_argument("ko_output_file", type=argparse.FileType('w'),
                    help="filepath to to the KO output file")
parser.add_argument("ec_output_file", type=argparse.FileType('w'),
                    help="filepath to the EC output file")
parser.add_argument("phylo_output_file", type=argparse.FileType('w'),
                    help="filepath to the phylo output file")
parser.add_argument("ko_ec_gff_output_file", type=argparse.FileType('w'),
                    help="filepath to the KO/EC GFF output file")
parser.add_argument("last_blasttabplus_file", nargs='?',
                    type=argparse.FileType('r'), default=sys.stdin,
                    help="""file path to the blasttab+ output (can also be
                    provided via STDIN)""")
parser.add_argument("-l", "--aln_len_ratio", metavar="aln_len_ratio",
                    type=float, default=0.7,
                    help="""the minimum required alignment length ratio""")
parser.add_argument("-m", "--min_ko_hits",
                    metavar="min_ko_hits", type=int, default=2,
                    help="""minimum number of hits with KO assignment required
                    to be in top hits""")
parser.add_argument("-n", "--number_of_top_hits", metavar="number_of_top_hits",
                    type=int, default=5,
                    help="number of top hits to look at")
parser.add_argument("-v", "--version", action="version",
                    version="%(prog)s 1.0.0")
args, input = parser.parse_known_args()


class Hit:
    def __init__(self, hit_fields):
        self.query_gene = hit_fields[0]
        self.subject_gene = hit_fields[1]
        self.percent_id = hit_fields[2]
        self.query_gene_start = int(hit_fields[6])
        self.query_gene_end = int(hit_fields[7])
        self.subject_gene_start= int(hit_fields[8])
        self.subject_gene_end = int(hit_fields[9])
        self.evalue = hit_fields[10]
        self.bitscore = hit_fields[11]
        self.query_gene_len = int(hit_fields[12])
        self.subject_gene_len = int(hit_fields[13])


    def has_ko(self, md5_lookup):
        if not hasattr(self, "ko_ids"):
            if "ko_ids" in md5_lookup[self.subject_gene]:
                self.ko_ids = md5_lookup[self.subject_gene]["ko_ids"]
            else:
                self.ko_ids = None
        return self.ko_ids is not None


    def meets_alignment_requirements(self, project_type, aln_len_ratio):
        if project_type == "isolate":
            q_aln_len = self.query_gene_end - self.query_gene_start + 1
            s_aln_len = self.subject_gene_end - self.subject_gene_start + 1
            if (q_aln_len / self.query_gene_len < aln_len_ratio or
                    s_aln_len / self.subject_gene_len < aln_len_ratio):
                return False
            self.alignment_length = q_aln_len
        else:
            if self.query_gene_len <= self.subject_gene_len:
                gene_len = self.query_gene_len
                aln_len = self.query_gene_end - self.query_gene_start + 1
            else:
                gene_len = self.subject_gene_len
                aln_len = self.subject_gene_end - self.subject_gene_start + 1
            if aln_len / gene_len < aln_len_ratio:
                return False
            self.alignment_length = (self.query_gene_end -
                                     self.query_gene_start + 1)
        return True


    def ko_ec_out(self, ko_file_handle, ec_file_handle,
                  gff_file_handle, md5_lookup, tool_version):
        ko_ec_fields = [self.query_gene]
        ko_ec_fields.append(md5_lookup[self.subject_gene]["gene_ids"][0])
        ko_ec_fields.append("")
        ko_ec_fields.append(self.percent_id)
        ko_ec_fields.append(str(self.query_gene_start))
        ko_ec_fields.append(str(self.query_gene_end))
        ko_ec_fields.append(str(self.subject_gene_start))
        ko_ec_fields.append(str(self.subject_gene_end))
        ko_ec_fields.append(self.evalue)
        ko_ec_fields.append(self.bitscore)
        ko_ec_fields.append(str(self.alignment_length))
        gff_fields = []
        gff_fields.append(self.query_gene)
        gff_fields.append(tool_version)
        gff_fields.append("KO:" + self.ko_ids.replace(":", "_KO:"))
        for ko_id in self.ko_ids.split(":"):
            ko_ec_fields[2] = "KO:" + ko_id
            ko_file_handle.write("\t".join(ko_ec_fields) + "\n")
        if "ec_ids" in md5_lookup[self.subject_gene]:
            gff_fields[-1] += ("__EC:" +
                               md5_lookup[self.subject_gene]["ec_ids"].replace(":",
                                                                               "_EC:"))
            for ec_id in md5_lookup[self.subject_gene]["ec_ids"].split(":"):
                ko_ec_fields[2] = "EC:" + ec_id
                ec_file_handle.write("\t".join(ko_ec_fields) + "\n")
        gff_fields.append(str(self.query_gene_start))
        gff_fields.append(str(self.query_gene_end))
        gff_fields.append(self.bitscore)
        gff_fields.append(".\t.")
        gff_fields.append("ID=" + self.query_gene + "_" +
                          str(self.query_gene_start) + "_" +
                          str(self.query_gene_end) + ";subject_gene_ids=" +
                          ",".join(md5_lookup[self.subject_gene]["gene_ids"]) +
                          ";subject_start=" + str(self.subject_gene_start) +
                          ";subject_end=" + str(self.subject_gene_end) +
                          ";evalue=" + self.evalue + ";percent_identity=" +
                          self.percent_id + ";alignment_length=" +
                          str(self.alignment_length) + ";query_gene_length=" +
                          str(self.query_gene_len) + ";subject_gene_length=" +
                          str(self.subject_gene_len))
        gff_file_handle.write("\t".join(gff_fields) + "\n")


    def phylo_out(self, file_handle, md5_lookup, phylogeny_lookup):
        taxon_id = md5_lookup[self.subject_gene]["taxon_id"]
        fields = [self.query_gene]
        fields.append(md5_lookup[self.subject_gene]["gene_ids"][0])
        fields.append(taxon_id)
        fields.append(self.percent_id)
        fields.append(phylogeny_lookup[taxon_id])
        file_handle.write("\t".join(fields) + "\n")



def create_md5_lookup():
    print(str(datetime.now()) + " - Creating MD5 lookup...")
    md5_lookup = {}
    for line in args.md5_mapping_file:
        line_split = line.rstrip().split("\t")
        md5_hash = line_split[0]
        md5_lookup[md5_hash] = {}
        md5_lookup[md5_hash]["gene_ids"] = line_split[1].split(":")
        if line_split[2] != "N/A":
            md5_lookup[md5_hash]["ko_ids"] = line_split[2]
        if line_split[3] != "N/A":
            md5_lookup[md5_hash]["ec_ids"] = line_split[3]
        md5_lookup[md5_hash]["taxon_id"] = line_split[4]

    return md5_lookup


def create_phylogeny_lookup():
    print(str(datetime.now()) + " - Creating phylogeny lookup...")
    phylogeny_lookup = {}
    for line in args.phylo_mapping_file:
        taxon_id, phylogeny = line.rstrip().split("\t")
        phylogeny_lookup[taxon_id] = phylogeny

    return phylogeny_lookup


def get_tool_and_version():
    last_version_line = args.last_blasttabplus_file.readline()
    if not last_version_line.startswith("# LAST version"):
        print("""ERROR: First line of blastab file does not contain the version
              information!""", file = sys.stderr)
        sys._exit(1)
    version = last_version_line.rstrip().split()[3]
    return "lastal " + version


def get_hit_to_report(ko_hits):
    for i in range(1, len(ko_hits)):
        if ko_hits[0].ko_ids != ko_hits[i].ko_ids:
            return None
    for hit in ko_hits:
        if hit.meets_alignment_requirements(args.project_type,
                                            args.aln_len_ratio):
            return hit

    return None


md5_lookup = create_md5_lookup()
phylogeny_lookup = create_phylogeny_lookup()
tool_version = get_tool_and_version()

current_gene = ""
best_ko_hit = None
ko_hits = []
hit_counter = 0
skip = False
print(str(datetime.now()) + " - Running over lastal's blasttab+ output now...")
for line in args.last_blasttabplus_file:
    if line.startswith("#"): # Ignore comment lines
        continue

    hit = Hit(line.rstrip().split())

    if hit.query_gene != current_gene:
        # New gene, check if there's a hit to report for previous gene
        if len(ko_hits) >= args.min_ko_hits:
            hit_to_report = get_hit_to_report(ko_hits)
            if hit_to_report:
                hit_to_report.ko_ec_out(args.ko_output_file,
                                        args.ec_output_file,
                                        args.ko_ec_gff_output_file,
                                        md5_lookup, tool_version)

        # Setup for new gene and checks
        current_gene = hit.query_gene
        ko_hits = []
        hit_counter = 0
        skip = False
        # Report phylogeny based on top hit
        hit.phylo_out(args.phylo_output_file, md5_lookup, phylogeny_lookup)

    if skip:
        continue

    hit_counter += 1

    if hit_counter > args.number_of_top_hits:
        skip = True
        continue

    if hit.has_ko(md5_lookup):
        ko_hits.append(hit)

if len(ko_hits) >= args.min_ko_hits:
    hit_to_report = get_hit_to_report(ko_hits)
    if hit_to_report:
        hit_to_report.ko_ec_out(args.ko_output_file, args.ec_output_file,
                                args.ko_ec_gff_output_file, md5_lookup,
                                tool_version)


args.ko_output_file.close()
args.ec_output_file.close()
args.ko_ec_gff_output_file.close()
args.phylo_output_file.close()
print(str(datetime.now()) + " - Done.")
# _exit() gets called to avoid the long running garbage collection.
os._exit(0)

