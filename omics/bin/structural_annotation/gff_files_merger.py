#!/usr/bin/env python3

import argparse
import sys
import os

parser = argparse.ArgumentParser(description="""This script merges two or more
                                 GFF files containing gene predictions.

                                 The order in which the paths to the GFF files
                                 are given as arguments also determines the
                                 priority in which their genes will be
                                 accepted. Given that all predictions in the
                                 first GFF will be kept, predictions within the
                                 subsequent GFF files will only be accepted if
                                 they don't overlap with any of the already
                                 stored predictions in a conflicting way.
                                 If additionally a GFF file is given via the -a
                                 option, all predictions in it will be kept, no
                                 matter if they overlap with other genes or
                                 not.

                                 The final (merged) GFF will get printed to
                                 stdout.""")
parser.add_argument("main_gff", metavar="highest_priority_gff",
                    help="file path to the gff with highest priority")
parser.add_argument("other_gffs", metavar="lower_priority_gff", nargs="+",
                    help="file path(s) to lower priority gff(s)")
parser.add_argument("-a", "--allowed_to_overlap_gff",
                    metavar="allowed_to_overlap_gff",
                    help="""file path to gff containing genes that are allowed
                    to overlap""")
parser.add_argument("-f", "--fna_file",
                    metavar="fna_file", type=argparse.FileType('r'),
                    help="""file path to the fasta file containing the
                    nucleotide sequences of the in the GFF file listed
                    contigs""")
parser.add_argument("-m", "--min_cds_len",
                    metavar="min_cds_len", default=75,
                    help="""minimum length a CDS must have [default: 75]""")
args = parser.parse_args()



class Alternate_Start_Finder:
    def __init__(self, fna_file_handle):
        self.fh = fna_file_handle
        self.sequence_lookup = {}
        self.seq_name_2_fh_pos = {}
        self.start_codons = frozenset(["ATG", "GTG", "TTG"])
        self.complements = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


    def read_contig_seq(self):
        seq = ""
        """tell() is disabled in a for loop (when next() is invoked).
        So we have to use while and readline()."""
        line = self.fh.readline()
        while not line.startswith(">"):
            if line == "":
                # Reached end of file
#                print("EOF, returning seq of len " + str(len(seq)),
#                      file=sys.stderr)
                return seq
            seq += line.rstrip()
            line = self.fh.readline()
        # We just read a contig name line
        self.fh.seek(self.fh.tell() - len(line))
        return seq


    def get_nucleotide_sequence(self, contig_name):
#        print("Trying to get sequence for contig: " + contig_name,
#              file=sys.stderr)
        """Check, if we pulled the sequence already."""
        if contig_name in self.sequence_lookup:
#            print("Found already stored sequence for contig: " + contig_name,
#                  file=sys.stderr)
            return self.sequence_lookup[contig_name]
        """Otherwise check, if we know where the sequence starts."""
        if contig_name in self.seq_name_2_fh_pos:
#            print("Found stored offset for contig: " + contig_name,
#                  file=sys.stderr)
            self.fh.seek(self.seq_name_2_fh_pos[contig_name])
            self.sequence_lookup[contig_name] = self.read_contig_seq()
            return self.sequence_lookup[contig_name]
        """If not, look for the sequence and remember all passed sequences and
        their positions in the file. (See above tell() comment.)"""
        line = self.fh.readline()
        while not line.startswith(">" + contig_name):
            if line == "":
                # We are at end of file
                print("Reached end of file and could not find contig: " +
                      contig_name + "\nAborting!", file=sys.stderr)
                os._exit(1)
            if line.startswith(">"):
                seq_name = line.rstrip().split("\s", 1)[0][1:]
#                print("Storing offset for contig: " +seq_name, file=sys.stderr)
                self.seq_name_2_fh_pos[seq_name] = self.fh.tell()
            line = self.fh.readline()
        if line.startswith(">" + contig_name):
            self.sequence_lookup[contig_name] = self.read_contig_seq()
            return self.sequence_lookup[contig_name]
        return None


    def get_alternate_start(self, contig_name, min_start, max_start, strand):
        alternate_start = -1
#        print("      min_start: " + str(min_start) + " , max_start: " +
#              str(max_start), file = sys.stderr)
        if min_start >= max_start:
            return alternate_start
        modulo_remainder = (max_start - min_start) % 3
#        print("      Modulo remainder: " + str(modulo_remainder),
#              file = sys.stderr)
        if modulo_remainder > 0:
            if strand == "+":
                min_start += modulo_remainder
            else:
                max_start -= modulo_remainder
        slice_start = min_start - 1 # -1 due to 0-based index
        slice_end = max_start
        if strand == "+":
            slice_end += 2 # To include codon starting at last possible start
        else:
            slice_start -= 2 # To include last possible codon

        contig = self.get_nucleotide_sequence(contig_name)
        searchspace = contig[slice_start:slice_end]
        if strand == "-": # Build reverse complement
            searchspace = "".join([self.complements[nt] for nt in
                                   searchspace[::-1]])
#        print("      s_start: " + str(slice_start) + " , s_end: " +
#              str(slice_end) + "\n      searchspace:\n" +
#              "\n".join([searchspace[0+i:i+75] for i in range(0,
#                                                              len(searchspace),
#                                                              75)]),
#              file=sys.stderr)
        for start_offset in range(0, len(searchspace)-2, 3):
            codon = searchspace[start_offset:start_offset+3]
            if codon in self.start_codons:
#                print("      Found " + codon + " at offset " +
#                      str(start_offset), file=sys.stderr)
                if strand == "+":
                    alternate_start = min_start + start_offset
                else:
                    alternate_start = max_start - start_offset
                break
        return alternate_start


# For debugging:
#    def get_original_gene_seq(self, contig_name, start, end, strand):
#        contig = self.get_nucleotide_sequence(contig_name)
#        gene_seq = contig[start-1:end]
#        if strand == "-":
#            gene_seq = "".join([self.complements[nt] for nt in gene_seq[::-1]])
#        return "\n".join([gene_seq[0+i:i+75]
#                          for i in range(0, len(gene_seq), 75)])



alternate_start_finder = None
if args.fna_file:
    alternate_start_finder = Alternate_Start_Finder(args.fna_file)



class Gene:
    def __init__(self, fields):
        self.seq_id = fields[0]
        self.source = fields[1]
        self.type = fields[2]
        self.start = int(fields[3])
        self.end = int(fields[4])
        self.score = fields[5]
        self.strand = fields[6]
        self.phase = fields[7]
        self.attributes = fields[8].rstrip()
# For debugging:
#        self.id = (self.seq_id + "_" + str(self.start) + "_" + str(self.end) +
#                   " (" + self.strand + " , " + self.type + " , " +
#                   self.source + ")")
#        print("\nCreated: " + self.id, file=sys.stderr)


    def add_child(self, gene):
        if not hasattr(self, "childs"):
            self.childs = []
        self.childs.append(gene)


    def get_overlap_length(self, gene):
        # First check if the overlaps overlap at all
        if self.start > gene.end or self.end < gene.start:
            return 0

        overlap_start = max(self.start, gene.start)
        overlap_end = min(self.end, gene.end)
        return overlap_end - overlap_start + 1


    def got_shortened(self, gene, allowed_overlap):
#        print("    Trying to shorten the following gene now: " + self.id,
#              file=sys.stderr)
        if self.strand == "+" and self.end > gene.end:
            min_start = gene.end - allowed_overlap + 1
            max_start = self.end - args.min_cds_len + 1
            alternate_start = alternate_start_finder.get_alternate_start(self.seq_id,
                                                                         min_start,
                                                                         max_start,
                                                                         self.strand)
            if alternate_start > -1:
#                print("The original gene sequence was:\n" +
#                      alternate_start_finder.get_original_gene_seq(self.seq_id,
#                                                                   self.start,
#                                                                   self.end,
#                                                                   self.strand),
#                      file=sys.stderr)
                attribute_fields = self.attributes.split(";")
                attribute_fields.append("shortened=original start " +
                                        str(self.start))
                self.start = alternate_start
                new_gene_id = (self.seq_id + "_" + str(self.start) + "_" +
                               str(self.end))
                attribute_fields[0] = "ID=" + new_gene_id
                self.attributes = ";".join(attribute_fields)
                self.end_shortened_during_current_iteration = True
#                print("      SHORTENED " + self.id + " (new start is " +
#                      str(alternate_start) + ")", file=sys.stderr)
#                print("      New Gene ID is: " + new_gene_id, file=sys.stderr)
#                id_split = self.id.split()
#                id_split[0] = new_gene_id
#                self.id = " ".join(id_split)
                return True
        elif self.strand == "-" and self.start < gene.start:
            max_end = self.start + args.min_cds_len - 1
            min_start = gene.start + allowed_overlap - 1
            alternate_start = alternate_start_finder.get_alternate_start(self.seq_id,
                                                                         max_end,
                                                                         min_start,
                                                                         self.strand)
            if alternate_start > -1:
#                print("The original gene sequence was:\n" +
#                      alternate_start_finder.get_original_gene_seq(self.seq_id,
#                                                                   self.start,
#                                                                   self.end,
#                                                                   self.strand),
#                      file=sys.stderr)
#                print("      SHORTENED " + self.id + " (new end is " +
#                      str(alternate_start) + ")\n", file=sys.stderr)
                attribute_fields = self.attributes.split(";")
                attribute_fields.append("shortened=original end " + str(self.end))
                self.end = alternate_start
                new_gene_id = (self.seq_id + "_" + str(self.start) + "_" +
                               str(self.end))
                attribute_fields[0] = "ID=" + new_gene_id
                self.attributes = ";".join(attribute_fields)
                self.end_shortened_during_current_iteration = True
#                print("      New Gene ID is: " + new_gene_id, file=sys.stderr)
#                id_split = self.id.split()
#                id_split[0] = new_gene_id
#                self.id = " ".join(id_split)
                return True
#        print("      Not able to shorten gene! Ignoring it!", file=sys.stderr)
        return False


    def has_no_inacceptable_overlap(self, gene, allowed_overlap_lookup):
#        print("Comparing: " + self.id + " and " + gene.id, file=sys.stderr)
        allowed_overlap = min(allowed_overlap_lookup[self.type],
                              allowed_overlap_lookup[gene.type])
        overlap = self.get_overlap_length(gene)
#        print("    Overlap: " + str(overlap) + " (allowed: " + str(allowed_overlap)
#              + ")", file=sys.stderr)
        if overlap > allowed_overlap:
#            print("    OVERLAP is over the threshold!", file=sys.stderr)
            """Check, if this gene is a CDS and can be shortened."""
            if self.type == "CDS":
                """ CDSs from the same method are always allowed to overlap. """
                if gene.type == "CDS" and self.source == gene.source:
#                    print("    Is CDS - CDS from same method. Allowed!",
#                          file=sys.stderr)
                    return True
                if alternate_start_finder is not None:
                    return self.got_shortened(gene, allowed_overlap)
            return False
        else:
            """ Make sure we do not add an identical CDS gene
                that's within the allowed overlap cutoff."""
            if (overlap > 0 and self.type == "CDS" and gene.type == "CDS"
                    and self.start == gene.start and self.end == gene.end):
#                print("    The CDSs are identical. Ignoring " + self.id,
#                      file=sys.stderr)
                return False
#        print("    HAS _NO_ INACCEPTABLE OVERLAP!", file=sys.stderr)
        return True


    def to_gff(self):
        print(self.seq_id + "\t" + self.source + "\t" + self.type + "\t" +
              str(self.start) + "\t" + str(self.end) + "\t" + self.score +
              "\t" + self.strand + "\t" + self.phase + "\t" + self.attributes)
        if hasattr(self, "childs"):
            for gene in self.childs:
                gene.to_gff()



class Contig_GFF:
    def __init__(self):
        self.genes = []
        self.last_used_idx = 0


    def reset_index_tracker(self):
        self.last_used_idx = 0


    def add_gene(self, gene, allowed_overlap_lookup,
                 allowed_to_significantly_overlap=False):
        """If it's a CDS, make sure it has the required minimum length."""
        if (gene.type == "CDS" and
              (gene.end - gene.start + 1) < args.min_cds_len):
            return False

        """If there are no genes stored yet, just add the gene."""
        if not self.genes:
            self.genes.append(gene)
            return True

        """Otherwise only add this gene if it is not conflicting with any of
        the already accepted genes.
        First check, if the gene can get added at the end."""
        if gene.end > self.genes[-1].end:
#            print("New gene has > end than last gene in list.",
#                  file=sys.stderr)
            if allowed_to_significantly_overlap:
                self.genes.append(gene)
                self.last_used_idx = len(self.genes) - 1
                return True

            if (gene.has_no_inacceptable_overlap(self.genes[-1],
                                                   allowed_overlap_lookup)):
                if (hasattr(gene, "end_shortened_during_current_iteration") and
                        gene.end_shortened_during_current_iteration == True):
                    self.last_used_idx -= 5
                    if self.last_used_idx < 0:
                        self.last_used_idx = 0
                    gene.end_shortened_during_current_iteration = False
                else:
                    self.genes.append(gene)
                    self.last_used_idx = len(self.genes) - 1
                    return True
            else:
                return False

        """If not, check if the gene fits between two genes within our overlap
        costraints."""
        loop_boundary = len(self.genes) #- 1
#        print("Loop boundary: " + str(loop_boundary), file=sys.stderr)
        while self.last_used_idx < loop_boundary:
#            print("While iteration. LUIDX: " + str(self.last_used_idx) +
#                  "\nGene end: " +
#                  str(gene.end) + " , index gene end: " +
#                  str(self.genes[self.last_used_idx].end), file=sys.stderr)
            if gene.end <= self.genes[self.last_used_idx].end:
                if allowed_to_significantly_overlap:
#                    print("ALLOWED TO OVERLAP. Instering at " +
#                          str(self.last_used_idx), file=sys.stderr)
                    self.genes.insert(self.last_used_idx, gene)
                    return True

                if gene.has_no_inacceptable_overlap(self.genes[self.last_used_idx],
                                                    allowed_overlap_lookup):
#                    print("    No ina overlap with gene at idx: " +
#                          str(self.last_used_idx) +
#                          "\n    Checking with previous now (if present).",
#                          file=sys.stderr)
                    if (hasattr(gene, "end_shortened_during_current_iteration")
                          and
                          gene.end_shortened_during_current_iteration == True):
                        self.last_used_idx -= 15
                        if self.last_used_idx < 0:
                            self.last_used_idx = 0
                        gene.end_shortened_during_current_iteration = False
#                        print("    The gene got shortened and has a new end " +
#                              "coordinate. Resetting LUIDX to: " +
#                              str(self.last_used_idx) + "\n    Resetting end " +
#                              "shortened in current iteration marker to: " +
#                              "False\n    Going to recompare against " +
#                              "previous 15 genes now.", file=sys.stderr)
                        continue
                    if (self.last_used_idx == 0
                          or
                          gene.has_no_inacceptable_overlap(self.genes[self.last_used_idx-1],
                                                           allowed_overlap_lookup)):
                        if (hasattr(gene, "end_shortened_during_current_iteration")
                              and
                              gene.end_shortened_during_current_iteration == True):
                            self.last_used_idx -= 15
                            if self.last_used_idx < 0:
                                self.last_used_idx = 0
                            gene.end_shortened_during_current_iteration = False
#                            print("    The gene got shortened and has a new end " +
#                                  "coordinate. Resetting LUIDX to: " +
#                                  str(self.last_used_idx) + "\n    Resetting end " +
#                                  "shortened in current iteration marker to: " +
#                                  "False\n    Going to recompare against " +
#                                  "previous 15 genes now.", file=sys.stderr)
                            continue
#                        print("    All good! Instering at " +
#                              str(self.last_used_idx), file=sys.stderr)
                        self.genes.insert(self.last_used_idx, gene)
                        return True
                return False
            self.last_used_idx += 1
#        print("    Adding gene to the end of the list.", file=sys.stderr)
        self.genes.append(gene)
#        return False
        return True


    def add_child_to_last_stored_gene(self, gene):
        self.genes[self.last_used_idx].add_child(gene)


    def to_gff(self):
        self.genes.sort(key=lambda g: g.start)
        for gene in self.genes:
            gene.to_gff()




class GFF_Files_Merger:
    def __init__(self, gff_file):
        self.final_gff_data = {}
        self.allowed_overlap_lookup = {}
        self.allowed_overlap_lookup["rRNA"] = 10
        self.allowed_overlap_lookup["tRNA"] = 10
        self.allowed_overlap_lookup["tmRNA"] = 10
        self.allowed_overlap_lookup["ncRNA"] = 10
        self.allowed_overlap_lookup["CRISPR"] = 10
        self.allowed_overlap_lookup["CDS"] = 90
        self.contig_sequence_lookup = None
        self.previous_gene_stored = False
        self.add_genes_from_file(gff_file)


    def add_genes_from_file(self, gff_file,
                            allowed_to_significantly_overlap=False):
        fr = open(gff_file, "r")
        seq_name = ""
        self.previous_gene_stored = False
        for line in fr:
            """Ignore comments or empty lines."""
            if not line or line[0] == "#" or line == "\n":
                continue

            """Create a Gene object from the new GFF line."""
            gene = Gene(line.split("\t"))

            """Reset variables if we've a new sequence."""
            if seq_name != gene.seq_id:
                seq_name = gene.seq_id
                if seq_name not in self.final_gff_data:
                    self.final_gff_data[seq_name] = Contig_GFF()
                else:
                    self.final_gff_data[seq_name].reset_index_tracker()
                self.previous_gene_stored = False

            """If we've a child line, add it to the previous gene (if we stored
            it)."""
            if "Parent=" in gene.attributes:
                if self.previous_gene_stored:
                    self.final_gff_data[seq_name].add_child_to_last_stored_gene(gene)
                continue

#            print("Seq: " + seq_name + " , LUIDX: " +
#                  str(self.final_gff_data[seq_name].last_used_idx),
#                  file=sys.stderr)
            """Add the gene if possible."""
            self.previous_gene_stored = self.final_gff_data[seq_name].add_gene(gene,
                                                        self.allowed_overlap_lookup,
                                                        allowed_to_significantly_overlap)
#            print(gene.id + " got stored: " + str(self.previous_gene_stored) +
#                  "\n", file=sys.stderr)
        fr.close()


    def print_final_gff(self):
        for seq_name in sorted(self.final_gff_data):
            self.final_gff_data[seq_name].to_gff()



#print("INSERTING " + args.main_gff, file=sys.stderr)
gff_files_merger = GFF_Files_Merger(args.main_gff)
"""Now add the remaining GFFs."""
for gff_file in args.other_gffs:
#    print("\nINSERTING " + gff_file, file=sys.stderr)
    gff_files_merger.add_genes_from_file(gff_file)

"""Add genes that can overlap completely with everything."""
if args.allowed_to_overlap_gff:
    allowed_to_significantly_overlap = True
#    print("\nINSERTING ALLOWED TO OVERLAP FILE " + args.allowed_to_overlap_gff,
#         file=sys.stderr)
    gff_files_merger.add_genes_from_file(args.allowed_to_overlap_gff,
                                         allowed_to_significantly_overlap)

"""Print out the final GFF"""
gff_files_merger.print_final_gff()

