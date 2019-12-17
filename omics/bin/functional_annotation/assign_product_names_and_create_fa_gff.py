#!/usr/bin/env python3

import argparse
import fileinput
import os
#import sys


parser = argparse.ArgumentParser(description="""This script will assign names
                                 to the given proteins.
                                 For the assigment the script looks at the
                                 functional annotation results in the following
                                 order:

                                 KO > TIGRfam > COG > Pfam

                                 If none of the above results have a hit for a
                                 given gene it'll be annotated with
                                 'hypothetical protein'.
                                 If a gene got assigned multiple IDs of the
                                 same type (KO, TIGRfam, COG or Pfam) then the
                                 corresponding terms will get concatenated to
                                 the product name via '/'. If a gene has
                                 multiple hits to the same ID its corresponding
                                 term will only get listed once.""")
parser.add_argument("mapping_dir",
                    help="""path to directory that contains the product names
                    mapping files""")
parser.add_argument("structural_annotation_gff",
                    help="""path to the GFF file containing the structural
                    annotations""")
fa_files = parser.add_argument_group("GFFs also used for product name assignments")
fa_files.add_argument("-k", "--ko",
                    metavar="ko_file", type=argparse.FileType('r'),
                    help="path to KO hits file")
fa_files.add_argument("-t", "--tigrfam", metavar="tigrfam_file",
                    type=argparse.FileType('r'),
                    help="path to TIGRfam hits file")
fa_files.add_argument("-c", "--cog", metavar="cog_file",
                    type=argparse.FileType('r'),
                    help="path to COG hits file")
fa_files.add_argument("-p", "--pfam", metavar="pfam_file",
                    type=argparse.FileType('r'),
                    help="path to Pfam hits file")
additional_fa_files = parser.add_argument_group("GFFs not used for the product name assignments")
additional_fa_files.add_argument("-f", "--cath_funfam",
                                 metavar="cath_funfam_file",
                                 type=argparse.FileType('r'),
                                 help="path to Cath-Funfam hits file")
additional_fa_files.add_argument("-s", "--smart", metavar="smart_file",
                                 type=argparse.FileType('r'),
                                 help="path to SMART hits file")
additional_fa_files.add_argument("-u", "--superfamily",
                                 metavar="superfamily_file",
                                 type=argparse.FileType('r'),
                                 help="path to SuperFamily hits file")
additional_fa_files.add_argument("-e", "--cleavage_site",
                                 metavar="cleavage_site_file",
                                 type=argparse.FileType('r'),
                                 help="path to cleavage sites file")
additional_fa_files.add_argument("-r", "--transmembrane_helix",
                                 metavar="transmembrane_helix_file",
                                 type=argparse.FileType('r'),
                                 help="path to transmembrane helices file")
parser.add_argument("-v", "--version", action="version",
                    version="%(prog)s 1.0.0")
args = parser.parse_args()
if not (args.ko or args.tigrfam or args.cog or args.pfam):
    parser.error("""No functional annotation GFF files provided. Please provide
                 at least one via the 'GFFs also used for product name
                 assignments' group arguments.""")
if not os.path.isdir(args.mapping_dir):
    parser.error(args.mapping_dir + " is not a directory.")
if not os.access(args.mapping_dir, os.R_OK):
    parser.error(args.mapping_dir + " is not readable.")


isology_ranking = {}
id_to_name = {}
tfam_to_isology = {}
gene_to_product_name = {}
gene_to_annotations = {}
product_names_fa_types = frozenset(["ko", "tigrfam", "cog", "pfam"])

def create_isology_ranking():
    isology_ranking["equivalog"] = 1
    isology_ranking["hypoth_equivalog"] = 2
    isology_ranking["paralog"] = 3
    isology_ranking["exception"] = 4
    isology_ranking["equivalog_domain"] = 5
    isology_ranking["hypoth_equivalog_domain"] = 6
    isology_ranking["paralog_domain"] = 7
    isology_ranking["subfamily"] = 8
    isology_ranking["superfamily"] = 9
    isology_ranking["subfamily_domain"] = 10
    isology_ranking["domain"] = 11
    isology_ranking["signature"] = 12
    isology_ranking["repeat"] = 13


def parse_and_add_to_dict(tsv_file, mapping_dict):
    for line in fileinput.input(tsv_file):
        line = line.rstrip()
        if line:
            fa_id, value = line.split("\t")
            mapping_dict[fa_id] = value


def create_lookups():
    tsv_files = [os.path.join(args.mapping_dir, ls_file) for ls_file in
                os.listdir(args.mapping_dir) if ls_file.endswith(".tsv")]
    for tsv_file in tsv_files:
        if tsv_file.endswith("tfam_isology.tsv"):
            parse_and_add_to_dict(tsv_file, tfam_to_isology)
        else:
            parse_and_add_to_dict(tsv_file, id_to_name)


def check_isology(product_names, sources):
    isology_mapping = {}
    for i in range(0, len(sources)):
        isology = tfam_to_isology[sources[i]]
        if isology in isology_mapping:
            isology_mapping[isology]["pns"].append(product_names[i])
            isology_mapping[isology]["srcs"].append(sources[i])
        else:
            new_product_names = [product_names[i]]
            new_sources = [sources[i]]
            isology_mapping[isology] = {}
            isology_mapping[isology]["pns"] = new_product_names
            isology_mapping[isology]["srcs"] = new_sources
    if len(isology_mapping.keys()) == 1:
        return product_names, sources
    else:
        isologies = list(isology_mapping.keys())
        best_ranking_isology = isologies[0]
        for i in range(1, len(isologies)):
            if (isology_ranking[isologies[i]] <
                isology_ranking[best_ranking_isology]):
                best_ranking_isology = isologies[i]
        return (isology_mapping[best_ranking_isology]["pns"],
                isology_mapping[best_ranking_isology]["srcs"])


def remove_ec_and_split_ids(ids):
    cleaned_ids = []
    for concat_id in ids:
        ko_ids_and_ec_numbers = concat_id.split("__")
        cleaned_ids.extend(ko_ids_and_ec_numbers[0].split("_"))
    return cleaned_ids


def make_pn_unique(product_names, sources):
    new_product_names = []
    new_sources = []
    hypothetical_present = False
    non_hypothetical_present = False
    for pn in product_names:
        if pn not in new_product_names:
            new_product_names.append(pn)
            if pn == "hypothetical protein":
                hypothetical_present = True
            else:
                non_hypothetical_present = True
    if hypothetical_present and non_hypothetical_present:
        new_product_names.remove("hypothetical protein")
    for fa_id in sources:
        if fa_id not in new_sources:
            new_sources.append(fa_id)
    return (new_product_names, new_sources)


def create_gene_mapping(gene, ids, fa_type):
    if fa_type == "ko":
        ids = remove_ec_and_split_ids(ids)
    product_names = []
    sources = []
    for fa_id in ids:
        if fa_id in id_to_name:
            product_names.append(id_to_name[fa_id])
            sources.append(fa_id)
    if len(product_names) > 0:
        if fa_type == "tigrfam":
            product_names, sources = check_isology(product_names, sources)
        product_names, sources = make_pn_unique(product_names, sources)
#        gene_to_product_name[gene] = {}
        gene_to_product_name[gene]["product_name"] = str.join("/",
                                                              product_names)
        gene_to_product_name[gene]["source"] = str.join("/", sources)


def add_gene_annotation(gene_id, fa_type, fa_id, rest_of_line):
    if gene_id not in gene_to_product_name:
        gene_to_product_name[gene_id] = {}
    if "annotations" not in gene_to_product_name[gene_id]:
        gene_to_product_name[gene_id]["annotations"] = {}
    if fa_type not in gene_to_product_name[gene_id]["annotations"]:
        if fa_type == "transmembrane_helix_parts":
            gene_to_product_name[gene_id]["annotations"][fa_type] = []
        else:
            gene_to_product_name[gene_id]["annotations"][fa_type] = set()
    if fa_type == "cleavage_site_network":
        fa_id = rest_of_line.rstrip().split("network=")[1].split(";")[0]
    if fa_type == "transmembrane_helix_parts":
        start, stop, ignore = rest_of_line.split("\t", 2)
        fa_id += "_" + start + "_" + stop
        gene_to_product_name[gene_id]["annotations"][fa_type].append(fa_id)
    else:
        gene_to_product_name[gene_id]["annotations"][fa_type].add(fa_id)
#        print("add_gene_annotation: Stored annotation = " +
#              gene_id + " - " + fa_type + " - " + fa_id, file=sys.stderr)


def parse_fa_file(file_handle, fa_type):
    prev_gene = ""
    ids = []
    for line in file_handle:
        gene_id, ignore, fa_id, rest = line.split("\t", 3)
        add_gene_annotation(gene_id, fa_type, fa_id, rest)
        if (fa_type not in product_names_fa_types
                or "product_name" in gene_to_product_name[gene_id]):
            continue
        if gene_id != prev_gene:
            create_gene_mapping(prev_gene, ids, fa_type)
            ids = []
            prev_gene = gene_id
        ids.append(fa_id)
    create_gene_mapping(prev_gene, ids, fa_type)


def process_fa_files():
    if args.ko:
        parse_fa_file(args.ko, "ko")
    if args.tigrfam:
        parse_fa_file(args.tigrfam, "tigrfam")
    if args.cog:
        parse_fa_file(args.cog, "cog")
    if args.pfam:
        parse_fa_file(args.pfam, "pfam")
    if args.cath_funfam:
        parse_fa_file(args.cath_funfam, "cath_funfam")
    if args.smart:
        parse_fa_file(args.smart, "smart")
    if args.superfamily:
        parse_fa_file(args.superfamily, "superfamily")
    if args.cleavage_site:
        parse_fa_file(args.cleavage_site, "cleavage_site_network")
    if args.transmembrane_helix:
        parse_fa_file(args.transmembrane_helix, "transmembrane_helix_parts")


def add_annotations(attributes_field, annotations_dict):
    for fa_type in sorted(list(annotations_dict.keys())):
#        print("\tadd_annotations: fa_type="+fa_type, file=sys.stderr)
        if fa_type == "ko":
            kos_and_ecs = annotations_dict["ko"].pop().split("__")
            attributes_field += (";ko=" + kos_and_ecs[0].replace("_", ","))
            if len(kos_and_ecs) == 2:
                attributes_field += (";ec_number=" +
                                     kos_and_ecs[1].replace("_", ","))
        elif fa_type == "transmembrane_helix_parts":
            attributes_field += (";" + fa_type + "=" +
                                 ",".join(annotations_dict[fa_type]))
        else:
            attributes_field += (";" + fa_type + "=" +
                                 ",".join(sorted(annotations_dict[fa_type])))
    return attributes_field


def get_product_name_and_gff_lines(fields):
    product_name_line = ""
    gene_id = fields[-1].split(";", 1)[0].split("=", 1)[1]
#    print("\nChecking " + gene_id + " (" + fields[2] + ")", file=sys.stderr)
    if fields[2] == "CDS":
        if gene_id in gene_to_product_name:
            if "product_name" in gene_to_product_name[gene_id]:
                product_name_line = (gene_id + "\t" +
                                     gene_to_product_name[gene_id]["product_name"]
                                     + "\t" +
                                     gene_to_product_name[gene_id]["source"])
                fields[-1] += (";product=" +
                               gene_to_product_name[gene_id]["product_name"] +
                               ";product_source=" +
                               gene_to_product_name[gene_id]["source"])
            else:
                product_name_line = (gene_id + "\thypothetical protein\t" +
                                     "Hypo-rule applied")
                fields[-1] += (";product=hypothetical protein" +
                               ";product_source=Hypo-rule applied")

            if "annotations" in gene_to_product_name[gene_id]:
                fields[-1] = add_annotations(fields[-1],
                                             gene_to_product_name[gene_id]["annotations"])
        else:
            product_name_line = (gene_id + "\thypothetical protein\t" +
                                 "Hypo-rule applied")
            fields[-1] += (";product=hypothetical protein" +
                           ";product_source=Hypo-rule applied")

    else:
        product = fields[-1].split("product=", 1)[1].split(";")[0]
        source = fields[2]
        if source == "rRNA":
            source += "_" + product.split(" ", 1)[0]
        product_name_line = (gene_id + "\t" + product + "\t" + source)

    gff_line = "\t".join(fields) + "\n"
    return product_name_line, gff_line



create_isology_ranking()
create_lookups()
process_fa_files()

# Create output files
functional_annotation_gff = args.structural_annotation_gff.replace("_structural_",
                                                                   "_functional_")
product_names_file = args.structural_annotation_gff[:args.structural_annotation_gff.rfind("structural")]
product_names_file += "product_names.tsv"
sa_gff_fh = open(args.structural_annotation_gff, "r")
fa_gff_fh = open(functional_annotation_gff, "w")
product_names_fh = open(product_names_file, "w")
# Run over structural annotation gff
for gff_line in sa_gff_fh:
    fields = gff_line.rstrip().split("\t")
    if fields[2] == "CDS" or fields[2].endswith("RNA"):
        product_name_line, gff_line = get_product_name_and_gff_lines(fields)
        product_names_fh.write(product_name_line + "\n")
    fa_gff_fh.write(gff_line)
sa_gff_fh.close()
fa_gff_fh.close()
product_names_fh.close()

