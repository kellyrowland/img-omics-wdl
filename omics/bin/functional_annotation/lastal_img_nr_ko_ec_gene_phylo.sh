#!/usr/bin/env bash

if [[ $# -lt 5 || $# -gt 9 ]]
then
	printf "Usage: $(basename $0) <project_type> <proteins_fasta> " >&2
	printf "<img_nr_last_db> <md5_lookup> <phylogeny_lookup> " >&2
	printf "<number_of_threads (default: 0, uses all available cores)> " >&2
	printf "<number_of_top_hits_to_consider (default 5)> " >&2
	printf "<min_ko_hits (default: 2)> " >&2
	printf "<aln_length_ratio (default 0.7)>\n" >&2
	exit 1
fi


project_type=$1
if [[ $project_type != "isolate" && $project_type != "metagenome" ]]
then
	printf "Unknown project type: $project_type" >&2
	exit 1
fi

proteins_fasta=$2
if [[ ! -f $proteins_fasta ]]
then
	echo "$proteins_fasta is not a file! Aborting!" >&2
	exit 1
fi

img_nr_last_db=$3

md5_lookup=$4

phylogeny_lookup=$5

number_of_threads=0
if [[ ! -z "$6" ]]
then
	number_of_threads=$6
fi

number_of_top_hits_to_consider=5
if [[ ! -z "$7" ]]
then
	no_assignment_cutoff=$7
fi

min_ko_hits=2
if [[ ! -z "$8" ]]
then
	min_ko_hits=$8
fi

aln_length_ratio=0.7
if [[ ! -z "$9" ]]
then
	aln_length_ratio=$9
fi

last_blasttab=${proteins_fasta%.*}.img_nr.last.blasttab
echo "$(date +%F_%T) - Running lastal against IMG-NR db now..."
lastal -f blasttab+ -P $number_of_threads $img_nr_last_db \
			$proteins_fasta 1> $last_blasttab
exit_code=$?
if [[ $exit_code -ne 0 ]]
then
	echo "$(date +%F_%T) - lastal failed! Aborting!" >&2
	exit $exit_code
fi


ko_tsv_out=${proteins_fasta%_*}_ko.tsv
ec_tsv_out=${proteins_fasta%_*}_ec.tsv
ko_ec_gff_out=${proteins_fasta%_*}_ko_ec.gff
phylo_tsv_out=${proteins_fasta%_*}_gene_phylogeny.tsv
echo "$(date +%F_%T) - Selecting KOs, EC numbers and gene phylogeny now..."
$(dirname $0)/lastal_img_nr_ko_ec_gene_phylo_hit_selector.py \
				-l $aln_length_ratio \
				-m $min_ko_hits \
				-n $number_of_top_hits_to_consider \
				$project_type \
				$md5_lookup \
				$phylogeny_lookup \
				$ko_tsv_out \
				$ec_tsv_out \
				$phylo_tsv_out \
				$ko_ec_gff_out < $last_blasttab
exit_code=$?
if [[ $exit_code -ne 0 ]]
then
	printf "$(date +%F_%T) - " >&2
	printf "lastal_img_nr_ko_ec_hit_selector.py failed! Aborting!\n" >&2
	exit $exit_code
fi

echo "$(date +%F_%T) - Done with KO, EC and gene phylogeny annotation."

