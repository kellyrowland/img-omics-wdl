#!/usr/bin/env bash

if [[ $# -ne 2 ]]
then
	printf "Usage: $(basename $0) <input_fasta> " >&2
	printf "<project_type (isolate|metagenome)>\n" >&2
	exit 1
fi


input_fasta=$1
if [[ ! -f $input_fasta ]]
then
	echo "$input_fasta is not a file! Aborting!" >&2
	exit 1
fi


project_type=$2


gff_out=${input_fasta%_*}_genemark.gff
genes_fasta=${input_fasta%_*}_genemark_genes.fna
proteins_fasta=${input_fasta%_*}_genemark_proteins.faa


if [[ $project_type == "isolate" ]]
then
	echo "$(date +%F_%T) - Running gmsn2.pl now..."
	gms2.pl --seq $input_fasta --genome-type auto --output $gff_out \
			--format gff --fnn $genes_fasta --faa $proteins_fasta
	exit_code=$?
	if [[ $exit_code -eq 0 ]]
	then
		mv log gms2.log
		rm group*
	else
		echo "gms2.pl failed! Aborting!." >&2
		exit $exit_code
	fi
elif [[ $project_type == "metagenome" ]]
then
	model_file="$(dirname $(readlink -f $(which gmhmmp2)))/mgm_11.mod"
	printf "$(date +%F_%T) - Will run gmhmmp2 with model file "
	printf "$(basename $model_file) now...\n"
	gmhmmp2 --Meta $model_file --incomplete_at_gaps 30 -o $gff_out \
			--format gff --NT $genes_fasta --AA $proteins_fasta \
			--seq $input_fasta
	exit_code=$?
	if [[ $exit_code -ne 0 ]]
	then
		echo "gmhmmp2 failed! Aborting!." >&2
		exit $exit_code
	fi
else
	echo "Unknown project type: $project_type" >&2
	echo "Valid project types are: 'isolate' and 'metagenome'" >&2
	exit 1
fi


# Switch * occurrences to X
echo "$(date +%F_%T) - Switching * to X now..."
sed -i 's/\*/X/g' $proteins_fasta


# Fix genemark version and gene names
echo "$(date +%F_%T) - Running unify_gene_ids.py now..."
$(dirname $0)/unify_gene_ids.py $gff_out $genes_fasta $proteins_fasta
exit_code=$?
if [[ $exit_code -ne 0 ]]
then
	echo "unify_gene_ids.py failed! Aborting!" >&2
	exit $exit_code
fi


echo "$(date +%F_%T) - Done with GeneMark module."
