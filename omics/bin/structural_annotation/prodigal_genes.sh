#!/usr/bin/env bash

if [[ $# -lt 2 || $# -gt 3 ]]
then
	printf "Usage: $(basename $0) <input_fasta> " >&2
	printf "<project_type (isolate|metagenome)> " >&2
	printf "<transation_table (default: 11)>\n" >&2
	exit 1
fi

input_fasta=$1
if [[ ! -f $input_fasta ]]
then
	echo "$input_fasta is not a file! Aborting!" >&2
	exit 1
fi

project_type=$2
if [[ $project_type != "isolate" && $project_type != "metagenome" ]]
then
	echo "Unknown project type: $project_type" >&2
	echo "Allowed types are 'isolate' and 'metagenome'. Aborting!" >&2
    exit 1
fi

translation_table=11
if [[ ! -z "$3" ]]
then
	translation_table=$3
fi

log=${input_fasta%_*}_prodigal.out
gff_out=${input_fasta%_*}_prodigal.gff
genes_fasta=${input_fasta%_*}_prodigal_genes.fna
proteins_fasta=${input_fasta%_*}_prodigal_proteins.faa

if [[ $project_type == "isolate" && $(grep -v '^>' $input_fasta | wc -m) -ge 20000 ]]
then
	training_file=${input_fasta%_*}_prodigal.trn
	prodigal -i $input_fasta -t $training_file -g \
			 $translation_table -q &>> $log
	exit_code=$?
	if [[ $exit_code -eq 0 ]]
	then
		echo "$(date +%F_%T) - Created training file: $training_file"
	else
		echo "Model creation failed! Aborting!" >&2
		exit $exit_code
	fi

	echo "$(date +%F_%T) - Running Prodigal now..."
	prodigal -f gff -g $translation_table -p single -m \
			 -i $input_fasta -t $training_file \
			 -o $gff_out -d $genes_fasta -a $proteins_fasta &>> $log
	if [[ $? -ne 0 ]]
	then
		echo "Prodigal failed! Aborting!" >&2
		exit $?
	fi
else
	echo "$(date +%F_%T) - Running Prodigal now..."
	prodigal -f gff -p meta -m \
			 -i $input_fasta -o $gff_out \
			 -d $genes_fasta -a $proteins_fasta &>> $log
	exit_code=$?
	if [[ $exit_code -ne 0 ]]
	then
		echo "Prodigal failed! Aborting!" >&2
		exit $exit_code
	fi
fi

# Remove trailing * of protein sequences and switch other * to X
echo "$(date +%F_%T) - Removing trailing * in protein seqs now..."
sed -i 's/\*$//g' $proteins_fasta
echo "$(date +%F_%T) - Switching all other * to X now..."
sed -i 's/\*/X/g' $proteins_fasta

# Fix gene IDs
echo "$(date +%F_%T) - Running unify_gene_ids.py now..."
$(dirname $0)/unify_gene_ids.py $gff_out $genes_fasta $proteins_fasta
exit_code=$?
if [[ $exit_code -ne 0 ]]
then
	echo "unify_gene_ids.py failed! Aborting!" >&2
	exit $exit_code
fi

echo "$(date +%F_%T) - Done with Prodigal module."
