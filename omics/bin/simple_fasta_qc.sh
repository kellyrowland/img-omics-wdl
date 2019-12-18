#!/usr/bin/env bash

if [[ $# -lt 3 || $# -gt 7 ]]
then
	printf "Usage: $(basename $0) <project_type (isolate|metagenome)> " >&2
	printf "<input_fasta> <project_id> " >&2
	printf "[<rename (yes|no)> <n_ratio_cutoff> " >&2
	printf "<seqs_per_million_bp_cutoff> <min_seq_length>]\n" >&2
	exit 1
fi

project_type=$1
if [[ $project_type != "isolate" && $project_type != "metagenome" ]]
then
	echo "Unknown project type: $project_type! Aborting!" >&2
	exit 1
fi

input_fasta=$2
if [[ ! -f $input_fasta ]]
then
	echo "$input_fasta is not a file! Aborting!" >&2
	exit 1
fi

project_id=$3
if [[ -z "$project_id" ]]
then
	echo "No project ID given! Aborting!" >&2
	exit 1
fi

rename=0
if [[ ! -z "$4" && $4=="yes" ]]
then
	rename=1
fi

n_ratio_cutoff=0.5
if [[ ! -z "$5" ]]
then
	n_ratio_cutoff=$5
fi

seqs_per_million_bp_cutoff=500
if [[ ! -z "$6" ]]
then
	seqs_per_million_bp_cutoff=$6
fi

min_seq_length=150
if [[ ! -z "$7" ]]
then
	min_seq_length=$7
fi

tmp_fasta="${input_fasta}.tmp"
qced_fasta="$(dirname $input_fasta)/${project_id}_contigs.fna"

# Ignore empty lines
grep -v '^\s*$' $input_fasta | \
# Delete carriage return characters
tr -d '\r' | \
# Remove whitespaces between > and sequence name
sed 's/^>[[:blank:]]*/>/g' > $tmp_fasta

# Check for ratio of non-ACGT letters (which would get replaced with Ns)
acgt_count=`grep -v '^>' $tmp_fasta | grep -o [acgtACGT] | wc -l`
n_count=`grep -v '^>' $tmp_fasta | grep -o '[^acgtACGT]' | wc -l`
n_ratio=`echo $n_count $acgt_count | awk '{printf "%f", $1 / $2}'`
if (( $(echo "$n_ratio >= $n_ratio_cutoff" | bc) ))
then
	printf -v n_ratio_percent "%.2f" $(echo "$n_ratio * 100" | bc)
	printf -v n_cutoff_percent "%.2f" $(echo "$n_ratio_cutoff * 100" | bc)
	printf "Amount of non-ACGT letters ($n_ratio_percent%%) crossed " >&2
	printf "the threshold of $n_cutoff_percent%% and thus the contigs " >&2
	printf "will not get processed! Aborting!\n" >&2
	rm $tmp_fasta
	exit 1
fi

if [[ $project_type == "isolate" ]]
then
	# Make sure we don't have more than 500 sequences per million base pair
	seq_count=`grep -c '^>' $tmp_fasta`
	bp_count=`grep -v '^>' $tmp_fasta | tr -d '\n' | wc -m`
	seqs_per_million_bp=$seq_count
	if (( $bp_count > 1000000 ))
	then
		divisor=$(echo $bp_count | awk '{printf "%.f", $1 / 1000000}')
		seqs_per_million_bp=$(echo $seq_count $divisor | \
							  awk '{printf "%.2f", $1 / $2}')
	fi
	if (( $(echo "$seqs_per_million_bp > $seqs_per_million_bp_cutoff" | bc) ))
	then
		printf "The submitted file contains $seqs_per_million_bp " >&2
		printf "sequences per million base pair (or less), which is " >&2
		printf "over the threshold of $seqs_per_million_bp_cutoff " >&2
		printf "sequences per million base pair. Thus the contigs " >&2
		printf "will _not_ get processed! Aborting!\n" >&2
		rm $tmp_fasta
		exit 1
	fi
fi

# Run through sanity script
fasta_sanity_cmd="$(dirname $0)/fasta_sanity.py $tmp_fasta $qced_fasta"
if [[ rename ]]
then
	fasta_sanity_cmd="$fasta_sanity_cmd -p $project_id"
fi
fasta_sanity_cmd="$fasta_sanity_cmd -l $min_seq_length"
$fasta_sanity_cmd
exit_code=$?
if [[ $exit_code -ne 0 ]]
then
	echo "fasta_sanity.py failed! Aborting!" >&2
	exit $exit_code
fi
rm $tmp_fasta

echo $qced_fasta

