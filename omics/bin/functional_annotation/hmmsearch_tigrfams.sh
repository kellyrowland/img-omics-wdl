#!/usr/bin/env bash

if [[ $# -lt 2 || $# -gt 5 ]]
then
	printf "Usage: $(basename $0) <proteins_fasta> <tigrfam_hmm_db> " >&2
	printf "<number_of_additional_threads (default: 0)>\n" >&2
	printf "<min_aln_length_ratio (default 0.7)> " >&2
	printf "<max_overlap_ratio (default 0.1)> " >&2
	exit 1
fi


proteins_fasta=$1
if [[ ! -f $proteins_fasta ]]
then
	echo "$proteins_fasta is not a file! Aborting!" >&2
	exit 1
fi

tigrfam_hmm_db=$2
if [[ ! -f $tigrfam_hmm_db ]]
then
	echo "$tigrfam_hmm_db is not a file! Aborting!" >&2
	exit 1
fi

number_of_additional_threads=0
if [[ ! -z "$3" ]]
then
	number_of_additional_threads=$3
fi

min_aln_length_ratio=0.7
if [[ ! -z "$4" ]]
then
	min_aln_length_ratio=$4
fi

max_overlap_ratio=0.1
if [[ ! -z "$5" ]]
then
	max_overlap_ratio=$5
fi


echo "$(date +%F_%T) - Calling hmmsearch to predict TIGRFAMs now..."
hmmsearch --notextw --cut_nc --cpu $number_of_additional_threads \
			--domtblout ${proteins_fasta%.*}.tigrfam.domtblout \
    		$tigrfam_hmm_db $proteins_fasta 1> /dev/null
exit_code=$?
if [[ $exit_code -ne 0 ]]
then
	echo "$(date +%F_%T) - hmmsearch failed! Aborting!" >&2
	exit $exit_code
fi
echo "$(date +%F_%T) - Getting hmmsearch version now..."
tool_and_version=$(hmmsearch -h | grep HMMER)
tool_and_version=${tool_and_version%%;*}
tool_and_version=${tool_and_version:2}

echo "$(date +%F_%T) - Sort, filter and select TIGRFAMs now..."
tigrfams_gff=${proteins_fasta%_*}_tigrfam.gff
grep -v '^#' ${proteins_fasta%.*}.tigrfam.domtblout | \
awk '{print $1,$3,$4,$6,$13,$14,$16,$17,$20,$21}' | \
sort -k1,1 -k6,6nr -k5,5n | \
$(dirname $0)/hmmsearch_hit_selector.py -a $min_aln_length_ratio \
										-o $max_overlap_ratio \
										"$hmmer_version" 1> $tigrfams_gff
exit_code=$?
if [[ $exit_code -ne 0 ]]
then
	echo "$(date +%F_%T) - One of the pipe commands failed! Aborting!" >&2
	exit $exit_code
fi

echo "$(date +%F_%T) - Done with TIGRFAM annotation."
