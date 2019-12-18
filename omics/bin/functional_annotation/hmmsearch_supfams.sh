#!/usr/bin/env bash

if [[ $# -lt 2 || $# -gt 6 ]]
then
	printf "Usage: $(basename $0) <proteins_fasta> <supfam_hmm_db> " >&2
	printf "<number_of_additional_threads (default: 0)>\n" >&2
	printf "<min_domain_evalue_cutoff (default 0.01)> " >&2
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

supfam_hmm_db=$2
if [[ ! -f $supfam_hmm_db ]]
then
	echo "$supfam_hmm_db is not a file! Aborting!" >&2
	exit 1
fi

number_of_additional_threads=0
if [[ ! -z "$3" ]]
then
	number_of_additional_threads=$3
fi

min_domain_evalue_cutoff=0.01
if [[ ! -z "$4" ]]
then
	min_domain_evalue_cutoff=$4
fi

min_aln_length_ratio=0.7
if [[ ! -z "$5" ]]
then
	min_aln_length_ratio=$5
fi

max_overlap_ratio=0.1
if [[ ! -z "$6" ]]
then
	max_overlap_ratio=$6
fi


echo "$(date +%F_%T) - Calling hmmsearch against the SuperFamily db now..."
hmmsearch --notextw --domE $min_domain_evalue_cutoff \
			--cpu $number_of_additional_threads \
			--domtblout ${proteins_fasta%.*}.supfam.domtblout \
    		$supfam_hmm_db $proteins_fasta 1> /dev/null
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

echo "$(date +%F_%T) - Filtering fragments and converting output to gff now..."
supfam_gff=${proteins_fasta%_*}_supfam.gff
grep -v '^#' ${proteins_fasta%.*}.supfam.domtblout | \
awk '{print $1,$3,$4,$5,$6,$7,$8,$13,$14,$16,$17,$20,$21}' | \
sort -k1,1 -k7,7nr -k6,6n | \
$(dirname $0)/hmmsearch_fragmented_hits_filter.py -a $min_aln_length_ratio \
												  -o $max_overlap_ratio \
											"$tool_and_version" > $supfam_gff
exit_code=$?
if [[ $exit_code -ne 0 ]]
then
	echo "$(date +%F_%T) - One of the pipe commands failed! Aborting!" >&2
	exit $exit_code
fi

echo "$(date +%F_%T) - Done with SuperFamily annotation."
