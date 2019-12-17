#!/usr/bin/env bash

if [[ $# -lt 3 || $# -gt 4 ]]
then
	printf "Usage: $(basename $0) <proteins_fasta> " >&2
	printf "<pfam_hmm_db> <pfam_claninfo_tsv> " >&2
	printf "<number_of_additional_threads (default: 0)>\n" >&2
	exit 1
fi


proteins_fasta=$1
if [[ ! -f $proteins_fasta ]]
then
	echo "$proteins_fasta is not a file! Aborting!" >&2
	exit 1
fi

pfam_hmm_db=$2
if [[ ! -f $pfam_hmm_db ]]
then
	echo "$pfam_hmm_db is not a file! Aborting!" >&2
	exit 1
fi

pfam_claninfo_tsv=$3
if [[ ! -f $pfam_claninfo_tsv ]]
then
	echo "$pfam_claninfo_tsv is not a file! Aborting!" >&2
	exit 1
fi

number_of_additional_threads=0
if [[ ! -z "$4" ]]
then
	number_of_additional_threads=$4
fi


echo "$(date +%F_%T) - Calling hmmsearch to predict Pfams now..."
hmmsearch --notextw --cut_tc \
			--cpu $number_of_additional_threads \
			--domtblout ${proteins_fasta%.*}.pfam.domtblout \
    		$pfam_hmm_db $proteins_fasta 1> /dev/null
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

echo "$(date +%F_%T) - Executing clan filtering and transforming to GFF now..."
pfams_gff=${proteins_fasta%_*}_pfam.gff
grep -v '^#' ${proteins_fasta%.*}.pfam.domtblout | \
awk '{print $1,$4,$5,$6,13,$14,$16,$17,$20,$21}' | \
sort -k1,1 -k6,6nr -k5,5n | \
$(dirname $0)/pfam_clan_filter.py "$tool_and_version" \
	$pfam_claninfo_tsv > $pfams_gff
exit_code=$?
if [[ $exit_code -ne 0 ]]
then
	echo "$(date +%F_%T) - The pipe command failed! Aborting!" >&2
	exit $exit_code
fi

echo "$(date +%F_%T) - Done with Pfam annotation."
