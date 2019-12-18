#!/usr/bin/env bash

if [[ $# -lt 4 || $# -gt 5 ]]
then
	printf "Usage: $(basename $0) <input_fasta> " >&2
	printf "<rfam_cm> <rfam_clan_info_tsv> <rfam_img_ncbi_feature_tsv>" >&2
	printf "<number_of_additional_threads (default: 0)>\n" >&2
	exit 1
fi

input_fasta=$1
if [[ ! -f $input_fasta ]]
then
	echo "$input_fasta is not a file! Aborting!" >&2
	exit 1
fi

rfam_cm=$2
if [[ ! -f $rfam_cm ]]
then
	echo "$rfam_cm is not a file! Aborting!" >&2
	exit 1
fi

rfam_clan_info_tsv=$3
if [[ ! -f $rfam_clan_info_tsv ]]
then
	echo "$rfam_clan_info_tsv is not a file! Aborting!" >&2
	exit 1
fi

rfam_img_ncbi_feature_tsv=$4
if [[ ! -f $rfam_img_ncbi_feature_tsv ]]
then
	echo "$rfam_img_ncbi_feature_tsv is not a file! Aborting!" >&2
	exit 1
fi

number_of_additional_threads=0
if [[ ! -z "$5" ]]
then
	number_of_additional_threads=$5
fi


echo "$(date +%F_%T) - Running cmsearch now..."
rfam_tbl_out="${input_fasta%_*}_rfam.tbl"
cmsearch --notextw --cut_tc --cpu $number_of_additional_threads \
  --tblout $rfam_tbl_out $rfam_cm $input_fasta 1> /dev/null
exit_code=$?
if [[ $exit_code -ne 0 ]]
then
	echo "cmsearch failed! Aborting!" >&2
	exit $exit_code
fi

echo "$(date +%F_%T) - Getting cmsearch version now..."
tool_and_version=$(cmsearch -h | grep INFERNAL)
tool_and_version=${tool_and_version:2}

echo "$(date +%F_%T) - Executing clan filtering and GFF transformation now..."
rfam_gff="${input_fasta%_*}_rfam.gff"
grep -v '^#' $rfam_tbl_out |\
awk '$17 == "!" {print $1,$3,$4,$6,$7,$8,$9,$10,$11,$15,$16}' |\
sort -k1,1 -k10,10nr -k11,11n |\
$(dirname $0)/rfam_clan_filter.py "$tool_and_version" \
	$rfam_clan_info_tsv $rfam_img_ncbi_feature_tsv > $rfam_gff
exit_code=$?
if [[ $exit_code -ne 0 ]]
then
	echo "The pipe command failed! Aborting!" >&2
	exit $exit_code
fi


printf "$(date +%F_%T) - Separating feature types now (to "
printf "accommodate the overlap resolution during the gff files merging)...\n"
misc_and_regulatory_gff="${input_fasta%_*}_rfam_misc_bind_misc_feature_regulatory.gff"
awk -F'\t' '$3 == "misc_bind" || $3 == "misc_feature" || \
			$3 == "regulatory" {print $0}' $rfam_gff > $misc_and_regulatory_gff
exit_code=$?
if [[ $exit_code -ne 0 ]]
then
	echo "Pulling out misc and regulatory failed! Aborting!" >&2
	exit $exit_code
fi

rrna_gff="${input_fasta%_*}_rfam_rrna.gff"
awk -F'\t' '$3 == "rRNA" {print $0}' $rfam_gff > $rrna_gff
exit_code=$?
if [[ $exit_code -ne 0 ]]
then
	echo "Pulling out rRNAs failed! Aborting!" >&2
	exit $exit_code
fi

ncrna_tmrna_gff="${input_fasta%_*}_rfam_ncrna_tmrna.gff"
awk -F'\t' '$3 == "ncRNA" || $3 == "tmRNA" {print $0}' \
		$rfam_gff > $ncrna_tmrna_gff
exit_code=$?
if [[ $exit_code -ne 0 ]]
then
	echo "Pulling out ncRNAs and tmRNAs failed! Aborting!" >&2
	exit $exit_code
fi


echo "$(date +%F_%T) - Done with Rfam annotation."
