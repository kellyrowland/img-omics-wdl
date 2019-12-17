#!/usr/bin/env bash

if [[ $# -lt 1 || $# -gt 2 ]]
then
	printf "Usage: $(basename $0) <proteins_fasta> <organism_type>\n" >&2
	exit 1
fi


proteins_fasta=$1
if [[ ! -f $proteins_fasta ]]
then
	echo "$proteins_fasta is not a file! Aborting!" >&2
	exit 1
fi

organism_type=$2
if [[ -z $organism_type ]]
then
	organism_type="gram-"
fi
if [[ ! -z $organism_type && \
		$organism_type != "gram+" && \
		$organism_type != "gram-" && \
		$organism_type != "euk" ]]
then
	echo "Unknown organism type: $organism_type! Aborting!" >&2
	exit 1
fi

gff_out=${proteins_fasta%_*}_cleavage_sites.gff
echo "$(date +%F_%T) - Getting SignalP version..."
signalp_version=$(signalp -V)
echo "$(date +%F_%T) - Calling SignalP now and creating gff output..."
signalp -t $organism_type -f short $proteins_fasta | \
	grep -v '^#' | \
	awk -v sv="$signalp_version" -v ot="$organism_type" \
		'$10 == "Y" {print $1"\t"sv"\tcleavage_site\t"$3-1"\t"$3"\t"$2\
		"\t.\t.\tD-score="$9";network="$12";organism_type="ot}' > $gff_out
exit_code=$?
if [[ $exit_code -ne 0 ]]
then
	echo "$(date +%F_%T) - Signalp or the parsing pipe failed! Aborting!" >&2
	exit $exit_code
fi

echo "$(date +%F_%T) - Done with cleavage site annotation."
