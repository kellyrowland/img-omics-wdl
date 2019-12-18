#!/usr/bin/env bash

if [[ $# != 2 ]]
then
	printf "Usage: $(basename $0) <proteins_fasta> <model_file>\n" >&2
	exit 1
fi


proteins_fasta=$1
if [[ ! -f $proteins_fasta ]]
then
	echo "$proteins_fasta is not a file! Aborting!" >&2
	exit 1
fi

model_file=$2
if [[ ! -f $model_file ]]
then
	echo "$model_file is not a file! Aborting!" >&2
	exit 1
fi


echo "$(date +%F_%T) - Getting decodeanhmm version now..."
tool_and_version=$(decodeanhmm -v 2>&1 | head -1)

echo "$(date +%F_%T) - Calling decodeanhmm now to predict TMHs and \
transforming output to gff..."
gff_out="${proteins_fasta%_*}_tmh.gff"
background="0.081 0.015 0.054 0.061 0.040 0.068 0.022 0.057 0.056 0.093 0.025"
background="$background 0.045 0.049 0.039 0.057 0.068 0.058 0.067 0.013 0.032"
sed 's/\*/X/g' $proteins_fasta | \
	decodeanhmm -N 1 -background $background -PrintNumbers \
	$model_file 2> /dev/null | \
		$(dirname $0)/decodeanhmm_parser.py "$tool_and_version" > $gff_out
exit_code=$?
if [[ $exit_code -ne 0 ]]
then
	echo "$(date +%F_%T) - Something in the pipe failed! Aborting!" >&2
	exit $exit_code
fi

echo "$(date +%F_%T) - Done with TMH annotation."
