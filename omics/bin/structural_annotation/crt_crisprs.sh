#!/usr/bin/env bash

if [[ $# -ne 1 ]]
then
	echo "Usage: $(basename $0) <input_fasta>" >&2
	exit 1
fi

input_fasta=$1
if [[ ! -f $input_fasta ]]
then
	echo "$input_fasta is not a file! Aborting!" >&2
	exit 1
fi

crt_out=${input_fasta%_*}_crt.out

echo "$(date +%F_%T) - Running CRT now..."
full_path_to_crt=$(which CRT-CLI.jar)
java -Xmx1536m -jar $full_path_to_crt $input_fasta $crt_out 1> /dev/null
exit_code=$?
if [[ $exit_code -ne 0 ]]
then
	echo "CRT failed! Aborting!" >&2
	exit $exit_code
fi

echo "$(date +%F_%T) - Transforming the CRT output now..."
tool_and_version=$(java -jar $full_path_to_crt -version | cut -d' ' -f1,6)
exit_code=$?
if [[ $exit_code -ne 0 ]]
then
	echo "Failed to get CRT version! Aborting!" >&2
	exit $exit_code
fi
$(dirname $0)/transform_crt_output.py $crt_out "$tool_and_version"
exit_code=$?
if [[ $exit_code -ne 0 ]]
then
	echo "Transforming CRT output failed! Aborting!" >&2
	exit $exit_code
fi

echo "$(date +%F_%T) - Done running CRT module."
