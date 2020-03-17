#!/usr/bin/env bash

if [[ $# -lt 2 || $# -gt 3 ]]
then
	printf "Usage: $(basename $0) <input_fasta> <project_type> " >&2
	printf "<threads (default: 1)>\n" >&2
	exit 1
fi

input_fasta=$1
if [[ ! -f $input_fasta ]]
then
	echo "$input_fasta is not a file! Aborting!" >&2
	exit 1
fi

project_type=$2
if [[ "$project_type" != "isolate" && "$project_type" != "metagenome" ]]
then
	echo "Unknown project type: $project_type" >&2
	echo "The project type has to be either 'isolate' or 'metagenome'." >&2
	exit 1
fi

threads=1
if [[ ! -z "$3" ]]
then
	threads=$3
fi


if [[ $threads -eq 1 ]]
then
	echo "$(date +%F_%T) - Running tRNAscan-SE in bacterial mode now..."
	bacterial_run_out="${input_fasta%_*}_trnascan_bacterial.out"
	tRNAscan-SE -B --thread 0 $input_fasta &> $bacterial_run_out
	exit_code=$?
	if [[ $exit_code -ne 0 ]]
	then
		echo "tRNAscan-SE bacterial run failed! Aborting!" >&2
		exit $exit_code
	fi

	echo "$(date +%F_%T) - Running tRNAscan-SE in archaeal mode now..."
	archaeal_run_out="${input_fasta%_*}_trnascan_archaeal.out"
	tRNAscan-SE -A --thread 0 $input_fasta &> $archaeal_run_out
	exit_code=$?
	if [[ $exit_code -ne 0 ]]
	then
		echo "tRNAscan-SE archaeal run failed! Aborting!" >&2
		exit $exit_code
	fi
else
	printf "$(date +%F_%T) - Splitting up input fasta into about $threads "
	printf "pieces now and then run tRNAscan-SE on them separately. "
	printf "First with the bacterial and then with the archaeal model...\n"
	tmp_dir=$(dirname $input_fasta)
	filesize=$(ls -l $input_fasta | awk '{print $5}')
    blocksize=$((($filesize / $threads) + 200000 )) # +200k to lower number of splits, since cori often gets
                                                    # hickups when almost all cores are used for GNU parallel

	# Use parallel to split up the input and
	# run tRNAscan-SE in parallel on those splits
	cat $input_fasta | parallel --pipe --recstart '>' \
							--blocksize $blocksize \
							'cat > '$tmp_dir'/tmp.$$.split.fna; \
							tRNAscan-SE -B --thread 0 \
							'$tmp_dir'/tmp.$$.split.fna &> \
							'$tmp_dir'/tmp.bacterial.$$.out; \
                            tRNAscan-SE -A --thread 0 \
                            '$tmp_dir'/tmp.$$.split.fna &> \
                            '$tmp_dir'/tmp.archaeal.$$.out;'
	exit_code=$?
	if [[ $exit_code -ne 0 ]]
	then
		echo "GNU parallel run failed! Aborting!" >&2
		exit $exit_code
	fi

	echo "$(date +%F_%T) - Concatenating bacterial split result files now..."
	bacterial_run_out="${input_fasta%_*}_trnascan_bacterial.out"
	cat $tmp_dir/tmp.bacterial.* > $bacterial_run_out
	exit_code=$?
	if [[ $exit_code -ne 0 ]]
	then
		echo "Concatenating tmp bacterial output failed! Aborting!" >&2
		exit $exit_code
	fi

	echo "$(date +%F_%T) - Concatenating archaeal split result files now..."
	archaeal_run_out="${input_fasta%_*}_trnascan_archaeal.out"
	cat $tmp_dir/tmp.archaeal.* > $archaeal_run_out
	exit_code=$?
	if [[ $exit_code -ne 0 ]]
	then
		echo "Concatenating tmp archaeal output failed! Aborting!" >&2
		exit $exit_code
	fi

	echo "$(date +%F_%T) - Deleting tmp files now..."
	rm $tmp_dir/tmp.*
fi

out_gff="${input_fasta%_*}_trna.gff"
echo "$(date +%F_%T) - Picking best results now and creating gff..."
$(dirname $0)/trna_pick_and_transform_to_gff.py \
			${input_fasta%_*}_trnascan_* > $out_gff
exit_code=$?
if [[ $exit_code -ne 0 ]]
then
	echo "The script trna_pick_and_transform_to_gff.py failed! Aborting!" >&2
	exit $exit_code
fi

echo "$(date +%F_%T) - Done with tRNA prediction."
