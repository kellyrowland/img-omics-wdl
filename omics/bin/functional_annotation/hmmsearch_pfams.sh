#!/usr/bin/env bash

if [[ $# -lt 3 || $# -gt 6 ]]
then
	printf "Usage: $(basename $0) <proteins_fasta> " >&2
	printf "<pfam_hmm_db> <pfam_claninfo_tsv> " >&2
	printf "<number_of_additional_threads (default: 0)> " >&2
	printf "<number_of_parallel_hmmsearch_instances (default: 0)> " >&2
	printf "<approximate_number_of_total_proteins (default: 0)>\n" >&2
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

number_of_parallel_hmmsearch_instances=0
if [[ ! -z "$5" ]]
then
	number_of_parallel_hmmsearch_instances=$5
fi

approximate_number_of_total_proteins=0
if [[ ! -z "$6" ]]
then
	approximate_number_of_total_proteins=$6
fi


if [[ $number_of_additional_threads -gt $number_of_parallel_hmmsearch_instances ]]
then
    hmmsearch_threads=$(echo $number_of_additional_threads / $number_of_parallel_hmmsearch_instances | bc)
    printf "$(date +%F_%T) - Splitting up proteins fasta into $number_of_parallel_hmmsearch_instances "
    printf "pieces now and then run hmmsearch on them separately with $hmmsearch_threads "
    printf "threads each against the Pfam db...\n"
    tmp_dir=$(dirname $proteins_fasta)
    filesize=$(ls -l $proteins_fasta | awk '{print $5}')
    blocksize=$((($filesize / $number_of_parallel_hmmsearch_instances) + 20000))

    hmmsearch_base_cmd="hmmsearch --notextw --cut_tc"
    if [[ $approximate_number_of_total_proteins -gt 0 ]]
    then
        hmmsearch_base_cmd="$hmmsearch_base_cmd -Z $approximate_number_of_total_proteins"
    fi
    hmmsearch_base_cmd="$hmmsearch_base_cmd --cpu $hmmsearch_threads "
    # Use parallel to split up the input and
    # run hmmsearch in parallel on those splits
    cat $proteins_fasta | parallel --pipe --recstart '>' \
                            --blocksize $blocksize \
                            'cat > '$tmp_dir'/tmp.$$.split.faa; ' \
                            $hmmsearch_base_cmd \
			                '--domtblout '$tmp_dir'/tmp.pfam.$$.domtblout' \
    		                $pfam_hmm_db $tmp_dir'/tmp.$$.split.faa 1> /dev/null;'
    exit_code=$?
    if [[ $exit_code -ne 0 ]]
    then
        echo "GNU parallel run failed! Aborting!" >&2
        exit $exit_code
    fi

    echo "$(date +%F_%T) - Concatenating split result files now..."
    cat $tmp_dir/tmp.pfam.* > ${proteins_fasta%.*}.pfam.domtblout
    exit_code=$?
    if [[ $exit_code -ne 0 ]]
    then
        echo "Concatenating split outputs failed! Aborting!" >&2
        exit $exit_code
    fi

    echo "$(date +%F_%T) - Deleting tmp files now..."
    rm $tmp_dir/tmp.*
else
    echo "$(date +%F_%T) - Calling hmmsearch to predict Pfams now..."
    hmmsearch_cmd="hmmsearch --notextw --cut_tc"
    if [[ $approximate_number_of_total_proteins -gt 0 ]]
    then
        hmmsearch_cmd="$hmmsearch_cmd -Z $approximate_number_of_total_proteins"
    fi
    hmmsearch_cmd="$hmmsearch_cmd --domtblout ${proteins_fasta%.*}.pfam.domtblout "
    hmmsearch_cmd="$hmmsearch_cmd $pfam_hmm_db $proteins_fasta 1> /dev/null"
    $hmmsearch_cmd
    exit_code=$?
    if [[ $exit_code -ne 0 ]]
    then
        echo "$(date +%F_%T) - hmmsearch failed! Aborting!" >&2
        exit $exit_code
    fi
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
