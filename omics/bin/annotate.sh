#!/usr/bin/env bash

if [[ $# -ne 1 ]]
then
	echo "Usage: $(basename $0) <annotation_config.yaml>\n" >&2
	exit 1
fi

imgap_bin_dir=$(dirname $0)

config_yaml=$1
source $imgap_bin_dir/set_env.sh $config_yaml

run_folder=$(dirname $imgap_input_fasta)

# Create associative array to store imgap run information
declare -A imgap_run_info=(["imgap_version"]=$(cat $imgap_bin_dir/../VERSION))

if [[ "$imgap_structural_annotation_execute" == "True" ]]
then
	sa_log=${imgap_input_fasta%_*}_structural_annotation.log

    if [[ ! -f $run_folder/SA_ALL_DONE ]]
    then
        if [[ "$imgap_structural_annotation_pre_qc_execute" == "True" ]]
        then
            if [[ ! -f $run_folder/PRE_QC_DONE ]]
            then
                echo "$(date +%F_%T) - Starting pre-annotation qc now..."
                echo "Pre-QC" >> $run_folder/started_modules.log
                qc_cmd="/usr/bin/time $imgap_bin_dir/qc/pre-annotation/simple_fasta_qc.sh "
                qc_cmd="$qc_cmd $imgap_project_type $imgap_input_fasta "
                qc_cmd="$qc_cmd $imgap_project_id yes"
                qced_fasta=`$qc_cmd`
                exit_code=$?
                if [[ $exit_code -ne 0 ]]
                then
                    echo "$(date +%F_%T) - The pre-annotation qc failed! Aborting!" >&2
                    echo "Pre-Annotation QC failed!" > $run_folder/ERROR
                    exit $exit_code
                fi
                export imgap_input_fasta=$qced_fasta
                touch $run_folder/PRE_QC_DONE
            else
                echo "$(date +%F_%T) - Pre-annotation qc appears to have been run earlier already. Skipping!"
            fi
        fi

        if [[ ! -f $run_folder/SA_DONE ]]
        then
            echo "$(date +%F_%T) - Starting structural annotation now..."
            /usr/bin/time $imgap_bin_dir/structural_annotation/structural_annotation.sh \
                            &> $sa_log
            exit_code=$?
            if [[ $exit_code -ne 0 ]]
            then
                echo "$(date +%F_%T) - The structural annotation failed! Aborting!" >&2
                echo "Structural annotation failed!" > $run_folder/ERROR
                exit $exit_code
            fi
            touch $run_folder/SA_DONE
        else
            printf "$(date +%F_%T) - "
            printf "The Structural annotation script appears to have been run earlier already. "
            printf "Skipping!\n"
        fi


        if [[ "$imgap_project_type" == "isolate" ]]
        then
            if [[ ! -f $run_folder/POST_QC_DONE ]]
            then
                echo "$(date +%F_%T) - Starting post-structural-annotation qc now..."
                echo "Post-QC" >> $run_folder/started_modules.log
                qc_cmd="/usr/bin/time $imgap_bin_dir/qc/post-annotation"
                qc_cmd="${qc_cmd}/genome_structural_annotation_sanity.py "
                qc_cmd="$qc_cmd $imgap_input_fasta "
                qc_cmd="$qc_cmd ${imgap_input_fasta%_*}_structural_annotation.gff"
                $qc_cmd
                exit_code=$?
                if [[ $exit_code -ne 0 ]]
                then
                    printf "$(date +%F_%T) - The post-structural-annotation " >&2
                    printf "qc script failed! Aborting!\n" >&2
                    echo "Post-Annotation QC failed!" > $run_folder/ERROR
                    exit $exit_code
                fi

                if [[ -f "$(dirname $imgap_input_fasta)/QC_NEEDED.txt" ]]
                then
                    printf "$(date +%F_%T) - A QC_NEEDED.txt file is present.\n" >&2
            #		printf "Stopping all further processing now!\n" >&2
                    printf "Please check $(dirname $imgap_input_fasta)/QC_NEEDED.txt for " >&2
                    printf "details!\n" >&2
            #		exit 1
                fi
                touch $run_folder/POST_QC_DONE
            else
                printf "$(date +%F_%T) - "
                printf "Post QC of teh structural annotation  appears to have been run earlier already. "
                printf "Skipping!\n"
            fi
        fi
        touch $run_folder/SA_ALL_DONE
    else
        printf "$(date +%F_%T) - "
        printf "All steps of the structural annotation appear to have been run earlier already. "
        printf "Skipping!\n"
    fi

    imgap_run_info["sa_bins"]="$(grep -A 1 '^Used Programs' $sa_log | tail -1)"
	imgap_run_info["sa_dbs"]="$(grep -A 1 '^Used DBs' $sa_log | tail -1)"
fi


if [[ "$imgap_functional_annotation_execute" == "True" ]]
then
    fa_log=${imgap_input_fasta%_*}_functional_annotation.log

	if [[ ! -f $run_folder/FA_DONE ]]
	then
        echo "$(date +%F_%T) - Starting functional annotation now..."
        if [[ "$imgap_structural_annotation_execute" == "True" ]]
        then
            export imgap_input_fasta=${imgap_input_fasta%_*}_proteins.faa
        fi
        /usr/bin/time $imgap_bin_dir/functional_annotation/functional_annotation.sh \
                        &> $fa_log
        exit_code=$?
        if [[ $exit_code -ne 0 ]]
        then
            echo "$(date +%F_%T) - The functional annotation failed! Aborting!" >&2
            echo "Functional annotation failed!" > $run_folder/ERROR
            exit $exit_code
        fi
        touch $run_folder/FA_DONE
    else
        echo "$(date +%F_%T) - The FA appears to have been run earlier already. Skipping!"
    fi

	imgap_run_info["fa_bins"]="$(grep -A 1 '^Used Programs' $fa_log | tail -1)"
	imgap_run_info["fa_dbs"]="$(grep -A 1 '^Used DBs' $fa_log | tail -1)"
fi

echo "ALL DONE" >> $run_folder/started_modules.log

imgap_run_info_file=$run_folder/${imgap_project_id}_imgap.info
echo "IMGAP Version: ${imgap_run_info["imgap_version"]}" > $imgap_run_info_file
if [[ -v imgap_run_info["sa_bins"] ]]
then
	printf "Structural Annotation Programs Used: " >> $imgap_run_info_file
	printf "${imgap_run_info["sa_bins"]}\n" >> $imgap_run_info_file
fi

if [[ -v imgap_run_info["sa_dbs"] ]]
then
	printf "Structural Annotation DBs Used: " >> $imgap_run_info_file
	printf "${imgap_run_info["sa_dbs"]}\n" >> $imgap_run_info_file
fi

if [[ -v imgap_run_info["fa_bins"] ]]
then
	printf "Functional Annotation Programs Used: " >> $imgap_run_info_file
	printf "${imgap_run_info["fa_bins"]}\n" >> $imgap_run_info_file
fi

if [[ -v imgap_run_info["fa_dbs"] ]]
then
	printf "Functional Annotation DBs Used: " >> $imgap_run_info_file
	printf "${imgap_run_info["fa_dbs"]}\n" >> $imgap_run_info_file
fi

done_file=$run_folder/DONE
printf "Structural Annotation: " >> $done_file
printf "$imgap_structural_annotation_execute\n" >> $done_file
printf "Functional Annotation: " >> $done_file
printf "$imgap_functional_annotation_execute\n" >> $done_file

echo "$(date +%F_%T) - Done with all annotation tasks."
echo

