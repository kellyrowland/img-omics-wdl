#!/usr/bin/env bash

if [[ $# -gt 1 ]]
then
	printf "Usage: $(basename $0) <annotation_config.yaml>\n" >&2
	exit 1
fi

sa_bin_dir=$(dirname $0)

config_yaml=$1
if [[ ! -z "$config_yaml" ]]
then
	# If an argument was given, it is assumed that it is the config yaml file.
	# It's used as argument to ../set_env.sh, which gets sourced to export
	# the environment variables.
	if [[ ! -f $config_yaml ]]
	then
		echo "$config_yaml is not a file! Aborting!" >&2
		exit 1
	fi
	source $sa_bin_dir/../set_env.sh $config_yaml
fi

# If no argument was given to the script, it is assumed that the environment
# variables are already set.


run_folder=$(dirname $imgap_input_fasta)
cd $run_folder


if [[ "$imgap_structural_annotation_trnascan_se_execute" == "True" ]]
then
    if [[ ! -f $run_folder/TRNA_DONE ]]
    then
        # Reducing the number of additional threads by 1/3,
        # since cori seems to have a lot of GNU parallel thread hickups
        # when the thread nomber is close to the number oc cores.
        trna_threads=$(echo "$imgap_additional_threads * 0.66" | bc)
        trna_threads=${trna_threads%.*}
        echo "$(date +%F_%T) - Predicting tRNAs now..."
        echo "TRNA" >> $run_folder/started_modules.log
        trna_log=${imgap_input_fasta%_*}_trna.log
        /usr/bin/time $sa_bin_dir/trnascan-se_trnas.sh \
                        $imgap_input_fasta $imgap_project_type \
                        $trna_threads &> $trna_log
        exit_code=$?
        if [[ $exit_code -ne 0 ]]
        then
            echo "predict_trnas.sh failed! Aborting!" >&2
            exit $exit_code
        fi
        touch $run_folder/TRNA_DONE
    else
        echo "$(date +%F_%T) - The tRNA prediction appears to have been run earlier already. Skipping!"
    fi
fi


if [[ "$imgap_structural_annotation_rfam_execute" == "True" ]]
then
    if [[ ! -f $run_folder/RFAM_DONE ]]
    then
        echo "$(date +%F_%T) - Predicting RNAs now..."
        echo "RFAM" >> $run_folder/started_modules.log
        rfam_log=${imgap_input_fasta%_*}_rfam.log
        /usr/bin/time $sa_bin_dir/cmsearch_rfams.sh $imgap_input_fasta \
                        $imgap_structural_annotation_rfam_cm \
                        $imgap_structural_annotation_rfam_claninfo_tsv \
                        $imgap_structural_annotation_rfam_feature_lookup_tsv \
                        $imgap_additional_threads &> $rfam_log
        exit_code=$?
        if [[ $exit_code -ne 0 ]]
        then
            echo "cmsearch_rfams.sh failed! Aborting!" >&2
            exit $exit_code
        fi
        touch $run_folder/RFAM_DONE
    else
        echo "$(date +%F_%T) - The RNA prediction appears to have been run earlier already. Skipping!"
    fi
fi


crt_log=${imgap_input_fasta%_*}_crt.log
prodigal_log=${imgap_input_fasta%_*}_prodigal.log
genemark_log=${imgap_input_fasta%_*}_genemark.log
# If we've can use at least three threads,
# we can run the last three tools in parallel
if [[ $imgap_additional_threads -ge 2 ]]
then
	PIDS=()
	TOOL_NAMES=()
	index_tracker=0

	if [[ "$imgap_structural_annotation_crt_execute" == "True" ]]
	then
        if [[ ! -f $run_folder/CRT_DONE ]]
        then
            index_tracker=$((index_tracker + 1))
            echo "$(date +%F_%T) - Starting CRISPRs prediction now..."
            echo "CRT" >> $run_folder/started_modules.log
            (/usr/bin/time $sa_bin_dir/crt_crisprs.sh \
                            $imgap_input_fasta &> $crt_log) &
            PIDS[$index_tracker]=$!
            TOOL_NAMES[$index_tracker]="CRT"
        else
            echo "$(date +%F_%T) - The CRISPR prediction appears to have been run earlier already. Skipping!"
        fi
	fi

	if [[ "$imgap_structural_annotation_prodigal_execute" == "True" ]]
	then
        if [[ ! -f $run_folder/Prodigal_DONE ]]
        then
            index_tracker=$((index_tracker + 1))
            echo "$(date +%F_%T) - Starting CDSs prediction via Prodigal now..."
            echo "Prodigal" >> $run_folder/started_modules.log
            (/usr/bin/time $sa_bin_dir/prodigal_genes.sh $imgap_input_fasta \
                                         $imgap_project_type &> $prodigal_log) &
            PIDS[$index_tracker]=$!
            TOOL_NAMES[$index_tracker]="Prodigal"
        else
            echo "$(date +%F_%T) - Prodigal appears to have been run earlier already. Skipping!"
        fi
	fi

	if [[ "$imgap_structural_annotation_genemark_execute" == "True" ]]
	then
        if [[ ! -f $run_folder/GeneMark_DONE ]]
        then
            index_tracker=$((index_tracker + 1))
            echo "$(date +%F_%T) - Starting CDSs predicition via GeneMark now..."
            echo "GeneMark" >> $run_folder/started_modules.log
            (/usr/bin/time $sa_bin_dir/genemark_genes.sh $imgap_input_fasta \
                                            $imgap_project_type &> $genemark_log) &
            PIDS[$index_tracker]=$!
            TOOL_NAMES[$index_tracker]="GeneMark"
        else
            echo "$(date +%F_%T) - GeneMark appears to have been run earlier already. Skipping!"
        fi
	fi

	if [[ $index_tracker -gt 0 ]]
	then
		for i in $(seq 1 $index_tracker)
		do
			wait ${PIDS[$i]}
			exit_code=$?
			if [[ $exit_code -ne 0 ]]
			then
				echo "${TOOL_NAMES[$i]} failed." >&2
				exit $exit_code
			fi
            touch $run_folder/${TOOL_NAMES[$i]}_DONE
		done
	fi
else
	if [[ "$imgap_structural_annotation_crt_execute" == "True" ]]
	then
        if [[ ! -f $run_folder/CRT_DONE ]]
        then
            echo "$(date +%F_%T) - Predicting CRISPRs now..."
            echo "CRT" >> $run_folder/started_modules.log
            /usr/bin/time $sa_bin_dir/crt_crisprs.sh \
                            $imgap_input_fasta &> $crt_log &
            exit_code=$?
            if [[ $exit_code -ne 0 ]]
            then
                echo "crt_crisprs.sh failed! Aborting!" >&2
                exit $exit_code
            fi
            touch $run_folder/CRT_DONE
        else
            echo "$(date +%F_%T) - The CRISPR prediction appears to have been run earlier already. Skipping!"
        fi
	fi

	if [[ "$imgap_structural_annotation_prodigal_execute" == "True" ]]
    then
        if [[ ! -f $run_folder/Prodigal_DONE ]]
        then
            echo "$(date +%F_%T) - Predicting CDSs via Prodigal now..."
            echo "Prodigal" >> $run_folder/started_modules.log
            /usr/bin/time $sa_bin_dir/prodigal_genes.sh \
                            $imgap_input_fasta $imgap_project_type &> $prodigal_log
            exit_code=$?
            if [[ $exit_code -ne 0 ]]
            then
                echo "prodigal_genes.sh failed! Aborting!" >&2
                exit $exit_code
            fi
            touch $run_folder/Prodigal_DONE
        else
            echo "$(date +%F_%T) - Prodigal appears to have been run earlier already. Skipping!"
        fi
	fi

	if [[ "$imgap_structural_annotation_genemark_execute" == "True" ]]
	then
        if [[ ! -f $run_folder/GeneMark_DONE ]]
        then
            echo "$(date +%F_%T) - Predicting CDSs via GeneMark now..."
            echo "GeneMark" >> $run_folder/started_modules.log
            /usr/bin/time $sa_bin_dir/genemark_genes.sh \
                            $imgap_input_fasta $imgap_project_type &> $genemark_log
            exit_code=$?
            if [[ $exit_code -ne 0 ]]
            then
                echo "genemark_genes.sh failed! Aborting!" >&2
                exit $exit_code
            fi
            touch $run_folder/GeneMark_DONE
        else
            echo "$(date +%F_%T) - GeneMark appears to have been run earlier already. Skipping!"
        fi
	fi
fi


final_gff=${imgap_input_fasta%_*}_structural_annotation.gff
if [[ ! -f $run_folder/GFF_MERGE_DONE ]]
then
    echo "$(date +%F_%T) - Merging gff files now..."
    echo "GFF Merging" >> $run_folder/started_modules.log
    merger_args="-f $imgap_input_fasta"
    if [[ "$imgap_structural_annotation_rfam_execute" == "True" ]]
    then
        misc_and_regulatory_gff=${imgap_input_fasta%_*}_rfam_misc_bind_misc_feature_regulatory.gff
        rrna_gff=${imgap_input_fasta%_*}_rfam_rrna.gff
        merger_args="$merger_args -a $misc_and_regulatory_gff $rrna_gff"

#        cmsearch_version=$(cmsearch -h | grep INFERNAL | cut -d' ' -f2,3)
#        used_programs+=("$cmsearch_version")
#        rfam_version=$(basename $(realpath $(dirname $imgap_structural_annotation_rfam_cm)))
#        used_dbs+=("Rfam $rfam_version")
    fi

    if [[ "$imgap_structural_annotation_trnascan_se_execute" == "True" ]]
    then
        trna_gff=${imgap_input_fasta%_*}_trna.gff
        merger_args="$merger_args $trna_gff"

        trnascan_version=$(tRNAscan-SE 2>&1 | head -2 | tail -1)
#        used_programs+=("$trnascan_version")
    fi

    if [[ "$imgap_structural_annotation_rfam_execute" == "True" ]]
    then
        ncrna_tmrna_gff=${imgap_input_fasta%_*}_rfam_ncrna_tmrna.gff
        merger_args="$merger_args $ncrna_tmrna_gff"
    fi

    if [[ "$imgap_structural_annotation_crt_execute" == "True" ]]
    then
        crt_gff=${imgap_input_fasta%_*}_crt.gff
        merger_args="$merger_args $crt_gff"

        crt_version=$(java -jar $(which CRT-CLI.jar) -version | cut -d' ' -f1,6)
#        used_programs+=("$crt_version")
    fi

    if [[ "$imgap_structural_annotation_genemark_execute" == "True" ]]
    then
        genemark_gff=${imgap_input_fasta%_*}_genemark.gff
        merger_args="$merger_args $genemark_gff"

        genemark_version=$(gmhmmp2 2>&1 | grep version)
#        used_programs+=("$(echo ${genemark_version/ersion /})")
    fi

    if [[ "$imgap_structural_annotation_prodigal_execute" == "True" ]]
    then
        prodigal_gff=${imgap_input_fasta%_*}_prodigal.gff
        merger_args="$merger_args $prodigal_gff"

        prodigal_version=$(prodigal -v 2>&1 | head -2 | tail -1 | cut -d':' -f1)
#        used_programs+=("$prodigal_version")
    fi


    $sa_bin_dir/gff_files_merger.py $merger_args 1> $final_gff
    exit_code=$?
    if [[ $exit_code -ne 0 ]]
    then
        echo "Merging gff files failed! Aborting!" >&2
        exit $exit_code
    fi

    touch $run_folder/GFF_MERGE_DONE
else
    echo "$(date +%F_%T) - GFF files merging appears to have been run earlier already. Skipping!"
fi

if [[ "$imgap_structural_annotation_prodigal_execute" == "True" ||
		"$imgap_structural_annotation_genemark_execute" == "True" ]]
then
    if [[ ! -f $run_folder/FASTA_MERGE_DONE ]]
    then
        final_genes=${imgap_input_fasta%_*}_genes.fna
        final_proteins=${imgap_input_fasta%_*}_proteins.faa
        ff_genes_args="$final_gff"
        ff_proteins_args="$final_gff"

        echo "Fasta Merging" >> $run_folder/started_modules.log

        echo "$(date +%F_%T) - Creating final genes fasta file(s) now..."
        # We can't ran the next two commands in parallel, since
        # the fna (genes) run also rewrites the gff file.
        if [[ "$imgap_structural_annotation_genemark_execute" == "True" ]]
        then
            genemark_genes=${imgap_input_fasta%_*}_genemark_genes.fna
            genemark_proteins=${imgap_input_fasta%_*}_genemark_proteins.faa
            ff_genes_args="$ff_genes_args $genemark_genes"
            ff_proteins_args="$ff_proteins_args $genemark_proteins"
        fi

        if [[ "$imgap_structural_annotation_prodigal_execute" == "True" ]]
        then
            prodigal_genes=${imgap_input_fasta%_*}_prodigal_genes.fna
            prodigal_proteins=${imgap_input_fasta%_*}_prodigal_proteins.faa
            ff_genes_args="$ff_genes_args $prodigal_genes"
            ff_proteins_args="$ff_proteins_args $prodigal_proteins"
        fi

        $sa_bin_dir/fasta_files_merger.py $ff_genes_args 1> $final_genes
        exit_code=$?
        if [[ $exit_code -ne 0 ]]
        then
                echo "Creating the final genes fasta file(s) failed." >&2
            exit $exit_code
        fi

        echo "$(date +%F_%T) - Creating the final proteins fasta file(s) now..."
        $sa_bin_dir/fasta_files_merger.py $ff_proteins_args 1> $final_proteins
        exit_code=$?
        if [[ $exit_code -ne 0 ]]
        then
            echo "Creating the final proteins fasta files failed." >&2
            exit $exit_code
        fi

        touch $run_folder/FASTA_MERGE_DONE
    else
        echo "$(date +%F_%T) - Fasta files merging appears to have been run earlier already. Skipping!"
    fi
fi


if [[ "$imgap_structural_annotation_gff_and_fasta_stats_execute" == "True" ]]
then
    echo "$(date +%F_%T) - Calculating gff and final fasta stats now..."
    echo "GFF and Fasta Stats" >> $run_folder/started_modules.log
    $sa_bin_dir/gff_and_final_fasta_stats.py $imgap_input_fasta $final_gff
    exit_code=$?
    if [[ $exit_code -ne 0 ]]
    then
        echo "Calculating GFF and final fasta stats failed." >&2
        exit $exit_code
    fi
fi


echo "$(date +%F_%T) - Done with structural annotation."


IFS=$'\n'
# Arrays to store tool and DB versions
used_programs=($(cut -f2 $final_gff | sort | uniq))

# Write out tools and DB version info
if [[ ${#used_programs[@]} -gt 0 ]]
then
	echo "Used Programs:"
	printf "${used_programs[0]}"
	for i in $(seq 1 $((${#used_programs[@]}-1)))
	do
		printf "; ${used_programs[$i]}"
	done
	printf "\n"
fi

if [[ "$imgap_structural_annotation_rfam_execute" == "True" ]]
then
    rfam_version=$(basename $(realpath $(dirname $imgap_structural_annotation_rfam_cm)))
	printf "Used DBs:\nRfam $rfam_version\n"
fi

