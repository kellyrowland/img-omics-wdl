#!/usr/bin/env bash

if [[ $# -gt 1 ]]
then
	printf "Usage: $(basename $0) <annotation_config.yaml>\n" >&2
	exit 1
fi

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
	source $(dirname $0)/../set_env.sh $config_yaml
fi


# Change into the folder where the input file is
run_folder=$(dirname $imgap_input_fasta)
cd $run_folder


# Create log filenames
ko_ec_log=${imgap_input_fasta%_*}_img_nr_ko_ec_gene_phylo.log
cath_funfam_log=${imgap_input_fasta%_*}_cath_funfam.log
pfam_log=${imgap_input_fasta%_*}_pfam.log
supfam_log=${imgap_input_fasta%_*}_supfam.log
cog_log=${imgap_input_fasta%_*}_cog.log
tigrfam_log=${imgap_input_fasta%_*}_tigrfam.log
smart_log=${imgap_input_fasta%_*}_smart.log
signalp_log=${imgap_input_fasta%_*}_signalp.log
tmh_log=${imgap_input_fasta%_*}_tmh.log


# Function to check if an element is not already in an array.
element_not_in_array ()
{
	local array="$1[@]"
	local seeking=$2
	local not_in=0

	for element in "${!array}"
	do
		if [[ $element == $seeking ]]
		then
			not_in=1
			break
		fi
	done

	return $not_in
}


# Arrays to store tool and db versions
used_programs=()
used_dbs=()
# Cmd args for the assign_product_names_and_create_fa_gff.py script,
# in case we are going to run it.
cmd_args=""


if [[ "$imgap_functional_annotation_ko_ec_execute" == "True" ]]
then
    if [[ ! -f $run_folder/KO_DONE ]]
    then
        printf "$(date +%F_%T) - "
        printf "Starting KO/EC and gene phylogeny prediction now...\n"
        echo "KO/EC" >> $run_folder/started_modules.log
        /usr/bin/time $(dirname $0)/lastal_img_nr_ko_ec_gene_phylo.sh \
                            $imgap_project_type $imgap_input_fasta \
                            $imgap_functional_annotation_ko_ec_img_nr_db \
                            $imgap_functional_annotation_ko_ec_md5_mapping \
                            $imgap_functional_annotation_ko_ec_taxon_to_phylo_mapping \
                            $imgap_additional_threads &> $ko_ec_log
        exit_code=$?
        if [[ $exit_code -ne 0 ]]
        then
            echo "KO/EC and gene phylogeny annotation failed! Aborting!" >&2
            exit $exit_code
        fi
        touch $run_folder/KO_DONE
    else
        printf "$(date +%F_%T) - "
        printf "KO/EC and gene phylogeny appears to have been run earlier already. Skipping!\n"
    fi

	cmd_args="-k ${imgap_input_fasta%_*}_ko_ec.gff"
	lastal_version=$(lastal -V)
	used_programs+=("$lastal_version")
	img_nr_version=$(echo $(basename $(realpath $(dirname \
						$imgap_functional_annotation_ko_ec_img_nr_db))))
	used_dbs+=("IMG-NR $img_nr_version")
fi


if [[ "$imgap_functional_annotation_smart_execute" == "True" ]]
then
    if [[ ! -f $run_folder/SMART_DONE ]]
    then
		echo "$(date +%F_%T) - Starting SMART prediction now..."
        echo "SMART" >> $run_folder/started_modules.log
		/usr/bin/time $(dirname $0)/hmmsearch_smart.sh $imgap_input_fasta \
									$imgap_functional_annotation_smart_db \
                                    $imgap_additional_threads \
									$imgap_functional_annotation_parallel_hmmsearch_instances \
                                    $imgap_functional_annotation_hmmsearch_Z_arg &> $smart_log
		exit_code=$?
		if [[ $exit_code -ne 0 ]]
		then
			echo "SMART annotation failed! Aborting!" >&2
			exit $exit_code
		fi
        touch $run_folder/SMART_DONE
    else
        printf "$(date +%F_%T) - "
        printf "SMART prediction appears to have been run earlier already. Skipping!\n"
	fi

	cmd_args="$cmd_args -s ${imgap_input_fasta%_*}_smart.gff"
	hmmsearch_version=$(hmmsearch -h | grep HMMER | cut -d' ' -f2,3)
	used_programs+=("$hmmsearch_version")
	smart_version=$(echo $(basename $(realpath $(dirname \
						$imgap_functional_annotation_smart_db))))
	used_dbs+=("SMART $smart_version")
fi


if [[ "$imgap_functional_annotation_cog_execute" == "True" ]]
then
    if [[ ! -f $run_folder/COG_DONE ]]
    then
		echo "$(date +%F_%T) - Starting COG prediction now..."
        echo "COG" >> $run_folder/started_modules.log
		/usr/bin/time $(dirname $0)/hmmsearch_cogs.sh $imgap_input_fasta \
									$imgap_functional_annotation_cog_db \
									$imgap_additional_threads \
									$imgap_functional_annotation_parallel_hmmsearch_instances \
                                    $imgap_functional_annotation_hmmsearch_Z_arg &> $cog_log
		exit_code=$?
		if [[ $exit_code -ne 0 ]]
		then
			echo "COG annotation failed! Aborting!" >&2
			exit $exit_code
		fi
        touch $run_folder/COG_DONE
    else
        printf "$(date +%F_%T) - "
        printf "COG prediction appears to have been run earlier already. Skipping!\n"
	fi

	cmd_args="$cmd_args -c ${imgap_input_fasta%_*}_cog.gff"
	hmmsearch_version=$(hmmsearch -h | grep HMMER | cut -d' ' -f2,3)
	if element_not_in_array used_programs "$hmmsearch_version"
	then
		used_programs+=("$hmmsearch_version")
	fi
	cog_version=$(echo $(basename $(realpath $(dirname \
						$imgap_functional_annotation_cog_db))))
	used_dbs+=("COG $cog_version")
fi


if [[ "$imgap_functional_annotation_tigrfam_execute" == "True" ]]
then
    if [[ ! -f $run_folder/TIGRFAM_DONE ]]
    then
		echo "$(date +%F_%T) - Starting TIGRFAM prediction now..."
        echo "TIGRFAM" >> $run_folder/started_modules.log
		/usr/bin/time $(dirname $0)/hmmsearch_tigrfams.sh $imgap_input_fasta \
									$imgap_functional_annotation_tigrfam_db \
									$imgap_additional_threads \
									$imgap_functional_annotation_parallel_hmmsearch_instances \
                                    $imgap_functional_annotation_hmmsearch_Z_arg &> $tigrfam_log
		exit_code=$?
		if [[ $exit_code -ne 0 ]]
		then
			echo "TIGRFAM annotation failed! Aborting!" >&2
			exit $exit_code
		fi
        touch $run_folder/TIGRFAM_DONE
    else
        printf "$(date +%F_%T) - "
        printf "TIGRFAM prediction appears to have been run earlier already. Skipping!\n"
	fi

	cmd_args="$cmd_args -t ${imgap_input_fasta%_*}_tigrfam.gff"
	hmmsearch_version=$(hmmsearch -h | grep HMMER | cut -d' ' -f2,3)
	if element_not_in_array used_programs "$hmmsearch_version"
	then
		used_programs+=("$hmmsearch_version")
	fi
	tigrfam_version=$(echo $(basename $(realpath $(dirname \
						$imgap_functional_annotation_tigrfam_db))))
	used_dbs+=("TIGRFAM $tigrfam_version")
fi


if [[ "$imgap_functional_annotation_superfamily_execute" == "True" ]]
then
    if [[ ! -f $run_folder/SUPERFAM_DONE ]]
    then
		echo "$(date +%F_%T) - Starting SUPERFAMILY prediction now..."
        echo "SUPERFAM" >> $run_folder/started_modules.log
		/usr/bin/time $(dirname $0)/hmmsearch_supfams.sh $imgap_input_fasta \
								$imgap_functional_annotation_superfamily_db \
								$imgap_additional_threads \
								$imgap_functional_annotation_parallel_hmmsearch_instances \
                                $imgap_functional_annotation_hmmsearch_Z_arg &> $supfam_log
		exit_code=$?
		if [[ $exit_code -ne 0 ]]
		then
			echo "SUPERFAMILY annotation failed! Aborting!" >&2
			exit $exit_code
		fi
        touch $run_folder/SUPERFAM_DONE
    else
        printf "$(date +%F_%T) - "
        printf "SuperFamily prediction appears to have been run earlier already. Skipping!\n"
	fi

	cmd_args="$cmd_args -u ${imgap_input_fasta%_*}_supfam.gff"
	hmmsearch_version=$(hmmsearch -h | grep HMMER | cut -d' ' -f2,3)
	if element_not_in_array used_programs "$hmmsearch_version"
	then
		used_programs+=("$hmmsearch_version")
	fi
	supfam_version=$(echo $(basename $(realpath $(dirname \
						$imgap_functional_annotation_superfamily_db))))
	used_dbs+=("SuperFamily $supfam_version")
fi


if [[ $imgap_functional_annotation_pfam_execute == "True" ]]
then
    if [[ ! -f $run_folder/PFAM_DONE ]]
    then
		echo "$(date +%F_%T) - Starting PFAM prediction now..."
        echo "PFAM" >> $run_folder/started_modules.log
		/usr/bin/time $(dirname $0)/hmmsearch_pfams.sh $imgap_input_fasta \
								$imgap_functional_annotation_pfam_db \
								$imgap_functional_annotation_pfam_claninfo_tsv \
								$imgap_additional_threads \
								$imgap_functional_annotation_parallel_hmmsearch_instances \
                                $imgap_functional_annotation_hmmsearch_Z_arg &> $pfam_log
		exit_code=$?
		if [[ $exit_code -ne 0 ]]
		then
			echo "PFAM annotation failed! Aborting!" >&2
			exit $exit_code
		fi
        touch $run_folder/PFAM_DONE
    else
        printf "$(date +%F_%T) - "
        printf "PFAM prediction appears to have been run earlier already. Skipping!\n"
	fi

	cmd_args="$cmd_args -p ${imgap_input_fasta%_*}_pfam.gff"
	hmmsearch_version=$(hmmsearch -h | grep HMMER | cut -d' ' -f2,3)
	if element_not_in_array used_programs "$hmmsearch_version"
	then
		used_programs+=("$hmmsearch_version")
	fi
	pfam_version=$(echo $(basename $(realpath $(dirname \
						$imgap_functional_annotation_pfam_db))))
	used_dbs+=("Pfam $pfam_version")
fi


if [[ $imgap_functional_annotation_cath_funfam_execute == "True" ]]
then
    if [[ ! -f $run_folder/CATH-FUNFAM_DONE ]]
    then
		echo "$(date +%F_%T) - Starting CATH-FUNFAM prediction now..."
        echo "CATH-FUNFAM" >> $run_folder/started_modules.log
		/usr/bin/time $(dirname $0)/hmmsearch_cath_funfams.sh \
					$imgap_input_fasta \
					$imgap_functional_annotation_cath_funfam_db \
					$imgap_additional_threads \
					$imgap_functional_annotation_parallel_hmmsearch_instances \
                    $imgap_functional_annotation_hmmsearch_Z_arg &> $cath_funfam_log
		exit_code=$?
		if [[ $exit_code -ne 0 ]]
		then
			echo "CATH-FUNFAM annotation failed! Aborting!" >&2
			exit $exit_code
		fi
        touch $run_folder/CATH-FUNFAM_DONE
    else
        printf "$(date +%F_%T) - "
        printf "Cath-FunFam prediction appears to have been run earlier already. Skipping!\n"
	fi

	cmd_args="$cmd_args -f ${imgap_input_fasta%_*}_cath_funfam.gff"
	hmmsearch_version=$(hmmsearch -h | grep HMMER | cut -d' ' -f2,3)
	if element_not_in_array used_programs "$hmmsearch_version"
	then
		used_programs+=("$hmmsearch_version")
	fi
	funfam_version=$(echo $(basename $(realpath $(dirname \
						$imgap_functional_annotation_cath_funfam_db))))
	used_dbs+=("Cath-Funfam $funfam_version")
fi

# SignalP and TMHMM will only get executed on isolates
if [[ $imgap_project_type == "isolate" ]]
then
	if [[ "$imgap_functional_annotation_signalp_execute" == "True" ]]
	then
        if [[ ! -f $run_folder/SIGNALP_DONE ]]
        then
			if [[ ! -z "$imgap_functional_annotation_signalp_gram_stain" ]]
			then
				imgap_functional_annotation_signalp_gram_stain="gram-"
			fi
			echo "$(date +%F_%T) - Starting cleavage sites prediction now..."
            echo "SignalP" >> $run_folder/started_modules.log
			/usr/bin/time $(dirname $0)/signal_peptide_cleavage_sites.sh \
							$imgap_input_fasta \
							$imgap_functional_annotation_signalp_gram_stain \
							&> $signalp_log
			exit_code=$?
			if [[ $exit_code -ne 0 ]]
			then
				echo "Cleavage sites annotation failed! Aborting!" >&2
				exit $exit_code
			fi
            touch $run_folder/SIGNALP_DONE
        else
            printf "$(date +%F_%T) - "
            printf "SignalP prediction appears to have been run earlier already. Skipping!\n"
		fi

		cmd_args="$cmd_args -e ${imgap_input_fasta%_*}_cleavage_sites.gff"
		signalp_version=$(signalp -V)
		used_programs+=("$signalp_version")
    fi


	if [[ "$imgap_functional_annotation_tmhmm_execute" == "True" ]]
	then
        if [[ ! -f $run_folder/TMHMM_DONE ]]
        then
			printf "$(date +%F_%T) - "
			printf "Starting transmembrane protein prediction now...\n"
            echo "TMHMM" >> $run_folder/started_modules.log
			/usr/bin/time $(dirname $0)/transmembrane_protein_topology.sh \
									$imgap_input_fasta \
									$imgap_functional_annotation_tmhmm_model \
									&> $tmh_log
			exit_code=$?
			if [[ $exit_code -ne 0 ]]
			then
				printf "Transmembrane protein topology annotation failed!" >&2
				printf " Aborting!\n" >&2
				exit $exit_code
			fi
            touch $run_folder/TMHMM_DONE
        else
            printf "$(date +%F_%T) - "
            printf "TMHMM prediction appears to have been run earlier already. Skipping!\n"
		fi

		cmd_args="$cmd_args -r ${imgap_input_fasta%_*}_tmh.gff"
		decodeanhmm_version=$(decodeanhmm 2>&1 | head -1)
		used_programs+=("$decodeanhmm_version")
	fi
fi


if [[ -f ${imgap_input_fasta%_*}_structural_annotation.gff ]]
then
	echo "$(date +%F_%T) - Assigning product names now..."
    echo "Product Name Assignment" >> $run_folder/started_modules.log
	cmd_args="$cmd_args $imgap_functional_annotation_product_names_mapping_dir"
	cmd_args="$cmd_args ${imgap_input_fasta%_*}_structural_annotation.gff"
	/usr/bin/time $(dirname $0)/assign_product_names_and_create_fa_gff.py \
						$cmd_args &> ${imgap_input_fasta%_*}_product_names.log
	exit_code=$?
	if [[ $exit_code -ne 0 ]]
	then
		echo "Product name assignment failed! Aborting!" >&2
		exit $exit_code
	fi
else
	printf "Could not find the structural annotation file "
	printf "(${imgap_input_fasta%_*}_structural_annotation.gff).\n"
	printf "Thus skipping product name assignment and "
	printf "functional annotation gff file creation.\n"
fi

echo "$(date +%F_%T) - Done with functional annotation."

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

if [[ ${#used_dbs[@]} -gt 0 ]]
then
	echo "Used DBs:"
	printf "${used_dbs[0]}"
	for i in $(seq 1 $((${#used_dbs[@]}-1)))
	do
		printf "; ${used_dbs[$i]}"
	done
	printf "\n"
fi

