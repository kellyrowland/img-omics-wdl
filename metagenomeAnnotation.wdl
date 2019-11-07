workflow metagenomeAnnotation {

  Int    num_splits
  String imgap_input_dir
  String imgap_input_fasta
  String imgap_project_id
  String imgap_project_type
  Int additional_threads
  # structural annotation
  Boolean sa_execute
  File    sa_bin
  Boolean sa_pre_qc_execute
  File    sa_pre_qc_bin
  String  sa_pre_qc_rename
  File    sa_post_qc_bin
  Boolean sa_trnascan_se
  Boolean sa_rfam_execute
  File    sa_rfam_cm
  File    sa_rfam_claninfo_tsv
  File    sa_rfam_feature_lookup_tsv
  Boolean sa_crt_execute
  Boolean sa_prodigal_execute
  Boolean sa_genemark_execute
  Boolean sa_gff_and_fasta_stats_execute
  # functional annotation
  Boolean fa_execute
  File    fa_bin
  String  fa_product_names_mapping_dir
  Boolean fa_ko_ec_execute
  File    fa_ko_ec_img_nr_db
  File    fa_ko_ec_md5_mapping
  File    fa_ko_ec_taxon_to_phylo_mapping
  Boolean fa_cath_funfam_execute
  File    fa_cath_funfam_db
  Boolean fa_pfam_execute
  File    fa_pfam_db
  File    fa_pfam_claninfo_tsv
  Boolean fa_superfamily_excute
  File    fa_superfamily_db
  Boolean fa_cog_execute
  File    fa_cog_db
  Boolean fa_tigrfam_execute
  File    fa_tigrfam_db
  Boolean fa_smart_execute
  File    fa_smart_db
  Boolean fa_signalp_execute
  String  fa_signalp_gram_stain
  Boolean fa_tmhmm_execute
  File    fa_tmhmm_model

  call setup {
    input:
      splits=num_splits
  }

  scatter(n in setup.n_splits) {
#    call annotate {
#      input:
#        split_index=split,
#    }
     call test {
       input:
         val=n
     }
  }

  # structural annotation
  if(sa_execute && sa_pre_qc_execute) {
 #   call pre_qc {
 #     input:
 #       qc_bin = sa_pre_qc_bin,
 #       project_type = imgap_project_type,
 #       input_fasta = imgap_input_fasta,
 #       project_id = imgap_project_id,
 #       rename = sa_pre_qc_rename
 #   }
  }

  if(sa_execute) {
 #   call structural_annotation {
 #     input:
 #       sa_bin = sa_bin
 #   }
  }

  if(sa_execute && imgap_project_type == "isolate") {
#    call post_qc {
#      input:
#        qc_bin = sa_post_qc_bin
#    }
  }

  if(fa_execute) {
#    call functional_annotation {
#      input:
#        fa_bin = fa_bin
#    }
  }
}

task test {
  Int val

  command {
    echo ${val}
  }
}

task setup {
  String conda_env="imgap"
  Int splits

  command {
    module load python/3.7-anaconda-2019.07 hmmer/3.1b2 parallel
    source activate ${conda_env}
    python -c 'for i in range(${splits}): print(i+1)'
  }
  output {
    Array[Int] n_splits = read_lines(stdout())
  }
}

task pre_qc {

  File qc_bin
  String project_type
  File input_fasta
  String project_id
  String rename = "yes"

  command {
    ${qc_bin} ${project_type} ${input_fasta} ${project_id} ${rename}
  }
  output {
    File fasta = "${project_id}_contigs.fna"
  }
}

task structural_annotation {

  File sa_bin 

  command {
    ${sa_bin}
  }
}

task post_qc {

  File qc_bin
  File input_fasta

  command {
    ${qc_bin} ${input_fasta} "${input_fasta}_structural_annotation.gff"
  }
  output {
    File out = "${input_fasta}_structural_annotation.gff"
  }
}

task functional_annotation {
  File    fa_bin
  String  product_names_mapping_dir
  Boolean ko_ec_execute
  File    ko_ec_img_nr_db
  File    ko_ec_md5_mapping
  File    ko_ec_taxon_to_phylo_mapping
  Boolean cath_funfam_execute
  File    cath_funfam_db
  Boolean pfam_execute
  File    pfam_db
  File    pfam_claninfo_tsv
  Boolean superfamily_excute
  File    superfamily_db
  Boolean cog_execute
  File    cog_db
  Boolean tigrfam_execute
  File    tigrfam_db
  Boolean smart_execute
  File    smart_db
  Boolean signalp_execute
  String  signalp_gram_stain
  Boolean tmhmm_execute
  File    tmhmm_model

  command {
    ${fa_bin}
  }
}

