import "structural-annotation.wdl" as sa

workflow metagenome_annotation {

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
  Boolean sa_trnascan_se_execute
  File    sa_trnascan_se_bin
  Boolean sa_rfam_execute
  File    sa_rfam_bin
  File    sa_rfam_cm
  File    sa_rfam_claninfo_tsv
  File    sa_rfam_feature_lookup_tsv
  Boolean sa_crt_execute
  File    sa_crt_bin
  Boolean sa_prodigal_execute
  File    sa_prodigal_bin
  Boolean sa_genemark_execute
  File    sa_genemark_bin
  File    sa_gff_merge_bin
  File    sa_fasta_merge_bin
  Boolean sa_gff_and_fasta_stats_execute
  File    sa_gff_and_fasta_stats_bin
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
      n_splits = num_splits,
      dir = imgap_input_dir,
      fasta = imgap_input_fasta
  }

  scatter(split in setup.splits) {

    if(sa_execute) {
      call sa.annotate {
        input:
          imgap_project_id = imgap_project_id,
          additional_threads = additional_threads,
          imgap_project_type = imgap_project_type,
          imgap_input_fasta = split,
          bin = sa_bin,
          pre_qc_execute = sa_pre_qc_execute,
          pre_qc_bin = sa_pre_qc_bin,
          pre_qc_rename = sa_pre_qc_rename,
          post_qc_bin = sa_post_qc_bin,
          trnascan_se_execute = sa_trnascan_se_execute,
          trnascan_se_bin = sa_trnascan_se_bin,
          rfam_execute = sa_rfam_execute,
          rfam_bin = sa_rfam_bin,
          rfam_cm = sa_rfam_cm,
          rfam_claninfo_tsv = sa_rfam_claninfo_tsv,
          rfam_feature_lookup_tsv = sa_rfam_feature_lookup_tsv,
          crt_execute = sa_crt_execute,
          crt_bin = sa_crt_bin,
          prodigal_execute = sa_prodigal_execute,
          prodigal_bin = sa_prodigal_bin,
          genemark_execute = sa_genemark_execute,
          genemark_bin = sa_genemark_bin,
          gff_merge_bin = sa_gff_merge_bin,
          fasta_merge_bin = sa_fasta_merge_bin,
          gff_and_fasta_stats_execute = sa_gff_and_fasta_stats_execute,
          gff_and_fasta_stats_bin = sa_gff_and_fasta_stats_bin
      }
    }

#  if(fa_execute) {
#    call functional_annotation {
#      input:
#        fa_bin = fa_bin
#    }
#  }
  }
}

task setup {
  String dir
  Int    n_splits
  String fasta

  command {
    python -c 'for i in range(${n_splits}): print("${dir}"+str(i+1)+"/${fasta}")'
  }
  output {
    Array[File] splits = read_lines(stdout())
  }
}
