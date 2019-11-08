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
      splits=num_splits
  }

  scatter(n in setup.n_splits) {
#     call test {
#       input:
#         val=n,
#         dir=imgap_input_dir,
#         fasta=imgap_input_fasta
#     }
    # structural annotation
    if(sa_execute && sa_pre_qc_execute) {
      call pre_qc {
        input:
          bin = sa_pre_qc_bin,
          project_type = imgap_project_type,
          val = n,
          dir = imgap_input_dir,
          input_fasta = imgap_input_fasta,
          project_id = imgap_project_id,
          rename = sa_pre_qc_rename
      }
    }
    if(sa_execute && sa_trnascan_se_execute) {
      call trnascan_se {
        input:
          bin = sa_trnascan_se_bin,
          input_fasta = imgap_input_fasta,
          project_type = imgap_project_type,
          threads = additional_threads
      }
    }
    if(sa_execute && sa_rfam_execute) {
      call rfam {
        input:
          bin = sa_rfam_bin,
          input_fasta = imgap_input_fasta,
          cm = sa_rfam_cm,
          claninfo_tsv = sa_rfam_claninfo_tsv,
          feature_lookup_tsv = sa_rfam_feature_lookup_tsv,
          threads = additional_threads
      }
    }
    if(sa_execute && sa_crt_execute) {
      call crt {
        input:
          bin = sa_crt_bin,
          input_fasta = imgap_input_fasta
      }
    }
    if(sa_execute && sa_prodigal_execute) {
      call prodigal {
        input:
          bin = sa_prodigal_bin,
          input_fasta = imgap_input_fasta,
          project_type = imgap_project_type,
      }
    }
    if(sa_execute && sa_genemark_execute) {
      call genemark {
        input:
          bin = sa_genemark_bin,
          input_fasta = imgap_input_fasta,
          project_type = imgap_project_type,
      }
    }
    if(sa_execute) {
      call gff_merge {
        input:
          bin = sa_gff_merge_bin,
          input_fasta = imgap_input_fasta
      }
    }
    if(sa_prodigal_execute || sa_genemark_execute) {
      call fasta_merge {
        input:
          bin = sa_fasta_merge_bin,
          input_fasta = imgap_input_fasta,
          final_gff = gff_merge.final_gff
      }
    }
    if(sa_execute && sa_gff_and_fasta_stats_execute) {
      call gff_and_fasta_stats {
        input:
          bin = sa_gff_and_fasta_stats_bin,
          input_fasta = imgap_input_fasta,
          final_gff = gff_merge.final_gff
      }
    }
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
  String dir
  String fasta

  command {
    head -n 5 ${dir}${val}/${fasta}
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

  File   bin
  String project_type
  Int    val
  String dir
  File   input_fasta
  String project_id
  String rename = "yes"

  command {
    ${bin} ${project_type} ${dir}${val}/${input_fasta} ${project_id} ${rename}
  }
  output {
    File fasta = "${project_id}_contigs.fna"
  }
}

task trnascan_se {

  File   bin
  File   input_fasta
  String project_type
  Int    threads

  command {
    ${bin} ${input_fasta} ${project_type} ${threads}
  }
  output {
    File trna_log = stdout()
  }
}

task rfam {

  File bin
  File input_fasta
  File cm
  File claninfo_tsv
  File feature_lookup_tsv
  Int  threads

  command {
    ${bin} ${input_fasta} ${cm} ${claninfo_tsv} ${feature_lookup_tsv}
  }
  output {
    File rfam_log = stdout()
  }
}

task crt {

  File bin
  File input_fasta

  command {
    ${bin} ${input_fasta}
  }
  output {
    File crt_log = stdout()
  }
}

task prodigal {

  File   bin
  File   input_fasta
  String project_type

  command {
    ${bin} ${input_fasta} ${project_type}
  }
  output {
    File prodigal_log = stdout()
  }
}

task genemark {

  File   bin
  File   input_fasta
  String project_type

  command {
    ${bin} ${input_fasta} ${project_type}
  }
  output {
    File genemark_log = stdout()
  }
}

task gff_merge {

  File  bin
  File  input_fasta
  File? misc_and_regulatory_gff
  File? rrna_gff
  File? trna_gff
  File? ncrna_tmrna_gff
  File? crt_gff
  File? genemark_gff
  File? prodigal_gff

  command {
    ${bin} -f ${input_fasta} ${"-a " + misc_and_regulatory_gff + " " + rrna_gff} \
    ${trna_gff} ${ncrna_tmrna_gff} ${crt_gff} ${genemark_gff} ${prodigal_gff}
  }
  output {
    File final_gff = stdout()
  }
}

task fasta_merge {

  File  bin
  File  input_fasta
  File  final_gff
  File? genemark_genes
  File? genemark_proteins
  File? prodigal_genes
  File? prodigal_proteins

  command {
    ${bin} ${final_gff} ${genemark_genes} ${prodigal_genes} 1> ${input_fasta}_genes.fna
    ${bin} ${final_gff} ${genemark_proteins} ${prodigal_proteins} 1> ${input_fasta}_proteins.faa
  }
  output {
    File final_genes = "basename(input_fasta)_genes.fna"
    File final_proteins = "basename(input_fasta)_proteins.faa"
  }
}

task gff_and_fasta_stats {

  File bin
  File input_fasta
  File final_gff

  command {
    ${bin} ${input_fasta} ${final_gff}
  }
}

task post_qc {

  File qc_bin
  File input_fasta

  command {
    ${qc_bin} ${input_fasta} "${input_fasta}_structural_annotation.gff"
  }
  output {
    File out = "basename(input_fasta)_structural_annotation.gff"
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

