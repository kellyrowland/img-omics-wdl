import "trnascan.wdl" as trnascan
import "rfam.wdl" as rfam
import "crt.wdl" as crt
import "prodigal.wdl" as prodigal
import "genemark.wdl" as genemark

workflow s_annotate {

  File    imgap_input_fasta
  String  imgap_project_id
  String  imgap_project_type
  String  output_dir
  Int     additional_threads
  Boolean pre_qc_execute
  File    pre_qc_bin
  String  pre_qc_rename
  File    post_qc_bin
  Boolean trnascan_se_execute
  File    trnascan_se_bin
  File    trnascan_pick_and_transform_to_gff_bin
  Boolean rfam_execute
  File    rfam_cmsearch_bin
  File    rfam_clan_filter_bin
  File    rfam_cm
  File    rfam_claninfo_tsv
  File    rfam_feature_lookup_tsv
  Boolean crt_execute
  File    crt_cli_jar
  File    crt_transform_bin
  Boolean prodigal_execute
  File    prodigal_bin
  File    unify_bin
  Boolean genemark_execute
  File    genemark_iso_bin
  File    genemark_meta_bin
  File    genemark_meta_model
  File    gff_merge_bin
  File    fasta_merge_bin
  Boolean gff_and_fasta_stats_execute
  File    gff_and_fasta_stats_bin

  if(pre_qc_execute) {
    call pre_qc {
      input:
        bin = pre_qc_bin,
        project_type = imgap_project_type,
        input_fasta = imgap_input_fasta,
        project_id = imgap_project_id,
        rename = pre_qc_rename
    }
  }
  if(trnascan_se_execute) {
    call trnascan.trnascan {
      input:
        trnascan_se_bin = trnascan_se_bin,
        pick_and_transform_to_gff_bin = trnascan_pick_and_transform_to_gff_bin,
        imgap_input_fasta = imgap_input_fasta,
        imgap_project_id = imgap_project_id,
        imgap_project_type = imgap_project_type,
        additional_threads = additional_threads,
        output_dir = output_dir
    }
  }
  if(rfam_execute) {
    call rfam.rfam {
      input:
        cmsearch_bin = rfam_cmsearch_bin,
        clan_filter_bin = rfam_clan_filter_bin,
        imgap_input_fasta = imgap_input_fasta,
        imgap_project_id = imgap_project_id,
        imgap_project_type = imgap_project_type,
        cm = rfam_cm,
        claninfo_tsv = rfam_claninfo_tsv,
        feature_lookup_tsv = rfam_feature_lookup_tsv,
        additional_threads = additional_threads,
        output_dir = output_dir
    }
  }
  if(crt_execute) {
    call crt.crt {
      input:
        crt_cli_jar = crt_cli_jar,
        crt_transform_bin = crt_transform_bin,
        imgap_input_fasta = imgap_input_fasta,
        imgap_project_id = imgap_project_id,
        output_dir = output_dir
    }
  }
  if(prodigal_execute) {
    call prodigal.prodigal {
      input:
        prodigal_bin = prodigal_bin,
        prodigal_unify_bin = unify_bin,
        imgap_input_fasta = imgap_input_fasta,
        imgap_project_id = imgap_project_id,
        imgap_project_type = imgap_project_type,
        output_dir = output_dir
    }
  }
  if(genemark_execute) {
    call genemark.genemark {
      input:
        genemark_iso_bin = genemark_iso_bin,
        genemark_meta_bin = genemark_meta_bin,
        genemark_meta_model = genemark_meta_model,
        genemark_unify_bin = unify_bin,
        imgap_input_fasta = imgap_input_fasta,
        imgap_project_id = imgap_project_id,
        imgap_project_type = imgap_project_type,
        output_dir = output_dir
    }
  }
  call gff_merge {
    input:
      bin = gff_merge_bin,
      input_fasta = imgap_input_fasta,
      project_id = imgap_project_id,
      misc_and_regulatory_gff = rfam.misc_bind_misc_feature_regulatory_gff,
      rrna_gff = rfam.rrna_gff,
      trna_gff = trnascan.gff,
      ncrna_tmrna_gff = rfam.ncrna_tmrna_gff,
      crt_gff = crt.gff, 
      genemark_gff = genemark.gff,
      prodigal_gff = prodigal.gff,
      output_dir = output_dir
  }
  if(prodigal_execute || genemark_execute) {
    call fasta_merge {
      input:
        bin = fasta_merge_bin,
        input_fasta = imgap_input_fasta,
        project_id = imgap_project_id,
        final_gff = gff_merge.final_gff,
        genemark_genes = genemark.genes,
        genemark_proteins = genemark.proteins,
        prodigal_genes = prodigal.genes,
        prodigal_proteins = prodigal.proteins,
        output_dir = output_dir
    }
  }
  if(gff_and_fasta_stats_execute) {
    call gff_and_fasta_stats {
      input:
        bin = gff_and_fasta_stats_bin,
        input_fasta = imgap_input_fasta,
        project_id = imgap_project_id,
        final_gff = gff_merge.final_gff
    }
  }
  if(imgap_project_type == "isolate") {
    call post_qc {
      input:
        qc_bin = post_qc_bin,
        input_fasta = imgap_input_fasta,
        project_id = imgap_project_id,
        output_dir = output_dir
    }
  }
  output {
    File  gff = gff_merge.final_gff
    File? proteins = fasta_merge.final_proteins
  }
}

task pre_qc {

  File   bin
  String project_type
  File   input_fasta
  String project_id
  String rename = "yes"

  command {
    ${bin} ${project_type} ${input_fasta} ${project_id} ${rename}
  }
  output {
    File fasta = "${project_id}_contigs.fna"
  }
}

task gff_merge {

  File   bin
  File   input_fasta
  String project_id
  File?  misc_and_regulatory_gff
  File?  rrna_gff
  File?  trna_gff
  File?  ncrna_tmrna_gff
  File?  crt_gff
  File?  genemark_gff
  File?  prodigal_gff
  String output_dir

  command {
    ${bin} -f ${input_fasta} ${"-a " + misc_and_regulatory_gff + " " + rrna_gff} \
    ${trna_gff} ${ncrna_tmrna_gff} ${crt_gff} \
    ${genemark_gff} ${prodigal_gff} 1> ${project_id}_structural_annotation.gff
    cp ./${project_id}_structural_annotation.gff ${output_dir}
  }
  output {
    File final_gff = "${project_id}_structural_annotation.gff"
  }
}

task fasta_merge {

  File   bin
  File   input_fasta
  String project_id
  File   final_gff
  File?  genemark_genes
  File?  genemark_proteins
  File?  prodigal_genes
  File?  prodigal_proteins
  String output_dir

  command {
    ${bin} ${final_gff} ${genemark_genes} ${prodigal_genes} 1> ${project_id}_genes.fna
    ${bin} ${final_gff} ${genemark_proteins} ${prodigal_proteins} 1> ${project_id}_proteins.faa
    cp ./${project_id}_genes.fna ./${project_id}_proteins.faa ${output_dir}
  }
  output {
    File final_genes = "${project_id}_genes.fna"
    File final_proteins = "${project_id}_proteins.faa"
  }
}

task gff_and_fasta_stats {

  File   bin
  File   input_fasta
  String project_id
  File   final_gff

  command {
    ${bin} ${input_fasta} ${final_gff}
  }
}

task post_qc {

  File   qc_bin
  File   input_fasta
  String project_id
  String output_dir

  command {
    ${qc_bin} ${input_fasta} "${project_id}_structural_annotation.gff"
    cp ./${project_id}_structural_annotation.gff ${output_dir}
  }
  output {
    File out = "${project_id}_structural_annotation.gff"
  }
}
