import "trnascan.wdl" as trnascan

workflow annotate {

  Boolean test = true
  File   imgap_input_fasta
  String imgap_project_id
  String imgap_project_type
  Int additional_threads
  # structural annotation
  File    bin
  Boolean pre_qc_execute
  File    pre_qc_bin
  String  pre_qc_rename
  File    post_qc_bin
  Boolean trnascan_se_execute
  File    trnascan_se_bin
  File    trnascan_pick_and_transform_to_gff_bin
  Boolean rfam_execute
  File    rfam_bin
  File    rfam_cm
  File    rfam_claninfo_tsv
  File    rfam_feature_lookup_tsv
  Boolean crt_execute
  File    crt_bin
  Boolean prodigal_execute
  File    prodigal_bin
  Boolean genemark_execute
  File    genemark_bin
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
        additional_threads = additional_threads
    }
  }
  if(rfam_execute) {
    call rfam {
      input:
        bin = rfam_bin,
        input_fasta = imgap_input_fasta,
        project_id = imgap_project_id,
        cm = rfam_cm,
        claninfo_tsv = rfam_claninfo_tsv,
        feature_lookup_tsv = rfam_feature_lookup_tsv,
        threads = additional_threads
    }
  }
  if(crt_execute) {
    call crt {
      input:
        bin = crt_bin,
        input_fasta = imgap_input_fasta,
        project_id = imgap_project_id
    }
  }
  if(prodigal_execute) {
    call prodigal {
      input:
        bin = prodigal_bin,
        input_fasta = imgap_input_fasta,
        project_id = imgap_project_id,
        project_type = imgap_project_type
    }
  }
  if(genemark_execute) {
    call genemark {
      input:
        bin = genemark_bin,
        input_fasta = imgap_input_fasta,
        project_id = imgap_project_id,
        project_type = imgap_project_type
    }
  }
  if(!test) {
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
      prodigal_gff = prodigal.gff
  }
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
        prodigal_proteins = prodigal.proteins
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
#    call post_qc {
#      input:
#        qc_bin = post_qc_bin
#    }
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

task rfam {

  File   bin
  File   input_fasta
  String project_id
  File   cm
  File   claninfo_tsv
  File   feature_lookup_tsv
  Int    threads

  command {
    ${bin} ${input_fasta} ${cm} ${claninfo_tsv} ${feature_lookup_tsv} ${threads} &> ${project_id}_rfam.log
  }
  output {
    File log = "${project_id}_rfam.log"
    File tbl = "${project_id}_rfam.tbl"
    File misc_bind_misc_feature_regulatory_gff = "${project_id}_rfam_misc_bind_misc_feature_regulatory.gff"
    File ncrna_tmrna_gff = "${project_id}_rfam_ncrna_tmrna.gff"
    File rrna_gff = "${project_id}_rfam_rrna.gff"
  }
}

task crt {

  File   bin
  File   input_fasta
  String project_id

  command {
    ${bin} ${input_fasta}  &> ${project_id}_crt.log
  }
  output {
    File log = "${project_id}_crt.log"
    File crisprs = "${project_id}_crt.crisprs"
    File gff = "${project_id}_crt.gff"
    File out = "${project_id}_crt.out"
  }
}

task prodigal {

  File   bin
  File   input_fasta
  String project_id
  String project_type

  command {
    ${bin} ${input_fasta} ${project_type} &> ${project_id}_prodigal.log
  }
  output {
    File log = "${project_id}_prodigal.log"
    File gff = "${project_id}_prodigal.gff"
    File out = "${project_id}_prodigal.out"
    File genes = "${project_id}_prodigal_genes.fna" 
    File proteins = "${project_id}_prodigal_proteins.faa" 
  }
}

task genemark {

  File   bin
  File   input_fasta
  String project_id
  String project_type

  command {
    ${bin} ${input_fasta} ${project_type} &> ${project_id}_genemark.log
  }
  output {
    File log = "${project_id}_genemark.log"
    File gff = "${project_id}_genemark.gff"
    File genes = "${project_id}_genemark_genes.fna" 
    File proteins = "${project_id}_genemark_proteins.faa" 
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

  command {
    ${bin} -f ${input_fasta} ${"-a " + misc_and_regulatory_gff + " " + rrna_gff} \
    ${trna_gff} ${ncrna_tmrna_gff} ${crt_gff} \
    ${genemark_gff} ${prodigal_gff} 1> ${project_id}_structural_annotation.gff
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

  command {
    ${bin} ${final_gff} ${genemark_genes} ${prodigal_genes} 1> ${project_id}_genes.fna
    ${bin} ${final_gff} ${genemark_proteins} ${prodigal_proteins} 1> ${project_id}_proteins.faa
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

  command {
    ${qc_bin} ${input_fasta} "${project_id}_structural_annotation.gff"
  }
  output {
    File out = "${project_id}_structural_annotation.gff"
  }
}
