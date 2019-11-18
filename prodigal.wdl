workflow prodigal {

  String imgap_input_fasta
  String imgap_project_id
  String imgap_project_type
  File   prodigal_bin
  File   prodigal_unify_bin

  if(imgap_project_type == "isolate") {
    call fasta_len {
      input:
        input_fasta = imgap_input_fasta
    }
  }
  if(fasta_len.wc >= 20000) {
    call iso_big {
      input:
        bin = prodigal_bin,
        input_fasta = imgap_input_fasta,
        project_id = imgap_project_id
    }
  }
  if(fasta_len.wc < 20000) {
    call iso_small {
      input:
        bin = prodigal_bin,
        input_fasta = imgap_input_fasta,
        project_id = imgap_project_id
    }
  }
  if(imgap_project_type == "metagenome") {
    call metag {
      input:
        bin = prodigal_bin,
        input_fasta = imgap_input_fasta,
        project_id = imgap_project_id
    }
  }

  call clean_and_unify {
    input:
      iso_big_proteins_fasta = iso_big.proteins,
      iso_small_proteins_fasta = iso_small.proteins,
      meta_proteins_fasta = metag.proteins,
      iso_big_genes_fasta = iso_big.genes,
      iso_small_genes_fasta = iso_small.genes,
      meta_genes_fasta = metag.genes,
      iso_big_gff = iso_big.gff,
      iso_small_gff = iso_small.gff,
      meta_gff = metag.gff,
      unify_bin = prodigal_unify_bin,
      project_id = imgap_project_id
  }

  output {
    File gff = clean_and_unify.gff
    File genes = clean_and_unify.genes
    File proteins = clean_and_unify.proteins
  }
}

task fasta_len {

  File input_fasta

  command {
    grep -v '^>' ${input_fasta} | wc -m
  }
  output {
    Int wc = read_int(stdout())
  }
}

task iso_big {

  File   bin
  File   input_fasta
  Int?   translation_table = 11
  String project_id
  File   train = "${project_id}_prodigal.trn"

  command {
    ${bin} -i ${input_fasta} -t ${train} -g ${translation_table} -q
    ${bin} -f gff -g ${translation_table} -p single -m -i ${input_fasta} \
    -t ${train} -o ${project_id}_prodigal.gff \
    -d ${project_id}_prodigal_genes.fna -a ${project_id}_prodigal_proteins.faa
  }
  output {
    File gff = "${project_id}_prodigal.gff"
    File genes = "${project_id}_prodigal_genes.fna"
    File proteins = "${project_id}_prodigal_proteins.faa"
  }
}

task iso_small {

  File   bin
  File   input_fasta
  String project_id

  command {
    ${bin} -f gff -p meta -m -i ${input_fasta} \
    -o ${project_id}_prodigal.gff -d ${project_id}_prodigal_genes.fna \
    -a ${project_id}_prodigal_proteins.faa
  }
  output {
    File gff = "${project_id}_prodigal.gff"
    File genes = "${project_id}_prodigal_genes.fna"
    File proteins = "${project_id}_prodigal_proteins.faa"
  }
}

task metag {

  File   bin
  File   input_fasta
  String project_id

  command {
    ${bin} -f gff -p meta -m -i ${input_fasta} \
    -o ${project_id}_prodigal.gff -d ${project_id}_prodigal_genes.fna \
    -a ${project_id}_prodigal_proteins.faa
  }
  output {
    File gff = "${project_id}_prodigal.gff"
    File genes = "${project_id}_prodigal_genes.fna"
    File proteins = "${project_id}_prodigal_proteins.faa"
  }
}

task clean_and_unify {
  
  File?  iso_big_proteins_fasta
  File?  iso_small_proteins_fasta
  File?  meta_proteins_fasta
  File?  iso_big_genes_fasta
  File?  iso_small_genes_fasta
  File?  meta_genes_fasta
  File?  iso_big_gff
  File?  iso_small_gff
  File?  meta_gff
  File   unify_bin
  String project_id

  command {
    sed -i 's/\*$//g' ${iso_big_proteins_fasta} ${iso_small_proteins_fasta} ${meta_proteins_fasta}
    sed -i 's/\*/X/g' ${iso_big_proteins_fasta} ${iso_small_proteins_fasta} ${meta_proteins_fasta}
    ${unify_bin} ${iso_big_gff} ${iso_small_gff} ${meta_gff} \
                 ${iso_big_genes_fasta} ${iso_small_genes_fasta} ${meta_genes_fasta} \
                 ${iso_big_proteins_fasta} ${iso_small_proteins_fasta} ${meta_proteins_fasta}
  }
  output {
    File gff = "${project_id}_prodigal.gff"
    File genes = "${project_id}_prodigal_genes.fna"
    File proteins = "${project_id}_prodigal_proteins.faa"
  }
}
