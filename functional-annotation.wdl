workflow f_annotate {

  String  imgap_project_id
  String  imgap_project_type
  Int     additional_threads
  File    input_fasta
  Boolean ko_ec_execute
  File    ko_ec_img_nr_db
  File    ko_ec_md5_mapping
  File    ko_ec_taxon_to_phylo_mapping
  File    lastal_bin
  File    selector_bin

  if(ko_ec_execute) {
    call ko_ec {
      input:
        project_id = imgap_project_id,
        project_type = imgap_project_type,
        input_fasta = input_fasta,
        threads = additional_threads,
        nr_db = ko_ec_img_nr_db,
        md5 = ko_ec_md5_mapping,
        phylo = ko_ec_taxon_to_phylo_mapping,
        lastal = lastal_bin,
        selector = selector_bin
    }
  }
}

task ko_ec {

  String project_id
  String project_type
  Int    threads = 0
  File   input_fasta
  File   nr_db
  File   md5
  File   phylo
  Int    top_hits = 5
  Int    min_ko_hits = 7
  Float  aln_length_ratio = 0.7
  File   lastal
  File   selector


  command {
    ${lastal} -f blasttab+ -P ${threads} ${nr_db} ${input_fasta} 1> ${project_id}_proteins.img_nr.last.blasttab
    ${selector} -l ${aln_length_ratio} -m ${min_ko_hits} -n ${top_hits} \
                ${project_type} ${md5} ${phylo} \
                ${project_id}_ko.tsv ${project_id}_ec.tsv \
                ${project_id}_gene_phylogeny.tsv ${project_id}_ko_ec.gff \
                < ${project_id}_proteins.img_nr.last.blasttab
  }
  output {
    File last_blasttab = "${project_id}_proteins.img_nr.last.blasttab"
  }
}