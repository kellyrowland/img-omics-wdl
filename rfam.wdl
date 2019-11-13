workflow rfam {

  String imgap_input_fasta
  String imgap_project_id
  String imgap_project_type
  Int    additional_threads
  File   cmsearch_bin
  File   cm
  File   claninfo_tsv
  File   feature_lookup_tsv
  File   clan_filter_bin


call cmsearch {
      input:
        bin = cmsearch_bin,
        input_fasta = imgap_input_fasta,
        project_id = imgap_project_id,
        cm = cm,
        claninfo_tsv = claninfo_tsv,
        feature_lookup_tsv = feature_lookup_tsv,
        threads = additional_threads
    }

call clan_filter {
      input:
        bin = clan_filter_bin,
        project_id = imgap_project_id
    }

  output {
    File rrna_gff = cmsearch.rrna_gff
    File ncrna_tmrna_gff = cmsearch.ncrna_tmrna_gff
    File misc_bind_misc_feature_regulatory_gff = clan_filter.misc_bind_misc_feature_regulatory_gff
  }
}

task cmsearch {

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
    File ncrna_tmrna_gff = "${project_id}_rfam_ncrna_tmrna.gff"
    File rrna_gff = "${project_id}_rfam_rrna.gff"
  }
}

task clan_filter {

  File   bin
  String project_id

  command {
    ${bin}
  }
  output {
    File misc_bind_misc_feature_regulatory_gff = "${project_id}_rfam_misc_bind_misc_feature_regulatory.gff"
  }
}

