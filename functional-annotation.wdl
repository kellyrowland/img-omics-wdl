workflow f_annotate {

  String  imgap_project_id
  String  imgap_project_type
  Int     additional_threads
  File    input_fasta
  Boolean ko_ec_execute
  String  ko_ec_img_nr_db
  File    ko_ec_md5_mapping
  File    ko_ec_taxon_to_phylo_mapping
  File    lastal_bin
  File    selector_bin
  Boolean smart_execute
  File    smart_db
  File    hmmsearch_bin
  File    frag_hits_filter_bin
  Boolean cog_execute
  File    cog_db
  Boolean tigrfam_execute
  File    tigrfam_db
  File    hit_selector_bin
  Boolean superfam_execute
  File    superfam_db
  Boolean pfam_execute
  File    pfam_db
  File    pfam_claninfo_tsv
  File    pfam_clan_filter
  Boolean cath_funfam_execute
  File    cath_funfam_db
  Boolean signalp_execute
  String  signalp_gram_stain
  File    signalp_bin
  Boolean tmhmm_execute
  File    tmhmm_model
  File    tmhmm_decode
  File    tmhmm_decode_parser

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
  if(smart_execute) {
    call smart {
      input:
        project_id = imgap_project_id,
        input_fasta = input_fasta,
        threads = additional_threads,
        smart_db = smart_db,
        hmmsearch = hmmsearch_bin,
        frag_hits_filter = frag_hits_filter_bin
    }
  }
  if(cog_execute) {
    call cog {
      input:
        project_id = imgap_project_id,
        input_fasta = input_fasta,
        threads = additional_threads,
        cog_db = cog_db,
        hmmsearch = hmmsearch_bin,
        frag_hits_filter = frag_hits_filter_bin
    }
  }
  if(tigrfam_execute) {
    call tigrfam {
      input:
        project_id = imgap_project_id,
        input_fasta = input_fasta,
        threads = additional_threads,
        tigrfam_db = tigrfam_db,
        hmmsearch = hmmsearch_bin,
        hit_selector = hit_selector_bin
    }
  }
  if(superfam_execute) {
    call superfam {
      input:
        project_id = imgap_project_id,
        input_fasta = input_fasta,
        threads = additional_threads,
        superfam_db = superfam_db,
        hmmsearch = hmmsearch_bin,
        frag_hits_filter = frag_hits_filter_bin
    }
  }
  if(pfam_execute) {
    call pfam {
      input:
        project_id = imgap_project_id,
        input_fasta = input_fasta,
        threads = additional_threads,
        pfam_db = pfam_db,
        pfam_claninfo_tsv = pfam_claninfo_tsv,
        pfam_clan_filter = pfam_clan_filter,
        hmmsearch = hmmsearch_bin
    }
  }
  if(cath_funfam_execute) {
    call cath_funfam {
      input:
        project_id = imgap_project_id,
        input_fasta = input_fasta,
        threads = additional_threads,
        cath_funfam_db = cath_funfam_db,
        hmmsearch = hmmsearch_bin,
        frag_hits_filter = frag_hits_filter_bin
    }
  }
  if(imgap_project_type == "isolate" && signalp_execute) {
    call signalp {
      input:
        project_id = imgap_project_id,
        input_fasta = input_fasta,
        gram_stain = signalp_gram_stain,
        signalp = signalp_bin
    }
  }
  if(imgap_project_type == "isolate" && tmhmm_execute) {
    call tmhmm {
      input:
        project_id = imgap_project_id,
        input_fasta = input_fasta,
        model = tmhmm_model,
        decode = tmhmm_decode,
        decode_parser = tmhmm_decode_parser
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
    File ko_tsv = "${project_id}_ko.tsv"
    File ec_tsv = "${project_id}_ec.tsv"
    File phylo_tsv = "${project_id}_gene_phylogeny.tsv"
    File ko_ec_gff = "${project_id}_ko_ec.gff"
  }
}

task smart {
  
  String project_id
  File   input_fasta
  File   smart_db
  Int    threads = 0
  Float  min_domain_eval_cutoff = 0.01
  Float  aln_length_ratio = 0.7
  Float  max_overlap_ratio = 0.1
  File   hmmsearch
  File   frag_hits_filter

  command <<<
    ${hmmsearch} --notextw --domE ${min_domain_eval_cutoff} --cpu ${threads} \
                 --domtblout ${project_id}_proteins.smart.domtblout \
                 ${smart_db} ${input_fasta}
    grep -v '^#' ${project_id}_proteins.smart.domtblout | \
    awk '{print $1,$3,$4,$5,$6,$7,$8,$13,$14,$16,$17,$20,$21}' | \
    sort -k1,1 -k7,7nr -k6,6n | \
    ${frag_hits_filter} -a ${aln_length_ratio} -o ${max_overlap_ratio} \
                        > ${project_id}_smart.gff
  >>>
  output {
    File gff = "${project_id}_smart.gff"
  }
}

task cog {
  
  String project_id
  File   input_fasta
  File   cog_db
  Int    threads = 0
  Float  min_domain_eval_cutoff = 0.01
  Float  aln_length_ratio = 0.7
  Float  max_overlap_ratio = 0.1
  File   hmmsearch
  File   frag_hits_filter

  command <<<
    ${hmmsearch} --notextw --domE ${min_domain_eval_cutoff} --cpu ${threads} \
                 --domtblout ${project_id}_proteins.cog.domtblout \
                 ${cog_db} ${input_fasta}
    grep -v '^#' ${project_id}_proteins.cog.domtblout | \
    awk '{print $1,$3,$4,$5,$6,$7,$8,$13,$14,$16,$17,$20,$21}' | \
    sort -k1,1 -k7,7nr -k6,6n | \
    ${frag_hits_filter} -a ${aln_length_ratio} -o ${max_overlap_ratio} \
                        > ${project_id}_cog.gff
  >>>
  output {
    File gff = "${project_id}_cog.gff"
  }
}

task tigrfam {
  
  String project_id
  File   input_fasta
  File   tigrfam_db
  Int    threads = 0
  Float  aln_length_ratio = 0.7
  Float  max_overlap_ratio = 0.1
  File   hmmsearch
  File   hit_selector

  command <<<
    ${hmmsearch} --notextw --cut_nc --cpu ${threads} \
                 --domtblout ${project_id}_proteins.tigrfam.domtblout \
                 ${tigrfam_db} ${input_fasta}
    grep -v '^#' ${project_id}_proteins.tigrfam.domtblout | \
    awk '{print $1,$3,$4,$6,$13,$14,$16,$17,$20,$21}' | \
    sort -k1,1 -k6,6nr -k5,5n | \
    ${hit_selector} -a ${aln_length_ratio} -o ${max_overlap_ratio} \
                        > ${project_id}_tigrfam.gff
  >>>
  output {
    File gff = "{project_id}_tigrfam.gff"
  }
}

task superfam {

  String project_id
  File   input_fasta
  File   superfam_db
  Int    threads = 0
  Float  min_domain_eval_cutoff = 0.01
  Float  aln_length_ratio = 0.7
  Float  max_overlap_ratio = 0.1
  File   hmmsearch
  File   frag_hits_filter

  command <<<
    ${hmmsearch} --notextw --domE ${min_domain_eval_cutoff} --cpu ${threads} \
                 --domtblout ${project_id}_proteins.supfam.domtblout \
                 ${superfam_db} ${input_fasta}
    grep -v '^#' ${project_id}_proteins.supfam.domtblout | \
    awk '{print $1,$3,$4,$5,$6,$7,$8,$13,$14,$16,$17,$20,$21}' | \
    sort -k1,1 -k7,7nr -k6,6n | \
    ${frag_hits_filter} -a ${aln_length_ratio} -o ${max_overlap_ratio} \
                        > ${project_id}_supfam.gff
  >>>
  output {
    File gff = "${project_id}_supfam.gff"
  }
}

task pfam {
  
  String project_id
  File   input_fasta
  File   pfam_db
  File   pfam_claninfo_tsv
  Int    threads = 0
  File   hmmsearch
  File   pfam_clan_filter

  command <<<
    ${hmmsearch} --notextw --cut_tc --cpu ${threads} \
                 --domtblout ${project_id}_proteins.pfam.domtblout \
                 ${pfam_db} ${input_fasta}
    grep -v '^#' ${project_id}_proteins.pfam.domtblout | \
    awk '{print $1,$3,$4,$6,$13,$14,$16,$17,$20,$21}' | \
    sort -k1,1 -k6,6nr -k5,5n | \
    ${pfam_clan_filter} ${pfam_claninfo_tsv} > ${project_id}_pfam.gff
  >>>
  output {
    File gff = "{project_id}_pfam.gff"
  }
}

task cath_funfam {
  
  String project_id
  File   input_fasta
  File   cath_funfam_db
  Int    threads = 0
  Float  min_domain_eval_cutoff = 0.01
  Float  aln_length_ratio = 0.7
  Float  max_overlap_ratio = 0.1
  File   hmmsearch
  File   frag_hits_filter

  command <<<
    ${hmmsearch} --notextw --domE ${min_domain_eval_cutoff} --cpu ${threads} \
                 --domtblout ${project_id}_proteins.cath_funfam.domtblout \
                 ${cath_funfam_db} ${input_fasta}
    grep -v '^#' ${project_id}_proteins.cath_funfam.domtblout | \
    awk '{print $1,$3,$4,$5,$6,$7,$8,$13,$14,$16,$17,$20,$21}' | \
    sort -k1,1 -k7,7nr -k6,6n | \
    ${frag_hits_filter} -a ${aln_length_ratio} -o ${max_overlap_ratio} \
                        > ${project_id}_cath_funfam.gff
  >>>
  output {
    File gff = "${project_id}_cath_funfam.gff"
  }
}

task signalp {
  
  String project_id
  File   input_fasta
  String gram_stain
  File   signalp

  command <<<
    signalp_version=$(${signalp} -V)
    ${signalp} -t ${gram_stain} -f short ${input_fasta} | \
    grep -v '^#' | \
    awk -v sv="$signalp_version" -v ot="${gram_stain}" \
        '$10 == "Y" {print $1"\t"sv"\tcleavage_site\t"$3-1"\t"$3"\t"$2\
        "\t.\t.\tD-score="$9";network="$12";organism_type="ot}' > ${project_id}_cleavage_sites.gff
  >>>
  output {
    File gff = "${project_id}_cleavage_sites.gff"
  }
}

task tmhmm {
  
  String project_id
  File   input_fasta
  File   model
  File   decode
  File   decode_parser

  command <<<
  tool_and_version=$(${decode} -v 2>&1 | head -1)
    background="0.081 0.015 0.054 0.061 0.040 0.068 0.022 0.057 0.056 0.093 0.025"
    background="$background 0.045 0.049 0.039 0.057 0.068 0.058 0.067 0.013 0.032"
    sed 's/\*/X/g' ${input_fasta} | \
    ${decode} -N 1 -background $background -PrintNumbers \
    ${model} 2> /dev/null | ${decode_parser} "$tool_and_version" > ${project_id}_tmh.gff
  >>>
  output {
    File gff = "${project_id}_tmh.gff"
  }
}
