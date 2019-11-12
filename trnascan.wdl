workflow trnascan {

  String imgap_input_fasta
  String imgap_project_id
  String imgap_project_type
  Int    additional_threads
  File   trnascan_se_bin
  File   pick_and_transform_to_gff_bin


  if(imgap_project_type == "isolate") {
    call trnascan_isolate {
      input:
        bin = trnascan_se_bin,
        input_fasta = imgap_input_fasta,
        project_id = imgap_project_id,
        threads = additional_threads
    }
  }
  if(imgap_project_type == "metagenome") {
    call trnascan_metagenome {
      input:
        bin = trnascan_se_bin,
        input_fasta = imgap_input_fasta,
        project_id = imgap_project_id,
        threads = additional_threads
    }
  }
  call pick_and_transform_to_gff {
    input:
      bin = pick_and_transform_to_gff_bin,
      project_id = imgap_project_id,
      bacterial_out = trnascan_isolate.bacterial_out,
      archaeal_out = trnascan_isolate.archaeal_out,
      metagenome_out = trnascan_metagenome.out
  }
  output {
    File gff = pick_and_transform_to_gff.gff
  }
}

task trnascan_isolate {

  File   bin
  String input_fasta
  String project_id
  Int    threads

  command {
    ${bin} -B --thread ${threads} ${input_fasta} &> ${project_id}_trnascan_bacterial.out
    ${bin} -A --thread ${threads} ${input_fasta} &> ${project_id}_trnascan_archaeal.out
  }
  output {
    File bacterial_out = "${project_id}_trnascan_bacterial.out"
    File archaeal_out = "${project_id}_trnascan_archaeal.out"
  }
}

task trnascan_metagenome {

  File   bin
  String input_fasta
  String project_id
  Int    threads

  command <<<
    filesize=$(ls -l ${input_fasta} | awk '{print $5}')
    blocksize=$(($filesize / ${threads-2}))
    cat ${input_fasta} | parallel --pipe --recstart '>' --blocksize $blocksize \
    'cat > 'tmp.$$.split.fna;' ${bin} -G --thread 0 tmp.$$.split.fna &> tmp.general.$$.out;'
    cat tmp.general.* > ${project_id}_trnascan_general.out
  >>>
  output {
    File out = "${project_id}_trnascan_general.out"
  }
}

task pick_and_transform_to_gff {

  File   bin
  String project_id
  File?  bacterial_out
  File?  archaeal_out
  File?  metagenome_out
  
  command {
    ${bin} ${bacterial_out} ${archaeal_out} ${metagenome_out} > ${project_id}_trna.gff
  }
  output {
    File gff = "${project_id}_trna.gff"
  }
}
