workflow s_annotate {

  File    imgap_input_fasta
  String  imgap_project_id

    call pre_qc {
      input:
        input_fasta = imgap_input_fasta,
        project_id = imgap_project_id
    }
}

task pre_qc {

  File   input_fasta
  String project_id

  command <<<
    tmp_fasta="${input_fasta}.tmp"
    qced_fasta="${project_id}_contigs.fna"
	echo Found it $qced_fasta and $tmp_fasta
  >>>

  runtime {
    time: "00:40:00"
    mem: "5G"
    poolname: "marcel_split1"
    node: 1
    nwpn: 1
    docker: "jfroula/img-omics:0.1.1"
  }
	
  output {
    #File fasta = "${project_id}_contigs.fna"
	String out = stdout()
  }
}
