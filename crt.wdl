workflow crt {

  String imgap_input_fasta
  String imgap_project_id
  String output_dir
  String crt_transform_bin

  call run {
    input:
      input_fasta = imgap_input_fasta,
      project_id = imgap_project_id,
      out_dir = output_dir
  }

  call transform {
    input:
      transform_bin = crt_transform_bin,
      project_id = imgap_project_id,
      crt_out = run.out,
      out_dir = output_dir
  }

  output {
    File crisprs = transform.crisprs
    File gff = transform.gff
  }
}

task run {

  File   input_fasta
  String project_id
  String out_dir

  command {
	jar="java -Xmx1536m -jar /opt/omics/bin/CRT-CLI.jar"
    $jar ${input_fasta} ${project_id}_crt.out
    #cp ./${project_id}_crt.out ${out_dir}
  }

  runtime {
	docker: "jfroula/img-omics:0.1.7"
    cluster: "cori"
    time: "1:00:00"
    mem: "86G"
    poolname: "justtest"
    shared: 1
    node: 1
    nwpn: 1
    constraint: "haswell"
  }

  output {
    File out = "${project_id}_crt.out"
  }
}

task transform {

  String transform_bin
  File   crt_out
  String project_id
  String crt_out_local = basename(crt_out)
  String out_dir

  command {
    mv ${crt_out} ./${crt_out_local}
	jar="java -Xmx1536m -jar /opt/omics/bin/CRT-CLI.jar"
    tool_and_version=$($jar -version | cut -d' ' -f1,6)
    ${transform_bin} ${crt_out_local} "$tool_and_version"
    #cp -r ./${project_id}_crt.* ${out_dir}
  }

  runtime {
	docker: "jfroula/img-omics:0.1.7"
    cluster: "cori"
    time: "1:00:00"
    mem: "86G"
    poolname: "justtest"
    shared: 1
    node: 1
    nwpn: 1
    constraint: "haswell"
  }

  output{
    File crisprs = "${project_id}_crt.crisprs"
    File gff = "${project_id}_crt.gff"
  }
}

