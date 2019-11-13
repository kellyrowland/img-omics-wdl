workflow crt {

  String imgap_input_fasta
  String imgap_project_id
  File   bin

  call crt_task {
    input:
      bin = bin,
      input_fasta = imgap_input_fasta,
      project_id = imgap_project_id
  }

  output {
    File crisprs = crt_task.crisprs
    File gff = crt_task.gff
  }
}

task crt_task {

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

