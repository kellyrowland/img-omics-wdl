import "sub.wdl" as sa

workflow annotation {

  Int     num_splits
  File    imgap_input_dir
  String  imgap_input_fasta
  String  imgap_project_id

  call setup {
    input:
      n_splits = num_splits,
      dir = imgap_input_dir
  }

  scatter(split in setup.splits) {

      call sa.s_annotate {
        input:
          imgap_project_id = imgap_project_id,
          imgap_input_fasta = "${split}"+"/"+"${imgap_input_fasta}"
      }
  }
}

task setup {
  File dir
  Int  n_splits

  command {
    python -c 'for i in range(${n_splits}): print("${dir}/"+str(i+1)+"/")'
  }

  runtime {
    time: "00:40:00"
    mem: "5G"
    poolname: "marcel_split1"
    node: 1
    nwpn: 1
    docker: "jfroula/img-omics:0.1.1"
  }

  output {
    Array[File] splits = read_lines(stdout())
  }
}
