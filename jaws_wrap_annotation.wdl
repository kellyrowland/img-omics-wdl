workflow img_cloud {
    File input_file
    File sub_json
    File databases="s3://bf-20191119-staging/databases.20201029.tar"
    Int mb_split
    String img_container="bfoster1/img-omics:0.1.5"
    call split {input: infile=input_file, container=img_container, blocksize=mb_split}
    scatter(pathname in split.files) {
        call img_annot{input: infile=pathname, all_infile=input_file, db=databases, container=img_container, sub_j=sub_json}
    }
}


task split{
     File infile
     String container
     String blocksize=50
     String tmp_dir="."
     String file_of_files = "splits_out.fof"
      runtime {
        docker: container
        poolname: "somethingunique"
        shared: 0
        node: 1
        nwpn: 1
        mem: "115G"
        time: "10:00:00"
    }
    
     command{
        git clone https://github.com/kellyrowland/img-omics-wdl.git -b cloud
	img-omics-wdl/split.py ${infile} ${blocksize} ${tmp_dir};
	ls 
     }

   output{
       Array[File] files = glob("*_Ga*_contigs.fna")
#       Array[File] files = [glob("1_Ga*_contigs.fna")[0],glob("2_Ga*_contigs.fna")[0]]
  }
}


task img_annot{
     File infile
     File all_infile
     File sub_j
     File db
     String container
     String outfile = "cromwell-executions.tar.gz"
     String outmetadata = "metadata.ouput.json"
      runtime {
        docker: container
        poolname: "somethingunique"
        shared: 0
        node: 1
        nwpn: 1
        mem: "115G"
        time: "10:00:00"
    }
     

    command{
	git clone https://github.com/kellyrowland/img-omics-wdl.git -b cloud &&	cd img-omics-wdl
	tar -xf ${db}
	mkdir -p splits/1
	cp ${infile} splits/1/GaXXXXXXX_contigs.fna
	cp ${all_infile} splits/1/GaXXXXXXX_allcontigs.fna
	java -Dconfig.file=local.conf -jar /opt/omics/bin/cromwell.jar run -i ${sub_j} annotation.wdl -m ${outmetadata} || exit 1
	find cromwell-executions -size +1G |  grep "/inputs/" |  xargs  rm 
	tar -cpvzf ${outfile} cromwell-executions

     }
     output{
	File fileout = outfile
	File metadataout = outmetadata
     }
}

