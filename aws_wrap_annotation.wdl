workflow img_cloud {
    File input_file
    File sub_json
    File databases="s3://bf-20191119-staging/databases.20201029.tar"
    String img_container="bfoster1/img-omics:0.5.0"
    call split {input: infile=input_file, container=img_container}
    scatter(pathname in split.files) {
        call img_annot{input: infile=pathname, all_infile=input_file, db=databases, container=img_container, sub_j=sub_json}
    }
}


task split{
     File infile
     String container
     String blocksize=20000000
     String tmp_dir="."
     runtime {
            docker: container
	    backend: "i3-120D-ceq"
	    memory: "120 GiB"
            cpu:  16
            maxRetries: 1
     }
    
     command{
           blocksize=${blocksize}
	   tmp_dir=${tmp_dir}
           cat ${infile} | parallel --pipe --recstart '>' \
                             --blocksize $blocksize \
                             'cat > '$tmp_dir'/tmp.$$.split.faa';
     }
     output {
#        Array[File] files = glob("*split.faa")
        Array[File] files = ["tmp.1.split.faa", "tmp.2.split.faa"]
  }
}


task img_annot{
     File infile
     File all_infile
     File sub_j
     File db
     String container
     String outfile = "/cromwell_root/img-omics-wdl/cromwell-executions.tar.gz"
     String outmetadata = "/cromwell_root/img-omics-wdl/metadata.ouput.json"
     String filename_resources ="resources.log"
     runtime {
            docker: container
	    backend: "i3-120D-ceq"
	    memory: "120 GiB"
            cpu:  16
            maxRetries: 1
     }
    command{
	git clone https://github.com/kellyrowland/img-omics-wdl.git -b cloud &&	cd img-omics-wdl
	ls /cromwell_root;
	cd /cromwell_root; ls ; tar -xf ${db} ;cd -
	mkdir -p splits/1
	cp ${infile} splits/1/GaXXXXXXX_contigs.fna
	cp ${all_infile} splits/1/GaXXXXXXX_allcontigs.fna
	java -Dconfig.file=local.conf -jar /opt/omics/bin/cromwell.jar run -i ${sub_j} annotation.wdl -m ${outmetadata}
	find /cromwell_root/img-omics-wdl/cromwell-executions -size +1G |  grep "/inputs/" |  xargs  rm 
	tar -cpvzf ${outfile} /cromwell_root/img-omics-wdl/cromwell-executions

     }
     output{
	File fileout = outfile
	File metadataout = outmetadata
     }
}

