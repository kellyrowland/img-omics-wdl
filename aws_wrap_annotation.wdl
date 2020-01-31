workflow img_cloud {
    File input_file
    File databases="s3://bf-20191119-staging/databases.tar"
    String img_container="bfoster1/img-omics:0.1.0"
    call split {input: infile=input_file, container=img_container}
    scatter(pathname in split.files) {
        call img_annot{input: infile=pathname,db=databases, container=img_container}
    }
}


task split{
     File infile
     String container
     String blocksize=30000000
     String tmp_dir="."
     runtime {
            docker: container
            backend: "i3-120-16c-spot-ceq"
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
        Array[File] files = glob("*split.faa")
  }
}


task img_annot{
     File infile
     File db
     String container
     String outfile = "/cromwell_root/img-omics-wdl/cromwell-executions.tar.gz"
     String outmetadata = "/cromwell_root/img-omics-wdl/metadata.ouput.json"
     String filename_resources ="resources.log"
     runtime {
            docker: container
            backend: "i3-120-16c-spot-ceq"
            memory: "120 GiB"
            cpu:  16
            maxRetries: 1
     }
    command{
	echo $(curl --fail --max-time 10 --silent http://169.254.169.254/latest/meta-data/public-hostname)
        touch ${filename_resources};
	curl --fail --max-time 10 --silent https://bitbucket.org/berkeleylab/jgi-meta/get/master.tar.gz | tar --wildcards -zxvf - "*/bin/resources.bash" && ./*/bin/resources.bash > ${filename_resources} &
        sleep 30
	
	git clone https://github.com/kellyrowland/img-omics-wdl.git -b cloud &&	cd img-omics-wdl
	ls /cromwell_root;
	cd /cromwell_root; ls ; tar -xf ${db} ;cd -
	mkdir -p splits/1
	cp ${infile} splits/1/GaXXXXXXX_contigs.fna
	java -jar /opt/omics/bin/cromwell.jar run -i inputs.cloud.json annotation.wdl -m ${outmetadata}
	tar -cpvzf ${outfile} /cromwell_root/img-omics-wdl/cromwell-executions

     }
     output{
	File fileout = outfile
	File metadataout = outmetadata
	File resourceout = filename_resources
     }
}

