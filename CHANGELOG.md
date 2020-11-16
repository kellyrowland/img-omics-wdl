# Changelog

### 11/16/2020 version 0.1.6  
I made the following changes to Brian’s “cloud” version to make it work on cori. I created a branch jeffs_version and PR: (https://github.com/kellyrowland/img-omics-wdl/pull/8).

1) This line was in the inputs.cloud.json but I had to hard-code it inside crt.wdl because JAWS sees strings with slashes in them '/' as file paths and tries to verify that the file/dir exists. Thus you can't include paths inside a command or paths that exist only inside the container.
"annotation.sa_crt_cli_jar": "java -Xmx1536m -jar /opt/omics/bin/CRT-CLI.jar"


2) removed all paths pointing to scripts and files inside the container (for the reason outlined above).


3) added docker to runtime block


4) added full path to /opt/omics/bin/${model} inside of the task gmmeta in genemark.wdl


5) There was a line in all the hmmsearch scripts where the hmmsearch command was getting built. This redirect: 1> /dev/null code was part of the steps when building the command. However, when running this command, the bash shell saw the 1> as another argument to hmmsearch since it had single quotes around it (i.e. hmmsearch -i <whatever> —dblout ‘>1’ /dev/null). Therefore, I had to move the redirect code to to where the command gets run (i.e on the next line:  $hmmsearchcmd 1> /dev/null).


6) created a new dockerfile bfoster1/img-omics:0.1.6 that included the changes in the hmmsearch scripts.

Files where I had to move the redirect code:
    hmmsearch_cath_funfams.sh
    hmmsearch_tigrfams.sh
    hmmsearch_smart.sh
    hmmsearch_supfams.sh
    hmmsearch_cogs.sh
7) This task (ko_ec) was turned off in the inputs.json since I got the following error:
    `information!ERROR: First line of blastab file does not contain the version information!`
