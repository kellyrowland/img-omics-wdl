# IMG Omics Annotation WDL

## Setup

### Requirements

* Cromwell (and associated Java dependencies)
* Input binaries, databases, etc. specified in the `inputs.json` file

### Input

Input options are listed in the `inputs.json` file and include locations of
input databases, locations of necessary executables, and Boolean variables
to toggle run options on and off.

The inputs have some executables' absolute paths within the docker container. 
The scripts (bin), libraries (lib) and supporting script (programs) are all 
under "omics".  And this whole directory is copied to the container.
The original source for these scripts was 
`/global/dna/projectdirs/microbial/omics/gbp/img-annotation-pipeline/bin`

The workflow is written to take in a pre-split fasta file of a metagenome or isolate genome. 
The directory structure looks like:

`/path/to/splits/<split number>/file.fna`

The splits should go from 1 to N, and N should be specified in `inputs.json`
as the variable `metagenome_annotation.num_splits`. The workflow assumes that
all FASTA splits have the same file name. So, for example, if we have an
example FASTA `example.fna` split into 3, the directory structure would look
like:

```
/path/splits/1/example.fna
/path/splits/2/example.fna
/path/splits/3/example.fna
```

where "example" should be set as the value of the variable 
`annotation.imgap_project_id` and `/path/splits/` should be set as
the value of the variable `annotation.imgap_input_dir`.

## Instructions

Edit `inputs.json` to pick whichever parts of the annotation pipeline(s) should
be run. Make sure that the locations of the input databases and executables are
correctly set as well.

For example, on cori the ref databases are
`/global/dna/shared/databases/jaws/img/`

And at lbl (i.e. LAB IT) they are here:
`/global/scratch/jaws/ref_data/img`

The annotation workflow is structured:

```
|annotation
|- setup
|- structural annotation
|-- pre-qc
|-- trnascan
|-- rfam
|-- crt
|-- prodigal
|-- genemark
|-- gff_merge
|-- fasta_merge
|-- gff_and_fasta_stats
|-- post-qc
|- functional annotation
|-- ko_ec
|-- smart
|-- cog
|-- tigrfam
|-- superfamily
|-- pfam
|-- cath_funfam
|-- signalp
|-- tmhmm
|-- product_name_assign
```

Run the workflow with the command:

# using cromwell.jar (for testing)
`java -jar <Cromwell> annotation.wdl -i inputs.json`

# or through jaws
`jaws submit annotation.wdl inputs.json`

in this directory.
