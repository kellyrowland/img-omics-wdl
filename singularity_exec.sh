#!/bin/bash
export SINGULARITY_CACHEDIR=$(pwd)/sif_files
export SINGULARITY_PULLFOLDER=$(pwd)/sif_files
export SINGULARITY_TMPDIR=$(pwd)/sif_files
export SINGULARITY_LOCALCACHEDIR=$(pwd)/sif_files
singularity exec --bind $1:$2 --bind /global/scratch/jaws/ref_data:/refdata docker://$3 $4 $5
