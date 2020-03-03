#!/bin/bash
<<<<<<< HEAD
export SINGULARITY_CACHEDIR=./sif_files
export SINGULARITY_PULLFOLDER=./sif_files
export SINGULARITY_TMPDIR=./sif_files
export SINGULARITY_LOCALCACHEDIR=./sif_files
singularity exec --bind $1:$2 --bind /global/scratch/jaws/ref_data:/refdata docker://$3 $4 $5
