#!/bin/bash
docker=$3
#sif_path=/global/scratch/jaws/sif_files
sif_path=/global/scratch/jfroula/JAWS/img-omics-wdl/sif_files
export SINGULARITY_CACHEDIR=${sif_path}
export SINGULARITY_PULLFOLDER=${sif_path}
export SINGULARITY_TMPDIR=${sif_path}
export SINGULARITY_LOCALCACHEDIR=${sif_path}

#DOCKER_NAME=$(sed -e 's/[^A-Za-z0-9._-]/_/g' <<< ${docker})
#IMAGE=${sif_path}/$DOCKER_NAME.sif
#if [ ! -f $IMAGE ]; then
#    singularity pull $IMAGE docker://${docker}
#fi

singularity exec --bind $1:$2 --bind /global/scratch/jaws/ref_data:/refdata docker://$3 $4 $5
#singularity exec --bind $1:$2 --bind /global/scratch/jaws/ref_data:/refdata $IMAGE $4 $5
