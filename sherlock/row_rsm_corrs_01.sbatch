#!/bin/bash
#
#set the job name (output file and error file)
#SBATCH --job-name=row_01_rh
#SBATCH --output=row_01_rh.out
#SBATCH --error=row_01_rh.error
#################
#set time needed
#SBATCH --time=24:00:00
#################
#SBATCH --nodes=1
#################
#request 256GB RAM on CPU
#SBATCH --mem=1024000

module load python/3.6.1
source ${OAK}/biac2/kgs/projects/Dawn/sherlock_venv/bin/activate 
srun python ${OAK}/biac2/kgs/projects/Dawn/NSD/code/streams/scripts/make_rsm_corrs_by_row.py --subjid "01" --hemi 'rh' --roi_name 'tessellate_1000'

