#!/bin/bash
#
#set the job name (output file and error file)
#SBATCH --job-name=rsm_corrs_06_rh
#SBATCH --output=rsm_corrs_06_rh.out
#SBATCH --error=rsm_corrs_06_rh.error
#################
#set time needed
#SBATCH --time=12:00:00
#################
#SBATCH --nodes=1
#################
#request 256GB RAM on CPU
#SBATCH --mem=256000

module load python/3.6.1
source ${OAK}/biac2/kgs/projects/Dawn/sherlock_venv/bin/activate 
srun python ${OAK}/biac2/kgs/projects/Dawn/NSD/code/streams/scripts/make_rsm_corrs.py --subjid "06" --hemi 'rh' --roi_name 'tessellate_300' --thresh 0.0

