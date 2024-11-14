#!/bin/bash
#SBATCH --job-name=fl52_job 
#SBATCH --account=commons 
#SBATCH --partition=commons  
#SBATCH --time=20:00:00 
#SBATCH --mail-user=fl52@rice.edu 
#SBATCH --mail-type=ALL
#SBATCH --array=0-100
#SBATCH --error=/scratch/fl52/log/job.err
#SBATCH --output=/scratch/fl52/log/job.out
echo "My job ran on:" echo $SLURM_NODELIST
if [[ -d $SHARED_SCRATCH/$USER && -w $SHARED_SCRATCH/$USER ]]; then
    cd $SHARED_SCRATCH/$USER
    srun pp_highMultGen $SLURM_ARRAY_TASK_ID	    
fi

