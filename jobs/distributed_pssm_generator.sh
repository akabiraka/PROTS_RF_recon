#!/usr/bin/sh

## this must be run from directory where run.py exists.
## --workdir is not used in this file.

#SBATCH --job-name=PROTS_RF_recon
#SBATCH --qos=csqos
#SBATCH --output=/scratch/akabir4/PROTS_RF_recon/outputs/argo_logs/distributed_pssm_generator-%N-%j.output
#SBATCH --error=/scratch/akabir4/PROTS_RF_recon/outputs/argo_logs/distributed_pssm_generator-%N-%j.error
#SBATCH --mail-user=<akabir4@gmu.edu>
#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --partition=all-HiPri
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000MB

##for training
#SBATCH --array=0-202
##here total 203 unique proteins, where each may have multiple mutations
##for testing
##SBATCH --array=0-36
##the array task is set in the environment variable $SLURM_ARRAY_TASK_ID in python you
##can scrape it with ID = int(os.environ["SLURM_ARRAY_TASK_ID"])

python data_gen/distributed_pssm_generator.py
