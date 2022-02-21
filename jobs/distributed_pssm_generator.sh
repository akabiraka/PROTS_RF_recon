#!/usr/bin/sh

## this must be run from directory where run.py exists.
## --workdir is not used in this file.

#SBATCH --job-name=wild
#SBATCH --output=/scratch/akabir4/PROTS_RF_recon/outputs/argo_logs/distributed_pssm_generator-%N-%j.output
#SBATCH --error=/scratch/akabir4/PROTS_RF_recon/outputs/argo_logs/distributed_pssm_generator-%N-%j.error
#SBATCH --mail-user=<akabir4@gmu.edu>
#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --partition=all-LoPri
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000MB

#SBATCH --array=0-36
python data_generators/distributed_pssm.py
