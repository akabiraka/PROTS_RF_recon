#!/usr/bin/sh

#SBATCH --job-name=prots_feature
#SBATCH --output=/scratch/akabir4/PROTS_RF_recon/outputs/argo_logs/prots_feature-%N-%j.output
#SBATCH --error=/scratch/akabir4/PROTS_RF_recon/outputs/argo_logs/prots_feature-%N-%j.error
#SBATCH --mail-user=<akabir4@gmu.edu>
#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --partition=all-LoPri
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000MB

#SBATCH --array=0-116279
##the array task is set in the environment variable $SLURM_ARRAY_TASK_ID in python you
##can scrape it with ID = int(os.environ["SLURM_ARRAY_TASK_ID"])

python data_generators/feature_prots.py
