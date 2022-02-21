import sys
sys.path.append("../PROTS_RF_recon")

import pandas as pd
import subprocess
import os
from objects.PDBData import PDBData
from objects.PSSM import PSSM
import data_generators.utils as Utils


# configurations
pdb_dir = "data/pdbs/"
pdbs_clean_dir = "data/pdbs_clean/"
fastas_dir = "data/fastas/"
CIF = "mmCif"
input_file_path = "data/dataset_5_train.csv"
# input_file_path = "data/dataset_5_test.csv"
     
# object initialization
PDBData = PDBData(pdb_dir=pdb_dir)
pssm = PSSM()

# data generation
dfs = pd.read_csv(input_file_path)

def run_command(command):
   return subprocess.getoutput(command)

def generate_wild_pssm(wild_fasta_file):
    pssm.set_up(wild_fasta_file)


def generate_mutant_pssm(mutant_fasta_file):
    pssm.set_up(mutant_fasta_file)


def generate_pssm_for_ith_wild_fasta(i):
    pdb_chain_ids = dfs["pdb_id"]+dfs["chain_id"]
    unique_pdb_chain_ids = pdb_chain_ids.drop_duplicates().to_list()
    unique_pdb_chain_ids.sort()
    print(len(unique_pdb_chain_ids))
    pdb_chain_id = unique_pdb_chain_ids[i]
    print(pdb_chain_id)
    wild_fasta_file = fastas_dir+pdb_chain_id+".fasta"
    generate_wild_pssm(wild_fasta_file)


def generate_pssm_for_ith_mutant_fasta(i):
    pdb_chain_mutation_ids = dfs["pdb_id"]+dfs["chain_id"]+ "_" + dfs["mutation"]
    pdb_chain_mutation_ids = pdb_chain_mutation_ids.drop_duplicates().to_list()
    pdb_chain_mutation_ids.sort()
    print(len(pdb_chain_mutation_ids))
    pdb_chain_mutation_id = pdb_chain_mutation_ids[i]
    print(pdb_chain_mutation_id)
    mutant_fasta_file = fastas_dir+pdb_chain_mutation_id+".fasta"
    generate_mutant_pssm(mutant_fasta_file)


def remove_ith_proteins_pssms(i):
    unique_pdb_ids = dfs["pdb_id"].drop_duplicates().to_list()
    unique_pdb_ids.sort()
    ith_pdb_id = unique_pdb_ids[i]
    ith_protein_mutation_dfs = dfs[dfs["pdb_id"]==ith_pdb_id]
    ith_protein_mutation_dfs = ith_protein_mutation_dfs.reset_index()
    print(ith_pdb_id, ith_protein_mutation_dfs.shape)
    print(ith_protein_mutation_dfs)
    print("removing all pssms for {} ...".format(ith_pdb_id))
    run_command("rm -rf data/pssms/{}*".format(ith_pdb_id))

# i = int(os.environ["SLURM_ARRAY_TASK_ID"]) 
i=11
# remove_ith_proteins_pssms(i)
# generate_pssm_for_ith_wild_fasta(i)
generate_pssm_for_ith_mutant_fasta(i)