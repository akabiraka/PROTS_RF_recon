import sys
sys.path.append("../PROTS_RF_recon")

import pandas as pd
import os
from objects.PDBData import PDBData
from features.PSSM import PSSM


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

    
i = int(os.environ["SLURM_ARRAY_TASK_ID"])    
# i=161
unique_pdb_ids = dfs["pdb_id"].drop_duplicates().to_list()
unique_pdb_ids.sort()
ith_pdb_id = unique_pdb_ids[i]
ith_protein_mutation_dfs = dfs[dfs["pdb_id"]==ith_pdb_id]
# print(ith_protein_mutation_dfs)

for i, row in ith_protein_mutation_dfs.iterrows():
    pdb_id, chain_id, mutation, mutation_site, wild_residue, mutant_residue = row["pdb_id"], row["chain_id"], row["mutation"], int(row["mutation_site"]), row["wild_residue"], row["mutant_residue"]
   
    wild_fasta_file = fastas_dir+pdb_id+chain_id+".fasta"
    mutant_fasta_file = fastas_dir+pdb_id+chain_id+"_"+mutation+".fasta"

    pssm.set_up(wild_fasta_file)
    pssm.set_up(mutant_fasta_file)

    with open("outputs/argo_logs/pssms_done.txt", "a") as f:
        f.write("{}, {}, {}\n".format(i, wild_fasta_file, mutant_fasta_file))
