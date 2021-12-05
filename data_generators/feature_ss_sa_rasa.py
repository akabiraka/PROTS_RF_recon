import sys
sys.path.append("../PROTS_RF_recon")

import numpy as np
import pandas as pd
from objects.PDBData import PDBData
from objects.MutationUtils import MutationUtils

# configs
pdbs_cln_dir="data/pdbs_clean/"
input_file_path="data/dataset_5_train.csv"
# input_file_path="data/dataset_5_test.csv"
out_dir="data/features/ss_sa/"
n_rows_to_skip = 0
n_rows_to_evalutate = 10

# objects 
pdbdata=PDBData(pdb_dir=pdbs_cln_dir)
mutation_utils=MutationUtils()
df = pd.read_csv(input_file_path)

# helix: H,G,I; sheet:B,E; coil:T,S,-
ss_dict={"H":"H", "G":"H", "I":"H", "B":"B", "E":"B", "T":"C", "S":"C", "-":"C"}
def get_ss_onehot(ss):
    letter = ss_dict.get(ss)
    return np.array([0.1 if char != letter else 0.9 for char in "HBC"], dtype=np.float32)

def get_sa_onehot(rasa):
    letter="E" if rasa>=0.25 else "B" #0=exposed, 0=buried
    return np.array([0.1 if char != letter else 0.9 for char in "EB"], dtype=np.float32)



for i, row in df.iterrows():
    if i+1 <= n_rows_to_skip: continue
    
    # extracting the data
    pdb_id, chain_id, mutation, mutation_site, wild_residue, mutant_residue, ddg = mutation_utils.get_row_items(row)
    
    cln_pdb_file = pdbs_cln_dir+pdb_id+chain_id+".pdb"
    
    ss, rasa = pdbdata.get_ss_and_rasa_at_residue(pdb_id, chain_id, cln_pdb_file, mutation_site)

    print("saving ss, sa and rasa: {}{}, mutation:{}, mutation_site:{}".format(pdb_id, chain_id, mutation, mutation_site))
    features=np.concatenate((get_ss_onehot(ss), get_sa_onehot(rasa), [rasa]))
    with open(out_dir+pdb_id+chain_id+".npy", 'wb') as f: np.save(f, features)
    # with open(out_dir+pdb_id+chain_id+".npy", 'rb') as f: print(np.load(f))

    print()
    if i+1 == n_rows_to_skip+n_rows_to_evalutate: 
        break
