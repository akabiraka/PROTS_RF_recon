import sys
sys.path.append("../PROTS_RF_recon")

import numpy as np
import pandas as pd

from features.SS_SA_RASA import SS_SA_RASA
from objects.MutationUtils import MutationUtils

# configs
pdbs_cln_dir="data/pdbs_clean/"
input_file_path="data/dataset_5_train.csv"
# input_file_path="data/dataset_5_test.csv"
out_dir="data/features/ss_sa_rasa/"
n_rows_to_skip = 0
n_rows_to_evalutate = 10

# objects 
ss_sa_rasa=SS_SA_RASA()
mutation_utils=MutationUtils()
df = pd.read_csv(input_file_path)


for i, row in df.iterrows():
    if i+1 <= n_rows_to_skip: continue
    
    # extracting the data
    pdb_id, chain_id, mutation, mutation_site, wild_residue, mutant_residue, ddg = mutation_utils.get_row_items(row)
    
    cln_pdb_file = pdbs_cln_dir+pdb_id+chain_id+".pdb"
    
    ss_onehot, sa_onehot, rasa = ss_sa_rasa.get_ss_sa_rasa(pdb_id, chain_id, cln_pdb_file, mutation_site)

    print("saving ss, sa and rasa: {}{}, mutation:{}, mutation_site:{}".format(pdb_id, chain_id, mutation, mutation_site))
    features=np.concatenate((ss_onehot, sa_onehot, [rasa]))
    with open(out_dir+pdb_id+chain_id+".npy", 'wb') as f: np.save(f, features)
    with open(out_dir+pdb_id+chain_id+".npy", 'rb') as f: print(np.load(f))

    print()
    if i+1 == n_rows_to_skip+n_rows_to_evalutate: 
        break
