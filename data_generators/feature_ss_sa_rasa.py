import sys
sys.path.append("../PROTS_RF_recon")

import numpy as np
import pandas as pd

from objects.SS_SA_RASA import SS_SA_RASA
import data_generators.utils as Utils

# configs
pdbs_cln_dir="data/pdbs_clean/"
input_file_path="data/dataset_5_train.csv"
# input_file_path="data/dataset_5_test.csv"
out_dir="data/features/ss_sa_rasa/"
n_rows_to_skip = 0
n_rows_to_evalutate = 100000

# objects 
ss_sa_rasa=SS_SA_RASA()
df = pd.read_csv(input_file_path)


for i, row in df.iterrows():
    if i+1 <= n_rows_to_skip: continue
    
    # extracting the data
    pdb_id, chain_id, mutation, mutation_site, wild_residue, mutant_residue, ddg = Utils.get_row_items(row)
    cln_pdb_file = pdbs_cln_dir+pdb_id+chain_id+".pdb"
    
    print("Row no:{}->computing ss, sa and rasa: {}{}, mutation:{}, mutation_site:{}".format(i+1, pdb_id, chain_id, mutation, mutation_site))
    ss_onehot, sa_onehot, rasa = ss_sa_rasa.get_ss_sa_rasa(pdb_id, chain_id, cln_pdb_file, mutation_site)

    # saving features
    features=np.concatenate((ss_onehot, sa_onehot, [rasa]))
    Utils.save_as_pickle(features, out_dir+pdb_id+chain_id+".pkl")
    # features = Utils.load_pickle(out_dir+pdb_id+chain_id+".pkl")
    print("saved features of shape: {}, {}".format(features.shape, features))

    print()
    if i+1 == n_rows_to_skip+n_rows_to_evalutate: 
        break
