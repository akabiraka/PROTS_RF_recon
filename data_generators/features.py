import sys
sys.path.append("../PROTS_RF_recon")

import numpy as np
import pandas as pd
import data_generators.utils as Utils

input_file_path = "data/dataset_5_train.csv"
# input_file_path = "data/dataset_5_test.csv"
n_rows_to_skip = 0
n_rows_to_evalutate = 1#00000

df = pd.read_csv(input_file_path)

for i, row in df.iterrows():
    if i+1 <= n_rows_to_skip: continue

    # extracting the data
    # pdb_id, chain_id, mutation, mutation_site, wild_residue, mutant_residue, ddg = Utils.get_row_items(row)
    pdb_id, chain_id, mutation = "1a0f", "A", "S_11_A"
    print("Row no:{}->{}{}, mutation:{}".format(i+1, pdb_id, chain_id, mutation))

    ss_sa_rasa_features = Utils.load_pickle("data/features/ss_sa_rasa/"+pdb_id+chain_id+".pkl")
    relative_diff_features = Utils.load_pickle("data/features/relative_diff/"+pdb_id+chain_id+"_"+mutation+".pkl")
    # evolutionary_features = Utils.load_pickle("data/features/evolutionary/"+pdb_id+chain_id+"_"+mutation+".pkl")

    features = np.concatenate((ss_sa_rasa_features, relative_diff_features))
    print(ss_sa_rasa_features.shape, relative_diff_features.shape, features.shape)

    print()
    if i+1 == n_rows_to_skip+n_rows_to_evalutate: 
        break