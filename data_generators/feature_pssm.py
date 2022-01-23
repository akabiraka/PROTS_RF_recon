import sys
sys.path.append("../PROTS_RF_recon")

import pandas as pd
import numpy as np

from objects.PSSM import PSSM
import data_generators.utils as Utils

# configs
pdbs_cln_dir="data/pdbs_clean/"
fastas_dir="data/fastas/"
input_file_path="data/dataset_5_train.csv"
# input_file_path="data/dataset_5_test.csv"
out_dir="data/features/pssms/"
n_rows_to_skip = 0
n_rows_to_evalutate = 10

# objects 
pssm=PSSM()
df = pd.read_csv(input_file_path)

for i, row in df.iterrows():
    if i+1 <= n_rows_to_skip: continue
    
    # extracting the data
    pdb_id, chain_id, mutation, mutation_site, wild_residue, mutant_residue, ddg = Utils.get_row_items(row)
    
    clean_pdb_file = pdbs_cln_dir+pdb_id+chain_id+".pdb"
    zero_based_mutation_site = Utils.get_zero_based_mutation_site(clean_pdb_file, chain_id, mutation_site)
    print("Row no:{}->{}{}, mutation:{}, mutation_site:{}, zero_based_mutation_site:{}".format(i+1, pdb_id, chain_id, mutation, mutation_site, zero_based_mutation_site))

    print("computing wildtype features ...")
    wild_fasta_file = fastas_dir+pdb_id+chain_id+".fasta"
    pssm.set_up(wild_fasta_file)
    wtlo = pssm.get_logodds(zero_based_mutation_site) #wildtype logodds
    wtsm = pssm.get_softmax(zero_based_mutation_site) #wildtype softmax

    wtlo5 = pssm.get_avg_logodds(zero_based_mutation_site, 5) #wildtype logodds of 5 neighbors
    wtlo9 = pssm.get_avg_logodds(zero_based_mutation_site, 9) #wildtype logodds of 9 neighbors
    wtlo15 = pssm.get_avg_logodds(zero_based_mutation_site, 15) #wildtype logodds of 15 neighbors

    wtsm5 = pssm.get_avg_softmax(zero_based_mutation_site, 5) #wildtype softmax of 5 neighbors
    wtsm9 = pssm.get_avg_softmax(zero_based_mutation_site, 9) #wildtype softmax of 9 neighbors
    wtsm15 = pssm.get_avg_softmax(zero_based_mutation_site, 15) #wildtype softmax of 15 neighbors

    print("computing variant features ...")
    mutant_fasta_file = fastas_dir+pdb_id+chain_id+"_"+mutation+".fasta"
    pssm.set_up(mutant_fasta_file)
    vlo = pssm.get_logodds(zero_based_mutation_site) #variant logodds
    vsm = pssm.get_softmax(zero_based_mutation_site) #variant softmax

    vlo5 = pssm.get_avg_softmax(zero_based_mutation_site, 5) #variant softmax of 5 neighbors
    vlo9 = pssm.get_avg_softmax(zero_based_mutation_site, 9) #variant softmax of 9 neighbors
    vlo15 = pssm.get_avg_softmax(zero_based_mutation_site, 15) #variant softmax of 15 neighbors

    vsm5 = pssm.get_avg_logodds(zero_based_mutation_site, 5) #variant logodds of 5 neighbors
    vsm9 = pssm.get_avg_logodds(zero_based_mutation_site, 9) #variant logodds of 9 neighbors
    vsm15 = pssm.get_avg_logodds(zero_based_mutation_site, 15) #variant logodds of 15 neighbors

    
    features = np.concatenate((wtlo, wtsm, wtsm5, wtsm9, wtsm15, wtlo5, wtlo9, wtlo15, vlo, vsm, vsm5, vsm9, vsm15, vlo5, vlo9, vlo15))
    with open(out_dir+pdb_id+chain_id+"_"+mutation+".npy", 'wb') as f: np.save(f, features)
    # with open(out_dir+pdb_id+chain_id+"_"+mutation+".npy", 'rb') as f: print(np.load(f))
    print("saved pssm features of shape: {}: {}".format(features.shape, features))

    print()
    if i+1 == n_rows_to_skip+n_rows_to_evalutate: 
        break