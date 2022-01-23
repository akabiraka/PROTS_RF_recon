import sys
sys.path.append("../PROTS_RF_recon")

import pandas as pd
import numpy as np

from objects.RelativeDifference import RelativeDifference
import data_generators.utils as Utils

# configs
pdbs_cln_dir="data/pdbs_clean/"
fastas_dir="data/fastas/"
input_file_path="data/dataset_5_train.csv"
# input_file_path="data/dataset_5_test.csv"
out_dir="data/features/relative_diff/"
n_rows_to_skip = 0
n_rows_to_evalutate = 1000000

# objects 
relative_diff=RelativeDifference()
df = pd.read_csv(input_file_path)

for i, row in df.iterrows():
    if i+1 <= n_rows_to_skip: continue
    
    # extracting the data
    pdb_id, chain_id, mutation, mutation_site, wild_residue, mutant_residue, ddg = Utils.get_row_items(row)
    print("Row no:{}->{}{}, mutation:{}".format(i+1, pdb_id, chain_id, mutation))

    wild_fasta_file = fastas_dir+pdb_id+chain_id+".fasta"
    _, wt_posi = relative_diff.count_positive_charged_residues(wild_fasta_file)
    _, wt_nega = relative_diff.count_negative_charged_residues(wild_fasta_file)
    _, wt_char = relative_diff.count_charged_residues(wild_fasta_file)
    _, wt_smal = relative_diff.count_small_residues(wild_fasta_file)
    _, wt_tiny = relative_diff.count_tiny_residues(wild_fasta_file)

    mutant_fasta_file = fastas_dir+pdb_id+chain_id+"_"+mutation+".fasta"
    _, mt_posi = relative_diff.count_positive_charged_residues(mutant_fasta_file)
    _, mt_nega = relative_diff.count_negative_charged_residues(mutant_fasta_file)
    _, mt_char = relative_diff.count_charged_residues(mutant_fasta_file)
    _, mt_smal = relative_diff.count_small_residues(mutant_fasta_file)
    _, mt_tiny = relative_diff.count_tiny_residues(mutant_fasta_file)

    posi = wt_posi - mt_posi
    nega = wt_nega - mt_nega
    char = wt_char - mt_char
    smal = wt_smal - mt_smal
    tiny = wt_tiny - mt_tiny


    features = np.array([posi, nega, char, smal, tiny])
    Utils.save_as_pickle(features, out_dir+pdb_id+chain_id+"_"+mutation+".pkl")
    # features = Utils.load_pickle(out_dir+pdb_id+chain_id+"_"+mutation+".pkl")
    print("saved features of shape: {}, {}".format(features.shape, features))

    print()
    if i+1 == n_rows_to_skip+n_rows_to_evalutate: break