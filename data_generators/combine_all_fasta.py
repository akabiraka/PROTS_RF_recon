import sys
sys.path.append("../PROTS_RF_recon")

import pandas as pd
import data_generators.utils as Utils

fastas_dir = "data/fastas/"
input_file_path = "data/dataset_5_train.csv"
# input_file_path = "data/dataset_5_test.csv"
n_rows_to_skip = 0
n_rows_to_evalutate = 100#00

df = pd.read_csv(input_file_path) 


out_file = open("all.fasta", "w")



def combine_wild():
    unique_pdb_ids = (df["pdb_id"]+df["chain_id"]).drop_duplicates().to_list()
    for pdb_id in unique_pdb_ids:
        wild_fasta_file = fastas_dir+pdb_id+".fasta"
        f = open(wild_fasta_file, "r") 
        for line in f:
            out_file.write(line)
        out_file.write("\n")



def combine_wild_and_variants():
    unique_pdb_ids = set()
    for i, row in df.iterrows():
        if i+1 <= n_rows_to_skip: continue

        pdb_id, chain_id, mutation, _, _, _, _ = Utils.get_row_items(row)
        print("Row no:{}->{}{}, mutation:{}".format(i+1, pdb_id, chain_id, mutation))

        # creating necessary file paths
        wild_fasta_file = fastas_dir+pdb_id+chain_id+".fasta"
        mutant_fasta_file = fastas_dir+pdb_id+chain_id+"_"+mutation+".fasta"   

        if pdb_id+chain_id not in unique_pdb_ids:
            f = open(wild_fasta_file, "r") 
            for line in f:
                out_file.write(line)
            out_file.write("\n")
            unique_pdb_ids.add(pdb_id+chain_id)
        
        f = open(mutant_fasta_file, "r") 
        for line in f:
            out_file.write(line)
        out_file.write("\n")

        if i+1 == n_rows_to_skip+n_rows_to_evalutate: 
            break


combine_wild()
# combine_wild_and_variants()


