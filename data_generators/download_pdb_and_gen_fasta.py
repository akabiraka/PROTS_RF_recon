import sys
sys.path.append("../PROTS_RF_recon")

import csv
import pandas as pd
import numpy as np
import torch
from objects.PDBData import PDBData
from objects.Selector import ChainAndAminoAcidSelect
import data_generators.utils as Utils

# configurations
pdb_dir = "data/pdbs/"
pdbs_clean_dir = "data/pdbs_clean/"
fastas_dir = "data/fastas/"
CIF = "mmCif"

# input_file_path = "data/dataset_5_train.csv"
input_file_path = "data/dataset_5_test.csv"
n_rows_to_skip = 0
n_rows_to_evalutate = 100000

# object initialization
pdbdata = PDBData(pdb_dir=pdb_dir)

df = pd.read_csv(input_file_path)

for i, row in df.iterrows():
    if i+1 <= n_rows_to_skip: continue

    # extracting the data
    pdb_id, chain_id, mutation, mutation_site, wild_residue, mutant_residue, ddg = Utils.get_row_items(row)
    
    # creating necessary file paths
    cln_pdb_file = pdbs_clean_dir+pdb_id+chain_id+".pdb"
    wild_fasta_file = fastas_dir+pdb_id+chain_id+".fasta"
    mutant_fasta_file = fastas_dir+pdb_id+chain_id+"_"+mutation+".fasta"
    
    # downloading and cleaning PDB structure
    pdbdata.download_structure(pdb_id=pdb_id)
    pdbdata.clean(pdb_id=pdb_id, chain_id=chain_id, selector=ChainAndAminoAcidSelect(chain_id))

    # getting 0-based mutation site
    zero_based_mutation_site = Utils.get_zero_based_mutation_site(cln_pdb_file, chain_id, mutation_site)
    print("Row no:{}->{}{}, mutation:{}, 0-based_mutaiton_site:{}".format(i+1, pdb_id, chain_id, mutation, zero_based_mutation_site))

    # generating wild and mutant fasta file
    pdbdata.generate_fasta_from_pdb(pdb_id=pdb_id, chain_id=chain_id, input_pdb_filepath=cln_pdb_file, save_as_fasta=True, output_fasta_dir=fastas_dir)
    pdbdata.create_mutant_fasta_file(wild_fasta_file=wild_fasta_file, mutant_fasta_file=mutant_fasta_file, zero_based_mutation_site=zero_based_mutation_site, mutant_residue=mutant_residue, mutation=mutation)
    
    print()
    if i+1 == n_rows_to_skip+n_rows_to_evalutate: 
        break

