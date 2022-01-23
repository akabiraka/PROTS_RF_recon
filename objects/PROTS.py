import sys
sys.path.append("../PROTS_RF_recon")

import pandas as pd
import data_generators.utils as Utils
from objects.PDBData import PDBData
from itertools import permutations
from Bio import SeqIO
from Bio.PDB import PDBParser, PDBIO, DSSP
import os
import re
import json

class PROTS(object):
    def __init__(self) -> None:
        super().__init__()

        self.pdbdata = PDBData()
        self.pdbio = PDBIO()
        self.pdbparser = PDBParser(QUIET=True)
        
        self.fastas_dir = "data/fastas/"
        self.pdbs_cln_dir = "data/pdbs_clean/"
        
        self.ss_dict={"H":"H", "G":"H", "I":"H", "B":"B", "E":"B", "T":"C", "S":"C", "-":"C"}
        
        self.out_file = "data/features/prots/fragment_vs_feature_occ_count_1.csv"
        self.json_out_file = "data/features/prots/fragment_vs_feature_occ_indices_per_protein.json"
        if not os.path.exists(self.out_file): 
            pd.DataFrame(columns=["fragment", "n_occurance", "n_helix", "n_sheet", "n_coil", "n_buried", "n_exposed", "n_intermediate"]).to_csv(self.out_file, index=False)
        


    def do_permutation(self, fragment_size=4, force=False):
        col_name = "fragment"
        if os.path.exists(self.out_file) and force==False:   
            print("Permutations already computed. To compute again use force=True")
            return (pd.read_csv(self.out_file))[col_name].to_list() 
        else:
            print("Computing permutations ...") 
            aminos_acids = "ARNDCEQGHILKMFPSTWYV"
            fragments=["".join(i) for i in permutations(aminos_acids, fragment_size)]
            out_df = pd.read_csv(self.out_file)
            out_df[col_name] = fragments
            out_df.to_csv(self.out_file, index=False)
            return fragments
            

    
    def get_occurance_of_fragment(self, fragment):
        df = pd.read_csv("data/dataset_5_train.csv")
        sum=0
        for i, row in df.iterrows():
            pdb_id, chain_id, mutation, _, _, _, _ = Utils.get_row_items(row)
            fasta_file = self.fastas_dir+pdb_id+chain_id+"_"+mutation+".fasta"
            fasta = next(SeqIO.parse(fasta_file, "fasta"))
            sum = sum + fasta.seq.count(fragment) # the num of time a fragment occurs in a seq
        return sum


    
    def get_sa_type(self, rasa):
        if rasa<0.25: sa="B" # Buried
        elif rasa>0.5: sa="E" # Exposed
        else: sa="I" # Intermediate 
        return sa

    def get_frag_occ_in_ss_sa(self, fragment, fasta_seq, ss, sa):
        n_helix, n_sheet, n_coil, n_buried, n_exposed, n_intermediate=0,0,0,0,0,0
        occ_indices = [m.start() for m in re.finditer(fragment, str(fasta_seq))]
        
        for i in occ_indices:
            if "H" in ss[i:i+len(fragment)]: n_helix=1
            if "B" in ss[i:i+len(fragment)]: n_sheet=1
            if "C" in ss[i:i+len(fragment)]: n_coil=1

            if "B" in sa[i:i+len(fragment)]: n_buried=1
            if "E" in sa[i:i+len(fragment)]: n_exposed=1
            if "I" in sa[i:i+len(fragment)]: n_intermediate=1
        return n_helix, n_sheet, n_coil, n_buried, n_exposed, n_intermediate, occ_indices

    def get_all_feat_occ_for_a_frag(self, fragment):
        df = pd.read_csv("data/dataset_5_train.csv")
        pdb_chain_ids = (df["pdb_id"]+df["chain_id"]).unique().tolist()

        n_helix, n_sheet, n_coil, n_buried, n_exposed, n_intermediate=0,0,0,0,0,0

        all_occurances = []
        for pdb_chain_id in pdb_chain_ids:
            pdb_id = pdb_chain_id[:4] # "1h7m" 
            chain_id = pdb_chain_id[4] # "A" 
            # pdb_chain_id="1h7mA"

            fasta_file = self.fastas_dir+pdb_chain_id+".fasta"
            fasta = next(SeqIO.parse(fasta_file, "fasta"))

            cln_pdb_file = self.pdbs_cln_dir+pdb_chain_id+".pdb"
            ss, sa, _ = self.pdbdata.get_full_ss_and_sa(pdb_id, chain_id, cln_pdb_file, self.ss_dict, self.get_sa_type)
            # print(ss)
            # print(sa)
            # print(fasta.seq)
            helix, beta, coil, buried, exposed, intr, occ_indices_in_a_seq = self.get_frag_occ_in_ss_sa(fragment, str(fasta.seq), ss, sa)
            n_helix+=helix; n_sheet+=beta; n_coil+=coil; n_buried+=buried; n_exposed+=exposed; n_intermediate+=intr

            print("    fragment {} occurred in {} for features: {}".format(fragment, pdb_chain_id, occ_indices_in_a_seq))
            all_occurances.append({pdb_chain_id: occ_indices_in_a_seq})
            
            # break
        return n_helix, n_sheet, n_coil, n_buried, n_exposed, n_intermediate, all_occurances

    def compute_all_feat_occ_for_all_frag(self, row_index):
        out_df = pd.read_csv(self.out_file)
        print(out_df.loc[row_index])

        all_frag_occ_indices = []
        # for i, row in out_df.iterrows():

        fragment = out_df.loc[row_index, "fragment"]
        print("Computing occ for: ", row_index, fragment)

        n_occurance = self.get_occurance_of_fragment(fragment)
        print("    fragment {} occurred in the dataset: {}".format(fragment, n_occurance))
        
        n_helix, n_sheet, n_coil, n_buried, n_exposed, n_intermediate, all_occurances = self.get_all_feat_occ_for_a_frag(fragment=fragment)
        print("    feature occ: ", n_helix, n_sheet, n_coil, n_buried, n_exposed, n_intermediate)
        
        
        out_df.loc[row_index, "n_occurance"] = n_occurance
        out_df.loc[row_index, "n_helix"] = n_helix
        out_df.loc[row_index, "n_sheet"] = n_sheet
        out_df.loc[row_index, "n_coil"] = n_coil
        out_df.loc[row_index, "n_buried"] = n_buried
        out_df.loc[row_index, "n_exposed"] = n_exposed
        out_df.loc[row_index, "n_intermediate"] = n_intermediate


        all_frag_occ_indices.append({fragment: all_occurances})
        # break

        print("    updating the occ and json file ....\n")
        out_df.to_csv(self.out_file, index=False)

        with open(self.json_out_file, "a") as f:
            json_out = json.dumps({fragment: all_occurances})
            f.write(json_out)
            f.write(",\n")
