import sys
sys.path.append("../PROTS_RF_recon")

from objects.MutationUtils import MutationUtils
from objects.PDBData import PDBData
from itertools import permutations
from Bio import SeqIO
from Bio.PDB import PDBParser, PDBIO, DSSP
import os
import re

class PROTS(object):
    def __init__(self) -> None:
        super().__init__()
        self.pdbdata = PDBData()
        self.mutation_utils = MutationUtils()
        self.pdbio = PDBIO()
        self.pdbparser = PDBParser(QUIET=True)
        self.fastas_dir = "data/fastas/"
        self.pdbs_cln_dir = "data/pdbs_clean/"
        self.out_file = "data/features/prots/xxx.csv"
        self.ss_dict={"H":"H", "G":"H", "I":"H", "B":"B", "E":"B", "T":"C", "S":"C", "-":"C"}


    def do_permutation(self, fragment_size=4, force=False):
        col_name = "fragments"
        if os.path.exists(self.out_file) and force==False:   
            print("Permutations already computed. To compute again use force=True")
            return (pd.read_csv(self.out_file))[col_name].to_list() 
        else:
            print("Computing permutations ...") 
            aminos_acids = "ARNDCEQGHILKMFPSTWYV"
            fragments=["".join(i) for i in permutations(aminos_acids, fragment_size)]
            df = pd.DataFrame(fragments, columns=[col_name])
            df.to_csv(self.out_file, index=False)
            return fragments
            

    
    def get_occurance_of_fragment(self, fragment, df):
        sum=0
        for i, row in df.iterrows():
            pdb_id, chain_id, mutation, mutation_site, wild_residue, mutant_residue, ddg = self.mutation_utils.get_row_items(row)
            fasta_file = self.fastas_dir+pdb_id+chain_id+"_"+mutation+".fasta"
            fasta = next(SeqIO.parse(fasta_file, "fasta"))
            sum = sum + fasta.seq.count(fragment) # the num of time a fragment occurs in a seq
        return sum

    
    def compute_occurance(self, fragments, df):
        for i, fragment in enumerate(fragments):
            print(i, fragment, self.get_occurance_of_fragment(fragment, df))
        return 0

    
    def get_sa_type(self, rasa):
        if rasa<0.25: sa="B" # Buried
        elif rasa>0.5: sa="E" # Exposed
        else: sa="I" # Intermediate 
        return sa

    def get_frag_occ_in_ss_sa(self, fragment, fasta_seq, ss, sa):
        n_helix, n_sheet, n_coil, n_buried, n_exposed, n_intermediate=0,0,0,0,0,0
        indices = [m.start() for m in re.finditer(fragment, str(fasta_seq))]
        print("fragment {} occurred in {}".format(fragment, indices))
        for i in indices:
            if "H" in ss[i:i+len(fragment)]: n_helix=1
            if "B" in ss[i:i+len(fragment)]: n_sheet=1
            if "C" in ss[i:i+len(fragment)]: n_coil=1

            if "B" in sa[i:i+len(fragment)]: n_buried=1
            if "E" in sa[i:i+len(fragment)]: n_exposed=1
            if "I" in sa[i:i+len(fragment)]: n_intermediate=1
        return n_helix, n_sheet, n_coil, n_buried, n_exposed, n_intermediate

    def compute_feature_occurance(self, pdb_chain_ids, fragment="PDGH"):
        n_helix, n_sheet, n_coil, n_buried, n_exposed, n_intermediate=0,0,0,0,0,0
        for pdb_chain_id in pdb_chain_ids:
            pdb_id = pdb_chain_id[:4]
            chain_id = pdb_chain_id[4]
            print(pdb_chain_id)

            fasta_file = self.fastas_dir+pdb_chain_id+".fasta"
            fasta = next(SeqIO.parse(fasta_file, "fasta"))

            cln_pdb_file = self.pdbs_cln_dir+pdb_chain_id+".pdb"
            ss, sa, _ = self.pdbdata.get_full_ss_and_sa(pdb_id, chain_id, cln_pdb_file, self.ss_dict, self.get_sa_type)
            # print(ss)
            # print(sa)
            # print(fasta.seq)
            helix, beta, coil, buried, exposed, intr = self.get_frag_occ_in_ss_sa(fragment, str(fasta.seq), ss, sa)
            n_helix+=helix; n_sheet+=beta; n_coil+=coil; n_buried+=buried; n_exposed+=exposed; n_intermediate+=intr
            
            # break
            print(n_helix, n_sheet, n_coil, n_buried, n_exposed, n_intermediate, "\n")
        return n_helix, n_sheet, n_coil, n_buried, n_exposed, n_intermediate

    def compute_feature_feagment_occurance(self):
        pass
        

import pandas as pd
df = pd.read_csv("data/dataset_5_train.csv")
# stabilizing_df, destabilizing_df = df[df["ddG"]>=0.0], df[df["ddG"]<0.0]
# print(df.shape, stabilizing_df.shape, destabilizing_df.shape)
pdb_chain_ids = (df["pdb_id"]+df["chain_id"]).unique().tolist()
# print(len(pdb_chain_ids))

prots = PROTS()
# fragments = prots.do_permutation(4, False)
# print(len(fragments))


prots.compute_feature_occurance(pdb_chain_ids)
# n_occurance_st = prots.compute_occurance(fragments, stabilizing_df)
# n_occurance_de = prots.compute_occurance(fragments, destabilizing_df)

# prots.compute_feature_occurance(stabilizing_df, feature="helix", fragment="ARND")