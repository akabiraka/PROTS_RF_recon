import sys
sys.path.append("../PROTS_RF_recon")

import os
import pandas as pd
import numpy as np
from scipy.special import softmax
from Bio.Blast.Applications import NcbipsiblastCommandline

class PSSM(object):
    """PSSM stands for position-specific socring-matrix.
    https://www.cs.rice.edu/~ogilvie/comp571/2018/09/11/pssm.html
    This url describes first 20 columns gives the log-odds value of each amino-acid at each position.
    The next 20 colums gives the frequency matrix of each amino-acid at each position. 
    """
    def __init__(self, db=None, output_dir=None) -> None:
        super().__init__()
        # self.db = "3rd_party_items/swissprot_db/swissprot" if db is None else db
        self.db = "3rd_party_items/rp_seqs_55/rp_seqs_55"
        self.output_dir = "data/pssms/" if output_dir is None else output_dir
        self.psiblast_exe = "3rd_party_items/ncbi-blast-2.12.0+/bin/psiblast"
        self.pdb_id = None

        
    def set_up(self, fasta_file, force=False):
        """This blast run the query sequence 3 iterations against a db using psiblast program,
        and save the output file in pssms directory. 

        Args:
            fasta_file (str): file path
            force (bool): whether to enforce PSSM set up from start
        """
        pdbid = fasta_file.split("/")[2].split(".")[0]
        output_file_path = self.output_dir + pdbid +".pssm"
        self.pdb_id = pdbid 

        if os.path.exists(output_file_path) and force==False: 
            print("PSSM is already set up for {}. To set-up again, set force=True.".format(pdbid))
            return
        else:
            print("Computing PSSM for {} using psi-blast ... ...".format(pdbid))    
            E_VALUE_TRESH = 10
            cline = NcbipsiblastCommandline(cmd=self.psiblast_exe, db=self.db, query=fasta_file,\
                                            evalue=E_VALUE_TRESH, outfmt=5, num_iterations=3,\
                                            save_pssm_after_last_round=True, out_ascii_pssm=output_file_path)# 
                                            # out = out_xml, out_pssm=out_pssm, out_ascii_pssm=output_file_path)
            cline()
        

    def __get_pssm_file(self):
        return self.output_dir + self.pdb_id + ".pssm"
    

    def __parse_pssm_output_file(self):
        """Parse PSSM raw file generated by PSI-BLAST.

        Args:
            pssm_file (str): a pssm file path

        Returns:
            dataframe: raw result as dataframe
        """
        pssm_file = self.__get_pssm_file()
        col_names = pd.read_csv(pssm_file, delim_whitespace=True, header=None, skiprows=[0, 1], nrows=1)
        col_names = col_names.loc[0, 0:19].tolist()
        residue_dict = {key: i for i, key in enumerate(col_names)}

        df = pd.read_csv(pssm_file, delim_whitespace=True, header=None, skiprows=[0, 1, 2]) # skipping 1st 2 rows
        df = df.head(-5) # removing last 5 rows
        
        return df, residue_dict

   
    def get_logodds(self, i):
        """First 20 columns gives the log-odds value of each amino-acid at each position.
        """
        df, residue_dict = self.__parse_pssm_output_file()
        # print(df.head())
        log_odds = np.array(df.loc[i, 2:21], dtype=np.float32)
        # print(log_odds)
        value = log_odds[residue_dict[df.loc[i, 1]]]
        return np.array([value])


    def get_avg_logodds(self, i, neighbors):
        df, residue_dict = self.__parse_pssm_output_file()
        seq_len = df.shape[0]
        # print(seq_len, df.tail())

        x=0
        if i-int(neighbors/2) < 0:
            for j in range(0, neighbors):
                # print(i+j)
                x=x+self.get_logodds(i+j)

        elif i+int(neighbors/2) > (seq_len-1):
            for j in range(seq_len-neighbors, seq_len):
                # print(j)
                x=x+self.get_logodds(j)

        else: 
            for j in range(-int(neighbors/2), int(neighbors/2)+1):
                # print(i+j)
                x=x+self.get_logodds(i+j)

        return x/neighbors

    def get_softmax(self, i):
        """Probability of i-th residue at i-th position.
        """
        df, residue_dict = self.__parse_pssm_output_file()
        # print(df.head())
        freq_mat = np.array(df.loc[i, 22:41], dtype=np.float32)
        probs = softmax(freq_mat)
        # print("should be 1: ", freq_mat, probs, probs.sum())
        value = probs[residue_dict[df.loc[i, 1]]]
        # print(value)
        return np.array([value])
        

    def get_avg_softmax(self, i, neighbors):
        df, residue_dict = self.__parse_pssm_output_file()
        seq_len = df.shape[0]
        # print(seq_len, df.tail())

        x=0
        if i-int(neighbors/2) < 0:
            for j in range(0, neighbors):
                # print(i+j)
                x=x+self.get_softmax(i+j)

        elif i+int(neighbors/2) > (seq_len-1):
            for j in range(seq_len-neighbors, seq_len):
                # print(j)
                x=x+self.get_softmax(j)

        else: 
            for j in range(-int(neighbors/2), int(neighbors/2)+1):
                # print(i+j)
                x=x+self.get_softmax(i+j)

        return x/neighbors
    
    
   
# # sample usage
# pssm = PSSM()

# variant features at mutation site
# fasta_file = "data/fastas/1a5eA_D_74_N.fasta"
# pssm.set_up(fasta_file)
# print(pssm.get_logodds(73))
# print(pssm.get_softmax(73))

# print(pssm.get_avg_softmax(73, 5))
# print(pssm.get_avg_softmax(73, 9))
# print(pssm.get_avg_softmax(73, 15))

# print(pssm.get_avg_logodds(73, 5))
# print(pssm.get_avg_logodds(73, 9))
# print(pssm.get_avg_logodds(73, 15))  

# wildtype features at mutation site
# fasta_file = "data/fastas/1a5eA.fasta"
# pssm.set_up(fasta_file)
# print(pssm.get_logodds(73))
# print(pssm.get_softmax(73))

# print(pssm.get_avg_softmax(73, 5))
# print(pssm.get_avg_softmax(73, 9))
# print(pssm.get_avg_softmax(73, 15))

# print(pssm.get_avg_logodds(73, 5))
# print(pssm.get_avg_logodds(73, 9))
# print(pssm.get_avg_logodds(73, 15))