import sys
sys.path.append("../PROTS_RF_recon")

import numpy as np
from Bio import SeqIO

class RelativeDifference(object):
    def __init__(self) -> None:
        super().__init__()


    def __get_cnt_and_normalized(self, fasta_file, residues):
        fasta = next(SeqIO.parse(fasta_file, "fasta"))
        cnt = np.sum([fasta.seq.count(res) for res in residues])
        cnt_normalized = cnt/len(fasta.seq) 
        return cnt, cnt_normalized


    def count_positive_charged_residues(self, fasta_file):
        pos_residues=["R", "K", "H"]
        cnt, cnt_normalized = self.__get_cnt_and_normalized(fasta_file, pos_residues)
        # print(cnt, cnt_normalized)
        return cnt, cnt_normalized

    def count_negative_charged_residues(self, fasta_file):
        neg_residues=["E", "D"]
        cnt, cnt_normalized = self.__get_cnt_and_normalized(fasta_file, neg_residues) 
        # print(cnt, cnt_normalized)
        return cnt, cnt_normalized

    def count_charged_residues(self, fasta_file):
        charged_residues=["R", "K", "H", "E", "D"]
        cnt, cnt_normalized = self.__get_cnt_and_normalized(fasta_file, charged_residues)  
        # print(cnt, cnt_normalized)
        return cnt, cnt_normalized

    def count_small_residues(self, fasta_file):
        sml_residues=["T", "D"]
        cnt, cnt_normalized = self.__get_cnt_and_normalized(fasta_file, sml_residues) 
        # print(cnt, cnt_normalized)
        return cnt, cnt_normalized

    def count_tiny_residues(self, fasta_file):
        tiny_residues=["A", "G", "P", "S"]
        cnt, cnt_normalized = self.__get_cnt_and_normalized(fasta_file, tiny_residues) 
        # print(cnt, cnt_normalized)
        return cnt, cnt_normalized


# relative_diff=RelativeDifference()
# fasta_file="data/fastas/1a5eA.fasta"
# relative_diff.count_positive_charged_residues(fasta_file)
# relative_diff.count_negative_charged_residues(fasta_file)
# relative_diff.count_charged_residues(fasta_file)
# relative_diff.count_small_residues(fasta_file)
# relative_diff.count_tiny_residues(fasta_file)

# fasta_file="data/fastas/1a5eA_D_74_N.fasta"
# relative_diff.count_positive_charged_residues(fasta_file)
# relative_diff.count_negative_charged_residues(fasta_file)
# relative_diff.count_charged_residues(fasta_file)
# relative_diff.count_small_residues(fasta_file)
# relative_diff.count_tiny_residues(fasta_file)
