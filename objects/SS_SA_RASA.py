import sys
sys.path.append("../PROTS_RF_recon")

import numpy as np
from objects.PDBData import PDBData

class SS_SA_RASA(object):
    def __init__(self) -> None:
        super().__init__()
        self.pdbdata=PDBData()

        # helix: H,G,I; sheet:B,E; coil:T,S,-
        self.ss_dict={"H":"H", "G":"H", "I":"H", "B":"B", "E":"B", "T":"C", "S":"C", "-":"C"}

    def get_ss_onehot(self, ss):
        letter = self.ss_dict.get(ss)
        return np.array([0.1 if char != letter else 0.9 for char in "HBC"], dtype=np.float32)

    def get_sa_onehot(self, rasa):
        letter="E" if rasa>=0.25 else "B" #0=exposed, 0=buried
        return np.array([0.1 if char != letter else 0.9 for char in "EB"], dtype=np.float32)

    def get_ss_sa_rasa(self, pdb_id, chain_id, cln_pdb_file, mutation_site):
        ss, rasa = self.pdbdata.get_ss_and_rasa_at_residue(pdb_id, chain_id, cln_pdb_file, mutation_site)
        ss_onehot, sa_onehot=self.get_ss_onehot(ss), self.get_sa_onehot(rasa)
        return ss_onehot, sa_onehot, rasa


    