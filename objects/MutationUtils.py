import sys
sys.path.append("../PROTS_RF_recon")

from objects.PDBData import PDBData


class MutationUtils(object):
    def __init__(self)->None:
        super().__init__()
        self.pdbdata = PDBData(pdb_dir="data/pdbs_clean/")

    def get_row_items(self, row):
        return row["pdb_id"].lower()[:4], row["chain_id"], row["mutation"], int(row["mutation_site"]), row["wild_residue"], row["mutant_residue"], row["ddG"]

    def get_zero_based_mutation_site(self, clean_pdb_file, chain_id, mutation_site):
        residue_ids_dict = self.pdbdata.get_residue_ids_dict(pdb_file=clean_pdb_file, chain_id=chain_id)
        return residue_ids_dict.get(mutation_site)