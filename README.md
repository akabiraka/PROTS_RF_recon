# PROTS-RF: reconstruction

Base paper title: PROTS-RF: A Robust Model for Predicting MutationInduced Protein Stability Changes.

#### Feature generation

* **Download and clean**: `python data_generators/download_pdb_and_gen_fasta.py`
* **Evolutionary features** (dimension=16):
  * `sbatch jobs/distributed_pssm_generator.sh`
  * `python data_generators/feature_pssm.py`
* **Secondary structure, solvent accessibility, relative accessible surface area** (dimension=6)
  * Install DSSP: `sudo apt install dssp`
    * from: https://ssbio.readthedocs.io/en/latest/instructions/dssp.html
    * This installs 3.0.0 but current version is 4.0
    * 4.0 is not compatible with Biopython DSSP module. It can be installed by following:
      * Home: https://swift.cmbi.umcn.nl/gv/dssp/
      * Github: https://github.com/PDB-REDO/dssp
  * `python data_generators/sa_sa_rasa.py`
* **Relative difference** (dimension=6)
