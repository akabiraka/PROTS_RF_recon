import imp
import sys
sys.path.append("../PROTS_RF_recon")

import os
from objects.PROTS import PROTS

prots = PROTS()
fragments = prots.do_permutation(4, False)
print("# of fragments: ", len(fragments))

# distributed fragment-feature occurance counter
# i=116279 #0-116279
i = int(os.environ["SLURM_ARRAY_TASK_ID"])   
prots.compute_all_feat_occ_for_all_frag(i)