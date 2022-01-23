import sys
sys.path.append("../PROTS_RF_recon")
from multiprocessing import Pool
import time
import os
from objects.PROTS import PROTS

prots = PROTS()
fragments = prots.do_permutation(4, True)
print("# of fragments: ", len(fragments))

# distributed fragment-feature occurance counter
# i=116279 #0-116279
# i = int(os.environ["SLUsRM_ARRAY_TASK_ID"])   
# prots.compute_all_feat_occ_for_all_frag(i)

start_time = time.time()
with Pool(100) as p:
    p.map(prots.compute_all_feat_occ_for_all_frag, range(0, 2000))

print("{} minutes".format((time.time() - start_time)/60))
    