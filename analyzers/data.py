import sys
sys.path.append("../PROTS_RF_recon")

import pandas as pd


def print_class_distribution(df):
    stabilizing = df[df["ddG"]>=0]
    destabilizing = df[df["ddG"]<0]
    proteins = df["pdb_id"].unique().tolist()
    print("#-stabilizing: ", stabilizing.shape, "#-destabilizing: ", destabilizing.shape, "#-proteins: ", len(proteins))

print_class_distribution(pd.read_csv("data/dataset_5_train.csv"))
print_class_distribution(pd.read_csv("data/dataset_5_test.csv"))