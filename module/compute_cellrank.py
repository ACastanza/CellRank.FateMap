import os, sys, re, argparse, shutil, warnings
from optparse import OptionParser

import anndata
import anndata as ad
import scanpy as sc
import scvelo as scv
import cellrank as cr
import numpy as np

__author__ = "Anthony S. Castanza"
__email__ = "acastanza@ucsd.edu"
__version__="1.0.0"

# warnings.simplefilter("ignore", category=UserWarning)
# warnings.simplefilter("ignore", category=FutureWarning)
# warnings.simplefilter("ignore", category=DeprecationWarning)

def main():
    usage = "%prog [options]" + "\n"
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input-file", action="store",
                    dest="input_file", help="h5ad (anndata) file.")
    ap.add_argument("-p", "--plot", default="png", action="store", dest="plot",
                    help="Save velocity plots as png or svg")
    ap.add_argument("-o", "--out", default="result", action="store",
                    dest="output", help="Output file basename")
    ap.add_argument("-j", "--cpu", default="1", action="store", dest="ncores",
                    help="CPU cores to use for transition dynamics calculation")

    options = ap.parse_args()

    adata = anndata.read_h5ad(options.input_file)

# Detect/create clustering
    if "clusters" in list(adata.obs):
        cluster_type = "clusters"
        cluster_out = "dataset_clusters"
        print("Found 'clusters' key in dataset. We'll use this for plots and any differential kinetics.\n")
    elif "clusters" not in list(adata.obs):
        if "leiden" in list(adata.obs):
            cluster_type = "leiden"
            cluster_out = "leiden_clusters"
            print("Found 'leiden' clustering in dataset. We'll use this for plots and any differential kinetics.\n")
        elif "leiden" not in list(adata.obs):
            if "walktrap" in list(adata.obs):
                cluster_type = "walktrap"
                cluster_out = "walktrap_clusters"
                print(
                    "Found 'walktrap' clustering in dataset. We'll use this for plots and any differential kinetics.\n")
            else:
                print("Didn't find any clustering in dataset, clustering data using method: 'leiden'.\nWe'll use this for plots and any differential kinetics.\n")
                sc.tl.leiden(adata)
                cluster_type = "leiden"
                cluster_out = "leiden_clusters"

    cr.tl.terminal_states(adata, cluster_key=cluster_type, weight_connectivities=0.2)
    cr.pl.terminal_states(adata, save=, save=options.output + "_CellRank_FateMap." + options.plot)

if __name__ == '__main__':
    main()
