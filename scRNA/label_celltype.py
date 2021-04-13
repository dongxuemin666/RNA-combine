import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-i','--object', help="Scanpy object",default='')
parser.add_argument('-c','--cell_type', help="cell types from cluster, starting from cluster1, example: -g B_cell T_cell",default='',nargs='+')
parser.add_argument('-o','--output', help="The directory of output",default='./')
args = parser.parse_args()
if len(args.object):
    adata = sc.read(args.object)
    new_cluster_names = args.cell_type
    adata.obs['cell_type']=adata.obs['leiden']
    adata.rename_categories('cell_type', new_cluster_names)
    sc.pl.umap(adata, color='cell_type', legend_loc='on data', title='', frameon=False)
    adata.write(args.output+'/'+"object.h5ad")
    plt.savefig(args.output+'/'+"cell_type_plot.png")
    plt.cla()

else:
	print('please input the Scanpy object')
	sys.exit(1)



