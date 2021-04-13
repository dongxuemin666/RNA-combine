import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import os

parser = argparse.ArgumentParser()

parser.add_argument('-i','--input', type=str,help="The scanpy object produced by scrna.py")
parser.add_argument('-s','--subset_celltype', type=str,help="Due to not all cell types have evolutionary relations,you should provide a list of cell types that indeed have evolutionary relations or may have",default='',nargs='+')
parser.add_argument('-o','--output', help="The directory of output",default='./')
parser.add_argument('-r','--root', type=str, help="The root cell types",default='')




args = parser.parse_args()


sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80, facecolor='white')


adata = sc.read(args.input)

if ('cell_type' in adata.obs.columns.values.tolist()):
    adata=adata[adata.obs['cell_type'].isin(args.subset_celltype)]
    sc.tl.paga(adata, groups='cell_type')
    sc.pl.paga_compare(adata, threshold=0.03, title='', right_margin=0.2, size=10, edge_width_scale=0.5,legend_fontsize=12, fontsize=12, frameon=False, edges=True, save=False)
    plt.savefig(args.output+"/trajectory_plot.png")
    plt.cla()    
    adata.uns['iroot'] = np.flatnonzero(adata.obs['cell_type'] == args.root)[0]
    sc.tl.dpt(adata)
    sc.pl.umap(adata, color=['cell_type', 'dpt_pseudotime'], legend_loc='on data')
    plt.savefig(args.output+"/psedotime_plot.png")
    plt.cla()
elif ('leiden' in adata.obs.columns.values.tolist()):
    adata=adata[adata.obs['leiden'].isin(args.subset_celltype)]
    sc.tl.paga(adata, groups='leiden')
    sc.pl.paga_compare(adata, threshold=0.03, title='', right_margin=0.2, size=10, edge_width_scale=0.5,legend_fontsize=12, fontsize=12, frameon=False, edges=True, save=False)
    plt.savefig(args.output+"/trajectory_plot.png")
    plt.cla() 
    adata.uns['iroot'] = np.flatnonzero(adata.obs['leiden']  == args.root)[0]
    sc.tl.dpt(adata)
    sc.pl.umap(adata, color=['leiden', 'dpt_pseudotime'], legend_loc='on data')
    plt.savefig(args.output+"/psedotime_plot.png")




