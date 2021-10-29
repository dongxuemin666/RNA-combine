import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-i','--object', help="Scanpy object",default='')
parser.add_argument('-g','--plot_gene_umap', help="genes to plot on umap, example: -g CD3E CD8A",default='',nargs='+')
parser.add_argument('-v','--plot_gene_violin', help="genes to plot in violin plot, example: -v CD3E CD8A",default='',nargs='+')
parser.add_argument('-o','--output', help="The directory of output",default='./')
args = parser.parse_args()
if len(args.object):
	adata = sc.read(args.object)
else:
	print('please input the Scanpy object after umap calculation')
	sys.exit(1)
if len(args.plot_gene_umap):
	for i in range(0,len(args.plot_gene_umap)):
		#print(type(args.plot_gene_umap[i]))
		sc.pl.umap(adata, color=args.plot_gene_umap[i])
		plt.savefig(args.output+'/'+args.plot_gene_umap[i]+"_umap_expression.png")
		plt.cla()


if len(args.plot_gene_violin):
	for i in range(0,len(args.plot_gene_violin)):
		sc.pl.violin(adata,args.plot_gene_violin[i] , groupby='leiden')
		plt.savefig(args.output+'/'+args.plot_gene_violin[i]+"_violin_expression.png")
		plt.cla()
