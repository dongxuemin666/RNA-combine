import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import os

parser = argparse.ArgumentParser()

parser.add_argument('-i','--input', type=str,help="The input matrix")
parser.add_argument('-o','--output', type=str,help="The directory of output",default='./')
parser.add_argument('-m','--organism', type=str,help="Human and mouse are supported now",default='Human')
parser.add_argument('-l','--min_genes', type=int,help="The lowest number of genes in cells", default=100)
parser.add_argument('-x','--max_genes', type=int,help="The highest number of genes in cells", default=2500)
parser.add_argument('-c','--mitochondrial_gene_percentage', type=int,help="The lowest percentage of mitochondrial genes", default=10)
parser.add_argument('-p','--num_of_PCs', type=int,help="Number of PCs to calculate cell neighbors ", default=40)

args = parser.parse_args()

dirc='mkdir'+' '+args.output+'/outs'
out=args.output+'/outs'
os.system(dirc)
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80, facecolor='white')


adata = sc.read(args.input,delimiter=',')   


adata.var_names_make_unique()


sc.pl.highest_expr_genes(adata, n_top=20, )
plt.savefig(out+"/highest_expressed_genes.png")
plt.cla()

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
if (args.organism=='Human'):
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
if (args.organism=='Mouse'):
    adata.var['mt'] = adata.var_names.str.startswith('mt-')
  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)

plt.savefig(out+"/violin_for_quality.png")
plt.cla()

#sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
#sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

adata = adata[adata.obs.n_genes_by_counts < args.max_genes, :]
adata = adata[adata.obs.n_genes_by_counts > args.min_genes, :]
adata = adata[adata.obs.pct_counts_mt < args.mitochondrial_gene_percentage, :]

sc.pp.normalize_total(adata, target_sum=1e4)

sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]

sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, svd_solver='arpack')

sc.pl.pca_variance_ratio(adata, log=True)
plt.savefig(out+"/PCA_variance_ratio.png")
plt.cla()

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=args.num_of_PCs)

sc.tl.umap(adata)

sc.tl.leiden(adata)

sc.pl.umap(adata, color=['leiden'])
plt.tight_layout()
plt.savefig(out+"/umap_plot_clusters.png")
plt.cla()


sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
maker_gene=pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']}).head(150)
maker_gene.to_csv(out+'/marker_genes_for_each_cluster.csv')


results_file = out+'/object.h5ad'
adata.write(results_file)
