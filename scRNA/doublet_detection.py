import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import os
import scrublet as scr
import scipy.io
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42


parser = argparse.ArgumentParser()

parser.add_argument('-i','--input', type=str,help="The directory of input matrix")
parser.add_argument('-o','--output', type=str,help="The directory of output",default='./')


args = parser.parse_args()



input_dir = args.input
os.system('gunzip '+input_dir+"/*" )

counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx').T.tocsc()
genes = np.array(scr.load_genes(input_dir + '/genes.tsv', delimiter='\t', column=1))

print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
print('Number of genes in gene list: {}'.format(len(genes)))

#os.system('gzip '+input_dir+"/*" )

scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)

doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2,min_cells=3,min_gene_variability_pctl=85,n_prin_comps=30)


print('Running UMAP...')
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
print('Done.')

scrub.plot_embedding('UMAP', order_points=True);
plt.savefig(args.output+"/umap_doublet_plot.png")
plt.cla()



len_array = np.array(range(0,len(predicted_doublets)))

result = list(len_array[list(predicted_doublets)])
print("There are "+str(len(result))+" doublets.")


adata = sc.read_10x_mtx(
    input_dir,  # the directory with the `.mtx` file
    var_names='gene_symbols'               # use gene symbols for the variable names (variables-axis index)
)  

print(adata)
barcodes=adata.obs.index.tolist()
genes=adata.var.index.tolist()
mat=adata.X.todense()
filtered=pd.DataFrame(data=mat,columns=genes,index=barcodes)
filter_name=list(filtered.iloc[result,:].index)
filtered=filtered.drop(filter_name,axis=0) 
filtered.to_csv(args.output+'/cell_gene_filtered_matrix.csv')
