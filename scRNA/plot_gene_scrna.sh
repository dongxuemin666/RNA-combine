#!/bin/bash


object=`sed '/^object_4=/!d;s/.*=//' conf_scRNA.txt`
plot_gene_umap=`sed '/^plot_gene_umap_4=/!d;s/.*=//' conf_scRNA.txt`
plot_gene_violin=`sed '/^plot_gene_violin_4=/!d;s/.*=//' conf_scRNA.txt`
output=`sed '/^output_4=/!d;s/.*=//' conf_scRNA.txt`


python plot_gene_scrna.py -i $object -g $plot_gene_umap -v $plot_gene_violin -o $output