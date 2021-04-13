#!/bin/bash

input=`sed '/^input_1=/!d;s/.*=//' conf_scRNA.txt`
output=`sed '/^output_1=/!d;s/.*=//' conf_scRNA.txt`
organism=`sed '/^organism_1=/!d;s/.*=//' conf_scRNA.txt`
min_genes=`sed '/^min_genes_1=/!d;s/.*=//' conf_scRNA.txt`
max_genes=`sed '/^max_genes_1=/!d;s/.*=//' conf_scRNA.txt`
mitochondrial_gene_percentage=`sed '/^mitochondrial_gene_percentage_1=/!d;s/.*=//' conf_scRNA.txt`
num_of_PCs=`sed '/^num_of_PCs_1=/!d;s/.*=//' conf_scRNA.txt`

python cell_clustering.py -i $input -o $output -m $organism -l $min_genes -x $max_genes -c $mitochondrial_gene_percentage -p $num_of_PCs