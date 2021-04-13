#!/bin/bash


input_significant_gene_file=`sed '/^input_significant_gene_file=/!d;s/.*=//' conf_function.txt`
gene_type=`sed '/^gene_type=/!d;s/.*=//' conf_function.txt`
number_of_significant_genes=`sed '/^number_of_significant_genes=/!d;s/.*=//' conf_function.txt`
output=`sed '/^output=/!d;s/.*=//' conf_function.txt`
function_type=`sed '/^function_type=/!d;s/.*=//' conf_function.txt`
organism=`sed '/^organism=/!d;s/.*=//' conf_function.txt`




Rscript function.R -i $input_significant_gene_file -g $gene_type -n $number_of_significant_genes -o $output -f $function_type -m $organism 