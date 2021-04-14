#!/bin/bash


normal_gene_matrix=`sed '/^normal_gene_matrix=/!d;s/.*=//' conf_ssn.txt`
case_gene_matrix=`sed '/^case_gene_matrix=/!d;s/.*=//' conf_ssn.txt`
pvalue=`sed '/^pvalue=/!d;s/.*=//' conf_ssn.txt`
out=`sed '/^out=/!d;s/.*=//' conf_ssn.txt`



python interaction.py -i $normal_gene_matrix -o $out

python compute_SSN.py -background="${out}/all_gene_pairs.txt" -pvalue=$pvalue -ref="$normal_gene_matrix" -sample="$case_gene_matrix" -out="$out"
