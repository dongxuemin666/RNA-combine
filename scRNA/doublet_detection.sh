#!/bin/bash

input=`sed '/^input_7=/!d;s/.*=//' conf_scRNA.txt`
output=`sed '/^output_7=/!d;s/.*=//' conf_scRNA.txt`



if [ -x "$input/features.tsv.gz"]; then
mv $input/features.tsv.gz $input/genes.tsv.gz
fi



python doublet_detection.py -i $input -o $output
