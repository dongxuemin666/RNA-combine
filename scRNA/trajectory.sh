#!/bin/bash


input=`sed '/^input_6=/!d;s/.*=//' conf_scRNA.txt`
subset_celltype=`sed '/^subset_celltype_6=/!d;s/.*=//' conf_scRNA.txt`
output=`sed '/^output_6=/!d;s/.*=//' conf_scRNA.txt`
root=`sed '/^root_6=/!d;s/.*=//' conf_scRNA.txt`


python trajectory.py -i $input -s $subset_celltype -o $output -r $root  