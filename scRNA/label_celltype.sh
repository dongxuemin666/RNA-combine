#!/bin/bash

object=`sed '/^object_3=/!d;s/.*=//' conf_scRNA.txt`
cell_type=`sed '/^cell_type_3=/!d;s/.*=//' conf_scRNA.txt`
# for example T_cell B_cell NK 
output=`sed '/^output_3=/!d;s/.*=//' conf_scRNA.txt`

python label_celltype.py -i $object -c $cell_type -o $output