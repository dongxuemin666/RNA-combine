#!/bin/bash


database=`sed '/^database_5=/!d;s/.*=//' conf_scRNA.txt`
genelist=`sed '/^genelist_5=/!d;s/.*=//' conf_scRNA.txt`
organism=`sed '/^organism_5=/!d;s/.*=//' conf_scRNA.txt`


python search_cell_type.py -m $database -g $genelist -r $organism 