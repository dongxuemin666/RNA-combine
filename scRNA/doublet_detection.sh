#!/bin/bash

input=`sed '/^input_7=/!d;s/.*=//' conf_scRNA.txt`
output=`sed '/^output_7=/!d;s/.*=//' conf_scRNA.txt`

python doublet_detection.py -i $input -o $output