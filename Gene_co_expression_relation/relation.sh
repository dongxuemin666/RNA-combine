#!/bin/bash

method=`sed '/^method=/!d;s/.*=//' conf_relation.txt`
input_matrix=`sed '/^input_matrix=/!d;s/.*=//' conf_relation.txt`
output=`sed '/^output=/!d;s/.*=//' conf_relation.txt`



Rscript relation.R -e $method -i $input_matrix -o $output