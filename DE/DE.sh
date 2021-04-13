#!/bin/bash



method=`sed '/^method=/!d;s/.*=//' conf_DE.txt`
input_matrix=`sed '/^input_matrix=/!d;s/.*=//' conf_DE.txt`
metadata=`sed '/^metadata=/!d;s/.*=//' conf_DE.txt`
output=`sed '/^output=/!d;s/.*=//' conf_DE.txt`
pvalue=`sed '/^pvalue=/!d;s/.*=//' conf_DE.txt`
logfoldchange=`sed '/^logfoldchange=/!d;s/.*=//' conf_DE.txt`



Rscript DE.R -e $method -i $input_matrix -m $metadata -o $output -p $pvalue -f $logfoldchange 