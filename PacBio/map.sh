#!/bin/bash

bam_file_2=`sed '/^bam_file_2=/!d;s/.*=//' conf_pacbio.txt`
reference_genome_2=`sed '/^reference_genome_2=/!d;s/.*=//' conf_pacbio.txt`
out_path_2=`sed '/^out_path_2=/!d;s/.*=//' conf_pacbio.txt`

cd $out_path_2
b=${bam_file_2##*/}
pre_name=${b%.bam}
blasr $bam_file_2 $reference_genome_2 --bam --out ${pre_name}.alignments.bam