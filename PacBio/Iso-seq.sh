#!/bin/bash

bam_file_1=`sed '/^bam_file_1=/!d;s/.*=//' conf_pacbio.txt`
primer_fasta_file_1=`sed '/^primer_fasta_file_1=/!d;s/.*=//' conf_pacbio.txt`
out_path_1=`sed '/^out_path_1=/!d;s/.*=//' conf_pacbio.txt`

cd $out_path_1
b=${bam_file_1##*/}
prefix=${b%.bam}

ccs $bam_file_1 ${prefix}.ccs.bam --min-rq 0.9

lima ${prefix}.ccs.bam $primer_fasta_file_1 ${prefix}.fl.bam --isoseq --peek-guess

file=$(ls | grep .subreadset.xml)
pre_name=${file%.subreadset.xml}

isoseq3 refine ${pre_name}.bam $primer_fasta_file_1 ${prefix}.flnc.bam

isoseq3 cluster ${prefix}.flnc.bam ${prefix}.clustered.bam --verbose --use-qvs