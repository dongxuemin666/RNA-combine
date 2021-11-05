#!/bin/bash

method_2=`sed '/^method_2=/!d;s/.*=//' conf_pacbio.txt`
bam_file_2=`sed '/^bam_file_2=/!d;s/.*=//' conf_pacbio.txt`
reference_genome_2=`sed '/^reference_genome_2=/!d;s/.*=//' conf_pacbio.txt`
out_path_2=`sed '/^out_path_2=/!d;s/.*=//' conf_pacbio.txt`

if [ "$method" == "blasr" ]; then
cd $out_path_2
b=${bam_file_2##*/}
pre_name=${b%.bam}
blasr $bam_file_2 $reference_genome_2 --bam --out ${pre_name}.alignments.bam
fi

if [ "$method" == "minimap2" ]; then
cd $out_path_2
b=${bam_file_2##*/}
pre_name=${b%.bam}
samtools sort -n $bam_file_2 -@ 1 -o ${pre_name}.sorted.bam
samtools bam2fq ${pre_name}.sorted.bam > ${pre_name}.sorted.fastq
minimap2 -ax map-pb $reference_genome_2 ${pre_name}.sorted.fastq > ${pre_name}.sam

fi

