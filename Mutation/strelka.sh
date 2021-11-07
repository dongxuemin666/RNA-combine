#!/bin/bash

normalsam=`sed '/^normalsam_3=/!d;s/.*=//' conf_mutation.txt`
tumorsam=`sed '/^tumorsam_3=/!d;s/.*=//' conf_mutation.txt`
genome=`sed '/^genome_3=/!d;s/.*=//' conf_mutation.txt`
run_dir=`sed '/^run_dir_3=/!d;s/.*=//' conf_mutation.txt`
threads=`sed '/^threads_3=/!d;s/.*=//' conf_mutation.txt`
region=`sed '/^region_3=/!d;s/.*=//' conf_mutation.txt`


cp $normalsam $run_dir/normal.sam
cp $tumorsam $run_dir/tumor.sam

sed -i '/^\@SQ.*\_/d' $run_dir/normal.sam
sed -i '/^\@SQ.*\_/d' $run_dir/tumor.sam

samtools view -bS $run_dir/normal.sam  >$run_dir/normal.bam
samtools view -bS $run_dir/tumor.sam    >$run_dir/tumor.bam
samtools sort $run_dir/normal.bam -o $run_dir/normal.sort.bam
samtools sort $run_dir/tumor.bam -o $run_dir/tumor.sort.bam

samtools index $run_dir/normal.sort.bam
samtools index $run_dir/tumor.sort.bam
samtools view -b -h $run_dir/normal.sort.bam $region >$run_dir/normal.filtered.bam
samtools view -b -h $run_dir/tumor.sort.bam $region >$run_dir/tumor.filtered.bam


#samtools sort $run_dir/normal.filtered.bam -o $run_dir/normal.filtered.sort.bam
#samtools sort $run_dir/tumor.filtered.bam -o $run_dir/tumor.filtered.sort.bam
samtools index $run_dir/normal.filtered.bam
samtools index $run_dir/tumor.filtered.bam


script=../scripts

raw_dir=$(pwd)

cd ${script}/strelka-2.9.2.centos6_x86_64/bin
strelka_script=$(pwd)
${strelka_script}/configureStrelkaSomaticWorkflow.py \
    --normalBam $run_dir/normal.filtered.bam \
    --tumorBam $run_dir/tumor.filtered.bam \
    --referenceFasta $genome \
    --exome \
    --runDir $run_dir
    #--callRegions $region

cd $run_dir
runWorkflow.py -m local -j $threads
