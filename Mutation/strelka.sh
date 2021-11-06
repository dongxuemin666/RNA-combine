#!/bin/bash

normalbam=`sed '/^normalbam_3=/!d;s/.*=//' conf_mutation.txt`
tumorbam=`sed '/^tumorbam_3=/!d;s/.*=//' conf_mutation.txt`
genome=`sed '/^genome_3=/!d;s/.*=//' conf_mutation.txt`
run_dir=`sed '/^run_dir_3=/!d;s/.*=//' conf_mutation.txt`
threads=`sed '/^threads_3=/!d;s/.*=//' conf_mutation.txt`





script=../scripts

raw_dir=$(pwd)

cd ${script}/strelka-2.9.2.centos6_x86_64/bin
strelka_script=$(pwd)
${strelka_script}/configureStrelkaSomaticWorkflow.py \
    --normalBam $normalbam \
    --tumorBam $tumorbam \
    --referenceFasta $genome \
    --runDir $run_dir

cd $run_dir
runWorkflow.py -m local -j $threads
