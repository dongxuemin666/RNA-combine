#!/bin/bash

cellranger_path=`sed '/^cellranger_path_2=/!d;s/.*=//' conf_scRNA.txt`
cellranger_ref_genome=`sed '/^cellranger_ref_genome_2=/!d;s/.*=//' conf_scRNA.txt`
out=`sed '/^out_2=/!d;s/.*=//' conf_scRNA.txt`
id=`sed '/^id_2=/!d;s/.*=//' conf_scRNA.txt`
fastq_path=`sed '/^fastq_path_2=/!d;s/.*=//' conf_scRNA.txt`
sample_name=`sed '/^sample_name_2=/!d;s/.*=//' conf_scRNA.txt`
expected_cell_number=`sed '/^expected_cell_number_2=/!d;s/.*=//' conf_scRNA.txt`
cores_to_run=`sed '/^cores_to_run_2=/!d;s/.*=//' conf_scRNA.txt`

cd $out

$cellranger_path/cellranger count --id=$id \
                   --transcriptome=$cellranger_ref_genome \
                   --fastqs=$fastq_path \
                   --sample=$sample_name \
                   --expect-cells=$expected_cell_number \
                   --localcores=$cores_to_run 
