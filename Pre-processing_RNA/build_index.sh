#!/bin/bash

sortmerna_index_path=`sed '/^sortmerna_index_path_0=/!d;s/.*=//' conf_prerna.txt`
hisat_index_path=`sed '/^hisat_index_path_0=/!d;s/.*=//' conf_prerna.txt`
hisat_index_prefix=`sed '/^hisat_index_prefix_0=/!d;s/.*=//' conf_prerna.txt`
genome=`sed '/^genome_0=/!d;s/.*=//' conf_prerna.txt`
#hisat2_path=`sed '/^hisat2_path_0=/!d;s/.*=//' conf_prerna.txt`

cd ../scripts
sortmerna_path=$(pwd)


cd $sortmerna_index_path || exit 1
mkdir sortmerna_index
cd sortmerna_index
$sortmerna_path/indexdb_rna --ref $sortmerna_path/rRNA_databases/silva-bac-16s-id90.fasta,./silva-bac-16s-db:\
$sortmerna_path/rRNA_databases/silva-bac-23s-id98.fasta,./silva-bac-23s-db:\
$sortmerna_path/rRNA_databases/silva-arc-16s-id95.fasta,./silva-arc-16s-db:\
$sortmerna_path/rRNA_databases/silva-arc-23s-id98.fasta,./silva-arc-23s-db:\
$sortmerna_path/rRNA_databases/silva-euk-18s-id95.fasta,./silva-euk-18s-db:\
$sortmerna_path/rRNA_databases/silva-euk-28s-id98.fasta,./silva-euk-28s:\
$sortmerna_path/rRNA_databases/rfam-5s-database-id98.fasta,./rfam-5s-db:\
$sortmerna_path/rRNA_databases/rfam-5.8s-database-id98.fasta,./rfam-5.8s-db

cd $hisat_index_path || exit 1
#$hisat2_path/
hisat2-build -p 4 $genome $hisat_index_prefix
