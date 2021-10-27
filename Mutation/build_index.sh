#!/bin/bash

genome=`sed '/^genome_1=/!d;s/.*=//' conf_mutation.txt`
#gatk_path=`sed '/^gatk_path_1=/!d;s/.*=//' conf_mutation.txt`
#samtools_path=`sed '/^samtools_path_1=/!d;s/.*=//' conf_mutation.txt`
known_snp=`sed '/^known_snp_1=/!d;s/.*=//' conf_mutation.txt`

#$samtools_path/
samtools faidx $genome

#$gatk_path/
gatk CreateSequenceDictionary -R $genome -O ${genome%.*}.dict

#ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/v4.0/00-All.vcf.gz
#$gatk_path/gatk IndexFeatureFile -I $known_snp





