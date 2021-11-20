#!/bin/bash

bam_path=`sed '/^bam_path_2=/!d;s/.*=//' conf_mutation.txt`
#samtools_path=`sed '/^samtools_path_2=/!d;s/.*=//' conf_mutation.txt`
#gatk_path=`sed '/^gatk_path_2=/!d;s/.*=//' conf_mutation.txt`
genome=`sed '/^genome_2=/!d;s/.*=//' conf_mutation.txt`
known_snp=`sed '/^known_snp_2=/!d;s/.*=//' conf_mutation.txt`
interval=`sed '/^interval_2=/!d;s/.*=//' conf_mutation.txt`
sample_amount=`sed '/^sample_amount_2=/!d;s/.*=//' conf_mutation.txt`
cd $bam_path || exit 1
 
mkdir outs

 for file in $(ls | grep .bam)
         do
 cd outs
 #bam file id
 pre_name=${file%.bam}
 
 
 
 
 
 
 echo $pre_name>>id.txt
#sort
#$samtools_path/
samtools sort -o ${pre_name}.sorted.bam $bam_path/${pre_name}.bam
#mark dupicates
time gatk MarkDuplicates -I ${pre_name}.sorted.bam -O ${pre_name}.markdup.bam -M ${pre_name}.metrics --CREATE_INDEX
# #build calibration model
# #time $gatk_path/gatk BaseRecalibrator -R $genome -I ${pre_name}.markdup.bam -O ${pre_name}.recal.table --known-sites $known_snp
# #calibration
# #time $gatk_path/gatk ApplyBQSR -R $genome -I ${pre_name}.markdup.bam -bqsr ${pre_name}.recal.table -O ${pre_name}.recal.bam
# #build intermediate files
# #time $gatk_path/gatk HaplotypeCaller -R $genome -I ${pre_name}.recal.bam -ERC GVCF --dbsnp $known_snp -O ${pre_name}.snps.indels.vcf


time gatk SplitNCigarReads -R $genome \
-I ${pre_name}.markdup.bam \
-O  ${pre_name}.markdup.split.bam










if [ "$sample_amount" == "multiple" ]; then

time gatk HaplotypeCaller -R $genome -I ${pre_name}.markdup.split.bam -ERC GVCF  -O ${pre_name}.snps.indels.vcf

fi

if [ "$sample_amount" == "single" ]; then

time gatk HaplotypeCaller -R $genome -I ${pre_name}.markdup.split.bam -O genotype.vcf.gz


fi



#ID for vcf files
 echo ${pre_name}.snps.indels.vcf>>vcf_files.txt
 cd ..
         done
 cd outs
 
 
 
 
# #For  GenomicsDBImport parameters
if [ "$sample_amount" == "multiple" ]; then

paste id.txt vcf_files.txt > vcf_map.txt
#merge
time gatk  GenomicsDBImport --genomicsdb-workspace-path my_database -L $interval \
--sample-name-map vcf_map.txt

time gatk GenotypeGVCFs -R $genome -V gendb://my_database -O genotype.vcf.gz


fi

#if [ "$sample_amount" == "single" ]; then

#time cp ${pre_name}.snps.indels.vcf genotype.vcf
#gzip genotype.vcf


#fi

# select SNPs
time gatk SelectVariants -select-type SNP -V genotype.vcf.gz -O genotype.snp.vcf.gz

#filter SNPs
time gatk VariantFiltration -V genotype.snp.vcf.gz \
    --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "Filter" \
    -O genotype.snp.filter.vcf.gz
# select Indels
time gatk SelectVariants -select-type INDEL \
    -V genotype.vcf.gz \
    -O genotype.indel.vcf.gz
#fileter Indels
time gatk VariantFiltration -V genotype.indel.vcf.gz \
--filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "Filter" \
    -O genotype.indel.filter.vcf.gz

# merge filtered SNP and Indel

time gatk MergeVcfs \
    -I genotype.snp.filter.vcf.gz -I genotype.indel.filter.vcf.gz \
    -O genotype.filter.snp.indel.vcf.gz
    
gatk SelectVariants -R $genome -V genotype.filter.snp.indel.vcf.gz \
-O genotype.pass.snp.indel.vcf.gz -select "vc.isNotFiltered()"

# rm no-usage files
 # rm -rf genotype.snp.vcf.gz genotype.snp.vcf.gz.tbi genotype.indel.vcf.gz genotype.indel.vcf.gz.tbi \
 # genotype.filter.snp.indel.vcf.gz genotype.filter.snp.indel.vcf.gz.tbi genotype.vcf.gz genotype.vcf.gz.tbi \
 # *txt \
 # *recal.bam *recal.bai *recal.table *metrics *markdup.bam *markdup.bai my_database *sorted.bam
