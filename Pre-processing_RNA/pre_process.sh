#!/bin/bash

fq_path=`sed '/^fq_path=/!d;s/.*=//' conf_prerna.txt`
threads=`sed '/^threads=/!d;s/.*=//' conf_prerna.txt`

#hisat2_path=`sed '/^hisat2_path=/!d;s/.*=//' conf_prerna.txt`
sortmerna_index_path=`sed '/^sortmerna_index_path=/!d;s/.*=//' conf_prerna.txt`
#samtools_path=`sed '/^samtools_path=/!d;s/.*=//' conf_prerna.txt`
#featureCounts_path=`sed '/^featureCounts_path=/!d;s/.*=//' conf_prerna.txt`

hisat_index_path=`sed '/^hisat_index_path=/!d;s/.*=//' conf_prerna.txt`
hisat_index_prefix=`sed '/^hisat_index_prefix=/!d;s/.*=//' conf_prerna.txt`
gtf=`sed '/^gtf=/!d;s/.*=//' conf_prerna.txt`
PE_or_SE=`sed '/^PE_or_SE=/!d;s/.*=//' conf_prerna.txt`
read_length=`sed '/^read_length=/!d;s/.*=//' conf_prerna.txt`

cd ../scripts
sortmerna_path=$(pwd)
trim_path=$(pwd)/Trimmomatic-0.38


flag=1

if [ "$PE_or_SE" == "PE" ]; then
cd $fq_path || exit 1
for file in $(ls | grep _1.fastq)
        do
        pre_name=${file%_1.fastq}
        fq1="${pre_name}_1.fastq"
        fq2="${pre_name}_2.fastq"

        # step 1
        # remove ribosome genes
        mkdir 1.rm_rrna
        cd 1.rm_rrna 
        $sortmerna_path/merge-paired-reads.sh ../$fq1 ../$fq2 ${pre_name}_merge.fastq

        $sortmerna_path/sortmerna --ref $sortmerna_path/rRNA_databases/silva-bac-16s-id90.fasta,$sortmerna_index_path/sortmerna_index/silva-bac-16s-db:\
$sortmerna_path/rRNA_databases/silva-bac-23s-id98.fasta,$sortmerna_index_path/sortmerna_index/silva-bac-23s-db:\
$sortmerna_path/rRNA_databases/silva-arc-16s-id95.fasta,$sortmerna_index_path/sortmerna_index/silva-arc-16s-db:\
$sortmerna_path/rRNA_databases/silva-arc-23s-id98.fasta,$sortmerna_index_path/sortmerna_index/silva-arc-23s-db:\
$sortmerna_path/rRNA_databases/silva-euk-18s-id95.fasta,$sortmerna_index_path/sortmerna_index/silva-euk-18s-db:\
$sortmerna_path/rRNA_databases/silva-euk-28s-id98.fasta,$sortmerna_index_path/sortmerna_index/silva-euk-28s:\
$sortmerna_path/rRNA_databases/rfam-5s-database-id98.fasta,$sortmerna_index_path/sortmerna_index/rfam-5s-db:\
$sortmerna_path/rRNA_databases/rfam-5.8s-database-id98.fasta,$sortmerna_index_path/sortmerna_index/rfam-5.8s-db --reads ${pre_name}_merge.fastq \
--num_alignments 1 --fastx --aligned ${pre_name}_rRNA --other ${pre_name}_non_rRNA --paired_in --log -v -a $threads
        # paired-in: If one of the paired-end reads is Aligned, put both reads into Aligned FASTA/Q file


      $sortmerna_path/unmerge-paired-reads.sh ${pre_name}_non_rRNA.fastq ${pre_name}_non_rRNA_1.fastq ${pre_name}_non_rRNA_2.fastq


        # Step 2
        # Trim reads using Trimmomatic
        cd ..
        mkdir 2.trim
        cd 2.trim
java -jar $trim_path/trimmomatic-0.38.jar PE -threads $threads \
../1.rm_rrna/${pre_name}_non_rRNA_1.fastq \
../1.rm_rrna/${pre_name}_non_rRNA_2.fastq \
-baseout ${pre_name}_nonrRNA_trimmed.fastq \
ILLUMINACLIP:$trim_path/adapters/TruSeq3-PE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$read_length \
    2> ${pre_name}_Trimmomatic.log

        # Step 3
        # Map reads to hg19 reference
        cd ..
        mkdir 3.align
        cd 3.align
        #$hisat2_path/
	hisat2 -t -p $threads -q -x $hisat_index_path/$hisat_index_prefix \
        --rg-id ${pre_name} \
        --rg SM:${pre_name} --rg LB:${pre_name} --rg PL:ILLUMINA \
        -1 ../2.trim/${pre_name}_nonrRNA_trimmed_1P.fastq \
        -2 ../2.trim/${pre_name}_nonrRNA_trimmed_2P.fastq \
        -S ${pre_name}.sam
                2> ${pre_name}_alignment.log

        # sam to bam
        #$samtools_path/
	samtools view -S ${pre_name}.sam -b > ${pre_name}.bam
        # sort default chromosome location
        #$samtools_path/samtools sort ${pre_name}.bam -o ${pre_name}_sorted.bam
        # build index
        #$samtools_path/
	#samtools index ${pre_name}.bam
        # sort by reads name
        #$samtools_path/
	samtools sort  ${pre_name}.bam -o ${pre_name}_nsorted.bam
    	samtools index ${pre_name}_nsorted.bam


        cd ..
        mkdir 4.count
        cd 4.count
        # Step 4
        # Count number of reads on genes
        #$featureCounts_path/
	featureCounts \
                -T $threads \
                -p \
                -t exon \
                -g gene_id \
                -a $gtf \
                -o ${pre_name}_featureCounts.txt \
                ../3.align/${pre_name}_nsorted.bam \
                2> ${pre_name}_featureCounts.log
        sed -i '1d' ${pre_name}_featureCounts.txt
        awk '{print $1,$7}' ${pre_name}_featureCounts.txt >${pre_name}.count
        
        if [ $flag == 1 ]; then
        awk '{print $1}' ${pre_name}_featureCounts.txt >gene_expression_matrix.txt
        join gene_expression_matrix.txt ${pre_name}.count --nocheck-order >gene_expression_matrix_${flag}.txt
        else
        flag_1=$((flag-1))
        join gene_expression_matrix_${flag_1}.txt ${pre_name}.count --nocheck-order >gene_expression_matrix_${flag}.txt
        fi
        
        cd ..
        flag=$((flag+1))
        done
                cd 4.count
        flag=$((flag-1))
        mv gene_expression_matrix_${flag}.txt expression_matrix.txt
        rm gene_expression_matrix* *count
fi





if [ "$PE_or_SE" == "SE" ]; then
cd $fq_path || exit 1
for file in $(ls | grep .fastq)
        do
        
        pre_name=${file%.fastq}

# 
# 
#         #step 1
#         #remove ribosome genes
        mkdir 1.rm_rrna
        cd 1.rm_rrna

        $sortmerna_path/sortmerna --ref $sortmerna_path/rRNA_databases/silva-bac-16s-id90.fasta,$sortmerna_index_path/sortmerna_index/silva-bac-16s-db:\
$sortmerna_path/rRNA_databases/silva-bac-23s-id98.fasta,$sortmerna_index_path/sortmerna_index/silva-bac-23s-db:\
$sortmerna_path/rRNA_databases/silva-arc-16s-id95.fasta,$sortmerna_index_path/sortmerna_index/silva-arc-16s-db:\
$sortmerna_path/rRNA_databases/silva-arc-23s-id98.fasta,$sortmerna_index_path/sortmerna_index/silva-arc-23s-db:\
$sortmerna_path/rRNA_databases/silva-euk-18s-id95.fasta,$sortmerna_index_path/sortmerna_index/silva-euk-18s-db:\
$sortmerna_path/rRNA_databases/silva-euk-28s-id98.fasta,$sortmerna_index_path/sortmerna_index/silva-euk-28s:\
$sortmerna_path/rRNA_databases/rfam-5s-database-id98.fasta,$sortmerna_index_path/sortmerna_index/rfam-5s-db:\
$sortmerna_path/rRNA_databases/rfam-5.8s-database-id98.fasta,$sortmerna_index_path/sortmerna_index/rfam-5.8s-db --reads $fq_path/${pre_name}.fastq \
--num_alignments 1 --fastx --aligned ${pre_name}_rRNA --other ${pre_name}_non_rRNA --paired_in --log -v -a $threads
        # paired-in: If one of the paired-end reads is Aligned, put both reads into Aligned FASTA/Q file

# 
#         #Step 2
#         #Trim reads using Trimmomatic
        cd ..
mkdir 2.trim
cd 2.trim
java -jar $trim_path/trimmomatic-0.38.jar SE -threads $threads \
../1.rm_rrna/${pre_name}_non_rRNA.fastq \
${pre_name}_nonrRNA_trimmed.fastq ILLUMINACLIP:$trim_path/adapters/TruSeq3-SE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$read_length 2> ${pre_name}_Trimmomatic.log

#         # Step 3
#         # Map reads to hg19 reference
        cd ..
        mkdir 3.align
        cd 3.align
        #$hisat2_path/
	hisat2 -t  -p $threads  -x $hisat_index_path/$hisat_index_prefix --rg-id ${pre_name} --rg SM:${pre_name} --rg LB:${pre_name} --rg PL:ILLUMINA -U ../2.trim/${pre_name}_nonrRNA_trimmed.fastq -S ${pre_name}.sam \
                2> ${pre_name}_alignment.log

        # sam to bam
        #$samtools_path/
	samtools view -S ${pre_name}.sam -b > ${pre_name}.bam
        # sort default chromosome location
        #$samtools_path/samtools sort ${pre_name}.bam -o ${pre_name}_sorted.bam
        # build index
        #$samtools_path/
	#samtools index ${pre_name}.bam
        # sort by reads name
        #$samtools_path/
	samtools sort  ${pre_name}.bam -o ${pre_name}_nsorted.bam
#     
	samtools index ${pre_name}_nsorted.bam
# 
# 
         cd ..
         mkdir 4.count
        cd 4.count
        # Step 4
        # Count number of reads on genes
        #$featureCounts_path/
	featureCounts \
                -T $threads \
                -p \
                -t exon \
                -g gene_id \
                -a $gtf \
                -o ${pre_name}_featureCounts.txt \
                ../3.align/${pre_name}_nsorted.bam \
                2> ${pre_name}_featureCounts.log

        sed -i '1d' ${pre_name}_featureCounts.txt
        awk '{print $1,$7}' ${pre_name}_featureCounts.txt >${pre_name}.count
        
        if [ $flag == 1 ]; then
        awk '{print $1}' ${pre_name}_featureCounts.txt >gene_expression_matrix.txt
        join gene_expression_matrix.txt ${pre_name}.count --nocheck-order >gene_expression_matrix_${flag}.txt
        else
        flag_1=$((flag-1))
        join gene_expression_matrix_${flag_1}.txt ${pre_name}.count --nocheck-order >gene_expression_matrix_${flag}.txt
        fi
        
        cd ..
        flag=$((flag+1))
        done
        cd 4.count
        flag=$((flag-1))
        mv gene_expression_matrix_${flag}.txt expression_matrix.txt
        rm gene_expression_matrix* *count
        
fi
