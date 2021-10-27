# RNA-combine
A toolkit for comprehensive analysis on multiple platform transcriptome data 

## Schematic
RNA-combine mainly contains three modules, NGS bulk RNA-seq data analysis, scRNA-seq data analysis, PacBio data analysis

![schematic](https://github.com/dongxuemin666/RNA-combine/blob/master/Schematic.svg)



## Usage
### 0 Environment building and overview

Download the codes to your directory and build the environment using below code  
```Bash
conda install --file=requirement.txt
```

To use function modules in each directories, first, you need to edit parameters for concerned modules in configuration file, such as the path of softwares, the path of raw sequencing data. This process is quite simple, because all parameters are explained in configration file. Then, run corresponding bash commands. 
### 1 Bulk NGS RNA-seq data analysis
#### 1.1 Pre-processing_RNA
Pre-process the bulf RNA-seq data.  
Inputs are fastq files.  
Outputs mainly include aligned bam files and gene count matrix, with intermediate files saved in 4 directories. *1.rm_rrna* save sequences after removing rRNA sequences,with log files recording rRNA information for each sample. *2.trim* save sequences after cutting adapters and low quality bases, with log files recording information of cutted sequences. *3.align* save aligned bam files, with log files recording alignment information. *4.count* save gene count files which annotate sequences to annoated reference gtf files, with log and summary files recording annotation details. 
 
To use it, first you need to edit the configuration file "conf_prerna.txt".  
Second, it is optional to get the index files for rRNA and reference genome, depending on whether you have index files or not. 


```Bash
./build_index.sh
```
Third, simply run pre_process.sh
```Bash
./pre_process.sh
```
Then all results will be outputted to same directory of fastq files.
#### 1.2 Differential（DE） analysis
Explore the differential expressed  genes in different conditions, utilizing the bulf RNA-seq data.  
Input is gene count file.  
For outputs, "sample_correlation_plot.png" shows the sample distances.  
"volcano_plot.png" illustrates the significantly differential genes in differential conditions.  
"de_gene_foldchange_above_0.csv" and "de_gene_foldchange_above_0.csv" give more details about DE genes. Be sure to read comparison mode in the output.  
To use it, first you need to edit the configuration file "conf_DE.txt".  
Second, simply run "DE.sh"
```Bash
./DE.sh
```
#### 1.3 Function annotation
Function enrichment based on significant DE genes.   
Input is a gene list.  
For outputs, *function_annotation.pdf* shows the enriched Go terms of KEGG pathways.  
To use it, first, you need to edit the configuration file "conf_function.txt".  
Second you simply run "function.sh"  
```Bash
./function.sh
```
Then all results will be outputted to output directory.
#### 1.4 Gene co-expression relation 
Explore the gene co-expression pattern based on Spearman, Pearson, PCA-PMI(Juan Zhao et al. PNAS, 2016).    
Input is a gene exprssion matrix, rows are genes, columns are samples.    
For outputs, pdf files illustrate gene co-expression pattern, csv files give more details about coefficients.   
To use it, first, you need to edit the configuration file "conf_relation.txt".  
Second you simply run "relation.sh"  
```Bash
./relation.sh
```
Then all results will be outputted to output directory.
#### 1.5 Differential gene co-expression pattern
Explore differential gene co-expression patterns in different conditions, mainly utilizing SSN(Xiaoping Liu et al. Nucleic Acids Research, 2016).  
For input, normal_gene_matrix is used to build background network, case_gene_matrix was added to background network to form perturbation network. csv format is needed. 
Then the differential gene relations(differential edges) were calculated between normal samples and case samples.  
For outputs, the coefficients of differential edges are given, including Pearson corelation deviation betwwen normal samples and case samples, p-value for differential edges.
To use it, first, you need to edit the configuration file "conf_ssn.txt".  
Second you simply run "SSN.sh"  
```Bash
./SSN.sh
```
Then all results will be outputted to output directory.
#### 1.6 Variant detection
Detect variants based on RNA-seq data, including SNP and INDEL.  
Inputs are bam files produced by pre_process.sh  
For outputs, *xxx.pass.snp.indel.vcf.gz* file give details of called SNPs and INDELs.  
To use it, first, you need to edit the configuration file "conf_mutation.txt".  
Second you simply run "mutaion.sh"  
```Bash
./mutaion.sh
```
Then all results will be outputted to the sample directory of bam files.
#### 1.7 RNA splicing analysis
RNA splicing analysis includes three kinds of methods, including exon-based method: DEXseq, transcript-based method(not only)：StringTie+ballgown, event-based method: rMATS  

For outputs of DEXseq, *Significant_differential_exons.csv* give details about differential exons in different conditions, while *exon_plot_for_one_gene.pdf* visualize differential exon information for specific gene. 
For outputs of StringTie+ballgown, *Significant_differential_transcripts.csv* give details about differential transcripts in different conditions, while *Transcipt_plot_for_one_gene.pdf* visualize differential transcript information for specific gene. 
For outputs of rMATS, txt files give details about differential splicing events in different conditions(Search rMATS for description for each file), while *event_plot* save visualization information for all differential splicing events. 

To use it, first, you need to edit the configuration file "conf_splicing.txt". Choose the method, then edit the corresponding module.  

Second you simply run "Splicing_pipeline.sh"  
```Bash
./Splicing_pipeline.sh
```
Then all results will be outputted to output directory.  
### 2 scRNA data analysis
For scRNA data analysis, 7 functions are introduced.


#### 2.1 cellranger.sh
Turn 10X fastq files to cell_gene matrix, only support data produced by 10X platform. This process mainly utilizes CellRanger software.  
The gene-cell matrix was at "outs/filtered_feature_bc_matrix" directory.  which could be the input of doublet_detection.sh module.  
To use it, first you need to edit cell_ranger module of the configuration file conf_scRNA.txt. 
Second you need simply run cell_clustering.sh 
```Bash
./cellranger.sh
```
Then the results will be outputted to output directory
#### 2.2 doublet_detection.sh
Detect doublets based on gene-cell matrix 
For outputs, umap_doublet_plot.png is UMAP plot of doublets and single cells. cell_gene_filtered_matrix.csv is gene-cell matrix under filtering of doublets, which could be input of cell_clustering.sh moudule  
To use it, first you need to edit doublet_detection.sh module of the configuration file conf_scRNA.txt. 
Second you need simply run doublet_detection.sh 
```Bash
./doublet_detection.sh
```
Then the results will be outputted to output directory 
#### 2.3 cell_clustering.sh 
Designed for quality control, data dimension reduction, cell clustering, marker gene detection for each cluster et. This process mainly utilizes scanpy python package.  
The input is 10X gene-cell matrix.  
For output, "highest_expressed_genes.png" illustrates the highest expressed genes in this sample.  
"PCA_variance_ratio.png" shows the variance of PCs,  PC1 to 'elbow' PC of this plot are commonly used for subsequent analysis, which is corresponding to "num_of_PCs_1" parameter in configurarion file.  
"violin_for_quality.png" shows the quality of your data includes the number of genes, the number of UMIs, the percentage of mitochondrial genes.  
"umap_plot_clusters.png" shows the cell clusters enriched.
"object.h5ad" is the scanpy object that contains all information of analyzed data, which could be input for othe analyses.  
To use it, first you need to edit cell_clustering module of the configuration file conf_scRNA.txt.  
Second you need simply run cell_clustering.sh 
```Bash
./cell_clustering.sh
```
Then the results will be outputted to output directory
#### 2.4 plot_gene_scrna.sh
Plot the gene expression in the form of violin plot or UMAP plot.  
First you need to edit plot_gene_scrna module of the configuration file conf_scRNA.txt .  
Second you need simply run plot_gene_scrna.sh
```Bash
./plot_gene_scrna.sh
```
Then the results will be outputted to output directory
#### 2.5 search_cell_type.sh 
Predict the cell type based on your provided significant gene(SG) list. For each SG, the program searches correpsonding cell type to it in database. Finally, the program summarize what cell type it is most like. Right now, cell marker datasets in CellMarker database and PanglaoDB database are downloaded and clearned, you can choose to use them.
To use it, first you need to edit search_cell_type module of the configuration file conf_scRNA.txt .  
Second you need simply run search_cell_type.sh 
```Bash
./search_cell_type.sh
```
Then the results will be outputted to output directory

#### 2.6 label_celltype.sh
After searching the cell type of each cell cluster, you need to label the cell type information in scRNA object.  
To use it, first you need to edit label_celltype module of the configuration file conf_scRNA.txt.  
Second you need simply run label_celltype.sh 
```Bash
./label_celltype.sh
```
Then the results will be outputted to output directory
#### 2.7 trajectory.sh
You may want to detect the evolutionary relations among different cell types, you can use this function to do so.
To use it, first you need to edittrajectory module of the configuration file conf_scRNA.txt.  
Second you need simply run trajectory.sh
```Bash
./trajectory.sh
```
Then the results will be outputted to output directory

### 3 PacBio RNA data analysis
For PacBio RNA data analysis, 2 functions are introduced.

Run install.sh before run pipelines for installing needed software.

#### 3.1 Iso-seq.sh
Turn raw read bam files to full-length, non-concatemer, unique transcript files. 
Fasta outputs are divided into HQ and LQ files, and only one concensus for one transcript cluster. 
To use it, first you need to edit Iso-seq.sh module of the configuration file conf_pacbio.txt. 
Second you need simply run Iso-seq.sh 
```Bash
./Iso-seq.sh
```
Then the results will be outputted to output directory. Among all of results, (prefix).hq.fasta.gz are isoforms with predicted accuracy ≥ 0.99, (prefix).lq.fasta.gz are isoforms with predicted accuracy < 0.99.
#### 3.2 map.sh
Align sub-reads to reference genome 
For output, ...alignments.bam file is aligned bam file 
To use it, first you need to edit map.sh module of the configuration file conf_pacbio.txt. 
Second you need simply run map.sh 
```Bash
./map.sh
```
Then the results will be outputted to output directory 

## Dependencies
|Methods|Description|
| ------------- |:-------------:| 
|Sortmerna (Kopylova E. et al., 2012)	|Removing rRNA sequences |
|Hisat2 (Kim D et al., 2015)	| Aligning reads to reference genome |
|samtools	| Converting file formats and auxiliary functionalities |
|featureCounts (Yang Liao et al., 2013)	| Counting number of reads on genes | 
|DESeq2(Love, M.I. et al., 2014), limma(Ritchie ME et al., 2015), edgeR(Robinson MD, 2010)|	 Differential gene analysis |
|ClusterProfiler (Yu G et al., 2012)     | Function enrichment|
|DEXseq(Anders S et al., 2012), StringTie(Pertea M et al., 2015), rMATS(Shihao Shen et al., 2014)	| Splicing site detection |
|GATK |	Variant calling |
|CellRanger	| Produce barcode-gene matrix based on scRNA fastq files |
|Scrublet (Samuel L. et al., 2019) |	Doublet detection in scRNA data |
|Scanpy (Wolf, F. et al., 2018) |	filtering, normalization, clustering, dimension reduction, trajectory analysis (and so on) of scRNA data |
|isoseq3 |	Produce full-length, non-concatemer, unique transcripts based on PacBio raw bam read files |
|blasr	| Map PacBio reads to reference genome |



## Maintainer
DONG Xuemin, Institute of Zoology, dongxuemin18@mails.ucas.ac.cn


## References
Anders S, Reyes A et al. (2012). Detecting differential usage of exons from RNA-seq data. Genome Research, 22, 2008-2017.  
Anthony M. Bolger, Marc Lohse and et al. (2014) Trimmomatic: a flexible trimmer for Illumina sequence data, Bioinformatics, 30, 2114–2120.  
Identification of Cell Doublets in Single-Cell Transcriptomic Data. Cell Systems, 8, 281-291.  
Kim D, Langmead B et al. (2015) HISAT: a fast spliced aligner with low memory requirements. Nature Methods, 12, 357–360.  
Kopylova E., Noé L. et al. (2012) SortMeRNA: Fast and accurate filtering of ribosomal RNAs in metatranscriptomic data, Bioinformatics, 28, 3211–3217.  
Love, M.I., Huber, W. et al. (2014) Moderated estimation of fold change and dispersion for RNA seq data with DESeq2. Genome Biology, 15,550.  
Pertea M, Pertea GM et al. (2015) StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. Nature Biotechnology, 33, 290–295.  
Ritchie ME, Phipson B et al. (2015) limma powers differential expression analyses for RNA sequencing and microarray studies. Nucleic Acids Research, 43, 47.  
Robinson MD, McCarthy DJ et al. (2010) edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics, 26, 139-140.  
Samuel L. Wolock, Romain Lopez et al. (2019) Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data. Cell Systems, 8, 281-291.  
Shihao Shen, Juw Won Park et al. (2014) rMATS: Robust and flexible detection of differential alternative splicing from replicate RNA-Seq data. PNAS, 111, 5593-5601.  
Wolf, F., Angerer, P. & Theis, F. (2018) SCANPY: large-scale single-cell gene expression data analysis. Genome Biol, 19, 15.  
Yang Liao, Gordon K Smyth et al. (2014) featureCounts: an efficient general-purpose program for assigning sequence reads to genomic features. Bioinformatics, 30, 923-930.  
Yu G, Wang L et al. (2012) clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology, 16, 284-287.  
