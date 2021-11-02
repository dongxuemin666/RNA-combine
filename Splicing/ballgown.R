if(!require("getopt")) install.packages("getopt",repos="http://cran.rstudio.com/")
library(getopt)
library(tools)
if(!require("ballgown")) BiocManager::install("ballgown")
if(!require("genefilter")) BiocManager::install("genefilter")
if(!require("dplyr")) install.packages("dplyr")
library(ballgown)
library(genefilter)
library(dplyr)


spec <- matrix(
  c(
    "metaData",  "m", 1, "character", "meta-data of the samples, introduce their conditions",
    "qvalue", "q", 1,"double","This is the highest q_value for detetcting significant differential transcripts",
    "output","o",1,"character","This is output directory",
    "Gene_to_plot","g",1,"character","Plot transcript expression for specific gene",
    "help",   "h", 0, "double",  "This is Help,input_matrix,metadata,pvalue,logfoldchange are mandatory parameters"
  ),
  byrow=TRUE, ncol=5)


opt <- getopt(spec=spec)
if( !is.null(opt$help) || is.null(opt$metaData) || is.null(opt$qvalue) || is.null(opt$output) ){
  cat(paste(getopt(spec=spec, usage = T), "\n"))
  quit()
}



if(file_ext(opt$metaData)=='csv')
{
  pheno_data <- read.csv(opt$metaData, header = TRUE)
}else if(file_ext(opt$metaData)=='txt')
{
  pheno_data <- read.table(opt$metaData, header = TRUE)
}else
{
  print('Only csv and txt files are supported!')
  quit()    
}

bg_chrX <- ballgown(dataDir = "ballgown",
                    samplePattern = "*",
                    pData = pheno_data)

bg_chrX_filt <- subset(bg_chrX, "rowVars(texpr(bg_chrX)) >1", genomesubset=TRUE)

results_transcripts <- stattest(bg_chrX_filt,
                                feature="transcript",
                                covariate="dex",
                                getFC=TRUE, meas="FPKM")

results_transcripts <- data.frame(geneNames = geneNames(bg_chrX_filt),
                                  geneIDs = geneIDs(bg_chrX_filt),
                                  results_transcripts)

results_de=results_transcripts %>% filter(qval < opt$qvalue)
gene_to_plot=as.character(results_de[results_de$geneNames==opt$Gene_to_plot,'geneIDs'])



write.csv(results_de,file='Significant_differential_transcripts.csv')
if(!require("ggplot2")) install.packages("ggplot2")
if(!require("cowplot")) install.packages("cowplot")
library(ggplot2)
library(cowplot)

results_transcripts$mean <- rowMeans(texpr(bg_chrX_filt))

pdf(file=paste(opt$output,'/volcano_plot.pdf',sep=''))
ggplot(results_transcripts, aes(log2(mean), log2(fc), colour = qval< opt$qvalue)) +
  scale_color_manual(values=c("#999999", "#FF0000")) +
  geom_point() +
  geom_hline(yintercept=0)
dev.off()

pdf(file=paste(opt$output,'/Transcipt_plot_for_one_gene.pdf',sep=''))
plotMeans(gene_to_plot, bg_chrX_filt,groupvar="dex",legend=FALSE)
dev.off()




