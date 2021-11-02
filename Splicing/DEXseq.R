
if(!require("getopt")) install.packages("getopt")
library(getopt)
library(tools)

spec <- matrix(
  c("count_dir",'c',1,"character","The directory of count matrix",
    "metaData",  "m", 1, "character", "meta-data of the samples, introduce their conditions",
    "pvalue", "p", 1,"double","This is the highest p_value for detetcting significant differential exons",
    "output","o",1,"character","This is output directory",
    "gene_to_plot","g",1,"character","Plot exon expression of specific gene",
    "help",   "h", 0, "double",  "This is Help,input_matrix,metadata,pvalue,logfoldchange are mandatory parameters"
  ),
  byrow=TRUE, ncol=5)

# 使用getopt方法
opt <- getopt(spec=spec)
if( !is.null(opt$help) || is.null(opt$count_dir) || is.null(opt$metaData) || is.null(opt$pvalue) || is.null(opt$output)){
  # ... 这里你也可以自定义一些东放在里面
  cat(paste(getopt(spec=spec, usage = T), "\n"))
  quit()
}

inDir = opt$count_dir
countFiles = list.files(inDir, pattern=".txt$", full.names=TRUE)
flattenedFile = list.files(inDir, pattern="gff$", full.names=TRUE)

if(file_ext(opt$metaData)=='csv')
{
  sampleTable <- read.csv(opt$metaData, header = TRUE,row.names = 1)
}else if(file_ext(opt$metaData)=='txt')
{
  sampleTable <- read.table(opt$metaData, header = TRUE,row.names = 1)
}else
{
  print('Only csv and txt files are supported!')
  quit()    
}


if(!require("DEXSeq")) BiocManager::install("DEXSeq")

library( "DEXSeq" )

dxd = DEXSeqDataSetFromHTSeq(
  countFiles,
  sampleData=sampleTable,
  design= ~ sample + exon + dex:exon,
  flattenedfile=flattenedFile )

dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd)
plotDispEsts( dxd )
dxd = testForDEU( dxd)
dxd = estimateExonFoldChanges( dxd, fitExpToVar="dex")

dxr1 = DEXSeqResults( dxd )

pdf(file=paste(opt$output,'/volcano_plot.pdf',sep=''))
plotMA( dxr1, cex=0.8 )
dev.off()

dxr2=data.frame(dxr1)
dxr3=na.omit(dxr2[dxr2$padj < opt$pvalue,])

write.csv(dxr3,file=paste(opt$output,'/Significant_differential_exons.csv',sep=''))

pdf(file=paste(opt$output,'/exon_plot_for_one_gene.pdf',sep=''))
plotDEXSeq( dxr1, opt$gene_to_plot, fitExpToVar="dex", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2,splicing=TRUE )
dev.off()
