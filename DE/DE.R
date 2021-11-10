
if(!require("getopt")) install.packages("getopt")
library(getopt)
library(tools)

spec <- matrix(
  c("method",'e',1,"character","Methods for DE analysis, DEseq2, edgeR, limma are supported",
    "input_matrix",  "i", 1, "character", "Input matrix file",
    "metadata", "m", 1, "character",  "Metadata that introduce the condition of samples",
    "output",  "o", 2, "character",  "The directory of output",
    "pvalue", "p", 1,"double","This is the highest p_value for detetcting significant differential genes",
    "logfoldchange","f",1,"double","This is the bound of logfoldchange(absolute value) for detetcting significant differential genes",
    "help",   "h", 0, "double",  "This is Help,input_matrix,metadata,pvalue,logfoldchange are mandatory parameters"
  ),
  byrow=TRUE, ncol=5)

# 使用getopt方法
opt <- getopt(spec=spec)
if( !is.null(opt$help) || is.null(opt$input_matrix) || is.null(opt$metadata) || is.null(opt$pvalue) || is.null(opt$logfoldchange) ){
  # ... 这里你也可以自定义一些东放在里面
  cat(paste(getopt(spec=spec, usage = T), "\n"))
  quit()
}

if(!is.null(opt$ouput) & (!file.exists(opt$output))){dir.create(opt$output)}

if(opt$method=='DESeq2')
{
  if(!require("DESeq2")) BiocManager::install("DESeq2")
  if(!require("ggplot2")) BiocManager::install("ggplot2")
  
  library(DESeq2)
  library(ggplot2)
  
  if(file_ext(opt$input_matrix)=='csv')
  {
    countData <- read.csv(opt$input_matrix, header = TRUE)
    print('This is the part of count data')
    head(countData)
  }else if(file_ext(opt$input_matrix)=='txt')
  {
    countData <- read.table(opt$input_matrix, header = TRUE)
    print('This is the part of count data')
    head(countData)
  }else
  {
      print('Only csv and txt files are supported!')
      quit()    
    }
    

  
  #metaDataName <- "http://bioconnector.org/workshops/data/airway_metadata.csv"
  #download.file(metaDataName, destfile = "airway_metadata.csv", method = "auto")
  if(file_ext(opt$metadata)=='csv')
  {
    metaData <- read.csv(opt$metadata, header = TRUE)
    print('This is the metadata')
    print(metaData)
  }else if(file_ext(opt$metadata)=='txt')
  {
    metaData <- read.table(opt$metadata, header = TRUE)
    print('This is the metadata')
    print(metaData)
  }else
  {
    print('Only csv and txt files are supported!')
    quit()  
    }

  
  dds <- DESeqDataSetFromMatrix(countData=countData, 
                                colData=metaData, 
                                design=~dex, tidy = TRUE)
  
  dds <- DESeq(dds)
  
  res <- results(dds)
  
  ok=summary(res)
  print(ok)
  res=na.omit(res)
  
  res <- res[order(res$pvalue),]
  pdf(file=paste(opt$output,'/volcano_plot.pdf',sep=''))
  #reset par
  par(mfrow=c(1,1))
  # Make a basic volcano plot
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot"))
  
  # Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
  with(subset(res, pvalue<(opt$pvalue) ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
  with(subset(res, pvalue<(opt$pvalue) & abs(log2FoldChange)>(opt$logfoldchange)), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
  
  dev.off()
  vsdata <- vst(dds, blind=FALSE)
  pdf(file=paste(opt$output,'/pca_plot.pdf',sep=''))
  plotPCA(vsdata, intgroup="dex") 
  dev.off()
  
  res_interest=res[(res$pvalue)<opt$pvalue & (res$log2FoldChange)>(opt$logfoldchange),]
  res_interest <- res_interest[order(res_interest$pvalue),]
  write.csv(res_interest,paste(opt$output,'/de_gene_foldchange_above_0.csv',sep=''))
  
  
  res_interest=res[(res$pvalue)<opt$pvalue & res$log2FoldChange<(-(opt$logfoldchange)),]
  res_interest <- res_interest[order(res_interest$pvalue),]
  write.csv(res_interest,paste(opt$output,'/de_gene_foldchange_below_0.csv',sep=''))
  
}


if(opt$method=='edgeR')
{
  if(!require("edgeR")) BiocManager::install("edgeR")
  if(!require("ggplot2")) BiocManager::install("ggplot2")
  library("edgeR")
  library(ggplot2)
  
  
  if(file_ext(opt$input_matrix)=='csv')
  {
    data_clean <- read.csv(opt$input_matrix, header = TRUE,row.names=1)
    print('This is the part of count data')
    head(data_clean)
  }else if(file_ext(opt$input_matrix)=='txt')
  {
    data_clean <- read.table(opt$input_matrix, header = TRUE, row.names=1)
    print('This is the part of count data')
    head(data_clean)
  }else
  {
    print('Only csv and txt files are supported!')
    quit()    
  }
  
  
  
  #metaDataName <- "http://bioconnector.org/workshops/data/airway_metadata.csv"
  #download.file(metaDataName, destfile = "airway_metadata.csv", method = "auto")
  if(file_ext(opt$metadata)=='csv')
  {
    metaData <- read.csv(opt$metadata, header = TRUE)
    print('This is the metadata')
    print(metaData)
  }else if(file_ext(opt$metadata)=='txt')
  {
    metaData <- read.table(opt$metadata, header = TRUE)
    print('This is the metadata')
    print(metaData)
  }else
  {
    print('Only csv and txt files are supported!')
    quit()  
  }
  
  

  cpm_log <- cpm(data_clean, log = TRUE)
  median_log2_cpm <- apply(cpm_log, 1, median)
  expr_cutoff <- -1
  
  data_clean <- data_clean[median_log2_cpm > expr_cutoff, ]
  
  
  group <-metaData[,2]
  group
  y <- DGEList(counts = data_clean, group = group)
  #y
  
  y <- calcNormFactors(y)
  
  y <- estimateDisp(y)
  print('Biological coefficient of varitation')
  sqrt(y$common.dispersion)
  
  et <- exactTest(y)
  results_edgeR <- topTags(et, n = nrow(data_clean), sort.by = "none")
  print('!!!This is so important, this shows the results for comparison mode')
  print(results_edgeR$comparison)
  
  res=results_edgeR$table
  res=na.omit(res)
  res <- res[order(res$PValue),]
  pdf(file=paste(opt$output,'/volcano_plot.pdf',sep=''))
  #reset par
  par(mfrow=c(1,1))
  # Make a basic volcano plot
  with(res, plot(logFC, -log10(PValue), pch=20, main="Volcano plot" ))
  
  # Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
  with(subset(res, PValue<(opt$pvalue) ), points(logFC, -log10(PValue), pch=20, col="blue"))
  with(subset(res, PValue<(opt$pvalue) & abs(logFC)>(opt$logfoldchange)), points(logFC, -log10(PValue), pch=20, col="red"))
  
  dev.off()
  
  cpm_log <- cpm(data_clean, log = TRUE)
  pdf(file=paste(opt$output,'/sample_corelation_plot.pdf',sep=''),width=10, height=10)
  heatmap(cor(cpm_log))
  dev.off()
  
  res_interest=res[(res$PValue)<opt$pvalue & (res$logFC)>(opt$logfoldchange),]
  res_interest <- res_interest[order(res_interest$PValue),]
  write.csv(res_interest,paste(opt$output,'/de_gene_foldchange_above_0.csv',sep=''))
  
  
  res_interest=res[(res$PValue)<opt$pvalue & res$logFC<(-(opt$logfoldchange)),]
  res_interest <- res_interest[order(res_interest$PValue),]
  write.csv(res_interest,paste(opt$output,'/de_gene_foldchange_below_0.csv',sep=''))
  
}

if(opt$method=='limma')
{
  if(!require("edgeR")) BiocManager::install("edgeR")
  if(!require("ggplot2")) BiocManager::install("ggplot2")
  
  #library(getopt)
  library("edgeR")
  library(ggplot2)
  if(file_ext(opt$input_matrix)=='csv')
  {
    data_clean <- read.csv(opt$input_matrix, header = TRUE,row.names=1)
    print('This is the part of count data')
    head(data_clean)
  }else if(file_ext(opt$input_matrix)=='txt')
  {
    data_clean <- read.table(opt$input_matrix, header = TRUE, row.names=1)
    print('This is the part of count data')
    head(data_clean)
  }else
  {
    print('Only csv and txt files are supported!')
    quit()    
  }
  
  
  
  #metaDataName <- "http://bioconnector.org/workshops/data/airway_metadata.csv"
  #download.file(metaDataName, destfile = "airway_metadata.csv", method = "auto")
  if(file_ext(opt$metadata)=='csv')
  {
    metaData <- read.csv(opt$metadata, header = TRUE)
    print('This is the metadata')
    print(metaData)
  }else if(file_ext(opt$metadata)=='txt')
  {
    metaData <- read.table(opt$metadata, header = TRUE)
    print('This is the metadata')
    print(metaData)
  }else
  {
    print('Only csv and txt files are supported!')
    quit()  
  }
  

  cpm_log <- cpm(data_clean, log = TRUE)
  median_log2_cpm <- apply(cpm_log, 1, median)
  expr_cutoff <- -1
  
  data_clean <- data_clean[median_log2_cpm > expr_cutoff, ]
  
  
  group <-metaData[,2]
  group
  y <- DGEList(counts = data_clean, group = group)
  #y
  
  y <- calcNormFactors(y)
  
  mm <- model.matrix(~0 + group)
  y <- voom(y, mm, plot = F)
  
  fit <- lmFit(y, mm)
  condition=paste(colnames(coef(fit))[2],colnames(coef(fit))[1],sep='-')
  contr <- makeContrasts(condition, levels = colnames(coef(fit)))
  print('!!!This is so important, this shows the comparison mode')
  print(contr)
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, n = Inf)
  head(top.table, 20)
  
  res=top.table
  res=na.omit(res)
  res <- res[order(res$adj.P.Val),]
  pdf(file=paste(opt$output,'/volcano_plot.pdf',sep=''))
  #reset par
  par(mfrow=c(1,1))
  # Make a basic volcano plot
  with(res, plot(logFC, -log10(adj.P.Val), pch=20, main="Volcano plot"))
  
  # Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
  with(subset(res,adj.P.Val<(opt$pvalue) ), points(logFC, -log10(adj.P.Val), pch=20, col="blue"))
  with(subset(res, adj.P.Val<(opt$pvalue) & abs(logFC)>(opt$logfoldchange)), points(logFC, -log10(adj.P.Val), pch=20, col="red"))
  
  dev.off()
  
  cpm_log <- cpm(data_clean, log = TRUE)
  pdf(file=paste(opt$output,'/sample_corelation_plot.pdf',sep=''))
  heatmap(cor(cpm_log))
  dev.off()
  
  res_interest=res[(res$adj.P.Val)<opt$pvalue & (res$logFC)>(opt$logfoldchange),]
  res_interest <- res_interest[order(res_interest$adj.P.Val),]
  write.csv(res_interest,paste(opt$output,'/de_gene_foldchange_above_0.csv',sep=''))
  
  
  res_interest=res[(res$adj.P.Val)<opt$pvalue & res$logFC<(-(opt$logfoldchange)),]
  res_interest <- res_interest[order(res_interest$adj.P.Val),]
  write.csv(res_interest,paste(opt$output,'/de_gene_foldchange_below_0.csv',sep=''))

}



if(opt$method=='T-test')
{

  
  if(file_ext(opt$input_matrix)=='csv')
  {
    countData <- read.csv(opt$input_matrix, header = TRUE,row.names=1)
    print('This is the part of count data')
    head(countData)
  }else if(file_ext(opt$input_matrix)=='txt')
  {
    countData <- read.table(opt$input_matrix, header = TRUE,row.names=1)
    print('This is the part of count data')
    head(countData)
  }else
  {
      print('Only csv and txt files are supported!')
      quit()    
    }
    

  if(file_ext(opt$metadata)=='csv')
  {
    metaData <- read.csv(opt$metadata, header = TRUE)
    print('This is the metadata')
    print(metaData)
  }else if(file_ext(opt$metadata)=='txt')
  {
    metaData <- read.table(opt$metadata, header = TRUE)
    print('This is the metadata')
    print(metaData)
  }else
  {
    print('Only csv and txt files are supported!')
    quit()  
    }  

normal=metaData[metaData['dex']=="normal"]
tumor=metaData[metaData['dex']!="normal"]


Pvalue<-c(rep(0,nrow(countData))) 
log2_FC<-c(rep(0,nrow(countData))) 
fdr=c(rep(0,nrow(countData)))


for(i in 1:nrow(countData)){
 if(sd(countData[i,colnames(countData) %in% normal])==0&&sd(countData[i,colnames(countData) %in% tumor])==0){
 Pvalue[i] <-NA
 log2_FC[i]<-NA
fdr[i]=NA
 }else{
 y=t.test(as.numeric(countData[i,colnames(countData) %in% normal]),as.numeric(countData[i,colnames(countData) %in% tumor]))
 Pvalue[i]<-y$p.value
 log2_FC[i]<-log2((mean(as.numeric(countData[i,colnames(countData) %in% normal]))+0.001)/(mean(as.numeric(countData[i,colnames(countData) %in% tumor]))+0.001)) 
fdr[i]=p.adjust(Pvalue[i], "BH") 
}
}

#fdr[i]=p.adjust(Pvalue[i], "BH") 

res<-data.frame(cbind(countData,log2_FC,Pvalue))
#res=as.data.frame(res)

#print(head(res))
res=na.omit(res) 
print("log2_FC, Normal/Case")
#print(head(res))

  res_interest=res[(res$Pvalue)<opt$pvalue & (res$log2_FC)>(opt$logfoldchange),]
  res_interest <- res_interest[order(res_interest$Pvalue),]
  write.csv(res_interest,paste(opt$output,'/de_gene_foldchange_above_0.csv',sep=''))
  
  
  res_interest=res[(res$Pvalue)<opt$pvalue & res$log2_FC<(-(opt$logfoldchange)),]
  res_interest <- res_interest[order(res_interest$Pvalue),]
  write.csv(res_interest,paste(opt$output,'/de_gene_foldchange_below_0.csv',sep=''))

  
  
  res <- res[order(res$Pvalue),]
  pdf(file=paste(opt$output,'/volcano_plot.pdf',sep=''))
  #reset par
  par(mfrow=c(1,1))
  # Make a basic volcano plot
  with(res, plot(log2_FC, -log10(Pvalue), pch=20, main="Volcano plot"))
  
  # Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
  with(subset(res, Pvalue<(opt$pvalue) ), points(log2_FC, -log10(Pvalue), pch=20, col="blue"))
  with(subset(res, Pvalue<(opt$pvalue) & abs(log2_FC)>(opt$logfoldchange)), points(log2_FC, -log10(Pvalue), pch=20, col="red"))
  
  dev.off()

  

  
}

