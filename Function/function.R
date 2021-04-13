#!/usr/bin/Rscript
if(!require("getopt")) install.packages("getopt")
library(getopt)


spec <- matrix(
  c("input_significant_gene_file",  "i", 1, "character", "Input the significant gene file that is produced by differential gene analysis",
    "gene_type", "g", 1, "character",  "The ID type of genes, SYMBOL, ENTREZID, ENSEMBL,ENSEMBLTRANS are supported",
    "number_of_significant_genes", "n", 2, "numeric",  "the number of significant genes you want to analyze in you list",
    "output",  "o", 2, "character",  "The directory of output",
    "function_type","f",1,"character","This is function type you want choose for function enrichment, GO and KEGG are supported",
    "organism", "m", 1,"character","This is the organism of your samples, right now only human and mouse are supported ",
    "help",   "h", 0, "double",  "This is Help,input_significant_gene_file,gene_type,function_type,organism are mandatory parameters"
  ),
  byrow=TRUE, ncol=5)

# 使用getopt方法
opt <- getopt(spec=spec)
#print(opt$first)
if( !is.null(opt$help) || is.null(opt$input_significant_gene_file) || is.null(opt$gene_type) || is.null(opt$organism) || is.null(opt$function_type)){
  # ... 这里你也可以自定义一些东放在里面
  cat(paste(getopt(spec=spec, usage = T), "\n"))
  quit()
}

if(!require("clusterProfiler")) BiocManager::install("clusterProfiler")

if(!require("DOSE")) BiocManager::install("DOSE")
if(!require("Rgraphviz")) BiocManager::install("Rgraphviz")

library(clusterProfiler)
library(DOSE)
library(Rgraphviz)

x=read.csv(opt$input_significant_gene_file)[1:opt$number_of_significant_genes,]
print(x)

print(x)

if(opt$organism=='Human')
{
  if(!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db")
  library(org.Hs.eg.db)
  eg = bitr(x, fromType=opt$gene_type, toType="ENTREZID",  OrgDb="org.Hs.eg.db")
  
  if(opt$function_type=='GO')
  {
    ego=enrichGO(eg$ENTREZID, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                 qvalueCutoff = 0.2,keyType = 'ENTREZID')
  }  
  if(opt$function_type=='KEGG')
  {
    ego=enrichKEGG(eg$ENTREZID, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', 
               minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
  }
}

if(opt$organism=='Mouse')
{
  if(!require("org.Mm.eg.db")) BiocManager::install("org.Mm.eg.db")
  library(org.Mm.eg.db)
  eg = bitr(x, fromType=opt$gene_type, toType="ENTREZID",  OrgDb="org.Mm.eg.db")
  
  if(opt$function_type=='GO')
  {
    ego=enrichGO(eg$ENTREZID, OrgDb = org.Mm.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                 qvalueCutoff = 0.2,keyType = 'ENTREZID')
  }  
  if(opt$function_type=='KEGG')
  {
    ego=enrichKEGG(eg$ENTREZID, organism = 'mouse', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', 
                   minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
  }
}



pdf(file=paste(opt$output,'/function_annotation.pdf',sep=''))
myplot=barplot(ego)
print(myplot)
dev.off()