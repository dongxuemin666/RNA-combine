#!/usr/bin/Rscript
if(!require("getopt")) install.packages("getopt")
if(!require("pheatmap")) install.packages("pheatmap")

library(getopt)
library(tools)
library(pheatmap)




spec <- matrix(
  c("method",'e',1,"character","Methods for relation calculation, Pearson, Spearman, PCA-PMI(Juan Zhao et al. PNAS, 2016) are supported",
    "input_matrix",  "i", 1, "character", "Input matrix file",
    "output",  "o", 2, "character",  "The directory of output",
    "help",   "h", 0, "double",  "This is help"
  ),
  byrow=TRUE, ncol=5)

# 使用getopt方法
opt <- getopt(spec=spec)
if( !is.null(opt$help) || is.null(opt$input_matrix) || is.null(opt$method) ){
  # ... 这里你也可以自定义一些东放在里面
  cat(paste(getopt(spec=spec, usage = T), "\n"))
  quit()
}







# ---------------------------PCA-CMI------------------------------
#  This is a R code for PCA-CMI
#  inferring gene regulatory networks from gene expression data by path consitency algorithm
#  based on conditional mutual information.
#
# Input:
#   'dat' is expression of variable,in which row is varible and column is the sample;
#   'lamda' is the parameter decide the dependence;
#   'order0' is the parameter to end the program when order=order0;
#
# Output:
#   'G' is the 0-1 network or graph after pc algorith
#   'Gval' is the network with strenthness of dependence;
#   'order' is the order of the pc algorithm, here is equal to order0;
#
#  If nargin==2,the algorithm will be terminated untill there is no higher order networks.
#
#  Juan Zhao et al. PNAS, 2016







# compute conditional mutual information of x and y
cmi <- function(v1, v2, nargin = 2, vcs = NULL) {
      if (nargin == 2) {


 
            c1 = var(v1)
            c2 = var(v2)
            c3 = det(cov(t(rbind(v1, v2))))
            cmiv = 0.5 * log(c1 * c2 / c3)
	


      } else if (nargin == 3) {





        #matlab

       #n1 = size(vcs,1);
       n1=nrow(vcs)
       
       
       n = n1 +2

       Cov1 = var(v1)
       Cov2 = var(v2)
       #Covm = cov([v1;v2;vcs]');
       Covm=cov(t(rbind(v1,v2,vcs)))
       #print(Covm)
       #Covm1 = cov(vcs');
       if (is.null(dim(vcs))){
       Covm1 = var(vcs)
       InvCovm1=1/Covm1
       #print(InvCovm1)
       }
       else{
       Covm1=cov(as.matrix(vcs))
       InvCovm1=solve(Covm1)
       print(InvCovm1)}
       #Covm2 =cov([v1;vcs]');
       Covm2=cov(t(rbind(v1,vcs)))
       #Covm3 = cov([v2;vcs]');
       Covm3=cov(t(rbind(v2,vcs)))

       InvCov1 = 1/Cov1
       InvCov2 = 1/Cov2
       #InvCovm = inv(Covm);
  
       #InvCovm1 = inv(Covm1);
       InvCovm = solve(Covm)

       #InvCovm2 = inv(Covm2);
       InvCovm2 = solve(Covm2)
       #print(InvCovm2)
       #InvCovm3 = inv(Covm3);
       InvCovm3 = solve(Covm3)

       #C11 = -InvCovm(1,2)*(1/(InvCovm(2,2)-InvCovm3(1,1)+InvCov2))*InvCovm(1,2)+InvCovm(1,1);
       C11 = -InvCovm[1,2]*(1/(InvCovm[2,2]-InvCovm3[1,1]+InvCov2))*InvCovm[1,2]+InvCovm[1,1]
       C12 = 0

       C22 = -InvCovm[1,2]*(1/(InvCovm[1,1]-InvCovm2[1,1]+InvCov1))*InvCovm[1,2]+InvCovm[2,2]


       #C33 = -(InvCovm(2,3:2+n1)-InvCovm3(1,2:1+n1))'*(1/(InvCovm(2,2)-InvCovm3(1,1)+InvCov2))*(InvCovm(2,3:2+n1)-InvCovm3(1,2:1+n1))-(InvCovm(1,3:2+n1)-InvCovm2(1,2:1+n1))'*(1/(InvCovm(1,1)-InvCovm2(1,1)+InvCov1))*(InvCovm(1,3:2+n1)-InvCovm2(1,2:1+n1))+(InvCovm(3:2+n1,3:2+n1)-InvCovm3(2:1+n1,2:1+n1))+(InvCovm(3:2+n1,3:2+n1)-InvCovm2(2:1+n1,2:1+n1))+InvCovm1;

       
       if (is.null(dim(vcs))) {

         
         C13 =-InvCovm[1,2]*(1/(InvCovm[2,2]-InvCovm3[1,1]+InvCov2))*(InvCovm[2,3]-InvCovm3[1,2])+(InvCovm[1,3])
         C23 = -InvCovm[1,2]*(1/(InvCovm[1,1]-InvCovm2[1,1]+InvCov1))*(InvCovm[1,3]-InvCovm2[1,2])+(InvCovm[2,3])
         C33 = -(InvCovm[2,3]-InvCovm3[1,2])*(1/(InvCovm[2,2]-InvCovm3[1,1]+InvCov2))*(InvCovm[2,3]-InvCovm3[1,2])-(InvCovm[1,3]-InvCovm2[1,2])*(1/(InvCovm[1,1]-InvCovm2[1,1]+InvCov1))*(InvCovm[1,3]-InvCovm2[1,2])+(InvCovm[3,3]-InvCovm3[2,2])+(InvCovm[3,3]-InvCovm2[2,2])+InvCovm1
         InvC=rbind(c(C11,C12,C13),c(C12,C22,C23),c(C13,C23,C33))
         #print(InvC)
        #C0 = (det(Covm)*det(Covm1)/(det(Covm2)*det(Covm3)))*Cov1*Cov2*(InvCovm(2,2)-InvCovm3(1,1)+InvCov2)*(InvCovm(1,1)-InvCovm2(1,1)+InvCov1);
         C0 = (det(Covm)*Covm1/(det(Covm2)*det(Covm3)))*Cov1*Cov2*(InvCovm[2,2]-InvCovm3[1,1]+InvCov2)*(InvCovm[1,1]-InvCovm2[1,1]+InvCov1)
        #cmiv = 0.5 * (trace(InvC*Covm)+log(C0)-n); 
         #print(C0)
         cmiv=0.5*(sum(diag(InvC%*%Covm))+log(C0)-3)
         #print(cmiv)
         }
       else{
       
       C13 = t(as.matrix(-InvCovm[1,2]*(1/(InvCovm[2,2]-InvCovm3[1,1]+InvCov2))*(InvCovm[2,3:2+n1]-InvCovm3[1,2:1+n1])+(InvCovm[1,3:2+n1])))
       C23 = t(as.matrix(-InvCovm[1,2]*(1/(InvCovm[1,1]-InvCovm2[1,1]+InvCov1))*(InvCovm[1,3:2+n1]-InvCovm2[1,2:1+n1])+(InvCovm[2,3:2+n1])))
       
       C33 = -(as.matrix(InvCovm[2,3:2+n1]-InvCovm3[1,2:1+n1])*(1/(InvCovm[2,2]-InvCovm3[1,1]+InvCov2)))%*%(t(InvCovm[2,3:2+n1]-InvCovm3[1,2:1+n1]))-(as.matrix(InvCovm[1,3:2+n1]-InvCovm2[1,2:1+n1])*(1/(InvCovm[1,1]-InvCovm2[1,1]+InvCov1)))%*%(t(InvCovm[1,3:2+n1]-InvCovm2[1,2:1+n1]))+(InvCovm[3:2+n1,3:2+n1]-InvCovm3[2:1+n1,2:1+n1])+(InvCovm[3:2+n1,3:2+n1]-InvCovm2[2:1+n1,2:1+n1])+InvCovm1
       InvC=rbind(c(C11,C12,C13),c(C12,C22,C23),cbind(t(C13),t(C23),C33))
              #C0 = (det(Covm)*det(Covm1)/(det(Covm2)*det(Covm3)))*Cov1*Cov2*(InvCovm(2,2)-InvCovm3(1,1)+InvCov2)*(InvCovm(1,1)-InvCovm2(1,1)+InvCov1);
       C0 = (det(Covm)*det(Covm1)/(det(Covm2)*det(Covm3)))*Cov1*Cov2*(InvCovm[2,2]-InvCovm3[1,1]+InvCov2)*(InvCovm[1,1]-InvCovm2[1,1]+InvCov1)
       #cmiv = 0.5 * (trace(InvC*Covm)+log(C0)-n); 
       cmiv=0.5*(sum(diag(InvC%*%Covm))+log(C0)-n)
       
       }

       #InvC = [[C11,C12,C13];[C12,C22,C23];[[C13',C23'],C33]];
       # C = inv(InvC);  
       




            
      }
      
      #if(nargin==2)
      #cat('\n values: ',c1,c2,c3,cmiv,sep=' ')
      #if(nargin==3)
      #cat('\n values: ',c1,c2,c3,c4,cmiv,sep=' ')
      
      if (cmiv == Inf) {
            cmiv = 0
      }
      return(cmiv)
}

  
# edgereduce is pca_cmi
#function [G,Gval,t]=edgereduce(G,Gval,order,data,t,lamda)
edgereduce <- function(G, Gval, order, dat, t, lambda) {
      if (order == 0) {
            for (i in 1:(nrow(G) - 1)) {
                  for (j in (i + 1):nrow(G)) {
                        if (G[i, j] != 0) {
                              cmiv = cmi(dat[i,], dat[j,])
                              Gval[i, j] = cmiv
                              Gval[j, i] = cmiv
                              
                              if (cmiv < lambda) {
                                    G[i, j] = 0
                                    G[j, i] = 0
                              }
                        }
                  }
            }
            t = t + 1
            
      } else {
            for (i in 1:(nrow(G) - 1)) {
                  for (j in (i + 1):nrow(G)) {
                        if (G[i, j] != 0) {
                              adj = numeric()
                              
                              for (k in 1:nrow(G)) {
                                    if (G[i, k] != 0 && G[j, k] != 0) {
                                          adj = c(adj, k)
                                    }
                              }
                              if (length(adj) >= order) {
                                    if (length(adj) == 1) {
                                          combntnslist = adj
                                          combntnsrow = 1
                                    } else {
                                          combntnslist = t(combn(adj, order))
                                          combntnsrow = nrow(combntnslist)
                                    }
                                    cmiv = 0
                                    v1 = dat[i,]
                                    v2 = dat[j,]
                                    
                                    if (combntnsrow == 1) {
                                          vcs = dat[combntnslist,]
                                          cmiv = cmi(v1, v2, 3, vcs)
                                    } else {
                                          for (k in 1:combntnsrow) {
                                                vcs = dat[combntnslist[k,],]
                                                
                                                a = cmi(v1,
                                                        v2,
                                                        3,
                                                        vcs)
                                                
                                                cmiv = max(cmiv, a)
                                                
                                          }
                                    }
                                    Gval[i, j] = cmiv
                                    Gval[j, i] = cmiv
                                    
                                    if (cmiv < lambda) {
                                          G[i, j] = 0
                                          G[j, i] = 0
                                    }
                                    t = t + 1
                              }
                        }
                  }
            }
      }
      return(list(G = G, Gval = Gval, t = t))
}

pca_cmi <- function(dat, lambda, order0 = 0, order0_key = 'off') {
    # function [G,Gval,order]=pca_cmi(data,lamda,order0)
    library(Matrix)
    
    n_gene = nrow (dat)
    G = matrix(1, n_gene, n_gene)
    G = t(as.matrix(tril(G, -1)))
    
    G = G + t(G)
    Gval = G
    order = -1
    t = 0
    
    while (t == 0) {
        order = order + 1
        
        if (order0_key == 'on') {
            if (order > order0) {
                G = t(as.matrix(tril(G, -1)))
                
                Gval = t(as.matrix(tril(G, -1)))
                
                order = order - 1
                # The value of order is the last order of pc algorith
                return(list(
                    G = G,
                    Gval = Gval,
                    order = order
                ))
            }
        }
        
        res = edgereduce(G, Gval, order, dat, t, lambda)
        
        G = res$G
        Gval = res$Gval
        t = res$t
        
        if (t == 0) {
            cat('\n No edge is reduce! Algorithm  finished!')
            break
        } else {
            t = 0
        }
    }
    

    
    order = order - 1
    # The value of order is the last order of pc algorith
    return(list(G = G,
                Gval = Gval,
                order = order))
}





if(file_ext(opt$input_matrix)=='csv'){
x0 <- read.csv(opt$input_matrix, header = TRUE,row.names=1)
print('This is the part of count data')
head(x0)
  } else if(file_ext(opt$input_matrix)=='txt'){
    x0 <- read.table(opt$input_matrix, header = TRUE,row.names=1)
    print('This is the part of count data')
    head(x0)
  } else{
      print('Only csv and txt files are supported!')
      quit()    
    }
  
    
    

if(opt$method=='Spearman'){  
   pdf(file=paste(opt$output,'/Spearman_plot.pdf',sep=''))
   plot=pheatmap(cor(t(x0),method="spearman"),cluster_rows = FALSE,cluster_columns = FALSE)
   print(plot)
    dev.off()
       
   write.csv(cor(t(x0),method="spearman"),file=paste(opt$output,'/Spearman_relations.csv',sep=''))

}

if( opt$method=='Pearson'){  
   pdf(file=paste(opt$output,'/Pearson_plot.pdf',sep=''))
   plot=pheatmap(cor(t(x0),method="pearson"), cluster_rows = FALSE,cluster_columns = FALSE)
   print(plot)
    dev.off()
       
   write.csv(cor(t(x0),method="pearson"),file=paste(opt$output,'/Pearson_relations.csv',sep=''))

}
    
if( opt$method=='PCA-PMI'){  
    
    x2=as.matrix(x0)
    pan=pca_cmi(x2, 0.4)


    row.names(pan$G)=rownames(x2)
    row.names(pan$Gval)=rownames(x2)
    colnames(pan$G)=rownames(x2)
    colnames(pan$Gval)=rownames(x2)
    
    
    pdf(file=paste(opt$output,'/PMI_plot.pdf',sep=''))
    plot=pheatmap(pan$Gval,cluster_rows = FALSE,cluster_columns = FALSE)
    print(plot)
    dev.off()
       
    write.csv(pan$Gval,file=paste(opt$output,'/PMI_relations.csv',sep=''))

}    
























