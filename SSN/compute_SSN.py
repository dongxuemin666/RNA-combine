import os
import scipy.stats as stat
import math
import time
import sys



#  Xiaoping Liu et al., 2016
begin=time.asctime()
print(begin)

param={}
for i in range(1,len(sys.argv)):
    t=sys.argv[i].split("=")
    param[t[0].lower()]=t[1]
    
help_msg="""
usage: python construct_single_network.py -pvalue=threshold_p_value | -threshold=threshold_value -background=background_network_file -ref=reference_sample_file  -sample=sample_data_file -out=results_output_fold
Options and arguments:
-pvalue : set the threshold of p-value [0..1], if the -pvalue set 1, all edges will be outputted to the SSN
-background : background network to calculate the deltaPCC of edges based on the network
-ref : the expression profile of reference samples
-sample : the expression profile for the sample to be constructed the SSN
-out : the directory to store the SSN
"""
   
if "-help" in param.keys() or "-h" in param.keys():
    print(help_msg)

if "-ref" not in param.keys() or "-sample" not in param.keys() or "-background" not in param.keys():
    print("Parameter missing!")
    print(help_msg)
    exit()

disease={}
normal={}

name_c=[]

reference_file=param["-ref"]
sample_file=param["-sample"]
background=param["-background"]


# read the giving pvalue in robust way

if "-pvalue" not in param.keys() :
    print("The parameter -pvalue need be set!")
    exit()

p_value=1
if "-pvalue" in param.keys():
    p_value=float(param["-pvalue"])
    if p_value < 0 or p_value > 1:
        print("Please set the correct threshold of p-value in -pvalue [0..1]")
        exit()

if "-out" not in param.keys():    
    fold="."
else:
    fold=param["-out"]





if not os.path.exists(fold):
    os.mkdir(fold)
# constrcut gene expression matrix for background
f=open(reference_file)
flag=0
n=0
for p in f: # read line by line
    flag+=1 #how many genes
    t=p.split(',')
    # splicing
    if flag==1:
        n=len(t)#how many samples
        name_c=t[1:]
        continue
    normal[t[0]]=[float(t[i]) for i in range(1,n)]
    #construct a list of expression data for background
f.close()



name_d=[]
# constrcut gene expression matrix for disease sample
f=open(sample_file)
flag=0
n=0
for p in f:
    flag+=1
    t=p.split(',')
    if flag==1:
        n=len(t)
        #here numer of samples is n-1
        name_d=t[1:]
        continue #next round of for
    disease[t[0]]=[float(t[i]) for i in range(1,n)]
    #construct a list of expression data for disease
f.close()




if "-pvalue" in param.keys():
    p_value=float(param["-pvalue"])
    for i in range(len(name_d)):
        print("sample",name_d[i])
        fw=open(fold+os.sep+"ssn_"+name_d[i]+".txt","w")
        fw.write("Gene1\tGene2\tbackground_PCC\tditurb_PCC_"+name_d[i]+"\tdeltaPCC_"+name_d[i]+"\tp-value_"+name_d[i]+"\n")           
        f=open(background)
        # gene pair!!!
        index=i
        for p in f:
            t=p.split()
            r1=stat.pearsonr(normal[t[0]],normal[t[1]])[0]
            # pearson relations of one gene pair in n samples
            r2=stat.pearsonr(normal[t[0]]+[disease[t[0]][i]],normal[t[1]]+[disease[t[1]][i]])[0]
            # pearson relations of one gene pair in n+1 samples(added one sample is changing)
            r=r2-r1
            #delta 
            if r1==1:
                r1=0.99999999
            if r1==-1:
                r1=-0.99999999
            z=r/((1-r1*r1)/(len(normal[t[0]])-1))
            # in case of infinate number
            pvalue=1-stat.norm.cdf(abs(z))
            #calculate pvalue based on z score
            if pvalue < p_value:    
                fw.write(t[0]+"\t"+t[1]+"\t"+str(r1)+"\t"+str(r2)+"\t"+str(r)+"\t"+str(pvalue)+"\n")

        f.close()
        fw.close()



print("Begin time: "+begin)
print("End time: "+time.asctime())





