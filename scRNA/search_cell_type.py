import argparse
import numpy as np
import pandas as pd
import re

parser = argparse.ArgumentParser()
parser.add_argument('-m','--database', help="Database used for annotating cell types, CellMarker and PanglaoDB are supported, CellMarker support human genes and mouse genes",default='CellMarker')
parser.add_argument('-g','--genelist', help="Significant genes for specific cell clusters",default='')
parser.add_argument('-r','--organism', help="Human and Mouse are supported",default='Human')
args = parser.parse_args()
if (not len(args.database)) or (not len(args.genelist)):
    print('please input the database and genelist')
    sys.exit(1)

if args.database=='CellMarker':
    genelist=pd.read_csv(args.genelist,header=None)
    test=pd.read_csv('markers_database/CellMarker_filtered.csv')
    test1=test.loc[(test['speciesType']==args.organism),:]
    celltype=[]
    marker=[]
    print('Gene list you provided')
    print(genelist)
    for i in range(0,len(genelist)):
        for j in range(0,len(test1)):
            if ((genelist.iloc[i,0]) in (test1.iloc[j,9])):
                celltype.append(test1.iloc[j,6])
                marker.append(genelist.iloc[i,0])
    statistic={}
    for i in celltype:
        if i not in statistic.keys():
            statistic[i]=1
        else:
            statistic[i]=statistic[i]+1
    statistic=sorted(statistic.items(), key=lambda e:e[1], reverse=True)
    print('Below shows the enriched cell type and confidance(How many records support this)')
    print(pd.DataFrame(statistic,columns=['CellType', 'Number_of_records_supported']))

if args.database=='PanlaoDB':
    genelist0=pd.read_csv(args.genelist,header=None)
    test=pd.read_csv('markers_database/PanglaoDB_markers_27_Mar_2020.csv')
    if args.organism=='Human':
        test1=test.loc[(test['species']=='Mm Hs') | (test['species']=='Hs'),:]
        genelist=genelist0
    if args.organism=='Mouse':
        test1=test.loc[(test['species']=='Mm Hs') | (test['species']=='Mm'),:]
        genelist=[]
        for i in genelist0.iloc[:,0]:
            genelist.append(i.upper())
        genelist=pd.DataFrame(genelist)
    celltype=[]
    marker=[]
    print('Gene list you provided')
    print(genelist)
    for i in range(0,len(genelist)):
        for j in range(0,len(test1)):
            if ((genelist.iloc[i,0]) in (test1.iloc[j,1])):
                celltype.append(test1.iloc[j,2])
                marker.append(genelist.iloc[i,0])
    statistic={}
    for i in celltype:
        if i not in statistic.keys():
            statistic[i]=1
        else:
            statistic[i]=statistic[i]+1
    statistic=sorted(statistic.items(), key=lambda e:e[1], reverse=True)
    print('Below shows the enriched cell type and confidance(How many records support this)')
    print(pd.DataFrame(statistic,columns=['CellType', 'Number_of_records_supported']))



