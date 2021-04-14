import pandas as pd
from itertools import combinations, permutations

import argparse


parser = argparse.ArgumentParser()

parser.add_argument('-i','--input', type=str,help="The gene expression matrix")
parser.add_argument('-o','--output', type=str,help="The gene pairs")

args = parser.parse_args()

normal=pd.read_csv(args.input)
normal_gene=normal.iloc[:,0]
x1=list(set(normal_gene))
x3=list(combinations(x1, 2))

filename=open(args.output+'/all_gene_pairs.txt','w')
for value in range(len(x3)):
	for i in x3[value]:
		filename.write("{}  ".format(i))
	filename.write('\n')
filename.close()


