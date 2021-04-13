import pandas as pd
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-i','--input', type=str,help="The event file")
parser.add_argument('-o','--output', type=str,help="The directory of output",default='./')
parser.add_argument('-n','--number_of_top_event', type=int,help="The number of top(sorted by pvalue) events to plot", default=10)

args = parser.parse_args()

all=pd.read_table(args.input)

all.sort_values(['PValue'], ascending = True).iloc[0:args.number_of_top_event,:].to_csv(args.output+'/top_events.txt', sep='\t', index=False)
