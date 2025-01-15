#!/usr/bin/env python3
import sys, os
import numpy as np
import pandas as pd
from pandas.io.parsers import read_csv
import itertools as it
from collections import Counter
from gtfparse import read_gtf
import networkx as nx
import pickle

try:
    sample = sys.argv[1]
    dirsubgtf = sys.argv[2]
except:
    sys.exit("Please, (1) import sample name (same before .fc_small.tsv and _small.csv) & (2) output root name (to import _small_annot.gtf)")

# read gtf
gtf = read_gtf(dirsubgtf+'_small_annot.gtf')
gtf = pd.DataFrame(gtf, columns = gtf.columns)
gtf = gtf[gtf['feature']=='transcript']

##########################
# Multi-loci assignation # ONLY FOR SMALL BIOTYPES
##########################

# counts
#cdf = read_csv(sample + '_small.csv', sep = '\t', comment = '#')
#fcdf = cdf[cdf[cdf.columns[-1]]>0]

# edges
fdf = read_csv(sample + '.fc_small.tsv', sep = '\t', header = None, names = ['read','feature'])
if len(fdf) == 0:
    G = nx.Graph()
    pickle.dump(G, open(sample+'_small_graph.pickle', 'wb'))
    sys.exit()

fdf['feature'] = fdf['feature'].apply(lambda x: x.rsplit(','))
aggfdf = fdf.groupby('read').agg('sum')
cnt_feat = Counter(aggfdf['feature'].sum())
faggfdf = aggfdf[aggfdf['feature'].apply(lambda x: len(x))>1]
if len(faggfdf) == 0:
    G = nx.Graph();
    for idx in cnt_feat:
        G.add_node(idx, size = cnt_feat[idx], biotype = gtf[gtf['transcript_name']==idx]['transcript_biotype'].iloc[0])
    pickle.dump(G, open(sample+'_small_graph.pickle', 'wb'))
    sys.exit()
edges = faggfdf['feature'].apply(lambda x: [y for y in it.combinations(x,2)]).sum()
cnt_edges = Counter(edges)

# create network
G = nx.Graph()
#for idx in fcdf.index:
#    G.add_node(fcdf.loc[idx,'Geneid'], size = fcdf.loc[idx,fcdf.columns[-1]], biotype = gtf[gtf['transcript_name']==fcdf.loc[idx,'Geneid']]['transcript_biotype'].iloc[0])
for idx in cnt_feat:
    G.add_node(idx, size = cnt_feat[idx], biotype = gtf[gtf['transcript_name']==idx]['transcript_biotype'].iloc[0])

for edge in cnt_edges:
    G.add_edge(edge[0], edge[1], weight = cnt_edges[edge])

pickle.dump(G, open(sample+'_small_graph.pickle', 'wb'))
