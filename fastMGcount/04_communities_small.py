#!/usr/bin/env python3
import sys, os
import numpy as np
import pandas as pd
from pandas.io.parsers import read_csv
import pickle
from infomap import Infomap # https://mapequation.github.io/infomap/python/
import networkx as nx
import igraph as ig
from collections import Counter

try:
    outdir = sys.argv[1]
    outname = sys.argv[2]
    inbamfile = sys.argv[3]
except:
    sys.exit("Please, provide (1) output directory; (2) list of bamfiles (tsv file, columns 1 are ref bamfiles, column 2 are sample name which _small_graph.pickle is based on")

samples = read_csv(inbamfile, sep = '\t', header = None, names = ['bam','sample'])
Gs = {s: pickle.load(open(outdir + '/' + s+'_small_graph.pickle','rb')) for s in samples['sample']}

# add all networks
G = nx.Graph()
for s in Gs:
    for node in Gs[s].nodes:
        try:
            G.nodes[node]['size'] += Gs[s].nodes[node]['size']
        except:
            G.add_node(node, size = Gs[s].nodes[node]['size'], biotype = Gs[s].nodes[node]['biotype'])
for s in Gs:
    for edge in Gs[s].edges:
        try:
            G.edges[edge]['weight'] += Gs[s].edges[edge]['weight']
        except:
            G.add_edge(edge[0], edge[1], weight = Gs[s].edges[edge]['weight'])

# Rosvall equation
btypeloop = set([G.nodes[n]['biotype'] for n in G.nodes if len(G.nodes[n])>0])
MGcommunities = pd.Series()
for btype in btypeloop:
    nodes = [n for n in G.nodes if len(G.nodes[n])>0 and G.nodes[n]['biotype']==btype]
    feat2mg = pd.Series()
    H = G.subgraph(nodes)
    feat2int = {n: i for i, n in enumerate(H.nodes)}
    int2feat = {feat2int[n]: n for n in feat2int}
    im = Infomap("--two-level --directed")
    i = 0
    for n in H.nodes:
        edges = [e for e in H.edges if e[0]==n]
        if len(edges) > 0:
            i += 1
            for e in edges:
                im.add_link(feat2int[e[0]], feat2int[e[1]], H.edges[e]['weight']/H.nodes[e[0]]['size'])
    if i > 0:
        im.run()
        for node in im.tree:
            if node.is_leaf:
                feat2mg.loc[int2feat[node.node_id]] = '-'.join([btype,str(node.module_id)])
    MGcommunities = pd.concat([MGcommunities,feat2mg])

    h = ig.Graph.from_networkx(H)
    layout = h.layout(layout='auto')
    pos = {n: np.array(layout.coords[i]) for i, n in enumerate(H.nodes)}
    i = 0
    for node in pos:
        H.nodes[node]['position'] = pos[node]
        if node in list(feat2mg.index):
            H.nodes[node]['id'] = feat2mg.loc[node]
        else:
            i += 1
            H.nodes[node]['id'] = '-'.join([btype,str(im.num_top_modules+i+1)])

    pickle.dump(H, open(outdir + '/' + btype + '_MGgraph.pickle', 'wb'))

    # save H :)
    '''
    fig, ax = plt.subplots()
    for node in H.nodes:
        ax.scatter(H.nodes[node]['position'][0], H.nodes[node]['position'][1], c = Colors.colors[int(H.nodes[node]['id'].rsplit('-')[-1])], s = H.nodes[node]['size'], zorder=2)
    for edge in H.edges:
        ax.plot([H.nodes[edge[0]]['position'][0],H.nodes[edge[1]]['position'][0]], [H.nodes[edge[0]]['position'][1],H.nodes[edge[1]]['position'][1]], c = 'silver', lw = 0.5*H.edges[edge]['weight'], zorder=1)
    ax.set_title(btype); ax.set_xticks([]); ax.set_yticks([])
    fig.savefig(outdir + '/' + btype + '_MGgraph.pdf', bbox_inches = 'tight')
    '''

for i,s in enumerate(samples['sample']):
    sdf = read_csv(outdir + '/' + s + '_small.csv', sep = '\t', comment='#', index_col = [0,1,2,3,4,5]); sdf.columns = [s]
    edf = read_csv(outdir + '/' + s + '_long_exons.csv', sep = '\t', comment='#', index_col = [0,1,2,3,4,5]); edf.columns = [s]
    idf = read_csv(outdir + '/' + s + '_long_introns.csv', sep = '\t', comment='#', index_col = [0,1,2,3,4,5]); idf.columns = [s]
    if i == 0:
        small_df = sdf
        exon_df = edf
        intron_df = idf
    else:
        small_df = small_df.merge(sdf, how = 'outer', left_index = True, right_index = True)
        exon_df = exon_df.merge(edf, how = 'outer', left_index = True, right_index = True)
        intron_df = intron_df.merge(idf, how = 'outer', left_index = True, right_index = True)

small_df['MGcomm'] = [MGcommunities.loc[idx[0]] if idx[0] in MGcommunities.index else idx[0] for idx in small_df.index]
agg_small_df = small_df.groupby('MGcomm').sum()

small_df.to_csv(outdir + '/' + outname + '_small_counts.tsv', sep = '\t')
exon_df.to_csv(outdir + '/' + outname + '_exonic_counts.tsv', sep = '\t')
intron_df.to_csv(outdir + '/' + outname + '_intronic_counts.tsv', sep = '\t')

exon_df.index = [idx[0] + '-exon' for idx in exon_df.index]
intron_df.index = [idx[0] + '-intron' for idx in intron_df.index]
agg_df = pd.concat([exon_df,intron_df,agg_small_df])
agg_df.to_csv(outdir + '/' + outname + '_aggregated_counts.tsv', sep = '\t')

