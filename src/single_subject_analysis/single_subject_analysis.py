from config import *
import pandas as pd
import numpy as np
import networkx as nx
import scipy.stats
from sklearn import metrics
import bct
import matplotlib.pyplot as plt

def get_adjmtx(corrmtx,density,verbose=False):
    assert density<=1
    cutoff=scipy.stats.scoreatpercentile(corrmtx[np.triu_indices_from(corrmtx,1)],
                                         100-(100*density))

    # scipy.stats.scoreatpercentile:
    # Calculate the score at a given percentile of the input sequence.

    # np.triu_indices_from
    # Return the indices for the upper-triangle of arr.

    if verbose:
        print('cutoff:%0.3f'%cutoff)
    adjmtx=(corrmtx>cutoff).astype('int')
    adjmtx[np.diag_indices_from(adjmtx)]=0
    return(adjmtx)

# now generate a graph using NetworkX
# created previously using get_yeo_assignments.py
labeldir   = '%s/references/HCP-MMP1/MMP_yeo2011_networks.csv'%(rootdir)
labeldata  = pd.read_csv(labeldir)
def gengraph(adjmtx):
    G=nx.from_numpy_array(adjmtx)

    # get giant component
    Gc = max(nx.connected_components(G), key=len)
    Gc=G.subgraph(Gc)
    print('Giant component includes %d out of %d total nodes'%(len(Gc.nodes),len(G.nodes)))
    labeldata_Gc=labeldata.loc[list(Gc.nodes)]

    cl={0:'black',1:'red',2:'yellow',3:'green',4:'blue',5:'orange',6:'gray',7:'magenta'}
    colors=[cl[labeldata['Yeo7'].iloc[i]] for i in Gc.nodes]
    degrees=np.array([Gc.degree(i) for i in Gc.nodes])
    layout=nx.spring_layout(Gc)
    nx.draw_networkx(Gc,pos=layout,with_labels=False,node_color=colors,
                  node_size=degrees)
    _=plt.axis('off')
    yeodict={0:'Undefined',1:'Visual',2:'Somatomotor',3:'DorsalAttention',
             4:'VentralAttention',5:'Limbic',
             6:'Frontoparietal',7:'Default'}

    for i in yeodict:
        print(cl[i],':',yeodict[i])

        
def comdetc(corrmtx, adjmtx, density):
    # get adj matrix for giant component
    G=nx.from_numpy_array(adjmtx)

    # get giant component
    Gc = max(nx.connected_components(G), key=len)
    Gc = G.subgraph(Gc)
    print('Giant component includes %d out of %d total nodes'%(len(Gc.nodes),len(G.nodes)))
    labeldata_Gc=labeldata.loc[list(Gc.nodes)]

    Gc_nodelist=list(Gc.nodes)
    tmp=corrmtx[Gc_nodelist,:]
    corrmtx_Gc=tmp[:,Gc_nodelist]
    adjmtx=get_adjmtx(corrmtx_Gc,density)

    mod_binary=bct.modularity_louvain_und(adjmtx)
    print('modularity:',mod_binary[1])

    print('Multilevel modularity optimization identifed %d communities'%len(np.unique(mod_binary[0])))
    ari=metrics.adjusted_rand_score(mod_binary[0],
                                        labeldata_Gc['Yeo7'])
    print('Adjusted Rand index compared to Yeo 7 networks: %0.3f'%ari)

    degrees=np.array([Gc.degree(i) for i in Gc.nodes])
    layout=nx.spring_layout(Gc)
    nx.draw_networkx(Gc,pos=layout,with_labels=False,
            node_color=[mod_binary[0][i] for i in range(len(Gc.nodes))],
            node_size=degrees,cmap='viridis')

    _=plt.axis('off')

    
def clk(G,k):
    """computes average clustering coefficient for nodes with degree k"""
    ls= list(G.degree(nx.nodes(G)))
    s=0
    c=0
    for i in ls:
        if i[1]==k:
            s=s+ nx.clustering(G, i[0])
            c=c+1
    return s/c


#small world
def ml(G,l):
    """
    it computes the average number of nodes within a distance less than or equal l to
    from any given vertex.
    """
    s=0
    for j in G.nodes:
        s=s+len(nx.single_source_shortest_path_length(G, j, cutoff =l))-1 #-1 becouse it counts distance(i,i)=0<cutoff
    return s/nx.number_of_nodes(G)
