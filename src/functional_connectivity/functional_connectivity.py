import numpy as np
import matplotlib.pyplot as plt
from config import *
import os,sys
import numpy as np
import pandas as pd
from pandas import DataFrame as df
import networkx as nx

def reorder_corrs(corrmtx,labeldata,labels):
    """
    reorder correlation matrix according to network labels
    """

    idx=np.lexsort(([i for i in range(labeldata.shape[0])],labeldata[labels]))
    tmp=corrmtx[:,idx]
    return(tmp[idx,:],labeldata.iloc[idx,:])


def plot_reordered_corrs(corrmtx,labeldata,labels,colorbar=True):
    """
    plot correlation matrix after reordering
    """

    corr_reord,labeldata_reord=reorder_corrs(corrmtx,labeldata,labels)
    plt.imshow(corr_reord)
    # find breakpoints and plot lines
    breaks=np.array([int(not i) for i in labeldata_reord[labels].values[:-1]==labeldata_reord[labels].values[1:]])
    breaklocs=np.hstack((np.where(breaks)[0],np.array(corrmtx.shape[0]-1)))
    for b in breaklocs:
        plt.plot([0,corrmtx.shape[0]-1],[b,b],color='w',linewidth=0.5)
        plt.plot([b,b],[0,corrmtx.shape[0]-1],color='w',linewidth=0.5)
    # find label locations
    # add a zero to help find label locations
    breaklocs2=np.hstack(([0],breaklocs))
    label_locs=np.mean(np.vstack((breaklocs,breaklocs2[:-1])),0)
    networks=labeldata_reord[labels].values[breaklocs]
    ax=plt.gca()
    ax.set_yticks(label_locs)
    ax.set_yticklabels(networks)
    if colorbar:
        plt.colorbar()
    
    
def get_pearson(sub):
    """
    """
    from src.data import subjects
    print("processing sub-%s"%sub)
    ts = {}
    for ds in denoising_strategies:
        ts[ds] = np.load("%s/ds-%s/sub-%s_ds-%s.npy"
                          %(subjects.time_seriesdir, ds, sub, ds))
    correlations = {}
    for ct in correlation_types:
        correlations[ct]={}
        for ds in denoising_strategies:
            correlations[ct][ds] = np.corrcoef(ts[ds][:,1:].T)  # drop the zero roi data
            # we end up with a few NAN values because of an empty ROI, for now just zero them out
            correlations[ct][ds][np.isnan(correlations[ct][ds])]=0
    for ds in denoising_strategies:
        for ct in correlations:
            np.save("%s/data/04_correlations/corr-%s/ds-%s/sub-%s_ds-%s_corr-%s"
                    %(rootdir,ct,ds,sub,ds,ct),correlations[ct][ds])
            Gr=nx.from_numpy_array(correlations[ct][ds])
            nx.write_gexf(Gr,"%s/data/04_correlations/corr-%s/ds-%s/sub-%s_ds-%s_corr-%s.gexf"
                          %(rootdir,ct,ds,sub,ds,ct))
     
    
def all_pearsons(subs):
    """
    """
    for sub in subs:
        get_pearson(sub)
        