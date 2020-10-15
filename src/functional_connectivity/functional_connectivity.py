import numpy as np
import matplotlib.pyplot as plt
from config import *
import os,sys
import pandas as pd
from pandas import DataFrame as df
import networkx as nx
import nilearn
from nilearn.connectome import ConnectivityMeasure
from sklearn import covariance, preprocessing

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
	"""Compute Pearson correlation
	sub : subjects ID
	"""
	from src.data import subjects
	print("processing sub-%s"%sub)
	for ds in denoising_strategies:
		ts = np.load("%s/ds-%s/sub-%s_ds-%s.npy"%(subjects.time_seriesdir, ds, sub, ds))
		correlation = np.corrcoef(ts[:,1:].T)
		# we end up with a few NAN values because of an empty ROI, for now just zero them out
		correlation[np.isnan(correlation)] = 0

		np.save("%s/data/04_correlations/corr-%s/ds-%s/sub-%s_ds-%s_corr-%s"
			%(rootdir,'pearson',ds,sub,ds,'pearson'),correlation)
		Gr=nx.from_numpy_array(correlation)
		nx.write_gexf(Gr,"%s/data/04_correlations/corr-%s/ds-%s/sub-%s_ds-%s_corr-%s.gexf"
			%(rootdir,'pearson',ds,sub,ds,'pearson'))


def all_pearsons(subs):
	"""Compute Pearson correlation for all subjects in list subs
	subs : list of subjects IDs
	"""
	print('Computing Pearson correlation')
	for sub in subs:
		get_pearson(sub)


def get_partial(sub):
	"""Compute partial correlation
	sub : subjects ID
	"""
	from src.data import subjects
	print("processing sub-%s"%sub)
	for ds in denoising_strategies:
		ts = np.load("%s/ds-%s/sub-%s_ds-%s.npy"%(subjects.time_seriesdir, ds, sub, ds))
		estimator=nilearn.connectome.ConnectivityMeasure(kind='partial correlation')
		ts = [ts[:,1:]]
		correlation = estimator.fit_transform(ts)  # drop the zero roi data
		correlation = correlation[0]
		correlation[np.diag_indices_from(correlation)]=0
		# we end up with a few NAN values because of an empty ROI, for now just zero them out
		correlation[np.isnan(correlation)]=0
		np.save("%s/data/04_correlations/corr-%s/ds-%s/sub-%s_ds-%s_corr-%s"
			%(rootdir,'partial',ds,sub,ds,'partial'),correlation)
		Gr=nx.from_numpy_array(correlation)
		nx.write_gexf(Gr,"%s/data/04_correlations/corr-%s/ds-%s/sub-%s_ds-%s_corr-%s.gexf"
			%(rootdir,'partial',ds,sub,ds,'partial'))


def all_partial(subs):
	"""Compute partial correlation for all subjects in list subs
	subs : list of subjects IDs
	"""
	print('Computing partial correlation')
	for sub in subs:
		get_partial(sub)


def fALPHA(ts, alphaRange, tol, sub, corr):
    """
    """
    import networkx as nx
    global argmax,ar;
    GCE = []
    alphaRange1 = alphaRange
    correlations = {}
    Scaler = preprocessing.StandardScaler()
    X = Scaler.fit_transform(ts[:,1:])
    emp_cov = covariance.empirical_covariance(X)
    shrunk_cov = covariance.shrunk_covariance(emp_cov, shrinkage=0.8)
    correlations = {}
    selected_alpha = alphaRange[0]
    for alpha in alphaRange:
        if corr == 'all':
            try:
                corrM= covariance.graphical_lasso(shrunk_cov, alpha,max_iter=700)[1]
                corrM[np.diag_indices_from(corrM)]=0
                corrG =nx.from_numpy_matrix(corrM)
                Gcc = sorted(nx.connected_components(corrG),key = len)
                G0 = corrG.subgraph(Gcc[-1])
                if len(G0.nodes)>=359:
                    selected_alpha = alpha
                else:
                    break
            except FloatingPointError:
                print("Failed at subject %s with alpha=%s"%(sub,alpha))

        elif corr == 'separated':
            try:
                corrM= covariance.graphical_lasso(shrunk_cov, alpha,max_iter=700)[1]
                corrM[np.diag_indices_from(corrM)]=0
                PcorrM = corrM.clip(min=0)
                NcorrM = corrM.clip(max=0)

                PcorrG =nx.from_numpy_matrix(PcorrM)
                NcorrG =nx.from_numpy_matrix(NcorrM)

                Gcc = sorted(nx.connected_components(PcorrG),key = len)
                PG0 = PcorrG.subgraph(Gcc[-1])
                Gcc = sorted(nx.connected_components(NcorrG),key = len)
                NG0 = NcorrG.subgraph(Gcc[-1])
                if len(PG0.nodes)>=359 and len(NG0.nodes)>=359:
                    selected_alpha = alpha
                else:
                    break
            except FloatingPointError:
                print("Failed at subject %s with alpha=%s"%(sub,alpha))


    print('Alpha:         ',alphaRange)
    print('Selected Alpha:',selected_alpha)
    l = (alphaRange[1] - alphaRange[0])
    NEWalphaRange = np.arange(selected_alpha,
                            selected_alpha+l,
                            l/10)

    if l>tol+tol/10:
        fALPHA(ts,NEWalphaRange,tol,sub,corr)
    print('==============================================')
    return selected_alpha

def get_glasso(ts, alphaRange, tol,sub,corr,ds):
    """
    """
    print('sub: %s, with %s denoising strategies'%(sub,ds))
    alpha = fALPHA(ts=ts, alphaRange=alphaRange,
                   tol=tol, sub=sub,corr=corr)
    if corr == 'all':
        Scaler = preprocessing.StandardScaler()
        X = Scaler.fit_transform(ts[:,1:])
        emp_cov = covariance.empirical_covariance(X)
        shrunk_cov = covariance.shrunk_covariance(emp_cov, shrinkage=0.8)
        correlations = covariance.graphical_lasso(shrunk_cov, alpha, max_iter=500)[1]
        correlations[np.diag_indices_from(correlations)]=0

    elif corr == 'separated':
        Scaler = preprocessing.StandardScaler()
        X = Scaler.fit_transform(ts[:,1:])
        emp_cov = covariance.empirical_covariance(X)
        shrunk_cov = covariance.shrunk_covariance(emp_cov, shrinkage=0.8)
        correlations = {'positive':{},'negative':{}}
        corrM = covariance.graphical_lasso(shrunk_cov, alpha, max_iter=500)[1]
        corrM[np.diag_indices_from(corrM)]=0
        correlations['positive'] = corrM.clip(min=0)
        correlations['negative'] = corrM.clip(max=0)

    return(correlations)

def all_glasso(timeseries, alphaRange, tol, sub_list, corr):
    for sub in sub_list:
        for ds in denoising_strategies:
            corrM = get_glasso(ts=timeseries[sub][ds],alphaRange=alphaRange,
                                      tol=tol,sub=sub,corr=corr,ds=ds)
        if corr=='all':
            os.system('mkdir -p %s/data/04_correlations/corr-glasso/all/ds-%s/'
                      %(rootdir,ds))
            nx.write_gexf(nx.from_numpy_array(corrM),"%s/data/\
04_correlations/corr-glasso/all/ds-%s/sub-%s_ds-%s_corr-glasso-all.gexf"
                          %(rootdir,ds,sub,ds))
            np.save("%s/data/04_correlations/corr-glasso/all/ds-%s/\
sub-%s_ds-%s_corr-glasso-all"
                    %(rootdir,ds,sub,ds), corrM)

        if corr=='separated':
            os.system('mkdir -p %s/data/04_correlations/corr-glasso/separated/ds-%s/'
                      %(rootdir,ds))
            nx.write_gexf(nx.from_numpy_array(corrM['positive']),"%s/data/\
04_correlations/corr-glasso/separated/ds-%s/sub-%s_ds-%s_corr-glasso-positive.gexf"
                          %(rootdir,ds,sub,ds))
            np.save("%s/data/04_correlations/corr-glasso/separated/ds-%s/\
sub-%s_ds-%s_corr-glasso-positive"
                    %(rootdir,ds,sub,ds), corrM['positive'])
            nx.write_gexf(nx.from_numpy_array(corrM['negative']),"%s/data/\
04_correlations/corr-glasso/separated/ds-%s/sub-%s_ds-%s_corr-glasso-negative.gexf"
                          %(rootdir,ds,sub,ds))
            np.save("%s/data/04_correlations/corr-glasso/separated/ds-%s/\
sub-%s_ds-%s_corr-glasso-negative"
                    %(rootdir,ds,sub,ds), corrM['negative'])
    return 0
