from config import *
import pandas as pd
import numpy as np
import glob

def mvfunc(measures, sub_list, strlist):
    """
    prints out the mean and variance of local measures among nodes for each subject.
    Inputs:
    measures = a variable including network measures for all subjects defined in notebook 08.
    sub-list = subjects list
    strlist = denoising_strategies list

    outputs:
    a dictionary including 2 dataframes of shape (n,m) corresponded to mean and variance.
    where n= number of local network measures
          m= number of subjects
    """
    mean_loc = {}
    var_loc  = {}
    for i in range(int(len(strlist)/2)):
        for sub in sub_list:
            m = measures[strlist[i]][sub].mean(axis = 0).to_frame()
            v = measures[strlist[i]][sub].var( axis = 0).to_frame()

            try:
                mean_loc[strlist[i]] = pd.concat([mean_loc[strlist[i]],m],axis=1)
                var_loc[strlist[i]]  = pd.concat([var_loc[ strlist[i]],v],axis=1)
            except KeyError:
                mean_loc[strlist[i]] = m
                var_loc[strlist[i]]  = v
        mean_loc[strlist[i]]  = mean_loc[strlist[i]]
        var_loc[strlist[i]]   = var_loc [strlist[i]]
        mean_loc[strlist[i]].columns = sub_list
        var_loc[strlist[i]].columns  = sub_list
    return {'mean': mean_loc, 'var': var_loc}

#########################################################################################

def mvfunc_total(measures,groups,strlist,ds):
    """prints out the mean and the variance of local measures
    Inputs:
    measures = includes network measures for all subjects, defined in notebook 08.
    groups = list of subjects groups
    strlist = denoising_strategies list

    outputs:
    n*m dataframe
        where:
        n = number of local network measures
        m = 2 * number of groups
    """
    global mvfunc
    totalm = {}
    nstr = int(len(strlist)/2)
    local_strlist = strlist[:nstr]
    for g in groups:
        totalm[g]={}
    for strategy in local_strlist:
        col = []
        for g in groups:
            totalm[g][strategy] = {}
            col.extend([g+'_mean',g+'_var'])

            d = measures[g][strategy]
            totalm[g][strategy]['nan'] = pd.DataFrame({})
            for sub in d:
                totalm[g][strategy]['nan'] = pd.concat([totalm[g][strategy]['nan'],d[sub]])

    mv_total= pd.DataFrame({})
    for g in groups:
        mv_total=pd.concat([mv_total,mvfunc(totalm[g],['nan'],strlist)['mean']['loc_'+ds]],axis=1)
        mv_total=pd.concat([mv_total,mvfunc(totalm[g],['nan'],strlist)['var']['loc_'+ds]],axis=1)

    mv_total.columns  = col
    return mv_total
    #return totalm

#########################################################################################

def mes_rel(measures,groups,strlist,ds):
    """
    computes covariance and correlation between local network measures

    inputs:
    measures = includes network measures for all subjects, defined in notebook 08.
    groups = list of subjects groups
    strlist = denoising_strategies list

    outputs
    an n*n relationship matrix where n is the number of local measures
    """
    subjects={}
    for g in groups:
        subjects[g] = list(measures[g]['loc_'+ds].keys())

    mes_cov  = {}
    mes_corr = {}
    for i in range(int(len(strlist)/2)):
        mes_cov[strlist[i]] = {}
        mes_corr[strlist[i]] = {}
        for g in groups:

            mes_cov [strlist[i]][g] = mvfunc(measures[g], subjects[g], strlist)['mean'][strlist[i]].transpose().cov()
            mes_corr[strlist[i]][g] = mvfunc(measures[g], subjects[g], strlist)['mean'][strlist[i]].transpose().corr()
    return {'cov': mes_cov, 'corr': mes_corr}

#########################################################################################

def sub_rel(measures,groups,strlist,ds):
    """
    computes covariance and correlation between subjects.

    inputs:
    measures = a variable including network measures for all subjects, defined in notebook 08.
    groups = groups of patients
    strlist = denoising_strategies list

    outputs
    an n*n relationship matrix where n is the number of subject
    """

    subjects={}
    for g in groups:
        subjects[g] = list(measures[g]['loc_'+ds].keys())
    sub_cov  = {}
    sub_corr = {}

    for i in range(int(len(strlist)/2)):
        sub_cov[strlist[i]] = {}
        sub_corr[strlist[i]] = {}
        for g in groups:
            sub_cov [strlist[i]][g] = mvfunc(measures[g], subjects[g], strlist)['mean'][strlist[i]].cov()
            sub_corr[strlist[i]][g] = mvfunc(measures[g], subjects[g], strlist)['mean'][strlist[i]].corr()

    return {'cov': sub_cov, 'corr': sub_corr}

#########################################################################################

def read_M(sub,ds,corr):
    """ read the correlation matrix for a subject with a given denoising strategy and correlation type
    """

    filesnp = glob.glob(rootdir + "/data/04_correlations/corr-%s/ds-%s/*%s*.npy"
                      %(corr,ds,sub))
    return np.load(filesnp[0])

#########################################################################################

def read_adj_M(sub,ds,corr,tm,negative):
    """ read the adjacency matrix for a subject with a given denoising strategy, correlation type, and thresholding method
    """
    if negative:
    	filesnp = glob.glob(rootdir + "/data/05_adjacency_matrices/negative/*%s*/corr-%s/ds-%s/*%s*.npy"
    	                  %(tm,corr,ds,sub))
    else:
        filesnp = glob.glob(rootdir + "/data/05_adjacency_matrices/positive/*%s*/corr-%s/ds-%s/*%s*.npy"
                          %(tm,corr,ds,sub))
    return np.load(filesnp[0])
