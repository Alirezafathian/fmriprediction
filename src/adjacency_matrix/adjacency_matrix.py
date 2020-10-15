from config import *
import numpy as np
import networkx as nx
import scipy.stats
import glob, os


def read_G(sub,ds,corr):
    """ read the correlation graph for a subject with a given denoising strategy
    """
    
    files = glob.glob(rootdir + "/data/04_correlations/corr-%s/ds-%s/*%s*.gexf"
                      %(corr,ds,sub))
    return nx.read_gexf(files[0])


def read_M(sub,ds,corr):
    """ read the correlation matrix for a subject with a given denoising strategy
    """
    
    filesnp = glob.glob(rootdir + "/data/04_correlations/corr-%s/ds-%s/*%s*.npy"
                      %(corr,ds,sub))
    return np.load(filesnp[0])


def get_adjmtx(corrmtx,density,verbose=False):
    """ generating adjacency matrix by thresholding the correlation matrix
    """
    
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


def fPSW(G, t):
    """
    """
    
    G1     = G.copy()
    e1     = G1.number_of_edges()
    remove = [(n1,n2) for (n1,n2,w) in list(G1.edges.data('weight')) if w<t]
    G1.remove_edges_from(remove)
    e2     = G1.number_of_edges()
    psw    = e2/e1
    return {'G': G1, 'psw': psw}


def fGCE(G, thresholds, tol):
    """
    """
    
    global argmax,tsh
    tsh  = thresholds
    G1   = G.copy()
    diag = [(i,i) for i in list(G1.nodes)]
    G1.remove_edges_from(diag)

    GCE = []
    for t in thresholds:
        psw = fPSW(G1,t)
        E   = nx.global_efficiency(psw['G'])
        gce = E-psw['psw']
        GCE.append(gce)

    GCE = np.asarray(GCE)
    print('Threshold: ',thresholds)
    print('GCE:       ',GCE)
    argmax = np.argmax(GCE)
    print('Maximum GCE:\nThreshold: %s, Global Efficiency: %s\n==============================================\n'%(thresholds[argmax],GCE[argmax]))

    GCE2 = GCE.copy()
    GCE2[np.argmax(GCE2)] = -100
    r = [thresholds[np.argmax(GCE)], thresholds[np.argmax(GCE2)]]
    thresholds1 = np.arange(min(r),max(r),(max(r)-min(r))/10)

    if thresholds[2]-thresholds[1]>tol+tol/10:
        fGCE(G,thresholds1,tol)
    print('==============================================')
    return fPSW(G1,tsh[argmax])['G']


def userdefined_adj(sub,correlation):
    """
    """

    for ds in denoising_strategies:
        for ct in correlation_types:
            for td in thresholding_values:
                                      
                print('sub-%s, ds-%s, corr-%s, td-%.3f: '%(sub,ds,ct,td))
                adjM = get_adjmtx(correlation[ds][ct]['M'],td,verbose=True)
                adjG = nx.from_numpy_matrix(adjM)
                    
                os.system('mkdir -p %s/data/05_adjacency_matrices/positive/tm-userdefined-%.3f/\
corr-%s/ds-%s/'%(rootdir,td,ct,ds))

                nx.write_gexf(adjG,"%s/data/05_adjacency_matrices/positive/tm-userdefined-%.3f/\
corr-%s/ds-%s/sub-%s_ds-%s_corr-%s_tm-userdefined_density-%.3f.gexf"
                              %(rootdir,td,ct,ds,sub,ds,ct,td))
                    
                np.save("%s/data/05_adjacency_matrices/positive/tm-userdefined-%.3f/corr-%s/\
ds-%s/sub-%s_ds-%s_corr-%s_tm-userdefined_density-%.3f"
                         %(rootdir,td,ct,ds,sub,ds,ct,td), adjM)
                
def all_userdefined(subs,correlations):
    for sub in subs:
        userdefined_adj(sub,correlations[sub])
                    
                    
def gce_adj(sub,correlation, tol):
    """
    """
    
    thresholds = np.arange(0,1,.1)
    for ds in denoising_strategies:
        for ct in correlation_types:
            print('sub-%s, ds-%s, corr-%s: '%(sub,ds,ct))
            adjG = fGCE(correlation[ds][ct]['G'], thresholds,tol)
            for n1, n2, d in adjG.edges(data=True):
                d.pop('weight', None)
                adjM = nx.to_numpy_array(adjG)
    
            os.system('mkdir -p %s/data/05_adjacency_matrices/positive/tm-gce/\
corr-%s/ds-%s/'%(rootdir,ct,ds))
                    
            nx.write_gexf(adjG,"%s/data/05_adjacency_matrices/positive/tm-gce/\
corr-%s/ds-%s/sub-%s_ds-%s_corr-%s_tm-gce.gexf"
                          %(rootdir,ct,ds,sub,ds,ct))
                    
            np.save("%s/data/05_adjacency_matrices/positive/tm-gce/corr-%s/\
ds-%s/sub-%s_ds-%s_corr-%s_tm-gce"
                    %(rootdir,ct,ds,sub,ds,ct), adjM)

            
def all_gce(subs,correlations,tol):
    for sub in subs:
        gce_adj(sub,correlations[sub],tol)



def Ngce_adj(sub,correlation, tol):
    """
    """
    
    thresholds = np.arange(0,1,.1)
    for ds in denoising_strategies:
        for ct in correlation_types:
            print('sub-%s, ds-%s, corr-%s: '%(sub,ds,ct))
            M = correlation[ds][ct]['M'].copy()
            M = -M
            Gr = nx.from_numpy_array(M)
            adjG = fGCE(Gr, thresholds,tol)
            for n1, n2, d in adjG.edges(data=True):
                d.pop('weight', None)
                adjM = nx.to_numpy_array(adjG)
    
            os.system('mkdir -p %s/data/05_adjacency_matrices/negative/tm-gce/\
corr-%s/ds-%s/'%(rootdir,ct,ds))
                    
            nx.write_gexf(adjG,"%s/data/05_adjacency_matrices/negative/tm-gce/\
corr-%s/ds-%s/sub-%s_ds-%s_corr-%s_tm-gce.gexf"
                          %(rootdir,ct,ds,sub,ds,ct))
                    
            np.save("%s/data/05_adjacency_matrices/negative/tm-gce/corr-%s/\
ds-%s/sub-%s_ds-%s_corr-%s_tm-gce"
                    %(rootdir,ct,ds,sub,ds,ct), adjM)

            
def all_Ngce(subs,correlations,tol):
    for sub in subs:
        Ngce_adj(sub,correlations[sub],tol)
