from config import *
import pandas as pd
import numpy as np
import networkx as nx
import glob, os
import bct
from sklearn import preprocessing


def normalize(df):
    '''
    normalize dataframe columns by mean
    '''
    min_max_scaler      = preprocessing.MinMaxScaler()
    x                   = df.values #returns a numpy array
    x_scaled            = min_max_scaler.fit_transform(x)
    normalized          = pd.DataFrame(x_scaled)
    normalized.columns  = df.columns
    return normalized

def loc_mes_disconnected(G,M):
    """
    computes and saves local measures for a disconnected graph G
    """
    m = {}
    m['ID'] = G.nodes
    m['degree'] = \
    bct.degrees_und(M)
    m['betweenness'] = \
    np.fromiter(nx.betweenness_centrality(G).values(), dtype=float)
    m['betweennessbct'] = \
    bct.betweenness_bin(M)
    m['eigenvector'] = \
    np.fromiter(nx.eigenvector_centrality(G, max_iter=500).values(), dtype=float)
    m['eigenvectorbct'] = \
    bct.eigenvector_centrality_und(M)
    m['katz'] = \
    np.fromiter(nx.katz_centrality_numpy(G).values(), dtype=float)
    m['closeness'] = \
    np.fromiter(nx.closeness_centrality(G).values(), dtype=float)
    m['load'] = \
    np.fromiter(nx.load_centrality(G).values(), dtype=float)
    m['clustering_coef'] = \
    np.fromiter(nx.clustering(G).values(), dtype=float)
    m['clustering_coefbct'] = \
    bct.clustering_coef_bu(M)
    m['pagerank'] = \
    np.fromiter(nx.pagerank(G).values(), dtype=float)
    m['pagerank_d85bct'] = \
    bct.pagerank_centrality(M, d = .85)
    m['subgraph'] = \
    np.fromiter(nx.subgraph_centrality(G).values(), dtype=float)
    m['subgraphbct'] = \
    bct.subgraph_centrality(M)
    m['harmonic'] = \
    np.fromiter(nx.harmonic_centrality(G).values(), dtype=float)

    return pd.DataFrame.from_dict(m)


def glob_mes_disconnected(G):
    """
    computes and saves global measures for a disconnected graph G
    """
    
    m = {}
    # transitivity: the fraction of all possible triangles present in G.
    m['transitivity'] = [nx.transitivity(G)]
    # average_clustering: Compute the average clustering coefficient
    m['average_clustering'] = [nx.average_clustering(G)]
    m['local_efficiency'] = [nx.local_efficiency(G)]
    m['global_efficiency'] = [nx.global_efficiency(G)]
    m['pearson_correlation'] = nx.degree_pearson_correlation_coefficient(G)

    return pd.DataFrame.from_dict(m)


def loc_mes_connected(G):
    """
    computes and saves local measures for a connected graph G
    """
    
    #ordered list of degrees
    lss=[]
    for i in list(G.degree(nx.nodes(G))):
        lss.append(i[1])
    k=list(set(lss))

    m = {}
    m['ID'] = G.nodes
    m['degree'] = lss
    m['eccentricity'] = \
    np.fromiter(nx.eccentricity(G).values(), dtype=float)
    m['betweenness'] = \
    np.fromiter(nx.betweenness_centrality(G).values(), dtype=float)
    m['com_betweenness'] = \
    np.fromiter(nx.communicability_betweenness_centrality(G).values(), dtype=float)
    m['eigenvector'] = \
    np.fromiter(nx.eigenvector_centrality(G, max_iter=500).values(), dtype=float)
    m['katz'] = \
    np.fromiter(nx.katz_centrality_numpy(G).values(), dtype=float)
    m['closeness'] = \
    np.fromiter(nx.closeness_centrality(G).values(), dtype=float)
    m['current_flow_closeness'] = \
    np.fromiter(nx.current_flow_closeness_centrality(G).values(), dtype=float)
    m['load'] = \
    np.fromiter(nx.load_centrality(G).values(), dtype=float)
    m['clustering_coef'] = \
    np.fromiter(nx.clustering(G).values(), dtype=float)
    m['pagerank'] = \
    np.fromiter(nx.pagerank(G).values(), dtype=float)
    m['subgraph'] = \
    np.fromiter(nx.subgraph_centrality(G).values(), dtype=float)
    m['harmonic'] = \
    np.fromiter(nx.harmonic_centrality(G).values(), dtype=float)

    return pd.DataFrame.from_dict(m)


def glob_mes_connected(G):
    """
    computes and saves global measures for a connected graph G
    """
    
    m = {}
    m['density'] = nx.density(G)
    m['average_shortest_path_length'] = nx.average_shortest_path_length(G)
    # transitivity: the fraction of all possible triangles present in G.
    m['transitivity'] = nx.transitivity(G)
    # average_clustering: Compute the average clustering coefficient
    m['average_clustering'] = nx.average_clustering(G)
    m['center'] = [nx.center(G)]
    m['diameter'] = nx.diameter(G)
    m['radius'] = nx.radius(G)
    m['periphery'] = [nx.periphery(G)]
    m['local_efficiency'] = nx.local_efficiency(G)
    m['global_efficiency'] = nx.global_efficiency(G)
    m['pearson_correlation'] = nx.degree_pearson_correlation_coefficient(G)
    # The small-world coefficient defined as: sigma = C/Cr / L/Lr
    #m['sigma'] = nx.sigma(G)
    # The small-world coefficient defined as: omega = Lr/L - C/Cl
    #m['omega'] = nx.omega(G)
    return pd.DataFrame.from_dict(m)


def compute_measures(subjects,
                    denoising_strategies,
                    correlation_types,
                    thresholding_methods,
                    thresholding_values):
    global rootdir
    for sub in subjects:
        for ds in denoising_strategies:
            for ct in correlation_types:
                for tm in thresholding_methods:
                    for tv in thresholding_values:
                        if tm=='userdefined':
                            tm = '%s-%.3f'%(tm,tv)
                        # read data
                        corrGdir = glob.glob(rootdir + "/data/04_correlations/corr-%s/ds-%s/*%s*.gexf"
                                             %(ct,ds,sub))
                        corrMdir = glob.glob(rootdir + "/data/04_correlations/corr-%s/ds-%s/*%s*.npy"
                                             %(ct,ds,sub))
                        adjGdir  = glob.glob(rootdir + "/data/05_adjacency_matrices/tm-*%s*/corr-%s/ds-%s/*%s*.gexf"%(tm,ct,ds,sub))
                        adjMdir  = glob.glob(rootdir + "/data/05_adjacency_matrices/tm-*%s*/corr-%s/ds-%s/*%s*.npy"%(tm,ct,ds,sub))
                        corrG    = nx.read_gexf(corrGdir[0])
                        corrM    = np.load(corrMdir[0])
                        adjG     = nx.read_gexf(adjGdir[0])
                        adjM     = np.load(adjMdir[0])

                        # compute the giant component
                        adjG_gc = [adjG.subgraph(c).copy() for c in nx.connected_components(adjG)]
                        adjG_gc = adjG_gc[0]

                        # compute measures
                        loc_mes_adj     = loc_mes_disconnected(adjG,adjM)
                        glob_mes_adj    = glob_mes_disconnected(adjG)
                        gc_loc_mes_adj  = loc_mes_connected(adjG_gc)
                        gc_glob_mes_adj = glob_mes_connected(adjG_gc)

                        # compute normalized local measures
                        loc_mes_adj_norm       = normalize(loc_mes_adj)
                        loc_mes_adj_norm.ID    = loc_mes_adj.ID
                        gc_loc_mes_adj_norm    = normalize(gc_loc_mes_adj)
                        gc_loc_mes_adj_norm.ID = gc_loc_mes_adj.ID

                        # save measures
                        dirc = rootdir + '/data/06_network_measures/tm-%s/corr-%s/ds-%s/sub-%s'%(tm,ct,ds,sub)
                        os.system('mkdir -p %s'%dirc)
                        loc_mes_adj.to_csv \
                        ('%s/sub-%s_ds-%s_corr-%s_tm-%s_local_measures.csv'
                         %(dirc,sub,ds,ct,tm), sep='\t')

                        glob_mes_adj.to_csv \
                        ('%s/sub-%s_ds-%s_corr-%s_tm-%s_global_measures.csv'
                         %(dirc,sub,ds,ct,tm), sep='\t')

                        gc_loc_mes_adj.to_csv \
                        ('%s/sub-%s_ds-%s_corr-%s_tm-%s_local_measures_giant_component.csv'
                         %(dirc,sub,ds,ct,tm), sep='\t')

                        gc_glob_mes_adj.to_csv \
                        ('%s/sub-%s_ds-%s_corr-%s_tm-%s_global_measures_giant_component.csv'
                         %(dirc,sub,ds,ct,tm), sep='\t')

                        loc_mes_adj_norm.to_csv \
                        ('%s/sub-%s_ds-%s_corr-%s_tm-%s_local_measures_norm.csv'
                         %(dirc,sub,ds,ct,tm),sep='\t')

                        gc_loc_mes_adj_norm.to_csv \
                        ('%s/sub-%s_ds-%s_corr-%s_tm-%s_local_measures_giant_component_norm.csv'
                         %(dirc,sub,ds,ct,tm), sep='\t')

                        print('✓ subject: %s, Denoising Strategy: %s, Correlation Type: %s, \
Thresholding Methods: %s'%(sub,ds,ct,tm))

def add_measure(gl,c,bw,
                mes,mes_name,
                subjects,
                denoising_strategies,
                smethods):
    '''
    adding new measures to the existing measure files for all subjects

    inputs:
    gl = It could be 'global' or 'local'.
    c  = It could be 'connected' or 'disconnected'.
    bw = It could be 'binary' or 'weighted'
    mes = a function for computing the measure values
          for local measures the output would be a scaler
          for global measures the  output is a dictionary
          with node IDs as key and the corresponding meacure as value.
    subjects = list of subject IDs.
    denoising_strategies = list of denoising strategies that you want to include.
    thresholding_methods = list of thresholding methods that you want to include.
    '''

    global rootdir
    for sub in subjects:
        for ds in denoising_strategies:
            for ct in correlation_types:
                for tm in thresholding_methods:
                    for tv in thresholding_values:
                        if tm=='userdefined':
                            tm = '%s-%.3f'%(tm,tv)
                        print(sub,ds,tm,ct)
                        corrGdir = glob.glob(rootdir + "/data/04_correlations/corr-%s/ds-%s/*%s*.gexf"
                                             %(ct,ds,sub))
                        adjGdir  = glob.glob(rootdir + "/data/05_adjacency_matrices/tm-*%s*/corr-%s/ds-%s/*%s*.gexf"%(tm,ct,ds,sub))
                        corrG    = nx.read_gexf(corrGdir[0])
                        adjG     = nx.read_gexf(adjGdir[0])

                        # compute the giant component
                        adjG_gc = [adjG.subgraph(c).copy() for c in nx.connected_components(adjG)]
                        adjG_gc = adjG_gc[0]

                        corrG_gc      = corrG.copy()
                        not_connected = set(corrG.nodes) - set(adjG_gc.nodes)
                        corrG_gc.remove_nodes_from(not_connected)
                        
                        dirc = '%s/data/06_network_measures/tm-%s/corr-%s/ds-%s/sub-%s'%(rootdir,tm,ct,ds,sub)
                        if gl=='local':
                            if c=='disconnected':
                                gc_local_dir   = glob.glob("%s/*local_measures_giant_component.csv"
                                                           %(dirc))
                                # Local Measures of the giant component:
                                gc_loc_mes_adj = pd.read_csv(gc_local_dir[0], sep='\t')
                                local_dir      = glob.glob("%s/*local_measures.csv"%(dirc))
                                # Local measures of the whole graph:
                                loc_mes_adj    = pd.read_csv(local_dir[0], sep='\t')

                                if mes_name in loc_mes_adj:
                                    print('measure already computed for subject %s network or its giant component with %s DS, %s CT and %s TM.'%(sub,ds,ct,tm))
                                    continue
                                del loc_mes_adj['Unnamed: 0']
                                del gc_loc_mes_adj['Unnamed: 0']
                                if bw == 'binary':
                                    loc_new_mes    = mes(adjG)
                                    gc_loc_new_mes = mes(adjG_gc)
                                if bw == 'weighted':
                                    loc_new_mes    = mes(corrG)
                                    gc_loc_new_mes = mes(corrG_gc)
                                loc_mes_adj[mes_name]    = loc_new_mes.values()
                                gc_loc_mes_adj[mes_name] = gc_loc_new_mes.values()
                                loc_mes_adj_norm         = normalize(loc_mes_adj)
                                loc_mes_adj_norm.ID      = loc_mes_adj.ID
                                gc_loc_mes_adj_norm      = normalize(gc_loc_mes_adj)
                                gc_loc_mes_adj_norm.ID   = gc_loc_mes_adj.ID

                                loc_mes_adj.to_csv \
                                ('%s/sub-%s_ds-%s_corr-%s_tm-%s_local_measures.csv'
                                 %(dirc,sub,ds,ct,tm), sep='\t')
                                gc_loc_mes_adj.to_csv \
                                ('%s/sub-%s_ds-%s_corr-%s_tm-%s_local_measures_giant_component.csv'
                                 %(dirc,sub,ds,ct,tm), sep='\t')
                                loc_mes_adj_norm.to_csv \
                                ('%s/sub-%s_ds-%s_corr-%s_tm-%s_local_measures_norm.csv'
                                 %(dirc,sub,ds,ct,tm), sep='\t')
                                gc_loc_mes_adj_norm.to_csv \
                                ('%s/sub-%s_ds-%s_corr-%s_tm-%s_local_measures_giant_component_norm.csv'
                                 %(dirc,sub,ds,ct,tm), sep='\t')

                            if c=='connected':
                                gc_local_dir    = glob.glob("%s/*local_measures_giant_component.csv"
                                                            %(dirc))
                                # Local Measures of the giant component:
                                gc_loc_mes_adj  = pd.read_csv(gc_local_dir[0], sep='\t')
                                if mes_name in gc_loc_mes_adj:
                                    print('measure already computed for subject %s with %s DS, %s CT and %s TM'
                                          %(sub,ds,ct,tm))
                                    continue
                                del gc_loc_mes_adj['Unnamed: 0']

                                if bw == 'binary':
                                    gc_loc_new_mes = mes(adjG_gc)
                                if bw == 'weighted':
                                    gc_loc_new_mes = mes(corrG_gc)

                                gc_loc_mes_adj[mes_name] = gc_loc_new_mes.values()
                                gc_loc_mes_adj_norm      = normalize(gc_loc_mes_adj)
                                gc_loc_mes_adj_norm.ID   = gc_loc_mes_adj.ID
                                gc_loc_mes_adj.to_csv \
                                ('%s/sub-%s_ds-%s_corr-%s_tm-%s_local_measures_giant_component.csv'
                                 %(dirc,sub,ds,ct,tm), sep='\t')
                                gc_loc_mes_adj_norm.to_csv \
                                ('%s/sub-%s_ds-%s_corr-%s_tm-%s_local_measures_giant_component_norm.csv'
                                 %(dirc,sub,ds,ct,tm), sep='\t')

                        elif gl=='global':
                            if c=='disconnected':
                                gc_global_dir = glob.glob("%s/*global_measures_giant_component.csv"
                                                          %(dirc))
                                global_dir    = glob.glob("%s/*global_measures.csv"%(dirc))
                                # Global measures of the whole graph:
                                glob_mes_adj    = pd.read_csv(global_dir[0], sep='\t')
                                # Global Measures of the giant component:
                                gc_glob_mes_adj = pd.read_csv(gc_global_dir[0], sep='\t')
                                if mes_name in gc_glob_mes_adj:
                                    print('measure already computed for subject %s network or its giant component with %s DS, %s CT and %s TM.'%(sub,ds,ct,tm))
                                    continue
                                del glob_mes_adj['Unnamed: 0']
                                del gc_glob_mes_adj['Unnamed: 0']

                                if bw == 'binary':
                                    glob_new_mes    = mes(adjG)
                                    gc_glob_new_mes = mes(adjG_gc)
                                if bw == 'weighted':
                                    glob_new_mes    = mes(corrG)
                                    gc_glob_new_mes = mes(corrG_gc)
                                glob_mes_adj[mes_name]    = glob_new_mes
                                gc_glob_mes_adj[mes_name] = gc_glob_new_mes

                                glob_mes_adj.to_csv \
                                ('%s/sub-%s_ds-%s_corr-%s_tm-%s_global_measures.csv'
                                 %(dirc,sub,ds,ct,tm), sep='\t')
                                gc_glob_mes_adj.to_csv \
                                ('%s/sub-%s_ds-%s_corr-%s_tm-%s_global_measures_giant_component.csv'
                                 %(dirc,sub,ds,ct,tm),sep='\t')

                            if c=='connected':
                                gc_global_dir = glob.glob("%s/*global_measures_giant_component.csv"
                                                          %(dirc))
                                # Global Measures of the giant component:
                                gc_glob_mes_adj = pd.read_csv(gc_global_dir[0], sep='\t')
                                if mes_name in gc_glob_mes_adj:
                                    print('measure already computed for subject %s with %s DS, %s CT and %s TM.'%(sub,ds,ct,tm))
                                    continue
                                del gc_glob_mes_adj['Unnamed: 0']

                                if bw == 'binary':
                                    gc_glob_new_mes = mes(adjG_gc)
                                if bw == 'weighted':
                                    gc_glob_new_mes = mes(corrG_gc)
                                gc_glob_mes_adj[mes_name] = gc_glob_new_mes
                                gc_glob_mes_adj.to_csv \
                                ('%s/sub-%s_ds-%s_corr-%s_tm-%s_global_measures_giant_component.csv'
                                 %(dirc,sub,ds,ct,tm), sep='\t')
