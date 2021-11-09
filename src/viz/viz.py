from config import *
from src.data import subjects
import numpy as np
import nibabel as nb
import pandas as pd
import glob,os
import random
import scipy.stats
from nilearn import surface,datasets,plotting
import matplotlib.pyplot as plt

coord = pd.read_csv('%s/references/HCP-MMP1/MMP_yeo2011_networks.csv'%rootdir)
coordinates = coord.T.loc[['X','Y','Z']].T.to_numpy()

colors=[(255/255, 0/255, 0/255), (249/255, 6/255, 6/255), (242/255, 13/255, 13/255), (236/255, 19/255, 19/255),
    (230/255, 25/255, 25/255), (223/255, 32/255, 32/255), (217/255, 38/255, 38/255), (210/255, 45/255, 45/255),
    (204/255, 51/255, 51/255), (198/255, 57/255, 57/255), (191/255, 64/255, 64/255), (185/255, 70/255, 70/255),
    (179/255, 77/255, 77/255), (172/255, 83/255, 83/255), (166/255, 89/255, 89/255), (159/255, 96/255, 96/255),
    (153/255, 102/255, 102/255), (147/255, 108/255, 108/255), (140/255, 115/255, 115/255),
    (134/255, 121/255, 121/255), (128/255, 128/255, 128/255)]
colors=colors[-1:0:-1]

def get_adjmtx(corrmtx,density,verbose=False):
    """
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

def brain_viz(groups, denoising_strategy,
              correlation_type,subs,density,
              method):
    """
    """
    global colors, coordinates
    from src.data import subjects
    all_subs = {}
    all_subs['all'] = subjects.subjects['all'].copy()
    for g in subjects_groups:
        all_subs[g] = all_subs['all'][all_subs['all'].group==g].participant_id.to_list()

    dirc           = rootdir + "/data/04_correlations/corr-%s/ds-%s"%(correlation_type,denoising_strategy)
    final_subjects = {}
    if type(subs) == int:
        subL={}
        for g in groups:
            subL[g]  = subjects.get_processed_list(dirc+'/*')
            subL[g]  = list(set(subL[g]) & set(all_subs[g]))
            try:
                final_subjects[g] = random.sample(subL[g], k=subs)
            except ValueError:
                print('there is not enough subjecsts available')

    if type(subs) == list:
        groups = (all_subs['all'][all_subs['all'].participant_id.isin(subs)].group).to_list()
        for g in groups:
            final_subjects[g]   = list(set(subs) & set(all_subs[g]))

    if method == '3d':
        for g in groups:
            for sub in final_subjects[g]:
                filesnp = glob.glob("%s/*%s*.npy"%(dirc,sub))
                corrM = np.load(filesnp[0])

                node_size = corrM.sum(axis=1)
                node_size-=np.mean(node_size)
                if np.std(node_size)!=0:
                    node_size/=np.std(node_size)
                node_size*=5.0
                node_size+=10
                view = plotting.view_connectome(get_adjmtx(corrM,density),
                                                coordinates,node_size=node_size,
                                                edge_threshold=.9, colorbar=False,
                                                title_fontsize=15,linewidth=4,
                                                title=
                                                'Top %.3f%% Edges, %s Subject: %s, DS: %s, CT: %s'
                                                %(100*density,g,sub,denoising_strategy,correlation_type))
                view.open_in_browser()
                view
    if method == '2d':
        for g in groups:
            for sub in final_subjects[g]:
                print('Top %.3f%% Edges,%s Subject: %s, DS: %s, CT: %s'
                      %(100*density,g,sub,denoising_strategy,correlation_type))

                filesnp = glob.glob("%s/*%s*.npy"%(dirc,sub))
                corrM = np.load(filesnp[0])

                node_size = corrM.sum(axis=1)
                node_size-=np.mean(node_size)
                if np.std(node_size)!=0:
                    node_size/=np.std(node_size)
                node_size*=15.0
                node_size+=36

                M = get_adjmtx(corrM,density)
                b=np.linspace(0,19,int(np.max(M.sum(axis=1)))+1,dtype=int)
                scolors = [colors[i] for i in b]
                sortedcolor = []
                for d in M.sum(axis=1):
                    sortedcolor.append(scolors[d])
                plt.figure(figsize=(len(b),1))
                for i in range(len(scolors)):
                    plt.plot([10 * i], [0],markersize=60,c = scolors[i], marker='s')
                    plt.text(10 * i, 0, i, horizontalalignment='center',
                         verticalalignment='center',color='white',fontsize=20)
                plt.axis('off')
                plt.show()

                view = plotting.plot_connectome(get_adjmtx(corrM,density),coordinates,edge_threshold=.9,
                    node_size=node_size,axes = (0, 0, 3, 3),edge_cmap='tab10',node_color = sortedcolor)
                view
    return


def brain_viz_from_path(title,density,method,outpath = "",path='',cm=[]):
    """
    """

    global coordinates

    if path !='':
        corrM = np.load(path)
    else:
        corrM = cm.copy()

    node_size = corrM.sum(axis=1)
    node_size-=np.mean(node_size)
    if np.std(node_size)!=0:
        node_size/=np.std(node_size)

    M = get_adjmtx(corrM,density)

    if method == '3d':
        node_size*=5.0
        node_size+=10
        if outpath!= "":
            view = plotting.view_connectome(M,
                                            coordinates,node_size=node_size,
                                            edge_threshold=.9, colorbar=False,
                                            title_fontsize=15,linewidth=4,
                                            title= title)
            # uncomment this to open the plot in a web browser:
            view.save_as_html(outpath)
            view
        else:
            view = plotting.view_connectome(M,
                                            coordinates,node_size=node_size,
                                            edge_threshold=.9, colorbar=False,
                                            title_fontsize=15,linewidth=4,
                                            title= title)
            view.open_in_browser()
            view
    if method == '2d':

        b=np.linspace(0,19,int(np.max(M.sum(axis=1)))+1,dtype=int)
        scolors = [colors[i] for i in b]
        sortedcolor = []
        for d in M.sum(axis=1):
            sortedcolor.append(scolors[d])
        node_size*=15.0
        node_size+=36
        if outpath!= "":
            print(title)
            view = plotting.plot_connectome(M,
                                            coordinates,edge_threshold=.9, node_size=node_size,
                                            axes = (0, 0, 3, 3),edge_cmap='tab10',node_color = sortedcolor,
                                            output_file = outpath)
            view
        else:

            plt.figure(figsize=(len(b),1))
            for i in range(len(scolors)):
                plt.plot([10 * i], [0],markersize=60,c = scolors[i], marker='s')
                plt.text(10 * i, 0, i, horizontalalignment='center',
                        verticalalignment='center',color='white',fontsize=20)
            plt.axis('off')
            plt.show()

            print(title)
            view = plotting.plot_connectome(M,
                                            coordinates,edge_threshold=.9, node_size=node_size,
                                            axes = (0, 0, 3, 3),edge_cmap='tab10',node_color = sortedcolor)
            view
    return


def brain_viz_nodes(nodes,title,method,node_size):
    """
    """
    global coordinates
    NEWcoordinates = coordinates[nodes,:]
    node_size = np.array(node_size).astype('float32')
    node_size-=np.mean(node_size)
    if np.std(node_size)!=0:
        node_size/=np.std(node_size)

    if method == '2d':
        print(title)
        node_size*=15.0
        node_size+=36
        view=plotting.plot_connectome(np.zeros((len(nodes),len(nodes))),
            NEWcoordinates,edge_threshold=0, node_size=node_size,
            node_color = 'red',axes = (0, 0, 3, 3))
        view
    elif method == '3d':
        node_size*=5.0
        node_size+=10
        view=plotting.view_markers(NEWcoordinates,title = title ,marker_size=node_size)
        view.open_in_browser()
        view
    return
