from config import rootdir,subjects_groups,denoising_strategies,normalize_measures
from src.data import subjects
from src.data.materials import *
import src.group_level_analysis.group_level_analysis as gla
from src.viz import viz
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob,os
from IPython.display import display
import seaborn as sns
import random
import hoggorm as ho
import networkx as nx
from matplotlib.collections import LineCollection
import itertools
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib import cm
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import scipy.stats
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.axes._axes import _log as matplotlib_axes_logger
from src.adjacency_matrix.adjacency_matrix import fGCE



def general(ds,ct,negative_corr,tm,tv,enodes_name,selected_set,nb,binM='load'):
    global subjects_groups,denoising_strategies,normalize_measures,rootdir
    #all_sub_list = subjects.to_group_level_analysis[0]
    all_sub_list = list(subjects.subjects['all']['participant_id'].values)
    all_subjects = subjects.subjects['all'][subjects.subjects['all']['participant_id'].isin(all_sub_list)]
#===============================================================
    subjects_list = {}
    for g in subjects_groups:
        subjects_list[g] = all_subjects.where(all_subjects.group==g).dropna()
        print('\n%d %s subjects are included'%(len(subjects_list[g]),g))
        display(subjects_list[g])
        subjects_list[g] = subjects_list[g].participant_id.tolist()
#===============================================================
    group_couples = [(a,b) for a in subjects_groups for b in
                     [i for i in subjects_groups[subjects_groups.index(a)+1:]]]
#===============================================================
    if tm=='userdefined':
        tm = '%s-%.3f'%(tm,tv)
    sgn = 'positive'
    if negative_corr:
        sgn = 'negative'
#===============================================================
    enodes = list(range(360))
#===============================================================

#===============================================================
    sg = subjects_groups.copy()
    sg.append('all')

#===============================================================
    str_loc    = ['loc_'  + ds for ds in denoising_strategies]
    str_glob   = ['glob_' + ds for ds in denoising_strategies]
    strlist    = str_loc.copy()
    strlist.extend(str_glob)

#===============================================================
    norm = ""
    if normalize_measures:
        norm="_norm"
#===============================================================
#===============================================================
#===============================================================
#===============================================================
    cols = subjects_list[subjects_groups[0]].copy()
    sep=[len(subjects_list[subjects_groups[0]])]
    for i in range(1,len(subjects_groups)):
        cols.extend(subjects_list[subjects_groups[i]])
        sep.append(len(subjects_list[subjects_groups[i]]) + sep[i-1])
    sep.insert(0,0)
#===============================================================
    subL = []
    subLRAW = []
    removed_subL = []
    if enodes_name!= 'all':
        removed_subL = [i for i in range(360) if i not in enodes]
    for s in cols:
        M = gla.read_M(s,ds,ct)
        #for i in removed_subL:
            #M[i,:]= np.zeros(np.shape(M)[0])
            #M = M.T
            #M[i,:]= np.zeros(np.shape(M)[0])
            #M = M.T
        subLRAW.append(M)
        M = M[enodes,:][:,enodes]
        subL.append(M - np.mean(M))

#===============================================================

#===============================================================
    print(removed_subL)
    os.system("mkdir -p %s/data/08_group_level/mean_network/tm-%s/corr-%s/ds-%s/"%(rootdir,tm,ct,ds))
    meanM = {}
    meanG = {}
    meanM_bin = {}
    meanG_bin = {}
    for g in range(len(subjects_groups)):
        con  = np.dstack(subLRAW[sep[g]:sep[g+1]])
        meanM[subjects_groups[g]] = np.mean(con, axis=2)
        for i in removed_subL:
            for j in removed_subL:
                meanM[subjects_groups[g]][i][j] = 0.0
        np.fill_diagonal(meanM[subjects_groups[g]], 0.0)
        meanG[subjects_groups[g]] = nx.from_numpy_matrix(meanM[subjects_groups[g]])
        if enodes_name=='all' and nb=='grouplevel':
            nx.write_gexf(meanG[subjects_groups[g]],"%s/data/08_group_level/mean_network/tm-%s/corr-%s/ds-%s/%s_%s_ds-%s_corr-%s_tm-%s.gexf"
                          %(rootdir,tm,ct,ds,subjects_groups[g],'corr',ds,ct,tm))
            np.save("%s/data/08_group_level/mean_network/tm-%s/corr-%s/ds-%s/%s_%s_ds-%s_corr-%s_tm-%s"
                          %(rootdir,tm,ct,ds,subjects_groups[g],'corr',ds,ct,tm), meanM[subjects_groups[g]])
        if binM == 'compute':
            if negative_corr:
                meanG_bin[subjects_groups[g]] = fGCE(nx.from_numpy_array(-meanM[subjects_groups[g]].clip(max=0)), thresholds = np.arange(0,1,.1),tol=0.01)
            else:
                meanG_bin[subjects_groups[g]] = fGCE(nx.from_numpy_array(meanM[subjects_groups[g]].clip(min=0)), thresholds = np.arange(0,1,.1),tol=0.01)
            meanM_bin[subjects_groups[g]] = nx.to_numpy_matrix(meanG_bin[subjects_groups[g]])
        if binM == 'load':
            meanM_bin[subjects_groups[g]] = np.load("%s/data/08_group_level/mean_network/tm-%s/corr-%s/ds-%s/%s_%s_ds-%s_bin-%s-corr-%s_tm-%s.npy"
                                                    %(rootdir,tm,ct,ds,subjects_groups[g],'corr',ds,sgn,ct,tm))
            meanG_bin[subjects_groups[g]] = nx.from_numpy_array(meanM_bin[subjects_groups[g]])
        if enodes_name=='all' and nb=='grouplevel':
            nx.write_gexf(meanG_bin[subjects_groups[g]],"%s/data/08_group_level/mean_network/tm-%s/corr-%s/ds-%s/%s_%s_ds-%s_bin-%s-corr-%s_tm-%s.gexf"
                          %(rootdir,tm,ct,ds,subjects_groups[g],'corr',ds,sgn,ct,tm))
            np.save("%s/data/08_group_level/mean_network/tm-%s/corr-%s/ds-%s/%s_%s_ds-%s_bin-%s-corr-%s_tm-%s"
                          %(rootdir,tm,ct,ds,subjects_groups[g],'corr',ds,sgn,ct,tm), meanM_bin[subjects_groups[g]])
    con = np.dstack(subLRAW[0:sep[-1]])
    meanM['all'] = np.mean(con, axis=2)
    for i in removed_subL:
        for j in removed_subL:
            meanM['all'][i][j] = 0.0
    np.fill_diagonal(meanM['all'], 0.0)
    meanG['all'] = nx.from_numpy_matrix(meanM['all'])
    if enodes_name=='all' and nb=='grouplevel':
        nx.write_gexf(meanG['all'],"%s/data/08_group_level/mean_network/tm-%s/corr-%s/ds-%s/%s_%s_ds-%s_corr-%s_tm-%s.gexf"
                      %(rootdir,tm,ct,ds,'all','corr',ds,ct,tm))
        np.save("%s/data/08_group_level/mean_network/tm-%s/corr-%s/ds-%s/%s_%s_ds-%s_corr-%s_tm-%s"
                %(rootdir,tm,ct,ds,'all','corr',ds,ct,tm), meanM['all'])

    if binM == 'compute':
        if negative_corr:
            mgt = fGCE(nx.from_numpy_array(-meanM['all'].clip(max=0)), thresholds = np.arange(0,1,.1),tol=0.01)
            meanG_bin['all'] = nx.from_numpy_array(-nx.to_numpy_array(mgt))
        else:
            meanG_bin['all'] = fGCE(nx.from_numpy_array(meanM['all'].clip(min=0)), thresholds = np.arange(0,1,.1),tol=0.01)
        meanM_bin['all'] = nx.to_numpy_matrix(meanG_bin['all'])
    if binM == 'load':
        if negative_corr:
            meanM_bin['all'] = -np.load("%s/data/08_group_level/mean_network/tm-%s/corr-%s/ds-%s/%s_%s_ds-%s_bin-%s-corr-%s_tm-%s.npy"
                                                    %(rootdir,tm,ct,ds,'all','corr',ds,sgn,ct,tm))
            meanG_bin['all'] = nx.from_numpy_array(meanM_bin['all'])
        else:
            meanM_bin['all'] = np.load("%s/data/08_group_level/mean_network/tm-%s/corr-%s/ds-%s/%s_%s_ds-%s_bin-%s-corr-%s_tm-%s.npy"
                                                    %(rootdir,tm,ct,ds,'all','corr',ds,sgn,ct,tm))
            meanG_bin['all'] = nx.read_gexf("%s/data/08_group_level/mean_network/tm-%s/corr-%s/ds-%s/%s_%s_ds-%s_bin-%s-corr-%s_tm-%s.gexf"
                                                         %(rootdir,tm,ct,ds,'all','corr',ds,sgn,ct,tm))

    if enodes_name=='all' and nb=='grouplevel':
        nx.write_gexf(meanG_bin['all'],"%s/data/08_group_level/mean_network/tm-%s/corr-%s/ds-%s/%s_%s_ds-%s_bin-%s-corr-%s_tm-%s.gexf"
                      %(rootdir,tm,ct,ds,'all','corr',ds,sgn,ct,tm))
        np.save("%s/data/08_group_level/mean_network/tm-%s/corr-%s/ds-%s/%s_%s_ds-%s_bin-%s-corr-%s_tm-%s"
                %(rootdir,tm,ct,ds,'all','corr',ds,sgn,ct,tm), meanM_bin['all'])

#===============================================================
    return group_couples,sep,all_sub_list,all_subjects,subjects_list,sg,str_loc,str_glob,strlist,tm,sgn,cols,subL,subLRAW,meanM,meanG,meanM_bin,meanG_bin
#============================================================================================================
coord = pd.read_csv('%s/references/HCP-MMP1/MMP_yeo2011_networks.csv'%rootdir)
coordinates = coord.T.loc[['X','Y','Z']].T.to_numpy()

hullL = ConvexHull(coordinates[0:180,0:2])
hullR = ConvexHull(coordinates[180:360,0:2])
hull = ConvexHull(coordinates[:,0:2])
#============================================================================================================
def get_adjmtx(corrmtx,density,verbose=False):
    """
    """
    assert density<=1
    cutoff=scipy.stats.scoreatpercentile(corrmtx[np.triu_indices_from(corrmtx,1)],
                                         100-(100*density))

    if verbose:
        print('cutoff:%0.3f'%cutoff)
    adjmtx=(corrmtx>cutoff).astype('int')
    adjmtx[np.diag_indices_from(adjmtx)]=0
    return(adjmtx)

#===============================================================
colors= [(255/255, 184/255, 188/255),(255/255,136/255,144/255),(255/255,89/255,99/255),(200/255,0,13/255)]
def multilayerplot(ii,rnk,sp,mes,
    snd,sep,all_sub_list,all_subjects,subjects_list,sg,str_loc,str_glob,
    strlist,tm,sgn,measures,RAWmeasures,intersect,removed,removed_subL,r,cols,subL,subLRAW,
    subL_adj,subLRAW_adj,meanM,meanG,meanM_bin,meanG_bin,meanM_adj,meanG_adj):
    for i in ii:
        fig = plt.figure(figsize=(8,6))
        matplotlib_axes_logger.setLevel('ERROR')
        ax = fig.add_subplot(111, projection='3d')
        M = get_adjmtx(meanM[subjects_groups[i]],density=.005)

        for j in range(len(sp[subjects_groups[i]])):
            rk = np.insert(rnk[subjects_groups[i]][str(sp[subjects_groups[i]][j])],r,
                           np.zeros(len(r)))
            cd1 =coordinates[rk>0]
            x1 =cd1[:,0]
            y1 =cd1[:,1]
            z1 =np.zeros(len(y1))+4*j
            ax.scatter(x1, y1,z1,s=10, marker='o',c=(0,19/255,29/255))

            cd2 =coordinates[rk==0]
            x0 =cd2[:,0]
            y0 =cd2[:,1]
            z0 =np.zeros(len(y0))+4*j
            ax.scatter(x0, y0,z0,s=4, marker='o',c=(0,19/255,29/255), alpha=0.3)

            for simplex in hull.simplices:
                plt.plot(coordinates[:,0:2][simplex, 0],
                         coordinates[:,0:2][simplex, 1],
                         np.zeros(2)+4*j,'k-',linewidth=1,
                         c=(85/255,136/255,163/255),alpha=0.5)
            for simplex in hullL.simplices:
                plt.plot(coordinates[0:180,0:2][simplex, 0],
                         coordinates[0:180,0:2][simplex, 1],
                         np.zeros(2)+4*j, 'k-',linewidth=1,
                         c=(85/255,136/255,163/255))#(20/255,83/255,116/255)
            for simplex in hullR.simplices:
                plt.plot(coordinates[180:360,0:2][simplex, 0],
                         coordinates[180:360,0:2][simplex, 1],
                         np.zeros(2)+4*j, 'k-',linewidth=1,
                         c=(85/255,136/255,163/255))

            ll =np.shape(coordinates[:,0:2][hull.vertices,0])[0]
            x = list(coordinates[:,0:2][hull.vertices,0])
            y = list(coordinates[:,0:2][hull.vertices,1])
            z = list(np.zeros(ll)+4*j)
            verts = [list(zip(x,y,z))]
            llL =np.shape(coordinates[0:180,0:2][hullL.vertices,0])[0]
            xL = list(coordinates[0:180,0:2][hullL.vertices,0])
            yL = list(coordinates[0:180,0:2][hullL.vertices,1])
            zL = list(np.zeros(llL)+4*j)
            vertsL = [list(zip(xL,yL,zL))]
            llR =np.shape(coordinates[180:360,0:2][hullR.vertices,0])[0]
            xR = list(coordinates[180:360,0:2][hullR.vertices,0])
            yR = list(coordinates[180:360,0:2][hullR.vertices,1])
            zR = list(np.zeros(llR)+4*j)
            vertsR = [list(zip(xR,yR,zR))]
            ax.add_collection3d(Poly3DCollection(verts,alpha=.1,
                                                 facecolors=(85/255,136/255,163/255)))
            ax.add_collection3d(Poly3DCollection(vertsL,alpha=.1,
                                                 facecolors=(85/255,136/255,163/255)))
            ax.add_collection3d(Poly3DCollection(vertsR,alpha=.1,
                                                 facecolors=(85/255,136/255,163/255)))

            M1 =M[rk>0,:][:,rk>0]
            inter_corr = []
            for ii in range(len(cd1)):
                for jj in range(ii+1,len(cd1)):
                    if M1[ii][jj] >0.0:
                        inter_corr.append([(cd1[ii,0],cd1[ii,1],4*j),
                                           (cd1[jj,0],cd1[jj,1],4*j)])
            ax.add_collection(Line3DCollection(inter_corr, colors=colors[j],
                                               linewidths=(j+1)/3))

        for a in range(len(sp[subjects_groups[i]])-1):
            for b in range(a+1,len(sp[subjects_groups[i]])):
                intra_corr = []
                rka = np.insert(rnk[subjects_groups[i]][str(sp[subjects_groups[i]][a])],r,
                                np.zeros(len(r)))
                rkb = np.insert(rnk[subjects_groups[i]][str(sp[subjects_groups[i]][b])],r,
                                np.zeros(len(r)))
                cd1 = coordinates[rka>0]
                cd2 = coordinates[rkb>0]
                M1 =M[rka>0,:][:,rkb>0]
                for aa in range(len(cd1)):
                    for bb in range(len(cd2)):
                        if M1[aa][bb] >0.0:
                            intra_corr.append([(cd1[aa,0],cd1[aa,1],4*a),
                                               (cd2[bb,0],cd2[bb,1],4*b)])
                ax.add_collection(Line3DCollection(intra_corr, colors=colors[max(a,b)],
                                                   linewidths=(max(a,b)+1)/3))
        ax.text(70, 0, 0, "Low %s nodes"%mes, color=(0,19/255,29/255))
        ax.text(70, 0, 4*(len(sp[subjects_groups[i]])-1), "High %s nodes"%mes, color=(0,19/255,29/255))

        plt.axis('off')
        matplotlib_axes_logger.setLevel('ERROR')
        plt.show()

#===============================================================
colors= [(255/255, 184/255, 188/255),(255/255,136/255,144/255),(255/255,89/255,99/255),(200/255,0,13/255)]
def twolayerplot(ii, rnk,sp,mes,
    snd,sep,all_sub_list,all_subjects,subjects_list,sg,str_loc,str_glob,
    strlist,tm,sgn,measures,RAWmeasures,intersect,removed,removed_subL,r,cols,subL,subLRAW,
    subL_adj,subLRAW_adj,meanM,meanG,meanM_bin,meanG_bin,meanM_adj,meanG_adj):
    global colors
    for i in ii:
        fig = plt.figure(figsize=(10,10))
        matplotlib_axes_logger.setLevel('ERROR')
        ax = fig.add_subplot(111, projection='3d')
        M = get_adjmtx(meanM_bin[subjects_groups[i]],density=.01)
        #M = meanM_bin[subjects_groups[i]]
        for j in [0,1]:
            #for simplex in hull.simplices:
                #plt.plot(coordinates[:,0:2][simplex, 0],
                #         coordinates[:,0:2][simplex, 1],
                #         np.zeros(2)+2*j-0.1,'k-',linewidth=1,
                #         c=(85/255,136/255,163/255),alpha=0.5)
            for simplex in hullL.simplices:
                plt.plot(coordinates[0:180,0:2][simplex, 0],
                         coordinates[0:180,0:2][simplex, 1],
                         np.zeros(2)+2*j-0.1, 'k-',linewidth=1,
                         c=(85/255,136/255,163/255))#(20/255,83/255,116/255)
            for simplex in hullR.simplices:
                plt.plot(coordinates[180:360,0:2][simplex, 0],
                         coordinates[180:360,0:2][simplex, 1],
                         np.zeros(2)+2*j-0.1, 'k-',linewidth=1,
                         c=(85/255,136/255,163/255))

            ll =np.shape(coordinates[:,0:2][hull.vertices,0])[0]
            x = list(coordinates[:,0:2][hull.vertices,0])
            y = list(coordinates[:,0:2][hull.vertices,1])
            z = list(np.zeros(ll)+2*j-0.1)
            verts = [list(zip(x,y,z))]
            llL =np.shape(coordinates[0:180,0:2][hullL.vertices,0])[0]
            xL = list(coordinates[0:180,0:2][hullL.vertices,0])
            yL = list(coordinates[0:180,0:2][hullL.vertices,1])
            zL = list(np.zeros(llL)+2*j-0.1)
            vertsL = [list(zip(xL,yL,zL))]
            llR =np.shape(coordinates[180:360,0:2][hullR.vertices,0])[0]
            xR = list(coordinates[180:360,0:2][hullR.vertices,0])
            yR = list(coordinates[180:360,0:2][hullR.vertices,1])
            zR = list(np.zeros(llR)+2*j-0.1)
            vertsR = [list(zip(xR,yR,zR))]
            #ax.add_collection3d(Poly3DCollection(verts,alpha=.3,
            #                                     facecolors=(85/255,136/255,163/255)))
            ax.add_collection3d(Poly3DCollection(vertsL,alpha=.5,
                                                 facecolors=(85/255,136/255,163/255)))
            ax.add_collection3d(Poly3DCollection(vertsR,alpha=.5,
                                                 facecolors=(85/255,136/255,163/255)))


        """
        Scaling is done from here...
        """
        x_scale=411*9/519
        y_scale=9
        z_scale=3

        scale=np.diag([x_scale, y_scale, z_scale, 1.0])
        scale=scale*(1.0/scale.max())
        scale[3,3]=1.0

        def short_proj():
          return np.dot(Axes3D.get_proj(ax), scale)

        ax.get_proj=short_proj
        """
        to here
        """

        for a in range(len(sp[subjects_groups[i]])-1):
            for b in range(a+1,len(sp[subjects_groups[i]])):
                intra_corr = []
                rka = np.insert(rnk[subjects_groups[i]][str(sp[subjects_groups[i]][a])],r,
                                np.zeros(len(r)))
                rkb = np.insert(rnk[subjects_groups[i]][str(sp[subjects_groups[i]][b])],r,
                                np.zeros(len(r)))
                cd1 = coordinates[rka>0]
                cd2 = coordinates[rkb>0]
                M1 =M[rka>0,:][:,rkb>0]
                for aa in range(len(cd1)):
                    for bb in range(len(cd2)):
                        if M1[aa][bb] >0.0:
                            if a in [0]:
                                if b in [0]:
                                    intra_corr.append([(cd1[aa,0],cd1[aa,1],2*a),
                                                       (cd2[bb,0],cd2[bb,1],2*b)])
                                    colo = colors[2]
                                    lw=.8
                                else:
                                    intra_corr.append([(cd1[aa,0],cd1[aa,1],2*a),
                                                       (cd2[bb,0],cd2[bb,1],2+.25*(b-1))])
                                    colo = colors[3]
                                    lw=1
                            else:
                                if b in [0]:
                                    intra_corr.append([(cd1[aa,0],cd1[aa,1],2+.25*(a-1)),
                                                       (cd2[bb,0],cd2[bb,1],2*b)])
                                else:
                                    intra_corr.append([(cd1[aa,0],cd1[aa,1],2+.25*(a-1)),
                                                       (cd2[bb,0],cd2[bb,1],2+.25*(b-1))])
                                colo = colors[3]
                                lw=1
                ax.add_collection(Line3DCollection(intra_corr, colors=colo,linewidths=lw))
        for j in range(len(sp[subjects_groups[i]])):
            rk = np.insert(rnk[subjects_groups[i]][str(sp[subjects_groups[i]][j])],r,
                           np.zeros(len(r)))
            cd1 =coordinates[rk>0]
            M1 =M[rk>0,:][:,rk>0]
            inter_corr = []
            for ii in range(len(cd1)):
                for jj in range(ii+1,len(cd1)):
                    if M1[ii][jj] >0.0:
                        if j in [0]:
                            inter_corr.append([(cd1[ii,0],cd1[ii,1],2*j),
                                               (cd1[jj,0],cd1[jj,1],2*j)])
                            colo = colors[2]
                            lw=.8
                        else:
                            inter_corr.append([(cd1[ii,0],cd1[ii,1],2+.25*(j-1)),
                                               (cd1[jj,0],cd1[jj,1],2+.25*(j-1))])
                            colo = colors[3]
                            lw=.8
            ax.add_collection(Line3DCollection(inter_corr, colors=colo,linewidths=lw))

            x1 =cd1[:,0]
            y1 =cd1[:,1]
            z1 =np.zeros(len(y1))+2+.25*(j-1)
            cc = (0,0/255,0/255)
            if j == 0:
                cc = (36/255,36/255,36/255)
            if j in [0,1]:
                z1 =np.zeros(len(y1))+2*j
            ax.scatter(x1, y1,z1,s=15+2*(3**j), marker='o',c=cc,
                       depthshade=False,alpha=1)

            if j in [0,1]:
                cd2 =coordinates[rk==0]
                x0 =cd2[:,0]
                y0 =cd2[:,1]
                z0 =np.zeros(len(y0))+2*j
                ax.scatter(x0, y0,z0,s=0, marker='o',c=(0,0/255,0/255), alpha=0.3)


        ax.text(-120, 0, -.5, "Low %s nodes"%mes, color=(0,19/255,29/255))
        ax.text(-120, 0, 2.5, "High %s nodes"%mes, color=(0,19/255,29/255))

        plt.axis('off')
        matplotlib_axes_logger.setLevel('ERROR')
        plt.show()

#===============================================================
def combrnk(rnk1,rnk2,subjects_groups):
    rnk= {'all':{}}
    sp = {}
    for i in range(len(subjects_groups)):
        rnk[subjects_groups[i]]={}
        for j in range(4):
            rnk[subjects_groups[i]][str(j)] = rnk1[subjects_groups[i]][list(rnk1[subjects_groups[i]].keys())[j]]+rnk2[subjects_groups[i]][list(rnk2[subjects_groups[i]].keys())[j]]
        for j in range(4):
            for k in range(j+1,4):
                nz = np.where(rnk[subjects_groups[i]][str(k)]!=0)[0]
                for w in nz:
                    rnk[subjects_groups[i]][str(j)][w]= 0.0
        sp[subjects_groups[i]] = list(range(4))
    return rnk,sp
#===============================================================
def mes_func(G,mes):
    if mes=='betweenness':
        return np.fromiter(nx.betweenness_centrality(G).values(), dtype=float)
    if mes=='degree':
        return [G.degree(j) for j in G.nodes]
    if mes=='strength':
        dic = dict(G.degree(weight= 'weight'))
        return np.fromiter(dic.values(), dtype=float)
    if mes=='clustering_coef':
        return np.fromiter(nx.clustering(G).values(), dtype=float)
    if mes=='eccentricity':
        return np.fromiter(nx.eccentricity(G).values(), dtype=float)
    if mes=='eigenvector':
        return np.fromiter(nx.eigenvector_centrality(G,max_iter=500).values(),
                           dtype=float)

#===============================================================
def frnk(method,mes,
    snd,sep,all_sub_list,all_subjects,subjects_list,sg,str_loc,str_glob,
    strlist,tm,sgn,measures,RAWmeasures,intersect,removed,removed_subL,r,cols,subL,subLRAW,
    subL_adj,subLRAW_adj,meanM,meanG,meanM_bin,meanG_bin,meanM_adj,meanG_adj,
    sp1={}):
    fig = plt.figure()
    rnk = {'all':{}}
    c= ['#9bdeac','#e1d97e','#e2979c','#e7305b']
    sp=sp1.copy()
    for i in range(len(subjects_groups)):
        ax = fig.add_subplot(2,2,i+1)
        mesM = []
        for s in cols[sep[i]:sep[i+1]]:
            ax.plot(np.arange(len(intersect['loc_36p'])),np.sort(RAWmeasures['all'][str_loc[0]][s].T.loc[mes].values),c=c[i])
            try:
                mesM = np.vstack((mesM,RAWmeasures['all'][str_loc[0]][s].T.loc[mes].values))
            except ValueError:
                mesM = RAWmeasures['all'][str_loc[0]][s].T.loc[mes].values
        if method == 'mean measure':
            rnk['all'][subjects_groups[i]] = np.mean(mesM,axis=0)
        plt.title(subjects_groups[i])
        if method == 'mean network':
            rnk['all'][subjects_groups[i]] = mes_func(meanG_bin[subjects_groups[i]],mes)
            rnk['all'][subjects_groups[i]] = np.delete(rnk['all'][subjects_groups[i]], removed)
        if sp1 == {}:
            sp[subjects_groups[i]]=[0]
            for k in [1,1.5,2]:
                try:
                    msd = np.mean(rnk['all'][subjects_groups[i]]) + k * np.std(rnk['all'][subjects_groups[i]])
                    sp[subjects_groups[i]].append(np.where(np.sort(rnk['all'][subjects_groups[i]])>=msd)[0][0])
                except IndexError:
                    sp[subjects_groups[i]].append(intersect['loc_36p'][-1]-1)
                    continue
        for ss in sp[subjects_groups[i]]:
            plt.axvline(x=ss)

    plt.tight_layout(pad=.3)
    plt.show()

    fig = plt.figure()
    for i in range(len(subjects_groups)):
        ax = fig.add_subplot(2,2,i+1)
        plt.hist(rnk['all'][subjects_groups[i]],color=c[i])
        plt.title(subjects_groups[i])
    plt.tight_layout(pad=.3)
    plt.show()

    fig = plt.figure(figsize=(4,3))
    ax = fig.add_subplot(1,1,1)

    for i in range(len(subjects_groups)):
        ax.plot(np.arange(len(rnk['all'][subjects_groups[i]])),np.sort(rnk['all'][subjects_groups[i]]),c=c[i])
    plt.legend(subjects_groups)
    plt.show()

    pd.DataFrame.from_dict(sp)

    for i in range(len(subjects_groups)):
        rnk[subjects_groups[i]] = {}
        for j in range(len(sp[subjects_groups[i]])):
            try:
                rnk[subjects_groups[i]][str(sp[subjects_groups[i]][j])]=rnk['all'][subjects_groups[i]]*(
                    (rnk['all'][subjects_groups[i]]>=
                     np.sort(rnk['all'][subjects_groups[i]])[sp[subjects_groups[i]][j]]).astype(int))*(
                    (rnk['all'][subjects_groups[i]]<
                     np.sort(rnk['all'][subjects_groups[i]])[sp[subjects_groups[i]][j+1]]).astype(int))
            except IndexError:
                rnk[subjects_groups[i]][str(sp[subjects_groups[i]][j])]=rnk['all'][subjects_groups[i]]*(
                    (rnk['all'][subjects_groups[i]]>=
                     np.sort(rnk['all'][subjects_groups[i]])[sp[subjects_groups[i]][j]]).astype(int))*(
                    (rnk['all'][subjects_groups[i]]<360).astype(int))

    return rnk,sp
