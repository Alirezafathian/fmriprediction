from config import *
from src.data import subjects
import numpy as np
import nibabel as nb
import glob,os
import random
import scipy.stats
from nilearn import surface,datasets,plotting

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

    from src.data import subjects
    all_subs = {}
    all_subs['all'] = subjects.subjects['all'].copy()
    for g in subjects_groups:
        all_subs[g] = all_subs['all'][all_subs['all'].group==g].participant_id.to_list()

    atlasdir=rootdir + '/references/HCP-MMP1'
    atlas={'left':'lh.HCP-MMP1.fsaverage5.gii','right':'rh.HCP-MMP1.fsaverage5.gii'}
    fsaverage = datasets.fetch_surf_fsaverage()

    atlaslabels = []
    coordinates = []

    for hemi in ['left', 'right']:
        atlaslabeltable=nb.load(os.path.join(atlasdir,atlas[hemi])).labeltable.labels
        atlaslabels.extend([i.label for i in atlaslabeltable[1:]])

        vert =nb.load(os.path.join(atlasdir,atlas[hemi])).darrays[0].data
        rr, _ = surface.load_surf_mesh(fsaverage['pial_%s' % hemi])
        for k, label in enumerate(atlaslabels):
            if "Unknown" not in str(label):  # Omit the Unknown label.
                # Compute mean location of vertices in label of index k
                coordinates.append(np.mean(rr[vert == k], axis=0))

    coordinates = np.array(coordinates)  # 3D coordinates of parcels
    coordinates = coordinates[~np.isnan(coordinates).any(axis=1)]
    # droping the first roi
    coordinates = coordinates[1:]
    dirc          = rootdir + "/data/04_correlations/corr-%s/ds-%s"%(correlation_type,denoising_strategy)
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
                view = plotting.view_connectome(get_adjmtx(corrM,density),
                                                coordinates,node_size=6,
                                                edge_threshold=.9, colorbar=False,
                                                title_fontsize=15,linewidth=4,
                                                title=
                                                'Top %.3f%% Edges, %s Subject: %s, DS: %s, CT: %s'
                                                %(100*density,g,sub,denoising_strategy,correlation_type))
                # uncomment this to open the plot in a web browser:
                view.open_in_browser()
                view
    if method == '2d':
        for g in groups:
            for sub in final_subjects[g]:
                print('Top %.3f%% Edges,%s Subject: %s, DS: %s, CT: %s'
                      %(100*density,g,sub,denoising_strategy,correlation_type))

                filesnp = glob.glob("%s/*%s*.npy"%(dirc,sub))
                corrM = np.load(filesnp[0])
                view = plotting.plot_connectome(get_adjmtx(corrM,density),
                                                coordinates,edge_threshold=.9, node_size=30
                                                )
                                                #output_file =
                                                #'HCP-MMP_atlas_sub-%s_%s.pdf'
                                                #%(sub,denoising_strategy),
                                                #)

                # uncomment this to open the plot in a web browser:
                view
    return


def brain_viz_from_path(path,title,density,method,outpath = ""):
    """
    """

    atlasdir=rootdir + '/references/HCP-MMP1'
    atlas={'left':'lh.HCP-MMP1.fsaverage5.gii','right':'rh.HCP-MMP1.fsaverage5.gii'}
    fsaverage = datasets.fetch_surf_fsaverage()

    atlaslabels = []
    coordinates = []

    for hemi in ['left', 'right']:
        atlaslabeltable=nb.load(os.path.join(atlasdir,atlas[hemi])).labeltable.labels
        atlaslabels.extend([i.label for i in atlaslabeltable[1:]])

        vert =nb.load(os.path.join(atlasdir,atlas[hemi])).darrays[0].data
        rr, _ = surface.load_surf_mesh(fsaverage['pial_%s' % hemi])
        for k, label in enumerate(atlaslabels):
            if "Unknown" not in str(label):  # Omit the Unknown label.
                # Compute mean location of vertices in label of index k
                coordinates.append(np.mean(rr[vert == k], axis=0))

    coordinates = np.array(coordinates)  # 3D coordinates of parcels
    coordinates = coordinates[~np.isnan(coordinates).any(axis=1)]
    # droping the first roi
    coordinates = coordinates[1:]
    corrM = np.load(path)
    if method == '3d':
        if outpath!= "":
            view = plotting.view_connectome(get_adjmtx(corrM,density),
                                            coordinates,node_size=6,
                                            edge_threshold=.9, colorbar=False,
                                            title_fontsize=15,linewidth=4,
                                            title= title)
            # uncomment this to open the plot in a web browser:
            view.open_in_browser()
            view
        else:
            view = plotting.view_connectome(get_adjmtx(corrM,density),
                                            coordinates,node_size=6,
                                            edge_threshold=.9, colorbar=False,
                                            title_fontsize=15,linewidth=4,
                                            title= title)
            view.save_as_html()
            view
    if method == '2d':
        if outpath!= "":
            print(title)
            view = plotting.plot_connectome(get_adjmtx(corrM,density),
                                            coordinates,edge_threshold=.9, node_size=30,
                                            output_file = outpath)
            view

        else:
            print(title)
            view = plotting.plot_connectome(get_adjmtx(corrM,density),
                                            coordinates,edge_threshold=.9, node_size=30)
            view

    return
