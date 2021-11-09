#%matplotlib notebook
from config import *
from src.data import subjects
import os,sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import nilearn.plotting
import nilearn.datasets
import nilearn.signal
import nibabel
import sklearn.preprocessing
import pkg_resources


def denoise(sub):
    """
    """
    ses = 1
    #DATA_PATH = pkg_resources.resource_filename('brainnetworks', 'data/')
    from src.data import subjects
    subdir=subjects.preprocesseddir + '/sub-%s'%sub
    assert os.path.exists(subdir)

    sesdir=os.path.join(subdir,'ses-%d/func'%ses)

    # get freesurfer data if we don't already have it
    fsaverage = nilearn.datasets.fetch_surf_fsaverage(mesh='fsaverage5')
    print('processing sub-%s'%sub)
    # number of timepoints
    ld1 = os.path.join(sesdir,'sub-%s_ses-%d_task-rest_space-fsaverage5_hemi-%s.func.gii'%(sub,ses,'L'))
    ld2 = os.path.join(sesdir,'sub-%s_ses-%d_task-rest_space-fsaverage5_hemi-%s_bold.func.gii'%(sub,ses,'L'))

    try: 
        ntp=len(nibabel.load(ld1).darrays)
    except FileNotFoundError:
        ntp=len(nibabel.load(ld2).darrays)

    # load the preprocessed fMRI data

    hemispheres=['L','R']
    bold_origfile={}
    nverts=10242  # number of vertices in fsaverage5 per hemisphere
    bolddata_orig=np.zeros((ntp,nverts*2))
    for i,h in enumerate(hemispheres):

        try:
            bold_origfile[h]=os.path.join(sesdir,'sub-%s_ses-%d_task-rest_space-fsaverage5_hemi-%s.func.gii'%(sub,ses,h))
            d=nibabel.load(bold_origfile[h]).darrays
        except FileNotFoundError:
            bold_origfile[h]=os.path.join(sesdir,'sub-%s_ses-%d_task-rest_space-fsaverage5_hemi-%s_bold.func.gii'%(sub,ses,h))
            d=nibabel.load(bold_origfile[h]).darrays
        for tp in range(len(d)):
            bolddata_orig[tp,(i*nverts):((i+1)*nverts)]=d[tp].data

    print('data shape:',bolddata_orig.shape)
   

    # load the confound data

    confounds=pd.read_csv(os.path.join(sesdir,
                        'sub-%s_ses-%d_task-rest_desc-confounds_regressors.tsv'%(sub,ses)),
                         sep='\t',na_values='n/a')
    confounds=confounds.replace(np.nan,0)

    # add temporal derivatives of motion estimates
    # motionvars=['X', 'Y', 'Z', 'RotX', 'RotY', 'RotZ']
    # estimated head-motion parameters:
    motionvars=['trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z']

    for v in motionvars:
        confounds['%s_derivative1'%v]=0
        confounds['%s_derivative1'%v].iloc[1:]=confounds[v].iloc[1:].values - confounds[v].iloc[:-1].values
    
    #First we'll load in our data and check the shape
    bolddata = bolddata_orig
    try:
        a = confounds['non_steady_state_outlier00'].index[confounds['non_steady_state_outlier00'] == 1].tolist()[0]
        bolddata = bolddata_orig[a+1:,:]
        confounds = confounds.loc[a+1:]
    except:    
        pass
    ntp = confounds.shape[0]
    
    #Creatin Confounds tables:
    confounds9p = ['csf', 'white_matter', 'global_signal','trans_x', 'trans_y', 'trans_z',
                   'rot_x', 'rot_y', 'rot_z',]
    confounds36p = confounds9p + ['csf_derivative1', 'csf_derivative1_power2', 'csf_power2',
                                    'white_matter_derivative1', 'white_matter_derivative1_power2', 'white_matter_power2',
                                    'global_signal_derivative1', 'global_signal_derivative1_power2', 'global_signal_power2',
                                    'trans_x_derivative1', 'trans_x_derivative1_power2', 'trans_x_power2',
                                    'trans_y_derivative1', 'trans_y_derivative1_power2', 'trans_y_power2',
                                    'trans_z_derivative1', 'trans_z_derivative1_power2', 'trans_z_power2',
                                    'rot_x_derivative1', 'rot_x_derivative1_power2', 'rot_x_power2',
                                    'rot_y_derivative1', 'rot_y_derivative1_power2', 'rot_y_power2',
                                    'rot_z_derivative1', 'rot_z_derivative1_power2', 'rot_z_power2',
                                  ]
    confounds9p = confounds[confounds9p]
    confounds36p = confounds[confounds36p]
    
    #First create a matrix of zeroes that matches our signals matrix
    bolddata9p = np.zeros_like(bolddata)
    bolddata36p = np.zeros_like(bolddata)

    #Apply only to brain voxels
    bolddata9p = nilearn.signal.clean(bolddata, confounds=confounds9p.values)
    bolddata36p = nilearn.signal.clean(bolddata, confounds=confounds36p.values)
    
    atlasdir=rootdir + '/references/HCP-MMP1'
    atlas={'L':'lh.HCP-MMP1.fsaverage5.gii','R':'rh.HCP-MMP1.fsaverage5.gii'}
    atlasdata={}
    atlaslabels={}
    for a in atlas:
       atlaslabeltable=nibabel.load(os.path.join(atlasdir,atlas[a])).labeltable.labels
       atlaslabels[a]=[i.label for i in atlaslabeltable[1:]]
       atlasdata[a]=nibabel.load(os.path.join(atlasdir,atlas[a])).darrays[0].data 
    allatlaslabels=atlaslabels['L']+atlaslabels['R']
    allatlasdata=np.hstack((atlasdata['L'],atlasdata['R']+180))  
    
    roidata    = np.zeros((ntp,361))
    roidata9p  = np.zeros((ntp,361))
    roidata36p = np.zeros((ntp,361))

    for region in range(361):
        if region == 0:
            continue
        regionverts = allatlasdata == region
        for tp in range(ntp):
            tmp = bolddata[tp,:]
            roidata[tp,region] = np.mean(tmp[regionverts])
            tmp = bolddata9p[tp,:]
            roidata9p[tp,region] = np.mean(tmp[regionverts])
            tmp = bolddata36p[tp,:]
            roidata36p[tp,region] = np.mean(tmp[regionverts])
    fd_thresh=0.5
    tps_exceeding_fd_thresh=np.where(confounds.framewise_displacement.values>fd_thresh)
    tswindow=10
    tsmask=np.ones(confounds.shape[0])
    for tp in tps_exceeding_fd_thresh[0]:
        tsmask[(tp-1):(tp+tswindow)]=0

    print('%d good timepoints remaining after scrubbing (%d removed)'%(np.sum(tsmask),
                                                                       ntp - np.sum(tsmask)))
    roidata36p_scrubbed=roidata36p[np.where(tsmask)[0],:]
    os.system('mkdir -p %s/data/03_time_series/ds-9p %s/data/03_time_series/ds-36p \
    %s/data/03_time_series/ds-36pscrubbed %s/data/03_time_series/ds-raw'
              %(rootdir,rootdir,rootdir,rootdir))

    np.save("%s/data/03_time_series/ds-9p/sub-%s_ds-9p"
            %(rootdir, sub), roidata9p)
    np.save("%s/data/03_time_series/ds-36p/sub-%s_ds-36p"
            %(rootdir, sub), roidata36p)
    np.save("%s/data/03_time_series/ds-36pscrubbed/sub-%s_ds-36pscrubbed"
            %(rootdir, sub), roidata36p_scrubbed)
    np.save("%s/data/03_time_series/ds-raw/sub-%s_ds-raw"
            %(rootdir, sub), roidata)
    
def denoise_all(subs):
    """
    """
   
    for sub in subs:
        denoise(sub)