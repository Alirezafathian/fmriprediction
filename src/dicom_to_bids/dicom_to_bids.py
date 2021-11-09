from config import *
import glob, os
from src.data import subjects
import pandas as pd


def dicom_to_bids(to_convert_to_BIDs):
    
    global dicom_dir
    participants_file = open("%s/data/01_bids/participants.tsv"%(rootdir),"a+")
    for sub in subjects.to_convert_to_BIDs:
        if os.path.isdir('%s/data/01_bids/sub-%s/ses-1/anat'%(rootdir,sub)) and\
    len(os.listdir('%s/data/01_bids/sub-%s/ses-1/anat'%(rootdir,sub))) != 0:
            print('sub-%s: already covnverted'%sub)
            continue
        dirc = dicom_dir + '/sub-'+sub
        if not os.path.isdir(dirc):
            print('sub-%s: could not find DICOM files'%sub)
            continue
        SubFiles = glob.glob(dirc + '/*')
        if len(SubFiles) < 6:
            os.system('find %s -maxdepth 1 -type f -exec rm -f {} \;'%dirc)
            os.system('dcm2niix -f %%f_%%p_%%t_%%S -p y -z y -ba n %s'%dirc)
        elif len(SubFiles) > 6:
            print('sub-%s: dcm2niix generated more than 4 files. conversion failed'%sub)
            continue

        mprage_json = glob.glob("%s/*MPRAGE*json"%(dirc))[0]
        mprage_nii  = glob.glob("%s/*MPRAGE*gz"%(dirc))[0]
        fmri_json   = glob.glob("%s/*json"%(dirc))
        fmri_nii    = glob.glob("%s/*gz"%(dirc))
        fmri_json   = list(set(fmri_json) - {mprage_json})[0]
        fmri_nii    = list(set(fmri_nii)  - {mprage_nii})[0]
        d           = rootdir + '/data/01_bids/sub-' + sub
        print('sub-%s: converted'%sub)
        participants_file.write("sub-%s\r\n" %(sub))
        os.system('mkdir -p %s/ses-1/anat;mkdir %s/ses-1/func'%(d,d))
        os.system('cp \'%s\' %s/ses-1/func/sub-%s_ses-1_task-rest_bold.nii.gz'%(fmri_nii,d,sub))
        os.system('cp \'%s\' %s/ses-1/func/sub-%s_ses-1_task-rest_bold.json'%(fmri_json,d,sub))
        os.system('cp \'%s\' %s/ses-1/anat/sub-%s_ses-1_T1w.nii.gz'%(mprage_nii,d,sub))
        os.system('cp \'%s\' %s/ses-1/anat/sub-%s_ses-1_T1w.json'%(mprage_json,d,sub))
    participants_file.close()
    
    p = pd.read_csv("%s/data/01_bids/participants.tsv"%(rootdir), sep='\t')
    p=p.drop_duplicates()
    p.to_csv("%s/data/01_bids/participants.tsv"%(rootdir),index=False, sep='\t')
