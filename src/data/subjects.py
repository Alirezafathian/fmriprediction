from config import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
from IPython.display import display
import re
#########################################################################################
# study subjects:
subjects = {}
subjects['all'] = pd.read_csv(subjects_list_dir, sep=',')
for sg in subjects_groups:
    subjects[sg] = subjects['all']['group']==sg
    subjects[sg] = subjects['all'][subjects[sg]]
    subjects[sg] = list(subjects[sg].participant_id)
#########################################################################################
#directories
BIDsdir           = rootdir + "/data/01_bids"
preprocesseddir   = rootdir + "/data/02_fmriprep/fmriprep"
time_seriesdir    = rootdir + "/data/03_time_series"
denoiseddir       = rootdir + "/data/04_correlations"
adjdir            = rootdir + "/data/05_adjacency_matrices"
#########################################################################################

def get_processed_list(directory):
    processedfiles = glob.glob(directory)
    subs = [i.split('sub-') for i in processedfiles]
    subs = [i for i in subs if len(i)!=1 ]
    subs = [i[1] for i in subs]
    subs = [i.split('_')[0] for i in subs]
    subs = [i.split('.')[0] for i in subs]
    subs = list(set(subs))
    return subs

DICOMs            = get_processed_list(dicom_dir+"/*")
BIDs              = get_processed_list(BIDsdir+'/*')
preprocessed      = get_processed_list(preprocesseddir+"/*.html")
denoised          = get_processed_list(time_seriesdir+"/*/*")
corr_constructed  = get_processed_list(denoiseddir+"/*/*/*")
adj_constructed   = get_processed_list(adjdir+"/*/*/*/*")
#########################################################################################
not_processed = []
processed     = []
not_complete  = []
complete      = []
for sub in adj_constructed:
    c = 0
    tmp1 = []
    tmp2 = []

    for ds in denoising_strategies:
        for ct in correlation_types:
            for tm in thresholding_methods:
                for tv in thresholding_values:
                    if tm=='userdefined':
                        tm = '%s-%.3f'%(tm,tv)
                    dirc = '%s/data/06_network_measures/tm-%s/corr-%s/ds-%s/sub-%s'%(rootdir,tm,ct,ds,sub)
                    files = glob.glob(dirc+'/*')
                    #print(files)
                    if len(files)==6:
                        c+=1
                        tmp2.append([sub,ds,tm])
                    else:
                        tmp1.append([sub,ds,tm])
    #print(c,len(files))
    if 'userdefined' in thresholding_methods:
        k = len(thresholding_methods) + len(thresholding_values) - 1
    else:
        k = len(thresholding_methods)
    l = len(denoising_strategies)*len(correlation_types)*k
    if c==l:
        processed.append(sub)
    elif c==0:
        not_processed.append(sub)
    else:
        not_complete.extend(tmp1)
        complete.extend(tmp2)
not_processed = list(set(not_processed))
processed = list(set(processed))
#########################################################################################
to_download               = list(set(subjects['all'].participant_id) - set(DICOMs) - set(BIDs))
to_convert_to_BIDs        = list(set(DICOMs) - set(BIDs))
to_preprocess             = list(set(BIDs) - set(preprocessed))
to_denoise                = list(set(preprocessed) - set(denoised))
to_construct_correlation  = list(set(denoised) - set(corr_constructed))
to_construct_adj          = list(set(corr_constructed) - set(adj_constructed))
to_compute_measures       = [not_processed,not_complete]
to_single_sub_analysis    = [processed,complete]
to_group_level_analysis   = [processed,complete]
#########################################################################################
results = """
=========================================================================================\n
List of subjects to be downloaded:\n%s
=========================================================================================\n
List of subjects to be converted into BIDs format:\n%s
=========================================================================================\n
List of subjects to be preprocessed:\n%s
=========================================================================================\n
List of subjects to be denoised:\n%s
=========================================================================================\n
List of subjects to construct correlation matrix:\n%s
=========================================================================================\n
List of subjects to construct adjacency matrix:\n%s
=========================================================================================\n
List of subjects to compute network measures:\n%s
=========================================================================================\n
List of subjects to do single subject network analysis:\n%s
=========================================================================================\n
List of subjects ready for group level analysis:\n%s
\nNumber of subjects ready for group level analysis: %s
=========================================================================================\n
"""%(to_download,to_convert_to_BIDs,to_preprocess,to_denoise,to_construct_correlation,
     to_construct_adj,to_compute_measures,to_group_level_analysis,to_group_level_analysis,
     len(to_group_level_analysis[0]))
#########################################################################################

def print_subs():
    print(results)
if __name__ == "__main__":
    print_subs( )
