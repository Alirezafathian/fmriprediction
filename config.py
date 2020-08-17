subjects_groups      = ['CN','EMCI','LMCI','AD']
denoising_strategies = ['9p', '36p', '36pscrubbed']
correlation_types    = ['pearson'] #'glasso'
thresholding_methods = ['gce', 'userdefined']
thresholding_values  = [0.05] # density of the adjacency matrix with userdefined thresholding method.
normalize_measures   = True   # normalize local measures by mean

#directories:
rootdir              = '/home/alireza/Thesis/fmriprediction'
subjects_list_dir    = rootdir + '/data/subjects_list.csv'
dicom_dir            = rootdir + '/data/00_dicom'