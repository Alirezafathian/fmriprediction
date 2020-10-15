subjects_groups      = ['CN','EMCI','LMCI','AD'] # CN:   control,
					         # EMCI: early mild cognitive impairment,
						 # LMCI: late mild cognitive impairment,
						 # AD:   Alzheimer’s disease
denoising_strategies = ['36p'] # or '36pscrubbed', '9p'
correlation_types    = ['pearson'] # 'pearson' or 'glasso' or 'partial'
thresholding_methods = ['gce'] # or 'userdefined'
thresholding_values  = [0.05] # density of the adjacency matrix with userdefined thresholding method.
normalize_measures   = True   # normalize local measures by mean
negative_corr        = True # if True it will compute adjacency matrix based on negative correlation as well as adjacency matrix based on positive correlation
#directories:
rootdir              = '/home/alireza/Thesis/fmriprediction' # the directory of fmriprediction
subjects_list_dir    = rootdir + '/data/subjects_list.csv' 
dicom_dir            = rootdir + '/data/00_dicom'
