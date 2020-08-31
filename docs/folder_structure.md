Folder Structure of the Preoject
==============================
```
.
├── config.py
├── data
│   ├── 00_dicom
│   │   ├── sub-002S4746
│   │   │   ├── MPRAGE
│   │   │   │   └── 2013-06-03_15_00_48.0
│   │   │   │       └── S191323
│   │   │   │           ├── ADNI_002_S_4746_MR_MPRAGE_br_raw_20130604142850776_2_S191323_I375098.dcm
│   │   │   │           ├── ADNI_002_S_4746_MR_MPRAGE_br_raw_20130604142853366_151_S191323_I375098.dcm
│   │   │   │           ├── ...
│   │   │   │           └── ADNI_002_S_4746_MR_MPRAGE_br_raw_20130604143657963_104_S191323_I375098.dcm
│   │   │   └── Resting_State_fMRI
│   │   │       └── 2013-06-03_15_00_48.0
│   │   │           └── S191322
│   │   │               ├── ADNI_002_S_4746_MR_Resting_State_fMRI_br_raw_.dcm
│   │   │               ├── ADNI_002_S_4746_MR_Resting_State_fMRI_br_raw.dcm
│   │   │               ├── ...
│   │   │               └── ADNI_037_S_4214_MR_Axial_rsfMRI__Eyes_Open__br_raw.dcm
│   │   └── ...
│   ├── 01_bids
│   │   ├── dataset_description.json
│   │   ├── participants.tsv
│   │   ├── sub-002S4171
│   │   │   └── ses-1
│   │   │       ├── anat
│   │   │       │   ├── sub-002S4171_ses-1_T1w.json
│   │   │       │   └── sub-002S4171_ses-1_T1w.nii.gz
│   │   │       └── func
│   │   │           ├── sub-002S4171_ses-1_task-rest_bold.json
│   │   │           └── sub-002S4171_ses-1_task-rest_bold.nii.gz
│   │   └── task-rest_bold.json
│   ├── 02_fmriprep
│   │   ├── fmriprep
│   │   │   ├── dataset_description.json
│   │   │   ├── desc-aparcaseg_dseg.tsv
│   │   │   ├── desc-aseg_dseg.tsv
│   │   │   ├── logs
│   │   │   │   └── CITATION.md
│   │   │   ├── sub-002S4171
│   │   │   │   ├── anat
│   │   │   │   │   └── ...
│   │   │   │   ├── figures
│   │   │   │   │   └── ...
│   │   │   │   └── ses-1
│   │   │   │       ├── figures
│   │   │   │       │   └── ...
│   │   │   │       └── func
│   │   │   │           ├── ...
│   │   │   │           ├── sub-002S4171_ses-1_task-rest_desc-confounds_regressors.tsv
│   │   │   │           ├── sub-002S4171_ses-1_task-rest_space-fsaverage5_hemi-L.func.gii
│   │   │   │           ├── sub-002S4171_ses-1_task-rest_space-fsaverage5_hemi-R.func.gii
│   │   │   │           └── ...
│   │   │   ├── sub-002S4171.html
│   │   │   └── ...
│   │   └── freesurfer
│   │       └── ...
│   ├── 03_time_series
│   │   ├── ds-36p
│   │   │   ├── sub-002S4171_ds-36p.npy
│   │   │   └── ...
│   │   ├── ds-36pscrubbed
│   │   │   ├── sub-002S4171_ds-36pscrubbed.npy
│   │   │   └── ...
│   │   ├── ds-9p
│   │   │   ├── sub-002S4171_ds-9p.npy
│   │   │   └── ...
│   │   └── ds-raw
│   │       ├── sub-002S4251_ds-raw.npy
│   │       └── ...
│   ├── 04_correlations
│   │   ├── corr-glasso
│   │   │   ├── ds-36p
│   │   │   ├── ds-36pscrubbed
│   │   │   └── ds-9p
│   │   └── corr-pearson
│   │       ├── ds-36p
│   │       │   ├── sub-002S4171_ds-36p_corr-pearson.gexf
│   │       │   ├── sub-002S4171_ds-36p_corr-pearson.npy
│   │       │   └── ...
│   │       ├── ds-36pscrubbed
│   │       │   ├── sub-002S4171_ds-36pscrubbed_corr-pearson.gexf
│   │       │   ├── sub-002S4171_ds-36pscrubbed_corr-pearson.npy
│   │       │   └── ...
│   │       └── ds-9p
│   │           ├── sub-002S4171_ds-9p_corr-pearson.gexf
│   │           ├── sub-002S4171_ds-9p_corr-pearson.npy
│   │           └── ...
│   ├── 05_adjacency_matrices
│   │   ├── tm-gce
│   │   │   ├── corr-glasso
│   │   │   │   ├── ds-36p
│   │   │   │   ├── ds-36pscrubbed
│   │   │   │   └── ds-9p
│   │   │   └── corr-pearson
│   │   │       ├── ds-36p
│   │   │       │   ├── sub-002S4171_ds-36p_corr-pearson_tm-gce.gexf
│   │   │       │   ├── sub-002S4171_ds-36p_corr-pearson_tm-gce.npy
│   │   │       │   └── ...
│   │   │       ├── ds-36pscrubbed
│   │   │       │   ├── sub-002S4171_ds-36pscrubbed_corr-pearson_tm-gce.gexf
│   │   │       │   ├── sub-002S4171_ds-36pscrubbed_corr-pearson_tm-gce.npy
│   │   │       │   └── ...
│   │   │       └── ds-9p
│   │   │           ├── sub-002S4171_ds-9p_corr-pearson_tm-gce.gexf
│   │   │           ├── sub-002S4171_ds-9p_corr-pearson_tm-gce.npy
│   │   │           └── ...
│   │   └── tm-userdefined-0.050
│   │       └── corr-pearson
│   │           ├── ds-36p
│   │           │   ├── sub-002S4171_ds-36p_corr-pearson_tm-userdefined_density-0.050.gexf
│   │           │   ├── sub-002S4171_ds-36p_corr-pearson_tm-userdefined_density-0.050.npy
│   │           │   └── ...
│   │           ├── ds-36pscrubbed
│   │           │   ├── sub-002S4171_ds-36pscrubbed_corr-pearson_tm-userdefined_density-0.050.gexf
│   │           │   ├── sub-002S4171_ds-36pscrubbed_corr-pearson_tm-userdefined_density-0.050.npy
│   │           │   └── ...
│   │           └── ds-9p
│   │               ├── sub-002S4171_ds-9p_corr-pearson_tm-userdefined_density-0.050.gexf
│   │               ├── sub-002S4171_ds-9p_corr-pearson_tm-userdefined_density-0.050.npy
│   │               └── ...
│   ├── 06_network_measures
│   │   ├── tm-gce
│   │   │   ├── corr-glasso
│   │   │   │   ├── ds-36p
│   │   │   │   ├── ds-36pscrubbed
│   │   │   │   └── ds-9p
│   │   │   └── corr-pearson
│   │   │       ├── ds-36p
│   │   │       │   ├── sub-002S4171
│   │   │       │   │   ├── sub-002S4171_ds-36p_corr-pearson_tm-gce_global_measures.csv
│   │   │       │   │   ├── sub-002S4171_ds-36p_corr-pearson_tm-gce_global_measures_giant_component.csv
│   │   │       │   │   ├── sub-002S4171_ds-36p_corr-pearson_tm-gce_local_measures.csv
│   │   │       │   │   ├── sub-002S4171_ds-36p_corr-pearson_tm-gce_local_measures_giant_component.csv
│   │   │       │   │   ├── sub-002S4171_ds-36p_corr-pearson_tm-gce_local_measures_giant_component_norm.csv
│   │   │       │   │   └── sub-002S4171_ds-36p_corr-pearson_tm-gce_local_measures_norm.csv
│   │   │       │   └── ...
│   │   │       ├── ds-36pscrubbed
│   │   │       │   ├── sub-002S4171
│   │   │       │   │   └── ...
│   │   │       │   └── ...
│   │   │       └── ds-9p
│   │   │           ├── sub-002S4171
│   │   │           │   └── ...
│   │   │           └── ...
│   │   └── tm-userdefined-0.050
│   │       ├── corr-glasso
│   │       │   └── ...
│   │       └── corr-pearson
│   │           ├── ds-36p
│   │           │   ├── sub-002S4171
│   │           │   │   └── ...
│   │           │   └── ...
│   │           ├── ds-36pscrubbed
│   │           │   ├── sub-002S4171
│   │           │   │   └── ...
│   │           └── ds-9p
│   │               ├── sub-002S4171
│   │               │   └── ...
│   │               └── ...
│   ├── 07_group_level
│   │   ├── mantelC.npy
│   │   └── mean_network
│   │       ├── ...
│   │       └── viz
│   │           └── ...
│   ├── 08_features
│   │   └── fsisher_fsf
│   │       ├── fsisher_fsf_features.csv
│   │       └── targets.csv
│   └── subjects_list.csv
├── docs
│   └── pipeline
│       └── ...
├── notebooks
│   ├── 00_subjects.ipynb
│   ├── 01_dicom_to_bids.ipynb
│   ├── 02_fMRIPrep_preprocessing.ipynb
│   ├── 03_denoising_and_extracting_time-Series.ipynb
│   ├── 04_estimating_functional_connectivity.ipynb
│   ├── 05_adjacency_matrix.ipynb
│   ├── 06_computing_network_measures.ipynb
│   ├── 07_single_subject_network_analysis.ipynb
│   ├── 08_group_level_network_analysis.ipynb
│   ├── 09_feature_selection.ipynb
│   ├── 10_classification.ipynb
│   └── 11_visualiser.ipynb
├── README.md
├── references
│   ├── fsaverage5
│   │   ├── __init__.py
│   │   ├── pial_inflated.left.gii
│   │   ├── pial_inflated.left.gii.gz
│   │   ├── pial_inflated.right.gii
│   │   ├── pial_inflated.right.gii.gz
│   │   ├── pial.left.gii
│   │   ├── pial.left.gii.gz
│   │   ├── pial.right.gii
│   │   ├── pial.right.gii.gz
│   │   ├── sulc.left.gii
│   │   ├── sulc.left.gii.gz
│   │   ├── sulc.right.gii
│   │   └── sulc.right.gii.gz
│   ├── FSlicense
│   │   └── license.txt
│   ├── graphics
│   │   └── ...
│   └──  HCP-MMP1
│       ├── lh.HCP-MMP1.annot
│       ├── lh.HCP-MMP1.fsaverage5.gii
│       ├── MMP_yeo2011_networks.csv
│       ├── rh.HCP-MMP1.annot
│       └── rh.HCP-MMP1.fsaverage5.gii
├── reports
│   └── figures
├── requirements.txt
└── src
    ├── adjacency_matrix
    │   ├── adjacency_matrix.py
    │   └── __init__.py
    ├── data
    │   ├── __init__.py
    │   └── subjects.py
    ├── denoising_timeseries
    │   ├── denoising_timeseries.py
    │   └── __init__.py
    ├── dicom_to_bids
    │   ├── dicom_to_bids.py
    │   └── __init__.py
    ├── features
    │   └── __init__.py
    ├── fMRIPrep_preprocessing
    │   ├── fMRIPrep_preprocessing.sh
    │   └── __init__.py
    ├── functional_connectivity
    │   ├── functional_connectivity.py
    │   └── __init__.py
    ├── group_level_analysis
    │   ├── group_level_analysis.py
    │   └── __init__.py
    ├── __init__.py
    ├── models
    │   └── __init__.py
    ├── network_measures
    │   ├── __init__.py
    │   └── network_measures.py
    ├── single_subject_analysis
    │   ├── __init__.py
    │   └── single_subject_analysis.py
    └── viz
        ├── __init__.py
        └── viz.py
```
