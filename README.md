# Intracranial_EEG_generalisation

This repository shares the Matlab scripts used for the following manuscript:

Karimi-Rouzbahani, Hamid, McGonigal, Aileen, "Generalisability of epileptiform patterns across time and patients", MedRxiv, 2023.
doi: https://doi.org/10.1101/2023.08.29.23294708

This is the readme file for the scripts used in the study.

The analyis scripts are the ones whose file names start with Cx, x referring the order of use from epoching the data to plotting the final figures.
The scripts have coments to help with their understanding. But, please feel free to contact Hamid Karimi-Rouzbahani at hkarimi265@gmail.com if you have qny questions.
The other scripts, which do not start with a C, are feature extraction scripts, which are called within the "C2_Alldsets_Seizure_all_feature_extraction.m". Each of them
has been explained in the manuscript and cited in the following references:

Karimi-Rouzbahani, H., Shahmohammadi, M., Vahab, E., Setayeshi, S., & Carlson, T. (2021). Temporal variabilities provide additional category-related information in object
category decoding: a systematic comparison of informative EEG features. Neural Computation, 33(11), 3027-3072. https://doi.org/10.1162/neco_a_01436
&
Karimi-Rouzbahani, H., & Woolgar, A. (2022). When the whole is less than the sum of its parts: maximum object category information and behavioral prediction in multiscale
activation patterns. Frontiers in Neuroscience, 16, 825746. https://doi.org/10.3389/fnins.2022.825746


I provide below a short discription of each of the custom scripts used in the current manuscript - the scripts with C as their initial letter in file name. More details can be found within each script.

**C1_ds4100_ictal_interictal_epoching.m**
% This code loads the interictal and ictal recordings from the original dataset
% checks the channels; excludes bad channels, ECG/EKG channels and saves
% the channel information as well as signal epochs
% INPUTS: data from the original dataset
% OUTPUTS: epoched data to be used by C2_Alldsets_Seizure_all_feature_extraction 


**C2_Alldsets_Seizure_all_feature_extraction.m**
% This code extracts features from the interictal and ictal recordings 
% and saves the extracted features in new files
% Some features are extracted using external functions which should be
% available in paths accissible to Matlab
% INPUTS: epoched and cleaned data from C1_ds4100_ictal_interictal_epoching
% OUTPUTS: features data to be used by C3_Separating_target_non_target_contacts_all_feats



**C3_Separating_target_non_target_contacts_all_feats.m**
% This code separates contacts within and outside of the resection zone
% down-samples them into 28 samples and prepares them 
% for clssification in the next script and saves them
% INPUTS: features data from C2_Alldsets_Seizure_all_feature_extraction
% OUTPUTS: features data separated by contacts to be used by C4_Classification files


**C4_Classification_gen_acrss_patnt_inter_or_ictl_cmb_feat_ovrsmp.m**
% This code classifies resected and non-resected contacts, but does it by
% generalising across patients: trains the classifiers on all-minus-one
% patient and tests it on the left out patient using Decision Tree
% Classifiers
% INPUTS: features data separated by contacts from C3_Separating_target_non_target_contacts_all_feats
% OUTPUTS: classification data to be permuted by C5_Permuting_all_classifciation_windows


**C4_Classification_gen_acrss_time_ict_to_inter_cmb_feat_ovrsmp.m**
% This code classifies resected and non-resected contacts, but does it by
% generalising across time (interictal to ictal or vice versa):
% trains the classifiers on data from interictal (or vice versa) data and tests them on
% the data from the inctal window (or vice versa) using Decision Tree
% Classifiers
% INPUTS: features data separated by contacts from C3_Separating_target_non_target_contacts_all_feats
% OUTPUTS: classification data to be permuted by C5_Permuting_all_classifciation_windows


**C4_Classification_within_patnt_intr_or_ictl_cmb_feat_undrsmp.m**
% This code classifies resected and non-resected contacts within each patient
% using Decision Tree Classifiers
% INPUTS: features data separated by contacts from C3_Separating_target_non_target_contacts_all_feats
% OUTPUTS: classification data to be permuted by C5_Permuting_all_classifciation_windows


**C4_Classification_within_patnt_intr_or_ictl_cmb_feat_undrsmp_3_windows.m**
% This code classifies resected and non-resected contacts within each patient
% using Decision Tree Classifiers; It does it for an early, middle and late
% time window of data separately
% INPUTS: features data separated by contacts from C3_Separating_target_non_target_contacts_all_feats
% OUTPUTS: classification data to be permuted by C5_Permuting_3_classifciation_windows


**C5_Permuting_3_classifciation_windows.m**
% This code shuffles contact labels and recalcualtes the classification
% performance (AUCs) 1000 times to generate the null distribution for statistical
% inference: this was done for three sub-windows (each 2s) of the following
% classification only
% Within-patient classification (interictal and ictal)
% INPUTS: classification data from C5_Permuting_3_classifciation_windows
% OUTPUTS: classification data + random data to be plotted by C6_Plotting_classification_3_windows


**C5_Permuting_all_classifciation_windows.m**
% This code shuffles contact labels and recalcualtes the classification
% performance (AUCs) 1000 times to generate the null distribution for statistical
% inference: this was done for all of the six classifications we performed including
% Within-patient classification (interictal and ictal)
% Across-patient generalisation (interictal and ictal)
% Across-time generalisation (interictal to ictal and vice versa)]
% INPUTS: classification data from C5_Permuting_all_classifciation_windows
% OUTPUTS: classification data + random data to be plotted by C6_Plotting_classification_and_feature_importance


**C6_Plotting_classification_3_windows.m**
% This code plots the AUCs obtained from the classifications
% of contacts in each of the 3 time sub-windows in
% Within-patient classification (interictal and ictal)
% INPUTS: data from C5_Permuting_3_classifciation_windows
% OUTPUTS: figures and numbers produced in command window


**C6_Plotting_classification_and_feature_importance.m**
% This code plots the AUCs obtained from the classifications
% of contacts in 6 classifications including:
% Within-patient classification (interictal and ictal)
% Across-patient generalisation (interictal and ictal)
% Across-time generalisation (interictal to ictal and vice versa)
% INPUTS: data from C5_Permuting_all_classifciation_windows
% OUTPUTS: figures and numbers produced in command window


**C7_Plotting_correlations_across_feats_and_subjects.m**
% This code plots the correlation plots obtained from the classifications
% of contacts in 6 classifications including:
% Within-patient classification (interictal and ictal)
% Across-patient generalisation (interictal and ictal)
% Across-time generalisation (interictal to ictal and vice versa)
% INPUTS: data from C4_Classification files
% OUTPUTS: figures and numbers produced in command window

