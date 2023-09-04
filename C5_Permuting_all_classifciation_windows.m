% This code shuffles contact labels and recalcualtes the classification
% performance (AUCs) 1000 times to generate the null distribution for statistical
% inference: this was done for all of the six classifications we performed including
% Within-patient classification (interictal and ictal)
% Across-patient generalisation (interictal and ictal)
% Across-time generalisation (interictal to ictal and vice versa)]

% INPUTS: classification data from C5_Permuting_all_classifciation_windows
% OUTPUTS: classification data + random data to be plotted by C6_Plotting_classification_and_feature_importance
%%
clc;
clear all;
close all;
%% 1 and 2: within-patient classification
ictal_or_inter='interictal';
load(['Within_subject_performance_sepnorm',ictal_or_inter,'_all_feats_crcted_feats_imp_100_10.mat'])
ff=1; % 1= all feaures included
iter_rand=1; % 1= number of iterations
for fld=1:size(Ground_truth,5) % folds
    for iter_equalis=1:size(Ground_truth,2) % iterations
        for p=1:size(Ground_truth,3) % patients
            Performance_within_inter(ff,iter_equalis,p,iter_rand,fld,1:7)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand,fld},Predictions{ff,iter_equalis,p,iter_rand,fld});
            [~,~,~,Performance_within_inter(ff,iter_equalis,p,iter_rand,fld,8)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand,fld},Predictions{ff,iter_equalis,p,iter_rand,fld},1);
            for rep=1:1000 % repititions
                Predictions_tmp=randsample(Predictions{ff,iter_equalis,p,iter_rand,fld},length(Predictions{ff,iter_equalis,p,iter_rand,fld}));
                Performance_rand_within_inter(ff,iter_equalis,p,iter_rand,fld,1:7,rep)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand,fld},Predictions_tmp);
                [~,~,~,Performance_rand_within_inter(ff,iter_equalis,p,iter_rand,fld,8,rep)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand,fld},Predictions_tmp,1);
            end
        end
    end
    [fld]
end
Performance_within_inter=squeeze(nanmean(nanmean(Performance_within_inter),5));
Performance_rand_within_inter=squeeze(nanmean(nanmean(Performance_rand_within_inter),5));
[1]

ictal_or_inter='ictal';
load(['Within_subject_performance_sepnorm',ictal_or_inter,'_all_feats_crcted_feats_imp_100_10.mat'])
ff=1; % 1= all feaures included
iter_rand=1; % 1= number of iterations
for fld=1:size(Ground_truth,5) % folds
    for iter_equalis=1:size(Ground_truth,2) % iterations
        for p=1:size(Ground_truth,3) % patients
            Performance_within_ictal(ff,iter_equalis,p,iter_rand,fld,1:7)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand,fld},Predictions{ff,iter_equalis,p,iter_rand,fld});
            [~,~,~,Performance_within_ictal(ff,iter_equalis,p,iter_rand,fld,8)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand,fld},Predictions{ff,iter_equalis,p,iter_rand,fld},1);
            for rep=1:1000 % repititions
                Predictions_tmp=randsample(Predictions{ff,iter_equalis,p,iter_rand,fld},length(Predictions{ff,iter_equalis,p,iter_rand,fld}));
                Performance_rand_within_ictal(ff,iter_equalis,p,iter_rand,fld,1:7,rep)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand,fld},Predictions_tmp);
                [~,~,~,Performance_rand_within_ictal(ff,iter_equalis,p,iter_rand,fld,8,rep)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand,fld},Predictions_tmp,1);
            end
        end
    end
    [fld]
end
Performance_within_ictal=squeeze(nanmean(nanmean(Performance_within_ictal),5));
Performance_rand_within_ictal=squeeze(nanmean(nanmean(Performance_rand_within_ictal),5));
[2]

%% 3 and 4: across-patient generalisation
ictal_or_inter='interictal';
load(['Generalisation_performance_across_subjects_',ictal_or_inter,'_all_feats_comb_crcted_feats_ovrsmp_imp_100_10.mat'])
classif=1; % 1= decision-tree classification
ff=1; % 1= all feaures included
iter_rand=1; % 1= number of iterations
for iter_equalis=1:size(Ground_truth,2) % iterations
    for p=1:size(Ground_truth,3) % patients
        Performance_acrss_subj_inter(ff,iter_equalis,p,iter_rand,1:7)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif});
        [~,~,~,Performance_acrss_subj_inter(ff,iter_equalis,p,iter_rand,8)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif},1);
        for rep=1:1000 % repititions
            Predictions_tmp=randsample(Predictions{ff,iter_equalis,p,iter_rand,classif},length(Predictions{ff,iter_equalis,p,iter_rand,classif}));
            Performance_rand_acrss_subj_inter(ff,iter_equalis,p,iter_rand,1:7,rep)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions_tmp);
            [~,~,~,Performance_rand_acrss_subj_inter(ff,iter_equalis,p,iter_rand,8,rep)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions_tmp,1);
        end
    end
end
Performance_acrss_subj_inter=squeeze(Performance_acrss_subj_inter);
Performance_rand_acrss_subj_inter=squeeze(Performance_rand_acrss_subj_inter);
[3]

ictal_or_inter='ictal';
load(['Generalisation_performance_across_subjects_',ictal_or_inter,'_all_feats_comb_crcted_feats_ovrsmp_imp_100_10.mat'])
classif=1; % 1= decision-tree classification
ff=1; % 1= all feaures included
iter_rand=1; % 1= number of iterations
for iter_equalis=1:size(Ground_truth,2) % iterations
    for p=1:size(Ground_truth,3) % patients
        Performance_acrss_subj_ictal(ff,iter_equalis,p,iter_rand,1:7)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif});
        [~,~,~,Performance_acrss_subj_ictal(ff,iter_equalis,p,iter_rand,8)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif},1);
        for rep=1:1000 % repititions
            Predictions_tmp=randsample(Predictions{ff,iter_equalis,p,iter_rand,classif},length(Predictions{ff,iter_equalis,p,iter_rand,classif}));
            Performance_rand_acrss_subj_ictal(ff,iter_equalis,p,iter_rand,1:7,rep)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions_tmp);
            [~,~,~,Performance_rand_acrss_subj_ictal(ff,iter_equalis,p,iter_rand,8,rep)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions_tmp,1);
        end
    end
end
Performance_acrss_subj_ictal=squeeze(Performance_acrss_subj_ictal);
Performance_rand_acrss_subj_ictal=squeeze(Performance_rand_acrss_subj_ictal);
[4]

%% 5 and 6: across-time generalisation
load(['Generalisation_performance_across_time_all_feats_comb_crcted_feats_ovrsmp_imp_100_10.mat'])
classif=1; % 1= decision-tree classification
ff=1; % 1= all feaures included
iter_rand=1; % 1= number of iterations
for iter_equalis=1:size(Ground_truth,2) % iterations
    for p=1:size(Ground_truth,3) % patients
        Performance_acrss_time_inter(ff,iter_equalis,p,iter_rand,1:7)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif});
        [~,~,~,Performance_acrss_time_inter(ff,iter_equalis,p,iter_rand,8)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif},1);
        for rep=1:1000 % repititions
            Predictions_tmp=randsample(Predictions{ff,iter_equalis,p,iter_rand,classif},length(Predictions{ff,iter_equalis,p,iter_rand,classif}));
            Performance_rand_acrss_time_inter(ff,iter_equalis,p,iter_rand,1:7,rep)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions_tmp);
            [~,~,~,Performance_rand_acrss_time_inter(ff,iter_equalis,p,iter_rand,8,rep)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions_tmp,1);
        end
    end
end
Performance_acrss_time_inter=squeeze(Performance_acrss_time_inter);
Performance_rand_acrss_time_inter=squeeze(Performance_rand_acrss_time_inter);
[5]

load(['Generalisation_performance_across_time_ict_to_int_all_feats_comb_crcted_feats_ovrsmp_imp_100_10.mat'])
classif=1; % 1= decision-tree classification
ff=1; % 1= all feaures included
iter_rand=1; % 1= number of iterations
for iter_equalis=1:size(Ground_truth,2) % iterations
    for p=1:size(Ground_truth,3) % patients
        Performance_acrss_time_ictal(ff,iter_equalis,p,iter_rand,1:7)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif});
        [~,~,~,Performance_acrss_time_ictal(ff,iter_equalis,p,iter_rand,8)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif},1);
        for rep=1:1000 % repititions
            Predictions_tmp=randsample(Predictions{ff,iter_equalis,p,iter_rand,classif},length(Predictions{ff,iter_equalis,p,iter_rand,classif}));
            Performance_rand_acrss_time_ictal(ff,iter_equalis,p,iter_rand,1:7,rep)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions_tmp);
            [~,~,~,Performance_rand_acrss_time_ictal(ff,iter_equalis,p,iter_rand,8,rep)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions_tmp,1);
        end
    end
end
Performance_acrss_time_ictal=squeeze(Performance_acrss_time_ictal);
Performance_rand_acrss_time_ictal=squeeze(Performance_rand_acrss_time_ictal);
[6]
%% Saving the data and the permuted null data for all 6 classifications
save('random_permutations.mat',...
    'Performance_acrss_subj_ictal','Performance_acrss_subj_inter',...
    'Performance_acrss_time_ictal','Performance_acrss_time_inter',...
    'Performance_rand_acrss_subj_ictal','Performance_rand_acrss_subj_inter',...
    'Performance_rand_acrss_time_ictal','Performance_rand_acrss_time_inter',...
    'Performance_rand_within_ictal','Performance_rand_within_inter',...
    'Performance_within_ictal','Performance_within_inter','-v7.3')
