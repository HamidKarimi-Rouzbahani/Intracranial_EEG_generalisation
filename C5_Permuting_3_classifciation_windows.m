% This code shuffles contact labels and recalcualtes the classification
% performance (AUCs) 1000 times to generate the null distribution for statistical
% inference: this was done for three sub-windows (each 2s) of the following
% classification only
% Within-patient classification (interictal and ictal)

% INPUTS: classification data from C5_Permuting_3_classifciation_windows
% OUTPUTS: classification data + random data to be plotted by C6_Plotting_classification_3_windows
%%
clc;
clear all;
close all;
%%
% Interictal
ictal_or_inter='interictal';
load(['Within_subject_performance_sepnorm_1_',ictal_or_inter,'_all_feats_crcted_feats_imp_100_10.mat'])
ff=1; % 1= all feaures included
iter_rand=1; % 1= number of iterations
for tw=1:3 % time windows
    for fld=1:size(Ground_truth,6) % folds
        for iter_equalis=1:size(Ground_truth,3) % iterations
            for p=1:size(Ground_truth,4) % patients
                Performance_within_inter(tw,ff,iter_equalis,p,iter_rand,fld,1:7)=Evaluate(Ground_truth{tw,ff,iter_equalis,p,iter_rand,fld},Predictions{tw,ff,iter_equalis,p,iter_rand,fld});
                [~,~,~,Performance_within_inter(tw,ff,iter_equalis,p,iter_rand,fld,8)]=perfcurve(Ground_truth{tw,ff,iter_equalis,p,iter_rand,fld},Predictions{tw,ff,iter_equalis,p,iter_rand,fld},1);
                for rep=1:1000 % repititions
                    Predictions_tmp=randsample(Predictions{tw,ff,iter_equalis,p,iter_rand,fld},length(Predictions{tw,ff,iter_equalis,p,iter_rand,fld}));
                    Performance_rand_within_inter(tw,ff,iter_equalis,p,iter_rand,fld,1:7,rep)=Evaluate(Ground_truth{tw,ff,iter_equalis,p,iter_rand,fld},Predictions_tmp);
                    [~,~,~,Performance_rand_within_inter(tw,ff,iter_equalis,p,iter_rand,fld,8,rep)]=perfcurve(Ground_truth{tw,ff,iter_equalis,p,iter_rand,fld},Predictions_tmp,1);
                end
            end
        end
        [tw fld]
    end
end
Performance_within_inter=squeeze(nanmean(nanmean(Performance_within_inter,2),6));
Performance_rand_within_inter=squeeze(nanmean(nanmean(Performance_rand_within_inter,2),6));
[1]

% Ictal
ictal_or_inter='ictal';
load(['Within_subject_performance_sepnorm_1_',ictal_or_inter,'_all_feats_crcted_feats_imp_100_10.mat'])
ff=1; % 1= all feaures included
iter_rand=1; % 1= number of iterations
for tw=1:3 % time windows
    for fld=1:size(Ground_truth,6) % folds
        for iter_equalis=1:size(Ground_truth,3) % iterations
            for p=1:size(Ground_truth,4) % patients
                Performance_within_ictal(tw,ff,iter_equalis,p,iter_rand,fld,1:7)=Evaluate(Ground_truth{tw,ff,iter_equalis,p,iter_rand,fld},Predictions{tw,ff,iter_equalis,p,iter_rand,fld});
                [~,~,~,Performance_within_ictal(tw,ff,iter_equalis,p,iter_rand,fld,8)]=perfcurve(Ground_truth{tw,ff,iter_equalis,p,iter_rand,fld},Predictions{tw,ff,iter_equalis,p,iter_rand,fld},1);
                for rep=1:1000 % repititions
                    Predictions_tmp=randsample(Predictions{tw,ff,iter_equalis,p,iter_rand,fld},length(Predictions{tw,ff,iter_equalis,p,iter_rand,fld}));
                    Performance_rand_within_ictal(tw,ff,iter_equalis,p,iter_rand,fld,1:7,rep)=Evaluate(Ground_truth{tw,ff,iter_equalis,p,iter_rand,fld},Predictions_tmp);
                    [~,~,~,Performance_rand_within_ictal(tw,ff,iter_equalis,p,iter_rand,fld,8,rep)]=perfcurve(Ground_truth{tw,ff,iter_equalis,p,iter_rand,fld},Predictions_tmp,1);
                end
            end
        end
        [tw fld]
    end
end
Performance_within_ictal=squeeze(nanmean(nanmean(Performance_within_ictal,2),6));
Performance_rand_within_ictal=squeeze(nanmean(nanmean(Performance_rand_within_ictal,2),6));
[2]

save('random_permutations_within_time_wind.mat',...
    'Performance_within_ictal','Performance_rand_within_ictal',...
    'Performance_within_inter','Performance_rand_within_inter',...
    '-v7.3')
