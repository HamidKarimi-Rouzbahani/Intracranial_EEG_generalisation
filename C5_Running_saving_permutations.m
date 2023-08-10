clc;
clear all;
close all;
ictal_or_inter='interictal';
load(['Within_subject_performance_',ictal_or_inter,'_all_feats_crcted_feats_imp_100_10.mat'])
classif=1; %
ff=1;
iter_rand=1;
for fld=1:size(Ground_truth,5)    
    for iter_equalis=1:size(Ground_truth,2)
        for p=1:size(Ground_truth,3)
            Performance_within_inter(ff,iter_equalis,p,iter_rand,fld,1:7)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand,fld},Predictions{ff,iter_equalis,p,iter_rand,fld});
            [~,~,~,Performance_within_inter(ff,iter_equalis,p,iter_rand,fld,8)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand,fld},Predictions{ff,iter_equalis,p,iter_rand,fld},1);
            for rep=1:1000
                Predictions_tmp=randsample(Predictions{ff,iter_equalis,p,iter_rand,fld},length(Predictions{ff,iter_equalis,p,iter_rand,fld}));
                Performance_rand_within_inter(ff,iter_equalis,p,iter_rand,fld,1:7,rep)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand,fld},Predictions_tmp);
                [~,~,~,Performance_rand_within_inter(ff,iter_equalis,p,iter_rand,fld,8,rep)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand,fld},Predictions_tmp,1);
            end
        end
    end
end
Performance_within_inter=squeeze(nanmean(nanmean(Performance_within_inter),5));
Performance_rand_within_inter=squeeze(nanmean(nanmean(Performance_rand_within_inter),5));
[1]

ictal_or_inter='ictal';
load(['Within_subject_performance_',ictal_or_inter,'_all_feats_crcted_feats_imp_100_10.mat'])
classif=1; %
ff=1;
iter_rand=1;
for fld=1:size(Ground_truth,5)
    for iter_equalis=1:size(Ground_truth,2)
        for p=1:size(Ground_truth,3)
            Performance_within_ictal(ff,iter_equalis,p,iter_rand,fld,1:7)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand,fld},Predictions{ff,iter_equalis,p,iter_rand,fld});
            [~,~,~,Performance_within_ictal(ff,iter_equalis,p,iter_rand,fld,8)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand,fld},Predictions{ff,iter_equalis,p,iter_rand,fld},1);
            for rep=1:1000
                Predictions_tmp=randsample(Predictions{ff,iter_equalis,p,iter_rand,fld},length(Predictions{ff,iter_equalis,p,iter_rand,fld}));
                Performance_rand_within_ictal(ff,iter_equalis,p,iter_rand,fld,1:7,rep)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand,fld},Predictions_tmp);
                [~,~,~,Performance_rand_within_ictal(ff,iter_equalis,p,iter_rand,fld,8,rep)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand,fld},Predictions_tmp,1);
            end
        end
    end
end
Performance_within_ictal=squeeze(nanmean(nanmean(Performance_within_ictal),5));
Performance_rand_within_ictal=squeeze(nanmean(nanmean(Performance_rand_within_ictal),5));
[2]

ictal_or_inter='interictal';
load(['Generalisation_performance_across_subjects_',ictal_or_inter,'_all_feats_comb_crcted_feats_ovrsmp_imp_100_10.mat'])
classif=1; %
ff=1;
iter_rand=1;
for iter_equalis=1:size(Ground_truth,2)
    for p=1:size(Ground_truth,3)
        Performance_acrss_subj_inter(ff,iter_equalis,p,iter_rand,1:7)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif});
        [~,~,~,Performance_acrss_subj_inter(ff,iter_equalis,p,iter_rand,8)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif},1);
        for rep=1:1000
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
classif=1; %
ff=1;
iter_rand=1;
for iter_equalis=1:size(Ground_truth,2)
    for p=1:size(Ground_truth,3)
        Performance_acrss_subj_ictal(ff,iter_equalis,p,iter_rand,1:7)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif});
        [~,~,~,Performance_acrss_subj_ictal(ff,iter_equalis,p,iter_rand,8)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif},1);
        for rep=1:1000
            Predictions_tmp=randsample(Predictions{ff,iter_equalis,p,iter_rand,classif},length(Predictions{ff,iter_equalis,p,iter_rand,classif}));
            Performance_rand_acrss_subj_ictal(ff,iter_equalis,p,iter_rand,1:7,rep)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions_tmp);
            [~,~,~,Performance_rand_acrss_subj_ictal(ff,iter_equalis,p,iter_rand,8,rep)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions_tmp,1);
        end
    end
end
Performance_acrss_subj_ictal=squeeze(Performance_acrss_subj_ictal);
Performance_rand_acrss_subj_ictal=squeeze(Performance_rand_acrss_subj_ictal);
[4]


load(['Generalisation_performance_across_time_all_feats_comb_crcted_feats_ovrsmp_imp_100_10.mat'])
classif=1; %
ff=1;
iter_rand=1;
for iter_equalis=1:size(Ground_truth,2)
    for p=1:size(Ground_truth,3)
        Performance_acrss_time_inter(ff,iter_equalis,p,iter_rand,1:7)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif});
        [~,~,~,Performance_acrss_time_inter(ff,iter_equalis,p,iter_rand,8)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif},1);
        for rep=1:1000
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
classif=1; %
ff=1;
iter_rand=1;
for iter_equalis=1:size(Ground_truth,2)
    for p=1:size(Ground_truth,3)
        Performance_acrss_time_ictal(ff,iter_equalis,p,iter_rand,1:7)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif});
        [~,~,~,Performance_acrss_time_ictal(ff,iter_equalis,p,iter_rand,8)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif},1);
        for rep=1:1000
            Predictions_tmp=randsample(Predictions{ff,iter_equalis,p,iter_rand,classif},length(Predictions{ff,iter_equalis,p,iter_rand,classif}));
            Performance_rand_acrss_time_ictal(ff,iter_equalis,p,iter_rand,1:7,rep)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions_tmp);
            [~,~,~,Performance_rand_acrss_time_ictal(ff,iter_equalis,p,iter_rand,8,rep)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions_tmp,1);
        end
    end
end
Performance_acrss_time_ictal=squeeze(Performance_acrss_time_ictal);
Performance_rand_acrss_time_ictal=squeeze(Performance_rand_acrss_time_ictal);
[6]

save('random_permutations.mat',...
    'Performance_acrss_subj_ictal','Performance_acrss_subj_inter',...
    'Performance_acrss_time_ictal','Performance_acrss_time_inter',...
    'Performance_rand_acrss_subj_ictal','Performance_rand_acrss_subj_inter',...
    'Performance_rand_acrss_time_ictal','Performance_rand_acrss_time_inter',...
    'Performance_rand_within_ictal','Performance_rand_within_inter',...
    'Performance_within_ictal','Performance_within_inter','-v7.3')
