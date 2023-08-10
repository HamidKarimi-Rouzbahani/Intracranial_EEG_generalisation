clc
clear all
baselined=1;
%% Classification: Within ictal or interictal
all_patients=[1:39 41:56];
ff=0;
% for features_used=[1:34]
features_used=[1:34];
    ff=ff+1;
    for iter_equalis=1
        p=0;
        for Patient=all_patients
            p=p+1;
            
            ictal_or_inter='ictal';
            load(['Effect_size_data_',ictal_or_inter,'_crcted_all_feats_100_to_10_points.mat'],'data_target_all','data_non_target_all','positive_negative_modulation','num_samples_keep')
            
            data_targ_res=[];
            data_non_res=[];
            f=0;
            for feats=features_used
                f=f+1;
                if size(data_non_target_all{Patient,feats},1)>=size(data_target_all{Patient,feats},1)
                    if strcmp(ictal_or_inter,'ictal')
                        if baselined==1
                            targ=squeeze(([(data_target_all{Patient,feats}(:,:,num_samples_keep+1:end))-nanmean(data_target_all{Patient,feats}(:,:,1:num_samples_keep),3)])./abs([nanmean(data_target_all{Patient,feats}(:,:,num_samples_keep+1:end),3)+nanmean(data_target_all{Patient,feats}(:,:,1:num_samples_keep),3)]));
                            data_targ=reshape(targ,[size(targ,1)*size(targ,2)],[]);
                            nont=squeeze(([(data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end))-nanmean(data_non_target_all{Patient,feats}(:,:,1:num_samples_keep),3)])./abs([nanmean(data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end),3)+nanmean(data_non_target_all{Patient,feats}(:,:,1:num_samples_keep),3)]));
                            data_non=reshape(nont,[size(nont,1)*size(nont,2)],[]);
                            samp=randsample([1:size(data_targ,1)],size(data_non,1)-size(data_targ,1),true);
                            data_targ=vertcat(data_targ,data_targ(samp,:));
                        else
                            targ=squeeze((data_target_all{Patient,feats}(:,:,num_samples_keep+1:end)));
                            data_targ=reshape(targ,[size(targ,1)*size(targ,2)],[]);
                            nont=squeeze((data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end)));
                            data_non=reshape(nont,[size(nont,1)*size(nont,2)],[]);
                            samp=randsample([1:size(data_targ,1)],size(data_non,1)-size(data_targ,1),true);
                            data_targ=vertcat(data_targ,data_targ(samp,:));
                        end
                    else
                        targ=squeeze(data_target_all{Patient,feats});
                        data_targ=reshape(targ,[size(targ,1)*size(targ,2)],[]);
                        nont=squeeze((data_non_target_all{Patient,feats}));
                        data_non=reshape(nont,[size(nont,1)*size(nont,2)],[]);
                        samp=randsample([1:size(data_targ,1)],size(data_non,1)-size(data_targ,1),true);
                        data_targ=vertcat(data_targ,data_targ(samp,:));
                    end
                else
                    if strcmp(ictal_or_inter,'ictal')
                        if baselined==1
                            nont=squeeze(([(data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end))-nanmean(data_non_target_all{Patient,feats}(:,:,1:num_samples_keep),3)])./abs([nanmean(data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end),3)+nanmean(data_non_target_all{Patient,feats}(:,:,1:num_samples_keep),3)]));
                            data_non=reshape(nont,[size(nont,1)*size(nont,2)],[]);
                            
                            targ=squeeze(([(data_target_all{Patient,feats}(:,:,num_samples_keep+1:end))-nanmean(data_target_all{Patient,feats}(:,:,1:num_samples_keep),3)])./abs([nanmean(data_target_all{Patient,feats}(:,:,num_samples_keep+1:end),3)+nanmean(data_target_all{Patient,feats}(:,:,1:num_samples_keep),3)]));
                            data_targ=reshape(targ,[size(targ,1)*size(targ,2)],[]);
                            
                            samp=randsample([1:size(data_non,1)],size(data_targ,1)-size(data_non,1),true);
                            data_non=vertcat(data_non,data_non(samp,:));
                        else
                            nont=squeeze((data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end)));
                            data_non=reshape(nont,[size(nont,1)*size(nont,2)],[]);
                            
                            targ=squeeze((data_target_all{Patient,feats}(:,:,num_samples_keep+1:end)));
                            data_targ=reshape(targ,[size(targ,1)*size(targ,2)],[]);
                            
                            samp=randsample([1:size(data_non,1)],size(data_targ,1)-size(data_non,1),true);
                            data_non=vertcat(data_non,data_non(samp,:));
                        end
                    else
                        nont=squeeze(data_non_target_all{Patient,feats});
                        data_non=reshape(nont,[size(nont,1)*size(nont,2)],[]);
                        
                        targ=squeeze((data_target_all{Patient,feats}));
                        data_targ=reshape(targ,[size(targ,1)*size(targ,2)],[]);
                        samp=randsample([1:size(data_non,1)],size(data_targ,1)-size(data_non,1),true);
                        data_non=vertcat(data_non,data_non(samp,:));
                    end
                end
                data_targ_res(:,f)=reshape(data_targ,[size(data_targ,1)*size(data_targ,2)],[]);
                data_non_res(:,f)=reshape(data_non,[size(data_non,1)*size(data_non,2)],[]);
            end
            ytrain=[ones(size(data_targ_res,1),1);zeros(size(data_non_res,1),1)];
            Xtrain=[data_targ_res;data_non_res];
            
            
            ictal_or_inter='interictal';
            load(['Effect_size_data_',ictal_or_inter,'_crcted_all_feats_100_to_10_points.mat'],'data_target_all','data_non_target_all','positive_negative_modulation','num_samples_keep')
            data_targ_res=[];
            data_non_res=[];
            f=0;
            for feats=features_used
                f=f+1;
                if size(data_non_target_all{Patient,feats},1)>=size(data_target_all{Patient,feats},1)
                    if strcmp(ictal_or_inter,'ictal')
                        if baselined==1
                            targ=squeeze(([(data_target_all{Patient,feats}(:,:,num_samples_keep+1:end))-nanmean(data_target_all{Patient,feats}(:,:,1:num_samples_keep),3)])./abs([nanmean(data_target_all{Patient,feats}(:,:,num_samples_keep+1:end),3)+nanmean(data_target_all{Patient,feats}(:,:,1:num_samples_keep),3)]));
                            data_targ=reshape(targ,[size(targ,1)*size(targ,2)],[]);
                            nont=squeeze(([(data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end))-nanmean(data_non_target_all{Patient,feats}(:,:,1:num_samples_keep),3)])./abs([nanmean(data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end),3)+nanmean(data_non_target_all{Patient,feats}(:,:,1:num_samples_keep),3)]));
                            data_non=reshape(nont,[size(nont,1)*size(nont,2)],[]);
                            samp=randsample([1:size(data_targ,1)],size(data_non,1)-size(data_targ,1),true);
                            data_targ=vertcat(data_targ,data_targ(samp,:));
                        else
                            targ=squeeze((data_target_all{Patient,feats}(:,:,num_samples_keep+1:end)));
                            data_targ=reshape(targ,[size(targ,1)*size(targ,2)],[]);
                            nont=squeeze((data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end)));
                            data_non=reshape(nont,[size(nont,1)*size(nont,2)],[]);
                            samp=randsample([1:size(data_targ,1)],size(data_non,1)-size(data_targ,1),true);
                            data_targ=vertcat(data_targ,data_targ(samp,:));
                        end
                    else
                        targ=squeeze(data_target_all{Patient,feats});
                        data_targ=reshape(targ,[size(targ,1)*size(targ,2)],[]);
                        nont=squeeze((data_non_target_all{Patient,feats}));
                        data_non=reshape(nont,[size(nont,1)*size(nont,2)],[]);
                        samp=randsample([1:size(data_targ,1)],size(data_non,1)-size(data_targ,1),true);
                        data_targ=vertcat(data_targ,data_targ(samp,:));
                    end
                else
                    if strcmp(ictal_or_inter,'ictal')
                        if baselined==1
                            nont=squeeze(([(data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end))-nanmean(data_non_target_all{Patient,feats}(:,:,1:num_samples_keep),3)])./abs([nanmean(data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end),3)+nanmean(data_non_target_all{Patient,feats}(:,:,1:num_samples_keep),3)]));
                            data_non=reshape(nont,[size(nont,1)*size(nont,2)],[]);
                            
                            targ=squeeze(([(data_target_all{Patient,feats}(:,:,num_samples_keep+1:end))-nanmean(data_target_all{Patient,feats}(:,:,1:num_samples_keep),3)])./abs([nanmean(data_target_all{Patient,feats}(:,:,num_samples_keep+1:end),3)+nanmean(data_target_all{Patient,feats}(:,:,1:num_samples_keep),3)]));
                            data_targ=reshape(targ,[size(targ,1)*size(targ,2)],[]);
                            
                            samp=randsample([1:size(data_non,1)],size(data_targ,1)-size(data_non,1),true);
                            data_non=vertcat(data_non,data_non(samp,:));
                            
                        else
                            nont=squeeze((data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end)));
                            data_non=reshape(nont,[size(nont,1)*size(nont,2)],[]);
                            
                            targ=squeeze((data_target_all{Patient,feats}(:,:,num_samples_keep+1:end)));
                            data_targ=reshape(targ,[size(targ,1)*size(targ,2)],[]);
                            
                            samp=randsample([1:size(data_non,1)],size(data_targ,1)-size(data_non,1),true);
                            data_non=vertcat(data_non,data_non(samp,:));
                        end
                    else
                        nont=squeeze(data_non_target_all{Patient,feats});
                        data_non=reshape(nont,[size(nont,1)*size(nont,2)],[]);
                        
                        targ=squeeze((data_target_all{Patient,feats}));
                        data_targ=reshape(targ,[size(targ,1)*size(targ,2)],[]);
                        samp=randsample([1:size(data_non,1)],size(data_targ,1)-size(data_non,1),true);
                        data_non=vertcat(data_non,data_non(samp,:));
                    end
                end
                data_targ_res(:,f)=reshape(data_targ,[size(data_targ,1)*size(data_targ,2)],[]);
                data_non_res(:,f)=reshape(data_non,[size(data_non,1)*size(data_non,2)],[]);
            end
            ytest=[ones(size(data_targ_res,1),1);zeros(size(data_non_res,1),1)];
            Xtest=([data_targ_res;data_non_res]);
            
            c=0;
            for i=1:size(Xtrain,2)
                if mean(~isnan(Xtrain(:,i)))>0.5 && mean(~isnan(Xtest(:,i)))>0.5
                    c=c+1;
                    %                     Xnewtrain(:,c)=Xtrain(:,i);
                    %                     Xnewtest(:,c)=Xtest(:,i);
                    Xnewtrain(:,c)=normalize(Xtrain(:,i));
                    Xnewtest(:,c)=normalize(Xtest(:,i));
                    included_feats{ff,iter_equalis,p}(c)=i;
                end
            end
            if exist('Xnewtrain','var') && exist('Xnewtest','var')
                Xtrain=(Xnewtrain);
                Xtest=(Xnewtest);
            
            %         Xtrain=normalize(Xnewtrain,1);
            %         Xtest=normalize(Xnewtest,1);
            
            for iter_rand=1:1
                if iter_rand~=1
                    ytrain_r=randsample(ytrain,length(ytrain));
                else
                    ytrain_r=ytrain;
                end
                
                            tic
                            Classifier_Model = TreeBagger(50,Xtrain,ytrain_r,...
                                Method="classification",...
                                OOBPrediction="on",OOBPredictorImportance="on");
                            impCART{ff,iter_equalis,p,iter_rand} = Classifier_Model.OOBPermutedPredictorDeltaError;
                            preds_tmp=predict(Classifier_Model,Xtest);
                            for i=1:length(preds_tmp)
                                preds(i,1)=str2num(preds_tmp{i});
                            end
                            Predictions{ff,iter_equalis,p,iter_rand,1}=preds;
                            clearvars preds
                            Classifiers={'RF'};
                            toc
                
%                 tic
%                 Classifier_Model = fitcdiscr(Xtrain,ytrain_r,'DiscrimType','pseudoLinear');
%                 Predictions{ff,iter_equalis,p,iter_rand,1}=predict(Classifier_Model,Xtest);
%                 toc
%                 Classifiers={'LDA'};
                                
                Ground_truth{ff,iter_equalis,p,iter_rand}=ytest;
            end
            else
                Predictions{ff,iter_equalis,p,1,1}=nan;
                Ground_truth{ff,iter_equalis,p,1}=nan;
                Classifiers={'Any'};
            end
            [ff iter_equalis p]
            clearvars Xtest Xtrain ytrain ytrain_r ytest Xnewtrain Xnewtest
            save(['Generalisation_performance_across_time_ict_to_int_all_feats_comb_crcted_feats_ovrsmp_imp_100_10.mat'],'Classifiers','Predictions','Ground_truth','impCART','included_feats','-v7.3')
        end
    end
ccc
%% Plotting
clc;
clear all;
close all;
load(['Generalisation_performance_across_time_ict_to_int_all_feats_comb_crcted_feats_ovrsmp_imp_100_10.mat'])
Measures={'accuracy','sensitivity','specificity','precision','recall','f1-measure','gmean'};
classif=1; %
Classifiers={'RF'};
ff=1;
iter_rand=1;
for iter_equalis=1:size(Ground_truth,2)
    for p=1:size(Ground_truth,3)
        perf(p,1:7,iter_equalis)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif});
        [~,~,~,auc(p,iter_equalis)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif},1);
    end
end
auc=squeeze(nanmean(auc,2));
perf=squeeze(nanmean(perf,3));
resect=bar([1:8],[nanmean(perf) nanmean(auc)],0.6,'k');
hold on;
errorbar([1:8],[nanmean(perf) nanmean(auc)],[nanstd(perf) nanstd(auc)]/sqrt(size(perf,1)),"LineStyle","none","Color","k")
grid on;
title(['Cross-time generalisation'])
xticks([1:8])
xticklabels([Measures 'auc'])
ylim([0 1])
%% Features importance
feats=[1:8 10:12 14:36];
Features_labels={'Mean','Median','Variance','Skewness','Kurtosis',...
    'LZ Comp','Higichi FD','Katz FD','Lyap Exp','Hurst Exp',...
    'Samp Ent','Apprx Ent','Autocor','Hjorth Comp','Hjorth Mob',...
    'Mean Freq','Med Freq','Avg Freq','SEF','Pow Med Freq',...
    'Phs Med Freq','Power','Energy ratio','Delta','Theta',...
    'Alpha','Beta','Gamma','H-Gamma','Correl','Delta Coh',...
    'Theta Coh','Alpha Coh','Beta Coh','Gamma Coh','HGamma Coh'};
iter_rand=1;
feat_imp=nan(1,21,1);
for iter_equalis=1
    for p=1:55
        feat_imp(p,included_feats{ff,iter_equalis,p,iter_rand},iter_equalis)=impCART{ff,iter_equalis,p,iter_rand};
    end
end
figure;
bar([1:length(feats)],[nanmean(feat_imp,1)],0.6,'k');
hold on;
errorbar([1:length(feats)],[nanmean(feat_imp,1)],[nanstd(feat_imp)]/sqrt(size(feat_imp,1)),"LineStyle","none","Color","k")
grid on;
ylabel('Feature importance')
xticks([1:length(feats)])
xticklabels([Features_labels(feats)])
%% Per feature
clc;
clear all;
close all;
Features_labels={'Mean','Median','Variance','Skewness','Kurtosis',...
    'LZ Comp','Higichi FD','Katz FD','Lyap Exp','Hurst Exp',...
    'Samp Ent','Apprx Ent','Autocor','Hjorth Comp','Hjorth Mob',...
    'Mean Freq','Med Freq','Avg Freq','SEF','Pow Med Freq',...
    'Phs Med Freq','Power','Energy ratio','Delta','Theta',...
    'Alpha','Beta','Gamma','H-Gamma','Correl','Delta Coh',...
    'Theta Coh','Alpha Coh','Beta Coh','Gamma Coh','HGamma Coh'};

Measures={'accuracy','sensitivity','specificity','precision','recall','f1-measure','gmean'};
load(['Generalisation_performance_across_time_each_feat_crcted_feats_ovrsmp_imp_100_10.mat'])
classif=1; %
feats=[1:8 10:12 14:36];
ff=1;
for ff=1:length(feats)
    iter_rand=1;
    for iter_equalis=1:1
        for p=1:55
            if isnan(Predictions{ff,iter_equalis,p,iter_rand,classif})
                perf(ff,p,1:7,iter_equalis)=nan;
                auc(ff,p,iter_equalis)=nan;
            else
                perf(ff,p,1:7,iter_equalis)=Evaluate(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif});
                [~,~,~,auc(ff,p,iter_equalis)]=perfcurve(Ground_truth{ff,iter_equalis,p,iter_rand},Predictions{ff,iter_equalis,p,iter_rand,classif},1);
            end
        end
    end
end

for metric=1:7
    
    figure;
    if metric<7
        perft=squeeze(perf(:,:,metric));
    else
        perft=auc;
    end
    resect=bar([1:length(feats)],[nanmean(perft,2)],0.6,'k');
    hold on;
    errorbar([1:length(feats)],[nanmean(perft,2)],[nanstd(perft')]/sqrt(size(perft,2)),"LineStyle","none","Color","k")
    grid on;
    if metric <7
        title(['Cross-time generalisation; ',Measures(metric)])
    else
        title(['Cross-time generalisation; AUC'])
    end
    xticks([1:length(feats)])
    xticklabels([Features_labels(feats)])
    ylim([0 1])
end