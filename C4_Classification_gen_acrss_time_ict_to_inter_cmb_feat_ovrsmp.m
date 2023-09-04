% This code classifies resected and non-resected contacts, but does it by
% generalising across time (interictal to ictal or vice versa):
% trains the classifiers on data from interictal (or vice versa) data and tests them on
% the data from the inctal window (or vice versa) using Decision Tree
% Classifiers

% INPUTS: features data separated by contacts from C3_Separating_target_non_target_contacts_all_feats
% OUTPUTS: classification data to be permuted by C5_Permuting_all_classifciation_windows
%%
clc
clear all
inter_to_ictal=1;% 1=interictal to ictal; 0=ictal to interictal
%% Classification: Within ictal or interictal
all_patients=[1:39 41:56]; % patient indices
ff=0;
features_used=[1:34]; % features to incorporate
ff=ff+1;
for iter_equalis=1:1 % how many iterations
    p=0;
    for Patient=all_patients % classificaiton/generalisation is performed
        % for each patient separately
        p=p+1;
        
        
        % Preparing the training set
        % from interictal to ictal or vice versa?
        if inter_to_ictal==1
            ictal_or_inter='interictal';
        else
            ictal_or_inter='ictal';
        end
        load(['Effect_size_data_',ictal_or_inter,'_crcted_all_feats_100_to_10_points.mat'],'data_target_all','data_non_target_all','positive_negative_modulation','num_samples_keep')
        
        data_targ_res=[];
        data_non_res=[];
        f=0;
        for feats=features_used
            f=f+1;
                % Upsample the data from the class with lower number
                % of contacts (usually resected) to equalise them with
                % the class with higher number of contacts (usually non-resected)
            if size(data_non_target_all{Patient,feats},1)>=size(data_target_all{Patient,feats},1)
                if strcmp(ictal_or_inter,'ictal')% ictal data
                    % Normalising the feature values in the post onset by mean of values in the pre-onset
                    targ=squeeze(([(data_target_all{Patient,feats}(:,:,num_samples_keep+1:end))-nanmean(data_target_all{Patient,feats}(:,:,1:num_samples_keep),3)])./abs([nanmean(data_target_all{Patient,feats}(:,:,num_samples_keep+1:end),3)+nanmean(data_target_all{Patient,feats}(:,:,1:num_samples_keep),3)]));
                    data_targ=reshape(targ,[size(targ,1)*size(targ,2)],[]);
                    
                    nont=squeeze(([(data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end))-nanmean(data_non_target_all{Patient,feats}(:,:,1:num_samples_keep),3)])./abs([nanmean(data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end),3)+nanmean(data_non_target_all{Patient,feats}(:,:,1:num_samples_keep),3)]));
                    data_non=reshape(nont,[size(nont,1)*size(nont,2)],[]);
                    
                    samp=randsample([1:size(data_targ,1)],size(data_non,1)-size(data_targ,1),true);
                    data_targ=vertcat(data_targ,data_targ(samp,:));
                else % interictal data
                    targ=squeeze(data_target_all{Patient,feats});
                    data_targ=reshape(targ,[size(targ,1)*size(targ,2)],[]);
                    
                    nont=squeeze((data_non_target_all{Patient,feats}));
                    data_non=reshape(nont,[size(nont,1)*size(nont,2)],[]);
                    
                    samp=randsample([1:size(data_targ,1)],size(data_non,1)-size(data_targ,1),true);
                    data_targ=vertcat(data_targ,data_targ(samp,:));
                end
            else
                if strcmp(ictal_or_inter,'ictal') % ictal data
                    % Normalising the feature values in the post onset by mean of values in the pre-onset
                    nont=squeeze(([(data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end))-nanmean(data_non_target_all{Patient,feats}(:,:,1:num_samples_keep),3)])./abs([nanmean(data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end),3)+nanmean(data_non_target_all{Patient,feats}(:,:,1:num_samples_keep),3)]));
                    data_non=reshape(nont,[size(nont,1)*size(nont,2)],[]);
                    
                    targ=squeeze(([(data_target_all{Patient,feats}(:,:,num_samples_keep+1:end))-nanmean(data_target_all{Patient,feats}(:,:,1:num_samples_keep),3)])./abs([nanmean(data_target_all{Patient,feats}(:,:,num_samples_keep+1:end),3)+nanmean(data_target_all{Patient,feats}(:,:,1:num_samples_keep),3)]));
                    data_targ=reshape(targ,[size(targ,1)*size(targ,2)],[]);
                    
                    samp=randsample([1:size(data_non,1)],size(data_targ,1)-size(data_non,1),true);
                    data_non=vertcat(data_non,data_non(samp,:));
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
        
        % Now prepare the testing set
        % Similar to the above
        
        % from interictal to ictal or vice versa?
        if inter_to_ictal==1
            ictal_or_inter='ictal';
        else
            ictal_or_inter='interictal';
        end
        
        load(['Effect_size_data_',ictal_or_inter,'_crcted_all_feats_100_to_10_points.mat'],'data_target_all','data_non_target_all','positive_negative_modulation','num_samples_keep')
        data_targ_res=[];
        data_non_res=[];
        f=0;
        for feats=features_used
            f=f+1;
            if size(data_non_target_all{Patient,feats},1)>=size(data_target_all{Patient,feats},1)
                if strcmp(ictal_or_inter,'ictal')
                    targ=squeeze(([(data_target_all{Patient,feats}(:,:,num_samples_keep+1:end))-nanmean(data_target_all{Patient,feats}(:,:,1:num_samples_keep),3)])./abs([nanmean(data_target_all{Patient,feats}(:,:,num_samples_keep+1:end),3)+nanmean(data_target_all{Patient,feats}(:,:,1:num_samples_keep),3)]));
                    data_targ=reshape(targ,[size(targ,1)*size(targ,2)],[]);
                    
                    nont=squeeze(([(data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end))-nanmean(data_non_target_all{Patient,feats}(:,:,1:num_samples_keep),3)])./abs([nanmean(data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end),3)+nanmean(data_non_target_all{Patient,feats}(:,:,1:num_samples_keep),3)]));
                    data_non=reshape(nont,[size(nont,1)*size(nont,2)],[]);
                    
                    samp=randsample([1:size(data_targ,1)],size(data_non,1)-size(data_targ,1),true);
                    data_targ=vertcat(data_targ,data_targ(samp,:));
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
                    nont=squeeze(([(data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end))-nanmean(data_non_target_all{Patient,feats}(:,:,1:num_samples_keep),3)])./abs([nanmean(data_non_target_all{Patient,feats}(:,:,num_samples_keep+1:end),3)+nanmean(data_non_target_all{Patient,feats}(:,:,1:num_samples_keep),3)]));
                    data_non=reshape(nont,[size(nont,1)*size(nont,2)],[]);
                    
                    targ=squeeze(([(data_target_all{Patient,feats}(:,:,num_samples_keep+1:end))-nanmean(data_target_all{Patient,feats}(:,:,1:num_samples_keep),3)])./abs([nanmean(data_target_all{Patient,feats}(:,:,num_samples_keep+1:end),3)+nanmean(data_target_all{Patient,feats}(:,:,1:num_samples_keep),3)]));
                    data_targ=reshape(targ,[size(targ,1)*size(targ,2)],[]);
                    
                    samp=randsample([1:size(data_non,1)],size(data_targ,1)-size(data_non,1),true);
                    data_non=vertcat(data_non,data_non(samp,:));
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
                Xnewtrain(:,c)=normalize(Xtrain(:,i));
                Xnewtest(:,c)=normalize(Xtest(:,i));
                included_feats{ff,iter_equalis,p}(c)=i;
            end
        end
        if exist('Xnewtrain','var') && exist('Xnewtest','var')
            Xtrain=(Xnewtrain);
            Xtest=(Xnewtest);
            
            % randomising class labels (to generate null distribution for statistical
            % testing)?: No, we do it in another file; so iter_rand=1
            for iter_rand=1:1
                if iter_rand~=1
                    ytrain_r=randsample(ytrain,length(ytrain));
                else
                    ytrain_r=ytrain;
                end
                
                % classify the contacts using the decision tree
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
                Ground_truth{ff,iter_equalis,p,iter_rand}=ytest;
            end
        else
            Predictions{ff,iter_equalis,p,1,1}=nan;
            Ground_truth{ff,iter_equalis,p,1}=nan;
            Classifiers={'Any'};
        end
        [ff iter_equalis p]
        clearvars Xtest Xtrain ytrain ytrain_r ytest Xnewtrain Xnewtest
        % saving the classification results + the ground truths data
        if inter_to_ictal==1 % from interictal to ictal or vice versa?
            save(['Generalisation_performance_across_time_all_feats_comb_crcted_feats_ovrsmp_imp_100_10.mat'],'Classifiers','Predictions','Ground_truth','impCART','included_feats','-v7.3')
        else
            save(['Generalisation_performance_across_time_ict_to_int_all_feats_comb_crcted_feats_ovrsmp_imp_100_10.mat'],'Classifiers','Predictions','Ground_truth','impCART','included_feats','-v7.3')
        end
    end
end
