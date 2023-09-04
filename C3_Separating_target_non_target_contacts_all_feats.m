% This code separates contacts within and outside of the resection zone
% down-samples them into 28 samples and prepares them 
% for clssification in the next script and saves them

% INPUTS: features data from C2_Alldsets_Seizure_all_feature_extraction
% OUTPUTS: features data separated by contacts to be used by C4_Classification files
%%
clc;
clear all;
close all;
% loading the data containing all features from each patient
for Patient=[1:39 41:56]
    clearvars -except data_non_target_all data_target_all sel_group Patient mean_feature Features_labels chars
    ictal_or_inter='interictal'; % 'ictal' or 'interictal'
        
    if Patient==1
        Patient_initials='060';
        modality='seeg';
        ictal_trials=[1:3];
        inter_ictal_trials=[1:2];
    elseif Patient==2
        Patient_initials='064';
        ictal_trials=[1];
        modality='ecog';
        inter_ictal_trials=[2];
    elseif Patient==3
        Patient_initials='065';
        modality='ecog';
        ictal_trials=[1:3];
        inter_ictal_trials=[1:2];
    elseif Patient==4
        Patient_initials='070';
        ictal_trials=[1:3];
        modality='ecog';
        inter_ictal_trials=[1:2];
    elseif Patient==5
        Patient_initials='074';
        ictal_trials=[1:3];
        modality='ecog';
        inter_ictal_trials=[1:2];
    elseif Patient==6
        Patient_initials='075';
        ictal_trials=[1];
        modality='ecog';
        inter_ictal_trials=[1:2];
    elseif Patient==7
        Patient_initials='080';
        ictal_trials=[1:4];
        modality='ecog';
        inter_ictal_trials=[1:2];
    elseif Patient==8
        Patient_initials='082';
        ictal_trials=[1:5];
        modality='ecog';
        inter_ictal_trials=[1:2];
    elseif Patient==9
        Patient_initials='086';
        ictal_trials=[1:2];
        modality='ecog';
        inter_ictal_trials=[1:2];
    elseif Patient==10
        Patient_initials='087';
        ictal_trials=[1:2];
        modality='ecog';
        inter_ictal_trials=[1:2];
    elseif Patient==11
        Patient_initials='088';
        ictal_trials=[1:3];
        modality='ecog';
        inter_ictal_trials=[1:2];
    elseif Patient==12
        Patient_initials='089';
        ictal_trials=[1:4];
        modality='ecog';
        inter_ictal_trials=[1:2];
    elseif Patient==13
        Patient_initials='094';
        ictal_trials=[1:3];
        modality='ecog';
        inter_ictal_trials=[1:2];
    elseif Patient==14
        Patient_initials='097';
        ictal_trials=[1:5];
        modality='ecog';
        inter_ictal_trials=[1:2];
    elseif Patient==15
        Patient_initials='105';
        ictal_trials=[1:2];
        modality='ecog';
        inter_ictal_trials=[1:2];
    elseif Patient==16
        Patient_initials='106';
        ictal_trials=[1:3];
        modality='ecog';
        inter_ictal_trials=[1:2];
    elseif Patient==17
        Patient_initials='107';
        ictal_trials=[1:5];
        modality='ecog';
        inter_ictal_trials=[1:2];
    elseif Patient==18
        Patient_initials='111';
        ictal_trials=[1:5];
        modality='ecog';
        inter_ictal_trials=[1:2];
    elseif Patient==19
        Patient_initials='112';
        ictal_trials=[1:5];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==20
        Patient_initials='114';
        ictal_trials=[1:4];
        modality='ecog';
        inter_ictal_trials=[1:2];
    elseif Patient==21
        Patient_initials='116';
        ictal_trials=[1:3];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==22
        Patient_initials='117';
        ictal_trials=[1:3];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==23
        Patient_initials='123';
        ictal_trials=[1:4];
        modality='ecog';
        inter_ictal_trials=[1:2];
    elseif Patient==24
        Patient_initials='126';
        ictal_trials=[1:4];
        modality='ecog';
        inter_ictal_trials=[1:2];
    elseif Patient==25
        Patient_initials='130';
        ictal_trials=[1:5];
        modality='seeg';
        inter_ictal_trials=[1:2];
        %   elseif Patient==25
        %         Patient_initials='132';
        %         ictal_trials=[1:5];
        %         modality='seeg';
        %         inter_ictal_trials=[1:2];
    elseif Patient==26
        Patient_initials='133';
        ictal_trials=[1:5];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==27
        Patient_initials='134';
        ictal_trials=[1];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==28
        Patient_initials='135';
        ictal_trials=[1:2];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==29
        Patient_initials='138';
        ictal_trials=[1:5];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==30
        Patient_initials='139';
        ictal_trials=[1:3];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==31
        Patient_initials='140';
        ictal_trials=[1:3];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==32
        Patient_initials='141';
        ictal_trials=[1:5];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==33
        Patient_initials='142';
        ictal_trials=[1:3];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==34
        Patient_initials='144';
        ictal_trials=[1:5];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==35
        Patient_initials='146';
        ictal_trials=[1:3];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==36
        Patient_initials='148';
        ictal_trials=[1:5];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==37
        Patient_initials='150';
        ictal_trials=[1:5];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==38
        Patient_initials='151';
        ictal_trials=[1:5];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==39
        Patient_initials='157';
        ictal_trials=[1:5];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==40
        Patient_initials='158';
        ictal_trials=[4];
        modality='seeg';
        inter_ictal_trials=[];
    elseif Patient==41
        Patient_initials='160';
        ictal_trials=[1:3];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==42
        Patient_initials='162';
        ictal_trials=[1:5];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==43
        Patient_initials='163';
        ictal_trials=[1:3];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==44
        Patient_initials='164';
        ictal_trials=[1:3];
        modality='seeg';
        inter_ictal_trials=[1:2];
        %     elseif Patient==45
        %         Patient_initials='165';
        %         seizures=[];
        %         modality='seeg';
        %         inter_ictal_trials=[1:2];
    elseif Patient==45
        Patient_initials='166';
        ictal_trials=[1:2];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==46
        Patient_initials='171';
        ictal_trials=[1:4];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==47
        Patient_initials='172';
        ictal_trials=[1:5];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==48
        Patient_initials='173';
        ictal_trials=[1:5];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==49
        Patient_initials='177';
        ictal_trials=[1:3];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==50
        Patient_initials='179';
        ictal_trials=[1:2];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==51
        Patient_initials='180';
        ictal_trials=[1:5];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==52
        Patient_initials='181';
        ictal_trials=[1:5];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==53
        Patient_initials='185';
        ictal_trials=[1:5];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==54
        Patient_initials='187';
        ictal_trials=[1:5];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==55
        Patient_initials='188';
        ictal_trials=[1:5];
        modality='seeg';
        inter_ictal_trials=[1:2];
    elseif Patient==56
        Patient_initials='190';
        ictal_trials=[1:3];
        modality='seeg';
        inter_ictal_trials=[1:2];
    end
    
    % laoding the data
    clearvars tmp_data tmp_data_baseline
    if strcmp(ictal_or_inter,'ictal')
        load([Patient_initials,'_All2_Features_Seizure',num2str(ictal_trials(1)),'_less_time_res.mat'])
    elseif strcmp(ictal_or_inter,'interictal')
        load([Patient_initials,'_All2_Features_Interictal',num2str(inter_ictal_trials(1)),'_less_time_res.mat'])
    end
    featuress=[1:8 10:12 14:36];
    chans_feats=repmat(1:length(channels),[length(featuress) 1]);
    if strcmp(ictal_or_inter,'ictal')
        trials=ictal_trials;
    elseif strcmp(ictal_or_inter,'interictal')
        trials=inter_ictal_trials;
    end
    % Categorising the two types of channels
    wind_length=1; % only evalauting 1 window size (2s)
    f=0;
    for feat=featuress % separately for each feature
        f=f+1;
        tr=0;
        for trial=trials % trials refer to recordings
            tr=tr+1;
            if strcmp(ictal_or_inter,'ictal')
                load([Patient_initials,'_All2_Features_Seizure',num2str(trial),chars,'.mat'])
                % How many time windows before and after onset are included
                
%                 if feat<24
%                     pre_wind=1:floor(size(features,3)/2);
%                     post_wind=floor(size(features,3)/2)+1:size(features,3);
%                 else
                    pre_wind=1:15;
                    post_wind=16:29;                    
%                 end
            elseif strcmp(ictal_or_inter,'interictal')
                load([Patient_initials,'_All2_Features_Interictal',num2str(trial),chars,'.mat'])
                % How many time windows ar included
                limit_wind=1:find(features(24,1,:)~=0,1,'last');
            end
            
            patient_data_address=['F:\RESEARCH\Hamid\Multicentre dataset\ds004100\sub-HUP',Patient_initials,'\ses-presurgery\ieeg\'];
            folder_address=[patient_data_address,'sub-HUP',Patient_initials,'_ses-presurgery_task-',ictal_or_inter,'_acq-',modality,'_run-',sprintf('%02d',trial)];
            
            %loading channel info file to determine chantacts within and
            %outside of resection
            channel_info=tdfread([folder_address,'_channels.tsv']);
            
            % indices of target channels
            c=0;
            Target_area='resect'; % or 'soz'
            for ch=1:size(channel_info.status_description,1)
                if contains(channel_info.status_description(ch,:),Target_area) && ~isempty(find(contains(channels_id,strtrim(channel_info.name(ch,:)))))
                    c=c+1;
                    Specific_targ_chans{c}=channels_id{1,find(contains(channels_id,strtrim(channel_info.name(ch,:))))};
                end
            end
            
            if exist('Specific_targ_chans','var')
                chans_tmp=[];
                for ch=1:size(Specific_targ_chans,2)
                    inds=find(strcmp(channels_id,Specific_targ_chans{ch}));
                    chans_tmp=horzcat(chans_tmp,chans_feats(:,inds));
                end
                chans_target=chans_tmp;
                included_chans=channels;
                included_chans(chans_tmp(1,:))=[];
                chans_nontarget=chans_feats(:,included_chans);
                
                tmp_data_Target=squeeze(nanmean(features(feat,chans_target(f,:),:,wind_length),4));
                tmp_data_NonTarget=squeeze(nanmean(features(feat,chans_nontarget(f,:),:,wind_length),4));
                if size(tmp_data_NonTarget,2)==1
                    tmp_data_NonTarget=tmp_data_NonTarget';
                elseif size(tmp_data_Target,2)==1
                    tmp_data_Target=tmp_data_Target';
                end
                % tmp_data_Target and tmp_data_NonTarget contain channels
                % with and outside of target area
                
                
                if strcmp(ictal_or_inter,'ictal')                    
                    num_samples_keep=14; % how many time windows are included after onset
                    
                    % Keep the 14 feature values in the post onset and
                    % the pre-onset 
                    for s=1:size(tmp_data_Target,1)
                        data_pre=tmp_data_Target(s,pre_wind);
                        data_post=tmp_data_Target(s,post_wind);
                        data_target(s,tr,:)=[resample(data_pre,num_samples_keep,length(data_pre)) resample(data_post,num_samples_keep,length(data_post))];
                        positive_negative_modulation(Patient,f,s,tr,1)=((nanmean(data_target(s,tr,1:num_samples_keep),3)-nanmean(data_target(s,tr,num_samples_keep:end),3))>0);
                    end
                    % Keep the 14 feature values in the post onset and
                    % the pre-onset                    
                    for n=1:size(tmp_data_NonTarget,1)
                        data_pre=tmp_data_NonTarget(n,pre_wind);
                        data_post=tmp_data_NonTarget(n,post_wind);
                        data_non_target(n,tr,:)=[resample(data_pre,num_samples_keep,length(data_pre)) resample(data_post,num_samples_keep,length(data_post))];
                        positive_negative_modulation(Patient,f,n,tr,2)=((nanmean(data_non_target(n,tr,1:num_samples_keep),3)-nanmean(data_non_target(n,tr,num_samples_keep:end),3))>0);
                    end
                elseif strcmp(ictal_or_inter,'interictal')
                    num_samples_keep=28; % how many time windows are included interictally
                    for s=1:size(tmp_data_Target,1)
                        data_target(s,tr,:)=resample(tmp_data_Target(s,limit_wind),num_samples_keep,size(tmp_data_Target(s,limit_wind),2));
                        positive_negative_modulation(Patient,f,s,tr,1)=(nanmean(data_target(s,tr,:),3)>0);
                    end
                    for n=1:size(tmp_data_NonTarget,1)
                        data_non_target(n,tr,:)=resample(tmp_data_NonTarget(n,limit_wind),num_samples_keep,size(tmp_data_NonTarget(n,limit_wind),2));
                        positive_negative_modulation(Patient,f,n,tr,2)=(nanmean(data_non_target(n,tr,:),3)>0);
                    end
                end
            end
        end
        data_target_all{Patient,f}=data_target;
        data_non_target_all{Patient,f}=data_non_target;
        % data_target_all and data_non_target_all contain the data within
        % and outside of target area
    end
    [Patient]
end
% Saving the data
save(['Effect_size_data_',ictal_or_inter,'_crcted_all_feats_100_to_10_points.mat'],'data_target_all','data_non_target_all','positive_negative_modulation','num_samples_keep','-v7.3')