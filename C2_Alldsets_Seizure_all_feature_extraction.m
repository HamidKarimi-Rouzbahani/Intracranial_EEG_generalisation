% This code extracts features from the interictal and ictal recordings 
% and saves the extracted features in new files
% Some features are extracted using external functions which should be
% available in paths accissible to Matlab

% INPUTS: epoched and cleaned data from C1_ds4100_ictal_interictal_epoching
% OUTPUTS: features data to be used by C3_Separating_target_non_target_contacts_all_feats
%% Feature extraction
clc
clear all;
close all;
% adding eeglab toolbox
addpath(genpath('F:\Toolbox\eeglab2021.1'))
ictal_or_inter='interictal'; % 'ictal' or 'interictal'

for Patient=1:56
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
    
    % Which features to be extracted?
    %     {'Mean','Median','Variance','Skewness','Kurtosis',...
    %     'LZ Comp','Higichi FD','Katz FD','Lyap Exp','Hurst Exp',...
    %     'Samp Ent','Apprx Ent','Autocor','Hjorth Comp','Hjorth Mob',...
    %     'Mean Freq','Med Freq','Avg Freq','SEF','Pow Med Freq',...
    %     'Phs Med Freq','Power','Energy ratio','Delta','Theta',...
    %     'Alpha','Beta','Gamma','H-Gamma','Correl','Delta Coh',...
    %     'Theta Coh','Alpha Coh','Beta Coh','Gamma Coh','HGamma Coh'};
    featuress=[1:8 10:12 14:36];
    
    
    
    % Loading the data file from C0
    if strcmp(ictal_or_inter,'ictal')
        trials=ictal_trials;
    elseif strcmp(ictal_or_inter,'interictal')
        trials=inter_ictal_trials;
    end
    % loading each trial/epoch/recodings
    for trial=trials
        tic
        if strcmp(ictal_or_inter,'ictal')
            load([Patient_initials,'_SEEG_Seizure',num2str(trial),'.mat'])
            Fs=header.Fs; % loading sampling frequeny
            
            % Some preferences: including step (jump in signal) size, the
            % length of pre and post onset and length of feature extrction windows
            % all are set according to Fs so that it is almost the same
            % across patients
            step_size=round(1000*(Fs/500)); % /FS ms; ~2000ms
            pre_siz_analysis=15000*(Fs/500);  %/FS ms; ~30s
            post_siz_analysis=15000*(Fs/500); %/FS ms; ~30s
            window_spans=round([1000]*(Fs/500)); %/FS; ~2000 ms
            windows=[1:length([-pre_siz_analysis:post_siz_analysis])];
            Time_windows=[ceil((max(window_spans)/step_size)/2):floor((length(windows)/step_size-(max(window_spans)/step_size)/2))];
            
            onset=event_info.sample(1);
            
            signal=signal(:,onset-pre_siz_analysis:onset+post_siz_analysis);
            
        elseif strcmp(ictal_or_inter,'interictal')
            load([Patient_initials,'_SEEG_InterIctal',num2str(trial),'.mat'])
            Fs=header.Fs; % loading sampling frequency
            
            % Some preferences: including step (jump in signal) size, the
            % length of pre and post onset and length of feature extrction windows
            % all are set according to Fs so that it is almost the same
            % across patients
            step_size=round(1000*(Fs/500)); % /FS ms; ~2000ms
            pre_siz_analysis=1;  %/FS ms; 0s
            post_siz_analysis=150000*(Fs/500); %/FS ms; ~300s
            window_spans=round([1000]*(Fs/500)); %/FS; ~2000ms
            try
                signal=signal(:,pre_siz_analysis:post_siz_analysis);
            catch
                clc
                [Patient trial]
                display('Signal is shorter than 5 mins!')
            end
            windows=[1:size(signal,2)];
            Time_windows=[ceil((max(window_spans)/step_size)/2):floor((length(windows)/step_size-(max(window_spans)/step_size)/2))];
        end
        % channel indices
        channels=[1:length(channels_id)];
        
        % An empty matrix for loading the features into
        features = nan(max(featuress),length(channels),length(Time_windows), length(window_spans));
        t=0;
        for time_windoww=Time_windows % shifting through signal and extracting features
            t=t+1;
            ws=0;
            for window_span=window_spans % We only checked one length of time windows (2s)
                ws=ws+1;
                wind=windows([time_windoww*step_size+1:time_windoww*step_size+window_span]-round(window_span/2));
                coh_done=0;
                for feature=featuress % features when evaluating time windows
                    ch=0;
                    cor_done=0;
                    for channel = channels
                        ch=ch+1;
                        % extracting one trial (or 2s time window of data)
                        trial_data=signal(channel,wind);
                        
                        % Extracting each feature: the details of each feature
                        % are explained in the manuscript
                        if feature==1
                            %% Time features: Sigal mean
                            features(feature,ch,t,ws) = mean(trial_data);
                            
                        elseif feature==2
                            %% Signal median
                            features(feature,ch,t,ws) = median(trial_data);
                            
                        elseif feature==3
                            %% Signal variance
                            features(feature,ch,t,ws) = var(trial_data);
                            
                        elseif feature==4
                            %% Signal skewness
                            features(feature,ch,t,ws) = skewness(trial_data);
                            
                        elseif feature==5
                            %% Signal Kurtosis
                            features(feature,ch,t,ws) = kurtosis(trial_data);
                            
                        elseif feature==6
                            %% LZ complexity
                            threshold = median(trial_data);
                            trial_data(trial_data>= threshold)=1;
                            trial_data(trial_data< threshold)=0;
                            [features(feature,ch,t,ws),~] = calc_lz_complexity(trial_data, 'exhaustive', 1);
                            
                        elseif feature==7
                            %% Higuchi fractal dimension
                            maxtime = length(trial_data);
                            Kmax = 10;
                            features(feature,ch,t,ws) = Higuchi_FD(trial_data,Kmax);
                            
                        elseif feature==8
                            %% Katz fractal dimensions
                            features(feature,ch,t,ws) = Katz_FD(trial_data);
                            
                        elseif feature==9
                            %% Lyapunov exponent (largest LLE)
                            [features(feature,ch,t,ws),~] = lyaprosen(trial_data,0,0);
                            
                        elseif feature==10
                            %% Hurst Exponent
                            features(feature,ch,t,ws) = estimate_hurst_exponent(trial_data);
                            
                        elseif feature==11
                            %% Sample entropy
                            trial_dataN=abs(trial_data/max(abs(trial_data)))';
                            features(feature,ch,t,ws) = entropy (trial_dataN);
                            
                        elseif feature==12
                            %% Approximate Entropy
                            features(feature,ch,t,ws) = ApEn (2,0.2.* std(trial_data),trial_data,1);
                            
                        elseif feature==13
                            %% Within-trial correlation
                            numLags=ceil(length(trial_data)./2);
                            [acf,lags,~] =autocorr(trial_data,numLags);
                            features(feature,ch,t,ws) = mean(acf(2:end));
                            
                        elseif feature==14
                            %% Hjorth complexity
                            % this finds spread of the spectrum and represents the change in frequency
                            % Hcomplexity
                            st_size=1./Fs;
                            data_prime=(diff(trial_data)./st_size);
                            data_second=(diff(data_prime)./st_size);
                            features(feature,ch,t,ws) = (std(data_second).*std(trial_data))./(std(trial_data)).^2;
                            
                        elseif feature==15
                            %% Hmobility
                            st_size=1./Fs;
                            data_prime=(diff(trial_data)./st_size);
                            data_second=(diff(data_prime)./st_size);
                            features(feature,ch,t,ws) = std(data_prime)./std(trial_data);
                            
                        elseif feature==16
                            %% Mean Freq
                            features(feature,ch,t,ws) = meanfreq(trial_data,Fs);
                            
                        elseif feature==17
                            %% Median Freq
                            features(feature,ch,t,ws) = medfreq(trial_data,Fs);
                            
                        elseif feature==18
                            %% Average Freq
                            zeroscount=0;
                            for i=2:length(trial_data)
                                if (trial_data(i)>0 && trial_data(i-1)<0) || (trial_data(i)<0 && trial_data(i-1)>0)
                                    zeroscount=zeroscount+1;
                                end
                            end
                            features(feature,ch,t,ws) = zeroscount.*(length(trial_data)./Fs);
                            
                        elseif feature==19
                            %% Spectral edge frequency 95%
                            if var(trial_data)==0
                                features(feature,ch,t,ws) =0;
                            else
                                Fourier = fft(trial_data)/length(trial_data);
                                Fouriers = (abs(Fourier));                           % Spectrum
                                Fv = linspace(0, 1, fix(length(trial_data)/2)+1)*Fs/2;        % Frequency Vector
                                Iv = 1:length(Fv);                                      % Index Vector
                                IntSpectrum = cumtrapz(Fv, Fouriers(Iv));               % Numeric Integration
                                try
                                    features(feature,ch,t,ws) = interp1(IntSpectrum, Fv, 0.95*IntSpectrum(end), 'linear');    % Interploate To Find ‘SEF’
                                catch
                                    features(feature,ch,t,ws) = Fs*2;
                                end
                            end
                            
                        elseif feature==20
                            %% Power at Median Freq
                            amp = 2*abs(fft(trial_data))/length(trial_data);
                            phs = angle(fft(trial_data));
                            Fv = linspace(0, 1, fix(length(trial_data)/2)+1)*Fs/2;  % Frequency Vector
                            Iv = 1:length(Fv);                                      % Index Vector
                            fr_des = medfreq(trial_data,Fs);                        % Desired Frequency
                            ampv = amp(Iv);                                         % Trim To Length Of ‘Fv’
                            phsv = phs(Iv);                                         % Trim To Length Of ‘Fv’
                            ap = [ampv(:) phsv(:)];                                 % Amplitude & Phase Matrix
                            ap_des = interp1(Fv(:), ap, fr_des, 'linear');
                            features(feature,ch,t,ws) = ap_des(1);
                            
                        elseif feature==21
                            %% Phase at Median Freq
                            amp = 2*abs(fft(trial_data))/length(trial_data);
                            phs = angle(fft(trial_data));
                            Fv = linspace(0, 1, fix(length(trial_data)/2)+1)*Fs/2;  % Frequency Vector
                            Iv = 1:length(Fv);                                      % Index Vector
                            fr_des = medfreq(trial_data,Fs);                        % Desired Frequency
                            ampv = amp(Iv);                                         % Trim To Length Of ‘Fv’
                            phsv = phs(Iv);                                         % Trim To Length Of ‘Fv’
                            ap = [ampv(:) phsv(:)];                                 % Amplitude & Phase Matrix
                            ap_des = interp1(Fv(:), ap, fr_des, 'linear');
                            features(feature,ch,t,ws) = ap_des(2);
                            
                        elseif feature==22
                            %% Signal power
                            features(feature,ch,t,ws) = bandpower(trial_data);
                            
                        elseif feature==23
                            %% H/L energy ratio
                            HB_top=97;
                            HB_bot=12.4;
                            LB_top=12.4;
                            LB_bot=3.5;
                            data_to_process.times=[0:length(trial_data)-1]*(1/Fs);
                            data_to_process.data=trial_data;
                            data_to_process.srate=Fs;
                            data_to_process.nbchan=1;
                            data_to_process.trials=1;
                            data_to_process.trials=1;
                            data_to_process.event=[];
                            data_to_process.pnts=length(trial_data);
                            tmp_filt=pop_eegfiltnew(data_to_process, HB_bot,HB_top,[],1);
                            HB_energy=mean((tmp_filt.data).^2);
                            tmp_filt=pop_eegfiltnew(data_to_process, LB_bot,LB_top,[],1);
                            LB_energy=mean((tmp_filt.data).^2);
                            features(feature,ch,t,ws) = HB_energy./LB_energy;
                            
                        elseif feature==24
                            %% Delta band power
                            freqrange=[0.5 4];
                            features(feature,ch,t,ws)=bandpower(trial_data,Fs,freqrange);
                            
                        elseif feature==25
                            %% Theta band power
                            freqrange=[4 8];
                            features(feature,ch,t,ws)=bandpower(trial_data,Fs,freqrange);
                            
                        elseif feature==26
                            %% Alpha band power
                            freqrange=[8 13];
                            features(feature,ch,t,ws)=bandpower(trial_data,Fs,freqrange);
                            
                        elseif feature==27
                            %% Beta band power
                            freqrange=[13 30];
                            features(feature,ch,t,ws)=bandpower(trial_data,Fs,freqrange);
                            
                        elseif feature==28
                            %% Gamma band power
                            freqrange=[30 90];
                            features(feature,ch,t,ws)=bandpower(trial_data,Fs,freqrange);
                            
                        elseif feature==29
                            %% High-Gamma band power
                            if Fs==256
                                freqrange=[90 128];
                                features(feature,ch,t,ws)=bandpower(trial_data,Fs,freqrange);
                            else
                                freqrange=[90 248];
                                features(feature,ch,t,ws)=bandpower(trial_data,Fs,freqrange);
                            end
                        elseif feature==30
                            %% Inter-channel zero-lag correlation
                            if cor_done==0
                                ICcor_tmp=nan(size(signal,1),size(signal,1));
                                for ch1=1:size(signal,1)
                                    for ch2=ch1+1:size(signal,1)
                                        trial_data1=signal(ch1,wind);
                                        trial_data2=signal(ch2,wind);
                                        ICcor_tmp(ch1,ch2)=corr(trial_data1',trial_data2');
                                    end
                                end
                                cor_done=1;
                            end
                            features(feature,ch,t,ws)=squeeze(nanmean(ICcor_tmp(ch,:),2));
                            
                        elseif feature>30
                            if coh_done==0
                                if Fs==256
                                    ICcoh_tmp=nan(size(signal,1),size(signal,1),33);
                                else
                                    ICcoh_tmp=nan(size(signal,1),size(signal,1),63);
                                end
                                for ch1=1:size(signal,1)
                                    for ch2=ch1+1:size(signal,1)
                                        trial_data1=signal(ch1,wind);
                                        trial_data2=signal(ch2,wind);
                                        [tmp_IC,a,b]=mscohere(trial_data1,trial_data2,hamming(round(length(trial_data1)/30)),round(length(trial_data1)/60),round(length(trial_data1)/8),Fs);
                                        if Fs==256
                                            ICcoh_tmp(ch1,ch2,:)=tmp_IC(1:33);
                                        else
                                            ICcoh_tmp(ch1,ch2,:)=tmp_IC(1:63);
                                        end
                                    end
                                end
                                coh_done=1;
                            end
                            ICcoh=squeeze(nanmean(ICcoh_tmp(ch,:,:),2));
                            %% Inter-channel coherence
                            if feature==31
                                % Delta band
                                features(feature,ch,t,ws)=squeeze(nanmean(ICcoh([1:2])));
                            elseif feature==32
                                % Theta band
                                features(feature,ch,t,ws)=squeeze(nanmean(ICcoh([2:3])));
                            elseif feature==33
                                % Alpha band
                                features(feature,ch,t,ws)=squeeze(nanmean(ICcoh([3:4])));
                            elseif feature==34
                                % Beta band
                                features(feature,ch,t,ws)=squeeze(nanmean(ICcoh([4:8])));
                            elseif feature==35
                                % Gamma band
                                features(feature,ch,t,ws)=squeeze(nanmean(ICcoh([9:23])));
                            elseif feature==36
                                % High gamma band
                                if Fs==256
                                    features(feature,ch,t,ws)=squeeze(nanmean(ICcoh([24:33])));
                                else
                                    features(feature,ch,t,ws)=squeeze(nanmean(ICcoh([24:63])));
                                end
                            end
                        end
                    end
                end
                [Patient trial t]
            end
        end
        % saving the features data
        if strcmp(ictal_or_inter,'ictal')
            save([Patient_initials,'_All2_Features_Seizure',num2str(trial),'_less_time_res.mat'],'features','channels_id',...
                'pre_siz_analysis','post_siz_analysis','window_spans','step_size','featuress',...
                'channels','windows','Time_windows','channel_info')
        elseif strcmp(ictal_or_inter,'interictal')
            save([Patient_initials,'_All2_Features_Interictal',num2str(trial),'_less_time_res.mat'],'features','channels_id',...
                'pre_siz_analysis','post_siz_analysis','window_spans','step_size','featuress',...
                'channels','windows','Time_windows','channel_info')
        end
        toc
    end
end