%% Feature extraction
clc
clear all;
close all;
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
    elseif Patient==57
        Patient_initials='jh103';
        modality='ecog';
        ictal_trials=[1:3];
        inter_ictal_trials=[1:2];
    elseif Patient==58
        Patient_initials='jh105';
        ictal_trials=[1:5];
        modality='ecog';
        inter_ictal_trials=[1:2];
    elseif Patient==59
        Patient_initials='pt1';
        modality='ecog';
        ictal_trials=[1:4];
        inter_ictal_trials=[1:4];
    elseif Patient==60
        Patient_initials='pt2';
        ictal_trials=[1:3];
        modality='ecog';
        inter_ictal_trials=[1:4];
    elseif Patient==61
        Patient_initials='pt3';
        ictal_trials=[1:2];
        modality='ecog';
        inter_ictal_trials=[1:3];
    elseif Patient==62
        Patient_initials='umf001';
        ictal_trials=[1];
        modality='ecog';
        inter_ictal_trials=[1:2];
    end
            
    chars='_less_time_res'; % '' or '_more_windows'
    
    featuress=[1:8 10:12 14:23];
    % Loading the data file
    if strcmp(ictal_or_inter,'ictal')
        trials=ictal_trials;
    elseif strcmp(ictal_or_inter,'interictal')
        trials=inter_ictal_trials;
    end
    
    for trial=trials
        tic
        if strcmp(ictal_or_inter,'ictal')
            load([Patient_initials,'_SEEG_Seizure',num2str(trial),'.mat'])
            Fs=header.Fs;
            step_size=round(100*(Fs/500)); % /FS ms; 200ms
            pre_siz_analysis=15000*(Fs/500);  %/FS ms; 30s
            post_siz_analysis=15000*(Fs/500); %/FS ms; 30s
%             window_spans=round([10 100 500]*(Fs/500)); %/FS ms; 20,200 and 1000 ms
            window_spans=round([1000]*(Fs/500)); %/FS ms; 20,200 and 1000 ms
            windows=[1:length([-pre_siz_analysis:post_siz_analysis])]; % e.g., 1:120,000 time points the whole signal
            Time_windows=[ceil((max(window_spans)/step_size)/2):floor((length(windows)/step_size-(max(window_spans)/step_size)/2))]; % e.g., 2:300 sampled time points
            
            if Patient<57
                onset=event_info.sample(1);
            elseif Patient>56 && Patient<59
                for i=1:size(event_info.trial_type,1)
                    if contains(event_info.trial_type(i,:),'SZ EVENT')
                        onset=event_info.sample(i,1);
                    end
                end
            elseif Patient>58 && Patient<62
                for i=1:size(event_info.trial_type,1)
                    if contains(event_info.trial_type(i,:),'onset')
                        onset=event_info.sample(i,1);
                    end
                end
            elseif Patient==62
                for i=1:size(event_info.trial_type,1)
                    if contains(event_info.trial_type(i,:),'sz start')
                        onset=event_info.sample(i,1);
                    end
                end
            end
            signal=signal(:,onset-pre_siz_analysis:onset+post_siz_analysis);
            
        elseif strcmp(ictal_or_inter,'interictal')
            load([Patient_initials,'_SEEG_InterIctal',num2str(trial),'.mat'])
            Fs=header.Fs;
            step_size=round(100*(Fs/500)); % /FS ms; 200ms
            pre_siz_analysis=1;  %/FS ms; 0s
            post_siz_analysis=150000*(Fs/500); %/FS ms; 300s
%             window_spans=round([100 500 1000]*(Fs/500)); %/FS ms; 200,1000 and 2000 ms
            window_spans=round([1000]*(Fs/500)); %/FS ms; 200,1000 and 2000 ms
            try
                signal=signal(:,pre_siz_analysis:post_siz_analysis);
            catch
                clc
                [Patient trial]
                display('Signal is shorter than 5 mins!')
            end            
            windows=[1:size(signal,2)]; % e.g., 1:120,000 time points the whole signal
            Time_windows=[ceil((max(window_spans)/step_size)/2):floor((length(windows)/step_size-(max(window_spans)/step_size)/2))]; % e.g., 2:300 sampled time points            
        end
        channels=[1:length(channels_id)];
        
        features = nan(max(featuress),length(channels),length(Time_windows), length(window_spans));

        t=0;
        for time_windoww=Time_windows
            t=t+1;
            ws=0;
            for window_span=window_spans
                ws=ws+1;
                wind=windows([time_windoww*step_size+1:time_windoww*step_size+window_span]-round(window_span/2));
                
                for feature=featuress % features when evaluating time windows
                    ch=0;
                    for channel = channels
                        ch=ch+1;
                        trial_data=signal(channel,wind);
                        
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
                            %         threshold = mean(data);
                            trial_data(trial_data>= threshold)=1;
                            trial_data(trial_data< threshold)=0;
                            [features(feature,ch,t,ws),~] = calc_lz_complexity(trial_data, 'exhaustive', 1);
                            
                        elseif feature==7
                            %% Higuchi fractal dimension
                            maxtime = length(trial_data);
                            %                             Kmax = floor(maxtime./2);
                            Kmax = 10;
                            features(feature,ch,t,ws) = Higuchi_FD(trial_data,Kmax);
                            %                 % second implementation
                            %                 [HFD2(channel),~,~,~]  = HFD_LCALC(data);
                            %                 % thrid implementation
                            %                 HFD3(channel) = hfd(data,Kmax);
                            
                        elseif feature==8
                            %% Katz fractal dimensions
                            features(feature,ch,t,ws) = Katz_FD(trial_data);
                            
                        elseif feature==9
                            %% Lyapunov exponent (largest LLE)
                            [features(feature,ch,t,ws),~] = lyaprosen(trial_data,0,0);
                            
                        elseif feature==10
                            %% Hurst Exponent
                            features(feature,ch,t,ws) = estimate_hurst_exponent(trial_data);
                            %                 HE2(channel) = genhurst(data);
                            %                 HE3(channel) = hurstCC(data);
                            
                        elseif feature==11
                            %% Sample entropy
                            trial_dataN=abs(trial_data/max(abs(trial_data)))';
                            features(feature,ch,t,ws) = entropy (trial_dataN);
                            %                 Ent2(channel) = SampEn (2,0.2.* std(data),data,1);
                            
                            %                         features(feature,channel,t)=sampen(trial_data);
                        elseif feature==12
                            %% Approximate Entropy
                            features(feature,ch,t,ws) = ApEn (2,0.2.* std(trial_data),trial_data,1);
                            %                 Ent4(channel) = approx_entropy(2,0.2.* std(data),data);
                            
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
                            %                             HB_power=bandpower(trial_data,Fs,[HB_bot HB_top]);
                            %                             LB_power=bandpower(trial_data,Fs,[LB_bot LB_top]);
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
                            %                     features(feature,ch,t,ws) = HB_power./(HB_power+LB_power);
                            
                            %                             elseif feature==28
                            %                                 %% Wavelet transform
                            %                                 [c,l] = wavedec(trial_data,5,'sym2');
                            %                                 [ca5] = appcoef(c,l,'sym2',5);
                            %                                 [cd1,cd2,cd3,cd4,cd5] = detcoef(c,l,[1 2 3 4 5]);
                            %                                 feature_chosen(channel,1:length([ca5,cd1,cd2,cd3,cd4,cd5]))=[ca5,cd1,cd2,cd3,cd4,cd5];
                            %
                            %                             elseif feature==29
                            %                                 %% Hilbert transform amplitude
                            %                                 Hilb = hilbert(trial_data);
                            %                                 feature_chosen(channel,:)=abs(Hilb);
                            %
                            %                             elseif feature==30
                            %                                 %% Hilbert transform phase
                            %                                 Hilb = hilbert(trial_data);
                            %                                 feature_chosen(channel,:)=angle(Hilb);
                            %
                            %                             elseif feature==31
                            %                                 %% Phase-Amplitude Coupling
                            %                                 data.trial{1,1}= trial_data;
                            %                                 data.time{1,1}= [1:length(trial_data)]./Fs;
                            %                                 data.trialinfo=[100]';
                            %                                 data.label{1,1}='SampleData';
                            %                                 toi=[0.001 1.0]; % time of interest
                            %                                 phase=[0.5 12];   % phase(1):2.5:phase(2)
                            %                                 ampl=[24 120];   % amp(1):19:amp(2)
                            %                                 diag = 'no'; %'yes' or 'no' to turn on or off diagrams during computation
                            %                                 surrogates = 'no'; %'yes' or 'no' to turn on or off surrogates during computation
                            %                                 approach = 'tort';%,'ozkort','canolty','PLV';
                            %                                 [MI_matrix_raw,~] = calc_MI(data,toi,phase,ampl,diag,surrogates,approach);
                            %                                 feature_chosen(channel,:)=reshape(MI_matrix_raw(1:6,1:5),[30 1]);
                            %
                            %                             elseif feature==32
                            %                                 %% Inter-channel correlation 31 feature per trial
                            %                                 ICC=zeros(size(signal,1),1);
                            %                                 for ch2=1:size(signal,1)
                            %                                     trial_data2=signal(ch2,wind);
                            %                                     if band>1
                            %                                         trial_data2=eegfilt(trial_data2,Fs,lowband,highband,0,floor(length(trial_data2)/3),0,'fir1');
                            %                                     end
                            %                                         ICC(ch2,1)=corr(trial_data',trial_data2');
                            %                                 end
                            %                                 feature_chosen(channel,:)=ICC;
                            %
                            %                             elseif feature==33
                            %                                 %% CNN
                            %                                 trial_data=(trial_data+abs(min(trial_data)));
                            %                                 trial_data=(trial_data).*(255./max(trial_data));
                            %
                            %                                 im_ = single(uint8(trial_data));
                            %                                 im_ = repmat(im_,[1 ceil(227*227./length(im_))]);
                            %                                 im_=im_(1,1:227*227);
                            %                                 im_=reshape(im_,[227 227]);
                            %                                 im_ = imresize(im_,net.meta.normalization.imageSize(1:2));
                            %                                 im_ = im_ - net.meta.normalization.averageImage;
                            %
                            %                                 res = vl_simplenn(net, im_);
                            %                                 feature_chosen(channel,:)=squeeze(res.x);
                            %
                            %                             elseif feature==34
                            %                                 %% Signal samples
                            %                                 feature_chosen(channel,:)=trial_data;
                            %
                            %                             elseif feature==35
                            %                                 %% Within-trial correlation
                            %                                 numLags=size(trial_data,2)-1;
                            %                                 [acf,lags,~] =autocorr(trial_data,numLags);
                            %                                 feature_chosen(channel,:)= acf(2:end);
                            %                             end
                        end
                    end
                end
                [Patient trial t]
            end
        end
%         if strcmp(ictal_or_inter,'ictal')
%             save([Patient_initials,'_ds4100_Features_Seizure',num2str(trial),chars,'.mat'],'features','channels_id',...
%                 'pre_siz_analysis','post_siz_analysis','window_spans','step_size','featuress',...
%                 'channels','windows','Time_windows','channel_info')
%         elseif strcmp(ictal_or_inter,'interictal')
%             save([Patient_initials,'_ds4100_Features_Interictal',num2str(trial),chars,'.mat'],'features','channels_id',...
%                 'pre_siz_analysis','post_siz_analysis','window_spans','step_size','featuress',...
%                 'channels','windows','Time_windows','channel_info')
%         end

        if strcmp(ictal_or_inter,'ictal')
            save([Patient_initials,'_All_Features_Seizure',num2str(trial),chars,'.mat'],'features','channels_id',...
                'pre_siz_analysis','post_siz_analysis','window_spans','step_size','featuress',...
                'channels','windows','Time_windows','channel_info')
        elseif strcmp(ictal_or_inter,'interictal')
            save([Patient_initials,'_All_Features_Interictal',num2str(trial),chars,'.mat'],'features','channels_id',...
                'pre_siz_analysis','post_siz_analysis','window_spans','step_size','featuress',...
                'channels','windows','Time_windows','channel_info')
        end
        toc
    end
end