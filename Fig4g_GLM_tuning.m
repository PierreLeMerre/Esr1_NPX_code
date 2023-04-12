%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script that runs a Generalized Linear Model on spike data from a nwb file
% using behavioral signals as predictors.
%
%
% Written by pierre Le Merre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cls
format long
set(0, 'DefaultFigureRenderer', 'painters')

%% Load
% Restrict GLM to a time period
GLM_restrict = 0; % 1 yes, 0 no
last_block = 4; % if yes indicate last block

%Chose genotype
genotype2plot = 'WT';

% Allen Brain Atlas regions to count the units
Regions = {'ACAd','PL','ILA','ORBm'};

% binning window
bin_sz = 0.1; % in seconds

% Min ISI violation rate
MIN_ISI = 0.01; % less than 1% of ISI violations

% Min FR
MIN_FR = 0.1; % mean session FR above 0.1Hz

% Task
Task = 'Aversion';

% Your path to the Database
Path2Data = '/Volumes/labs/dmclab/Pierre/NPX_Database/mPFC/Aversion/';
D = dir([Path2Data '*.nwb']);

% Your path to analysis folder
BasePath = '/Users/pierre/Documents/MATLAB/';
Path2Ana = [BasePath 'neuropixelPFC/Matlab/analysis/'];

% Randomization
nboot = 1000; % nb of surrogate data

% files to plot
%Chose genotype
Genotype = genotype2plot;
if strcmp(Genotype,'NPY')
    f_start = 1;
    f_stop = 5;
elseif strcmp(Genotype,'VGlut2')
    f_start = 6;
    f_stop = 10;
elseif strcmp(Genotype,'WT')
    f_start = 11;
    f_stop = 15;
elseif strcmp(Genotype,'Esr-retro')
    f_start = 22;
    f_stop = 26;
elseif strcmp(Genotype,'Esr-antero')
    f_start = 27;
    f_stop = 31;
elseif strcmp(Genotype,'Esr')
    f_start = 22;
    f_stop = 31;
elseif strcmp(Genotype,'All')
    f_start = 1;
    f_stop = 31;
end



for f = f_start : f_stop
    
    if strcmp(Genotype,'All') && 15<f && f<22
        
    else
        
        % Read nwb file
        nwb = nwbRead([Path2Data D(f).name]);
        % get id and recording duration from meta file
        id = strsplit(nwb.general_session_id,'_');
        id2 = strsplit(id{2},'-');
        disp(['Computing GLM for mouse ' id{1} ', session ' id{2} '...'])
        
        % Get PFC regions
        ede_regions = nwb.general_extracellular_ephys_electrodes.vectordata.get('location').data.load();
        % Get unit main channel
        main_ch = nwb.units.electrodes.data.load();
        
        isregion_mtrx = zeros(numel(main_ch),numel(Regions));
        % loop over regions
        for R = 1:numel(Regions)
            region=char(Regions{R});
            isregion = zeros(1,numel(main_ch));
            
            % loop over units
            for i = 1 : numel(main_ch)
                
                % Make logical region vector
                if (strfind(char(ede_regions{main_ch(i)+1}),region)) > 0 % main_ch is 0 indexed!
                    isregion(i) = 1;
                end
                
            end
            isregion_mtrx(:,R) = isregion;
        end
        
        
        % Get unit spiketimes, Load jagged arrays
        unit_times_data = nwb.units.spike_times.data.load();
        unit_times_idx = nwb.units.spike_times_index.data.load();
        unit_ids = nwb.units.id.data.load(); % array of unit ids
        % Initialize times Map containers indexed by unit_ids
        unit_times = containers.Map('KeyType',class(unit_ids),'ValueType','any');
        last_idx = 0;
        for u = 1:length(unit_ids)
            unit_id = unit_ids(u);
            s_idx = last_idx + 1;
            e_idx = unit_times_idx(u);
            unit_times(unit_id) = unit_times_data(s_idx:e_idx);
            last_idx = e_idx;
        end
        
        % Get quality metrics
        spk_number = nan(1,length(unit_ids));
        spk_number(1) = unit_times_idx(1);
        for u = 2:length(unit_ids)
            spk_number(u) = unit_times_idx(u)-unit_times_idx(u-1);
        end
        qual = nwb.units.vectordata.get('quality').data.load;
        isi = qual(1,:);
        isi_perc = isi./spk_number;
        
        % Get mean session FR
        load([BasePath 'neuropixelPFC/Matlab/analysis/Mean_FR/' D(f).name '_fr.mat'])
        mean_session_fr = firing_rates;
        
        % behavioral signals for the task
        
        load([Path2Ana 'Binned_FR/' nwb.general_session_id  '_binned_fr.mat'])
        load([Path2Ana 'Binned_sounds/' nwb.general_session_id  '_binned_sounds_pure.mat'])
        load([Path2Ana 'Binned_sounds/' nwb.general_session_id  '_binned_sounds_blue.mat'])
        load([Path2Ana 'Binned_airpuff/' nwb.general_session_id  '_binned_airpuff_long.mat'])
        load([Path2Ana 'Binned_opto/' nwb.general_session_id  '_binned_opto.mat'])
        
        
        
%         %% Restrict GLM
    rec_dur = size(spk_cnt,2)*bin_sz;
    t = linspace(1,rec_dur,round(rec_dur/bin_sz));
%         
%         if GLM_restrict
%             
%             blocks = nwb.intervals_trials.vectordata.get('Block').data.load();
%             bb = find(blocks==last_block+1,1,'first');
%             trial_ts = nwb.intervals_trials.start_time.data.load();
%             stop_time = trial_ts(bb)-2; % 2s before next block
%             stop_bin = round(stop_time/bin_sz);
%             
%             % Restrict
%             t = t(:,1:stop_bin);
%             spk_cnt = spk_cnt(:,1:stop_bin);
%             sound_bin = sound_bin(1:stop_bin);
%             %sound_bin2 = sound_bin2(1:stop_bin);
%             opto_bin = opto_bin(1:stop_bin);
%             air_bin = air_bin(1:stop_bin);
%             
%         end
        
        %% Previsualize
        fig1 = figure;
        
        ax1 = subplot(6,1,1:2);
        imagesc(t,'CData',spk_cnt)
        caxis([0, 10])
        title('firing rate')
        colormap(flipud(gray))
        ylabel('units')
        
        ax2 = subplot(6,1,3);
        plot(mean(spk_cnt),'r')
        title('population rate')
        
        ax3 = subplot(6,1,4);
        plot(sound_bin,'k')
        hold on
        plot(sound_bin,'b')
        title('Sounds')
        
        ax4 = subplot(6,1,5);
        plot(air_bin,'k')
        title('Air Puff')
        
        ax5 = subplot(6,1,6);
        plot(opto_bin,'k')
        title('Opto')
        
        linkaxes([ax1,ax2,ax3,ax4,ax5],'x')
        
        sgtitle('Data and predictors')
        
        
        %% Run GLM
                
        %GLM parameters
        opts = statset('glmfit');
        opts.MaxIter = 1000; % default value for glmfit is 100, we increase it here to 1000.
        
        n_predictors = 4;
        X = [sound_bin' opto_bin' sound_bin2' air_bin']; %predictors n by p of p predictors 
        
        
        distr = 'poisson'; % Poisson link function
        t = nan(nboot+1,size(spk_cnt,1),n_predictors+1);
        p = nan(nboot+1,size(spk_cnt,1),n_predictors+1);
        beta = nan(nboot+1,size(spk_cnt,1),n_predictors+1);
        
        if ~(exist([Path2Ana 'GLM2/' nwb.general_session_id  '.mat'],'file')==2)
            
            for ii = 1:nboot+1
                if mod(ii,100)==0
                    disp(['Reached iteration number ' num2str(ii)])
                end
                if ii==1
                    % Run GLM
                    for k=1:size(spk_cnt,1)
                        
                        if sum(isregion_mtrx(k,:)>=1) == 1 && isi_perc(k)<MIN_ISI && mean_session_fr(k)>MIN_FR
                            if mean(spk_cnt(k,:)./0.2)>1 % mean FR > 0.2Hz
                            y = spk_cnt(k,:)'; % n by 1 observed response
                            [B,DEV,STATS] = glmfit(X, y, distr,'options', opts);
                            beta(ii,k,:) = STATS.beta;
                            t(ii,k,:) = STATS.t;
                            p(ii,k,:) = STATS.p;
                            end
                        end
                    end
                    
                    
                else
                    % Run GLM for circular binned spike data. If spk_cnt = [1:nbins], spk_cnt_circshift = [[T:x] [1:T]]
                    T = round(rand*size(spk_cnt,2));
                    spk_cnt_circshift = [spk_cnt(:,T+1:end) spk_cnt(:,1:T)];
                    for k=1:size(spk_cnt,1)
                        if mean(spk_cnt(k,:)./0.2)>1 % mean FR > 0.2Hz
                            y = spk_cnt_circshift(k,:)'; % n by 1 observed response
                            [B,DEV,STATS] = glmfit(X, y, distr,'options', opts);
                            beta(ii,k,:) = STATS.beta;
                            t(ii,k,:) = STATS.t;
                            p(ii,k,:) = STATS.p;

                        end
                    end
                    
                end
            end
            
            % Save GLM output
            save([Path2Ana 'GLM2/' nwb.general_session_id  '.mat'],'beta','t','p')
            disp(['GLM output saved in ' Path2Ana])
            
        else
            load([Path2Ana 'GLM2/' nwb.general_session_id  '.mat'])
            
        end
        
        %% PLot
        
        
        
          figure;
        for ii = 1:10+1
            sound_t = t(ii,:,2);
            opto_t = t(ii,:,3);
            sound2_t = t(ii,:,4);
            air_t = t(ii,:,5);
            
            ax1 = subplot(1,4,1);
            histogram(sound_t,'BinLimits',[-50 50],'BinEdges',-50:1:50)
            hold on
            xlim([-50 50])
            title('Sound 1')
            
            ax2 = subplot(1,4,2);
            histogram(opto_t,'BinLimits',[-20 20],'BinEdges',-50:1:50)
            hold on
            xlim([-50 50])
            title('Opto')
            
            ax3 = subplot(1,4,3);
            histogram(sound2_t,'BinLimits',[-50 50],'BinEdges',-50:1:50)
            hold on
            xlim([-50 50])
            title('Sound 2')
            
            ax4 = subplot(1,4,4);
            histogram(air_t,'BinLimits',[-20 20],'BinEdges',-20:1:20)
            hold on
            xlim([-20 20])
            title('Air Puff')
            
            
        end
        
        suptitle('t scores')
        set(gcf,'units','points','position',[600,1200,1920,450])
        
        
    end
end

