%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script that plot the PSTH of the first blocks per PFC region.
%
% Written by pierre Le Merre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cls
format long
set(0, 'DefaultFigureRenderer', 'painters')

%% Load
% Task

CNT = 1;
P = [];
G = [];
Brain_region = 'mPFC'; %mPFC Aud

% Genotype
genotype2plot = 'All'; %'All', 'WT', 'VGlut2', 'NPY', 'Esr'

% Restrict unit type
Unit2plot = 'All'; % 'FS','RS' or 'All' 

% Allen Brain Atlas regions to count the units
if strcmp(Brain_region,'mPFC')
    Tasks = {'Passive','Context','Aversion','Opto'}; %,'Detection'
    %Regions = {'ACAd','PL','ILA','ORBm'};
    %Regions = {'MOs','ACAd','ACAv','PL','ILA','ORBm','ORBvl','ORBl'};
    Regions = {'MOs','ACAd','ACAv','PL','ILA','ORBm','ORBvl','ORBl','STR','CP','ACB','LS','SH','AON','OLF','DP','TT'};
    HexColors = {'#99B898','#FECEA8','#FECEA8','#FF847C','#6C5B7B','#355C7D','#355C7D','#355C7D','#111111','#111111','#111111','#111111','#111111','#111111','#111111','#111111','#111111'};
    %HexColors = {'#FECEA8','#FF847C','#6C5B7B','#355C7D'};
elseif strcmp(Brain_region,'Aud')
    Tasks = {'Passive','Context','Detection'}; %
    %Regions = {'AUDp','AUDd','AUDpo','AUDv'}; %
    Regions = {'AUDd','AUDp','AUDpo','AUDv','CA1','CA2','CA3','DG','ECT','EPd','LP','MRN','PERI','PIR','POL','ProS','SSp','SSs','TEa','VISal','VISli'};
    HexColors = {'#1261A0','#1261A0','#1261A0','#1261A0','#111111','#111111','#111111','#111111','#111111','#111111','#111111','#111111','#111111','#111111','#111111','#111111','#111111','#111111','#111111','#111111','#111111'};
   
    %Regions = {'AUDd1','AUDd2/3','AUDd4','AUDd5','AUDd6a','AUDd6b'}; %
end

%% Compute PSTH

for T = 1 : numel(Tasks)
    
    %T = 3;
    Task = Tasks{T}; %'Passive','Context','Aversion'
        
    % Your path to the Database
    %Path2Data = '/Volumes/dmclab/Pierre/NPX_Database/';
    %Path2Data = '/Volumes/Maxtor/NWB/';
    Path2Data = ['/Volumes/T7/LeMerre_dataset/' Brain_region '/'];
    Path = [Path2Data Task filesep];
    D = dir([Path '*.nwb']);
    
    % Your path to analysis folder
    Path2Ana = '/Users/pierre/Documents/MATLAB/neuropixelPFC/Matlab/analysis/';
    
    % Load colormaps
    load('/Users/pierre/Documents/MATLAB/neuropixelPFC/Matlab/Utilities/Colormaps/sd_colormap.mat')
    
    % Min ISI violation rate
    MIN_ISI = 0.01; % less than 1% of ISI violations
    
    % Min FR
    MIN_FR = 0.1; % mean session FR above 0.1Hz
        
    % number of block to plot
    Block_nb = 1;
        
    % Save figure
    Fig_save = 1;
    FigPath = '/Users/pierre/OneDrive - KI.SE/Pierre_Shared/LHA-LHb-PFC/NatureSubmission';
    
    % Trial per condition
    Trial_per_condition = 50;
    
    % Task modulated
    task_mod = 0; % Restrict plot to task modulated units only
    
    % files to plot
    if strcmp(Brain_region,'mPFC')
        if strcmp(Task,'Aversion')
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
                f_start = 16;
                f_stop = 20;
            elseif strcmp(Genotype,'Esr-antero')
                f_start = 21;
                f_stop = 25;
            elseif strcmp(Genotype,'Esr')
                f_start = 16;
                f_stop = 25;
            elseif strcmp(Genotype,'All')
                f_start = 1;
                f_stop = 25;
            end
        elseif strcmp(Task,'Detection')
            f_start = 1;
            f_stop = numel(D);
        else
            f_start = 1;
            f_stop = numel(D);
        end
    else

            f_start = 1;
            f_stop = numel(D);

    end
    
    % PSTH heatmap
    %initialize global vars
    Ede_regions = [];
    Main_ch = [];
    ISregion_mtrx = [];
    
   
    for f = f_start : f_stop % Loop through mice
    
    % Read nwb file
    nwb = nwbRead([Path D(f).name]);
    % get id and recording duration from meta file
    id = strsplit(nwb.general_session_id,'_');
    
    disp(['Get units from mouse ' id{1} ', session ' id{2}])
    
    
    
    % Get PFC regions
    ede_regions = nwb.general_extracellular_ephys_electrodes.vectordata.get('location').data.load;
    %Ede_regions = [Ede_regions; ede_regions];
    % Get unit main channel
    main_ch = nwb.units.electrodes.data.load;
    %Main_ch = [Main_ch; main_ch];
    
    %FS and RS
    if strcmp(Unit2plot,'FS')
    fs_rs_bolean = nwb.units.vectordata.get('FS').data.load;
    elseif strcmp(Unit2plot,'RS')
    fs_rs_bolean = nwb.units.vectordata.get('RS').data.load;
    else
    fs_rs_bolean = ones(1,numel(nwb.units.vectordata.get('RS').data.load));    
    end
    
    isregion_mtrx = zeros(numel(main_ch),numel(Regions));
    % loop over regions
    for R = 1:numel(Regions)
        region=char(Regions{R});
        isregion = zeros(1,numel(main_ch));
        
        % loop over units
        for i = 1 : numel(main_ch)            
            % Make logical region vector
            if (strfind(ede_regions{main_ch(i)+1,:},region)) > 0 % main_ch is 0 indexed!
                isregion(i) = 1;
            end            
        end
        isregion_mtrx(:,R) = isregion;
    end
    %ISregion_mtrx = [ISregion_mtrx; isregion_mtrx];
    
    
    % PSTH
    % activity bin parameters
    bin_sz = 0.01; % in seconds
    sr = 30000; % sampling rate
    start_time = -1; % in seconds
    stop_time = 3;
    time_window = stop_time-start_time; % in seconds
    
    % Gaussian Kernel 50ms
    N = 0.05/bin_sz;
    g = gausswin(N);
    w=g./trapz(g); % normalize area = 1
    
    % Get unit spiketimes, Load jagged arrays
    unit_times_data = nwb.units.spike_times.data.load;
    unit_times_idx = nwb.units.spike_times_index.data.load;
    unit_ids = nwb.units.id.data.load; % array of unit ids
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
    
    % Get PSTH Events
    if strcmp(Task,'Aversion')
        trial_ts = nwb.intervals_trials.start_time.data.load;
        if f==3
            blocks = nwb.intervals_trials.vectordata.get('Block').data.load;
            blocks = blocks(1:573);
        elseif f==4
            blocks = nwb.intervals_trials.vectordata.get('Block').data.load;
            blocks = blocks(1:600);
        elseif f==20
            blocks = nwb.intervals_trials.vectordata.get('Block').data.load;
            blocks = blocks(51:end);
        else
            blocks = nwb.intervals_trials.vectordata.get('Block').data.load;
        end
        
    elseif strcmp(Task,'Context')
        trial_ts = nwb.intervals_trials.start_time.data.load;
        blocks = nwb.intervals_trials.vectordata.get('Block').data.load;
        if Block_nb==1
            blocks = blocks(1:Trial_per_condition);
        end
        
    elseif strcmp(Task,'Passive')
        trial_ts = nwb.stimulus_presentation.get('passive_10kHz').timestamps.load;
        blocks = nwb.intervals_trials.vectordata.get('Block').data.load;
        amp = nwb.stimulus_presentation.get('passive_10kHz').data.load;
        trial_ts = trial_ts(amp==1 | amp==2 | amp==3);
        blocks = blocks(1:numel(trial_ts));
        
    elseif strcmp(Task,'Opto')
        trial_ts = nwb.stimulus_presentation.get('AudStim_10kHz').timestamps.load;
        opto_ts = nwb.stimulus_presentation.get('Optogenetics').timestamps.load;
        trial_ts = [trial_ts; opto_ts];
        blocks = nwb.intervals_trials.vectordata.get('Block').data.load;
        blocks = blocks(1:numel(trial_ts));
        amp = nwb.stimulus_presentation.get('AudStim_10kHz').data.load;
        %trial_ts = trial_ts(amp==1 | amp==2 | amp==3);
        %blocks = blocks(1:numel(trial_ts));
        trial_ts = trial_ts(1:numel(blocks));
        
    elseif strcmp(Task,'Detection')
        trial_ts = nwb.stimulus_presentation.get('AudStim_10kHz').timestamps.load;
        blocks = nwb.intervals_trials.vectordata.get('Block').data.load;
        amp = nwb.stimulus_presentation.get('AudStim_10kHz').data.load;
        if strcmp(Brain_region,'Aud') % block number is weird
        blocks = ones(1,numel(trial_ts));
        end
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
    load(['/Users/pierre/Documents/MATLAB/neuropixelPFC/Matlab/analysis/Mean_FR/' D(f).name '_fr.mat'])
    mean_session_fr = firing_rates;
    
    % PSTH Loop
    PSTH_B = nan(Block_nb,numel(Regions),numel(main_ch),time_window/bin_sz);
    GK_B = nan(Block_nb,numel(Regions),numel(main_ch),time_window/bin_sz);
    
    
    % loop over blocks
    for b = 1:Block_nb
        
        if strcmp(id{1},'PL066') && (b==4 || b==5)
        eval(['bloc' num2str(b) '_ts = [];'])      
        elseif strcmp(id{1},'PL066') && (b==6 || b==7 || b==8)
        eval(['bloc' num2str(b) '_ts = trial_ts(blocks==' num2str(b+1) ');'])    
        else    
        eval(['bloc' num2str(b) '_ts = trial_ts(blocks==' num2str(b) ');'])
        end
        
        if b==2 && strcmp(Task,'Aversion')
            bloc2_ts = bloc2_ts(1:2:end); % one trial over 2
        elseif b==1 && strcmp(Task,'Opto')
            bloc1_ts = bloc1_ts(1:5:end);  % one trial over 5     
        end
        eval(['bb = bloc' num2str(b) '_ts;'])
        
        if numel(bb)>Trial_per_condition % make sure to take 50 trials
            eval(['bloc' num2str(b) '_ts = bloc' num2str(b) '_ts(1:' num2str(Trial_per_condition) ');'])
        end
        eval(['bb = bloc' num2str(b) '_ts;'])
        
        % get saturation epochs
        sat_start = nwb.analysis.get('spike saturation start').timestamps.load;
        sat_stop = nwb.analysis.get('spike saturation stop').timestamps.load;
        
        
        PSTH = nan(numel(Regions),numel(main_ch),time_window/bin_sz);
        GK = nan(numel(Regions),numel(main_ch),time_window/bin_sz);
        
        % Loop through regions
        for R = 1:numel(Regions)
            
            psth = nan(numel(main_ch),time_window/bin_sz);
            gk = nan(numel(main_ch),time_window/bin_sz);
            fr_base = nan(1,numel(bb));
            fr_cs = nan(1,numel(bb));
            fr_us = nan(1,numel(bb));
            cnt = 1;
            
            % Loop through units
            if task_mod
                load([Path2Ana 'Task_modulated_units/' id{1} '_' id{2} '_tmu.mat'])
            end
            
            for i = 1 : numel(unit_ids)
                
                if isregion_mtrx(i,R) == 1 && isi_perc(i)<MIN_ISI && mean_session_fr(i)>MIN_FR && fs_rs_bolean(i)==1
                    
                    if task_mod==1
                        if task_modulated(i) == 1
                            
                            unit = i;
                            vecSpikeTimes = unit_times(unit-1); % -1 indexed in unit_times!!!
                            
                            % exclude spikes in saturation epochs
                            for j = 1:numel(sat_start)
                                s1 = sat_start(numel(sat_start)-j+1);
                                s2 = sat_stop(numel(sat_start)-j+1);
                                idx = find(0<vecSpikeTimes-s1 &  vecSpikeTimes-s2<0);
                                if ~isempty(idx)
                                    vecSpikeTimes(idx) = [];
                                end
                            end
                            
                            %exclude trial_ts in saturation epochs
                            eval(['vecEventTimes = bloc' num2str(b) '_ts;'])
                            idx2 = find(0<vecEventTimes-s1 &  vecEventTimes-s2<0);
                            if ~isempty(idx2)
                                vecEventTimes(idx2) = [];
                            end
                            
                            % calculate simple PSTH
                            [SpikeRates,WindowCenters] = PSTH_Simple(vecSpikeTimes.*sr, vecEventTimes.*sr, start_time, stop_time, sr, bin_sz);
                            psth(cnt,:) = mean(SpikeRates,2)';
                            gk(cnt,:)=conv(squeeze(psth(cnt,:)),w,'same'); % Gaussian smoothing
                            %increase counter
                            cnt = cnt +1;
                            
                        end
                        
                        
                    else
                        
                        unit = i;
                        vecSpikeTimes = unit_times(unit-1); % -1 indexed in unit_times!!!
                        
                        
                        
                        % exclude spikes in saturation epochs
                        for j = 1:numel(sat_start)
                            s1 = sat_start(numel(sat_start)-j+1);
                            s2 = sat_stop(numel(sat_start)-j+1);
                            idx = find(0<vecSpikeTimes-s1 &  vecSpikeTimes-s2<0);
                            if ~isempty(idx)
                                vecSpikeTimes(idx) = [];
                            end
                        end
                        
                        %exclude trial_ts in saturation epochs
                        eval(['vecEventTimes = bloc' num2str(b) '_ts;'])
                        idx2 = find(0<vecEventTimes-s1 &  vecEventTimes-s2<0);
                        if ~isempty(idx2)
                            vecEventTimes(idx2) = [];
                        end
                        
                        % calculate simple PSTH
                        [SpikeRates,WindowCenters] = PSTH_Simple(vecSpikeTimes.*sr, vecEventTimes.*sr, start_time, stop_time, sr, bin_sz);
                        psth(cnt,:) = mean(SpikeRates,2)';
                        gk(cnt,:)=conv(squeeze(psth(cnt,:)),w,'same'); % Gaussian smoothing
                        %increase counter
                        cnt = cnt +1;
                        
                    end
                    
                    
                    
                end
            end % end of unit loop
            PSTH(R,:,:) = psth;
            GK(R,:,:) = gk;
            
        end % end of region loop
        
        PSTH_B(b,:,:,:) = PSTH;
        GK_B(b,:,:,:) = GK;
        
        
    end % end of block loop
    
    eval(['P.mouse' num2str(CNT) '= PSTH_B;'])
    eval(['G.mouse' num2str(CNT) '= GK_B;'])
    CNT = CNT + 1;
    end % end of mouse loop
    
end %end of Task loop

%% Plotting
%Heatmaps
time = linspace(start_time, stop_time, time_window/bin_sz);
h_delim = nan(1,numel(Regions));
BM = nan(numel(Regions),1);
bm = [];
BS = nan(numel(Regions),1);
P1M = nan(numel(Regions),1);
p1m = [];
P1S = nan(numel(Regions),1);
P2M = nan(numel(Regions),1);
p2m = [];
P2S = nan(numel(Regions),1);

fig1=figure;
for b = 1:Block_nb % loop over blocks
    subplot(1,Block_nb,b);
    temp3 = [];
    Z_temp = [];
    for R = 1:numel(Regions) % loop over regions
        temp2 = [];
        for f = 1 : CNT-1 % loop over mice
            
            eval(['pp = G.mouse' num2str(f) ';'])
            temp = squeeze(pp(b,R,:,:));
            temp(isnan(temp(:,1)),:) = [];
            if size(temp2,2)==1
                temp2 = temp2';
            end
            temp2 = [temp2; temp];
        end
        % reorganize rows by peak activity 500ms after sound
        if b==1
            % reorder by peak
            loc1 = nan(1,size(temp2,1));
            for ii = 1:size(temp2,1)
                [~,loc1(ii)] = max(temp2(ii,101:150));
            end
            [~,order1] = sort(loc1);
            eval(['fr_order' num2str(R) '=order1;'])
        else
        end
        eval(['order = fr_order' num2str(R) ';'])
        if ~isempty(temp2) && b~=4 && b~=5
        temp2 = temp2(order,:);
        end
        temp3 = [temp3; temp2];
        h_delim(R) = size(temp3,1);
    end
    Z_temp = (temp3 - mean(temp3(:,1:100),2))./std(temp3(:,1:100),[],2); % manual baseline Z-score
    for i = 1:size(Z_temp,1) % loop to remove nan values (no spikes during the baseline)
        if isnan(Z_temp(i,1))
            Z_temp(i,:) = zeros(1,time_window/bin_sz);
        end
    end
    
    % plot
    imagesc('XData',time,'CData',flip(Z_temp));
    xlim([start_time stop_time]);
    if ~isempty(temp3)
    ylim([0 size(temp3,1)]);
    end
    title(['Block ' num2str(b)]);
    
    if strcmp(Task,'Opto') || strcmp(Task,'Context') || strcmp(Task,'Passive')
        line([0 0],[0 size(temp3,1)],'LineWidth',1,'Color',[0 0 0])
    elseif strcmp(Task,'Aversion')
        if b==1 || b==2 || b==3 || b==4
            line([0 0],[0 size(temp3,1)],'LineWidth',1,'Color',[0 0 0])
        end
        if b==2 || b==4 || b==5 || b==6 || b==7 || b==8
            line([0.7 0.7],[0 size(temp3,1)],'LineWidth',1,'Color',[0 0 0])
        end
        
    end
    for i=1:numel(h_delim)-1
        line([start_time stop_time],[size(Z_temp,1)-h_delim(i) size(Z_temp,1)-h_delim(i)],'LineWidth',1,'Color',[0 0 0])
    end
    
    colormap(sd_colormap);
    c=colorbar('southoutside');
    c.Label.String = 'Firing rate (sd)';
    caxis([-3.5, 3.5])
    set(gcf,'units','points','position',[0,900,Block_nb*200,900])
end

% Count units per region
units_per_region = zeros(1,numel(Regions));
for i = 1:numel(Regions)
    if i ==1
        units_per_region(i) =  h_delim(i);
    else
        units_per_region(i) =  h_delim(i)-h_delim(i-1);
    end
end

if strcmp(Task,'Aversion')
    sgtitle(['PSTH, Aversion ' Genotype ' (n = ' num2str(sum(units_per_region)) ' ' Unit2plot ' units)'])
else
    sgtitle(['PSTH, ' Task ' (n = ' num2str(sum(units_per_region)) ' ' Unit2plot ' units)'])
end

total_per_region = sum(ISregion_mtrx,1);
if task_mod==0
    if Fig_save==1
        saveas(fig1,[FigPath '/' Genotype '_' Unit2plot '_units_HM' ],'epsc');
    end
elseif task_mod==1
    if Fig_save==1
        saveas(fig1,[FigPath '/' Genotype '_' Unit2plot '_units_HM_modulated_units'],'epsc');
    end
end


% psth
fig2=figure;
for b = 1:Block_nb % loop over blocks
    subplot(1,Block_nb,b);
    temp3 = [];
    Z_temp = [];
    offset = numel(Regions)*10;
    
    for R = 1:numel(Regions) % loop over regions
        temp2 = [];
        
        for f = 1 : CNT-1 % loop over mice
            
            eval(['pp = G.mouse' num2str(f) ';'])
            temp = squeeze(pp(b,R,:,:));
            temp(isnan(temp(:,1)),:) = [];
             if size(temp2,2)==1
                temp2 = temp2';
            end
            temp2 = [temp2; temp];
        end
        if ~isempty(temp2)
            disp(['Adding ' num2str(numel(mean(temp2,2))) ' units'])
            mean2plot = nanmean(temp2-mean(temp2(:,1:100),2))+offset; % baseline correction before mouse average i.e on the units       
            sem2plot = nanstd(temp2)/sqrt(size(temp2,2));
            %Mean Baseline
            %BM(R) = nanmean(temp2,2)-nanmean(temp2(:,1:100),2);
            BM(R) = nanmean(nanmean(temp2(:,1:100),2));
            bm = [bm; mean(temp2(:,1:100),2)];
            %bm = [bm; mean(temp2(:,1:100),2)];
            BS(R) = nanstd(nanmean(temp2(:,1:100),2));
            %Mean Peak 1
            %P1M(R) = nanmean(temp2,2)-nanmean(temp2(:,101:150),2);
            P1M(R) = nanmean(nanmean(temp2(:,101:150),2));
            p1m = [p1m; mean(temp2(:,101:150),2)-mean(temp2(:,1:100),2)];
            %p1m = [p1m; mean(temp2(:,101:150),2)];
            P1S(R) = nanstd(nanmean(temp2(:,101:150),2));
            %Mean Peak 2
            %P2M(R) = nanmean(temp2,2)-nanmean(temp2(:,171:230),2);
            P2M(R) = nanmean(nanmean(temp2(:,171:230),2));
            p2m = [p2m; mean(temp2(:,171:230),2)-mean(temp2(:,1:100),2)];
            %p2m = [p2m; mean(temp2(:,171:230),2)];
            P2S(R) = nanstd(nanmean(temp2(:,171:230),2));
            % save single unit matrix
            Mean_Mtrx = [bm p1m p2m];
            save([FigPath '/stats/BLock' num2str(b) '_BcorrOnMice_' Genotype '_' Unit2plot '_units_MeanFR.mat' ],'Mean_Mtrx');
            boundedline(time,mean2plot,sem2plot,'r','alpha','cmap',hex2rgb(HexColors{R}));
            hold on
            xlim([start_time stop_time]);
            ylim([0 60]);
            title(['Block ' num2str(b)]);
        end
        offset = offset - 10;
        
        if strcmp(Task,'Opto')  || strcmp(Task,'Context') || strcmp(Task,'Passive')
            line([0 0],ylim,'LineWidth',1,'Color',[0 0 0])
        elseif strcmp(Task,'Aversion')
            if b==1 || b==2 || b==3 || b==4
                line([0 0],ylim,'LineWidth',1,'Color',[0 0 0])
            end
            if b==2 || b==4 || b==5 || b==6 || b==7 || b==8
                line([0.7 0.7],ylim,'LineWidth',1,'Color',[0 0 0])
            end
        end
        
    end
    
    ylabel('Firing rate (sd)');
    xlabel('time (s)')
    set(gcf,'units','points','position',[800,900,Block_nb*200,900])
end
if strcmp(Task,'Aversion')
    sgtitle(['PSTH, Aversion ' Genotype ' (n = ' num2str(sum(units_per_region)) ' ' Unit2plot ' units)'])
else
    sgtitle(['PSTH, ' Task ' (n = ' num2str(sum(units_per_region)) ' ' Unit2plot ' units)'])
end

total_per_region = sum(ISregion_mtrx,1);
if task_mod==0
    if Fig_save==1
        saveas(fig2,[FigPath '/' Genotype '_' Unit2plot '_units_PSTH' ],'epsc');
    end
elseif task_mod==1
    if Fig_save==1
        saveas(fig2,[FigPath '/' Genotype '_' Unit2plot '_units_PSTH_modulated_units'],'epsc');
    end
end

% make table with all values
% PFC_Regions =  Regions;
% T = table(PFC_Regions,BM,BS,P1M,P1S,P2M,P2S);
%   if Fig_save==1
%         save([FigPath '/' Genotype '_' Unit2plot '_units_table.mat' ],'T');
%   end

%% Task modulated
fig3=figure;
bar(units_per_region)
set(gca,'XTick',1:numel(Regions))
set(gca,'XTickLabel',Regions);
xtickangle(45)
xlabel('Regions');
ylabel('Number of units');
if task_mod==0
    title(['Units per PFC region (n = ' num2str(sum(units_per_region)) '_' Unit2plot ' units)'])
elseif task_mod==1
    title(['Task-modulated units per PFC region (n = ' num2str(sum(units_per_region)) '_' Unit2plot ' units)'])
end
set(gcf,'units','points','position',[1600,1000,550,420])
if Fig_save==1
    if task_mod==0
        saveas(fig3,[FigPath '/' Genotype '_' Unit2plot '_unit-count'],'epsc');
    elseif task_mod==1
        saveas(fig3,[FigPath '/' Genotype '_' Unit2plot '_units_task-modulated'],'epsc');
    end
end

% fraction
if task_mod==1
    fig3=figure;
    units_per_region = zeros(1,numel(Regions));
    for i = 1:numel(Regions)
        if i ==1
            units_per_region(i) =  h_delim(i);
        else
            units_per_region(i) =  h_delim(i)-h_delim(i-1);
        end
    end
    bar(units_per_region./total_per_region)
    set(gca,'XTick',1:numel(Regions))
    set(gca,'XTickLabel',Regions);
    xtickangle(45)
    ylim([0 1])
    xlabel('Regions');
    ylabel('fraction of units');
    title(['Fraction of task-modulated units per PFC region (n = ' num2str(sum(units_per_region)) ' units)'])
    set(gcf,'units','points','position',[1600,445,550,420])
    if Fig_save==1
        saveas(fig3,[FigPath '/' id{1} '_' id{2} '_' Unit2plot '_units_fraction-task-modulated'],'epsc');
    end
end

%% Latency close up

% psth
fig2=figure;
for b = 1:Block_nb % loop over blocks
    subplot(1,Block_nb,b);
    temp3 = [];
    Z_temp = [];
    %offset = numel(Regions)*10;
    
    for R = 1:numel(Regions) % loop over regions
        temp2 = [];
        
        for f = 1 : CNT-1 % loop over mice
            
            eval(['pp = G.mouse' num2str(f) ';'])
            temp = squeeze(pp(b,R,:,:));
            temp(isnan(temp(:,1)),:) = [];
               if size(temp2,2)==1
                temp2 = temp2';
            end
            temp2 = [temp2; temp];
        end
        if ~isempty(temp2)
            mean2plot = nanmean(temp2)-nanmean(nanmean(temp2(:,1:100),1));%+offset; % baseline correction
            sem2plot = nanstd(temp2)/sqrt(size(temp2,1));
            boundedline(time,mean2plot,sem2plot,'r','alpha','cmap',hex2rgb(HexColors{R}));
            hold on
            xlim([start_time stop_time]);
            %ylim([0 20]);
            title(['Block ' num2str(b)]);
        end
        %offset = offset - 10;
        
        if strcmp(Task,'Opto')  || strcmp(Task,'Context') || strcmp(Task,'Passive')
            line([0 0],ylim,'LineWidth',1,'Color',[0 0 0])
        elseif strcmp(Task,'Aversion')
            if b==1 || b==2 || b==3 || b==4
                line([0 0],ylim,'LineWidth',1,'Color',[0 0 0])
            end
            if b==2 || b==4 || b==5 || b==6 || b==7 || b==8
                line([0.7 0.7],ylim,'LineWidth',1,'Color',[0 0 0])
            end
        end
        
    end
    
    ylabel('Firing rate (sd)');
    xlabel('time (s)')
    set(gcf,'units','points','position',[0,0,Block_nb*400,450])
end
if strcmp(Task,'Aversion')
    sgtitle(['PSTH, Aversion ' Genotype ' (n = ' num2str(sum(units_per_region)) ' units)'])
else
    sgtitle(['PSTH, ' Task ' (n = ' num2str(sum(units_per_region)) ' units)'])
end


%end

