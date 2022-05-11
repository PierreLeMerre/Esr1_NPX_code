%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script that compute the activity modes from nwb file.
% Inspired and adapted from Allen et al., Science, 2019
%
% Written by pielem Le Merre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cls
format long
set(0, 'DefaultFigureRenderer', 'painters')

%% Parameters
% Tasks
Tasks = {'Passive','Context','Aversion'};

Brain_region = 'mPFC';

% Your path to the Database
Path2Data = ['/Volumes/NPX/LeMerre_dataset/new/' Brain_region '/'];
%Path2Data = '/Volumes/Maxtor/NWB/';

% Your path to analysis
Path2Ana = '/Users/pielem/Documents/MATLAB/neuropixelPFC/Matlab/analysis/';

Task = Tasks{3}; % Aversion

% Allen Brain Atlas regions to count the units
Regions = {'ACAd','PL','ILA','ORBm'};
%Regions = {'PL'};


% Min ISI violation rate
MIN_ISI = 0.01; % less than 1% of ISI violations

% Min FR
MIN_FR = 0.1; % mean session FR above 0.1Hz

% Subset to plot
task_mod = 0;
qc = 1; %Dataset passing quality control
subset = 'tmu_OR';
Block_nb = 2; % Block which the subset is from
Block_nb2 = 4; % Combine with another block


% Get nwb files
Path = [Path2Data Task filesep];
D = dir([Path '*.nwb']);

% activity bin parameters
bin_sz = 0.05; % in seconds
sr = 30000; % sampling rate
pre = -2; % in seconds
post = 5; % in seconds
state_window = 2; %in seconds
cue_window = 0.5; %in seconds
aversive_window = 1; %in seconds
time_window = post-pre; % in seconds

% files to plot

Genotype = 'All';

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

%% Compute activity modes
% preallocate coding direction
cd = nan(3,f_stop-f_start,4,round(time_window/bin_sz));
cnt2 = 1;
cd2 = nan(3,f_stop-f_start,4,round(time_window/bin_sz));



for f = f_start : f_stop

    
    
    % nwb file
    nwb = nwbRead([Path D(f).name]);
    % get id and recording duration from meta file
    id = strsplit(nwb.general_session_id,'_');
    disp(['Computing activity modes for mouse ' id{1} ', session ' id{2} '...'])
    
    % Get PFC regions
    ede_regions = nwb.general_extracellular_ephys_electrodes.vectordata.get('location').data.load();
    % Get unit main channel
    main_ch = nwb.units.electrodes.data.load();
    
    isregion = zeros(1,numel(main_ch));
    % loop over regions
    for R = 1:numel(Regions)
        region=char(Regions{R});
        % loop over units
        for i = 1 : numel(main_ch)
            % Make logical region vector
            if (strfind(ede_regions(main_ch(i)+1,:),region)) > 0 % main_ch is 0 indexed!
                isregion(i) = 1;
            end
        end
    end
    
    
    if task_mod == 1
        % Load subset
        load([Path2Ana 'Task_modulated_units/Block' num2str(Block_nb) '/' id{1} '_' id{2} '_' subset '.mat'])
        tmu1 = tmu;
        load([Path2Ana 'Task_modulated_units/Block' num2str(Block_nb2) '/' id{1} '_' id{2} '_' subset '.mat'])
        tmu2 = tmu;
        tmu = unique([tmu1 tmu2]);
    elseif qc == 1
        % Load subset
        load([Path2Ana 'Task_modulated_units/' id{1} '_' id{2} '_qcu.mat'])
        tmu = qcu;
        tmu(isnan(tmu))=[];
    end
    
    if ~isempty(tmu)
    
    % get control trial timestamps (one second before sound)
    trial_ts = nwb.intervals_trials.start_time.data.load();
    if f==3
        blocks = nwb.intervals_trials.vectordata.get('Block').data.load();
        blocks = blocks(1:573);
    elseif f==4
        blocks = nwb.intervals_trials.vectordata.get('Block').data.load();
        blocks = blocks(1:600);
    elseif f==20
        blocks = nwb.intervals_trials.vectordata.get('Block').data.load();
        blocks = blocks(51:end);
    else
        blocks = nwb.intervals_trials.vectordata.get('Block').data.load();
    end
    
    trials1 = trial_ts(blocks==1); % to plot
    trials2 = trial_ts(blocks==2); % to plot
    trials2 = trials2(1:2:end); %one over 2, to have the same dimension as the others
    trials3 = trial_ts(blocks==3); % to plot
    trials4 = trial_ts(blocks==4); % to plot
    

    % state mode
    base0 = trial_ts(blocks==1); % control baseline
    base1 = trial_ts(blocks==4); % aversive state
    
    % cue mode
    cue0 = trial_ts(blocks==1 | blocks==2 | blocks==3 | blocks==4); % 
    cue1 = trial_ts(blocks==1 | blocks==2 | blocks==3 | blocks==4); % 

    % Aversive signal mode
    opto = trial_ts(blocks==2); % internal
    opto = opto(1:2:end); %one over 2, to have the same dimension as the puffs
    %opto = trial_ts(blocks==1);
    puff = trial_ts(blocks==1); % external
   
    
    % Get unit spiketimes
    % Load jagged arrays
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
    qual = nwb.units.vectordata.get('quality').data.load();
    isi = qual(1,:);
    isi_perc = isi./spk_number;
    
    %load('/Users/pielem/Documents/MATLAB/neuropixelPFC/Matlab/analysis/Task_modulated_units/128514_20191215_qc_bolean.mat')
    
    % Get mean session FR
    load(['/Users/pielem/Documents/MATLAB/neuropixelPFC/Matlab/analysis/Mean_FR/' D(f).name '_fr.mat'])
    mean_session_fr = firing_rates;
    
    %% Construct Activity matrices for ROC
    % units X trials x activity bins
    
    % load unit mean firing rate through the session to substract
    load([Path2Ana 'Mean_FR/' id{1} '_' id{2} '_fr.mat'])
    
    %units in PFC
    if task_mod==0 && qc ==0
        units_in_PFC = sum(isregion);
    else
        units_in_PFC = numel(tmu);
    end
    
    genotypes = nan(1,units_in_PFC);
    
    % Full trials
    %bloc1
    T1 = nan(units_in_PFC,numel(trials1),round(time_window/bin_sz));
    cnt = 1;

    % loop over units
    for i = 1 : numel(unit_ids)
        if task_mod==1 || qc ==1
            if ~isempty(tmu)
                if any(tmu==i)
                    spk_ts = unit_times(i-1); % -1 indexed in unit_times!!!
                    [SpikeRates,WindowCenters] = PSTH_Simple(spk_ts.*sr, trials1.*sr, pre, post, sr, bin_sz);
                    T1(cnt,:,:) = SpikeRates';
                    % Make Genotypes array
                    if f<6
                    genotypes(cnt) = 3; % NPY
                    elseif 5<f && f<11
                    genotypes(cnt) = 2;% VGLut2
                    elseif 11<f && f<16
                    genotypes(cnt) = 1;  % WT
                    elseif 15<f
                    genotypes(cnt) = 4;  %Esr
                    end
                    cnt = cnt + 1;                    
                end
            end
        else
            if  isregion(i)==1 && isi_perc(i)<MIN_ISI && mean_session_fr(i)>MIN_FR
                spk_ts = unit_times(i-1); % -1 indexed in unit_times!!!
                [SpikeRates,WindowCenters] = PSTH_Simple(spk_ts.*sr, trials1.*sr, pre, post, sr, bin_sz);
                T1(cnt,:,:) = SpikeRates';
                cnt = cnt + 1;
                U_ID.unit(cnt2) = unit_ids(i);  %% 0 indexed!!
                U_ID.file{cnt2} = [id{1} '_' id{2}];
                cnt2 = cnt2 + 1;
            end
        end
    end
    
    
    %bloc2
    T2 = nan(units_in_PFC,numel(trials2),round(time_window/bin_sz));
    cnt = 1;
    % loop over units
    for i = 1 : numel(unit_ids)
        
        if task_mod==1 || qc ==1
            if ~isempty(tmu)
                if any(tmu==i)
                    spk_ts = unit_times(i-1); % -1 indexed in unit_times!!!
                    [SpikeRates,WindowCenters] = PSTH_Simple(spk_ts.*sr, trials2.*sr, pre, post, sr, bin_sz);
                    T2(cnt,:,:) = SpikeRates';
                    cnt = cnt + 1;
                end
            end
        else
            if  isregion(i)==1 && isi_perc(i)<MIN_ISI && mean_session_fr(i)>MIN_FR
                spk_ts = unit_times(i-1); % -1 indexed in unit_times!!!
                [SpikeRates,WindowCenters] = PSTH_Simple(spk_ts.*sr, trials2.*sr, pre, post, sr, bin_sz);
                T2(cnt,:,:) = SpikeRates';
                cnt = cnt + 1;
            end
        end
    end
    
    
    %bloc3
    T3 = nan(units_in_PFC,numel(trials3),round(time_window/bin_sz));
    cnt = 1;
    % loop over units
    for i = 1 : numel(unit_ids)
        
        if task_mod==1 || qc ==1
            if ~isempty(tmu)
                if any(tmu==i)
                    spk_ts = unit_times(i-1); % -1 indexed in unit_times!!!
                    [SpikeRates,WindowCenters] = PSTH_Simple(spk_ts.*sr, trials3.*sr, pre, post, sr, bin_sz);
                    T3(cnt,:,:) = SpikeRates';
                    cnt = cnt + 1;
                end
            end
        else
            if  isregion(i)==1 && isi_perc(i)<MIN_ISI && mean_session_fr(i)>MIN_FR
                spk_ts = unit_times(i-1); % -1 indexed in unit_times!!!
                [SpikeRates,WindowCenters] = PSTH_Simple(spk_ts.*sr, trials3.*sr, pre, post, sr, bin_sz);
                T3(cnt,:,:) = SpikeRates';
                cnt = cnt + 1;
            end
        end
    end
    
    
    %bloc4
    T4 = nan(units_in_PFC,numel(trials4),round(time_window/bin_sz));
    cnt = 1;
    % loop over units
    for i = 1 : numel(unit_ids)
        if task_mod==1  || qc ==1
        if ~isempty(tmu)
            if any(tmu==i)
                spk_ts = unit_times(i-1); % -1 indexed in unit_times!!!
                [SpikeRates,WindowCenters] = PSTH_Simple(spk_ts.*sr, trials4.*sr, pre, post, sr, bin_sz);
                T4(cnt,:,:) = SpikeRates';
                cnt = cnt + 1;
            end
        else
            if  isregion(i)==1 && isi_perc(i)<MIN_ISI && mean_session_fr(i)>MIN_FR
                spk_ts = unit_times(i-1); % -1 indexed in unit_times!!!
                [SpikeRates,WindowCenters] = PSTH_Simple(spk_ts.*sr, trials4.*sr, pre, post, sr, bin_sz);
                T4(cnt,:,:) = SpikeRates';
                cnt = cnt + 1;
            end
        end
        end
    end
    
    
    
    % Baseline Activity
    % control
    B1 = nan(units_in_PFC,numel(base0),round(state_window/bin_sz));
    B1_labels = zeros(units_in_PFC,numel(base0));
    cnt = 1;
    % loop over units
    for i = 1 : numel(unit_ids)
        if task_mod==1  || qc ==1
        if ~isempty(tmu)
            if any(tmu==i)
                    spk_ts = unit_times(i-1); % -1 indexed in unit_times!!!
                    [SpikeRates,WindowCenters] = PSTH_Simple(spk_ts.*sr, base0.*sr, -state_window, 0, sr, bin_sz);
                    B1(cnt,:,:) = SpikeRates';
                    cnt = cnt + 1;
            end
        else
            if  isregion(i)==1 && isi_perc(i)<MIN_ISI && mean_session_fr(i)>MIN_FR
                spk_ts = unit_times(i-1); % -1 indexed in unit_times!!!
                [SpikeRates,WindowCenters] = PSTH_Simple(spk_ts.*sr, base0.*sr, -state_window, 0, sr, bin_sz);
                B1(cnt,:,:) = SpikeRates';
                cnt = cnt + 1;
            end
        end
        end
    end
    
    % aversive state
    B2 = nan(units_in_PFC,numel(base1),round(state_window/bin_sz));
    B2_labels = ones(units_in_PFC,numel(base1));
    cnt = 1;
    % loop over units
    for i = 1 : numel(unit_ids)
        if task_mod==1  || qc ==1
        if ~isempty(tmu)
            if any(tmu==i)
                    spk_ts = unit_times(i-1); % -1 indexed in unit_times!!!
                    [SpikeRates,WindowCenters] = PSTH_Simple(spk_ts.*sr, base1.*sr, -state_window, 0, sr, bin_sz);
                    B2(cnt,:,:) = SpikeRates';
                    cnt = cnt + 1;
            end
        else
            if  isregion(i)==1 && isi_perc(i)<MIN_ISI && mean_session_fr(i)>MIN_FR
                spk_ts = unit_times(i-1); % -1 indexed in unit_times!!!
                [SpikeRates,WindowCenters] = PSTH_Simple(spk_ts.*sr, base1.*sr, -state_window, 0, sr, bin_sz);
                B2(cnt,:,:) = SpikeRates';
                cnt = cnt + 1;
            end
        end
        end
    end
    
    % Cue Activity
    %control
    C1 = nan(units_in_PFC,numel(cue0),round(cue_window/bin_sz));
    C1_labels = zeros(units_in_PFC,numel(cue0));
    cnt = 1;
    % loop over units
    for i = 1 : numel(unit_ids)
        if task_mod==1  || qc ==1
        if ~isempty(tmu)
            if any(tmu==i)
                    spk_ts = unit_times(i-1); % -1 indexed in unit_times!!!
                    [SpikeRates,WindowCenters] = PSTH_Simple(spk_ts.*sr, cue0.*sr, -cue_window, 0, sr, bin_sz);
                    C1(cnt,:,:) = SpikeRates';
                    cnt = cnt + 1;
            end
        else
            if  isregion(i)==1 && isi_perc(i)<MIN_ISI && mean_session_fr(i)>MIN_FR
                spk_ts = unit_times(i-1); % -1 indexed in unit_times!!!
                [SpikeRates,WindowCenters] = PSTH_Simple(spk_ts.*sr, cue0.*sr, -cue_window, 0, sr, bin_sz);
                C1(cnt,:,:) = SpikeRates';
                cnt = cnt + 1;
            end
        end
        end
    end
    
    %cues
    C2 = nan(units_in_PFC,numel(cue1),round(cue_window/bin_sz));
    C2_labels = ones(units_in_PFC,numel(cue1));
    cnt = 1;
    % loop over units
    for i = 1 : numel(unit_ids)
        if task_mod==1  || qc ==1
        if ~isempty(tmu)
            if any(tmu==i)
                    spk_ts = unit_times(i-1); % -1 indexed in unit_times!!!
                    [SpikeRates,WindowCenters] = PSTH_Simple(spk_ts.*sr, cue1.*sr, 0, cue_window, sr, bin_sz);
                    C2(cnt,:,:) = SpikeRates';
                    cnt = cnt + 1;
            end
        else
            if  isregion(i)==1 && isi_perc(i)<MIN_ISI && mean_session_fr(i)>MIN_FR
                spk_ts = unit_times(i-1); % -1 indexed in unit_times!!!
                [SpikeRates,WindowCenters] = PSTH_Simple(spk_ts.*sr, cue1.*sr, 0, cue_window, sr, bin_sz);
                C2(cnt,:,:) = SpikeRates';
                cnt = cnt + 1;
            end
        end
        end
    end
    
    % Aversive signal mode
    % Opto Activity (internal)
    O = nan(units_in_PFC,numel(opto),round(aversive_window/bin_sz));
    O_labels = zeros(units_in_PFC,numel(opto));
    cnt = 1;
    % loop over units
    for i = 1 : numel(unit_ids)
        if task_mod==1  || qc ==1
        if ~isempty(tmu)
            if any(tmu==i)
                    spk_ts = unit_times(i-1); % -1 indexed in unit_times!!!
                    [SpikeRates,WindowCenters] = PSTH_Simple(spk_ts.*sr, opto.*sr, 0.7, aversive_window+0.7 ,sr, bin_sz);
                    O(cnt,:,:) = SpikeRates';
                    cnt = cnt + 1;
            end
        else
            if  isregion(i)==1 && isi_perc(i)<MIN_ISI && mean_session_fr(i)>MIN_FR
                spk_ts = unit_times(i-1); % -1 indexed in unit_times!!!
                [SpikeRates,WindowCenters] = PSTH_Simple(spk_ts.*sr, opto.*sr, 0.7, aversive_window+0.7 ,sr, bin_sz);
                O(cnt,:,:) = SpikeRates';
                cnt = cnt + 1;
            end
        end
        end
    end
    
    
    % Puff Activity (external)
    P = nan(units_in_PFC,numel(puff),round(aversive_window/bin_sz)); %
    P_labels = ones(units_in_PFC,numel(puff));
    cnt = 1;
    % loop over units
    for i = 1 : numel(unit_ids)
        if task_mod==1  || qc ==1
        if ~isempty(tmu)
            if any(tmu==i)
                    spk_ts = unit_times(i-1); % -1 indexed in unit_times!!!
                    [SpikeRates,WindowCenters] = PSTH_Simple(spk_ts.*sr, puff.*sr, 0.7, aversive_window+0.7, sr, bin_sz);
                    P(cnt,:,:) = SpikeRates';
                    cnt = cnt + 1;
            end
        else
            if  isregion(i)==1 && isi_perc(i)<MIN_ISI && mean_session_fr(i)>MIN_FR
                spk_ts = unit_times(i-1); % -1 indexed in unit_times!!!
                [SpikeRates,WindowCenters] = PSTH_Simple(spk_ts.*sr, puff.*sr, 0.7, aversive_window+0.7, sr, bin_sz);
                P(cnt,:,:) = SpikeRates';
                cnt = cnt + 1;
            end
        end
        end
    end
    
    
    
    %% Obtain AUC form ROC
    % State mode
    AUC = nan(1,units_in_PFC);
    for i = 1 : units_in_PFC
        labels = [B1_labels(i,:) B2_labels(i,:)];
        scores = [(mean(squeeze(B1(i,:,:)),2)-firing_rates(i))' (mean(squeeze(B2(i,:,:)),2)-firing_rates(i))'];
        [~,~,~,AUC(i)] = perfcurve(labels,scores,0);
    end
    AUC = AUC-0.5.*2; % To get negative and positive values :[-1 1] interval
    W1 = AUC;
    
    
    % Cue mode
    AUC = nan(1,units_in_PFC);
    for i = 1 : units_in_PFC
        labels = [C1_labels(i,:) C2_labels(i,:)];
        scores = [(mean(squeeze(C1(i,:,:)),2)-firing_rates(i))' (mean(squeeze(C2(i,:,:)),2)-firing_rates(i))'];
        [~,~,~,AUC(i)] = perfcurve(labels,scores,0);
    end
    AUC = AUC-0.5.*2; % To get negative and positive values :[-1 1] interval
    W2 = AUC;
    
    % Aversive signal mode
    AUC = nan(1,units_in_PFC);
    for i = 1 : units_in_PFC
        labels = [O_labels(i,:) P_labels(i,:)];
        scores = [(mean(squeeze(O(i,:,:)),2)-firing_rates(i))' (mean(squeeze(P(i,:,:)),2)-firing_rates(i))'];
        [~,~,~,AUC(i)] = perfcurve(labels,scores,0);
    end
    AUC = AUC-0.5.*2; % To get negative and positive values :[-1 1] interval
    W3 = AUC;
    
    W = [W1; W2; W3];
    
    %% Obtain weights form SVM
    rng(42); %for reproducibility
    % State mode
    labels = [B1_labels(1,:) B2_labels(1,:)];
    nF = [nanmean(B1,3) nanmean(B2,3)]';
    nF(:,isnan(nF(1,:))) = [];
    Md = fitclinear(nF,labels,'Prior','Uniform');
    W_SVM_1 = Md.Beta';
    
    % Cue mode
    labels = [C1_labels(1,:) C2_labels(1,:)];
    nF = [nanmean(C1,3) nanmean(C2,3)]';
    nF(:,isnan(nF(1,:))) = [];
    Md = fitclinear(nF,labels,'Prior','Uniform');
    W_SVM_2 = Md.Beta';
    
    labels = [O_labels(1,:) P_labels(1,:)];
    nF = [nanmean(O,3) nanmean(P,3)]';
    nF(:,isnan(nF(1,:))) = [];
    Md = fitclinear(nF,labels,'Prior','Uniform');
    W_SVM_3 = Md.Beta';
    
    W_SVM = [W_SVM_1; W_SVM_2; W_SVM_3];
    
    %% Obtain Weightd form average diff
    % State mode
    AVG_Diff = nan(1,units_in_PFC);
    for i = 1 : units_in_PFC
        AVG_Diff(i) = mean(mean(squeeze(B1(i,:,:)),1)-mean(squeeze(B2(i,:,:)),1));
    end
    W1_diff = AVG_Diff;
    W1_diff = W1_diff/norm(W1_diff); % divide by L2 norm
    
    % Cue mode
    AVG_Diff = nan(1,units_in_PFC);
    for i = 1 : units_in_PFC
        AVG_Diff(i) = mean(mean(squeeze(C1(i,:,:)),1)-mean(squeeze(C2(i,:,:)),1));
    end
    W2_diff = AVG_Diff;
    W2_diff = W2_diff/norm(W2_diff); % divide by L2 norm
    
    % Aversive signal mode
    AVG_Diff = nan(1,units_in_PFC);
    for i = 1 : units_in_PFC
        AVG_Diff(i) = mean(mean(squeeze(O(i,:,:)),1)-mean(squeeze(P(i,:,:)),1));
    end
    W3_diff = AVG_Diff;
    W3_diff = W3_diff/norm(W3_diff); % divide by L2 norm
    
    W_diff = [W1_diff; W2_diff; W3_diff];
    
    %% Make Activity Matrix of averaged trial for each neuron
    
    % Gaussian Kernel 50ms for smoothing
    %N = 0.05/bin_sz;
    N = 5;
    g = gausswin(N);
    w=g./trapz(g); % normalize area = 1
    
    X = nan(units_in_PFC,round(time_window/bin_sz*4));
    for i = 1 : units_in_PFC
        X(i,:) = [conv(mean(squeeze(T1(i,:,:)),1)-firing_rates(i),w,'same') conv(mean(squeeze(T2(i,:,:)),1)-firing_rates(i),w,'same')... % 
            conv(mean(squeeze(T3(i,:,:)),1)-firing_rates(i),w,'same') conv(mean(squeeze(T4(i,:,:)),1)-firing_rates(i),w,'same')]; %-firing_rates(i)
    end
    X(isnan(X(:,1)),:) = [];
    
    X_all = [X_all X];
    %% QR decomposition
    % Perform a QR decomposition the m-by-n matrix W such that W = Q*R.
    % The factor R is an m-by-n upper-triangular matrix, and the factor Q is an m-by-m orthogonal matrix.
    

    [Q,R] = qr(W_diff);
    
    [Q2,R2] = qr(W_SVM);

    
    %% Coding direction
    
    % Get Global activity space
    %load('/Users/pielem/Documents/MATLAB/neuropixelPFC/Matlab/analysis/Activity_modes/activity_space_SVM.mat')
    %load('/Users/pielem/Documents/MATLAB/neuropixelPFC/Matlab/analysis/Activity_modes/activity_space_AVG.mat')
    
    CD = R*X;    
    CD2 = R2*X;
    
    
    for dim = 1:3
        cd(dim,f,1,:) = CD(dim,1:time_window/bin_sz);
        cd(dim,f,2,:) = CD(dim,time_window/bin_sz+1:2*time_window/bin_sz);
        cd(dim,f,3,:) = CD(dim,2*time_window/bin_sz+1:3*time_window/bin_sz);
        cd(dim,f,4,:) = CD(dim,3*time_window/bin_sz+1:4*time_window/bin_sz);
        
        cd2(dim,f,1,:) = CD2(dim,1:time_window/bin_sz);
        cd2(dim,f,2,:) = CD2(dim,time_window/bin_sz+1:2*time_window/bin_sz);
        cd2(dim,f,3,:) = CD2(dim,2*time_window/bin_sz+1:3*time_window/bin_sz);
        cd2(dim,f,4,:) = CD2(dim,3*time_window/bin_sz+1:4*time_window/bin_sz);
    end
    

    
   
    end    
    
end


%% 3D Plot
f1 = figure;
plot3(squeeze(nanmean(cd(2,:,1,:),2)),squeeze(nanmean(cd(3,:,1,:),2)),squeeze(nanmean(cd(1,:,1,:),2)),'k');
hold on
plot3(squeeze(nanmean(cd(2,:,2,:),2)),squeeze(nanmean(cd(3,:,2,:),2)),squeeze(nanmean(cd(1,:,2,:),2)),'b');
plot3(squeeze(nanmean(cd(2,:,3,:),2)),squeeze(nanmean(cd(3,:,3,:),2)),squeeze(nanmean(cd(1,:,3,:),2)),'g');
plot3(squeeze(nanmean(cd(2,:,4,:),2)),squeeze(nanmean(cd(3,:,4,:),2)),squeeze(nanmean(cd(1,:,4,:),2)),'r');
xlabel('Cue (A.U.)');
xlim([-40 120])
ylabel('Aversive signal (A.U.)');
ylim([-40 120])
zlabel('State (A.U.)');
zlim([-40 120])
title(['Activity modes (AVG diff), ' Genotype])
grid on
set(gcf,'units','points','position',[0,900,500,450])
saveas(f1,['/Users/pielem/OneDrive - KI.SE/Shared/LHA-LHb-PFC/Meeting11/3D_Activity_modes_(AVG_diff)_'  Genotype '_PL'],'epsc');

f2 =figure;
plot3(squeeze(nanmean(cd2(2,:,1,:),2)),squeeze(nanmean(cd2(3,:,1,:),2)),squeeze(nanmean(cd2(1,:,1,:),2)),'k');
hold on
plot3(squeeze(nanmean(cd2(2,:,2,:),2)),squeeze(nanmean(cd2(3,:,2,:),2)),squeeze(nanmean(cd2(1,:,2,:),2)),'b');
plot3(squeeze(nanmean(cd2(2,:,3,:),2)),squeeze(nanmean(cd2(3,:,3,:),2)),squeeze(nanmean(cd2(1,:,3,:),2)),'g');
plot3(squeeze(nanmean(cd2(2,:,4,:),2)),squeeze(nanmean(cd2(3,:,4,:),2)),squeeze(nanmean(cd2(1,:,4,:),2)),'r');
xlabel('Cue (A.U.)');
xlim([-3.5 3.5])
ylabel('Aversive signal (A.U.)');
ylim([-3.5 3.5])
zlabel('State (A.U.)');
zlim([-3.5 3.5])
title(['Activity modes (SVM), ' Genotype])
grid on
set(gcf,'units','points','position',[0,370,500,450])
saveas(f2,['/Users/pielem/OneDrive - KI.SE/Shared/LHA-LHb-PFC/Meeting11/3D_Activity_modes_(SVM)_'  Genotype '_PL'],'epsc');

%% Indivudual modes
t = linspace(pre,post,round(time_window/bin_sz));
f3 = figure;
subplot(1,3,1)
plot(t,squeeze(nanmean(cd(2,:,1,:),2)),'k')
hold on
plot(t,squeeze(nanmean(cd(2,:,2,:),2)),'b')
plot(t,squeeze(nanmean(cd(2,:,3,:),2)),'g')
plot(t,squeeze(nanmean(cd(2,:,4,:),2)),'r')
title('Cue mode')
xlim([-1 3])
ylim([-40 120])

subplot(1,3,2)
plot(t,squeeze(nanmean(cd(3,:,1,:),2)),'k')
hold on
plot(t,squeeze(nanmean(cd(3,:,2,:),2)),'b')
plot(t,squeeze(nanmean(cd(3,:,3,:),2)),'g')
plot(t,squeeze(nanmean(cd(3,:,4,:),2)),'r')
title('Aversive signal mode')
xlim([-1 3])
ylim([-40 120])

subplot(1,3,3)
plot(t,squeeze(nanmean(cd(1,:,1,:),2)),'k')
hold on
plot(t,squeeze(nanmean(cd(1,:,2,:),2)),'b')
plot(t,squeeze(nanmean(cd(1,:,3,:),2)),'g')
plot(t,squeeze(nanmean(cd(1,:,4,:),2)),'r')
title('State mode')
xlim([-1 3])
ylim([-40 120])

suptitle(['Activity modes (AVG diff), ' Genotype])
set(gcf,'units','points','position',[500,900,2000,450])
saveas(f3,['/Users/pielem/OneDrive - KI.SE/Shared/LHA-LHb-PFC/Meeting11/Activity_modes_(AVG_diff)_'  Genotype '_PL'],'epsc');

f4 = figure;
subplot(1,3,1)
plot(t,squeeze(nanmean(cd2(2,:,1,:),2)),'k')
hold on
plot(t,squeeze(nanmean(cd2(2,:,2,:),2)),'b')
plot(t,squeeze(nanmean(cd2(2,:,3,:),2)),'g')
plot(t,squeeze(nanmean(cd2(2,:,4,:),2)),'r')
title('Cue mode')
xlim([-1 3])
ylim([-3.5 3.5])

subplot(1,3,2)
plot(t,squeeze(nanmean(cd2(3,:,1,:),2)),'k')
hold on
plot(t,squeeze(nanmean(cd2(3,:,2,:),2)),'b')
plot(t,squeeze(nanmean(cd2(3,:,3,:),2)),'g')
plot(t,squeeze(nanmean(cd2(3,:,4,:),2)),'r')
title('Aversive signal mode')
xlim([-1 3])
ylim([-3.5 3.5])

subplot(1,3,3)
plot(t,squeeze(nanmean(cd2(1,:,1,:),2)),'k')
hold on
plot(t,squeeze(nanmean(cd2(1,:,2,:),2)),'b')
plot(t,squeeze(nanmean(cd2(1,:,3,:),2)),'g')
plot(t,squeeze(nanmean(cd2(1,:,4,:),2)),'r')
title('State mode')
xlim([-1 3])
ylim([-3.5 3.5])

suptitle(['Activity modes (SVM), ' Genotype])
set(gcf,'units','points','position',[500,370,2000,450])
%saveas(f4,['/Users/pielem/OneDrive - KI.SE/Shared/LHA-LHb-PFC/Meeting11/Activity_modes_(SVM)_'  Genotype '_PL'],'epsc');
