%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script that compute PCA and hierachical clustering of trial
% activity from nwb file.
%
% Written by pierre Le Merre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cls
format long
set(0, 'DefaultFigureRenderer', 'painters')

%% Parameters

% Your path to the Database
Path2Data = '/Volumes/labs-2/dmclab/Pierre/NPX_Database/mPFC/Aversion/';


% Your path to analysis folder
BasePath = '/Users/pierre/Documents/MATLAB/';
Path2Ana = [BasePath 'Code_10122021/neuropixelPFC/Matlab/analysis/'];

% Load colormaps
load([BasePath 'Esr1_NPX_code/utilities/Colormaps/sd_colormap.mat'])
load([BasePath 'Esr1_NPX_code/utilities/Colormaps/cmap_genotype2.mat'])
load([BasePath 'Esr1_NPX_code/utilities/Colormaps/cmap_allen.mat'])

Task = 'Aversion';

% Allen Brain Atlas regions to count the units
Regions = {'ACAd','PL','ILA','ORBm'};


% Get nwb files
D = dir([Path2Data '*.nwb']);

% Subset to plot
plotsubset = 1; % bolean to plot a subset
glm = 1; % the subset is glm fitted units

% Min ISI violation rate
MIN_ISI = 0.01; % less than 1% of ISI violations

% Min FR
MIN_FR = 0.1; % mean session FR above 0.1Hz

% activity bin parameters
bin_sz = 0.01; % in seconds
sr = 30000; % sampling rate
pre = -1; % in seconds
post = 3; % in seconds
state_window = 1; %in seconds
cue_window = 0.5; %in seconds
aversive_window = 1; %in seconds
time_window = post-pre; % in seconds

% Gaussian kernel size for spatial smoothing in um
win_sz=10;

% files to plot
%Chose genotype
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


%% Preallocate

% initialize
P = nan(3,f_stop-f_start-5,4,round(time_window/bin_sz));

% Preallocation
disp('preallocating...')
cnt2 = 1;
Tot_units = 0;
U_ID = []; %Store unit and file origin

for f = f_start : f_stop
    if strcmp(Genotype,'All') && 15<f && f<22
        
    else
        % nwb file
        nwb = nwbRead([Path2Data D(f).name]);
        % get id and recording duration from meta file
        id = strsplit(nwb.general_session_id,'_');
        
        if plotsubset==1
            % GLM subset
            if glm ==1
                load([Path2Ana 'Task_modulated_GLM/' id{1} '_' id{2} '_tmu.mat'])
            else
                disp('Subset cannot be found')
            end
        end
        
        for i = 1 : numel(tmu)
            U_ID.unit{cnt2} = tmu(i); %% 0 indexed!!
            U_ID.file{cnt2} = [id{1} '_' id{2}];
            cnt2 = cnt2 + 1;
        end
        
        Tot_units = Tot_units + numel(tmu);
    end
end

%% Run
bin_nb = (post-pre)/bin_sz*4;
X_all = nan(Tot_units,bin_nb);
FR_all = nan(Tot_units,bin_nb);
genotype_all = nan(1,Tot_units);
modul_all = nan(4,Tot_units);
modul_all2 = nan(4,Tot_units);
Active_bins = nan(1,Tot_units);
AP = nan(1,Tot_units);
ML = nan(1,Tot_units);
DV = nan(1,Tot_units);
RG = cell(1,Tot_units);
in = 1;
out = 0; % gets updated in the loop

for f = f_start : f_stop
    
    if strcmp(Genotype,'All') && 15<f && f<22
        
    else
        % nwb file
        nwb = nwbRead([Path2Data D(f).name]);
        % get id and recording duration from meta file
        id = strsplit(nwb.general_session_id,'_');
        
        disp(['preprocessing mouse ' id{1} ', session ' id{2} '...'])
        
        
        % Load subset
        
        if plotsubset==1
            
            % GLM subset
            if glm ==1
                load([Path2Ana 'Task_modulated_GLM/' id{1} '_' id{2} '_tmu.mat'])
                
            else
                disp('Subset cannot be found')
            end
            
        end
        
        
        % Get PFC regions
        ede_regions = nwb.general_extracellular_ephys_electrodes.vectordata.get('location').data.load();
        
        % Get unit main channel
        main_ch = nwb.units.electrodes.data.load();
        
        % Get Sterotaxic coordinates
        ap = nwb.general_extracellular_ephys_electrodes.vectordata.get('AP').data.load();
        ml = nwb.general_extracellular_ephys_electrodes.vectordata.get('ML').data.load();
        dv = nwb.general_extracellular_ephys_electrodes.vectordata.get('DV').data.load();
        
        
        isregion = zeros(1,numel(main_ch));
        % loop over regions
        for R = 1:numel(Regions)
            region=char(Regions{R});
            % loop over units
            for i = 1 : numel(main_ch)
                % Make logical region vector
                if (strfind(ede_regions{main_ch(i)+1,:},region)) > 0 % main_ch is 0 indexed!
                    isregion(i) = 1;
                end
            end
        end
        
        % get control trial timestamps (one second before sound)
        trial_ts = nwb.intervals_trials.start_time.data.load();
        if f==3
            blocks = nwb.intervals_trials.vectordata.get('Block').data.load();
            blocks = blocks(1:573);
        elseif f==4
            blocks = nwb.intervals_trials.vectordata.get('Block').data.load();
            blocks = blocks(1:600);
        elseif f==27
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
        
        %Genotype
        if 0<f && f<6 % NPY
            gen = 3;
        elseif 5<f && f<11 % VGlut2
            gen = 2;
        elseif 10<f && f<16 % WT
            gen = 1;
        elseif 22<f && f<32 % Esr
            gen = 4;
        end
        
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
        
        % Get mean session FR
        load([Path2Ana 'Mean_FR/' id{1} '_' id{2} '.nwb_fr.mat'])
        mean_session_fr = firing_rates;
        
        if plotsubset==1
            if glm ==1
                load([Path2Ana '/Task_modulated_GLM/' id{1} '_' id{2} '_tuning_scores.mat'])
                load([Path2Ana '/Task_modulated_GLM/' id{1} '_' id{2} '_tuned_bol.mat'])
                maxloc1 = zeros(1,size(tun_scores,2));
                maxloc2 = zeros(1,size(tun_scores,2));
                for i = 1 : size(tun_scores,2)
                    if any(tun_bol(:,i)==1 | tun_bol(:,i)==-1)
                        [~,maxloc1(i)] = max(tun_scores(:,i));
                        tun_scores(maxloc1(i),i) = -Inf;
                        [~,maxloc2(i)] = max(tun_scores(:,i));
                    end
                end
                bol_s1 = maxloc1==1;
                bol_o = maxloc1==2;
                bol_s2 = maxloc1==3;
                bol_a = maxloc1==4;
                
                bol2_s1 = maxloc2==1;
                bol2_o = maxloc2==2;
                bol2_s2 = maxloc2==3;
                bol2_a = maxloc2==4;
                
            end
        end
        
        
        % Construct Activity matrices for PCA
        % units X activity bins (4 trials)
        
        % load unit mean firing rate through the session to substract
        
        %units in PFC
        if task_mod==0
            units_in_PFC = sum(isregion);
        elseif task_mod==1
            units_in_PFC = numel(tmu);
        end
        
        disp([id{1} '_' id{2} ': ' num2str(units_in_PFC)])
        % update out counter
        out = out + units_in_PFC;
        
        % Full trials
        %bloc1
        T1 = nan(units_in_PFC,numel(trials1),round(time_window/bin_sz));
        cnt = 1;
        
        % loop over units
        for i = 1 : numel(unit_ids)
            if ~isempty(tmu)
                %for j = 1 : numel(tmu)
                if any(tmu==i) %unit_ids(i) == tmu(j)+1
                    spk_ts = unit_times(i-1); % -1 indexed in unit_times!!!
                    [SpikeRates,WindowCenters] = PSTH_Simple(spk_ts.*sr, trials1.*sr, pre, post, sr, bin_sz);
                    T1(cnt,:,:) = SpikeRates';
                    cnt = cnt + 1;
                end
                %end
            end
        end
        
        
        %bloc2
        T2 = nan(units_in_PFC,numel(trials2),round(time_window/bin_sz));
        cnt = 1;
        % loop over units
        for i = 1 : numel(unit_ids)
            if ~isempty(tmu)
                %for j = 1 : numel(tmu)
                if any(tmu==i) %unit_ids(i) == tmu(j)+1
                    spk_ts = unit_times(i-1); % -1 indexed in unit_times!!!
                    [SpikeRates,WindowCenters] = PSTH_Simple(spk_ts.*sr, trials2.*sr, pre, post, sr, bin_sz);
                    T2(cnt,:,:) = SpikeRates';
                    cnt = cnt + 1;
                end
                %end
            end
        end
        
        
        %bloc3
        T3 = nan(units_in_PFC,numel(trials3),round(time_window/bin_sz));
        cnt = 1;
        % loop over units
        for i = 1 : numel(unit_ids)
            
            if ~isempty(tmu)
                %for j = 1 : numel(tmu)
                if any(tmu==i) %unit_ids(i) == tmu(j)+1
                    spk_ts = unit_times(i-1); % -1 indexed in unit_times!!!
                    [SpikeRates,WindowCenters] = PSTH_Simple(spk_ts.*sr, trials3.*sr, pre, post, sr, bin_sz);
                    T3(cnt,:,:) = SpikeRates';
                    cnt = cnt + 1;
                end
                %end
            end
        end
        
        
        %bloc4
        T4 = nan(units_in_PFC,numel(trials4),round(time_window/bin_sz));
        cnt = 1;
        % loop over units
        for i = 1 : numel(unit_ids)
            if ~isempty(tmu)
                %for j = 1 : numel(tmu)
                if any(tmu==i) %unit_ids(i) == tmu(j)+1
                    spk_ts = unit_times(i-1); % -1 indexed in unit_times!!!
                    [SpikeRates,WindowCenters] = PSTH_Simple(spk_ts.*sr, trials4.*sr, pre, post, sr, bin_sz);
                    T4(cnt,:,:) = SpikeRates';
                    cnt = cnt + 1;
                end
                %end
            end
        end
        
        
        unit_genotype = nan(1,units_in_PFC);
        unit_modul = nan(4,units_in_PFC);
        unit_modul2 = nan(4,units_in_PFC);
        aap = nan(1,units_in_PFC);
        mml = nan(1,units_in_PFC);
        ddv = nan(1,units_in_PFC);
        reg = cell(1,units_in_PFC);
        cnt = 1;
        % loop over units
        for i = 1 : numel(unit_ids)
            if ~isempty(tmu)
                for j = 1 : numel(tmu)
                    if i == tmu(j)
                        %genotype
                        unit_genotype(cnt) = gen;
                        if plotsubset==1
                            % GLM tuning score
                            if glm == 1
                                for k = 1 : size(tun_scores,2)
                                    % primary tuning
                                    if i == tmu(k)+1 && bol_s1(k)==1 %% tmu is 0 indexed here!!
                                        unit_modul(1,cnt) = 1;
                                    end
                                    if i == tmu(k)+1 && bol_o(k)==1
                                        unit_modul(2,cnt) = 2;
                                    end
                                    if i == tmu(k)+1 && bol_s2(k)==1
                                        unit_modul(3,cnt) = 3;
                                    end
                                    if i == tmu(k)+1 && bol_a(k)==1
                                        unit_modul(4,cnt) = 4;
                                    end
                                    % secondary tuning
                                    if i == tmu(k)+1 && bol2_s1(k)==1 %% tmu is 0 indexed here!!
                                        unit_modul2(1,cnt) = 1;
                                    end
                                    if i == tmu(k)+1 && bol2_o(k)==1
                                        unit_modul2(2,cnt) = 2;
                                    end
                                    if i == tmu(k)+1 && bol2_s2(k)==1
                                        unit_modul2(3,cnt) = 3;
                                    end
                                    if i == tmu(k)+1 && bol2_a(k)==1
                                        unit_modul2(4,cnt) = 4;
                                    end
                                end
                            end
                        end
                        % stereotaxic cordinates
                        aap(cnt) = ap(main_ch(i)+1);
                        mml(cnt) = ml(main_ch(i)+1);
                        ddv(cnt) = dv(main_ch(i)+1);
                        reg{cnt} = ede_regions(main_ch(i)+1,:);
                        cnt = cnt + 1;
                    end
                end             
            end
        end
        
        
        % Make Activity Matrix of averaged trial for each neuron
        % Gaussian Kernel 50ms for smoothing
        N = 5;
        g = gausswin(N);
        w=g./trapz(g); % normalize area = 1
        
        X = nan(units_in_PFC,round(time_window/bin_sz*4));
        FR = nan(units_in_PFC,round(time_window/bin_sz*4));
        active_bins = nan(1,units_in_PFC);
        for i = 1 : units_in_PFC
            X(i,:) = [conv(mean(squeeze(T1(i,:,:)),1),w,'same') conv(mean(squeeze(T2(i,:,:)),1),w,'same')...
                conv(mean(squeeze(T3(i,:,:)),1),w,'same') conv(mean(squeeze(T4(i,:,:)),1),w,'same')];
            % Nb of active bins
            active_bins(i) = sum(abs(X(i,:))>0.5);  % 0.5Hz as threshold
            % Store FR
            FR(i,:) = X(i,:);
            % Z-score the data
            X(i,:) = zscore(X(i,:));
            
        end
        
        if ~isempty(tmu)
            X_all(in:out,:) = X;
            FR_all(in:out,:) = FR;
            Active_bins(in:out) = active_bins;
            genotype_all(in:out) = unit_genotype;
            modul_all(:,in:out) = unit_modul;
            modul_all2(:,in:out) = unit_modul2;
            AP(in:out) = aap;
            ML(in:out) = mml;
            DV(in:out) = ddv;
            for ii = in : out
                RG{ii} = reg{ii-in+1};
            end
            
            in = in + units_in_PFC;
            
            
        end
        
        
    end
end

%% Quality Control
% histogram of activity bin
%histogram(Active_bins,100);
% Set cutoff for the data more than 200 active bins over 1600.

disp('Quality control...')

X_all = X_all(Active_bins>200,:);
FR_all = FR_all(Active_bins>200,:);
genotype_all = genotype_all(Active_bins>200);
modul_all = modul_all(:,Active_bins>200);
modul_all2 = modul_all2(:,Active_bins>200);
AP = AP(Active_bins>200);
ML = ML(Active_bins>200);
DV = DV(Active_bins>200);
RG = RG(Active_bins>200);
G_ID.unit = U_ID.unit(Active_bins>200);
G_ID.file = U_ID.file(Active_bins>200);

%remove nan
FR_all(:,isnan(X_all(:,1))) = [];
genotype_all(isnan(X_all(:,1))) = [];
modul_all2(:,isnan(X_all(:,1))) = [];
AP(isnan(X_all(:,1))) = [];
ML(isnan(X_all(:,1))) = [];
DV(isnan(X_all(:,1))) = [];
RG(isnan(X_all(:,1))) = [];
G_ID.unit(isnan(X_all(:,1))) = [];
G_ID.file(isnan(X_all(:,1))) = [];
X_all(isnan(X_all(:,1))) = [];

%% Dimensionality reduction: PCA

disp('PCA...')

[coef,score,latent,tsquared,explained,mu] = pca(X_all');

for dim = 1:20
    P(dim,f,1,:) = score(1:time_window/bin_sz,dim)';
    P(dim,f,2,:) = score(time_window/bin_sz+1:2*time_window/bin_sz,dim)';
    P(dim,f,3,:) = score(2*time_window/bin_sz+1:3*time_window/bin_sz,dim)';
    P(dim,f,4,:) = score(3*time_window/bin_sz+1:4*time_window/bin_sz,dim)';
end

% PLot activity space
figure;
plot3(squeeze(nanmean(P(1,:,1,:),2)),squeeze(nanmean(P(2,:,1,:),2)),squeeze(nanmean(P(3,:,1,:),2)),'k');
hold on
plot3(squeeze(nanmean(P(1,:,2,:),2)),squeeze(nanmean(P(2,:,2,:),2)),squeeze(nanmean(P(3,:,2,:),2)),'b');
plot3(squeeze(nanmean(P(1,:,3,:),2)),squeeze(nanmean(P(2,:,3,:),2)),squeeze(nanmean(P(3,:,3,:),2)),'g');
plot3(squeeze(nanmean(P(1,:,4,:),2)),squeeze(nanmean(P(2,:,4,:),2)),squeeze(nanmean(P(3,:,4,:),2)),'r');
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
title(['PCA, ' Genotype])
grid on

%explained variance, elbow plot
figure;
plot(explained(1:20),'ko-')
hold on
plot(cumsum(explained),'go:')
xlim([1 20])
set(gca,'XTick',1:length(explained(1:20)))
xlabel('PCs');
ylabel('Explained variance');

Feat = coef(:,1:20);

%% Hierarchical Clustering
%Find the similarity or dissimilarity between every pair of objects in the data set
CLU_NB = 15;

disp('Hierarchical clustering...')

%rng(42);  % For reproducibility
Y = pdist(Feat);

%linkage: Group the objects into a binary, hierarchical cluster tree.
Z = linkage(Y,'ward');

figure;
[H,T,outperm] = dendrogram(Z,CLU_NB); %


%sort clusters by number of units
nb_of_units_in_clu = zeros(1,numel(unique(T)));
for i = 1 : numel(unique(T))
    nb_of_units_in_clu(i) = sum(T==i);
end
[~,order] = sort(nb_of_units_in_clu);
clu_list = 1:numel(unique(T));
clu_list = flip(clu_list(order));

% reorder units
[~,or]=sort(T);
figure; imagesc(X_all(or,:));
h = 0;
for i=1:numel(nb_of_units_in_clu)-1
    h = h + nb_of_units_in_clu(i);
    line([0 1600],[h h],'LineWidth',1,'Color',[0 0 0])
end
colormap(sd_colormap);
caxis([-2, 2])
c=colorbar('southoutside');
c.Label.String = 'Firing rate (sd)';


%% Reorganize Matrix by genotypes c(v2, no tree)

outperm2 = [11 4 10 5 14 1 8 6 12 7 2 9 3 13 15];

t = linspace(pre,post,time_window/bin_sz);
if strcmp(Genotype,'All')
    ConcX = [];
    ConcG = [];
    ConcM = [];
    ConcM2 = [];
    modul_frac = nan(CLU_NB,4);
    comodul_frac = nan(CLU_NB,1);
    sz_array = nan(CLU_NB,4);
    
    array4stats = nan(CLU_NB,14);
    
    for i = 1 : CLU_NB
        %subplot(10,4,1+cnt:3+cnt)
        x2plot = X_all(T==outperm2(i),:);
        fr2plot = FR_all(T==outperm2(i),:);
        gen2plot = genotype_all(T==outperm2(i))';
        modul2plot = modul_all(:,T==outperm2(i))';
        modul2plot2 = modul_all2(:,T==outperm2(i))';
        modul_frac(i,:) = [sum(modul_all(1,T==outperm2(i))==1) sum(modul_all(2,T==outperm2(i))==2)...
            sum(modul_all(3,T==outperm2(i))==3) sum(modul_all(4,T==outperm2(i))==4)];
        comodul_frac(i,1) = sum(modul_all(2,T==outperm2(i))==2 & modul_all(4,T==outperm2(i))==4);
        reg2plot = RG(T==outperm2(i));
        [~,order]=sort(gen2plot);
        ConcX = [ConcX; x2plot(order,:)];
        %imagesc(x2plot(order,:))
        ylabel(['Cluster ' num2str(outperm2(i))] )
        %subplot(10,4,4+cnt)
        ConcG = [ConcG; gen2plot(order)];
        ConcM = [ConcM; modul2plot(order,:)];
        ConcM2 = [ConcM2; modul2plot2(order,:)];
        %imagesc(gen2plot(order))
        %cnt = cnt + 4;
        
        % single cluster plots
        fg = figure;
        in = 1;
        out = time_window/bin_sz;
        step = time_window/bin_sz;
        baseline = nan(4,size(fr2plot,1));
        in2 = 1;
        out2 = time_window/bin_sz/4;
        step2 = time_window/bin_sz;
        cnt = 0;
        for j = 1:4
            eval(['ax(' num2str(j) ') = subplot(2,12,1+cnt:2+cnt);'])
            imagesc('XData',t,'CData',x2plot(order,in:out))
            %imagesc('XData',t,'CData',x2plot(order,in2:out2))
            colormap(sd_colormap);
            caxis([-3, 3])
            ylim([0 size(x2plot(order,:),1)])
            xlim([-1 3])
            title(['Block ' num2str(j)])
            line([0 0],[0 size(x2plot(order,:),1)],'LineWidth',1,'Color',[0 0 0])
            colorbar
            %freezeColors
            
            eval(['ax(' num2str(j+12) ') = subplot(2,12,13+cnt:14+cnt);'])
            boundedline(t,nanmean(fr2plot(order,in:out)),nanstd(fr2plot(order,in:out))./sqrt(size(fr2plot,1)),'k')
            ylim([-1 max(nanmean(fr2plot(order,:)))+1])
            xlim([-1 3])
            line([0 0],[0 max(nanmean(fr2plot(order,:)))+1],'LineWidth',1,'Color',[0 0 0])
            if j == 2 || j ==4
                line([0.7 0.7],[0 max(nanmean(fr2plot(order,:)))+1],'LineWidth',1,'Color',[0 0 0])
            end
            baseline(j,:) = nanmean(fr2plot(order,in2:out2),2);
            
            in = in + step;
            out = out + step;
            in2 = in2 + step2;
            out2 = out2 + step2;
            cnt = cnt + 2;
            
        end
        subplot(2,12,9);
        imagesc(flip(gen2plot(order,:)));
        if numel(unique(gen2plot(order,:)))==4
            colormap(cmap_genotype);
        elseif numel(unique(gen2plot(order,:)))==3
            colormap(cmap_genotype(2:4,:));
        end
        title('genotype')
        %freezeColors
        
        subplot(2,12,10);
        imagesc(flip(modul2plot(order,:)));
        colormap(inferno);
        caxis([0 4])
        title('primary tuning')
        set(gca,'XTick',1:4)
        xticklabels({'Sound 1', 'Opto', 'Sound2', 'Airpuff'})
        xtickangle(45)
        %freezeColors
        
        subplot(2,12,11);
        imagesc(flip(modul2plot2(order,:)));
        colormap(inferno);
        caxis([0 4])
        title('secondary tuning')
        set(gca,'XTick',1:4)
        xticklabels({'Sound 1', 'Opto', 'Sound2', 'Airpuff'})
        xtickangle(45)
        %freezeColors
        
        subplot(2,12,12);
        nonan = modul2plot(order,:);
        nonan(isnan(nonan))=0;
        nonan2 = modul2plot2(order,:);
        nonan2(isnan(nonan2))=0;
        sumnan = nonan + nonan2;
        comodul2plot = zeros(size(modul2plot(order,:),1),1);
        %comodul2plot(sum(sumnan,2)<6)=0;
        comodul2plot(sum(sumnan,2)==6)=2;
        %comodul2plot(sum(sumnan,2)>6)=0;
        imagesc(flip(comodul2plot));
        title('comodulation')
        set(gca,'XTick',1)
        xticklabels({'Opto-Airpuff'})
        xtickangle(45)
        caxis([0 6])
        colormap(inferno);
        
        subplot(2,12,21:24)
        %Wilcoxon signed rank
        p1 = ranksum(baseline(1,:),baseline(2,:))*6;
        p2 = ranksum(baseline(1,:),baseline(3,:))*6;
        p3 = ranksum(baseline(1,:),baseline(4,:))*6;
        p4 = ranksum(baseline(2,:),baseline(3,:))*6;
        p5 = ranksum(baseline(2,:),baseline(4,:))*6;
        p6 = ranksum(baseline(3,:),baseline(4,:))*6;
        
        array4stats(i,:) = [nanmean(baseline,2)' nanstd(baseline,[],2)' p1 p2 p3 p4 p5 p6];
        %array4stats = array4stats(:,[1 5 2 6 9 3 7 10 11 4 8 12 13 14]);
        
        bar(nanmean(baseline,2))
        hold on
        errorbar(nanmean(baseline,2),nanstd(baseline,[],2),'ko')
        ylabel('firing rate (Hz)')
        ylim([-5 20])
        xlim([0.5 4.5])
        set(gca,'XTick',1:4)
        set(gca,'XTickLabel',{'1', '2', '3', '4'});
        xlabel('blocks')
        title(['Baseline FR, pvalues: ' num2str(round(p1,3)) ', ' num2str(round(p2,3)) ', ' num2str(round(p3,3))...
            ', ' num2str(round(p4,3)) ', ' num2str(round(p5,3)) ', ' num2str(round(p6,3))])
        
        sgtitle(['Cluster ' num2str(outperm2(i))])
        set(gcf,'units','points','position',[0,800,1600,500])
        %saveas(fg,['/Users/pierre/OneDrive - KI.SE/Pierre_Shared/LHA-LHb-PFC/NatureSubmission/clusters/cluster' num2str(i)],'epsc')
        
        % Array for region dot plot
        max_Sz = numel(reg2plot);
        isregion2 = zeros(numel(Regions),max_Sz);
        % loop over regions
        for R = 1:numel(Regions)
            region=char(Regions{R});
            % loop over units
            for iii = 1 : max_Sz
                % Make logical region vector
                if ~isempty((strfind(reg2plot{iii}{1,1},region))) % main_ch is 0 indexed!
                    isregion2(R,iii) = 1;
                end
            end
        end
        sz_array(i,:) = (sum(isregion2,2)./max_Sz);
        
    end
    
    % global plot
    fig = figure;
    in = 1;
    out = time_window/bin_sz;
    step = time_window/bin_sz;
    cnt = 0;
    for j = 1:4
        eval(['ax(' num2str(j) ') = subplot(1,12,1+cnt:2+cnt);'])
        imagesc('XData',t,'CData',ConcX(:,in:out))
        h = 0;
        for i=1:numel(nb_of_units_in_clu)-1
            h = h + nb_of_units_in_clu(outperm2(i));
            line([-1 3],[h h],'LineWidth',1,'Color',[0 0 0])
        end
        line([0 0],[0 size(ConcX,1)],'LineWidth',1,'Color',[0 0 0])
        colormap(sd_colormap);
        caxis([-3, 3])
        ylim([0 size(ConcX,1)])
        xlim([-1 3])
        title(['Block ' num2str(j)])
        in = in + step;
        out = out + step;
        cnt = cnt + 2;
        freezeColors
    end
    
    
    ax5 = subplot(1,12,9);
    imagesc(flip(ConcG));
    h = 0;
    for i=1:numel(nb_of_units_in_clu)-1
        h = h + nb_of_units_in_clu(outperm2(numel(nb_of_units_in_clu)-i+1));
        line([0.5 1.5],[h h],'LineWidth',1,'Color',[1 1 1])
    end
    colormap(cmap_genotype);
    set(gca,'XTick',1)
    xticklabels({''})
    title('genotype')
    freezeColors
    
    ax6 = subplot(1,12,10);
    imagesc(flip(ConcM));
    h = 0;
    for i=1:numel(nb_of_units_in_clu)-1
        h = h + nb_of_units_in_clu(outperm2(numel(nb_of_units_in_clu)-i+1));
        line([0.5 4.5],[h h],'LineWidth',1,'Color',[1 1 1])
    end
    colormap(inferno);
    caxis([0 4])
    set(gca,'XTick',1:4)
    xticklabels({'Sound 1', 'Opto', 'Sound2', 'Airpuff'})
    xtickangle(45)
    title('primary tuning')
    freezeColors
    
    ax7 = subplot(1,12,11);
    imagesc(flip(ConcM2));
    h = 0;
    for i=1:numel(nb_of_units_in_clu)-1
        h = h + nb_of_units_in_clu(outperm2(numel(nb_of_units_in_clu)-i+1));
        line([0.5 4.5],[h h],'LineWidth',1,'Color',[1 1 1])
    end
    colormap(inferno);
    caxis([0 4])
    set(gca,'XTick',1:4)
    xticklabels({'Sound 1', 'Opto', 'Sound2', 'Airpuff'})
    xtickangle(45)
    title('secondary tuning')
    freezeColors
    
    ax8 = subplot(1,12,12);
    nonan = ConcM;
    nonan(isnan(nonan))=0;
    nonan2 = ConcM2;
    nonan2(isnan(nonan2))=0;
    sumnan = nonan + nonan2;
    comodul2plot = zeros(size(ConcM,1),1);
    %comodul2plot(sum(sumnan,2)<6)=0;
    comodul2plot(sum(sumnan,2)==6)=2;
    %comodul2plot(sum(sumnan,2)>6)=0;
    imagesc(flip(comodul2plot));
    h = 0;
    for i=1:numel(nb_of_units_in_clu)-1
        h = h + nb_of_units_in_clu(outperm2(numel(nb_of_units_in_clu)-i+1));
        line([0.5 1.5],[h h],'LineWidth',1,'Color',[1 1 1])
    end
    title('comodulation')
    set(gca,'XTick',1)
    xticklabels({'Opto-Airpuff'})
    xtickangle(45)
    caxis([0 6])
    colormap(inferno);
    
end

set(gcf,'units','points','position',[0,800,1200,1600])

%%
figure;
for i = 1 : CLU_NB
    subplot(CLU_NB,1,CLU_NB-i+1)
    bar([modul_frac(i,:)./nb_of_units_in_clu(outperm2(i))])
    xticklabels({'Sound 1', 'Opto', 'Sound 2', 'Airpuff'})
    ylim([0 0.6])
end
sgtitle('Primary tunings per cluster')
set(gcf,'units','points','position',[1200,800,400,1600])

%%
%remove zeros
sz_array(sz_array==0) = NaN;
figure;
for i = 1 : CLU_NB
    subplot(CLU_NB,1,CLU_NB-i+1)
    scatter([1 2 3 4],[1 1 1 1],sz_array(i,:).*1000,cmap_genotype,'filled')
    set(gca,'XTick',[1 2 3 4]);
    xticklabels({'ACAd', 'PL', 'ILA', 'ORBm'})
    xlim([0 4.5])
end
sgtitle('Region per cluster')
set(gcf,'units','points','position',[1200,800,400,1600])


%% Save Unit matrix

% if glm == 1
%     G_ID.cluster = T';
%     G_ID.outperm = outperm;
%     G_ID.ap = AP;
%     G_ID.ml = ML;
%     G_ID.dv = DV;
%     G_ID.region = RG;
%     G_ID.genotype = genotype_all;
%     save([Path2Ana 'Clustering/G_ID.mat'],'G_ID')
% else
%     U_ID.cluster = T';
%     U_ID.outperm = outperm;
%     U_ID.ap = AP;
%     U_ID.ml = ML;
%     U_ID.dv = DV;
%     U_ID.region = RG;
%     U_ID.genotype = genotype_all;
%     save([Path2Ana 'Clustering/U_ID.mat'],'U_ID')
% end




%% plot anat for cluster

% if plot_anat
%     
%     [av,tv,st] = get_Allen_Data;
%     bregma = allenCCFbregma;
%     bregma=[bregma(1) bregma(3) bregma(2)];
%     
%     figure;
%     [im] = plotAPSlice(bregma(1)-220,av,tv,st);
%     hh1 = imagesc(im);
%     colormap(cmap_allen)
%     freezeColors
%     g=get(hh1);
%     imwrite (g.CData, '/Users/pierre/OneDrive - KI.SE/Shared/LHA-LHb-PFC/Meeting6/coronal.png');
%     close
%     
%     
%     figure;
%     [imm] = plotMLSlice(bregma(2)+30,av,tv,st);
%     hh3 = imagesc(imm);
%     colormap(cmap_allen)
%     freezeColors
%     g2=get(hh3);
%     imwrite (g2.CData, '/Users/pierre/OneDrive - KI.SE/Shared/LHA-LHb-PFC/Meeting6/sagital.png');
%     close
%     
%     
%     
%     for cc = 1 : 10
%         %cc = 1
%         figure
%         clu2plot = outperm(cc);
%         cnt = 1;
%         
%         for gg = 1:4
%             log1 = genotype_all==gg;
%             log2 = T==clu2plot;
%             log2 = log2';
%             ap2plot = bregma(1)+AP(log2 & log1).*100;
%             ml2plot = bregma(2)-ML(log2 & log1).*100;
%             dv2plot = abs(bregma(3)-DV(log2 & log1).*100);
%             if gg==1
%                 load('/Users/pierre/Documents/MATLAB/neuropixelPFC/Matlab/Utilities/Colormaps/cmap_WT.mat');
%             elseif gg==2
%                 load('/Users/pierre/Documents/MATLAB/neuropixelPFC/Matlab/Utilities/Colormaps/cmap_VGlut2.mat');
%             elseif gg==3
%                 load('/Users/pierre/Documents/MATLAB/neuropixelPFC/Matlab/Utilities/Colormaps/cmap_NPY.mat');
%             elseif gg==4
%                 load('/Users/pierre/Documents/MATLAB/neuropixelPFC/Matlab/Utilities/Colormaps/cmap_Esr.mat');
%             end
%             
%             % Make heat map
%             %coronal
%             x=zeros(size(av,2),size(av,3));
%             for j = 1:size(av,3)
%                 for i = 1:size(av,2)
%                     temp_ij=find((round(ml2plot)==j) & (round(dv2plot)==i));
%                     if numel(temp_ij)==1
%                         x(i,j)=1;
%                     elseif numel(temp_ij)>1
%                         x(i,j)=1;
%                     end
%                 end
%             end
%             
%             %sagital
%             y=zeros(size(av,2),size(av,1));
%             for j = 1:size(av,1)
%                 for i = 1:size(av,2)
%                     temp_ij=find((round(ap2plot)==j) & (round(dv2plot)==i));
%                     if numel(temp_ij)==1
%                         y(i,j)=1;
%                     elseif numel(temp_ij)>1
%                         y(i,j)=1;
%                     end
%                 end
%             end
%             
%             
%             hm = imgaussfilt(x,win_sz);
%             figure;
%             hh2 = imagesc(hm);
%             
%             colormap(cmap)
%             c=colorbar('southoutside');
%             c.Label.String = 'activity densiy';
%             caxis([0, 0.005])
%             freezeColors
%             g=get(hh2);
%             imwrite (g.CData, ['/Users/pierre/OneDrive - KI.SE/Shared/LHA-LHb-PFC/Meeting6/cluster' num2str(clu2plot) '-hm' num2str(gg) '.png']);
%             close
%             
%             hm2 = imgaussfilt(y,win_sz);
%             figure;
%             hh4 = imagesc(hm2);
%             
%             colormap(cmap)
%             c=colorbar('southoutside');
%             c.Label.String = 'activity densiy';
%             caxis([0, 0.005])
%             freezeColors
%             g=get(hh4);
%             imwrite (g.CData, ['/Users/pierre/OneDrive - KI.SE/Shared/LHA-LHb-PFC/Meeting6/cluster' num2str(clu2plot) '-hm-sagital' num2str(gg) '.png']);
%             close
%             
%             
%             subplot(4,2,cnt)
%             atlas = imread('/Users/pierre/OneDrive - KI.SE/Shared/LHA-LHb-PFC/Meeting6/coronal.png');
%             heatmap = imread(['/Users/pierre/OneDrive - KI.SE/Shared/LHA-LHb-PFC/Meeting6/cluster' num2str(clu2plot) '-hm' num2str(gg) '.png']);
%             imagesc(imfuse(atlas,heatmap,'blend'));
%             colormap(cmap);
%             set(gca,'dataaspectratio',[1 1 1]);
%             set(gca,'XTick',[70 170 270 370 470 570 670 770 870 970 1070]);
%             set(gca,'XTickLabel',[-5 -4 -3 -2 -1 0 1 2 3 4 5]);
%             xlabel('ML (mm)')
%             set(gca,'YTick',[65 165 265 365 465 565 665 765]);
%             set(gca,'YTickLabel',[0 -1 -2 -3 -4 -5 -6 -7]);
%             ylabel('DV (mm)')
%             title(['Coronal Projection, cluster ' num2str(clu2plot)])
%             
%             subplot(4,2,cnt+1)
%             atlas2 = imread('/Users/pierre/OneDrive - KI.SE/Shared/LHA-LHb-PFC/Meeting6/sagital.png');
%             heatmap2 = imread(['/Users/pierre/OneDrive - KI.SE/Shared/LHA-LHb-PFC/Meeting6/cluster' num2str(clu2plot) '-hm-sagital' num2str(gg) '.png']);
%             imagesc(imfuse(atlas2,heatmap2,'blend'));
%             colormap(cmap);
%             set(gca,'dataaspectratio',[1 1 1]);
%             set(gca,'XTick',[40 140 240 340 440 540 640 740 840 940 1040 1140 1240 1340]);
%             set(gca,'XTickLabel',[-5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8]);
%             xlabel('AP (mm)')
%             set(gca,'YTick',[65 165 265 365 465 565 665 765]);
%             set(gca,'YTickLabel',[0 -1 -2 -3 -4 -5 -6 -7]);
%             ylabel('DV (mm)')
%             title(['Sagital Projection, cluster ' num2str(clu2plot)])
%             freezeColors
%             cnt = cnt + 2;
%         end
%         
%         set(gcf,'units','points','position',[0,800,1000,1600])
%         
%     end
%     
% end