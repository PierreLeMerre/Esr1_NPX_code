%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script that plot the PSTH and the raster of the selected units.
%
% Written by Pierre Le Merre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cls
format long
set(0, 'DefaultFigureRenderer', 'painters')

%% Load
% Task
Task = 'Aversion'; %'Passive','Context','Aversion', 'Opto'

% Allen Brain Atlas regions to count the units
Regions = {'ACAd','PL','ILA','ORBm'};
%Regions = {'MOs','ACAd','ACAv','PL','ILA','ORBm','ORBvl','ORBl'};

% Regions = {'ACAd2/3','ACAd5','ACAd6a','ACAd6b',...
%     'PL2/3','PL5','PL6a','PL6b',...
%     'ILA2/3','ILA5','ILA6a','ILA6b',...
%     'ORBm2/3','ORBm5','ORBm6a','ORBm6b'};

% Your path to the Database
%Path2Data = '/Volumes/dmclab/pielem/NPX_Database/';
Path2Data = '/Volumes/T7/NPX_Database/mPFC/Aversion/';
D = dir([Path2Data '*.nwb']);

% Your path to analysis folder
BasePath = '/Users/pielem/Documents/MATLAB/';

% Your path to analysis folder
Path2Ana = '/Users/pielem/Documents/MATLAB/neuropixelPFC/Matlab/analysis/';

% Subset to plot
task_mod = 0;
plotsubset = 0; % bolean to plot a subset
absolute = 0;
relative = 0; % relative change dataset
clustering = 0; % plot from clustering folder
conditioning = 0; % trace conditioning dataset
glm = 0;
subset = 'tmu_OR'; %'tmu_opto_neg' or just a cluster number;
Block_nb = 2; % Block which the subset is from
Block_nb2 = 4; % Combine with another block

% Min ISI violation rate
MIN_ISI = 0.01; % less than 1% of ISI violations

% Min FR
MIN_FR = 0.1; % mean session FR above 0.1Hz

% number of block to plot
Block_plot = 4;

% Trial per condition
Trial_per_condition = 50;

% Time window to plot
window = [-1 3];

% Psth bin size
psthBinSize = 0.01;

% Gaussian Kernel 50ms
N = 0.05/psthBinSize;
g = gausswin(N);
w=g./trapz(g); % normalize area = 1

% files to plot
Genotypes ={'NPY','VGlut2','WT','Esr'};

%for G = 1:numel(Genotypes) % GENOTYPE LOOP

Genotype = Genotypes{1};

if strcmp(Genotype,'NPY')
    f_start = 1;
    f_stop = 5;
elseif strcmp(Genotype,'VGlut2')
    f_start = 6;
    f_stop = 10;
elseif strcmp(Genotype,'WT')
    f_start = 11; t 
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


%initialize global vars
Ede_regions = [];
Main_ch = [];
ISregion_mtrx = [];
UNIT = [];

for f = f_start : f_stop % Loop through mice


    % Read nwb file
    nwb = nwbRead([Path2Data D(f).name]);
    % get id and recording duration from meta file
    id = strsplit(nwb.general_session_id,'_');

    id2 = strsplit(id{2},'-');

    disp(['Getting units from mouse ' id{1} ', session ' id{2} '...'])

    % Get PFC regions
    ede_regions = nwb.general_extracellular_ephys_electrodes.vectordata.get('location').data.load();
    Ede_regions = [Ede_regions ede_regions];
    % Get unit main channel
    main_ch = nwb.units.electrodes.data.load();
    Main_ch = [Main_ch; main_ch];

    isregion_mtrx = zeros(numel(main_ch),numel(Regions));
    % loop over regions
    for R = 1:numel(Regions)
        region=char(Regions{R});
        isregion = zeros(1,numel(main_ch));

        % loop over units
        for i = 1 : numel(main_ch)

            % Make logical region vector
            if (strfind(ede_regions(main_ch(i)+1,:),region)) > 0 % main_ch is 0 indexed!
                isregion(i) = 1;
            end

        end
        isregion_mtrx(:,R) = isregion;
    end
    ISregion_mtrx = [ISregion_mtrx; isregion_mtrx];

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

    % Get PSTH Events
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

    % Get mean session FR
    load([BasePath 'Esr1_NPX_code/analysis/Mean_FR/' D(f).name '_fr.mat'])
    mean_session_fr = firing_rates;

    % Get quality metrics
    spk_number = nan(1,length(unit_ids));
    spk_number(1) = unit_times_idx(1);
    for u = 2:length(unit_ids)
        spk_number(u) = unit_times_idx(u)-unit_times_idx(u-1);
    end
    qual = nwb.units.vectordata.get('quality').data.load;
    isi = qual(1,:);
    isi_perc = isi./spk_number;

    %preallocate
    PSTH = nan(Block_plot,numel(Regions),numel(unit_ids),diff(window/psthBinSize)+1);
    RASTERX = cell(Block_plot,numel(Regions));
    RASTERY = cell(Block_plot,numel(Regions));


    % loop over blocks
    for b = 1: Block_plot


        eval(['bloc' num2str(b) '_ts = trial_ts(blocks==' num2str(b) ');'])


        if b==2 && strcmp(Task,'Aversion')
            bloc2_ts = bloc2_ts(1:2:end); % one trial over 2
        end
        eval(['bb = bloc' num2str(b) '_ts;'])

        if numel(bb)>Trial_per_condition % make sure to take 50 trials
            eval(['bloc' num2str(b) '_ts = bloc' num2str(b) '_ts(1:' num2str(Trial_per_condition) ');'])
        end
        eval(['bb = bloc' num2str(b) '_ts;'])

        % get saturation epochs
        sat_start = nwb.analysis.get('spike saturation start').timestamps.load();
        sat_stop = nwb.analysis.get('spike saturation stop').timestamps.load();

        % Loop through regions
        for R = 1:numel(Regions)

            cnt = 1;

            % Loop through units
            for i = 1 : numel(unit_ids)

                if isregion_mtrx(i,R) == 1 && isi_perc(i)<MIN_ISI && mean_session_fr(i)>MIN_FR

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
                    [psth, bins, rasterX, rasterY, spikeCounts] = psthRasterAndCounts(vecSpikeTimes, vecEventTimes, window, psthBinSize);
                    gk=conv(squeeze(psth),w,'same'); % Gaussian smoothing
                    %
                    PSTH(b,R,cnt,:) = gk;
                    RASTERX{b,R}{1,cnt} = rasterX;
                    RASTERY{b,R}{1,cnt} = rasterY;
                    UNIT(cnt) = unit_ids(i);
                    %increase counter
                    cnt = cnt +1;

                end

            end % End of loop through units

        end % End of loop through regions

    end % End loop over blocks




    %% Plot

    figure;
    wvs = nwb.units.waveform_mean.data.load();

    for R = 1:numel(Regions) % Loop over regions
        disp(['Plotting region ' Regions{R}])
        cnt = 1;

        % Loop through units
        for i = 1 : numel(unit_ids)

            unit = i;

            if isregion_mtrx(unit,R) == 1 && isi_perc(i)<MIN_ISI && mean_session_fr(i)>MIN_FR &&  cnt<=size(RASTERX{b,R},2)
                
                % plot waveform
                subplot(2,5,1)
                plot(wvs(:,unit),'k')

                % loop over blocks
                for b = 1:Block_plot
                    subplot(2,5,b+1)
                    plot(bins,squeeze(PSTH(b,R,cnt,:)),'k')
                    xlim(window)
                    line([0 0],[0 max(max(squeeze(PSTH(:,R,cnt,:))))+5],'LineWidth',1,'Color',[0 0 0])
                    ylabel('Firing rate (Hz)')
                    title(['Block ' num2str(b)])

                    subplot(2,5,4+b+2)
                    if ~isempty(RASTERX{b,R}{1,cnt})
                        plot(RASTERX{b,R}{1,cnt},RASTERY{b,R}{1,cnt},'k')
                        xlim(window)
                        ylim([0 Trial_per_condition])
                        line([0 0],[0 Trial_per_condition],'LineWidth',1,'Color',[0 0 0])
                        if b==2 || b==4
                            line([0.7 0.7],[0 Trial_per_condition],'LineWidth',1,'Color',[0 0 0])
                            xlabel('Time (s)')
                            ylabel('Trials')
                        end
                    end
                end

                % title and figure
                sgtitle(['Genotype: ' Genotype ', Mouse: ' id{1} ', Session: ' id{2} ', Unit: ' num2str(unit) ', Region: ' Regions{R}])
                set(gcf,'units','points','position',[0,900,1000,600])
                saveas(gcf,['/Users/pielem/Library/CloudStorage/OneDrive-KI.SE/Pierre_Shared/LHA-LHb-PFC/Nature_Neuroscience_Revisions/All_Units/' Genotype '_' id{1} '_' id{2} '_Unit_' num2str(unit) '_' Regions{R} '.png'])

                close()

                %increase counter
                cnt = cnt +1;

            end % end of if loop

        end % End of the loop over units

    end % End of loop over regions


end % End of Loop through mice

%end % END OF GENOTYPE LOOP