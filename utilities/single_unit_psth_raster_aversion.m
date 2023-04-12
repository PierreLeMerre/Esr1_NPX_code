function  [psth, bins, rasterX, rasterY, spikeCounts,ede_region,main_ch] = ...
    single_unit_psth_raster_aversion(Path, Mouse, Session, unit, Block2plot, Trials2plot, timewindow, BinSize)

% Read nwb file
nwb = nwbRead([Path Mouse '_' Session '.nwb']);
% Get PFC regions
ede_regions = nwb.general_extracellular_ephys_electrodes.vectordata.get('location').data.load();

% Get unit main channel
main_ch = nwb.units.electrodes.data.load(unit);
ede_region = ede_regions(main_ch);

% Get unit spiketimes, Load jagged arrays
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

vecSpikeTimes = unit_times(unit-1); % -1 indexed in unit_times!!!


% Get PSTH Events
trial_ts = nwb.intervals_trials.start_time.data.load();
blocks = nwb.intervals_trials.vectordata.get('Block').data.load();
trial_ts = trial_ts(blocks == Block2plot);
trial_ts = trial_ts(1:Trials2plot);

% get saturation epochs
sat_start = nwb.analysis.get('spike saturation start').timestamps.load();
sat_stop = nwb.analysis.get('spike saturation stop').timestamps.load();

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
vecEventTimes = trial_ts;
idx2 = find(0<vecEventTimes-s1 &  vecEventTimes-s2<0);
if ~isempty(idx2)
    vecEventTimes(idx2) = [];
end

% Gaussian Kernel 50ms
%N = 0.05/psthBinSize;
N = 10;
g = gausswin(N);
w=g./trapz(g); % normalize area = 1
%preallocate
PSTH = nan(Trials2plot/10,1,diff(timewindow/BinSize)+1); %,numel(Regions)
RASTERX = cell(Trials2plot/10,1);
RASTERY = cell(Trials2plot/10,1);
pp =0;
for p = 1 : 10 : Trials2plot
    pp = pp + 1;
    vec =  vecEventTimes(p:p+10-1);
    
    % calculate simple PSTH
    [psth, bins, rasterX, rasterY, spikeCounts] = psthRasterAndCounts(vecSpikeTimes, vec, timewindow, BinSize);
     gk=conv(squeeze(psth),w,'same'); % Gaussian smoothing
    PSTH(pp,:) = gk ;%- mean(gk(1:100)); % baseline correction
    RASTERX{pp,1}{1,1} = rasterX;
    RASTERY{pp,1}{1,1} = rasterY;
end


[cmap] = cbrewer('seq', 'Reds',Trials2plot/10);
figure;
subplot(2,1,1)
pp=0;
for p = 1 : 10 : Trials2plot
    pp = pp + 1;
    hold on
    plot(bins,PSTH(pp,:),'Color',cmap(pp,:))
end
xlim(timewindow)
ylim([0 10]) %max(psth)+2
line([0 0],[0 max(psth)+2],'LineWidth',1,'Color',[0 0 0])
if Block2plot==2 || Block2plot==4
    line([0.7 0.7],[0 max(psth)+2],'LineWidth',1,'Color',[0 0 0])
    xlabel('Time (s)')
    ylabel('Trials')
end


subplot(2,1,2)
inc = 0;
pp=0;
for p = 1 : 10 : Trials2plot
    pp = pp + 1;
    hold on
plot(RASTERX{pp,1}{1,1},RASTERY{pp,1}{1,1}+inc,'Color',cmap(pp,:))
inc = inc +10;
end
xlim(timewindow)
ylim([0 Trials2plot])
line([0 0],[0 Trials2plot],'LineWidth',1,'Color',[0 0 0])
if Block2plot==2 || Block2plot==4
    line([0.7 0.7],[0 Trials2plot],'LineWidth',1,'Color',[0 0 0])
    xlabel('Time (s)')
    ylabel('Trials')
end

suptitle(['Mouse: ' Mouse ', Session: ' Session ', Unit: ' num2str(unit_ids(unit)) ', Region: ' ede_region{1}])

end