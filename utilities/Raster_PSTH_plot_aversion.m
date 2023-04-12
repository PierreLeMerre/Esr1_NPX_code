function Raster_PSTH_plot_aversion(nwb_filepath,neuron_id,SR,nb_of_blocks)

unit = neuron_id; % number of the unit to plot
fp = nwb_filepath; % nwb from where to plot the unit
sr = SR; % sampling rate
b_nb = nb_of_blocks; % number of block to plot (integer)

% nwb file
nwb = nwbRead(fp);
% get id and recording duration from meta file
id = strsplit(nwb.general_session_id,'_');
% get PFC regions
ede_regions = nwb.general_extracellular_ephys_electrodes.vectordata.get('location').data.load();
% get unit main channel
main_ch = nwb.units.electrodes.data.load();

unit_region = (char(ede_regions{main_ch(unit)+1}));


% get control trial timestamps (one second before sound)
trial_ts = nwb.intervals_trials.start_time.data.load();
blocks = nwb.intervals_trials.vectordata.get('Block').data.load();


% Get unit spiketimes
% Load jagged arrays
unit_times_data = nwb.units.spike_times.data.load();
unit_times_idx = nwb.units.spike_times_index.data.load();
unit_ids = nwb.units.id.data.load(); % array of unit ids
% Initialize times Map containers indexed by unit_ids
unit_times = containers.Map('KeyType',class(unit_ids),'ValueType','any');
last_idx = 0;
unit_id = unit_ids(unit);
s_idx = last_idx + 1;
e_idx = unit_times_idx(unit);
unit_times(unit_id) = unit_times_data(s_idx:e_idx);

% Parameters
spk_ts = unit_times(unit-1).*sr; % -1 indexed in unit_times!!!
bin_sz = 0.01; % in seconds
pre = -1; % in seconds
post = 3; % in seconds

figure
for b = 1 : b_nb
trials1 = trial_ts(blocks==b).*sr; % to plot
% Raster plot
ax1 = subplot(2,b_nb,b);
spikeTimes = cell(length(trials1),1);
for Trial=1:length(trials1)
    timestamps=(spk_ts(spk_ts>trials1(Trial)+pre*sr & spk_ts<=trials1(Trial)+post*sr)'-trials1(Trial));
    if isempty(timestamps)==0
    spikeTimes{Trial} = timestamps./sr;
    end
end
plotSpikeRaster(spikeTimes,'PlotType','scatter');
line([0 0],[0 length(trials1)],'LineWidth',1,'Color',[0 0 0])
set(gca,'XTickLabel',[]);
set(gca,'XLabel',[]);
ylabel('Trial number')
title(['Block' num2str(b)])


% PSTH
ax2 = subplot(2,b_nb,b+b_nb);
% WindowCenters is the center of the bins (xaxis)
WindowStep=bin_sz*sr;
NumberOfTrials = numel(trials1);
RelativeEdges = pre*sr:WindowStep:post*sr;
NumberOfEdges = numel(RelativeEdges);
NumberOfWindows = NumberOfEdges - 1;
% bin the spikes in windows to have one column per trial
SpikeCounts=zeros(NumberOfWindows,NumberOfTrials);
for AnchorTimeIndex = 1:NumberOfTrials
    BinEdges = trials1(AnchorTimeIndex) + RelativeEdges;
    [SpikeCounts(:,AnchorTimeIndex),~]=histcounts(spk_ts, BinEdges);
end
SpikeRates = SpikeCounts / bin_sz;
WindowCenters = (RelativeEdges(1:NumberOfWindows) + bin_sz/2) / sr;

psth = nanmean(SpikeRates,2);
plot(WindowCenters,psth,'k')
line([0 0],[0 max(psth)+5],'LineWidth',1,'Color',[0 0 0])
xlabel('Time (s)')
ylim([0 max(psth)+5]);
ylabel('firing rate (Hz)')
linkaxes([ax1 ax2],'x')


end
suptitle(['Mouse ' id{1} ', session ' id{2} ', unit ' num2str(unit) ' (' unit_region ')'])

end