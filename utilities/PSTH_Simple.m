%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to compute a simple PSTH
%
% Written by Vahid Esmaeili and adapted by Pierre Le Merre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [SpikeRates,WindowCenters]= PSTH_Simple(SpikeTimes, AnchorTimes, PreTime, PostTime,sr, WindowSize)

% SpikeTimes in samples
% Anchortimes in samples
% PreTime in seconds
% PostTime in seconds
% sr is the smpling rte
% WindowSize in seconds (usually people use 10-50ms)

% WindowCenters is the center of the bins (xaxis)

WindowStep=WindowSize*sr;
NumberOfTrials = numel(AnchorTimes);
RelativeEdges = PreTime*sr:WindowStep:PostTime*sr;
NumberOfEdges = numel(RelativeEdges);
NumberOfWindows = NumberOfEdges - 1;

% bin the spikes in windows to have one column per trial
SpikeCounts=zeros(NumberOfWindows,NumberOfTrials);
for AnchorTimeIndex = 1:NumberOfTrials
    BinEdges = AnchorTimes(AnchorTimeIndex) + RelativeEdges;
    [SpikeCounts(:,AnchorTimeIndex),~]=histcounts(SpikeTimes, BinEdges);
end

SpikeRates = SpikeCounts / WindowSize;
WindowCenters = (RelativeEdges(1:NumberOfWindows) + WindowSize/2 );

