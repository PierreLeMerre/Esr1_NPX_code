# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 14:56:56 2021

@author: Marina Slashcheva
"""
from pynwb import NWBFile, NWBHDF5IO

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import os.path
import sys
import scipy.io
import h5py



#%%
def SlidingWindow(spikeTimes,eventTimes, binsize, step, window):
    ''' 
    9.07.2020 MS
    
    supplementary function called within 'getOverlappingSpikeCounts'
    bin_num_NI, time_NI, binArray_NI = SlidingWindow(spikeTimes, nat_im_events, psthBinSize, step, window_NI)
    spikes, times, trials = getOverlappingSpikeCounts(Sess, binsize=bin_size, step=step)
        
    '''    
    
    if len(eventTimes) == 0:
        return 0, np.empty(0), np.empty(0)
    
    bin_num = len(np.arange(window[0], window[1]+step-binsize, step))
    binArray = np.zeros((len(eventTimes), bin_num))
    
    timeVec_tmp = np.zeros((len(eventTimes), bin_num))
    for ev in range(len(eventTimes)):
        start = window[0] + eventTimes[ev]
        end = window[1] + eventTimes[ev]
    
        spikes = spikeTimes[(spikeTimes > start) & (spikeTimes < end)]
        
        bin_starts = np.arange(start, end+step-binsize, step) #####
        bin_ends = bin_starts + binsize
        last_index = np.searchsorted(spikes, bin_ends, side='right')
        first_index = np.searchsorted(spikes, bin_starts, side='left')
    
        binArray[ev, :] = (last_index - first_index)[:bin_num]
        timeVec_tmp[ev, :] = bin_starts[:bin_num] + 1/2*step
    
    timeVec_tmp = np.reshape(timeVec_tmp, (1, np.size(timeVec_tmp)))[0]
    return  bin_num, timeVec_tmp, binArray


def SpikeCount_and_PSTH_per_unit(spikeTimes, eventTimes, psthBinSize, window):
    
    ''' 
    9.07.2020 MS
    Function to compute spike counts and PSTH for 1 unit over all trials given
    
    INPUTS
    spikeTimes - vector of spike times
    eventTimes - vector with start time for all stimuli to compute
    psthBinSize - size of the bin in sec
    window - [-0.2, 0.8] - time window around eventTimes
    
    OUTPUTS
    psth - average spiking of the unit across trials in HZ
    spikeCounts - total number of spikes in each bin across all trials
    binborders - specific borders
    binArray - array of size nEvents x nBins with spike counts
    
    USE
    psth, binBorders, spikeCounts, binArray = psth_per_unit(spikeTimes, eventTimes, psthBinSize, window)
    
    '''
    if len(eventTimes) == 0:
        return np.empty(0), np.empty(0), np.empty(0), np.empty(0)
    
    spikeTimes = spikeTimes[(spikeTimes > np.min(eventTimes)+window[0]) & (spikeTimes < np.max(eventTimes) + window[1])]
    
    binBorders = np.linspace(window[0], window[1], round((window[1] - window[0])/psthBinSize)+1)
    numBins = len(binBorders) - 1   

    binArray = np.zeros((len(eventTimes), numBins))

    for ev in range(len(eventTimes)):
        counts, binBordRealtime = np.histogram(spikeTimes, bins = binBorders+eventTimes[ev])
        binArray[ev, :] = counts
        
    spikeCounts = np.sum(binArray, 0) # total spike counts over all trials
    psth = np.mean(binArray/psthBinSize, 0) # same but in Hz
    
    #plt.plot(binBorders[:-1], psth)
    return psth, binBorders, spikeCounts, binArray


def SpikeCount_per_unit(spikeTimes, eventTimes, psthBinSize, window):
    
    ''' 
    9.07.2020 MS
    Function to compute spike counts for 1 unit over all trials given
    
    INPUTS
    spikeTimes - vector of spike times
    eventTimes - vector with start time for all stimuli to compute
    psthBinSize - size of the bin in sec
    window - [-0.2, 0.8] - time window around eventTimes
    
    OUTPUTS
    spikeCounts - total number of spikes in each bin across all trials
    binborders - specific borders
    binArray - array of size nEvents x nBins with spike counts
    
    USE
    binBorders, spikeCounts, binArray = SpikeCount_per_unit(spikeTimes, eventTimes, psthBinSize, window)
    
    '''
    if len(eventTimes) == 0:
        return np.empty(0), np.empty(0), np.empty(0), np.empty(0)
    
    spikeTimes = spikeTimes[(spikeTimes > np.min(eventTimes)+window[0]) & (spikeTimes < np.max(eventTimes) + window[1])]
    
    binBorders = np.linspace(window[0], window[1], round((window[1] - window[0])/psthBinSize)+1)
    numBins = len(binBorders) - 1   

    binArray = np.zeros((len(eventTimes), numBins))

    for ev in range(len(eventTimes)):
        counts, binBordRealtime = np.histogram(spikeTimes, bins = binBorders+eventTimes[ev])
        binArray[ev, :] = counts
        
    #spikeCounts = np.sum(binArray, 0) # total spike counts over all trials
    return binBorders, binArray

# %%

def RemoveSaturationSpikes(start, stop, spike_vec):
    '''
    Parameters
    ----------
    start : numpy.ndarray, start of saturation epoch
    stop : numpy.ndarray, start of saturation epoch
    spike_vec : numpy.ndarray, spikes of single unit

    Returns: clean spike vector, number of removed spikes
    -------
    TYPE
    Remove spikes of saturation epoch (+- 1 sec around saturation)
    '''
    ind_to_remove = []
    
    for sat in range(len(start)):
        if spike_vec[0] > stop[sat] or spike_vec[-1] < start[sat]:
            continue
        #print(start[sat], stop[sat])
        sat_start_ind = np.where(spike_vec > start[sat]-1)[0][0]
        sat_stop_ind = np.where(spike_vec < stop[sat]+1)[0][-1] 
        #print(sat_start_ind, sat_stop_ind)
        
        # Periods of saturation where sat_start_ind > sat_stop_ind do not have spikes
        # Periods of saturation where sat_start_ind = sat_stop_ind have one spike
        if sat_start_ind == sat_stop_ind:
            ind_to_remove.append([sat_start_ind])
        if sat_start_ind < sat_stop_ind:
            ind_to_remove.append(list(range(sat_start_ind, sat_stop_ind+1)))
    
    if len(ind_to_remove) != 0:
        unique_to_remove = np.unique(np.concatenate(ind_to_remove))
        spike_vec_clean = np.delete(spike_vec, unique_to_remove)
        # print(len(unique_to_remove))
        return spike_vec_clean, len(unique_to_remove)
    if len(ind_to_remove) == 0:
        return spike_vec, 0














