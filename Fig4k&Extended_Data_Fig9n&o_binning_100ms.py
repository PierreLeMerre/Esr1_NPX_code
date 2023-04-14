# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 16:59:53 2021

@author: admin
Copy of 'Binning_ITI' - final version
"""
from pynwb import NWBFile, NWBHDF5IO

import numpy as np
import pandas as pd
import os
import os.path
import sys
import scipy.io
import h5py
import re

# Import my custom functions
sys.path.append('C:\\Users\\admin\\OneDrive - KI.SE\\Skrivbordet\\General\\Projects\\Aversion\\neuropixelPFC\\Code_python\\')
from MS_suppl_functions import *

#%%
binsize = 0.1
window = [-2.1, -0.1]

# set up working folder
os.chdir('C:\\Users\\admin\\OneDrive - KI.SE\\Skrivbordet\\General\\Projects\\Aversion\\neuropixelPFC\\Code_python\\')
save_folder = 'C:\\Users\\admin\\OneDrive - KI.SE\\Skrivbordet\\General\\Projects\\Aversion\\neuropixelPFC\\Code_python\\Decoding\\'
session_folder = 'L:\\dmclab\\Pierre\\NPX_Database\\mPFC\\Aversion'

# Upload the table with session description - take genotypes from there
sess_descrip = pd.read_csv(os.path.join("Aversion_session_description.tsv"),  sep="\t")

all_genotypes = ['NPY-Cre', 'C57BL-6J', 'Vglut2-cre', 'Esr1-Cre-antero', 'Esr1-Cre-retro']
pfc_areas = ['ACAd', 'PL', 'ILA', 'ORBm'] # b'AI'  is out

all_sessions_output = []
all_sess_genotype = []
all_units_gt = []
all_units_sess = []
all_units_area = []
all_units_area_short = []
all_units_ids = []

for sess in sess_descrip['Name'][:]:
    print(sess)
    
    f = NWBHDF5IO((session_folder + '\\' + sess), 'r')
    data_nwb = f.read()
    sess_genotype = sess_descrip['Genotype'][sess_descrip['Name'] == sess].values[0]
    print(sess_genotype)   
    
    # I have to open the session as hdf file because pynwb does not see all fields for some reason (units/quality)
    hfile = h5py.File(session_folder + '\\' + sess)
    ISIs = hfile['units/quality'][()][:, 0]
    spike_counts = np.diff(hfile['units/spike_times_index'][()], prepend=0)
    hfile.close()
    
    # Upload computed mean FR for all session units
    # Github path to precalculated FRs: https://github.com/PierreLeMerre/Esr1_NPX_code/tree/main/analysis/Mean_FR

    mean_FR_path = 'C:\\Users\\admin\\OneDrive - KI.SE\\Skrivbordet\\General\\Projects\\Aversion\\neuropixelPFC\\Matlab\\analysis\\Mean_FR\\'
    mean_FR_fullpath = mean_FR_path + sess + '_fr' + '.mat'
    ff = scipy.io.loadmat(mean_FR_fullpath)
    mean_fr = ff['firing_rates'][0]
    
    ## QUALITY CONTROL IN 3 STEPS
    
    # QC 1 - select units in PFC - ['ACAd', 'PL', 'IL', 'ORBm']
    all_names = [x.decode('UTF8') for x in list(data_nwb.units.electrodes.to_dataframe().reset_index()['location'])] # list of all units areas
    all_names_short = [re.findall("[a-zA-Z]+", x)[0] for x in all_names] # short areas without a specified layer
    if_pfc = np.array([(x in pfc_areas) for x in all_names_short]) # binary value

    # QC 2 - ISI violation - if does not violate ISI - put 1 (ISI should be less than 0.01)
    if_ISI = (ISIs / spike_counts) < 0.01
        
    # QC 3 - mean FR
    if_FR = mean_fr > 0.1 # binary value, FR > 0.1
        
    ## GETTING THE FINAL SUBSET OF UNITS 
    unit_QC = (np.sum(np.vstack((if_pfc, if_ISI, if_FR)), axis=0) == 3)
    unit_ids_QC = np.array(range(len(data_nwb.units)))[unit_QC] # 0 based indexing!!!
    print('Units in PFC: ', len(unit_ids_QC))
    # np.save(save_folder + sess[:-11] + '_QC_boolean.npy', unit_QC)
    
    # Areas where only QC units belong - saved later
    unit_area_QC = np.array(all_names)[unit_QC]
    unit_area_short_QC = np.array(all_names_short)[unit_QC]
    

    # Timestamps of the beginning of the trials in each of 4 blocks
    sound_pre = data_nwb.trials[:]['start_time'].values[(data_nwb.trials[:]['Block'] == 1)]
    opto  = data_nwb.trials[:]['start_time'].values[(data_nwb.trials[:]['Block'] == 2)]
    sound_post = data_nwb.trials[:]['start_time'].values[(data_nwb.trials[:]['Block'] == 3)]
    airpuff = data_nwb.trials[:]['start_time'].values[(data_nwb.trials[:]['Block'] == 4)]
    
    # Take 49 trials in each block, for opto block with 100 trials take every 2nd to sample uniformly 
    eventTimes = np.concatenate((sound_pre[:49], opto[::2][:49], sound_post[:49], airpuff[:49]), axis=0)
    
    # preallocate the matrix for binned spikes
    bin_num = len(np.arange(window[0], window[1], binsize))
    binned_spikes = np.zeros((len(unit_ids_QC), bin_num * len(eventTimes))) # n neurons x n bins x n trials
    
    
    # Time frames of the saturation epochs where spikes have to be removed (for the function RemoveSaturationSpikes)
    saturated_start = data_nwb.analysis['spike saturation start'].timestamps[:]
    saturated_stop = data_nwb.analysis['spike saturation stop'].timestamps[:]
    total_spikes_removed = 0
    print(len(saturated_start), ' saturation periods')
    
     
    for current_un in range(len(unit_ids_QC)):        

        if np.isnan(saturated_start[0]):
            spike_vec = data_nwb.units[current_un]['spike_times'].values[0]
        else: 
            spike_vec_full = data_nwb.units[current_un]['spike_times'].values[0]
            spike_vec, spikes_out = RemoveSaturationSpikes(saturated_start, saturated_stop, spike_vec_full)
            total_spikes_removed = total_spikes_removed  + spikes_out
            
        # binArray is n trials x n bins, i transpose it later for simplicity
        binBorders, binArray = SpikeCount_per_unit(spike_vec, eventTimes, binsize, window) # sound_pre[:49]  opto[::2][:49] airpuff[:49]
                
        binned_spikes[current_un, :] = binArray.flatten()
    
    print(total_spikes_removed, ' spikes removed in the session due to saturation')
    print(np.unique(unit_area_QC), '\n')

    all_sess_genotype.append(sess_genotype)
    all_sessions_output.append(binned_spikes)
    
    all_units_gt.append([sess_genotype] * len(unit_ids_QC))
    all_units_sess.append([sess] * len(unit_ids_QC))
    all_units_area.append(unit_area_QC)
    all_units_area_short.append(unit_area_short_QC)

    f.close()

# Remove sessions with low quality based on excel file    
all_units_gt_filt = np.array(all_units_gt)[np.array(sess_descrip['Analyze'].values[:], dtype=bool)]
all_units_area_filt = np.array(all_units_area)[np.array(sess_descrip['Analyze'].values[:], dtype=bool)]
all_units_area_short_filt = np.array(all_units_area_short)[np.array(sess_descrip['Analyze'].values[:], dtype=bool)]
all_units_sess_filt = np.array(all_units_sess)[np.array(sess_descrip['Analyze'].values[:], dtype=bool)]
all_sessions_output_filt = np.array(all_sessions_output)[np.array(sess_descrip['Analyze'].values[:], dtype=bool)]


# PCA for virtual mice - all PFC units
genotype_per_unit = np.concatenate(all_units_gt_filt)
area_per_unit = np.concatenate(all_units_area_filt)
sess_per_unit = np.concatenate(all_units_sess_filt)
area_per_unit_short = np.concatenate(all_units_area_short_filt)

all_units_data = np.concatenate(all_sessions_output_filt)   

# Saving data to upload it faster
# np.save(save_folder + 'all_units_ITI_100ms_QC_all_blocks.npy', all_units_data)

# np.save(save_folder + 'area_per_unit_ITI_QC.npy', area_per_unit)
# np.save(save_folder + 'area_per_unit_short_ITI_QC.npy', area_per_unit_short)
# np.save(save_folder + 'genotype_per_unit_ITI_QC.npy', genotype_per_unit)
# np.save(save_folder + 'sess_per_unit_ITI_QC.npy', sess_per_unit)



