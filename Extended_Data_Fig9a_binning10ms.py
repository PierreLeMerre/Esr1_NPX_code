# -*- coding: utf-8 -*-

from pynwb import NWBFile, NWBHDF5IO

import numpy as np
import pandas as pd
import os
import os.path
import sys
import scipy.io
import h5py
import re
# import seaborn as sns; sns.set()

 
# Import my custom functions
sys.path.append('C:\\Users\\admin\\OneDrive - KI.SE\\Skrivbordet\\General\\Projects\\Aversion\\neuropixelPFC\\Code_python\\')
from MS_suppl_functions import *

# data_nwb.units.electrodes.to_dataframe().reset_index()['location']


# Settings for spike binning  
binsize = 0.01
step = 0.01
window = [-0.1, 3]


# set up working folder
os.chdir('C:\\Users\\admin\\OneDrive - KI.SE\\Skrivbordet\\General\\Projects\\Aversion\\neuropixelPFC\\Code_python\\')
save_folder = 'C:\\Users\\admin\\OneDrive - KI.SE\\Skrivbordet\\General\\Projects\\Aversion\\neuropixelPFC\\Code_python\\PCA_virtual_mice_zscore\\'

session_folder = 'L:\\dmclab\\Pierre\\NPX_Database\\mPFC\\Aversion'
sess_list = os.listdir(session_folder)

# Upload the table with session description - take genotypes from there
sess_descrip = pd.read_csv(os.path.join("Aversion_session_description.tsv"),  sep="\t")

pfc_areas = ['ACAd', 'PL', 'ILA', 'ORBm'] # b'AI'  is out
all_genotypes = ['NPY-Cre', 'C57BL-6J', 'Vglut2-cre', 'Esr1-Cre-antero', 'Esr1-Cre-retro']

# Basic description of blocks
trials_blocks = [0, 1, 2, 3]
trials_labels = ['sound_pre', 'sound + opto', 'sound_post', 'sound + air-puff']
trials_colors = ['green', 'red', 'yellow', 'blue']

all_sessions_output = []
all_sess_genotype = []
all_units_sess = []
all_units_area = []
all_units_area_short = []
all_units_ids = []

for sess in sess_descrip['Name'][:]:
    print(sess)
    # if sess == '225759_20200521-probe0.nwb':
    #     continue
    
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
    
    all_names = [x.decode('UTF8') for x in list(data_nwb.units.electrodes.to_dataframe().reset_index()['location'])] # list of all units areas
    all_names_short = [re.findall("[a-zA-Z]+", x)[0] for x in all_names] # short areas without a specified layer
    
        
    # QC 1 - select units in PFC - ['ACAd', 'PL', 'IL', 'ORBm']
    if_pfc = np.array([(x in pfc_areas) for x in all_names_short]) # binary value

    # QC 2 - ISI violation - if does not violate ISI - put 1 (ISI should be less than 0.01)
    if_ISI = (ISIs / spike_counts) < 0.01
        
    # QC 3 - mean FR
    if_FR = mean_fr > 0.1 # binary value, FR > 0.1
        
    ## GETTING THE FINAL SUBSET OF UNITS 
    unit_QC = (np.sum(np.vstack((if_pfc, if_ISI, if_FR)), axis=0) == 3)
    unit_ids_QC = np.array(range(len(data_nwb.units)))[unit_QC] # 0 based indexing!!!
    print('Units in PFC: ', len(unit_ids_QC))
    np.save(save_folder + sess[:-11] + '_QC_boolean.npy', unit_QC)
    
    # Areas where only QC units belong - saved later
    unit_area_QC = np.array(all_names)[unit_QC]
    unit_area_short_QC = np.array(all_names_short)[unit_QC]
    
    
    sound_pre = data_nwb.trials[:]['start_time'].values[(data_nwb.trials[:]['Block'] == 1)]
    opto  = data_nwb.trials[:]['start_time'].values[(data_nwb.trials[:]['Block'] == 2)]
    sound_post = data_nwb.trials[:]['start_time'].values[(data_nwb.trials[:]['Block'] == 3)]
    airpuff = data_nwb.trials[:]['start_time'].values[(data_nwb.trials[:]['Block'] == 4)]


    trials_array = []
    for ttype, bl in zip([sound_pre, opto, sound_post, airpuff], trials_blocks):
        trials_array = trials_array + [str(bl)]*len(ttype)
    
    eventTimes = np.concatenate((sound_pre, opto, sound_post, airpuff), axis=0)


    # preallocate the matrix for binned spikes
    bin_num = len(np.arange(window[0], window[1]+step-binsize, step))
    binned_spikes = np.zeros((len(unit_ids_QC), bin_num, len(trials_array))) # n neurons x n bins x n trials
    
    # Time frames of the saturation epochs where spikes have to be removed (for the function RemoveSaturationSpikes)
    saturated_start = data_nwb.analysis['spike saturation start'].timestamps[:]
    saturated_stop = data_nwb.analysis['spike saturation stop'].timestamps[:]
    total_spikes_removed = 0
    # print(len(saturated_start), ' saturation periods')
    
    total_dur = 0
    for g in range(len(saturated_stop)):
        dur = saturated_stop[g] - saturated_start[g]
        total_dur = total_dur + dur
        # print(dur)
    # print('total', total_dur, ' s')
    
    
    for i in range(len(unit_ids_QC)):

        # after reset_index() units get id numbers starting from 0
        current_un = unit_ids_QC[i]
        #incl_area_list.append((data_nwb.units.electrodes.to_dataframe().reset_index()['location'][current_un]).decode('utf-8'))        
        
        if np.isnan(saturated_start[0]):
            spike_vec = data_nwb.units[current_un]['spike_times'].values[0]
        else: 
            spike_vec_full = data_nwb.units[current_un]['spike_times'].values[0]
            spike_vec, spikes_out = RemoveSaturationSpikes(saturated_start, saturated_stop, spike_vec_full)
            total_spikes_removed = total_spikes_removed  + spikes_out
            
        # binArray is n trials x n bins, i transpose it later for simplicity
        bin_num, time, binArray = SlidingWindow(spike_vec, eventTimes, binsize, step, window)
        binned_spikes[i, :, :] = binArray.T
        
    print(total_spikes_removed, ' spikes removed in the session due to saturation')
    print(np.unique(unit_area_QC), '\n')
    

    av_binned_sound_pre = np.mean((binned_spikes[:, :, np.in1d(trials_array, '0')]), axis = 2)
    av_binned_opto = np.mean((binned_spikes[:, :, np.in1d(trials_array, '1')]), axis = 2)
    av_binned_sound_post = np.mean((binned_spikes[:, :, np.in1d(trials_array, '2')]), axis = 2)
    av_binned_aiprpuff = np.mean((binned_spikes[:, :, np.in1d(trials_array, '3')]), axis = 2)
    
    av_concat = np.concatenate((av_binned_sound_pre, av_binned_opto, av_binned_sound_post, av_binned_aiprpuff), axis = 1)
    
    all_sessions_output.append(av_concat)
    all_sess_genotype.append(sess_genotype)
    all_units_sess.append([sess_genotype] * int(len(unit_ids_QC)))
    all_units_area.append(unit_area_QC)
    all_units_area_short.append(unit_area_short_QC)
    all_units_ids.append(unit_ids_QC)
    
    f.close()


# Remove sessions with low quality based on excel file    
all_units_sess_filt = np.array(all_units_sess)[np.array(sess_descrip['Analyze'].values[:], dtype=bool)]
all_units_area_filt = np.array(all_units_area)[np.array(sess_descrip['Analyze'].values[:], dtype=bool)]
all_units_area_short_filt = np.array(all_units_area_short)[np.array(sess_descrip['Analyze'].values[:], dtype=bool)]
all_units_ids_filt = np.array(all_units_ids)[np.array(sess_descrip['Analyze'].values[:], dtype=bool)]

all_sessions_output_filt = np.array(all_sessions_output)[np.array(sess_descrip['Analyze'].values[:], dtype=bool)]


# PCA for virtual mice - all PFC units
genotype_per_unit = np.concatenate(all_units_sess_filt)
area_per_unit = np.concatenate(all_units_area_filt)
area_per_unit_short = np.concatenate(all_units_area_short_filt)
id_per_unit = np.concatenate(all_units_ids_filt)

all_units_data = np.concatenate(all_sessions_output_filt)   

np.unique(genotype_per_unit,  return_counts=True)

# Saving data to upload it faster next time
np.save(save_folder + 'all_units_binned_10_QC.npy', all_units_data)
np.save(save_folder + 'area_per_unit_QC.npy', area_per_unit)
np.save(save_folder + 'area_per_unit_short_QC.npy', area_per_unit_short)
np.save(save_folder + 'genotype_per_unit_QC.npy', genotype_per_unit)
np.save(save_folder + 'id_per_unit_QC.npy', id_per_unit)




f.close()

# %% 
