# %% This script requires output from 'Fig5k&Extended_Data_Fig9n&o_binning_100ms.py'

from pynwb import NWBFile, NWBHDF5IO

import numpy as np
import matplotlib.pyplot as plt
import os
import os.path
import sys
from scipy import stats as st # z-scoring

# Import my custom functions
sys.path.append('C:\\Users\\admin\\OneDrive - KI.SE\\Skrivbordet\\General\\Projects\\Aversion\\neuropixelPFC\\Code_python\\')
from MS_suppl_functions import *

# K-Fold cross-validation
from sklearn.model_selection import RepeatedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, roc_auc_score


#%% 
# Indicate to blocks to decode against each other, for example [1,2] to decode state in block 1 (sounds) against block 2 (sound + opto)
to_decode = [1,2]

# Loading data - for all units that passed QC!
# ITI of 49 trials in block1 or block2 - [-2.1, -0.1], 100ms bins

# set up working folder
os.chdir('C:\\Users\\admin\\OneDrive - KI.SE\\Skrivbordet\\General\\Projects\\Aversion\\neuropixelPFC\\Code_python\\')
save_folder = 'C:\\Users\\admin\\OneDrive - KI.SE\\Skrivbordet\\General\\Projects\\Aversion\\neuropixelPFC\\Code_python\\Decoding\\'

data_all_blocks = np.load(save_folder + 'all_units_ITI_100ms_QC_all_blocks.npy')

# Choosing the trials of selected blocks in 'to_decode'
block_subset = np.zeros((data_all_blocks.shape[1], ), dtype=bool)
points_in_block = data_all_blocks.shape[1] // 4
for b in to_decode:
    print(b)
    [start, stop] = [points_in_block * (b -1),  points_in_block * b]
    print(start, stop)
    block_subset[start:stop] = True
    

# Z-score 2 blocks together
data_conc_norm = st.zscore(data_all_blocks, axis=1)
# delete deurons with no spikes (NaNs)
nan_values = ~np.isnan(np.sum(data_conc_norm, axis=1))
data_conc_norm = data_conc_norm[nan_values, :]

# Load information to identify each neuron's genotype, session etc
area_all = np.load(save_folder + 'area_per_unit_short_ITI_QC.npy')[nan_values]
sess_all = np.load(save_folder + 'sess_per_unit_ITI_QC.npy')[nan_values]
genotype_all= np.load(save_folder + 'genotype_per_unit_ITI_QC.npy')[nan_values]
genotype_all= np.array([i[:3] for i in genotype_all])

# Apply the array with block selection to only take 2 blocks for decoding
data_conc_norm = data_conc_norm[:, block_subset]

# Assign binary values to the blocks, first block will be 0s, second = 1s
y = np.zeros((np.shape(data_conc_norm)[1]))
y[points_in_block:] = 1 # change to points_in_block


# %% This is to calculate and plot for different numbers of cells - per genotype

# Cell numbers to decode 
num_points = np.arange(5, 155, 5)
metric = 'Accuracy' # or 'AUC'

# This is to split data always into 50-50%
rcv = RepeatedKFold(n_splits=2, n_repeats=25, random_state=1)
classifier = LogisticRegression(solver='lbfgs', max_iter = 1000) #, penalty='l2', C=0.001)


all_gt = []
all_ar = []
all_mean_val = []
all_sd_val = []

for gt in range(4):
    print('Decoding for ', np.unique(genotype_all)[gt])
    gt_ind = np.where(genotype_all == np.unique(genotype_all)[gt])[0]
    
    # Loop through areas
    gt_areas_all = np.unique(area_all[gt_ind], return_counts=True)
    gt_areas = gt_areas_all[0][np.where(gt_areas_all[1] > 20)]
    gt_areas_num = gt_areas_all[1][np.where(gt_areas_all[1] > 20)]
    
    auc_per_area = np.zeros((len(gt_areas), len(num_points), 50))
    #auc_per_area_shuff = np.zeros((len(gt_areas), 50))

    for ar in range(len(gt_areas)):
        
        # To accumulate all data 
        all_gt.append(np.unique(genotype_all)[gt])
        all_ar.append(gt_areas[ar])

        area_all_ind = np.where(area_all == gt_areas[ar])[0]
        area_gt_ind = np.intersect1d(gt_ind, area_all_ind)
        print(gt_areas[ar],', ', gt_areas_num[ar], ' neurons')

        # Now iterate over different number of units for predictions: 
        for n in range(len(num_points)):
            if len(area_gt_ind) >= num_points[n]:
    
                spl = 0
                for train_index, test_index in rcv.split(data_conc_norm.T):   
                    num_ind =  np.random.choice(area_gt_ind, size=num_points[n], replace=False)
                    X = data_conc_norm[num_ind, :].T
        
                    trainX, testX = X[train_index], X[test_index]
                    trainy, testy = y[train_index], y[test_index]
                    trainy_shuff = np.random.permutation(trainy)
                        
                    # fit a model - real data
                    classifier.fit(trainX, trainy)
                    if metric == 'AUC':
                        # predict probabilities for the positive outcome only
                        lr_probs = classifier.predict_proba(testX)[:, 1]
                        # calculate scores
                        auc_per_area[ar, n, spl] = roc_auc_score(testy, lr_probs)
                    if metric == 'Accuracy':
                    # This is to predict accuracy, not AUC
                        y_pred = classifier.predict(testX)
                        auc_per_area[ar, n, spl] = accuracy_score(testy, y_pred)
    
                    # # fit a model - shuffled data
                    # classifier.fit(trainX, trainy_shuff)
                    # shuff_probs = classifier.predict_proba(testX)[:, 1]
                    # auc_per_area_shuff[ar, spl] = roc_auc_score(testy, shuff_probs)
                    
                    spl = spl + 1
            
    auc_mean = np.mean(auc_per_area, axis=2)
    auc_sd = np.std(auc_per_area, axis=2)
    
    all_mean_val.append(auc_mean)
    all_sd_val.append(auc_sd)

# %%
# Now using all data plot per area

cat_gts = np.array(['C57', 'Esr', 'NPY', 'Vgl'])
#cat_colors = np.array(['black', 'darkorchid', 'deepskyblue', 'darkkhaki'])
cat_colors = np.array(['#414341', '#E38F8B', '#9472A3', '#4972A0'])

all_gt = np.array(all_gt)
all_ar = np.array(all_ar)

all_mean_val_conc = np.concatenate((all_mean_val))
all_sd_val_conc = np.concatenate((all_sd_val))

fig, axs = plt.subplots(1,4, figsize=(14, 4.5))
axs = axs.ravel()

#areas_to_plot = np.unique(all_ar, return_counts=True)
areas_to_plot = ['ACAd', 'PL', 'ILA', 'ORBm']
for a in range(len(areas_to_plot)):
    ar_ind = np.where(all_ar == areas_to_plot[a])[0]
    gt_present = all_gt[ar_ind]
    for i in range(len(ar_ind)):
        col_line = cat_colors[np.where(gt_present[i] == cat_gts)[0][0]]
        
        x_axis = num_points[np.where(all_mean_val_conc[ar_ind[i], :] != 0)]
        y_axis = all_mean_val_conc[ar_ind[i], np.where(all_mean_val_conc[ar_ind[i], :] != 0)][0]
        y_err =  all_sd_val_conc[ar_ind[i], np.where(all_mean_val_conc[ar_ind[i], :] != 0)][0]

        #axs[a].scatter(x_axis, y_axis, label = gt_present[i], color=col_line)
        #axs[a].errorbar(x_axis, y_axis, yerr = y_err, ecolor=col_line, fmt='None')
        axs[a].fill_between(x_axis, y_axis-y_err, y_axis+y_err, facecolor = col_line, alpha=0.4)
        axs[a].plot(x_axis, y_axis, color=col_line)

    axs[a].set_ylim(0.4,1)
    axs[a].set_xlim(0, np.max(num_points))
    axs[a].set_ylabel(metric)
    axs[a].set_xlabel('Number of units')
    if a == 3:
        axs[a].legend(cat_gts)
    axs[a].set_title(areas_to_plot[a])
    axs[a].axhline(0.5, color='black', lw=1, ls = '--', alpha=0.3)
    
fig.tight_layout()
figname = os.path.join(save_folder, 'LogReg_per_area_Accuracy_bl{0}-bl{1}.svg'.format(str(to_decode[0]), str(to_decode[1])))
#fig.savefig(figname)

