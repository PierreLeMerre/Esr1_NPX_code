# -*- coding: utf-8 -*-

# This script requires the output from 'Extended_Data_Fig9a_binning10ms.py'

from pynwb import NWBFile, NWBHDF5IO

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import os.path
import sys
import scipy.io
import h5py
from scipy import stats as st # z-scoring
from sklearn.preprocessing import normalize # L1 and L2 norm
# import seaborn as sns; sns.set()

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from scipy.ndimage.filters import gaussian_filter1d
from mpl_toolkits.mplot3d import Axes3D
 
# Import my custom functions
sys.path.append('C:\\Users\\admin\\OneDrive - KI.SE\\Skrivbordet\\General\\Projects\\Aversion\\neuropixelPFC\\Code_python\\')
from MS_suppl_functions import *

## ## TO DISPLAY 3d PLOT IN BROWSER
import plotly
import plotly.io as pio
from plotly.subplots import make_subplots
import plotly.graph_objs as go

import plotly.express as px

pio.renderers.default='browser'

# switch back to spyder
# pio.renderers.default='svg'

#%% 
all_genotypes = ['NPY-Cre', 'C57BL-6J', 'Vglut2-cre', 'Esr1-Cre-antero', 'Esr1-Cre-retro']
# Select the blocks that you want to plot
blocks = np.array((1, 1, 1, 1), dtype='bool')

save_folder = 'C:\\Users\\admin\\OneDrive - KI.SE\\Skrivbordet\\General\\Projects\\Aversion\\neuropixelPFC\\Code_python\\PCA_virtual_mice_zscore\\'
all_units_data = np.load(save_folder + 'all_units_binned_10_QC.npy')
area_per_unit = np.load(save_folder + 'area_per_unit_short_QC.npy')
genotype_per_unit = np.load(save_folder + 'genotype_per_unit_QC.npy')

nb = np.shape(all_units_data)[1]//len(blocks)
use_trials = [[bl] * nb for bl in blocks]
use_trials = np.concatenate((use_trials))

trials_labels = np.array(['sound_pre', 'sound + opto', 'sound_post', 'sound + air-puff'])
#trials_colors = np.array(['forestgreen', 'crimson', 'gold', 'navy'])
#colorscales = np.array(['Greens', 'Reds', 'YlOrBr', 'Blues'])

print(np.shape(all_units_data))

#%%

# z-score of each unit
data_mean = np.mean(all_units_data, axis=1)
data_std = np.std(all_units_data, axis=1)

all_units_data_zs = np.nan_to_num((all_units_data - data_mean[:, None]) / data_std[:, None])

# smoothing of each unit with gaussian kernel separately for blocks and then concatenating them back

# Plotting of smoothed traces of 1 unit
fig, axs = plt.subplots(2,2, figsize=(14, 10)) ##
axs = axs.ravel()##

all_data_smoothed = []
for bl in range(len(blocks)):
   
    blocks_to_use = np.array((0,0,0,0), dtype='bool')
    blocks_to_use[bl] = True
    use_trials = np.concatenate(([[bl] * nb for bl in blocks_to_use]))
    
    block_smoothed = gaussian_filter1d(all_units_data_zs[:, use_trials], sigma = 2, axis=1)
    all_data_smoothed.append(block_smoothed)
    
    axs[bl].plot(all_units_data_zs[0, use_trials])
    axs[bl].plot(block_smoothed[0, :])
    
    
all_units_data = np.hstack(np.array(all_data_smoothed))
print(np.shape(all_units_data), '\n')

print(np.unique(area_per_unit, return_counts=True))
gen_short = [i[:3] for i in genotype_per_unit]
print(np.unique(gen_short, return_counts=True))

#%%

# Plot all genotypes in the same PCA space for each genotype

# Subsample 500+ random cell from each genotype and concatenate (4 x 1000) x T (timepoints)
# The number of cells in the subset is defined as 90% of random cells in the genotype with the least number of cells
un_counts = np.unique(gen_short, return_counts=True)
unit_sub = int(np.min(un_counts[1]) - np.min(un_counts[1])*0.1)

gen_short_labels = ['NPY-Cre', 'C57BL-6J', 'Vglut2-cre', 'Esr1-Cre']
gt_colors = np.array(['#9472A3', '#414341', '#4972A0', '#E38F8B'])
block_labels = ['Block1', 'Block2', 'Block3', 'Block4']
#colorscales = np.array(['Greens', 'Reds', 'YlOrBr', 'Blues'])

# Join Esr1 antero- and retrograde tracing into 1 genotype

all_gen_short = ['NPY', 'C57', 'Vgl', 'Esr']
gen_per_unit_short = [x[:3] for x in genotype_per_unit]


# Plotting and saving virtual mice separately (also for ESR1 antero and retrograde)
# Labels of the subplots (genotypes) are in the wrong order here
fig = make_subplots(
      rows=2, cols=2,
      specs=[[{'type': 'scene'}, {'type': 'scene'}],
             [{'type': 'scene'}, {'type': 'scene'}]],
      subplot_titles = block_labels)

r = 1
c = 1
for bl in range(len(blocks)):
    
    blocks_to_use = np.array((0,0,0,0), dtype='bool')
    blocks_to_use[bl] = True
    use_trials = [[bl] * nb for bl in blocks_to_use]
    use_trials = np.concatenate((use_trials))
    
    all_genotype_pca = []
    for gt in range(len(all_gen_short)):
        # select all units from the big table belonging to the same genotype
        gt_units = all_units_data[np.in1d(gen_per_unit_short, all_gen_short[gt]), :]
        gt_areas = area_per_unit[np.in1d(gen_per_unit_short, all_gen_short[gt])]
        #print(np.unique(gt_areas))
    
        ind_array = np.zeros((len(gt_units)))
        ind_array[np.random.choice(len(gt_units), unit_sub, replace=False)] = 1
    
        xx = gt_units[np.array(ind_array, dtype='bool'), :]
        all_genotype_pca.append(xx[:, use_trials])

    all_genotype_pca = np.vstack(np.array(all_genotype_pca))
    print(np.shape(all_genotype_pca))   

    # fitting PCA on trials averaged per block (different genotypes) and concatenated
    pca2 = PCA(n_components=10)

    PC_av_per_cond_all = pca2.fit_transform(all_genotype_pca.T)
    #PC_av_per_cond = PC_av_per_cond_all[:, :3]
    var2=np.cumsum(np.round(pca2.explained_variance_ratio_, decimals=3)*100)
    print(var2)
    
    PC_loadings = pca2.components_

    PC_av_per_cond = np.zeros((nb, len(all_gen_short), 10))
    g_ind = 0
    for i in range(len(all_gen_short)):
        PC_av_per_cond[:, i, :] = np.matmul(all_genotype_pca[g_ind:g_ind+unit_sub].T, pca2.components_.T[g_ind:g_ind+unit_sub])
        g_ind = g_ind + unit_sub
       
    #PC_av_per_cond_plot = np.reshape(PC_av_per_cond, (nb, len(all_gen_short), 3), order='F')
    PC_av_per_cond_plot = PC_av_per_cond[:, :, :3]
    print(np.shape(PC_av_per_cond_plot))

    for cnd in range(len(gen_short_labels)):
    
        fig.add_trace(go.Scatter3d(x=PC_av_per_cond_plot[:, cnd, 0], 
                                   y=PC_av_per_cond_plot[:, cnd, 1], 
                                   z=PC_av_per_cond_plot[:, cnd, 2], 
                                   mode='lines+markers',
                                   name=gen_short_labels[cnd],
                                   #marker=dict(size=4, color=np.linspace(0, nb), colorscale=colorscales[blocks][cnd]), # one of plotly colorscales),
                                   marker=dict(size=4, color=gt_colors[cnd], opacity=0.7),
                                   line=dict(color=gt_colors[cnd], width = 10),
                                   opacity = 0.7),
                     row=r, col=c)

    fig.update_layout(
        title =  str(unit_sub) + ' units',
        #scene = dict(xaxis_title='PC1; ' + str(var2[0])[:4] + '% cv',
                     #yaxis_title='PC2; ' + str(var2[1])[:4] + '% cv',
                     #zaxis_title='PC3; ' + str(var2[2])[:4] + '% cv'),
        scene = dict(xaxis_title='PC1',
                     yaxis_title='PC2',
                     zaxis_title='PC3'),
        #margin=dict(l=0, r=0, t=0, b=0),
        #width = 2000,
        #height = 1600,
        width = 2000,
        height = 1600,
        legend_title="Genotype",
        template="plotly_white")
    
    print(r, c)
    
    if bl == 0:
        c += 1
    if bl == 1:
        c = c - 1
        r += 1
    if bl == 2:
        c += 1
    
fig.show()
# enable saving if you need
#fig.write_html(save_folder + 'All_virt_10ms_QC' + ".html")