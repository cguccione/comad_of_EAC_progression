'''
Helper functions used to create selection plots
'''

#Basic imports
import numpy as np
import pandas as pd
from biom import load_table
import os

from python_scripts.neufit import *

def read_npz(file_path):
    # Load the .npz file
    npz_data = np.load(file_path)
    
    print(file_path.split('/')[-1].split('.npz')[0])

    # Access and print individual arrays
    for key in npz_data.keys():
        if key not in ['p', 'gR', 'dR', 'timeseries_time', 'timeseries_data', 'P0_nneutral', 'mean_freq_nneutral']:
            print(f"{key}:", npz_data[key])
    
    print ('HP_inital', npz_data['p'][0])
    return(npz_data)

def selection_plot(selection_output, save_figure=False, color_HP=True, color_t15=False):
    
    selection_fn = 'outputs/simulation_output/' + selection_output + '.npz'
    
    npz_data = read_npz(selection_fn)
    
    #Extract useful info out of npz_data
    n_reads = int(npz_data['N'])
    n_samples = int(npz_data['S'])
    m_value = int(npz_data['m'])
    
    #Remove 0 elements 
    ma_v = list(npz_data['mean_freq_nneutral'])
    o_v = list(npz_data['P0_nneutral'])
    ma = []
    o = []
    for i in range(0, len(ma_v)):
        if str(ma_v[i]) != '0.0':
            ma.append(ma_v[i])
            o.append(o_v[i])
    
    occurr_freqs = pd.DataFrame()
    occurr_freqs['mean_abundance'] = ma
    occurr_freqs['occurrence'] = o
    n_samples = len(ma)
    
    current_index = ['Unknown'] * len(occurr_freqs.index)
    
    if color_HP == True:
        #Add H.pylori to taxonomy
        current_index[0] = 'k__Bacteria;p__Proteobacteria;c__Epsilonproteobacteria;o__Campylobacterales;f__Helicobacteraceae;g__Helicobacter;s__Helicobacter pylori'
    
    '''
    else:
        taxa_df = pd.read_csv('outputs/occur_freqs/EAC_ICGC_wo-RI_species_WOL_scrubbed.tsv', sep='\t')

        current_index = list(taxa_df['full_taxonomy'])
    '''
    
    if color_t15 == True:
        for i in range(15):
            current_index[i] = 'k__Bacteria;p__Proteobacteria;c__Epsilonproteobacteria;o__Campylobacterales;f__Helicobacteraceae;g__Helicobacter;s__Helicobacter pylori'
    
    occurr_freqs['full_taxonomy'] = current_index
    
    beta_fit, r_square = neufit_calc(occurr_freqs, n_reads, n_samples, print_non_netural=False, print_report=True)
    
    neufit_plot(occurr_freqs, beta_fit, n_samples, n_reads, r_square, selection_output, sim_input=False, save_occur=False, from_simulation=True, save_plot=save_figure)
