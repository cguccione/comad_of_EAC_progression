'''
Helper functions used to calculate basic
alpha and beta diversity
'''

import pandas as pd
import os
from qiime2 import Artifact, Metadata

#Plotting
import matplotlib.pyplot as plt 
import seaborn as sns

#Alpha
from scipy.stats import mannwhitneyu
from qiime2.plugins.feature_table.methods import rarefy
from qiime2.plugins.diversity.pipelines import alpha

#Beta
import subprocess
from qiime2.plugins.emperor.visualizers import biplot
from qiime2.plugins import diversity

def alpha_plot(df, metric, fn):
    
    #Export everything for future use
    export_filepath = os.getcwd() + '/outputs/alpha_plots/'
    
    metric_list = list(set(df[metric]))
    print(metric_list)
    
    #Print helpful metrics
    print('Number of samples per group:')
    print(df.loc[:][metric].value_counts())
    print()
    
    #Observed p-values
    list_0 = list(df.loc[df[metric] ==  metric_list[0]]['observed_features'])
    list_1 = list(df.loc[df[metric] ==  metric_list[1]]['observed_features'])
    
    #Observed p-values
    _, pnorm = mannwhitneyu(list_0, list_1)
    opv = 'Mann Whiteny U rank test: p-value =' + str(pnorm)
    print(opv)
    
    #Observed features plotting
    sns.set_style('whitegrid')
    plt.figure(figsize=(8,6))
    sns.set_context(font_scale=1.5)
    sns.violinplot(y = 'observed_features', x = metric, data = df)
    sns.stripplot(y = 'observed_features', x = metric, data = df, color='black')
    sns.boxplot(y = 'observed_features', x = metric, data = df, width = 0.2, color='grey', showfliers = False)
    plt.xlabel('Sample Type')#metric)
    plt.ylabel('Observed Features')
    plt.gcf().text(0.3, 0.9, opv, fontsize=14, fontweight="bold")
    plt.rcParams['svg.fonttype'] = 'none'
    if fn != False:
        plt.savefig(export_filepath + str(fn) + '_observed.svg', format = 'svg')
    
    plt.show()
    
    plt.clf()
    
    #Shannon p-values
    list_0 = list(df.loc[df[metric] ==  metric_list[0]]['shannon_entropy'])
    list_1 = list(df.loc[df[metric] ==  metric_list[1]]['shannon_entropy'])
    
    #Shannon p-values
    _, pnorm = mannwhitneyu(list_0, list_1)
    spv = 'Mann Whiteny U rank test: p-value =' + str(pnorm)
    print(spv)
    
    #Shannon plotting
    sns.set_style('whitegrid')
    plt.figure(figsize=(8,6))
    sns.set_context(font_scale=1.5)
    sns.violinplot(y = 'shannon_entropy', x = metric, data = df)
    sns.stripplot(y = 'shannon_entropy', x = metric, data = df, color='black')
    sns.boxplot(y = 'shannon_entropy', x = metric, data = df, width = 0.2, color='grey', showfliers = False)
    plt.xlabel('Sample Type') #metric)
    plt.ylabel('Shannon Entropy')
    plt.gcf().text(0.3, 0.9, spv, fontsize=14, fontweight="bold")
    plt.rcParams['svg.fonttype'] = 'none'
    if fn != False:
        plt.savefig(export_filepath + str(fn) + '_shannon.svg', format = 'svg')
    

def alpha_diversity(table_list, rarefaction, metric, fn):
    
    ##would make this into function eventally because used in beta too
    table = pd.DataFrame()
    meta = pd.DataFrame(columns = ['sample_type'])
    for i in table_list:
        #Import biom table in pandas df/tsv form for current dataset
        biom_df = pd.read_csv('processed_data/pandas_df/' + i + '.tsv', sep = '\t', index_col=0)
        biom_df = biom_df.astype(int)
        table = pd.concat([table, biom_df], axis = 1)

        #Import metadata for current dataset
        meta_temp = pd.read_csv('processed_data/metadata/metadata_' + i + '.tsv', sep = '\t')
        meta_temp['sample_type'] = i
        meta = pd.merge(meta, meta_temp, how='outer')
    
    table = table.fillna(0)
    
    #Convert into q2 object / correct formatting
    table = table.T
    ft = Artifact.import_data("FeatureTable[Frequency]", table)
    meta = meta.set_index('sample_name', drop = True)
    
    #(Optional) Rarefaction
    if rarefaction != False:
        ft = rarefy(table=ft, sampling_depth = rarefaction)
        ft = ft.rarefied_table

    #Alpha diverstiy, observed
    alpha_result_o = alpha(table=ft, metric='observed_features')
    alpha_diversity_o = alpha_result_o.alpha_diversity
    alpha_series_o = alpha_diversity_o.view(pd.Series)
    df = pd.merge(meta, alpha_series_o, right_index = True, left_index = True)

    #Alpha diverstiy, Shannon
    alpha_result_s = alpha(table=ft, metric='shannon')
    alpha_diversity_s = alpha_result_s.alpha_diversity
    alpha_series_s = alpha_diversity_s.view(pd.Series)
    df = pd.merge(df, alpha_series_s, right_index = True, left_index = True)
    
    '''
    if tree != False:
        #Alpha diversity, Faiths
        alpha_diversity_f, = faith_pd(table=ft, phylogeny = tree)
        alpha_series_f = alpha_diversity_f.view(pd.Series)
        df = pd.merge(df, alpha_series_f, right_index = True, left_index = True)  
    '''
    
    alpha_plot(df, metric, fn)
    
    
def deicode_beta_diversity(table_list, metric, fn):
    
    ##would make this into function eventally because used in beta too
    table = pd.DataFrame()
    meta = pd.DataFrame(columns = ['sample_type'])
    for i in table_list:
        #Import biom table in pandas df/tsv form for current dataset
        biom_df = pd.read_csv('processed_data/pandas_df/' + i + '.tsv', sep = '\t', index_col=0)
        biom_df = biom_df.astype(int)
        table = pd.concat([table, biom_df], axis = 1)

        #Import metadata for current dataset
        meta_temp = pd.read_csv('processed_data/metadata/metadata_' + i + '.tsv', sep = '\t')
        meta_temp['sample_type'] = i
        meta = pd.merge(meta, meta_temp, how='outer')
    
    table = table.fillna(0)
    
    #Convert table into q2 object and export to qza
    table = table.T
    ft = Artifact.import_data("FeatureTable[Frequency]", table)
    ft.save('deicode_processing.qza')
    
    #Put sample data into correct formatting for q2
    meta = meta.set_index('sample_name')
    
    #Turn all true/false into strings
    mask = meta.applymap(type) != bool
    d = {True: 'True', False: 'False'}
    meta = meta.where(mask, meta.replace(d))
    
    sample_meta = Metadata(meta)
    
    #Run deicode from command line using qza created above
    subprocess.run(["qiime", "deicode", "rpca", "--i-table", "deicode_processing.qza", "--o-biplot", "deicode_biplot.qza", "--o-distance-matrix", "deicode_distance_test.qza"])
        
    #Import biplot back into python
    rpca_biplot = Artifact.load('deicode_biplot.qza')
    
    #Import biplot back into python
    rpca_distance_matrix = Artifact.load('deicode_distance_test.qza')
    
    #Create emperor visualization from the biplot result
    rpca_biplot_emperor = biplot(biplot = rpca_biplot, sample_metadata = sample_meta)
    
    #Calculate permanova 
    beta_result_o = diversity.actions.beta_group_significance(distance_matrix=rpca_distance_matrix, metadata=sample_meta.get_column(metric), method = 'permanova', pairwise = True)
    
    return(rpca_biplot_emperor, beta_result_o)
    