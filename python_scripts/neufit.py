'''
Helper functions used to run Neufit
'''

#Basic imports
import numpy as np
import pandas as pd
from biom import load_table
import os

#Neufit 
from lmfit import Parameters, Model, fit_report
from scipy.stats import beta
from statsmodels.stats.proportion import proportion_confint
from matplotlib import pyplot
from math import log10

#Rarefaction calc
from qiime2 import Metadata, Artifact
from qiime2.plugins import feature_table, diversity
from qiime2.plugins.feature_table.visualizers import summarize
from qiime2.plugins.diversity.visualizers import alpha_rarefaction

from matplotlib.ticker import MultipleLocator

#Import taxonomy (WOL)
#tax_in = pd.read_csv('/projects/wol/qiyun/wol2/taxonomy/lineages.txt', sep='\t', header=None, index_col=0)
tax_in = pd.read_csv('/Users/cguccion/Dropbox/Storage/HelpfulLabDocs/taxonomy_trees/WOL/lineages.txt', sep='\t', header=None, index_col=0)
tax_dict = tax_in[1].to_dict() 

pyrlori_name='k__Bacteria;p__Proteobacteria;c__Epsilonproteobacteria;o__Campylobacterales;f__Helicobacteraceae;g__Helicobacter;s__Helicobacter pylori'

#Strangly high gotus in ICGC run2
bad_gotu = ['G001408515', 'G001458175', 'G000601485', 'G000824785', 'G001748365', 'G001404055', 'G001049355', 'G000824825', 'G001458375', 'G001748465', 'G001458355', 'G000751555', 'G001049735', 'G000349545', 'G001458315', 'G000784355', 'G000310185', 'G001458335', 'G000724545', 'G900095615', 'G001261715', 'G001404515', 'G001403735', 'G000338055', 'G001458215', 'G000827595', 'G000613185', 'G000611675', 'G000333235', 'G001403795', 'G001258055', 'G000410875', 'G000582705', 'G000443165', 'G001282945', 'G001417815', 'G000497205', 'G000160455', 'G001458435', 'G000829675', 'G001283065', 'G000335475', 'G001458395', 'G001458415']
#G001273835, G001373415, G000417585

def calc_abdundance(biom_df, ignore_level):
    abundances = biom_df[biom_df.sum(1) > ignore_level]
    print ('dataset contains ' + str(abundances.shape[1]) + ' samples (sample_id, reads):')
    print (abundances.sum(0), '\n')
    print(abundances.sum(0)[:50])
    print(abundances.sum(0)[:51-80])
    return(abundances)

def subsample(counts, depth):
    # Subsamples counts to uniform depth, dropping all samples without enough depth
    for sample in counts:
        if counts[sample].sum() >= depth:
            flattened = np.repeat(np.arange(counts[sample].size), counts[sample])
            subsample = np.random.choice(flattened, depth, replace=False)
            counts[sample] = np.bincount(subsample, minlength=counts[sample].size)
        else:
            print ('dropping sample ' + sample + ' with ' + str(counts[sample].sum()) + ' reads < ' + str(depth))
            counts = counts.drop(sample, axis=1)
    return counts
    

def rarefraction(biom_df, ignore_level, rarefaction_level):
    abundances = calc_abdundance(biom_df, ignore_level)
    
    if rarefaction_level == 0 or rarefaction_level > max(abundances.sum(0)):
        rarefaction_level = min(abundances.sum(0))
        print('rarefying to highest possible uniform read depth')
    else:
        print ('rarefying to custom rarefaction level')
    print ('(' + str(rarefaction_level) + ' reads per sample)')
    
    if not all(n_reads == rarefaction_level for n_reads in abundances.sum(0)):
        abundances = subsample(abundances, rarefaction_level)
        abundances = abundances[abundances.sum(1) > 0]
    
    n_otus, n_samples = abundances.shape
    n_reads = rarefaction_level

    print ('\nfitting neutral expectation to dataset with ' + str(n_samples) + ' samples and ' + str(n_otus) + ' otus')
    
    '''
    print('abundances')
    display(abundances)
    abundances.to_csv('tmp.tsv', sep='\t')
    temp_df = abundances.loc['k__Bacteria;p__Bacteroidetes;c__Cytophagia;o__Cytophagales;f__Hymenobacteraceae;g__Hymenobacter;s__Hymenobacter sp. IS2118']
    display(temp_df[:50])
    '''
    
    return(abundances, n_samples, n_reads)

''' Will have to update this to not just be WOL like we have here
def add_taxonomy(occurr_freqs, taxonomy_file = 'neufit_inputs/lineages.txt'):
    taxonomy = pd.read_table(taxonomy_file, header=0, index_col=0, sep='\t')
    taxonomy.index.name = 'otu_id'
    occurr_freqs = occurr_freqs.join(taxonomy)
    return(occurr_freqs)
'''

def calc_mean_realtive_abudnace(abundances, n_samples, n_reads, rarefaction_level, add_tax=False):
    mean_relative_abundance = (1.0*abundances.sum(1))/n_reads/n_samples
    
    occurrence_frequency = (1.0*np.count_nonzero(abundances, axis=1))/n_samples

    occurr_freqs = pd.DataFrame(mean_relative_abundance, columns=['mean_abundance'])
    occurr_freqs['occurrence'] = occurrence_frequency
    #occurr_freqs.index.name = 'otu_id'
    
    occurr_freqs = occurr_freqs.sort_values(by=['mean_abundance'])
    display(occurr_freqs)

    return(occurr_freqs)

def beta_cdf(p, N, m):
    # Expected long term distribution under the neutral model (truncated cumulative beta-distribution)
    return beta.cdf(1.0, N*m*p, N*m*(1.0-p)) - beta.cdf(1.0/N, N*m*p, N*m*(1.0-p))

def neufit_calc(occurr_freqs, n_reads, n_samples, print_non_netural=False, print_report=True):
    
    #Fit the model
    params = Parameters()
    params.add('N', value=n_reads, vary=False)
    params.add('m', value=0.5, min=0.0, max=1.0)
    beta_model = Model(beta_cdf)
    beta_fit = beta_model.fit(occurr_freqs['occurrence'], params, p=occurr_freqs['mean_abundance'])
    
    r_square = 1.0 - np.sum(np.square(occurr_freqs['occurrence'] - beta_fit.best_fit))/np.sum(np.square(occurr_freqs['occurrence'] - np.mean(occurr_freqs['occurrence'])))
    
    if print_report == True:
        print(fit_report(beta_fit))
        print ('R^2 = ' + '{:1.2f}'.format(r_square))
    
    occurr_freqs['predicted_occurrence'] = beta_fit.best_fit
    occurr_freqs['lower_conf_int'], occurr_freqs['upper_conf_int'] = proportion_confint(occurr_freqs['predicted_occurrence']*n_samples, n_samples, alpha=0.05, method='wilson')
    
    if print_non_netural == True:
        above = occurr_freqs[occurr_freqs['occurrence'] > occurr_freqs['upper_conf_int']]
        below = occurr_freqs[occurr_freqs['occurrence'] < occurr_freqs['lower_conf_int']]
        print('Above uppper conf interval:', above)
        print()
        print('Below lower conf interval', below)
        
    return(beta_fit, r_square)

def neufit_plot(occurr_freqs, beta_fit, n_samples, n_reads, r_square, fn, HP_color = False, save_plot=True, sim_input=True, save_occur=True, from_simulation=False, non_color=False):
    pyplot.cla() #Clears previous plot - to avoid double keys
    
    font_size=18#5

    #Prepare results plot
    pyplot.xlabel('Mean relative abundance across samples', fontsize=font_size)#18
    pyplot.xscale('log')
    x_range = np.logspace(log10(min(occurr_freqs['mean_abundance'])/10), 0, 1000)
    pyplot.xlim(min(x_range), max(x_range))
    pyplot.xticks(fontsize=font_size)#18
    
    '''
    #pyplot.xticks([10**i for i in range(-6, 1)], fontsize=font_size)#18
    #pyplot.gca().xaxis.set_minor_locator(MultipleLocator(0.1))
    
    major_ticks = [10**i for i in range(-6, 1)]
    pyplot.xticks(major_ticks, fontsize=font_size)
    
    minor_ticks = []
    for i in range(len(major_ticks) - 1):
        minor_ticks.extend(np.linspace(major_ticks[i], major_ticks[i+1], num=9, endpoint=False)[1:])
    pyplot.gca().xaxis.set_minor_locator(MultipleLocator())  # Ensure minor ticks are spaced evenly
    '''

    pyplot.ylabel('Occurrence frequency in samples', fontsize=font_size)#18
    pyplot.ylim(-0.05, 1.05)
    pyplot.yticks(fontsize=font_size)#16

    #>>Calculate lower and upper range
    lower, upper = proportion_confint(beta_cdf(x_range, n_reads, beta_fit.best_values['m'])*n_samples, n_samples, alpha=0.05, method='wilson')

    #>>Calculate/plot the main fit line
    #Orginal plotting colors
    pyplot.plot(x_range, beta_cdf(x_range, n_reads, beta_fit.best_values['m']), '-', lw=5, color='darkred')
    pyplot.plot(x_range, lower, '--', lw=2, color='darkred')
    pyplot.plot(x_range, upper, '--', lw=2, color='darkred')
    pyplot.fill_between(x_range, lower, upper, color='lightgrey')
    pyplot.plot([1e-4, 1e-4], [-0.05, 1.05], color='black', linestyle='-', linewidth=1)
    
    #Set Marker size
    markersize = 5

    # Plot data points
    pyplot.plot(occurr_freqs['mean_abundance'], occurr_freqs['occurrence'], 'o', markersize=markersize, fillstyle='full', color='black')
    
    #Extract m value
    m = fit_report(beta_fit).split('(fixed)')[1]
    m = 'm = ' + str(round(float(m.strip('\n').split('+')[0].strip(' ').strip('m: ')),3))

    #Plot R^2 and m values
    #pyplot.text(0.05, 0.9, m, fontsize=16, transform=pyplot.gca().transAxes)
    pyplot.text(0.05, 0.9, '$R^2 = ' + '{:1.2f}'.format(r_square) + '$', fontsize=font_size, transform=pyplot.gca().transAxes)
    
    if 'full_taxonomy' not in occurr_freqs.columns:
        occurr_freqs['full_taxonomy'] = occurr_freqs.index
    
    #print(occurr_freqs['full_taxonomy'].iloc[0])
    print('-------')
    
    if 'p__' not in occurr_freqs['full_taxonomy'].iloc[0]:
        #Add full taxonomy instead of just GOTU
        occurr_freqs['full_taxonomy'] = occurr_freqs.index.map(lambda x: tax_dict.get(x, 'Unknown'))
        occurr_freqs['full_taxonomy'] = occurr_freqs.index
    
    display(occurr_freqs)
    
    HP_color_SAM = False
    EC_color = False
    color_key_out = 'full_taxonomy'
    #Adding phlya coloring
    for index, col in occurr_freqs.iterrows():
        if not isinstance(col['full_taxonomy'], str):
            continue
        if 's__Helicobacter pylori' in col['full_taxonomy']:
            HP_color = True
            print(col)
        if 'Helicobacter_pylori' in col['full_taxonomy']:
            HP_color_SAM = True
            print(col)
        if non_color == True:
            continue
        if 's__Escherichia coli' in col['full_taxonomy']:
            #EC_color = True
            continue
        if 'p__Firmicutes' in col['full_taxonomy']:
            pyplot.plot(col['mean_abundance'], col['occurrence'], 'o', 
                markersize=markersize, fillstyle='full', color='blue')
        elif 'p__Proteobacteria' in col['full_taxonomy']:
            pyplot.plot(col['mean_abundance'], col['occurrence'], 'o', 
                markersize=markersize, fillstyle='full', color='orange')
        elif 'p__Fusobacteria' in col['full_taxonomy']:
            pyplot.plot(col['mean_abundance'], col['occurrence'], 'o', 
                markersize=markersize, fillstyle='full', color='lightblue')
        elif 'p__Actinobacteria' in col['full_taxonomy']:
            pyplot.plot(col['mean_abundance'], col['occurrence'], 'o', 
                markersize=markersize, fillstyle='full', color='yellow')
        elif 'p__Bacteroidetes' in col['full_taxonomy']:
            pyplot.plot(col['mean_abundance'], col['occurrence'], 'o', 
                markersize=markersize, fillstyle='full', color='grey')
        elif col['full_taxonomy'] in bad_gotu:
            pyplot.plot(col['mean_abundance'], col['occurrence'], 'o', 
                markersize=markersize, fillstyle='full', color='yellow')
            
    
    #Plot SAM HP coloring
    if HP_color_SAM != False:
        hp_occurr_freqs = occurr_freqs[occurr_freqs['full_taxonomy'].str.contains('Helicobacter_pylori', case=False)].copy()
        pyplot.plot(hp_occurr_freqs['mean_abundance'], hp_occurr_freqs['occurrence'], 'o', markersize=markersize, fillstyle='full', color='magenta')    
    
    #Plot HP last so we can see coloring
    if HP_color != False:
        hp_occurr_freqs = occurr_freqs[occurr_freqs['full_taxonomy'].str.contains('s__Helicobacter pylori', case=False)].copy()
        pyplot.plot(hp_occurr_freqs['mean_abundance'], hp_occurr_freqs['occurrence'], 'o', markersize=markersize, fillstyle='full', color='magenta')
        
    if EC_color != False:
        ec_occurr_freqs = occurr_freqs[occurr_freqs['full_taxonomy'].str.contains('s__Escherichia coli', case=False)].copy()
        pyplot.plot(ec_occurr_freqs['mean_abundance'], ec_occurr_freqs['occurrence'], 'o', markersize=markersize, fillstyle='full', color='green')

    #Run standout microbes (optional)
    standout_microbes(occurr_freqs, fn, save_plot=save_plot)
    
    if save_plot != False:
        #pyplot.rcParams.update({'font.size': 5})
        
        if from_simulation == False:
            #Save plot
            pyplot.tight_layout()
            pyplot.gcf().set_size_inches(6,4)#(1.4, 0.98)#(7, 5)
            output='outputs/neufit_plots/'
            pyplot.savefig(output)
            pyplot.savefig(output + save_plot + '.png')
            pyplot.savefig(output + save_plot + '.pdf')
        else:
            #Save plot
            pyplot.tight_layout()
            pyplot.gcf().set_size_inches(6,4)#(1.4, 0.98)#(7, 5)
            output='outputs/simulation_neufit_plots/'
            pyplot.savefig(output + save_plot + '.png')
            pyplot.savefig(output + save_plot + '.pdf')
        
        #Save in svg format
        pyplot.rcParams['font.family'] = 'sans-serif'
        pyplot.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
        pyplot.rcParams['svg.fonttype'] = 'none'  # render SVG text as text, not curves
        pyplot.savefig(output + save_plot + '.svg')
        
    #Save df 
    if save_occur != False:
        output='outputs/occur_freqs/'
        occurr_freqs.to_csv(output + save_occur + '.tsv', sep='\t')
    
    #Save simulation input df
    if HP_color != False:
        #Move H.pylori to top row of dataframe
        occurr_freqs['hp'] = range(1,len(occurr_freqs)+1)
        occurr_freqs.loc[pyrlori_name, 'hp']=0
        occurr_freqs = occurr_freqs.sort_values('hp').drop('hp', axis = 1)
    
    if sim_input != False:
        sim='outputs/simulation_input/'
        occurr_freqs[['mean_abundance']].T.to_csv(sim + sim_input + '.csv', sep=',', index=False, header=False)

    pyplot.show()
    
def standout_microbes(occurr_freqs, fn, save_plot, threshold=0.1):
    #Create dataframe
    standoutMicrobes = pd.DataFrame(columns = ('Difference off Neutral Model',
                                               'Taxonomy', 'full_taxonomy', 'mean_abundance', 'occurrence'))
    
    #Loop and find most non-neutral microbes based upon threshold
    row_count = 0
    for index, col in occurr_freqs.iterrows():
        diff = abs(col['occurrence'] - col['predicted_occurrence'])
        if diff > threshold:
            standoutMicrobes.loc[row_count] = [diff, index, col['full_taxonomy'], col['mean_abundance'], col['occurrence']]
            row_count+=1
    
    standoutMicrobes = standoutMicrobes.sort_values(by =['Difference off Neutral Model'],
                                                    ascending=False)
    
    #Display and export non-neutral microbes as csv
    print("\nTop NonNeutral Microbes")
    display(standoutMicrobes)
    output='outputs/non_neutral/'
    
    if save_plot != False:
        standoutMicrobes.to_csv(output + save_plot + '.tsv', sep = '\t', index = False)
    
def calculate_summary_table(fn):
    
    #Import metadata and biom table into qiime2 format
    meta = Metadata.load(str('processed_data/metadata/metadata_' + fn + '.tsv'))
    ft = Artifact.load(str('processed_data/qza/' + fn + '.qza'))
    
    #Display table summary to see which rarefaction level to select 
    table_summary = feature_table.actions.summarize(table = ft,
                                                    sample_metadata = meta)
    return(table_summary)
    
def calculate_rarefaction(max_depth_rare, fn, steps=10):
    
    #Import metadata and biom table into qiime2 format
    meta = Metadata.load(str('processed_data/metadata/metadata_' + fn + '.tsv'))
    ft = Artifact.load(str('processed_data/qza/' + fn + '.qza'))
    
    alpha_rare = diversity.actions.alpha_rarefaction(table = ft,
                                                     max_depth = max_depth_rare,
                                                     metadata = meta,
                                                     steps = steps,
                                                     metrics = {'observed_features'})
    
    print('ft is here')
    display(ft)
    print('')
    
    return(alpha_rare)
    
def neufit_main(rarefaction_level, fn, ignore_level=0, taxonomy= None, non_color=False, save=False):
    
    #Import biom table in pandas df/tsv form for current dataset
    biom_df = pd.read_csv('processed_data/pandas_df/' + fn + '.tsv', sep = '\t', index_col=0)
    biom_df = biom_df.astype(int)
    display(biom_df[:3])

    #Import metadata for current dataset
    meta = pd.read_csv('processed_data/metadata/metadata_' + fn + '.tsv', sep = '\t')
    
    #Import data and calculate rarefaction
    ab, n_samples, n_reads = rarefraction(biom_df, ignore_level, rarefaction_level)
    
    #Calculate mean relative abundance
    occurr_freqs = calc_mean_realtive_abudnace(ab, n_samples, n_reads, rarefaction_level)
    
    #Neufit Calculation 
    beta_fit, r_square = neufit_calc(occurr_freqs, n_reads, n_samples)
    
    if save != False:
        save= save + '_r' + str(rarefaction_level)
    
    #Neufit Plotting
    if non_color==True:
        neufit_plot(occurr_freqs, beta_fit, n_samples, n_reads, r_square, fn + '_nonColor' + '_r' + str(rarefaction_level), non_color=True, save_plot=save, sim_input=save, save_occur=save)
        
    else:
        neufit_plot(occurr_freqs, beta_fit, n_samples, n_reads, r_square, fn + '_r' + str(rarefaction_level), save_plot=save, sim_input=save, save_occur=save)#, True) #Save file here'experimental_outputs/Hutch_combined_WOL_neufit.png')

    ''' Taxonomy file only for WOL rn .. I don't think this is used later but should check
    Will probs want to update this to [WOL, REP200... ect and then have all the files stored
    here for use
    
    #Import taxonomy file
    taxonomy_file = 'neufit_inputs/lineages.txt'
    
    Option to save occurr_freqs
    #occurr_freqs['mean_abundance'].to_csv('experimental_outputs/Hutch_combined_WOL.tsv', sep = '\t')
    '''