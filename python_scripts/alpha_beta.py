import pandas as pd
import qiime2 as q2
from qiime2 import Artifact, Metadata, Visualization
from biom import load_table
import seaborn as sns
import matplotlib.pyplot as plt
import os
import subprocess
import upsetplot
from matplotlib_venn import venn2
from scipy.stats import mannwhitneyu
import itertools
from scipy.stats import kruskal
import numpy as np
import scikit_posthocs as sp

from qiime2.plugins.taxa.visualizers import barplot
from qiime2.plugins.diversity.pipelines import alpha
from qiime2.plugins.diversity.visualizers import alpha_group_significance
from qiime2.plugins.diversity.visualizers import alpha_rarefaction
from qiime2.plugins.feature_table.methods import merge
from qiime2.plugins.emperor.visualizers import biplot
from qiime2.plugins import diversity
from qiime2.plugins.diversity_lib.methods import weighted_unifrac

#Color Palette
color_palette = ['#7e2954', '#ED6677', '#F26122', '#f0e442', '#56b3e9', '#009e74', '#483fa1', '#dddddd']

#Taxonomy
wol_taxonomy = Artifact.import_data('FeatureData[Taxonomy]', '/Users/cguccion/Dropbox/Storage/HelpfulLabDocs/taxonomy_trees/WOL/lineages_nonSpace.txt','HeaderlessTSVTaxonomyFormat')

def combine_qza_meta(fn_list, combo_meta_col=['sample_name', 'study_label']):
    '''Combine all imported filenames so they can be used in group analysis '''
    
    ft = Artifact.load(os.getcwd() + '/processed_data/qza/' + fn_list[0] + '.qza')
    
    meta = pd.read_csv(os.getcwd() + '/processed_data/metadata/metadata_' + fn_list[0] + '.tsv', sep='\t', index_col=0)
    meta['study_label'] = fn_list[0]
    
    count=0
    for fn in fn_list[1:]:
        
        #Merge taxonomy tables
        current_ft = Artifact.load(os.getcwd() + '/processed_data/qza/' + fn + '.qza')
        ft, = merge([ft, current_ft])
        
        #Merge metadata
        current_meta = pd.read_csv(os.getcwd() + '/processed_data/metadata/metadata_' + fn + '.tsv', sep='\t', index_col=0)
        current_meta['study_label'] = fn
        meta = pd.merge(meta, current_meta, how = 'outer', on= combo_meta_col, suffixes=('_' + str(count), '_' + str(count+1)))
    
        count+=1
        
    meta = meta[['study_label']]
    
    meta = q2.Metadata(meta)
    
    return(ft, meta)

def multi_alpha_rare(fn_list, name, combo_meta_col =['sample_name', 'study_label'], max_depth=200000):
    '''Alpha Rarefaction Graphs - allows us to see which rarefaction to chose'''
    
    ft, meta = combine_qza_meta(fn_list, combo_meta_col=combo_meta_col)
    
    alpha_rare = alpha_rarefaction(table= ft, max_depth= max_depth, metadata= meta)
    alpha_rare_v = alpha_rare.visualization
    alpha_rare_v.save(os.getcwd() + '/outputs/alpha_rare/alphaRareCurve_' + name + '.qzv')
    
    return(alpha_rare_v)

def multi_taxonomy_barplot(fn_list, name, taxonomy=wol_taxonomy, combo_meta_col =['sample_name', 'study_label']):
    '''Create taxonomy column barpolots '''
    
    ft, meta = combine_qza_meta(fn_list, combo_meta_col=combo_meta_col)
    
    barplot_output = barplot(table=ft, taxonomy=taxonomy, metadata=meta)
    barplot_output_v = barplot_output.visualization
    barplot_output_v.save(os.getcwd() + '/outputs/taxonomy_barplots/barplot_' + name + '.qzv')
    
    return(barplot_output_v)
        
def multi_basic_alpha(fn_list, name, metric = 'observed_features', combo_meta_col =['sample_name', 'study_label']):
    '''Shannon and Observed Alpha Diversity Graphs - creates the actual qza '''
    
    ft, meta = combine_qza_meta(fn_list, combo_meta_col=combo_meta_col)
    
    alpha_vector, = alpha(table= ft, metric=metric)
    alpha_output = alpha_group_significance(alpha_diversity = alpha_vector, metadata = meta)
    alpha_output_v = alpha_output.visualization
    alpha_output_v.save(os.getcwd() + '/outputs/alpha_plots/alpha_' + name + '_' + metric + '.qzv')
    
    return(alpha_output)

def alpha_plot(name, metric='observed_features'):
    '''Acutal pretty printing of Shannon and Observed Alpha Diversity Graphs'''
    
    df = pd.read_csv(os.getcwd() + '/outputs/alpha_plots/alpha_' + name + '_' + metric + '/data/metadata.tsv',
                     sep='\t', skiprows=[1])
    
    #Not sure if will use this, but all the p-values
    p_df = pd.read_csv(os.getcwd() + '/outputs/alpha_plots/alpha_' + name + '_' + metric + '/data/kruskal-wallis-pairwise-study_label.csv')
    
    sns.set_theme(style="white")
    
    if 'shannon' in metric:
        metric='shannon_entropy'
    
    sns.set_palette(color_palette) #d55e00 e69f00 
    plt.figure(figsize=(20, 6))
    # Box plot
    sns.boxplot(x='study_label', y=metric, data=df, showfliers=False)  
        
    # Strip plot
    sns.stripplot(x='study_label', y=metric, data=df, color='black', size=4, jitter=True)  
   
    plt.tight_layout()
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
    plt.rcParams['svg.fonttype'] = 'none'  # render SVG text as text, not curves
    plt.xticks(rotation=45, ha='right') 
    
    #Save in svg format
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
    plt.rcParams['svg.fonttype'] = 'none'  # render SVG text as text, not curves
    plt.savefig(os.getcwd() + '/outputs/alpha_plots/' + name + '_alpha_' + metric + '.svg')
    
    plt.show()
    
def alpha_stats(name, metric='observed_features'):
    '''Stats used in the publication for Alpha Diversity'''
    
    df = pd.read_csv(os.getcwd() + '/outputs/alpha_plots/alpha_' + name + '_' + metric + '/data/metadata.tsv',
                     sep='\t', skiprows=[1])
    
    display(df)
    

def beta_decoide(fn_list, name, metric='study_label', permutations=9999):
    
    '''
    Calculates beta diversity for unweighted 
    & weighted unifraq

    Parameters
    ---------
    ft: q2 FeatureTable[Frequency] object
        Qiime2 taxonomy object
    
    meta: pd df
        metadata
        
    metric: str
        The column in sample_meta used to
        calculate beta diversity across 
        samples
        
    return_biplot: bool (False)
        Returns the biplot instead 
        of the permaonva ouputs 
    
    Returns
    -------
    beta_result_d: q2 Visulization object
        Qiime2 visulaization object of 
        decoide
    
    Notes
    -----
    '''
    
    ft, meta = combine_qza_meta(fn_list)
    
    #Save ft for command line processing
    ft.save('deicode_processing.qza')
    
    #Run deicode from command line using qza created above
    ##!qiime deicode rpca --i-table deicode_processing.qza --o-biplot deicode_biplot.qza --o-distance-matrix deicode_distance_test.qza
    command = [
        'qiime',
        'deicode',
        'rpca',
        '--i-table',
        'deicode_processing.qza',
        '--o-biplot',
        'deicode_biplot.qza',
        '--o-distance-matrix',
        'deicode_distance_test.qza']
    
    # Run the command
    subprocess.run(command)
    
    #Import biplot back into python
    rpca_biplot = Artifact.load('deicode_biplot.qza')
    
    #Import biplot back into python
    rpca_distance_matrix = Artifact.load('deicode_distance_test.qza')
    
    #Create emperor visualization from the biplot result
    rpca_biplot_emperor = biplot(biplot = rpca_biplot, sample_metadata = meta)
    
    rpca_biplot_emperor_v = rpca_biplot_emperor.visualization
    rpca_biplot_emperor_v.save(os.getcwd() + '/outputs/RPCA/RPCA_' + name + '.qzv')
    
    #Calculate permanova 
    beta_result_d = diversity.actions.beta_group_significance(distance_matrix=rpca_distance_matrix,
                                                              metadata=meta.get_column(metric),
                                                              method = 'permanova', pairwise = True,
                                                              permutations=permutations)
    
    
    beta_result_v = beta_result_d.visualization
    beta_result_v.save(os.getcwd() + '/outputs/RPCA/permanova_' + name + '.qzv')
    
    return(rpca_biplot_emperor_v, beta_result_v)


def venn_digram_plot_(fn_list, name):
    
    fn_taxa_dict={}
    correct_fn={'normal': 'Healthy    ', 'GERD': 'GERD    ', 'BE_NP_T1':'BE NCO    ',
               'BE_Prog_T1': 'BE CO    ', 'EAC_Ross':'EAC cohort 2    ', 'EAC_ICGC':'EAC cohort 1    '}
    color_fn={'normal': '#743555', 'GERD': '#dc7784', 'BE_NP_T1':'#d96c3c',
               'BE_Prog_T1': '#dad15a', 'EAC_Ross':'#168a6b', 'EAC_ICGC':'#68add7'}
    
    for fn in fn_list:
        current_ft = Artifact.load(os.getcwd() + '/processed_data/qza/' + fn + '.qza')
        current_taxa = set(current_ft.view(pd.DataFrame).columns)
        fn = fn.split('_species')[0]
        fn_taxa_dict[fn] = current_taxa
    
    # Create all pairs of comparisons
    keys = list(fn_taxa_dict.keys())
    pairs = [(keys[i], keys[j]) for i in range(len(keys)) for j in range(i+1, len(keys))]

    # Plot
    plt.figure(figsize=(7, 7))
    for idx, (key1, key2) in enumerate(pairs, 1):
        plt.subplot(5, 3, idx)  # Adjust rows and columns based on the number of pairs
        venn = venn2([fn_taxa_dict[key1], fn_taxa_dict[key2]], set_labels=(correct_fn[key1], correct_fn[key2]))

        # Set colors and transparency
        venn.get_patch_by_id('11').set_color('#d4dcdd')
        venn.get_patch_by_id('10').set_color(color_fn[key1])
        venn.get_patch_by_id('01').set_color(color_fn[key2])

        # Set font sizes and line style
        for text in venn.set_labels:
            if text:
                text.set_fontsize(7)
                text.set_fontfamily('Arial')
        for subset in venn.subset_labels:
            if subset:
                subset.set_fontsize(6)
                subset.set_fontfamily('Arial')

    plt.rcParams['svg.fonttype'] = 'none' 
    plt.tight_layout()
    plt.savefig(os.getcwd() + '/outputs/venn_diagrams/venn_' + name + '.svg')
    plt.show()
    
def create_upset_(fn_list, name):
    
    fn_taxa_dict={}
    correct_fn={'normal': 'Healthy', 'GERD': 'GERD', 'BE_NP_T1':'BE NCO',
               'BE_Prog_T1': 'BE CO', 'EAC_Ross':'EAC early stage', 'EAC_ICGC':'EAC multi stage'}
    
    for fn in fn_list:
        current_ft = Artifact.load(os.getcwd() + '/processed_data/qza/' + fn + '.qza')
        current_taxa = set(current_ft.view(pd.DataFrame).columns)
        fn = fn.split('_species')[0]
        fn_taxa_dict[correct_fn[fn]] = current_taxa
        
    all_elems = fn_taxa_dict['Healthy'].union(fn_taxa_dict['GERD']).union(fn_taxa_dict['BE NCO'].union(fn_taxa_dict['BE CO']).union(fn_taxa_dict['EAC early stage']).union(fn_taxa_dict['EAC multi stage']))
    df = pd.DataFrame([[e in fn_taxa_dict['Healthy'], e in fn_taxa_dict['GERD'],
                        e in fn_taxa_dict['BE NCO'], e in fn_taxa_dict['BE CO'],
                       e in fn_taxa_dict['EAC early stage'], e in fn_taxa_dict['EAC multi stage']] for e in all_elems], columns = ['Healthy', 'GERD', 'BE NCO', 'BE CO', 'EAC early stage', 'EAC multi stage'])
    df_up = df.groupby(['Healthy', 'GERD', 'BE NCO', 'BE CO', 'EAC early stage', 'EAC multi stage']).size()
    
    upsetplot.plot(df_up, orientation='horizontal')    
    
    #Leave these high because shrink down figure in illustrater
    plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 18,  # Adjust font size as needed
    'axes.labelsize': 22,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18
    })
    plt.rcParams['svg.fonttype'] = 'none'  # render SVG text as text, not curves
    
    plt.savefig(os.getcwd() + '/outputs/venn_diagrams/upset_' + name + '.svg')
    plt.show()
    
def alpha_mannWhit(name, metric='observed_features'):
    '''Stats used for Mann Whitney'''
    
    df = pd.read_csv(os.getcwd() + '/outputs/alpha_plots/alpha_' + name + '_' + metric + '/data/metadata.tsv',
                     sep='\t', skiprows=[1])
    
    if metric == 'shannon':
        metric = 'shannon_entropy'
    
    #Calculate Mann Whitney U p-values across lables
    # Get the unique labels in the 'study_label' column
    study_labels = df['study_label'].unique()

    # Store the results
    results = []

    # Iterate over all pairwise combinations of study labels
    for label1, label2 in itertools.combinations(study_labels, 2):
        # Subset the data for each pair of study labels
        group1 = df[df['study_label'] == label1][metric]
        group2 = df[df['study_label'] == label2][metric]

        # Perform Mann-Whitney U test
        u_stat, p_value = mannwhitneyu(group1, group2, alternative='two-sided')

        # Append the result to the results list
        results.append({
            'Group 1': f"{label1.split('_species')[0]} (n= {len(df[df['study_label'] == label1])})",
            'Group 2': f"{label2.split('_species')[0]} (n= {len(df[df['study_label'] == label2])})",
            'u_statistic (' + metric + ')': u_stat,
            'p_value (' + metric + ')': p_value
        })

    # Convert the results into a DataFrame
    results_df = pd.DataFrame(results)

    # Display the results    
    display(results_df)

def alpha_KW(name, metric='observed_features'):
    '''Stats used in the publication for Kruskal'''
    
    df = pd.read_csv(os.getcwd() + '/outputs/alpha_plots/alpha_' + name + '_' + metric + '/data/metadata.tsv',
                     sep='\t', skiprows=[1])
    
    if metric == 'shannon':
        metric = 'shannon_entropy'
    
    # Group the data by study_label and get the shannon_entropy values for each group
    groups = [group[metric].values for _, group in df.groupby('study_label')]

    # Perform the Kruskal-Wallis test
    h_stat, p_value = kruskal(*groups)
    
    print(metric)
    print("Kruskal-Wallis H-statistic:", h_stat)
    print("p-value:", p_value)
    print()
    

def format_values(x):
    if isinstance(x, (float, int)):
        if x > 0.001:
            return f"{x:.4f}"  # Show 4 decimal places
        else:
            return f"{x:.1e}"  # Use scientific notation
    return x 

def alpha_KW_dunn(name, metric='observed_features'):
    '''Stats used in the publication for Kruskal with Dunn's test added'''
    
    df = pd.read_csv(os.getcwd() + '/outputs/alpha_plots/alpha_' + name + '_' + metric + '/data/metadata.tsv',
                     sep='\t', skiprows=[1])
    
    if metric == 'shannon':
        metric = 'shannon_entropy'
    
    # Group the data by study_label and get the shannon_entropy values for each group
    groups = [group[metric].values for _, group in df.groupby('study_label')]

    # Perform the Kruskal-Wallis test
    h_stat, p_value = kruskal(*groups)
    
    print(metric)
    print("Kruskal-Wallis H-statistic:", h_stat)
    print("p-value:", p_value)
    print()
    
    #Dunns PostHoc Calculations
    posthoc_results = sp.posthoc_dunn(df, val_col=metric, group_col='study_label', p_adjust='holm')
    posthoc_results = posthoc_results.applymap(format_values)
    
    #Print Dunn Result
    already_print=[]
    for row, col_data in posthoc_results.iterrows():
        for col, value in col_data.items():
            if pd.notnull(value):  # Ignore NaN values
                if row != col:
                    if str(row + col) not in already_print:
                        if float(value) < 0.05:
                            print(f"{row.split('_species')[0]} vs {col.split('_species')[0]}: {value}")
                        else:
                             print(f"{row.split('_species')[0]} vs {col.split('_species')[0]}: {value}" + "- NOT SIGNFIGANT")
                        already_print.append(str(col + row))

    display(posthoc_results)

    
def calc_med_alpha(name, metric='observed_features'):
    '''Stats used in the publication for median and Q1/Q3'''
    
    df = pd.read_csv(os.getcwd() + '/outputs/alpha_plots/alpha_' + name + '_' + metric + '/data/metadata.tsv',
                     sep='\t', skiprows=[1])
    
    if metric == 'shannon':
        metric = 'shannon_entropy'
    
    summary_df = df.groupby("study_label")[metric].agg(
        Median="median",
        Q1=lambda x: np.percentile(x, 25),
        Q3=lambda x: np.percentile(x, 75)
    ).reset_index()

    # Display the resulting DataFrame
    print(metric)
    display(summary_df)
    

def alphaPaperStats(name):
    '''Stats used in the publication for Alpha Diversity'''
    alpha_stats(name)
    #alpha_KW(name)
    calc_med_alpha(name)
    alpha_mannWhit(name)
    alpha_KW_dunn(name)
    print('------+++++++------+++++++------+++++++------+++++++------+++++++------+++++++------+++++++')
    alpha_stats(name, metric='shannon')
    #alpha_KW(name,  metric='shannon')
    calc_med_alpha(name,  metric='shannon')
    alpha_mannWhit(name, metric='shannon')
    alpha_KW_dunn(name, metric='shannon')