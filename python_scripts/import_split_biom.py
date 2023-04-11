'''
Helper functions used for pre-processing
and disease type splitting of biom tables and 
metadata. 
'''

import pandas as pd
from biom import load_table
from qiime2 import Artifact

def gerd_normal_data_split(disease_type, biom, fn):
    '''To split biom table from 14458, into a normal and biom table'''
    
    user_path = '/Users/cguccion/Dropbox/current_Spring23/ISMEJ_EsophModeling_Paper/Figures/comad_of_EAC_progression/'
    
    #Import metadata from Qiita
    meta = pd.read_csv(user_path + 'qiita_downloads/qiita14458_NormalGerdEsoph/sample_information_from_prep_14487.tsv', sep = '\t')
    
    #Create custom metadata for disease of intrest (just selected disease type)
    meta_custom =  meta[(meta.diagnosis == disease_type)]
    meta_custom = meta_custom.reset_index(drop = True)
    display(meta_custom[:3])
    
    #Create custom df of biom table (just selected disease type)
    table = load_table(biom).to_dataframe()
    df = table[table.columns.intersection(meta_custom['sample_name'].tolist())]
    print('Total number of samples in', fn, ' = ', df.shape[1])
    
    #Create custom qza table (just selected disease type)
    custom_qza = Artifact.import_data("FeatureTable[Frequency]", table.T)

    #Export everything for future use
    export_filepath = user_path + 'processed_data/'
    
    #Export metadata
    meta_filename = export_filepath + 'metadata/' +'metadata_' + fn + '.tsv'
    meta_custom.to_csv(meta_filename, sep = '\t', index = False)
         
    #Export biom table in the form of pandas df 
    filename = export_filepath + 'pandas_df/' + fn + '.tsv'
    df.to_csv(filename, sep = '\t')
    
    #Export qza
    custom_qza.save(export_filepath + 'qza/' + fn + '.qza')
    
    #Save qza as biom table
    Artifact.export_data(custom_qza, export_filepath + 'biom/' + fn + '.biom')
    
    return()