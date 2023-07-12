'''
Helper functions used for pre-processing
and disease type splitting of biom tables and 
metadata. 
'''

import pandas as pd
from biom import load_table
from qiime2 import Artifact
import os

def data_split_helper(biom, meta, fn):
    'Given metadata + biom, create all datatypes needed'
    
    #Create custom df of biom table (just selected disease type)
    table = load_table(biom).to_dataframe()
    df = table[table.columns.intersection(meta['sample_name'].tolist())]
    print('Total number of samples in', fn, ' = ', df.shape[1])
    
    #Create custom qza table (just selected disease type)
    custom_qza = Artifact.import_data("FeatureTable[Frequency]", df.T)

    #Export everything for future use
    export_filepath = os.getcwd() + '/processed_data/'
    
    #Export metadata
    meta_filename = export_filepath + 'metadata/' +'metadata_' + fn + '.tsv'
    meta.to_csv(meta_filename, sep = '\t', index = False)
         
    #Export biom table in the form of pandas df 
    filename = export_filepath + 'pandas_df/' + fn + '.tsv'
    df.to_csv(filename, sep = '\t')
    
    #Export qza
    custom_qza.save(export_filepath + 'qza/' + fn + '.qza')
    
    #Save qza as biom table
    Artifact.export_data(custom_qza, export_filepath + 'biom/' + fn + '.biom')
    
    return()

def data_split_helper_csv(meta, progression, timepoint, csv, fn):
    'Given metadata + csv, create other files - Particular to Sam input'
    
    #Convert filenames from Sam's file '23341' into filenames from ours '14598.PCGA.1123.BE.01'
    sam_df = pd.read_csv(csv)
    
    #Drop all columns which have normal samples (since these are blood)
    #Along with a few others that weren't in our data
    columns_to_drop = list(sam_df.columns[sam_df.columns.str.contains('N')])
    columns_to_drop = columns_to_drop + ['23761', '24858', '24879']
    sam_df = sam_df.drop(columns_to_drop, axis=1)
    
    #Switch column names to match those in our other data
    converted_cols=['']
    for orginal_sample_name in sam_df.columns.tolist()[1:]:
        sample_name = meta.loc[meta['orginal_sampleid'] == int(orginal_sample_name), 'sample_name'].values[0]
        converted_cols.append(sample_name)
    sam_df.columns = converted_cols
        
    #Create metadata
    meta_custom =  meta[(meta.sampletype.isin(timepoint)) & 
                        (meta.progressionstatus.isin(progression))]
    meta_custom = meta_custom.reset_index(drop = True)
    display(meta_custom[:3])
    
    #Subset to just selected disease type
    sam_df = sam_df.set_index(sam_df.columns[0])
    df = sam_df[sam_df.columns.intersection(meta_custom['sample_name'].tolist())]
    print('Total number of samples in', fn, ' = ', df.shape[1])
    
    #Since we have fractions here, multiply every cell by 100,000
    df = df * 100000
    
    #Create custom qza table (just selected disease type)
    custom_qza = Artifact.import_data("FeatureTable[Frequency]", df.T)

    #Export everything for future use
    export_filepath = os.getcwd() + '/processed_data/'
    
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
    
    
def Sam_BE_data_split(progression, timepoint, csv, fn):
    '''To split csv table from BE Samples from Sam using MetaPhlan4 '''
    
    #Import metadata from equivalent study on Qiita
    meta = pd.read_csv(os.getcwd() + '/qiita_downloads/qiita14598_BE_Esoph/sample_information_from_prep_14498.tsv', sep = '\t')
    
    #Convert from Sam to our format all files based on imported csv/meta_custom
    data_split_helper_csv(meta, progression, timepoint, csv, fn) 
    
    
def be_data_split(progression, timepoint, biom, fn):
    '''To split biom table from 14598, into various timepoint/progression types'''
    
    #Import metadata from Qiita
    meta = pd.read_csv(os.getcwd() + '/qiita_downloads/qiita14598_BE_Esoph/sample_information_from_prep_14498.tsv', sep = '\t')
    
    #Create metadata
    meta_custom =  meta[(meta.sampletype.isin(timepoint)) & 
                        (meta.progressionstatus.isin(progression))]
    meta_custom = meta_custom.reset_index(drop = True)
    meta_custom['sampletype_progressionstatus'] = meta_custom['sampletype'] + '_' + meta_custom['progressionstatus']
    display(meta_custom[:3])
    
    #Create all files based on biom/meta_custom
    data_split_helper(biom, meta_custom, fn)
    
    return()
    

def gerd_normal_data_split(disease_type, biom, fn):
    '''To split biom table from 14458, into a normal and gerd biom table'''
    
    #Import metadata from Qiita
    meta = pd.read_csv(os.getcwd() + '/qiita_downloads/qiita14458_NormalGerdEsoph/sample_information_from_prep_14487.tsv', sep = '\t')
    
    #Create custom metadata for disease of intrest (just selected disease type)
    meta_custom =  meta[(meta.diagnosis == disease_type)]
    meta_custom = meta_custom.reset_index(drop = True)
    display(meta_custom[:3])
    
    #Create all files based on biom/meta_custom
    data_split_helper(biom, meta_custom, fn)
    
    return()

def eac_data_split(biom, fn): 
    
    #Import metadata from Qiita
    meta = pd.read_csv(os.getcwd() + '/qiita_downloads/qiita14521_EAC_ICGC/sample_info_14521_20220509-203030.txt', sep = '\t')
    
    #Create all files based on biom/meta_custom
    data_split_helper(biom, meta, fn)
    display(meta[:3])
   
    return()

def normal_be_eac_data_split(disease_type, biom, fn):
    '''To split biom table from TBD, into a normal, BE, EAC biom table'''
    
    #Import metadata from Qiita
    meta = pd.read_csv(os.getcwd() + '/qiita_downloads/qiitaTBD_Norm_BE_EAC/qiita_sample_info.txt', sep = '\t')
    
    #Remove all the samples from the 'AHM1051'
    meta = meta[(meta.Patient != 'AHM1051')]
    
    #Create custom metadata for disease of intrest (just selected disease type)
    meta_custom = meta[(meta.Sample_type.isin(disease_type))]
    meta_custom = meta_custom.reset_index(drop = True)
    display(meta_custom[:3])
    
    #Create all files based on biom/meta_custom
    data_split_helper(biom, meta_custom, fn)
    
    return()

def HCC_data_split(tumor_type, host_sample_type, biom, fn):
    '''Used for HCC Amir Aim 3'''
    
    #Import metadata from Qiita
    meta = pd.read_csv(os.getcwd() + '/qiita_downloads/Aim3/qiita13756_HCC/qiita_13756_short_wblanks.csv', sep = ',')
    
    #Create metadata
    meta_custom =  meta[(meta.tumor_type.isin(tumor_type)) & 
                        (meta.host_sample_type.isin(host_sample_type))]
    meta_custom = meta_custom.reset_index(drop = True)
    meta_custom['sample_name'] = '13756.' + meta_custom['sample_name']
    display(meta_custom[:3])
    
    #Create all files based on biom/meta_custom
    data_split_helper(biom, meta_custom, fn)
    
    return()

def TCGA_subset_data_split(disease_type, biom, fn):
    '''Used for TCGA subset with just colon/liver Aim 3'''
    
    #Import metadata from Qiita
    meta = pd.read_csv(os.getcwd() + '/qiita_downloads/Aim3/qiita15027_TCGA_colonLiver/sample_information_from_prep_14766.tsv', sep = '\t')
    
    #Create metadata
    meta_custom =  meta[(meta.disease_type.isin(disease_type))]
    meta_custom = meta_custom.reset_index(drop = True)
    display(meta_custom[:3])

    #Create all files based on biom/meta_custom
    data_split_helper(biom, meta_custom, fn)
    
    return()