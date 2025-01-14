'''
Helper functions used for pre-processing
and disease type splitting of biom tables and 
metadata. 
'''

import pandas as pd
from biom import load_table
from qiime2 import Artifact
from qiime2.plugins import taxa
from qiime2.plugins.feature_table.methods import rarefy
import os

#This is the list of icgc_donor_ids that are in the Ross-Innes data. These are the ones that specifically need to be removed in ICGC
Ginny_donor = ['DO50326', 'DO234426', 'DO50436', 'DO50362', 'DO50381', 'DO234279', 'DO50331', 'DO234371', 'DO50318', 'DO50409',
              'DO234161', 'DO234232', 'DO234344', 'DO234406', 'DO50383', 'DO234226', 'DO50384', 'DO234408', 'DO50345', 
              'DO50445', 'DO50448', 'DO50374', 'DO50319']

#This allows for the subsetting of JUST the samples in the physical figures of the paper (not all extra ones we have that are listed in the supplemental) 
#This is used to subset the Ross-Innes paper 
Ginny_oac_bo = ['LP6007591', 'LP6007592', 'LP6005690-DNA_H01', 'LP6005691-DNA_H01', 'LP6005500-DNA_D01', 'LP6005501-DNA_D01', 'LP6005690-DNA_D01', 'LP6005691-DNA_E01', 'LP6005690-DNA_A02', 'LP6005691-DNA_C01', 'LP6005500-DNA_C01', 'LP6005501-DNA_C01', 'LP6007594', 'LP6007595', 'LP6005690-DNA_G01', 'LP6005691-DNA_G01', 'LP6007404-DNA_A01', 'LP6007405-DNA_A01', 'LP6007520-DNA_A01', 'LP6007521-DNA_A01', 'LP6005334-DNA_A03', 'LP6005335-DNA_A01', 'LP6007409-DNA_A01', 'LP6007410-DNA_A01', 'LP6005690-DNA_F01', 'LP6005691-DNA_F01', 'LP6005500-DNA_E01', 'LP6005501-DNA_E01', 'LP6005690-DNA_C01', 'LP6005691-DNA_A01', 'LP6005690-DNA_F03', 'LP6005691-DNA_A02', 'LP6005690-DNA_B02', 'LP6005691-DNA_D01', 'LP6005500-DNA_F01', 'LP6005501-DNA_F01', 'LP2000104-DNA_A01', 'LP2000110-DNA_A01', 'LP6005500-DNA_A01', 'LP6005501-DNA_A01', 'LP6005500-DNA_B01', 'LP6005501-DNA_B01', 'LP6005690-DNA_E01', 'LP6005691-DNA_B01', 'LP6007597', 'LP6007598']

def remove_reagents(df):
    bad_reagents = pd.read_csv('reagent_lab_contam_list.tsv', delimiter='\t')
    bad_r_list = bad_reagents['genus'].tolist()[1:]
    
    # Create a boolean mask using str.contains
    mask = df.index.to_series().apply(lambda x: any(substring in x for substring in bad_r_list))

    # Filter the DataFrame to exclude rows containing the strings
    df_new = df[~mask]

    return(df_new)

def keep_human(df):
    human_asso = pd.read_csv('human_asso_microbes.txt', delimiter='\t')
    human_asso = human_asso['species'].tolist()[1:]
    
    # Create a boolean mask using str.contains
    mask = df.index.to_series().apply(lambda x: any(substring in x for substring in human_asso))

    # Filter the DataFrame to include rows containing the strings
    df_new = df[mask]

    return(df_new)    

def data_split_helper(biom, meta, fn, micov_filter = False, zebra=False, rare_level=2000, zebra_threshold=0.01, zebra_rare=False, remove_bad_microbes=False, remove_reag=False, keep_hu=False, remove_bad_gotu=False, remove_bad_gotu_micov=False):
    'Given metadata + biom, create all datatypes needed'
    
    #Create custom df of biom table (just selected disease type)
    table = load_table(biom).to_dataframe()
    df = table[table.columns.intersection(meta['sample_name'].tolist())]
    
    if remove_bad_microbes == True:
        #Remove 'bad' microbes
        df = df.drop(bad_microbes_list)
        
    if remove_bad_gotu == True:
        #Remove 'bad' gotus that only apppear in high abudances in second ICGC run
        if 'species' in fn:
            df = df.drop(bad_species, errors='ignore')
        else:
            df = df.drop(bad_gotu)
            
    if remove_bad_gotu_micov == True:
        #Remove 'bad' gotus that only apppear in high abudances in second ICGC run + MICOV
        if 'species' in fn:
            df = df.drop(bad_species_micov, errors='ignore')
        else:
            df = df.drop(bad_gotu_micov)
            
    if micov_filter != False:
            micov_list_path = os.getcwd() + '/micov_results/micov_filtered_lists/'
            
            pass_= list( pd.read_csv(micov_list_path + 'PASSmicov_' + micov_filter + '.csv', sep='\t')['Pass'])
            fail_= list( pd.read_csv(micov_list_path + 'FAILmicov_' + micov_filter + '.csv', sep='\t')['Fail'])
            
            for i in df.index:
                
                #Species level
                if ';' in i:
                    taxdf = pd.read_csv(os.getcwd() + '/qiita_downloads/WOL_lineages.txt', sep='\t', names=['GOTU', 'Taxa'])
                    pass_species = list(taxdf[taxdf['GOTU'].isin(pass_)]['Taxa'])
                    pass_species = list(map(lambda s: s.replace('; ', ';'), pass_species))
                    fail_species = list(taxdf[taxdf['GOTU'].isin(fail_)]['Taxa'])
                    fail_species = list(map(lambda s: s.replace('; ', ';'), fail_species))
                #GOTU level
                else:
                    pass_species = pass_.copy()
                    fail_species = fail_.copy()

                if i not in pass_species and i not in fail_species:
                    print('WARNING: This GOTU had error durring micov and will not be dropped: ', i)

            df = df.drop(fail_species)
    
    print('Total number of samples in', fn, ' = ', df.shape[1])
    
    #Removes any microbes in the reagents list
    if remove_reag == True:
        df = remove_reagents(df)
    
    #Keeps only microbes which are in the human assoicated list
    if keep_hu == True:
        df = keep_human(df)
        
    display(df)    
        
    #Create custom qza table (just selected disease type)
    custom_qza = Artifact.import_data("FeatureTable[Frequency]", df.T)
    
    #Create rarified table off qza table
    rare_qza, = rarefy(table=custom_qza, sampling_depth = rare_level)
    print(rare_qza)
    
    #Convert rare table into pandas df
    rare_df = rare_qza.view(pd.DataFrame)

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
    
    print('DF len', len(df))
        
    ###------------ Save Rare Info ###------------
    
    #Export metadata -- I know this is the same but easier if they all follow a pattern
    meta_filename = export_filepath + 'metadata/' +'metadata_' + fn + '_r' + str(rare_level) + '.tsv'
    meta.to_csv(meta_filename, sep = '\t', index = False)
    
    #(RARE) Export biom table in the form of pandas df
    filename_r = export_filepath + 'pandas_df/' + fn + '_r' + str(rare_level) + '.tsv'
    rare_df.T.to_csv(filename_r, sep = '\t')
    
    #(RARE) Export qza
    rare_qza.save(export_filepath + 'qza/' + fn + '_r' + str(rare_level) + '.qza')
    
    #Save qza as biom table
    Artifact.export_data(rare_qza, export_filepath + 'biom/' + fn + '_r' + str(rare_level) + '.biom')
    
    if zebra_rare != False:
        #Convert from genome to species and output that as well
        t_fn = fn + '_r' + str(rare_level)
        species_from_zebra(rare_qza, export_filepath, t_fn, meta)
    
    return()

def data_split_helper_csv(meta, progression, timepoint, csv, fn , counts=False, rare_level=800):
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
    
    if counts == False:
        #Since we have fractions here, multiply every cell by 100,000
        df = df * 100000
    
    #Create custom qza table (just selected disease type)
    custom_qza = Artifact.import_data("FeatureTable[Frequency]", df.T)
    
    #Create rarified table off qza table
    rare_qza, = rarefy(table=custom_qza, sampling_depth = rare_level)
    print(rare_qza)
    
    #Convert rare table into pandas df
    rare_df = rare_qza.view(pd.DataFrame)

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
    
    print('DF len', len(df))
        
    ###------------ Save Rare Info ###------------
    
    #Export metadata -- I know this is the same but easier if they all follow a pattern
    meta_filename = export_filepath + 'metadata/' +'metadata_' + fn + '_r' + str(rare_level) + '.tsv'
    meta_custom.to_csv(meta_filename, sep = '\t', index = False)
    
    #(RARE) Export biom table in the form of pandas df
    filename_r = export_filepath + 'pandas_df/' + fn + '_r' + str(rare_level) + '.tsv'
    rare_df.T.to_csv(filename_r, sep = '\t')
    
    #(RARE) Export qza
    rare_qza.save(export_filepath + 'qza/' + fn + '_r' + str(rare_level) + '.qza')
    
    #Save qza as biom table
    Artifact.export_data(rare_qza, export_filepath + 'biom/' + fn + '_r' + str(rare_level) + '.biom')
    
    return()
    
    
def Sam_BE_data_split(progression, timepoint, csv, fn, counts=False, rare_level=800):
    '''To split csv table from BE Samples from Sam using MetaPhlan4 '''
    #Leave counts=False if using relative abudances data, else change to counts=True
    
    #Import metadata from equivalent study on Qiita
    meta = pd.read_csv(os.getcwd() + '/qiita_downloads/qiita14598_BE_Esoph/sample_information_from_prep_14498.tsv', sep = '\t')
    
    #Convert from Sam to our format all files based on imported csv/meta_custom
    data_split_helper_csv(meta, progression, timepoint, csv, fn, counts=counts, rare_level=rare_level) 
    
    
def be_data_split(progression, timepoint, biom, fn, micov_filter=False, zebra=False, exact_pairs= False, trim=False, remove_reag=False, keep_hu=False, remove_bad_gotu=False, rare_level=2000, sex_split=False):
    '''To split biom table from 14598, into various timepoint/progression types'''
    
    #Import metadata from Qiita
    meta = pd.read_csv(os.getcwd() + '/qiita_downloads/qiita14598_BE_Esoph/sample_information_from_prep_14498.tsv', sep = '\t')
    
    #Create metadata
    meta_custom =  meta[(meta.sampletype.isin(timepoint)) & 
                        (meta.progressionstatus.isin(progression))]
    meta_custom = meta_custom.reset_index(drop = True)
    meta_custom['sampletype_progressionstatus'] = meta_custom['sampletype'] + '_' + meta_custom['progressionstatus']
    
    if exact_pairs == True:
        new_meta = pd.DataFrame(columns = meta_custom.columns)
        for index, row in meta_custom.iterrows():
            if any((new_meta['pairednormal'] == row['pairednormal']) & (new_meta['sampletype'] == row['sampletype'])) == False:
                new_meta.loc[len(new_meta)] = row
        meta_custom = new_meta
    
    if trim!= False:
        meta_custom = meta_custom.sample(trim)
    display(meta_custom[:3])
    
    #Allows M/F split
    if sex_split != False:
        meta_custom = meta_custom[(meta_custom.gender.isin(sex_split))]
    
    #Create all files based on biom/meta_custom
    data_split_helper(biom, meta_custom, fn, micov_filter=micov_filter, zebra=zebra, remove_reag=remove_reag, keep_hu=keep_hu, remove_bad_gotu=remove_bad_gotu, rare_level=rare_level)
    
    return()
    

def gerd_normal_data_split(disease_type, biom, fn, micov_filter = False, zebra=False, trim=False, shotgun=True, remove_reag=False, keep_hu=False, remove_bad_gotu=False, rare_level=2000, sex_split=False):
    '''To split biom table from 14458, into a normal and gerd biom table'''
    
    #Import metadata from Qiita
    meta = pd.read_csv(os.getcwd() + '/qiita_downloads/qiita14458_NormalGerdEsoph/sample_information_from_prep_14487.tsv', sep = '\t')
    if shotgun == False:
        meta = pd.read_csv(os.getcwd() + '/qiita_downloads/qiita14325_16S_NormalGerdEsoph/matched_w_shotgun_sample_information_from_prep_12316.tsv', sep = '\t')
    
    #Create custom metadata for disease of intrest (just selected disease type)
    meta_custom =  meta[(meta.diagnosis == disease_type)]
    meta_custom = meta_custom.reset_index(drop = True)
    
    #Allows M/F split
    if sex_split != False:
        meta_custom = meta_custom[(meta_custom.gender.isin(sex_split))]
    
    if trim!= False:
        meta_custom = meta_custom.sample(trim)
    display(meta_custom[:3])
    
    #Create all files based on biom/meta_custom
    data_split_helper(biom, meta_custom, fn, micov_filter= micov_filter, zebra=zebra, remove_reag=remove_reag, keep_hu=keep_hu, remove_bad_gotu=remove_bad_gotu, rare_level=rare_level)
    
    return()

def eac_14857_data_split(biom, fn, micov_filter=False, trim=False, remove_bad_microbes=False, Ginny_list=True, remove_reag=False, keep_hu=False, zebra=False, remove_bad_gotu=False, remove_bad_gotu_micov=False, zebra_threshold=0.01, rare_level=2000, sex_split=False): 
    ''' Correct ICGC dataset'''
    
    #Import metadata from Qiita
    meta = pd.read_csv(os.getcwd() + '/qiita_downloads/qiita14857_EAC_ICGC/sample_information_from_ESAD_ICGC_14857_prep_13977.txt', sep = '\t')
    
    #This removes the samples that are in both the ICGC dataset and in Ginny's ... typically we would want these removed
    if Ginny_list == True:
        meta = meta[~meta.isin(Ginny_donor).any(axis=1)]
    
    if trim!= False:
        meta = meta.sample(trim)

    #Allows M/F split
    if sex_split != False:
        meta = meta[(meta.donor_sex.isin(sex_split))]
    
    display(meta[:3])
    
    #Create all files based on biom/meta_custom
    data_split_helper(biom, meta, fn, micov_filter=micov_filter, zebra=zebra, zebra_threshold=zebra_threshold, remove_bad_microbes=remove_bad_microbes, remove_reag=remove_reag, keep_hu=keep_hu, remove_bad_gotu=remove_bad_gotu, remove_bad_gotu_micov=remove_bad_gotu_micov, rare_level=rare_level)
   
    return()

def normal_be_eac_data_split(disease_type, biom, fn, micov_filter = False, exact_pairs=False, zebra=False, remove_bad_microbes=False, Ginny_list=True, remove_reag=False, keep_hu=False, remove_bad_gotu=False, rare_level=2000, sex_split=False):
    '''To split biom table from TBD, into a BE, EAC biom table from Ross-Innes'''
    #exact pairs just takes a single BE and EAC sample from each patient
    
    #Import metadata from Qiita -- We have updated the metadata to more accurate since Maria updated it  
    #meta = pd.read_csv(os.getcwd() + '/qiita_downloads/qiitaTBD_Norm_BE_EAC/qiita_sample_info.txt', sep = '\t')
    meta = pd.read_csv(os.getcwd() + '/compare_rossICGC/corrected_meta.tsv', sep = '\t')
    
    #Remove all the samples from the 'AHM1051'
    meta = meta[(meta.Patient != 'AHM1051')]
    
    #Create custom metadata for disease of intrest (just selected disease type)
    meta_custom = meta[(meta.Sample_type.isin(disease_type))]
    
    #Allows M/F split
    if sex_split != False:
        meta_custom = meta_custom[(meta_custom.Gender.isin(sex_split))]
    
    if exact_pairs == True:
        new_meta = pd.DataFrame(columns = meta_custom.columns)
        for index, row in meta_custom.iterrows():
            if any((new_meta['Patient'] == row['Patient']) & (new_meta['Sample_type'] == row['Sample_type'])) == False:
                new_meta.loc[len(new_meta)] = row
        meta_custom = new_meta
    
    #I think this is subsetting to JUST the samples that are acutally in the Ross-Innes paper ... there was a large amount of extra samples they provided so it removes those
    if Ginny_list == True:
        meta_custom = meta_custom[meta_custom.isin(Ginny_oac_bo).any(axis=1)]
    
    meta_custom = meta_custom.reset_index(drop = True)
    display(meta_custom[:3])
    
    #Create all files based on biom/meta_custom
    data_split_helper(biom, meta_custom, fn, micov_filter=micov_filter, zebra=zebra, remove_bad_microbes=remove_bad_microbes, remove_reag=remove_reag, keep_hu=keep_hu, remove_bad_gotu=remove_bad_gotu, rare_level=rare_level)
    
    return()
