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

Ginny_donor = ['DO50326', 'DO234426', 'DO50436', 'DO50362', 'DO50381', 'DO234279', 'DO50331', 'DO234371', 'DO50318', 'DO50409',
              'DO234161', 'DO234232', 'DO234344', 'DO234406', 'DO50383', 'DO234226', 'DO50384', 'DO234408', 'DO50345', 
              'DO50445', 'DO50448', 'DO50374', 'DO50319']

Ginny_oac_bo = ['LP6007591', 'LP6007592', 'LP6005690-DNA_H01', 'LP6005691-DNA_H01', 'LP6005500-DNA_D01', 'LP6005501-DNA_D01', 'LP6005690-DNA_D01', 'LP6005691-DNA_E01', 'LP6005690-DNA_A02', 'LP6005691-DNA_C01', 'LP6005500-DNA_C01', 'LP6005501-DNA_C01', 'LP6007594', 'LP6007595', 'LP6005690-DNA_G01', 'LP6005691-DNA_G01', 'LP6007404-DNA_A01', 'LP6007405-DNA_A01', 'LP6007520-DNA_A01', 'LP6007521-DNA_A01', 'LP6005334-DNA_A03', 'LP6005335-DNA_A01', 'LP6007409-DNA_A01', 'LP6007410-DNA_A01', 'LP6005690-DNA_F01', 'LP6005691-DNA_F01', 'LP6005500-DNA_E01', 'LP6005501-DNA_E01', 'LP6005690-DNA_C01', 'LP6005691-DNA_A01', 'LP6005690-DNA_F03', 'LP6005691-DNA_A02', 'LP6005690-DNA_B02', 'LP6005691-DNA_D01', 'LP6005500-DNA_F01', 'LP6005501-DNA_F01', 'LP2000104-DNA_A01', 'LP2000110-DNA_A01', 'LP6005500-DNA_A01', 'LP6005501-DNA_A01', 'LP6005500-DNA_B01', 'LP6005501-DNA_B01', 'LP6005690-DNA_E01', 'LP6005691-DNA_B01', 'LP6007597', 'LP6007598']

bad_microbes_list = ['k__Bacteria;p__Bacteroidetes;c__Cytophagia;o__Cytophagales;f__Hymenobacteraceae;g__Hymenobacter;s__Hymenobacter sp. IS2118',
                     'k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Micrococcales;f__Micrococcaceae;g__Kocuria;s__Kocuria marina',
                     'k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Corynebacteriales;f__Nocardiaceae;g__Nocardia;s__Nocardia veterana',
                     'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rickettsiales;f__;g__Candidatus Arcanobacter;s__Candidatus Arcanobacter lacustris', 
                     'k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Corynebacteriales;f__Nocardiaceae;g__Nocardia;s__Nocardia otitidiscaviarum', 
                     'k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Corynebacteriales;f__Nocardiaceae;g__Nocardia;s__Nocardia jiangxiensis', 
                     'k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Corynebacteriales;f__Nocardiaceae;g__Nocardia;s__Nocardia niigatensis', 
                     'k__Bacteria;p__Verrucomicrobia;c__;o__;f__;g__;s__Verrucomicrobia bacterium SCGC AAA164-L15',
                     'k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__;f__;g__;s__Gammaproteobacteria bacterium MFB021', 
                     'k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Mangrovibacter;s__Mangrovibacter sp. MFB070',
                     'k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Ralstonia;s__Ralstonia sp. PBA', 
                     'k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Shigella;s__Shigella dysenteriae', 
                     'k__Bacteria;p__Proteobacteria;c__Deltaproteobacteria;o__Myxococcales;f__;g__Enhygromyxa;s__Enhygromyxa salina', 
                     'k__Bacteria;p__Bacteroidetes;c__Cytophagia;o__Cytophagales;f__Cyclobacteriaceae;g__Lunatimonas;s__Lunatimonas lonarensis',
                     'k__Bacteria;p__Verrucomicrobia;c__;o__;f__;g__;s__Verrucomicrobia bacterium SCGC AAA164-O14', 
                     'k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Corynebacteriales;f__Nocardiaceae;g__Nocardia;s__Nocardia pneumoniae', 
                     'k__Bacteria;p__Acidobacteria;c__Acidobacteriia;o__Acidobacteriales;f__Acidobacteriaceae;g__Terracidiphilus;s__Terracidiphilus gabretensis',
                     'k__Bacteria;p__Acidobacteria;c__Acidobacteriia;o__Acidobacteriales;f__Acidobacteriaceae;g__Silvibacterium;s__Silvibacterium bohemicum',
                     'k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Lysinibacillus;s__Lysinibacillus xylanilyticus']

def filter_zebra(df, zebra):
    if zebra == 'Ross_WOL':
        zebra_df = pd.read_csv(os.getcwd() + '/zebra/ross_innes_wol_zebra_out.txt', sep ='\t')
    elif zebra == 'BE_WOL':
        zebra_df = pd.read_csv(os.getcwd() + '/zebra/BE_wol_zebra_out.txt', sep ='\t')
    elif zebra == 'EAC_WOL':
        zebra_df = pd.read_csv(os.getcwd() + '/zebra/ICGC_wol_zebra_out.txt', sep = '\t')
    elif zebra == 'EAC_TCGA_WOL':
        zebra_df = pd.read_csv(os.getcwd() + '/zebra/TCGA_EAC_wol_zebra_out.txt', sep = '\t')
    elif zebra == 'norm_GERD_WOL':
        zebra_df = pd.read_csv(os.getcwd() + '/zebra/mixed_Esoph_wol_zebra_out.txt', sep = '\t')
    
    #Subset zebra to at least 1% coverage
    zebra_df = zebra_df[zebra_df['coverage_ratio'] >= 0.01]
    
    #Subset biom with zebra
    df = df[df.index.isin(zebra_df['gotu'])]
    
    return(df)

def species_from_zebra(custom_qza, export_filepath, fn, meta):
    #Create species table from zebra
    
    #Convert from genome to species and output that as well
    taxdf = pd.read_csv(os.getcwd() + '/qiita_downloads/WOL_lineages.txt', sep='\t', index_col=0, header=None)
    taxdf.columns = ['Taxon']
    taxdf.index.name = 'Feature ID'
    taxonomy = Artifact.import_data('FeatureData[Taxonomy]', taxdf)
    species_qza = taxa.methods.collapse(table=custom_qza, taxonomy=taxonomy, level=7).collapsed_table
    
    #Convert to pandas df
    df = species_qza.view(pd.DataFrame)
    df = df.T
    
    #Change fn from genome to species
    fn = fn.replace('genome', 'species')
    
    #Export metadata
    meta_filename = export_filepath + 'metadata/' +'metadata_' + fn + '.tsv'
    meta.to_csv(meta_filename, sep = '\t', index = False)
         
    #Export biom table in the form of pandas df 
    filename = export_filepath + 'pandas_df/' + fn + '.tsv'
    df.to_csv(filename, sep = '\t')
    
    #Export qza
    species_qza.save(export_filepath + 'qza/' + fn + '.qza')
    
    #Save qza as biom table
    Artifact.export_data(species_qza, export_filepath + 'biom/' + fn + '.biom')
    
    return()

def data_split_helper(biom, meta, fn, zebra=False, rare_level=25000, zebra_rare=False, remove_bad_microbes=False):
    'Given metadata + biom, create all datatypes needed'
    
    #Create custom df of biom table (just selected disease type)
    table = load_table(biom).to_dataframe()
    df = table[table.columns.intersection(meta['sample_name'].tolist())]
    
    if remove_bad_microbes == True:
        #Remove 'bad' microbes
        df = df.drop(bad_microbes_list)
    
    print('Total number of samples in', fn, ' = ', df.shape[1])
    
    #Create custom qza table (just selected disease type)
    custom_qza = Artifact.import_data("FeatureTable[Frequency]", df.T)
    
    #Create rarified table off qza table
    rare_qza, = rarefy(table=custom_qza, sampling_depth = rare_level)
    print(rare_qza)
    
    #Convert rare table into pandas df
    rare_df = rare_qza.view(pd.DataFrame)
    
    if zebra != False:
        print('Zebra - non-Rare')
        df = filter_zebra(df, zebra)
        
    if zebra_rare != False:
        print('Zebra - Rare')
        rare_df = filter_zebra(rare_df, zebra_rare)

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
    
    if zebra != False:
        #Convert from genome to species and output that as well
        species_from_zebra(custom_qza, export_filepath, fn, meta)
        
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

def data_split_helper_csv(meta, progression, timepoint, csv, fn , counts=False):
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
    
    
def Sam_BE_data_split(progression, timepoint, csv, fn, counts=False):
    '''To split csv table from BE Samples from Sam using MetaPhlan4 '''
    #Leave counts=False if using relative abudances data, else change to counts=True
    
    #Import metadata from equivalent study on Qiita
    meta = pd.read_csv(os.getcwd() + '/qiita_downloads/qiita14598_BE_Esoph/sample_information_from_prep_14498.tsv', sep = '\t')
    
    #Convert from Sam to our format all files based on imported csv/meta_custom
    data_split_helper_csv(meta, progression, timepoint, csv, fn, counts=counts) 
    
    
def be_data_split(progression, timepoint, biom, fn, zebra=False, exact_pairs= False, trim=False):
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
    
    #Create all files based on biom/meta_custom
    data_split_helper(biom, meta_custom, fn, zebra=zebra)
    
    return()
    

def gerd_normal_data_split(disease_type, biom, fn, zebra=False, trim=False, shotgun=True):
    '''To split biom table from 14458, into a normal and gerd biom table'''
    
    #Import metadata from Qiita
    meta = pd.read_csv(os.getcwd() + '/qiita_downloads/qiita14458_NormalGerdEsoph/sample_information_from_prep_14487.tsv', sep = '\t')
    if shotgun == False:
        meta = pd.read_csv(os.getcwd() + '/qiita_downloads/qiita14325_16S_NormalGerdEsoph/matched_w_shotgun_sample_information_from_prep_12316.tsv', sep = '\t')
    
    #Create custom metadata for disease of intrest (just selected disease type)
    meta_custom =  meta[(meta.diagnosis == disease_type)]
    meta_custom = meta_custom.reset_index(drop = True)
    
    if trim!= False:
        meta_custom = meta_custom.sample(trim)
    display(meta_custom[:3])
    
    #Create all files based on biom/meta_custom
    data_split_helper(biom, meta_custom, fn, zebra=zebra)
    
    return()

def eac_14857_data_split(biom, fn, trim=False, remove_bad_microbes=False, Ginny_list=False): 
    
    #Import metadata from Qiita
    meta = pd.read_csv(os.getcwd() + '/qiita_downloads/qiita14857_EAC_ICGC/sample_information_from_ESAD_ICGC_14857_prep_13977.txt', sep = '\t')
    
    if Ginny_list == True:
        meta = meta[~meta.isin(Ginny_donor).any(axis=1)]
    
    if trim!= False:
        meta = meta.sample(trim)
    
    display(meta[:3])
    
    #Create all files based on biom/meta_custom
    data_split_helper(biom, meta, fn, zebra=False, remove_bad_microbes=remove_bad_microbes)
   
    return()

def eac_data_split(biom, fn, zebra=False, trim=False, stage=False): 
    
    #Import metadata from Qiita
    meta = pd.read_csv(os.getcwd() + '/qiita_downloads/qiita14521_EAC_ICGC/sample_info_14521_20220509-203030.txt', sep = '\t')
    
    if stage!= False:
        meta = pd.read_csv(os.getcwd() + '/qiita_downloads/qiita14521_EAC_ICGC/sample_info_14521_20220509-203030_extended.txt',
                           sep = '\t')
        if stage == 'low': #1 or 2
            meta = meta[(meta.tumour_grade.isin([1,2]))]
        elif stage == 'high': #3
            meta = meta[(meta.tumour_grade.isin([3]))]  
        meta = meta.reset_index(drop = True)
    
    if trim!= False:
        meta = meta.sample(trim)
    
    display(meta[:3])
    
    #Create all files based on biom/meta_custom
    data_split_helper(biom, meta, fn, zebra=zebra)
   
    return()

def eac_tcga_data_split(biom, fn, zebra=False, trim=False): 
    
    #Import metadata from Qiita
    meta = pd.read_csv(os.getcwd() + '/qiita_downloads/qiita14990_EAC_TCGA/sample_information_from_prep_14621.tsv', sep = '\t')
    
    if trim!= False:
        meta = meta.sample(trim)
    
    display(meta[:3])
    
    #Create all files based on biom/meta_custom
    data_split_helper(biom, meta, fn, zebra=zebra)
   
    return()

def normal_be_eac_data_split(disease_type, biom, fn, exact_pairs=False, zebra=False, remove_bad_microbes=False, Ginny_list=False):
    '''To split biom table from TBD, into a normal, BE, EAC biom table'''
    #exact pairs just takes a single BE and EAC sample from each patient
    
    #Import metadata from Qiita
    #meta = pd.read_csv(os.getcwd() + '/qiita_downloads/qiitaTBD_Norm_BE_EAC/qiita_sample_info.txt', sep = '\t')
    meta = pd.read_csv(os.getcwd() + '/compare_rossICGC/corrected_meta.tsv', sep = '\t')
    
    #Remove all the samples from the 'AHM1051'
    meta = meta[(meta.Patient != 'AHM1051')]
    
    #Create custom metadata for disease of intrest (just selected disease type)
    meta_custom = meta[(meta.Sample_type.isin(disease_type))]
    
    if exact_pairs == True:
        new_meta = pd.DataFrame(columns = meta_custom.columns)
        for index, row in meta_custom.iterrows():
            if any((new_meta['Patient'] == row['Patient']) & (new_meta['Sample_type'] == row['Sample_type'])) == False:
                new_meta.loc[len(new_meta)] = row
        meta_custom = new_meta
        
    if Ginny_list == True:
        meta_custom = meta_custom[meta_custom.isin(Ginny_oac_bo).any(axis=1)]
    
    meta_custom = meta_custom.reset_index(drop = True)
    display(meta_custom[:3])
    
    #Create all files based on biom/meta_custom
    data_split_helper(biom, meta_custom, fn, zebra=zebra, remove_bad_microbes=remove_bad_microbes)
    
    return()

############################ AIM 3 ########################

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