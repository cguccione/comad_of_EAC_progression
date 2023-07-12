'''
Helper functions used to run machine learning models
'''

#Basic imports
import numpy as np
import pandas as pd
import os

#ML imports
import matplotlib.pyplot as plt 
import seaborn as sns
from scipy.stats import mannwhitneyu, normaltest

from sklearn.feature_selection import chi2, VarianceThreshold, SelectKBest, f_classif, mutual_info_classif, RFE, SelectFromModel
from sklearn.model_selection import train_test_split, StratifiedKFold, LeaveOneOut, ShuffleSplit
from sklearn.linear_model import RidgeCV, LassoCV, Ridge, Lasso, LogisticRegression
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.metrics import RocCurveDisplay, auc
from sklearn import metrics

from sklearn.model_selection import ShuffleSplit
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import RocCurveDisplay, auc
from sklearn import metrics
import matplotlib.pyplot as plt

from statistics import mean
import warnings

def gradientBoosting_ML_cv(table, meta, metric, gb_max_depth=2, gb_n_estimators=100, cross_validation_splits=5):
    
    #Create list with metric to perfrom
    #ML on and convert list into True/False
    #instead of tumor/Background, 
    #and create list of all microbes
    
    X = table.T
    X = X.sort_index(ascending=True)
    X = X.to_numpy()
    meta = meta.sort_values('sample_name', ascending=True)
    y = meta[metric].tolist()
    y = [0 if i == y[0] else 1 for i in y]
    y = np.array(y)
    
    print('Cross validation testing size =', len(y)/cross_validation_splits)
    
    #Resources used to compile code below: 
    #https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc_crossval.html
    
    #Provides train/test indices to split data in train/test sets
    cross_val = StratifiedKFold(n_splits=cross_validation_splits)
    
    #Create gradient boosting classifier model
    gbc = GradientBoostingClassifier(max_depth = gb_max_depth, n_estimators = gb_n_estimators)
    
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)
    fig, ax = plt.subplots()
    
    #Loop through each cross validation scenario
    for i, (train, test) in enumerate(cross_val.split(X, y)):
        gbc.fit(X[train], y[train])
        viz = RocCurveDisplay.from_estimator(
            gbc,
            X[test],
            y[test],
            name="ROC fold {}".format(i),
            alpha=0.3,
            lw=1,
            ax=ax,
        )
        interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
        interp_tpr[0] = 0.0
        tprs.append(interp_tpr)
        aucs.append(viz.roc_auc)
        
    ax.plot([0, 1], [0, 1], linestyle="--", lw=2, color="r", label="Chance", alpha=0.8)

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    ax.plot(
        mean_fpr,
        mean_tpr,
        color="b",
        label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
        lw=2,
        alpha=0.8,
    )
    
    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    ax.fill_between(
        mean_fpr,
        tprs_lower,
        tprs_upper,
        color="grey",
        alpha=0.2,
        label=r"$\pm$ 1 std. dev.",
    )
    
    ax.set(
        xlim=[-0.05, 1.05],
        ylim=[-0.05, 1.05],
        title="Receiver operating characteristic example",
    )
    
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    plt.show()

def gradientBoosting_ML_ss(table, meta, metric, gb_max_depth=2, gb_n_estimators=100, cross_validation_splits=5, test_size = 0.3):
    
    #Create list with metric to perfrom
    #ML on and convert list into True/False
    #instead of tumor/Background, 
    #and create list of all microbes
    
    X = table.T
    X = X.sort_index(ascending=True)
    X = X.to_numpy()
    meta = meta.sort_values('sample_name', ascending=True)
    y = meta[metric].tolist()
    y = [0 if i == y[0] else 1 for i in y]
    y = np.array(y)
    
    print('Cross validation testing size =', (len(y) * (1- test_size))/cross_validation_splits)
    
    #Resources used to compile code below: 
    #https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc_crossval.html
    
    #Provides train/test indices to split data in train/test sets
    shuffle = ShuffleSplit(n_splits=cross_validation_splits, test_size=test_size)
    
    #Create gradient boosting classifier model
    gbc = GradientBoostingClassifier(max_depth = gb_max_depth, n_estimators = gb_n_estimators)
    
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)
    fig, ax = plt.subplots()
    
    #Loop through each cross validation scenario
    for i, (train, test) in enumerate(shuffle.split(X, y)):
        gbc.fit(X[train], y[train])
        viz = RocCurveDisplay.from_estimator(
            gbc,
            X[test],
            y[test],
            name="ROC fold {}".format(i),
            alpha=0.3,
            lw=1,
            ax=ax,
        )
        interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
        interp_tpr[0] = 0.0
        tprs.append(interp_tpr)
        aucs.append(viz.roc_auc)
        
    ax.plot([0, 1], [0, 1], linestyle="--", lw=2, color="r", label="Chance", alpha=0.8)

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    ax.plot(
        mean_fpr,
        mean_tpr,
        color="b",
        label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
        lw=2,
        alpha=0.8,
    )
    
    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    ax.fill_between(
        mean_fpr,
        tprs_lower,
        tprs_upper,
        color="grey",
        alpha=0.2,
        label=r"$\pm$ 1 std. dev.",
    )
    
    ax.set(
        xlim=[-0.05, 1.05],
        ylim=[-0.05, 1.05],
        title="Receiver operating characteristic example",
    )
    
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    plt.show()
    
def ML(fn, metric, cross_validation_splits=10, use_non_neutral=False):
    
    #Import biom table in pandas df/tsv form for current dataset
    biom_df = pd.read_csv('processed_data/pandas_df/' + fn + '.tsv', sep = '\t', index_col=0)
    biom_df = biom_df.astype(int)
    
    if use_non_neutral != False:
        #This uses only the non-neutral microbes for the provided input
        nn = pd.read_csv('outputs/non_neutral/' + use_non_neutral + '.tsv', sep = '\t')
        biom_df = biom_df.loc[nn['Taxonomy']]
        print('Using only top', len(biom_df), 'non-neutral microbes from', use_non_neutral)
    
    #Import metadata for current dataset
    meta = pd.read_csv('processed_data/metadata/metadata_' + fn + '.tsv', sep = '\t')
    #Subset meta to only include those in the df below
    meta = meta[meta['sample_name'].isin(list(biom_df.columns))]
    
    meta_temp = meta.sort_values('sample_name', ascending=True)
    all_cats = meta_temp[metric].tolist()
    print(set(all_cats))
    
    #Run machine learning stratified-split
    print('Stratfied-split')
    gradientBoosting_ML_ss(biom_df, meta, metric, cross_validation_splits=cross_validation_splits)
    
    print('Cross Validation')
    gradientBoosting_ML_cv(biom_df, meta, metric, cross_validation_splits=cross_validation_splits)

def ML_combine(fn, metric, table_list, cross_validation_splits=10):
    
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
        
    #Convert table to all ints
    table = table.fillna(0)
    biom_df = table.astype(int)
    
    #Subset meta to only include those in the df below
    meta = meta[meta['sample_name'].isin(list(biom_df.columns))]
    
    #Run machine learning stratified-split
    print('Stratfied-split')
    gradientBoosting_ML_ss(biom_df, meta, metric, cross_validation_splits=cross_validation_splits)
    
    print('Cross Validation')
    gradientBoosting_ML_cv(biom_df, meta, metric, cross_validation_splits=cross_validation_splits)
    
def ML_voom(fn, metric, cross_validation_splits=10):
    
    #Import biom table in pandas df/tsv form for current dataset
    biom_df = pd.read_csv('voom_output/' + fn + '.tsv', sep = '\t', index_col=0)
    biom_df = biom_df.astype(int)
    biom_df = biom_df.T
    
    #Import metadata for current dataset
    meta = pd.read_csv('voom_input/meta_' + fn + '.tsv', sep = '\t')
    #Add 'X' that matches voom
    meta['sample_name'] = 'X' + meta['sample_name'].astype(str)
    #Subset meta to only include those in the df below
    meta = meta[meta['sample_name'].isin(list(biom_df.columns))]
    
    #Run machine learning stratified-split
    print('Stratfied-split')
    gradientBoosting_ML_ss(biom_df, meta, metric, cross_validation_splits=cross_validation_splits)
    
    print('Cross Validation')
    gradientBoosting_ML_cv(biom_df, meta, metric, cross_validation_splits=cross_validation_splits)

def ML_voom_Sunday(fn, primary_cancer_type, host_body_site, sample_type, metric, cross_validation_splits=10, special_Case=False):
    
    #Import biom table in pandas df/tsv form for current dataset
    biom_df = pd.read_csv('voom_output/TCGAprimay_AmirTumor_genome_RS210-working.tsv', sep = '\t', index_col=0)
    biom_df = biom_df.astype(int)
    
    #Import metadata for current dataset
    meta = pd.read_csv('voom_input/meta_TCGAprimay_AmirTumor_genome_RS210.csv', sep = ',')
    #Add 'X' that matches voom
    meta['sample_name'] = 'X' + meta['sample_name'].astype(str)
    
    #Subset meta to only include what you want
    meta = meta[(meta.primary_cancer_type.isin(primary_cancer_type)) & 
                        (meta.host_body_site.isin(host_body_site)) &
                (meta.sample_type.isin(sample_type))]
    
    if special_Case == True:
        #SPECIAL CASE
        meta = meta[meta['sample_name'].str.contains('HCC')==False]
        
    
    #Subset biom_df based on meta
    biom_df = biom_df[biom_df.index.isin(meta['sample_name'])]
    biom_df = biom_df.T
    
    #Run machine learning stratified-split
    print('Stratfied-split')
    gradientBoosting_ML_ss(biom_df, meta, metric, cross_validation_splits=cross_validation_splits)
    
    print('Cross Validation')
    gradientBoosting_ML_cv(biom_df, meta, metric, cross_validation_splits=cross_validation_splits)
    
    
def voom_data_split(primary_cancer_type, host_body_site, sample_type, special_Case=False):
    
    #Import biom table in pandas df/tsv form for current dataset
    #biom_df = pd.read_csv('voom_output/TCGAprimay_AmirTumor_genome_RS210-working.tsv', sep = '\t', index_col=0)
    biom_df = pd.read_csv('voom_output/TCGAprimay_AmirTumor_genome_RS210.tsv', sep = '\t', index_col=0)
    biom_df = biom_df.astype(float)
    
    #Import metadata for current dataset
    meta = pd.read_csv('voom_input/meta_TCGAprimay_AmirTumor_genome_RS210.csv', sep = ',')
    #Add 'X' that matches voom
    #meta['sample_name'] = 'X' + meta['sample_name'].astype(str)
    
    #Subset meta to only include what you want
    meta = meta[(meta.primary_cancer_type.isin(primary_cancer_type)) & 
                        (meta.host_body_site.isin(host_body_site)) &
                (meta.sample_type.isin(sample_type))]
    
    if special_Case == True:
        #SPECIAL CASE
        meta = meta[meta['sample_name'].str.contains('HCC')==False]
        
    #Subset biom_df based on meta
    biom_df = biom_df[biom_df.index.isin(meta['sample_name'])]
    biom_df = biom_df.T
        
    return(biom_df, meta.reset_index(drop=True))
    