"""
util.py
=======
Author - Daniel Monyak
6-14-24
=======


Provides
    Helper functions needed by Select_fCpGs notebook

"""

import pandas as pd
import numpy as np
import os
import sys
from sklearn.cluster import KMeans
import EpiClockNBL.src.util as pc_util
pc_consts = pc_util.consts

## Defines acceptable ranges (inclusive) for each value
CLOCK_CRITERIA = {
    'normal':{
        'beta_nan_fracs':(0, 0.05),
        'beta_means':(0.4, 0.6)
    },
    'tumor':{
        'beta_nan_fracs':(0, 0.05),
        'beta_means':(0.4, 0.6)
    }
}

def getNeutralDNACpGs(chip_27k=False):
    """
    Returns
    -------
    List of "neutral" CpGs
        i.e. not associated with any gene or regulatory feature
    
    Notes
    -----
    Checks the manifests of both the 450K and 850K array
    """
    
    chip_annot_dir = os.path.join(pc_consts['BRCA_official_indir'], 'TCGA', 'chip_annots')
    manifest_450K = pd.read_table(os.path.join(chip_annot_dir, 'chip450_annot_cleaned.txt'), index_col=0, low_memory=False)
    manifest_850K = pd.read_table(os.path.join(chip_annot_dir, 'chip850_annot_cleaned.txt'), index_col=0, low_memory=False)
    
    # Requirements for neutral CpGs
    neutral_filters = [
        manifest_450K['Regulatory_Feature_Group'].isna(),
        manifest_450K['UCSC_RefGene_Name'].isna(),
        manifest_850K['Regulatory_Feature_Group'].isna(),
        manifest_850K['UCSC_RefGene_Name'].isna(),
    ]
    
    if chip_27k:
        neutral_filters.append(manifest_450K['Methyl27_Loci'])
    
    # Combine the filters to select neutral CpGs
    neutral_DNA_CpG_list = manifest_450K.index[pc_util.combineFilters(neutral_filters)].values
    
    # Sanity checks
    neutral_manifest_450K = manifest_450K.loc[neutral_DNA_CpG_list]
    neutral_manifest_850K = manifest_850K.loc[neutral_DNA_CpG_list]
    assert neutral_manifest_450K['Regulatory_Feature_Group'].isna().all()
    assert neutral_manifest_450K['UCSC_RefGene_Name'].isna().all()
    assert neutral_manifest_850K['Regulatory_Feature_Group'].isna().all()
    assert neutral_manifest_850K['UCSC_RefGene_Name'].isna().all()

    return neutral_DNA_CpG_list

def getDataDict(data_paths, filter_tum_samps=False, assert_same_sites=True):
    """
    Return a dictionary that holds data for the the TCGA and normal cohorts
    
    For each cohort, dictionary will hold
        beta_values - # samples x # CpGs
        pureSamples - list of samples after filtering for purity
            this step is already done in TCGA
    
    Parameters
    ----------
    data_paths : dictionary with cohort names ("tumor" and "normal") as keys and the path to the beta value table file as the values
    
    Returns
    -------
    data : dictionary
        Dictionary with 2 levels
        First level is between 'tumor' and 'normal'
            
    Notes
    -----
    Importing data can take up to 10 minutes once file is downloaded (if on an external file source)
        But should be shorter than this
    """
    
    # Data dictionary
    data = {'cohorts':[]}
    
    for cohort in data_paths:
        data[cohort] = {'beta_values':pd.read_table(data_paths[cohort], index_col=0)}
        data['cohorts'].append(cohort)
    
    idx_tumor = np.sort(data['tumor']['beta_values'].index.values)
    data['allCpGs'] = idx_tumor
    data['tumor']['beta_values'] = data['tumor']['beta_values'].loc[data['allCpGs']]
    
    # Clip unnecessary accessions suffixes from sample IDs
    data['tumor']['beta_values'] = data['tumor']['beta_values'].rename(lambda x:'-'.join(x.split('-')[:4]), axis='columns')
    
    # All TCGA samples have already been filtered for purity
    if filter_tum_samps:
        data['tumor']['purity'] = pc_util.getLUMP_values(data['tumor']['beta_values'])
        data['tumor']['pureSamples'] = data['tumor']['purity'].index[
            data['tumor']['purity'] >= pc_consts['LUMP_threshold']
        ].values
    else:
        data['tumor']['pureSamples'] = data['tumor']['beta_values'].columns.values
    
    if 'normal' in data:
        
        # Get pure samples for normals
        # Filter by LUMP value
        data['normal']['purity'] = pc_util.getLUMP_values(data['normal']['beta_values'])
        data['normal']['pureSamples'] = data['normal']['purity'].index[
            data['normal']['purity'] >= pc_consts['LUMP_threshold']
        ].values
        
        
        idx_normal = np.sort(data['normal']['beta_values'].index.values)
        
        if assert_same_sites:
            # Check that CpG Sites are the same in each DF
            assert np.all(idx_tumor == idx_normal)
        else:
            data['allCpGs'] = np.sort(np.intersect1d(idx_tumor, idx_normal)) # shared sites
            data['tumor']['beta_values'] = data['tumor']['beta_values'].loc[data['allCpGs']]
            

        # Use the same order of CpGs in the index
        data['normal']['beta_values'] = data['normal']['beta_values'].loc[data['allCpGs']]

    return data

def addMeanStdsNans(data):
    """
    Adds data objects to input dictionary
    
    For the following statistics for each CpG across all samples in the data
        mean
        standard deviation
        # NaN values
    Adds beta_values_SELECTION to the dictionary only if it's not there already
    
    Parameters
    ----------
    data : dictionary returned by getDataDict
    """
    
    for cohort in data['cohorts']:
        if 'beta_values_SELECTION' not in data[cohort]:
            data[cohort]['beta_values_SELECTION'] = data[cohort]['beta_values'][data[cohort]['pureSamples']]

        data[cohort]['beta_means'] = data[cohort]['beta_values_SELECTION'].mean(axis=1)
        data[cohort]['beta_stds'] = data[cohort]['beta_values_SELECTION'].std(axis=1)
        data[cohort]['beta_nan_fracs'] = data[cohort]['beta_values_SELECTION'].isna().mean(axis=1)


def getCpG_list(data, criteria, starting_CpG_list=None, n_select=None, sample_rank='tumor', stat_rank='beta_stds', good_end='higher', n_select_frac=None):
    """
    Return a set of CpGs that fit certain criteria
    
    Parameters
    ----------
    data : dictionary returned by getDataDict
    criteria : dictionary that defines acceptable ranges (inclusive) for each data type in each cohort
    starting_CpG_list : CpGs to consider initially
    n_select : number of CpGs to select; if None, return all CpGs that fit criteria
    sample_rank : If n_select is not None, use stat_rank from this cohort to rank the tumors
    stat_rank : If n_select is not None, use this stat from the sample_rank cohort to rank the tumors
    good_end : 'higher' or 'lower' - which end of stat_rank should be selected if n_select is not None
    
    Returns
    -------
    CpG_list : ndarray of names of CpGs that fit criteria
               only include those in starting_CpG_list
               limiting to n_select tumors, ranking using stat_rank in the sample cohort in sample_rank
    """
    
    if 'allCpGs' not in data:
        data['allCpGs'] = data['tumor']['beta_values'].index.values
    if starting_CpG_list is None:
        starting_CpG_list = data['allCpGs']
    
    siteFilters = []
    for sample in criteria.keys():
        samp_criteria = criteria[sample]
        for stat in samp_criteria.keys():
            siteFilters.append(data[sample][stat] >= samp_criteria[stat][0])
            siteFilters.append(data[sample][stat] <= samp_criteria[stat][1])
            
    combined_filter = pc_util.combineFilters(siteFilters)
    criteria_CpGs = np.intersect1d(data['allCpGs'][combined_filter], starting_CpG_list)
    if n_select is None:
        if n_select_frac is None:
            return criteria_CpGs
        
        n_select = int(n_select_frac * len(criteria_CpGs))
    
    print(f'Selecting from {len(criteria_CpGs)} balanced sites')
    
    return data[sample_rank][stat_rank].loc[criteria_CpGs].sort_values(ascending=(good_end == 'lower')).index[:n_select].values


def gen_CpG_set(data, neutral_DNA_CpG_list, n_select=None, criteria=None, n_select_frac=None):
    """
    Return the Clock set of CpGs
    
    Parameters
    ----------
    data : dictionary returned by getDataDict
    neutral_DNA_CpG_list : list of neutral CpGs returned by getNeutralDNACpGs
    n_select : how many CpGs to return
    
    Returns
    -------
    Clock_CpGs : ndarray of names of CpGs that fit CLOCK_CRITERIA
    """
    
    for cohort in data['cohorts']:
        if cohort in data:
            data[cohort]['beta_values_SELECTION'] = data[cohort]['beta_values'][data[cohort]['pureSamples']]
    
    print('Selecting CpGs with:')
    
    n_tumors = data['tumor']['beta_values_SELECTION'].shape[1]
    print(f'\t{n_tumors} TCGA samples')
    if 'normal' in data:
        n_normals = data['normal']['beta_values_SELECTION'].shape[1]
        print(f'\t{n_normals} normal samples')
    
    
    # Add beta_means, beta_stds, beta_nan_fracs to data dict
    addMeanStdsNans(data)
    
    if 'allCpGs' not in data:
        data['allCpGs'] = data['tumor']['beta_values'].index.values
    
    if criteria is None:
        criteria = CLOCK_CRITERIA
    
    Clock_CpGs = getCpG_list(data, criteria, starting_CpG_list=neutral_DNA_CpG_list, n_select=n_select, n_select_frac=n_select_frac)
    
    return Clock_CpGs


def clusteringWeights(km_beta_values, random_state=None):
    """
    Return the clustering weights w_s of each CpG site s
    See "Calculating clustering weight" notebook in "Select_fCpGs" directory
    
    Parameters
    ----------
    km_beta_values : DataFrame of beta values (# CpGs x # tumors)
    random_state : random_state variable to pass to KMeans for deterministic output
    
    Returns
    -------
    w : Series with clustering weight of each site
    """
    km = KMeans(n_clusters=4, random_state=random_state).fit(km_beta_values.T)
    
    n = pd.Series(km.labels_).value_counts().sort_index()
    k_dict = {}
    for i in range(n.shape[0]-1):
        for j in range(i+1, n.shape[0]):
            k_dict[(i, j)] = n.loc[i] * n.loc[j]
    k = pd.Series(data=k_dict)
    k_frac = k / k.sum()
    
    w = pd.Series(index=km_beta_values.index, data=0)
    n_clusters = km.cluster_centers_.shape[0]
    for i in range(n_clusters-1):
        for j in range(i+1, n_clusters):
            centroid_diff = pd.Series(
                data=np.abs(km.cluster_centers_[i] - km.cluster_centers_[j]),
                index=km_beta_values.index)
            w += centroid_diff * k_frac.loc[(i, j)]
    
    return w

def getBinEdges(start, stop, binwidth, hardStop=True):
    """
    Generate a list of bin edges given a start, stop, and bin width.

    Parameters:
    -----------
    start : float
        The starting edge of the binning range.
    stop : float
        The stopping edge of the binning range.
    binwidth : float
        The width of each bin.
    hardStop : bool, default=True
        If True, ensures that the last bin edge is exactly `stop` and works backwards to `start`.
        If False, starts from `start` and increments bin edges by `binwidth` until it passes `stop`.

    Returns:
    --------
    list of float
        A list of bin edges starting from `start` (or ending at `stop`, depending on `hardStop`).
    """
    if hardStop:
        bin_edges = [stop]
        # Work backwards from `stop` to `start`, subtracting `binwidth` each time
        while bin_edges[-1] > start:
            bin_edges.append(bin_edges[-1] - binwidth)
        return bin_edges[::-1]  # Return in ascending order
    else:
        bin_edges = [start]
        # Work forwards from `start` to `stop`, adding `binwidth` each time
        while bin_edges[-1] < stop:
            bin_edges.append(bin_edges[-1] + binwidth)
        return bin_edges
