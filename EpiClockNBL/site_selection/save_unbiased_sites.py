import os
import copy
import time
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import EpiClockNBL.util as nbl_util
nbl_consts = nbl_util.consts
from .util import CLOCK_CRITERIA, getDataDict, gen_CpG_set, getBinEdges

# IO directories

proj_dir = os.path.join(nbl_consts['official_indir'], 'TARGET')

figure_outdir = 'figures'
output_dir = 'outputs'
outfile_dir = os.path.join(output_dir, 'outfiles')
outdir = proj_dir
outfile_path = os.path.join(proj_dir, 'beta_values_unbiased_sites.txt')

def pipeline(verbose=True, make_figures=False):
    
    if verbose:
        print(f'Running save-unbiased-sites pipeline with verbose={verbose}, make_figures={make_figures}.')
        time.sleep(1)
        if make_figures:
            print(f'Saving figures to {os.path.abspath(figure_outdir)}')
        else:
            print('Set make_figures=True to generate the accompanying figures.')
        print('#'*100)
        print('#'*100)
        time.sleep(1)

    # Directory setup
    if verbose:
        print('\nDirectory setup...', end=' ')
        time.sleep(1)
        
    os.makedirs(figure_outdir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(outfile_dir, exist_ok=True)
    os.makedirs(outdir, exist_ok=True)
    
    if verbose:
        print('DONE')

    ####################################################################################
    ####################################################################################
    
    # Prepare data
    if verbose:
        print('\nImporting and preparing data, will take 2-30 min...', end=' ')
        time.sleep(1)

    path_dict = {}
    path_dict['tumor'] = os.path.join(proj_dir, 'cohort1.methyl.tsv')

    data = getDataDict(path_dict, filter_tum_samps=False, assert_same_sites=False)
    data['tumor']['beta_values'] = data['tumor']['beta_values'].dropna(axis=1, how='all')
    data['tumor']['pureSamples'] = data['tumor']['pureSamples'][np.isin(data['tumor']['pureSamples'], data['tumor']['beta_values'].columns)]
    
    if verbose:
        print('DONE')
    
    ####################################################################################
    ####################################################################################
    
    # Get unbiased sites
    if verbose:
        print('\nSelect unbiased sites...', end='\n\n')
        time.sleep(1)

    criteria = copy.deepcopy(CLOCK_CRITERIA)
    
    # NO NORMAL SAMPLE
    if 'normal' not in data['cohorts']:
        del criteria['normal']

    unbiased_sites = gen_CpG_set(data, criteria=criteria, neutral_DNA_CpG_list=None)
    data['tumor']['beta_values_SELECTION'].loc[unbiased_sites].to_csv(outfile_path, sep='\t')

    if verbose:
        print('\nDONE')
    
    ####################################################################################
    ####################################################################################
    
    if not make_figures:
        return
    #
    #
    #
    # Figure
    #
    #
    #
    
    if verbose:
        print('\nGenerating figure...', end='\n\n')
        time.sleep(1)

    ## Configure graph
    sf = nbl_consts['sf']
    figsize = np.array([9, 7])
    sns.set(rc={"savefig.bbox":'tight', 'axes.linewidth':sf}, font_scale=1, style='ticks')

    mean_only_criteria = copy.deepcopy(CLOCK_CRITERIA)
    del mean_only_criteria['normal']
    del mean_only_criteria['tumor']['beta_nan_fracs']

    mean_only_sites = gen_CpG_set(data, criteria=mean_only_criteria, neutral_DNA_CpG_list=None)

    ## Inter-tumor mean
    # Plot inter-tumor mean beta for each unbiased site

    # Create plot
    fig, ax = plt.subplots(1, 1, figsize=figsize * sf)

    binwidth=0.01

    bins = np.concatenate([
        getBinEdges(data['tumor']['beta_means'].min(), data['tumor']['beta_means'].loc[mean_only_sites].min(), binwidth, hardStop=True),
        getBinEdges(data['tumor']['beta_means'].loc[mean_only_sites].max(), data['tumor']['beta_means'].max(), binwidth, hardStop=False)
    ])
    sns.histplot(ax=ax,      # Biased sites
                x=data['tumor']['beta_means'].loc[list(set(data['allCpGs']).difference(set(mean_only_sites)))],
                color=nbl_consts['palette_jco'][2], alpha=nbl_consts['opacity'],
                bins = bins
                )
    sns.histplot(ax=ax,      # Unbiased sites
                x=data['tumor']['beta_means'].loc[mean_only_sites],
                color=nbl_consts['palette_jco'][1], alpha=nbl_consts['opacity'], binwidth=binwidth
                )

    # Customize figure
    ax.set_xlabel(r'Average $\beta$ in TARGET', fontsize=nbl_consts['labelfontsize'] * sf)
    ax.set_ylabel(ax.get_ylabel(), fontsize=nbl_consts['labelfontsize'] * sf)
    ax.set_title(f'All CpGs (n = {data["tumor"]["beta_means"].shape[0]:,})', fontsize=nbl_consts['labelfontsize'] * sf)
    ax.tick_params(axis='both', labelsize=nbl_consts['ticksfontsize'] * sf, width=sf, length=8 * sf)

    # Save figure
    fig.savefig(os.path.join(figure_outdir, 'pick_unbiased_sites.svg'), format='svg', pad_inches=0.1)
    fig.show()

    if verbose:
        print('\nDONE')