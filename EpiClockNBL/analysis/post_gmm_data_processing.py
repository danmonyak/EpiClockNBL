"""
Source code for processing the TARGET and Henrich cohort data after GMM fitting

"""

import pandas as pd
import numpy as np
import os
import time
import EpiClockNBL.util as nbl_util
nbl_consts = nbl_util.consts

# IO directories
proj_dir_TARGET = os.path.join(nbl_consts['official_indir'], 'TARGET')
proj_dir_Henrich = os.path.join(nbl_consts['official_indir'], 'Henrich')

def pipeline(verbose=True):
    
    if verbose:
        print(f'Running post_gmm_data_processing pipeline with verbose={verbose}.')
        time.sleep(1)
        print('#'*100)
        print('#'*100)
        time.sleep(1)

    ####################################################################################
    ####################################################################################
    
    # Load clinical tables
    if verbose:
        print('Loading clinical tables...', end=' ')
        time.sleep(1)

    clinical = {}
    clinical['TARGET'] = pd.read_table(os.path.join(proj_dir_TARGET, 'clinical.processed.tsv')).set_index('submitter_id')
    clinical['Henrich'] = pd.read_table(os.path.join(proj_dir_Henrich, 'clinical.processed.tsv')).set_index('Sample_geo_accession')

    if verbose:
        print('DONE')

    ####################################################################################
    ####################################################################################
    
    # Load GMM results
    if verbose:
        print('Loading GMM results...', end=' ')
        time.sleep(1)

    gmm_results_target = pd.read_csv(os.path.join(nbl_consts['repo_dir'], 'Gaussian Mixture Model', 'TARGET.GMM_results.csv'), index_col=0)
    gmm_results_target.index = gmm_results_target.index.map(lambda x:x.replace('.', '-'))

    gmm_results_henrich = pd.read_csv(os.path.join(nbl_consts['repo_dir'], 'Gaussian Mixture Model', 'Henrich.GMM_results.csv'), index_col=0)

    # Sensitivity analysis
    gmm_results_target_sensitivity_split1 = pd.read_csv(os.path.join(nbl_consts['repo_dir'], 'Gaussian Mixture Model', 'SENSITIVITY_SPLIT1.TARGET.GMM_results.csv'), index_col=0)
    gmm_results_target_sensitivity_split1.index = gmm_results_target_sensitivity_split1.index.map(lambda x:'-'.join(x.replace('.', '-').split('-')[:-1]))

    gmm_results_target_sensitivity_split2 = pd.read_csv(os.path.join(nbl_consts['repo_dir'], 'Gaussian Mixture Model', 'SENSITIVITY_SPLIT2.TARGET.GMM_results.csv'), index_col=0)
    gmm_results_target_sensitivity_split2.index = gmm_results_target_sensitivity_split2.index.map(lambda x:'-'.join(x.replace('.', '-').split('-')[:-1]))

    if verbose:
        print('DONE')

    ####################################################################################
    ####################################################################################

    # Add mitotic age estimate to clinical tables
    if verbose:
        print('Adding mitotic age estimate to clinical tables...', end=' ')
        time.sleep(1)

    # Add GMM phi
    clinical['TARGET'] = clinical['TARGET'].merge(gmm_results_target['phi.mean'].rename('phi'), left_on='sampleID', right_index=True, validate="one_to_one")
    clinical['Henrich'] = clinical['Henrich'].merge(gmm_results_henrich['phi.mean'].rename('phi'), left_index=True, right_index=True, validate="one_to_one")

    # Sensitivity analysis
    clinical['TARGET'] = clinical['TARGET'].merge(gmm_results_target_sensitivity_split1['phi.mean'].rename('phi.split1'), left_on='sampleID', right_index=True, validate="one_to_one")
    clinical['TARGET'] = clinical['TARGET'].merge(gmm_results_target_sensitivity_split2['phi.mean'].rename('phi.split2'), left_on='sampleID', right_index=True, validate="one_to_one")

    if verbose:
        print('DONE')
    
    ####################################################################################
    ####################################################################################

    # Save clinical tables
    if verbose:
        print('Saving clinical tables...', end='\n\n')
        time.sleep(1)

    nbl_util.writeWithoutOverwrite(
        os.path.join(proj_dir_TARGET, 'clinical.processed.withAges.tsv'),
        clinical['TARGET'],
        lambda filepath, data: data.reset_index(names='submitter_id').to_csv(filepath, sep='\t', index=False),
        lambda filepath: pd.read_table(filepath).set_index('submitter_id'),
        compareFunc=lambda old, new: (old.round(5).fillna(0) == new.round(5).fillna(0)).all(axis=None, skipna=True)
    )
    nbl_util.writeWithoutOverwrite(
        os.path.join(proj_dir_Henrich, 'clinical.processed.withAges.tsv'),
        clinical['Henrich'],
        lambda filepath, data: data.reset_index(names='Sample_geo_accession').to_csv(filepath, sep='\t', index=False),
        lambda filepath: pd.read_table(filepath).set_index('Sample_geo_accession'),
        compareFunc=lambda old, new: (old.round(5).fillna(0) == new.round(5).fillna(0)).all(axis=None, skipna=True)
    )

    if verbose:
        print('\nDONE')
