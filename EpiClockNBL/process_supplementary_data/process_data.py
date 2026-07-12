"""
Source code for processing the Henrich cohort data

1. Process clinical file
2. Import beta values
3. Filter samples

"""

import pandas as pd
import numpy as np
import os
import time
import EpiClockNBL.util as nbl_util
nbl_consts = nbl_util.consts

# IO directories
proj_dir = os.path.join(nbl_consts['official_indir'], 'Henrich')

def pipeline(verbose=True):
    
    if verbose:
        print(f'Running process_data pipeline with verbose={verbose}.')
        time.sleep(1)
        print('#'*100)
        print('#'*100)
        time.sleep(1)

    ####################################################################################
    ####################################################################################
    
    # Load beta values
    if verbose:
        print('Loading beta values...', end=' ')
        time.sleep(1)


    ## Beta values
    beta_values = pd.read_table(os.path.join(proj_dir, 'beta_values.txt'), index_col=0)

    if verbose:
        print('DONE')

    
    ####################################################################################
    ####################################################################################
    
    # Select subset with Clock CpGs
    if verbose:
        print('Select intersecting Clock CpGs of beta values...', end=' ')
        time.sleep(1)


    Clock_CpGs = np.loadtxt(os.path.join(nbl_consts['repo_dir'], '3. Select fCpGs', 'outputs', 'NBL_Clock_CpGs.txt'), dtype=str)

    Clock_CpGs_intersect = np.intersect1d(Clock_CpGs, beta_values.index)
    beta_values_Clock = beta_values.loc[Clock_CpGs_intersect]

    if verbose:
        print('DONE')

    ####################################################################################
    ####################################################################################
    
    # Save Clock beta values file
    if verbose:
        print('Saving Clock beta values file...', end='\n\n')
        time.sleep(1)


    nbl_util.writeWithoutOverwrite(
        filepath=os.path.join(nbl_consts['official_indir'], 'Henrich', 'Henrich.methyl.TARGET_Clock_sites.tsv'),
        data=beta_values_Clock,
        writeFunc=lambda filepath, df: df.to_csv(filepath, sep='\t'),
        readFunc=lambda filepath: pd.read_table(filepath, index_col=0),
        compareFunc=lambda data, existing_data:(data == existing_data).all().all(),
        verbose=True)

    if verbose:
        print('\nDONE')

    ####################################################################################
    ####################################################################################

    # Load clinical table
    if verbose:
        print('Loading and processing sample annotations...', end=' ')
        time.sleep(1)

    
    sample_annotations = pd.read_table(os.path.join(proj_dir, 'sample_annotations.txt'), index_col=0).T

    # Process all columns
    col_list = []
    for i in range(sample_annotations.shape[1]):
        col = sample_annotations.iloc[:, i]
        if col.name == 'Sample_characteristics_ch1':
            col_split = col.str.split(': ', expand=True)
            assert col_split[0].unique().shape[0] == 1
            col_list.append(col_split[1].rename(col_split.iloc[0, 0]))
        else:
            col_list.append(col)
    clinical = pd.concat(col_list, axis=1).set_index('Sample_geo_accession').drop('Sample_supplementary_file', axis=1)
    # clinical = pd.read_table(
    #     os.path.join(nbl_consts['official_indir'], 'Henrich', 'sample_annotations_clean.txt'), index_col='Sample_geo_accession'
    # )

    # Reformat strings
    clinical['inss stage'] = 'Stage ' + clinical['inss stage'].astype(str)
    clinical['inss stage'] = clinical['inss stage'].map(lambda x:'Stage 2' if x in['Stage 2.1', 'Stage 2.2'] else x)
    clinical['age at diagnosis'] = clinical['age at diagnosis'].map({'>=1.5 years':'Age >= 1.5 years', '<1.5 years':'Age < 1.5 years'})
    clinical['current risk category'] = clinical['current risk category'].map(lambda x:x.replace('-', ' ').capitalize())
    clinical['Stage grouped'] = clinical['inss stage'].map(nbl_util.stage_mapper)
    clinical['Risk grouped'] = clinical['current risk category'].map(nbl_util.risk_mapper)
    
    for col in ['Sample_channel_count', 'Sample_taxid_ch1', 'Sample_contact_zip/postal_code', 'Sample_data_row_count']:
        clinical[col] = clinical[col].astype(int)
    
    # Add dummy flag
    clinical['in_analysis_dataset'] = True

    if verbose:
        print('DONE')

    
    ####################################################################################
    ####################################################################################
    
    # Save processd clinical data
    if verbose:
        print('\nSave processed clinical data...', end=' ')
        time.sleep(1)

    nbl_util.writeWithoutOverwrite(
        os.path.join(proj_dir, 'clinical.processed.tsv'),
        clinical,
        lambda filepath, data: data.to_csv(filepath, sep='\t'),
        lambda filepath: pd.read_table(filepath).set_index('Sample_geo_accession'),
        compareFunc=lambda old, new: (old.round(5).fillna(0) == new.round(5).fillna(0)).all(axis=None, skipna=True)
    )

    if verbose:
        print('DONE')


