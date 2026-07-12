import pandas as pd
import numpy as np
import os
import time
import EpiClockNBL.util as nbl_util
nbl_consts = nbl_util.consts


# IO directories
proj_dir = os.path.join(nbl_consts['official_indir'], 'TARGET')

def pipeline(verbose=True):

    if verbose:
        print(f'\nRunning TARGET data processing pipeline with verbose={verbose}.')
        time.sleep(1)
        print('#'*100)
        print('#'*100)
        time.sleep(1)

    ####################################################################################
    ####################################################################################

    # Import data
    if verbose:
        print('\nImporting DNA methylation and clinical data...', end=' ')
        time.sleep(1)

    beta_values = pd.read_table(os.path.join(proj_dir, 'cohort1.methyl.tsv'))
    clinical = pd.read_table(os.path.join(proj_dir, 'cohort1.clinical.tsv'))
    clinical = clinical.set_index('submitter_id')

    if verbose:
        print('DONE')

    ####################################################################################
    ####################################################################################
    
    # Calculate LUMP purity values
    if verbose:
        print('\nCalculating LUMP purity values...', end=' ')
        time.sleep(1)

    lump_purity = nbl_util.getLUMP_values(beta_values)
    lump_purity.index = lump_purity.index.to_series().apply(nbl_util.getSampleID).values

    if verbose:
        print('DONE')

    ####################################################################################
    ####################################################################################
    
    # Process clinical data
    if verbose:
        print('\nProcessing clinical data...', end='\n')
        time.sleep(1)

    # Add sample IDs
    sample_ids = lump_purity.index.values
    patient_to_sample_ids = pd.Series(data=sample_ids, index=[nbl_util.sampleToPatientID(x) for x in sample_ids]).rename('sampleID')
    clinical = clinical.merge(patient_to_sample_ids, left_index=True, right_index=True, how='left', validate='one_to_one')
    clinical = clinical.dropna(subset=['sampleID'])

    # Create analysis filter: LUMP >= 0.7 and Neuroblastoma diagnosis
    clinical = clinical.merge(lump_purity.rename('LUMP'), left_on='sampleID', right_index=True, how='inner', validate='one_to_one')
    clinical['LUMP_pure'] = clinical['LUMP'] >= 0.7
    clinical['Neuroblastoma_diagnosis'] = clinical['primary_diagnosis'] == 'Neuroblastoma, NOS'
    clinical['in_analysis_dataset'] = clinical['LUMP_pure'] & clinical['Neuroblastoma_diagnosis']

    print(clinical['primary_diagnosis'].value_counts())
    print(f'{(~clinical['Neuroblastoma_diagnosis']).sum()} tumors were ganglioneuroblastomas')
    print(f'{(clinical['Neuroblastoma_diagnosis'] & ~clinical['LUMP_pure']).sum()} neuroblastoma tumors had LUMP < 0.7')

    # Process certain columns
    clinical['Age'] = clinical['age_at_diagnosis'] / 365
    # assert not clinical['Age'].astype(float).isna().any()
    clinical['Age >= 1.5'] = (clinical['Age'] >= 1.5).map(lambda x:'Age >= 1.5 years' if x else 'Age < 1.5 years')
    # clinical['Age >= 1.5'] = clinical['Age'].map(ageStratify)
    clinical['cog_neuroblastoma_risk_group'] = clinical['cog_neuroblastoma_risk_group'].map(lambda x:x.capitalize() if type(x) is str else x)

    if verbose:
        print('\nDONE')


    ####################################################################################
    ####################################################################################
    
    # Process clinical data
    if verbose:
        print('\nRetrieving and joining in MYCN amplication data from supplementary clinical data...', end=' ')
        time.sleep(1)

    supp_tbl = pd.read_excel(
        os.path.join(proj_dir, 'gdc_data', 'TARGET-NBL', 'Clinical', 'Clinical_Supplement', 'f5bfc8a5-b903-4418-b633-62234388a635', 'TARGET_NBL_ClinicalData_Discovery_20230523.xlsx'),
        sheet_name='Clinical Data', index_col=0
        )

    clinical = clinical.drop('MYCN status', axis=1, errors='ignore')
    clinical = clinical.merge(supp_tbl['MYCN status'], left_index=True, right_index=True, validate='one_to_one')
    clinical['MYCN status'] = clinical['MYCN status'].map({'Amplified':'MYCN-amplified', 'Not Amplified':'MYCN-nonamplified'})
    clinical['Stage grouped'] = clinical['inss_stage'].map(nbl_util.stage_mapper)
    clinical['Risk grouped'] = clinical['cog_neuroblastoma_risk_group'].map(nbl_util.risk_mapper)

    if verbose:
        print('DONE')

    ####################################################################################
    ####################################################################################
    
    # Process clinical data
    if verbose:
        print('\nRetrieving and joining in follow-up (survival) data from supplementary clinical data...', end=' ')
        time.sleep(1)

    # follow_up = pd.read_table(os.path.join(proj_dir, 'follow_up.tsv'), na_values=["'--"])
    # follow_up_selected = follow_up.loc[follow_up['case_submitter_id'].isin(clinical.index) & (follow_up['timepoint_category']=='Last Contact'),
    #                                ['case_submitter_id', 'days_to_follow_up', 'year_of_follow_up']].set_index('case_submitter_id')
    # clinical = clinical.merge(follow_up_selected, left_index=True, right_index=True, how='left', validate='one_to_one')

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
        lambda filepath, data: data.reset_index(names='submitter_id').to_csv(filepath, sep='\t', index=False),
        lambda filepath: pd.read_table(filepath).set_index('submitter_id'),
        compareFunc=lambda old, new: (old.round(5).fillna(0) == new.round(5).fillna(0)).all(axis=None, skipna=True)
    )

    if verbose:
        print('DONE')


