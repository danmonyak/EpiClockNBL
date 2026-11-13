"""
Source code for building the blood cell-type methylation lookup table used for the purity adjustment

The table holds, for each CpG site, the average beta value of each purified blood cell type
in GEO series GSE35069 (Reinius et al., 2012, PLoS One 7:e41361): 6 healthy male donors,
with whole blood, PBMC, granulocytes and seven isolated cell types measured on the 450K array.

"""

import os
import gzip
import time
import urllib.request
import pandas as pd
import EpiClockNBL.util as nbl_util
nbl_consts = nbl_util.consts

# IO directories
other_dir = os.path.join(nbl_consts['official_indir'], 'Other')

SERIES_MATRIX_URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE35nnn/GSE35069/matrix/GSE35069_series_matrix.txt.gz'
SERIES_MATRIX_FILENAME = 'GSE35069_series_matrix.txt.gz'
OUTPUT_FILENAME = 'GSE35069_cell_type_mean_beta.tsv'

# Columns of the lookup table, named by the sample source names in the series matrix
CELL_TYPES = ['Whole blood', 'PBMC', 'Neutrophils', 'Granulocytes',
              'CD4+ T cells', 'CD8+ T cells', 'CD14+ Monocytes', 'CD19+ B cells',
              'CD56+ NK cells', 'Eosinophils']
N_DONORS = 6
N_PROBES = 485_577


def pipeline(verbose=True):

    if verbose:
        print(f'Running wbc_lookup pipeline with verbose={verbose}.')
        time.sleep(1)
        print('#'*100)
        print('#'*100)
        time.sleep(1)

    os.makedirs(other_dir, exist_ok=True)
    series_matrix_path = os.path.join(other_dir, SERIES_MATRIX_FILENAME)

    ####################################################################################
    ####################################################################################

    # Download the series matrix file
    if os.path.exists(series_matrix_path):
        if verbose:
            print('Series matrix file already downloaded.')
    else:
        if verbose:
            print('Downloading series matrix file (about 120 MB)...', end=' ')
        urllib.request.urlretrieve(SERIES_MATRIX_URL, series_matrix_path)
        if verbose:
            print('DONE')

    ####################################################################################
    ####################################################################################

    # Read the sample cell types and the beta values
    if verbose:
        print('Reading series matrix file...', end=' ')

    source_names = None
    with gzip.open(series_matrix_path, 'rt', encoding='utf-8') as f:
        # Header lines, until the start of the data table
        for line in f:
            if line.startswith('!Sample_source_name_ch1'):
                source_names = [s.strip('"') for s in line.rstrip('\n').split('\t')[1:]]
            if line.startswith('!series_matrix_table_begin'):
                break
        # Data table, with one column of beta values per sample
        beta_values = pd.read_table(f, index_col=0, low_memory=False)

    # The last line marks the end of the table
    beta_values = beta_values.loc[beta_values.index != '!series_matrix_table_end'].apply(pd.to_numeric)
    cell_types = pd.Series(source_names, index=beta_values.columns)

    assert beta_values.shape == (N_PROBES, len(CELL_TYPES) * N_DONORS)
    assert (cell_types.value_counts().reindex(CELL_TYPES) == N_DONORS).all()

    if verbose:
        print('DONE')

    ####################################################################################
    ####################################################################################

    # Average the beta values of each cell type across the donors
    # Missing values are skipped, so the average uses the donors measured at that site
    if verbose:
        print('Averaging beta values of each cell type...', end=' ')

    lookup = pd.DataFrame({
        cell_type: beta_values.loc[:, cell_types == cell_type].mean(axis=1, skipna=True)
        for cell_type in CELL_TYPES
    })
    # Average of the two neutrophil-rich populations
    lookup['neutrophils ave'] = (lookup['Neutrophils'] + lookup['Granulocytes']) / 2
    lookup.index.name = 'ID_REF'

    if verbose:
        print('DONE')

    ####################################################################################
    ####################################################################################

    # Save the lookup table, don't overwrite an existing file with different values
    if verbose:
        print('Saving lookup table...', end='\n\n')

    nbl_util.writeWithoutOverwrite(
        os.path.join(other_dir, OUTPUT_FILENAME),
        lookup,
        lambda filepath, data: data.to_csv(filepath, sep='\t'),
        lambda filepath: pd.read_table(filepath, index_col=0),
        compareFunc=lambda new, old: (old.round(10).fillna(-1) == new.round(10).fillna(-1)).all(axis=None)
    )

    if verbose:
        print('\nDONE')
