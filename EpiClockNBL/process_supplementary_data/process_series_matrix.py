"""
Source code for processing the Henrich series matrix file

"""

import os
import time
from tqdm import tqdm
import EpiClockNBL.util as nbl_util
nbl_consts = nbl_util.consts

# IO directories
proj_dir = os.path.join(nbl_consts['official_indir'], 'Henrich')

def pipeline(verbose=True):
    
    if verbose:
        print(f'Running process_series_matrix pipeline with verbose={verbose}.')
        time.sleep(1)
        print('#'*100)
        print('#'*100)
        time.sleep(1)

    ####################################################################################
    ####################################################################################
    
    # Load beta values
    if verbose:
        print('Loading and processing series matrix file...', end='\n\n')
        time.sleep(1)

    
    header_lines = []
    beta_values_lines = []

    reading_header = False
    i = float('-Inf')
    with open(os.path.join(proj_dir, 'GSE73515_series_matrix.txt'), 'r') as f:
        for l in tqdm(f, total=382_261):

            if (i < 0) and (not reading_header) and l.startswith('!Sample_title'):
                reading_header = True
            elif (i < 0) and reading_header and l.startswith('!series_matrix_table_begin'):
                reading_header = False
                i = 0

            if reading_header:
                header_lines.append(
                    l.lstrip('!')
                )
            elif i > 0:
                beta_values_lines.append(
                    l
                )

            i += 1

    beta_values_lines.pop()

    if verbose:
        print(f'Read {len(header_lines)} lines from header, {len(beta_values_lines):,} lines from DNAm data section.')
        print('\nDONE')

    
    ####################################################################################
    ####################################################################################
    
    # Select subset with Clock CpGs
    if verbose:
        print('Writing header and DNAm data to files...', end=' ')
        time.sleep(1)

    with open(os.path.join(proj_dir, 'sample_annotations.txt'), 'w') as f:
        f.write(''.join(header_lines))

    with open(os.path.join(proj_dir, 'beta_values.txt'), 'w') as f:
        f.write(''.join(beta_values_lines))

    if verbose:
        print('DONE')



