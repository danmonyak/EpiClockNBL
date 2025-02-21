"""
consts.py
=======
Author - Daniel Monyak
6-14-24
=======


Provides
    Directories and variables used by other notebooks and source files

"""

import os
import sys

consts = {}

#! Set to parent directory of repo
labdir = os.getenv('lab')
if labdir is None:
    sys.exit('Must set $lab environmental variable to parent directory of MolecularClocks repository...')
    
consts['boxdir'] = os.path.join(os.getenv("HOME"), 'Library/CloudStorage/Box-Box/PROJECT 06023: MolClocks/Neuroblastoma_Paper')

#! Set to directory with data
consts['official_indir'] = os.path.join(os.getenv("HOME"), 'Library/CloudStorage/Box-Box/PROJECT 06023: MolClocks/Neuroblastoma_Paper/Datasets')
consts['repo_dir'] = os.path.join(labdir, 'PanCancerClock')

consts['CPE_threshold'] = 0.6
consts['LUMP_threshold'] = 0.7

consts['palette_jco'] = ["#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#7AA6DCFF", "#003C67FF", "#8F7700FF", "#3B3B3BFF", "#A73030FF", "#4A6990FF"]


