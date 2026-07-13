import numpy as np
import os

def arr2dToString(arr):
    return '\n'.join(['\t'.join(map(lambda x: format(x, '.18e'), x)) for x in arr]) + '\n'

def writeBetaValues(beta_list, beta_values_outfilepath, should_clear=True):
    if len(beta_list) == 0:
        return
    with open(beta_values_outfilepath, 'a') as f:
        beta_arr = np.stack(beta_list, axis=0)
        beta_arr_str = arr2dToString(beta_arr)
        f.writelines(beta_arr_str)
    if should_clear:
        beta_list.clear()

def writeNcells(n_cells_list, n_cells_outfilepath, should_clear=True):
    if len(n_cells_list) == 0:
        return
    with open(n_cells_outfilepath, 'a') as f:
        n_cells_arr = np.stack(n_cells_list, axis=0)
        n_cells_arr_str = arr2dToString(n_cells_arr)
        f.writelines(n_cells_arr_str)
        # n_cells_list_str = '\n'.join(map(str, n_cells_list)) + '\n'
        # f.writelines(n_cells_list_str)
    if should_clear:
        n_cells_list.clear()

def resetObjsDeleteFiles(ensmbl, lists, files):
    ensmbl.reInit()
    for l in lists:
        l.clear()
    for fi in files:
        if os.path.exists(fi):
            os.remove(fi)

def writeLine(filepath, text):
    with open(filepath, 'a') as f:
        f.write(text + '\n')