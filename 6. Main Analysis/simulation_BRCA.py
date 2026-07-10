"""
simulation.py
=======
Author - Daniel Monyak
9-3-24
=======

Source code for simulating fCpG sites in a growing population of cancer cells
Population is modeled as a 2d array of states:
    # cells x # CpGs
    4 possible states (two alleles)
Time is handled discretely
    Pass one day at a time
    Exponential RVs (time to next flip/division/death) are approximated as Bernoulli RVs (probability event occurs in a day)
    Expo(theta) ----> Bernouli(theta)
    
"""

from time import process_time, sleep, time
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import random
from scipy.stats import linregress
import datetime
import pickle
import os
import json
from EpiClockInvasiveBRCA.src.simulation_util import *

MAX_CELLS = int(1e8)


def betasFromStates(state_arr):
    """
    Returns array of beta values
    
    Parameters
    ----------
    state_arr : 2d array of states of cells x sites
    
    Returns
    -------
    beta_values : ndarray of beta values of fCpGs
    
    """
    meth_alleles = np.ceil(state_arr / 2)
    beta_values = (meth_alleles/2).mean(axis=0)
    return beta_values

def betaCorr(beta1, beta2):
    """
    Returns correlation of beta values
    
    Parameters
    ----------
    beta1, beta2 : arrays of beta values
    
    Returns
    -------
    rvalue : R-value of Pearson Correlation between beta1 and beta2 arrays
    
    """
    rvalue = linregress(beta1, beta2).rvalue
    return rvalue

    

def loadFromDirectory(indir=None, suffix=None):
    if indir is None:
        if suffix is None:
            raise Exception('Either indir or suffix must be passed...')

        indir = f'simulation_state_{suffix}'
        
    with open(os.path.join(indir, 'obj_info.json'), 'r') as f:
        obj_info = json.load(f)

    
    with open(os.path.join(indir, 'gen_pickled.bin'), 'rb') as f:
        obj_info['gen'] = pickle.load(f)
    obj_info['re_init'] = False

    print('loading')
    state_arr_alive = np.load(os.path.join(indir, 'state_arr.npy'))
    print('done loading')
    n_cells_alive, n_CpGs = state_arr_alive.shape

    print('reallocating memory')
    state_arr_zero = np.zeros([obj_info['max_cells'] - n_cells_alive, n_CpGs], dtype='byte')
    state_arr = np.concatenate([state_arr_alive, state_arr_zero], axis=0)
    assert (state_arr.shape[0] == obj_info['max_cells']) and (state_arr.shape[1] == n_CpGs)
    
    available_cells = list(range(obj_info['max_cells']))[::-1]
    living_cells = [available_cells.pop() for i in range(n_cells_alive)]        # Initialize with one cell

    new_ens = Ensemble(**obj_info)
    new_ens.available_cells = available_cells
    new_ens.living_cells = living_cells
    new_ens.state_arr = state_arr

    return new_ens

class EnsembleContainer:
    def __init__(self):
        self.ensmbl_list = []
    def addEnsemble(self, ensmbl):
        self.ensmbl_list.append(ensmbl)
    def getCombinedActiveStateArr(self):
        return np.concatenate([ensmbl.state_arr[ensmbl.living_cells] for ensmbl in self.ensmbl_list], axis=0)
    def getBetaValues(self, time_obj):
        time_obj['getBetaValues - line1'] -= time()
        combined_state_arr = self.getCombinedActiveStateArr()
        time_obj['getBetaValues - line1'] += time()
        time_obj['getBetaValues - line2'] -= time()
        betas = betasFromStates(combined_state_arr)
        time_obj['getBetaValues - line2'] += time()
        return betas
    def getNumCells(self):
        return sum([ensmbl.getNumCells() for ensmbl in self.ensmbl_list])
    def passDay(self):
        new_ensmbls = []
        del_ensmbls_idxs = []
        for i, ensmbl in enumerate(self.ensmbl_list):
            response = ensmbl.passDay()
            if response and (response['result'] == 'success'): # passDay was successful
                pass
            elif response and (response['result'] == 'split'):
                new_ensmbls.extend(response['data'])
                del_ensmbls_idxs.append(i)
            else:
                return False

        if len(del_ensmbls_idxs) > 0:
            for i, ensmbl in enumerate(self.ensmbl_list):
                if i not in del_ensmbls_idxs:
                    new_ensmbls.append(ensmbl)
            self.ensmbl_list = new_ensmbls
        
        return True
    
    def reInit(self):
        self.ensmbl_list = self.ensmbl_list[:1]
        self.ensmbl_list[0].reInit()
    def atCapacity(self):
        return any([ensmbl.atCapacity() for ensmbl in self.ensmbl_list])
    def getNumEnsembles(self):
        return len(self.ensmbl_list)


class Ensemble:
    """
    Contains 2d array of states
    Keeps pointers to
        initial parameters
        random generator object
    Represents an ensemble of fCpG sites in a growing population of cells
    
    Attributes
    ----------
    init_params : dict
        Initial parameters of simulation
    gen : Numpy random Generator
        Used to simulate random events
    max_cells: int
        Maximum number of cells in population
    available_cells: list
        List of cells (represented as ints) that are available for use
    living_cells: list
        List of living cells
    state_arr: 2d ndarray
        # cells x # CpGs
        4 possible states
            0: both alleles unmethylated
            1: first allele methylated, second allele unmethylated
            2: first allele unmethylated, second allele methylated
            3: both alleles methylated
            
            These state definitions derive their binary representations
                0 = 00b
                1 = 01b
                2 = 10b
                3 = 11b
                
                0 = unmethylated, 1 = methylated
                
    at_capacity: bool
        All possible cells have been exhausted
        
    Methods
    -------
    
    """
    
    def __init__(self, init_params, gen, max_cells=MAX_CELLS, re_init=True, split=False, split_limit=None, n_split=5, day=0,
                 save_info=None):
        """
        Parameters
        ----------
        init_params : dict
            Initial parameters of simulation
        gen : Numpy random Generator
            Used to simulate random events
        max_cells: int
            Maximum number of cells in population
            
        """
        
        self.init_params = init_params
        self.gen = gen
        if split_limit is not None:
            self.max_cells = split_limit * 5
        else:
            self.max_cells = max_cells
        
        # self.available_cells = list(range(self.max_cells))[::-1]
        # self.living_cells = [self.available_cells.pop()]        # Initialize with one cell
        
        # Determine # CpG sites
        if 'n_CpGs' in self.init_params:
            n_CpGs = self.init_params['n_CpGs']
        else:
            n_CpGs = sum(self.init_params['init_site_state_counts'])
        
        
        self.state_arr = np.zeros([self.max_cells, n_CpGs], dtype='byte')
        print('Allocated memory!')

        self.day = day
        if re_init:
            self.reInit()
        self.at_capacity = False

        ################################################################
        ################################################################
        
        self.split = split
        self.split_limit = split_limit
        self.n_split = n_split
        self.save_info = save_info

    def returnEnsembleSplits(self, split_ens_max_cells=None, save_not_return=False):
        n_cells = self.getNumCells()
        
        slice_idxs = list(map(int, np.linspace(0, n_cells, self.n_split + 1)))
        # print(f'slice_idxs: {slice_idxs}')
        new_ens_list = []
        for i in range(len(slice_idxs) - 1):
            start = slice_idxs[i]
            stop = slice_idxs[i+1]
            n_cells_new = stop - start
            # print(l[slice(start, stop)])

            if n_cells_new == 0:
                continue

            new_gen = np.random.default_rng(int(1000*self.gen.random()))
            if split_ens_max_cells is None:
                if save_not_return:
                    split_ens_max_cells = n_cells_new
                else:
                    split_ens_max_cells = self.max_cells
            new_ens = Ensemble(self.init_params, new_gen, split_ens_max_cells,
                               re_init=False,
                               split=True,
                               split_limit=self.split_limit, 
                               n_split=self.n_split,
                               day=self.day
                              )
            new_ens.available_cells = list(range(new_ens.max_cells))[::-1]
            new_ens.living_cells = [new_ens.available_cells.pop() for k in range(n_cells_new)]
            new_ens.state_arr[new_ens.living_cells] = self.state_arr[self.living_cells[slice(start, stop)]]

            # 
            if save_not_return:
                split_name = f'split_{i}'
                new_ens.saveAsDirectory(
                    parent_outdir=self.save_info['saveAsDirectory_parent_outdir'],
                    outdir=split_name
                )
                writeLine(self.save_info['split_jobs_filepath'], split_name)
                del new_ens
            else:
                new_ens_list.append(new_ens)

        if not save_not_return:
            return new_ens_list

    def saveAsDirectory(self, suffix=None, parent_outdir=None, outdir=None):
        if outdir is None:
            if suffix is None:
                now = datetime.datetime.now()
                suffix = now.strftime("%m-%d-%Y-%H-%M-%S")
        
            outdir = f'simulation_state_{suffix}'

        if parent_outdir is not None:
            os.makedirs(parent_outdir, exist_ok=True)
            outdir = os.path.join(parent_outdir, outdir)
        
        if os.path.exists(outdir):
            for fi in os.listdir(outdir):
                os.remove(os.path.join(outdir, fi))
        else:
            os.mkdir(outdir)
        
        # Directory for simple fields
        obj_info = {}
        for attr in ['init_params', 'max_cells', 'split', 'split_limit', 'n_split', 'day']:
            obj_info[attr] = getattr(self, attr)
        
        with open(os.path.join(outdir, 'obj_info.json'), 'w') as f:
            f.write(json.dumps(obj_info))
        
        gen_pickled = pickle.dumps(self.gen, protocol=pickle.HIGHEST_PROTOCOL)
        with open(os.path.join(outdir, 'gen_pickled.bin'), 'wb') as f:
            f.write(gen_pickled)
        
        np.save(os.path.join(outdir, 'state_arr.npy'), self.state_arr[self.living_cells])

    
    def copy(self, new_gen):
        """
        Create copy of fCpG ensemble
        Requires passing of a new generator object, so that copied ensembles do not generate the same random events
        
        Parameters
        ----------
        new_gen : Numpy random Generator
            New generator object
            
        Returns
        -------
        new_ens : Ensemble
            Copy of this Ensemble object
            
        """
        
        new_ens = Ensemble(self.init_params, new_gen, self.max_cells)
        new_ens.available_cells = self.available_cells.copy()
        new_ens.living_cells = self.living_cells.copy()
        new_ens.state_arr = self.state_arr.copy()
        return new_ens
    
    def getRandCell(self, max_cells=None):
        """
        Create Ensemble representing a metastasis of the current population of cells
        Create a new Ensemble and copy the state of a randomly picked living cell
        
        Parameters
        ----------
        max_cells : int
            maximum # cells in new Ensemble
            
        Returns
        -------
        new_ens : Ensemble
            New metastasis Ensemble
            
        """
        if max_cells is None:
            max_cells = self.max_cells
        new_ens = Ensemble(self.init_params, np.random.default_rng(), max_cells)
        cell_i = np.random.choice(self.living_cells, replace=False)
        new_ens.state_arr[new_ens.living_cells[0]] = self.state_arr[cell_i]
        return new_ens
    
    def getNumCells(self):
        """
        Return the number of living cells
        """
        return len(self.living_cells)
    
    def atCapacity(self):
        """
        Return True iff Ensemble has reached max number of living cells
        """
        return self.at_capacity
    
    def reInit(self):
        """
        Initialize or reinitialize (if all cells died) Ensemble 
        Only initialize one cells (that's how many are alive)
        
        If n_CpGs and init_site_state_probs are both in self.init_params
            initialize first cell with that number of CpGs and randomly distribute it among
            the 4 states according to the probabilities in init_site_state_probs
        If not, but init_site_state_counts is in self.init_params
            init_site_state_counts holds desired number of sites in each state in first cell
        """
        self.day = 0
        self.available_cells = list(range(self.max_cells))[::-1]
        self.living_cells = [self.available_cells.pop()]        # Initialize with one cell
        
        if ('n_CpGs' in self.init_params) and ('init_site_state_probs' in self.init_params):
            self.state_arr[self.living_cells[0]] = self.gen.choice(4, p=self.init_params['init_site_state_probs'], size=[1, self.init_params['n_CpGs']])
        elif 'init_site_state_counts' in self.init_params:
            count = 0
            for st in range(4):
                next_count = count + self.init_params['init_site_state_counts'][st]
                self.state_arr[self.living_cells[0], count:next_count] = st
                count = next_count
#             self.state_arr = np.reshape(np.concatenate([[st for i in range(self.init_params['init_site_state_counts'][st])] for st in range(4)], dtype=int, casting='unsafe'), [1, -1])
        else:
            sys.exit('Must provide a init_site_state_counts argument...')
    
    def getBetaValues(self):
        """
        Return beta values of all fCpGs, using all living cells
            
        Returns
        -------
        beta_values : ndarray
            Beta values of all fCpG sites
            
        """
        beta_values = betasFromStates(self.state_arr[self.living_cells])
        return beta_values
    
    def getRandCellsBetaValues(self, n=1):
        """
        Return beta values of n randomly picked living cells
            
        Returns
        -------
        beta_value_list : List of ndarrays
            Each ndarray is the array of beta values of a single cell
            
        """
        
        n = min(n, self.getNumCells())
        beta_value_list = []
        for cell in np.random.choice(self.living_cells, size=n, replace=False):
            beta_value_list.append(betasFromStates(np.reshape(self.state_arr[cell], [1, -1])))
        
        return beta_value_list
            
    def passDay(self, die=False):
        """
        Pass one time unit (one day)
        
        Pick number of each event from multinomial model
            Cumulative probabilities for events
            First range - divide
            Second range - die
            Third range - nothing
        If tumor is eliminated, reset the Ensemble and return False
        Dying cells - return cells to self.available_cells
        Dividing cells
            - pick cells from self.available_cells and copy states
            - determine new state of each new cell using flip rate to determine probabilities
        
        Parameters
        ----------
        die : bool
            if True, it forces the Ensemble has to die
        
        Returns
        -------
           bool : False if the tumor was eliminated or it ran out of available cells, else True
        
        """

        # Split cell population for efficiency
        if self.split and (self.getNumCells() > self.split_limit):
            save_not_return = self.save_info is not None
            return {'result':'split', 'data':self.returnEnsembleSplits(save_not_return=save_not_return)}
            
        ###########################################################################
        ###########################################################################
        ###########################################################################
            
        
        n_divide, n_die, n_nothing = self.gen.multinomial(self.getNumCells(), pvals=[self.init_params['prolif_rate'], self.init_params['death_rate'], 1 - self.init_params['prolif_rate'] - self.init_params['death_rate']])
        
        if die or (n_divide + n_nothing == 0):   # Tumor was eliminated
            print('Tumor died')
            self.available_cells.extend(self.living_cells)
            self.living_cells = [self.available_cells.pop()]
            return False
        
        # Shuffle order of living cells so that the same cells don't live forever
        self.gen.shuffle(self.living_cells)
        
        # Return dead cells to self.available_cells
        self.available_cells.extend(self.living_cells[n_divide + n_nothing:])
        self.living_cells = self.living_cells[:n_divide + n_nothing]
        
        # If the Ensemble will run out of cells
        # Reset cell lists and return False
        if n_divide > len(self.available_cells):
            print('Simulation reached capacity')
            self.at_capacity = True
            self.available_cells.extend(self.living_cells)
            self.living_cells = [self.available_cells.pop()]
            return False
        
        # Create new cells from divisions
        new_cell_states = self.state_arr[self.living_cells[:n_divide]]
        new_cell_idxs = [self.available_cells.pop() for i in range(n_divide)]
        self.living_cells.extend(new_cell_idxs)
        
        # Probability of each event
        #       Neither allele flips, first allele flips, second allele flips, both alleles flip
        mu = self.init_params['flip_rate']
        pvals = [(1 - mu)**2, mu*(1 - mu), mu*(1 - mu), mu**2]
        
        # Holds each state change for each new cell
        flip_event_arr = self.gen.choice(4, size=new_cell_states.shape, replace=True, p=pvals)
        
        # Bitwse XOR operation correctly performs state transitions
        #      0 (unmethylated) ^ 0 (no flip) = 0 (unmethylated)
        #      0 (unmethylated) ^ 1 (flip) = 1 (methylated)
        #      1 (methylated) ^ 0 (no flip) = 1 (methylated)
        #      1 (methylated) ^ 1 (flip) = 0 (unmethylated)
        self.state_arr[new_cell_idxs] = new_cell_states ^ flip_event_arr
#         new_cell_states ^ flip_event_arr
#         self.state_arr = np.vstack([self.state_arr, new_cell_states])
    
        # return True    # tumor is still alive
        self.day += 1
        return {'result':'success'}


def plotBetaValues(ax, beta_values=None, binwidth=None, color=None, opacity=None,
                   labelfontsize=None, ticksfontsize=None, sf=1, bins=None):
    """
    Plot histogram of beta values

    Parameters
    ----------
    ax : Matplotlib Axes object
        Axes to be used for plot
    beta_values : ndarray
        Beta values to be plotted
    binwidth : number or pair of numbers
    color : str or tuple of floats
    opacity : float
    labelfontsize : float
    ticksfontsize : float
    sf : float
        scale factor of figure
    bins : str, number, vector, or a pair of such values
        See Seaborn histplot documentation

    Returns
    -------
       bool : False if the tumor was eliminated or it ran out of available cells, else True

    """
    
    sns.histplot(ax=ax, x=beta_values, bins=bins, stat='proportion', color=color, alpha=opacity)
    ax.set_xlabel('β', fontsize=labelfontsize * sf)
    ax.set_ylabel('Proportion', fontsize=labelfontsize * sf)
    ax.set_xticks([0, 0.5, 1])
    ax.set_xlim(-0.05, 1.05)
    ax.tick_params(axis='both', labelsize=ticksfontsize * sf, width=sf, length=8 * sf)
