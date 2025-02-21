import numpy as np
import pandas as pd
import os
from math import ceil, isnan
from itertools import product, accumulate
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress, ranksums
import json

subdir_list = os.getcwd().split(os.sep)
while True:
    if subdir_list[-1] == 'EpiClockNBL':
        break
    
    subdir_list.pop()
        
repo_dir = os.path.join(os.sep, *subdir_list)

consts = json.loads(''.join(open(os.path.join(repo_dir, 'src', 'consts.json'), 'r').readlines()))
# consts = json.loads(''.join(open(os.path.join('/Users/danielmonyak/Documents/Duke/lab/EpiClockInvasiveBRCA', 'src', 'consts.json'), 'r').readlines()))

sampleToPatientID = lambda x: '-'.join(x.split('-')[:3])
getSampleID = lambda x: '-'.join(x.split('-')[:4])
isNaVec = np.vectorize(lambda x:(x is None) or ((type(x) is not str) and isnan(x)))

def combineFilters(filters):
    return list(accumulate(filters, lambda x,y:x&y))[-1]

def pearsonCorrelation(ser1, ser2, get_n_used=False):
    use_mask = ~(ser1.isna() | ser2.isna())
    res = linregress(ser1[use_mask], ser2[use_mask])
    if get_n_used:
        return res, use_mask.sum()
    else:
        return res

def wilcoxonRankSums(ser1, ser2, get_n_used=False):
    use_mask = ~(ser1.isna() | ser2.isna())
    res = ranksums(ser1[use_mask], ser2[use_mask])
    if get_n_used:
        return res, use_mask.sum()
    else:
        return res


def saveBoxPlotNew(sample_annotations, var_cat, var_y='c_beta', ax=None, restrict=True, use_groups=None,
                outdir='.', outfile=True, title=False, custom_title=None, xlabel=None, ylabel=None, palette=None,
                plot_ymax_mult=0.25, signif_bar_heights=0.03,
                   signif_fontsize=14, ylim=None,
                   figsize=(10, 10), labelfontsize=20, ticksfontsize=10, linewidth=1, fliersize=1, sf=1):
    """
    Create box plot from a sample annotations DataFrame
    Plot some numerical variable and stratfy by some categorical variable
    Can include significance bars (Wilcoxon rank-sum)
    
    Parameters
    ----------
    sample_annotations : Pandas Dataframe
        has columns var_cat and var_y
    var_cat : str
        categorical variable to stratify by
    var_y : str
        numerical variable to measure
    ax : matplotlib.axes.Axes
        Axes instance used for plot
        if None, it will create a new Axes instance
    restrict : boolean
        True iff we should restrict to samples with in_analysis_dataset==True
    use_groups : list of strs
        specific values of var_cat to compare
        if None, use all unique values of var_cat
    outdir : str
        output directory of figure file generated
    outfile : boolean
        True iff we want to save the figure to a file
    title : boolean
        use title
    custom_title : str
        title to use other than default (sample_annotations.name)
    xlabel : str
        custom x-label
    ylabel : str
        custom y-label
    palette : palette name, list, or dict
        color palette
    plot_ymax_mult : float
        proportionally changes y top limit
    signif_bar_heights : float
        height of vertical part of significance bars
        set to None if you don't want significance bars
    signif_fontsize : float
        fontsize of signficances * symbols
    ylim : tuple of floats
        ylim that is used if not None and signif_bar_heights is None
    figsize : tuple of floats
        figure size
    labelfontsize : float
        fontsize of x-label, y-label, x-ticks (categories), and title
    ticksfontsize : float
        fontsize of y-ticks
    linewidth : float
        linewidth argument passed to sns.boxplot
    fliersize : float
        fliersize argument passed to sns.boxplot
    sf : float
        scale factor
        change to alter the size of a figure - scales everything proportionally
    
    Returns
    -------
    ax : matplotlib.axes.Axes
        Axes instance used for plot

    """
    
    ########################
    ##### Select data
    ########################
    
    # Select samples wth var_cat in use_groups
    if use_groups is None:
        use_groups = sample_annotations[var_cat].unique()
        use_groups = np.sort(use_groups[~isNaVec(use_groups)])
    
    # Boolean mask of samples to use
    use_samples_mask = sample_annotations[var_cat].isin(use_groups)

    # Add to mask and set outfile name
    if restrict:
        use_samples_mask &= sample_annotations['in_analysis_dataset']
        if outfile:
            outfile_name = f'{sample_annotations.name}-{var_cat}-{var_y}-restrict.pdf'
    elif outfile:
        outfile_name = f'{sample_annotations.name}-{var_cat}-{var_y}.pdf'

    # Exclude samples with missing data
    use_samples_mask &= ~sample_annotations[var_y].isna()
    
    # Define final DataFrame for plotting
    ###### DELETE but check first
#     plot_data = sample_annotations.loc[use_samples_mask, [var_cat, var_y]].dropna()
    plot_data = sample_annotations.loc[use_samples_mask, [var_cat, var_y]]
    
    # Sanity check
    assert plot_data.shape[0] == plot_data.dropna().shape[0]
    
    ########################
    ##### Create plot
    ########################
    
    if ax is None:
        fig, ax = plt.subplots(figsize=np.array(figsize) * sf)
    sns.boxplot(ax=ax, data=plot_data, x=var_cat, y=var_y,
                order=use_groups, palette=palette,
                linewidth=linewidth * sf, fliersize=fliersize * sf)

    ######################################################
    ##### Customize plot labels, axes, ticks, ticklabels
    ######################################################
    
    if var_y == 'c_beta':
        ax.set_ylabel('$c_β$', fontsize=labelfontsize * sf)
    elif var_y == 'c_beta_adj1':
        ax.set_ylabel('$c_β^a$', fontsize=labelfontsize * sf)
    elif ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=labelfontsize * sf)

    ax.tick_params(axis='x', labelsize=labelfontsize * sf, width=sf, length=8 * sf)
    ax.tick_params(axis='y', labelsize=ticksfontsize * sf, width=sf, length=8 * sf)

    # Add (n = ...) under each x-tick (category)
    ax.set_xticks(ax.get_xticks(),
                  [group + f'\n(n = {(~(plot_data[var_y].isna()) & (plot_data[var_cat] == group)).sum()})' for group in use_groups])
    
    # Set title
    # Default = sample_annotations.name
    # By default, xlabel is erased from plot
    if title:
#         if label is None:
        ax.set_title(sample_annotations.name, fontsize=labelfontsize * sf)
#         else:
#             ax.set_title(f'{label} ({sample_annotations.name})', fontsize=labelfontsize * sf)
        ax.set_xlabel('', fontsize=labelfontsize * sf)
    elif custom_title is not None:
        ax.set_title(custom_title, fontsize=labelfontsize * sf)
    elif xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=labelfontsize * sf)

    ######################################################
    ##### Add significance bars
    ######################################################
    
    if signif_bar_heights is not None:
        max_y = plot_data[var_y].max()
        min_y = plot_data[var_y].min()

        # Create list of indexes
        idxs_1 = list(range(1, len(use_groups) + 1))
        
        # Create list of combinations of indexes
        pair_list = [(idxs_1[x], idxs_1[x + y]) for y in reversed(idxs_1) for x in range((len(idxs_1) - y))]

        # Hold combinations with a signifcant p-value
        significant_pairs = []
        for combo in pair_list:

            ser1 = plot_data.loc[plot_data[var_cat] == use_groups[combo[0]-1], var_y]
            ser2 = plot_data.loc[plot_data[var_cat] == use_groups[combo[1]-1], var_y]

            pvalue = wilcoxonRankSums(ser1, ser2).pvalue
            if pvalue < 0.05:
                significant_pairs.append([combo, pvalue])

        plot_ymax = plot_ymax_mult * 0.8 * max(1, len(significant_pairs)) * (max_y - min_y) + max_y
        ax.set_ylim(top=plot_ymax)

        # Get the y-axis limits
        bottom, top = ax.get_ylim()
        # top = plot_ymax
        y_range = top - bottom

        # Significance bars
        for i, combo in enumerate(significant_pairs):
            # Columns corresponding to the datasets of interest
            x1 = combo[0][0] - 1
            x2 = combo[0][1] - 1
            
            level = len(significant_pairs) - i
            
            # Plot the bar
            bar_height = max_y + (signif_bar_heights * level)

            bar_tips = bar_height - (y_range * 0.02)
            ax.plot(
                [x1, x1, x2, x2],
                [bar_tips, bar_height, bar_height, bar_tips], lw=linewidth * sf, c='k'
            )
            # Significance level
            pvalue = combo[1]
            if pvalue < 0.001:
                sig_symbol = '***'
            elif pvalue < 0.01:
                sig_symbol = '**'
            elif pvalue < 0.05:
                sig_symbol = '*'
            text_height = bar_height
            ax.text((x1 + x2) * 0.5, text_height, sig_symbol,
                    ha='center', va='bottom', c='k', fontsize=signif_fontsize * sf)
    elif ylim is not None:
        ax.set_ylim(ylim)
    
    # Save plot
    if outfile:
        fig.savefig(os.path.join(outdir, outfile_name), format='pdf', pad_inches=0.1)

    return ax


lump_CpGs = np.loadtxt(os.path.join(consts['repo_datadir'], 'lump-CpGs-44.txt'), dtype=str)
def getLUMP_values(beta_values):
    included_lump = np.intersect1d(lump_CpGs, beta_values.index)
    if included_lump.shape[0] != lump_CpGs.shape[0]:
        print(f'Only {included_lump.shape[0]} LUMP sites were available...')
    raw_LUMPs = (beta_values.loc[included_lump].mean(axis=0) / 0.85).to_frame()
    raw_LUMPs[1] = 1
    return raw_LUMPs.min(axis=1)

def getCorrelation(sample_annotations, var_x, var_y='methStdev', only_pure=False, use_samples=None, get_n_used=False):
    mask = ~sample_annotations[var_x].isna()
    if only_pure:
        mask &= sample_annotations['pure']
    if use_samples is not None:
        mask &= sample_annotations.index.isin(use_samples)
    df = sample_annotations.loc[mask]
    ser1 = df[var_y]
    ser2 = df[var_x]
    return pearsonCorrelation(ser1, ser2, get_n_used)


def saveCorrelationPlot(sample_annotations, var_y, var_x='c_beta', restrict=True, use_samples=None,
                        outdir='.', outfile=True, ax=None, text_x=0, text_y=0,
                        scatter_kws={}, line_kws={}, custom_title=None, xlabel=None, ylabel=None, bbox_dict=None, color='blue',
                       figsize=(10, 10), labelfontsize=20, ticksfontsize=10, s=1, sf=1):
    """
    Create scatter plot with a line of best fit from a sample annotations DataFrame
    Plot two numerical variables against each other
    
    Parameters
    ----------
    sample_annotations : Pandas Dataframe
    var_y : str
        numerical variable to plot on y-axis
    var_x : str
        numerical variable to plot on x-axis
    restrict : boolean
        True iff we should restrict to samples with in_analysis_dataset==True
    use_samples : list or ndarray of strs
        samples to restrict correlation to
    outdir : str
        output directory of figure file generated
    outfile : boolean
        True iff we want to save the figure to a file
    ax : matplotlib.axes.Axes
        Axes instance used for plot
        if None, it will create a new Axes instance
    text_x : float
        x coordinate of "R = ..." text
    text_y : float
        y coordinate of "R = ..." text
    scatter_kws : dict
        scatter plot keyword arguments -- see sns.regplot documentation
    line_kws : dict
        line plot keyword arguments -- see sns.regplot documentation
    xlabel : str
        custom x-label
    ylabel : str
        custom y-label
        default is var_y
    color : matplotlib color
        color of all plot elements
    figsize : tuple of floats
        figure size
    labelfontsize : float
        fontsize of x-label, y-label, and title
    ticksfontsize : float
        fontsize of x and y-ticks
    s : float
        size of scatter points
    sf : float
        scale factor
        change to alter the size of a figure - scales everything proportionally
        
    Notes
    -----
    sample_annotations.name must be set
    i.e. sample_annotations.name = 'dataset'

    """
    
    ########################
    ##### Select data
    ########################
    
    # Boolean mask of samples to use
    use_samples_mask = np.ones(shape=sample_annotations.shape[0], dtype=bool)
    
    # Add to mask and set outfile name
    if restrict:
        use_samples_mask &= sample_annotations['in_analysis_dataset']
        outfile_name = f'{sample_annotations.name}-{var_x}-{var_y}-pure.pdf'
    else:
        outfile_name = f'{sample_annotations.name}-{var_x}-{var_y}.pdf'
    
    # Set use_samples if it is None
    if use_samples is None:
        use_samples = sample_annotations.index[use_samples_mask]
    else:
        use_samples = np.intersect1d(use_samples, sample_annotations.index[use_samples_mask])
    
    # Define final DataFrame for plotting
    plot_data = sample_annotations.loc[use_samples, [var_x, var_y]].dropna()
    
    # Calculate Pearson correlation between variables
    res = getCorrelation(plot_data, var_x=var_x, var_y=var_y)
    
    ########################
    ##### Create plot
    ########################
    
    if ax is None:
        fig, ax = plt.subplots(figsize=np.array(figsize) * sf)
    
    scatter_kws['s'] = s * sf**2
    sns.regplot(ax=ax, data=plot_data, x=var_x, y=var_y, scatter_kws=scatter_kws, color=color,
               line_kws=line_kws)

    ######################################################
    ##### Customize plot labels, axes, ticks, ticklabels
    ######################################################

    # Y-label by default is var_y, split by _ and capitalized
    if ylabel is None:
        ax.set_ylabel(" ".join(var_y.split("_")).capitalize(), fontsize=labelfontsize * sf)
    else:
        ax.set_ylabel(ylabel, fontsize=labelfontsize * sf)

    ax.set_title(sample_annotations.name, fontsize=labelfontsize * sf)
    ax.tick_params(axis='both', labelsize=ticksfontsize * sf, width=sf, length=8 * sf)
    
    # Put (R = ...) on the plot
    ax.text(text_x, text_y, f'R = {res.rvalue:.2f}',
                        ha="center", va="bottom",
                        fontfamily='sans-serif', fontsize=0.8 * labelfontsize * sf, bbox=bbox_dict)
    
    if var_x == 'c_beta':
        ax.set_xlabel('$c_β$', fontsize=labelfontsize * sf)
    elif var_x == 'c_beta_adj1':
        ax.set_xlabel('$c_β^a$', fontsize=labelfontsize * sf)
    elif xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=labelfontsize * sf)
    else:
        ax.set_xlabel(ax.get_xlabel(), fontsize=labelfontsize * sf)
        
    # Save plot
    if outfile:
        fig.savefig(os.path.join(outdir, outfile_name), format='pdf', pad_inches=0.1)
        
        
def plotTumorWise(data, CpG_list=None, sample_type='tumor', sample_list=None, n_samps=30, ncols=3, suptitle='random pick of samples', extra_titles=None, title_formats=None, xlabel='Beta', random_seed=None, data_obj='beta_values_SELECTION',
                 outfile=False, outfile_name=None, outdir='images', choose_random=True, color='blue', ylim=None, bins='auto', figsize=None, text_fontsize=None, ticksfontsize=None, opacity=None, sf=1, tight_layout_pad=1):
    if sample_list is None:
        sample_list = data[sample_type]['pureSamples']
    
    n_samps = min(n_samps, len(sample_list))
    nrows = ceil(n_samps / ncols)
    if CpG_list is None:
        if data_obj == "beta_decomp":
            CpG_list = data.index
        else:
            CpG_list = data[sample_type][data_obj].index
        
    
    if choose_random:
        np.random.seed(random_seed)
        samples_randSamp = np.random.choice(sample_list, n_samps, replace=False)
    else:
        samples_randSamp = sample_list[:n_samps]
    
    if figsize is None:
        fig, axes = plt.subplots(nrows, ncols, figsize=(20 * sf, (3 + 3*nrows) * sf))
    else:
        fig, axes = plt.subplots(nrows, ncols, figsize=np.array(figsize) * sf)
    fig.suptitle(suptitle, y=0.99, fontsize=20, fontweight='bold')
    fig.tight_layout(pad=tight_layout_pad)
    
    for i, samp in enumerate(samples_randSamp):
        col = i % ncols
        if nrows > 1:
            row = i // ncols
            ax = axes[row, col]
        else:
            if ncols == 1:
                ax = axes
            else:
                ax = axes[col]
        
        if type(color) is list:
            cur_color = color[min(i, n_samps-1)]
        else:
            cur_color = color
        
        if data_obj == "beta_decomp":
            plot_data = data.loc[CpG_list, samp]
        else:
            plot_data = data[sample_type][data_obj].loc[CpG_list, samp]
        plot = sns.histplot(ax=ax, data=plot_data,
                            stat='proportion', binrange=(0, 1),
                           color=cur_color, bins=bins, alpha=opacity)
        if extra_titles is not None:
            title = extra_titles[i]
        
        if title_formats is None:
            title = samp
        else:
            if data_obj == "beta_decomp":
                title = f"{samp} {title_formats[samp]}"
            else:
                title = title_formats[i].format(samp)
        
        if sample_type == 'normal':
            title += f', age = {age_mapper[samp]}'
        
        ax.set_title(title, fontsize=text_fontsize * sf)
        ax.set_xlabel(xlabel, fontsize=text_fontsize * sf)
        if col == 0:
            ax.set_ylabel('Proportion', fontsize=text_fontsize * sf)
        else:
            ax.set_ylabel('')
        ax.tick_params(axis='both', labelsize=ticksfontsize * sf, width=sf, length=8 * sf)
        
        if ylim is not None:
            ax.set_ylim(ylim)
    
    if not outfile and not ncols == 1:
        fig.show()
    elif outfile_name is None:
        print('Provide a file name...')
    else:
        fig.savefig(os.path.join(outdir, outfile_name), format='pdf', pad_inches=0.1)    
    
    if ncols == 1:
        return fig

