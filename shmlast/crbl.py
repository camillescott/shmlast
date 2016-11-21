#/usr/bin/env python3

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from ficus import FigureManager
import numpy as np
from os import path
import pandas as pd
import seaborn as sns

from .hits import BestHits
from .last import MafParser

def get_reciprocal_best_last_translated(query_maf, database_maf):
    '''Perform Reciprocal Best Hits between the given MAF files.

    Args:
        query_maf (str): The query MAF file.
        database_maf (str): The translated datbase MAF file.
    Returns:
        tuple: DataFrames with the RBH's, query vs database, and database vs
            query hits.
    '''
    bh = BestHits(comparison_cols=['E', 'EG2'])
    qvd_df = MafParser(query_maf).read()
    qvd_df[['qg_name', 'q_frame']] = qvd_df.q_name.str.partition('_')[[0,2]]
    qvd_df.rename(columns={'q_name': 'translated_q_name',
                           'qg_name': 'q_name'},
                  inplace=True)
    qvd_df['ID'] = qvd_df.index

    dvq_df = MafParser(database_maf).read()
    dvq_df[['sg_name', 'frame']] = dvq_df.s_name.str.partition('_')[[0,2]]
    dvq_df.rename(columns={'s_name': 'translated_s_name',
                           'sg_name': 's_name'},
                  inplace=True)
    dvq_df['ID'] = dvq_df.index
    
    return bh.reciprocal_best_hits(qvd_df, dvq_df), qvd_df, dvq_df


def backmap_names(results_df, q_names, d_names):
    '''Map names from translated RBH's to original query and database names.

    Args:
        results_df (pandas.DataFrame): The results to backmap.
        q_names (pandas.DataFrame): Query name map.
        d_names (pandas.DataFrame): Database name map.
    Returns:
        pandas.DataFrame: Reference to results_df.
    '''

    results_df = pd.merge(results_df, 
                          q_names, 
                          left_on='q_name',
                          right_on='new_name')
    results_df['q_name'] = results_df['old_name']
    del results_df['old_name']

    results_df = pd.merge(results_df, 
                          d_names, 
                          left_on='s_name',
                          right_on='new_name')
    results_df['s_name'] = results_df['old_name']
    del results_df['old_name']
    del results_df['new_name_x']
    del results_df['new_name_y']

    return results_df


def scale_evalues(df, name='E', inplace=False):
    '''Log scale the evalue column specified by name.

    Args:
        df (pandas.DataFrame): The data.
        name (str): Column name with the evalues.
        inplace (bool): Perform the scaling inplace.
    Returns:
        tuple: The scaled DataFrame and the new column name of the scaled
            values.
    '''

    scaled_col_name = name + '_scaled'
    if inplace is False:
        df = df.copy()
    df[scaled_col_name] = df[name]
    df.loc[df[scaled_col_name] == 0.0, scaled_col_name] = 1e-300
    df[scaled_col_name] = -np.log10(df[scaled_col_name])
    return df, scaled_col_name


def fit_crbh_model(rbh_df, length_col='s_aln_len', feature_col='E'):
    '''Build the CRBH model on the given RBH's.

    Args:
        rbh_df (pandas.DataFrame): DataFrame with RBH's.
        length_col (str): The column with the subject lengths.
        feature_col (str): Score column to train on.
    Returns:
        pandas.DataFrame: The model.
    '''

    data = rbh_df[[length_col, feature_col]].rename(columns={length_col:'length'})
    data.sort_values('length', inplace=True)
    _, feature_col = scale_evalues(data, name=feature_col, inplace=True)

    # create a DataFrame for the model, staring with the alignment lengths
    fit = pd.DataFrame(np.arange(10, data['length'].max()), 
                       columns=['center'], dtype=int)
    
    # create the bins
    fit['size'] = fit['center'] * 0.1
    fit.loc[fit['size'] < 5, 'size'] = 5
    fit['size'] = fit['size'].astype(int)
    fit['left'] = fit['center'] - fit['size']
    fit['right'] = fit['center'] + fit['size']
    
    # do the fitting: it's just a sliding window with an increasing size
    def bin_mean(fit_row, df):
        hits = df[(df['length'] >= fit_row.left) & (df['length'] <= fit_row.right)]
        return hits[feature_col].mean()
    fit['fit'] = fit.apply(bin_mean, args=(data,), axis=1)
    model_df = fit.dropna()

    return model_df


def filter_hits_from_model(model_df, rbh_df, hits_df, feature_col='E',
                           id_col='ID', length_col='s_aln_len'):
    '''Filter a DataFrame of LAST best hits using the CRBH model.

    Args:
        model_df (pandas.DataFrame): The CRBH model.
        rbh_df (pandas.DataFrame): The RBH's.
        hits_df (pandas.DataFrame): The query vs database hits.
        feature_col (str): Column name of scores.
        id_col (str): Column with unique ID of hits.
        length_col (str): Column name to use for length.
    Returns:
        pandas.DataFrame: The CRBH's.
    '''

    hits_df, _ = scale_evalues(hits_df, name=feature_col, inplace=False)
    rbh_df, scaled_feature_col = scale_evalues(rbh_df, name=feature_col, inplace=False)

    # Merge the model into the subset of the hits which aren't in RBH
    comp_df = pd.merge(hits_df[hits_df[id_col].isin(rbh_df[id_col]) == False], 
                       model_df, left_on=length_col, right_on='center')

    crbl_df = comp_df[comp_df[scaled_feature_col] >= comp_df['fit']]

    del crbl_df['center']
    del crbl_df['left']
    del crbl_df['right']
    del crbl_df['fit']
    del crbl_df['size']

    return crbl_df


def plot_crbh_fit(model_df, hits_df, model_plot_fn, show=False,
                  figsize=(10,10), feature_col='E', length_col='s_aln_len'):

    plt.style.use('seaborn-ticks')

    with FigureManager(model_plot_fn, show=show, 
                       figsize=figsize) as (fig, ax):

        scatter_kws = {'s': 10, 'alpha':0.7}
        scatter_kws['c'] = sns.xkcd_rgb['ruby']
        scatter_kws['marker'] = 'o'
        line_kws = {'c': sns.xkcd_rgb['red wine'], 
                    'label':'Query Hits Regression'}
        sample_size = min(len(hits_df), 5000)
        hits_df, scaled_col = scale_evalues(hits_df, name=feature_col,
                                            inplace=False)
        sns.regplot(length_col, scaled_col, hits_df.sample(sample_size), order=1, 
                    label='Query Hits', scatter_kws=scatter_kws, 
                    line_kws=line_kws, color=scatter_kws['c'], ax=ax)

        scatter_kws['c'] = sns.xkcd_rgb['twilight blue']
        scatter_kws['marker'] = 's'
        sns.regplot('center', 'fit', model_df, 
                    fit_reg=False, x_jitter=True, y_jitter=True, ax=ax,
                    label='CRBL Fit', scatter_kws=scatter_kws, line_kws=line_kws)

        leg = ax.legend(fontsize='medium', scatterpoints=3, frameon=True)
        leg.get_frame().set_linewidth(1.0)

        ax.set_xlim(model_df['center'].min(), model_df['center'].max())
        ax.set_ylim(0, max(model_df['fit'].max(), hits_df[scaled_col].max()) + 50)
        ax.set_ylabel('Score ({0})'.format(feature_col))
        ax.set_xlabel('Alignment Length')


