#/usr/bin/env python3

from doit.tools import run_once, create_folder
from doit.task import clean_targets, dict_to_task
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from ficus import FigureManager
import numpy as np
from os.path import basename as base
import pandas as pd
import seaborn as sns
import sys

from .hits import BestHits
from .last import lastdb_task, lastal_task, MafParser
from .transeq import transeq_task, rename_task
from .util import ShortenedPythonAction, title

class RBL(object):

    def __init__(self, transcriptome_fn, database_fn, output_fn=None,
                 cutoff=.00001, n_threads=1, n_nodes=None, use_existing_db=None):

        self.transcriptome_fn = transcriptome_fn
        self.renamed_fn = self.transcriptome_fn + '.renamed'
        self.name_map_fn = self.transcriptome_fn + '.names.csv'
        self.translated_fn = self.renamed_fn + '.pep'

        self.database_fn = database_fn
        self.renamed_database_fn = self.database_fn + '.renamed'
        self.database_name_map_fn = self.database_fn + '.names.csv'
        self.n_threads = n_threads
        self.cutoff = cutoff
        self.use_existing_db = use_existing_db

        self.db_x_translated_fn = '{0}.x.{1}.maf'.format(base(self.renamed_database_fn),
                                                         base(self.translated_fn))
        self.translated_x_db_fn = '{0}.x.{1}.maf'.format(base(self.translated_fn),
                                                         base(self.renamed_database_fn))
        self.output_fn = output_fn
        if self.output_fn is None:
            self.output_fn = '{0}.x.{1}.rbl.csv'.format(base(self.transcriptome_fn),
                                                        base(self.database_fn))

        self.bh = BestHits(comparison_cols=['E', 'EG2'])

    @staticmethod
    def get_reciprocals(query, database, besthits):
        qvd_df = MafParser(query).read()
        qvd_df[['qg_name', 'q_frame']] = qvd_df.q_name.str.partition('_')[[0,2]]
        qvd_df.rename(columns={'q_name': 'translated_q_name',
                               'qg_name': 'q_name'},
                      inplace=True)
        qvd_df['ID'] = qvd_df.index

        dvq_df = MafParser(database).read()
        dvq_df[['sg_name', 'frame']] = dvq_df.s_name.str.partition('_')[[0,2]]
        dvq_df.rename(columns={'s_name': 'translated_s_name',
                               'sg_name': 's_name'},
                      inplace=True)
        dvq_df['ID'] = dvq_df.index
        
        return besthits.reciprocal_best_hits(qvd_df, dvq_df), qvd_df, dvq_df

    def reciprocal_best_last_task(self):
       
        def do_reciprocals():
            rbh_df, q_df, d_df = RBL.get_reciprocals(self.translated_x_db_fn,
                                                     self.db_x_translated_fn,
                                                     self.bh)
            rbh_df.to_csv(self.output_fn, index=False)
            qvd_df.to_csv(self.translated_x_db_fn + '.csv', index=False)
            dvq_df.to_csv(self.db_x_translated_fn + '.csv', index=False)

        td = {'name': 'reciprocal_best_last',
              'title': title,
              'actions': [ShortenedPythonAction(do_reciprocals)],
              'file_dep': [self.translated_x_db_fn,
                           self.db_x_translated_fn],
              'targets': [self.output_fn,
                          self.translated_x_db_fn + '.csv',
                          self.db_x_translated_fn + '.csv'],
              'clean': [clean_targets]}
        
        return dict_to_task(td)

    def rename_transcriptome_task(self):
        return rename_task(self.transcriptome_fn,
                           self.renamed_fn,
                           name_map_fn=self.name_map_fn)

    def rename_database_task(self):
        return rename_task(self.database_fn,
                           self.renamed_database_fn,
                           prefix='db',
                           name_map_fn=self.database_name_map_fn)

    def translate_task(self):
        return transeq_task(self.renamed_fn,
                            self.translated_fn)

    def format_transcriptome_task(self):
        return lastdb_task(self.translated_fn,
                           prot=True)

    def format_database_task(self):
        return lastdb_task(self.renamed_database_fn,
                           prot=True,
                           use_existing=self.use_existing_db)

    def align_transcriptome_task(self):
        return lastal_task(self.translated_fn,
                           self.renamed_database_fn + '.lastdb',
                           self.translated_x_db_fn,
                           translate=False, 
                           cutoff=self.cutoff,
                           n_threads=self.n_threads)

    def align_database_task(self):
        return lastal_task(self.renamed_database_fn,
                           self.translated_fn + '.lastdb',
                           self.db_x_translated_fn,
                           translate=False, 
                           cutoff=self.cutoff,
                           n_threads=self.n_threads)

    def tasks(self):
        yield self.rename_transcriptome_task()
        yield self.rename_database_task()
        yield self.translate_task()
        yield self.format_transcriptome_task()
        yield self.format_database_task()
        yield self.align_database_task()
        yield self.align_transcriptome_task()
        yield self.reciprocal_best_last_task()


class CRBL(RBL):

    def __init__(self, transcriptome_fn, database_fn, output_fn=None,
                 model_fn=None, cutoff=.00001, n_threads=1, use_existing_db=True):

        self.crbl_output_fn = output_fn
        self.crbl_output_prefix = output_fn
        if output_fn is None:
            self.crbl_output_prefix = '{0}.x.{1}.crbl'.format(base(transcriptome_fn),
                                                            base(database_fn))
            self.crbl_output_fn = self.crbl_output_prefix + '.csv'

        self.model_fn = model_fn
        if model_fn is None:
            self.model_fn = self.crbl_output_prefix + '.model.csv'
        self.model_plot_fn = self.model_fn + '.plot.pdf'

        super(CRBL, self).__init__(transcriptome_fn,
                                    database_fn,
                                    output_fn=None,
                                    cutoff=cutoff,
                                    n_threads=n_threads,
                                    use_existing_db=use_existing_db)

    @staticmethod
    def scale_evalue(df, name='E'):
        df['E_s'] = df[name]
        df.loc[df['E_s'] == 0.0, 'E_s'] = 1e-300
        df['E_s'] = -np.log10(df['E_s'])

    @staticmethod
    def fit_crbh_model(rbh_df):

        data = rbh_df[['s_aln_len', 'E']].rename(columns={'s_aln_len':'length'})
        data.sort_values('length', inplace=True)
        CRBL.scale_evalue(data)

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
            return hits['E_s'].mean()
        fit['fit'] = fit.apply(bin_mean, args=(data,), axis=1)
        model_df = fit.dropna()

        return model_df

    @staticmethod
    def filter_from_model(model_df, rbh_df, hits_df):

        CRBL.scale_evalue(hits_df)
        comp_df = pd.merge(hits_df[hits_df['ID'].isin(rbh_df['ID']) == False], 
                           model_df, left_on='s_aln_len', right_on='center')
    
        crbl_df = comp_df[comp_df['E_s'] >= comp_df['fit']]

        del crbl_df['center']
        del crbl_df['left']
        del crbl_df['right']
        del crbl_df['fit']
        del crbl_df['size']
        del crbl_df['translated_q_name']

        return crbl_df

    @staticmethod
    def plot_crbl_fit(model_df, rbh_df, hits_df, model_plot_fn, show=False,
                     figsize=(10,10)):

        plt.style.use('seaborn-ticks')

        with FigureManager(model_plot_fn, show=show, 
                           figsize=figsize) as (fig, ax):

            scatter_kws = {'s': 10, 'alpha':0.7}
            scatter_kws['c'] = sns.xkcd_rgb['ruby']
            scatter_kws['marker'] = 'o'
            line_kws = {'c': sns.xkcd_rgb['red wine'], 
                        'label':'Query Hits Regression'}
            sample_size = min(len(hits_df), 10000)
            sns.regplot('s_aln_len', 'E_s', hits_df.sample(sample_size), order=1, 
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
            ax.set_ylim(0, max(model_df['fit'].max(), hits_df['E'].max()) + 50)
            ax.set_title('CRBL Fit')

    def crbl_fit_task(self):

        def do_crbl_fit():
            rbh_df = pd.read_csv(self.output_fn)
            model_df = self.fit_crbh_model(rbh_df)
            model_df.to_csv(self.model_fn, index=False)

        td = {'name': 'fit_crbl_model',
              'title': title,
              'actions': [ShortenedPythonAction(do_crbl_fit)],
              'file_dep': [self.output_fn, 
                           self.translated_x_db_fn + '.csv',
                           self.db_x_translated_fn + '.csv'],
              'targets': [self.model_fn]}

        return dict_to_task(td)

    def crbl_filter_task(self):

        def do_crbl_filter():
            model_df = pd.read_csv(self.model_fn)
            rbh_df = pd.read_csv(self.output_fn)
            hits_df = pd.read_csv(self.translated_x_db_fn + '.csv')

            filtered_df = self.filter_from_model(model_df, rbh_df, hits_df)
            results = pd.concat([rbh_df, filtered_df], axis=0)
            results.to_csv(self.crbl_output_fn, index=False)

            self.plot_crbl_fit(model_df, rbh_df, hits_df, self.model_plot_fn)

        td = {'name': 'filter_crbl_hits',
              'title': title,
              'actions': [ShortenedPythonAction(do_crbl_filter)],
              'file_dep': [self.output_fn, 
                           self.translated_x_db_fn + '.csv',
                           self.db_x_translated_fn + '.csv',
                           self.model_fn],
              'targets': [self.crbl_output_fn, 
                          self.model_plot_fn]}

        return dict_to_task(td)

    def tasks(self):
        for tsk in super(CRBL, self).tasks():
            yield tsk
        yield self.crbl_fit_task()
        yield self.crbl_filter_task()
        
