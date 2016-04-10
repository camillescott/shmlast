#/usr/bin/env python3

from doit.tools import run_once, create_folder, title_with_actions
from doit.task import clean_targets, dict_to_task
import sys

from .hits import BestHits
from .last import lastdb_task, lastal_task, MafParser
from .transeq import transeq_task, rename_task

class ReciprocalBestLAST(object):

    def __init__(self, transcriptome_fn, database_fn, output_fn=None,
                 cutoff=.00001, n_threads=1):

        self.transcriptome_fn = transcriptome_fn
        self.renamed_fn = self.transcriptome_fn + '.renamed'
        self.name_map_fn = self.transcriptome_fn + '.names.csv'
        self.translated_fn = self.renamed_fn + '.pep'
        self.database_fn = database_fn
        self.n_threads = n_threads
        self.cutoff = cutoff

        self.db_x_translated_fn = '{0}.x.{1}.maf'.format(self.database_fn,
                                                            self.translated_fn)
        self.translated_x_db_fn = '{0}.x.{1}.maf'.format(self.translated_fn,
                                                            self.database_fn)
        self.output_fn = output_fn
        if self.output_fn is None:
            self.output_fn = '{0}.rbl.{1}.csv'.format(self.transcriptome_fn,
                                                      self.database_fn)

        self.bh = BestHits(comparison_cols=['E', 'q_aln_len'])


    def reciprocal_best_last_task(self):
        
        def cmd():
            qvd_df = MafParser(self.translated_x_db_fn).read()
            qvd_df[['qg_name', 'q_frame']] = qvd_df.q_name.str.partition('_')[[0,2]]
            qvd_df.rename(columns={'q_name': 'translated_q_name',
                                   'qg_name': 'q_name'},
                          inplace=True)
            qvd_df['ID'] = qvd_df.index

            dvq_df = MafParser(self.db_x_translated_fn).read()
            dvq_df[['sg_name', 'frame']] = dvq_df.s_name.str.partition('_')[[0,2]]
            dvq_df.rename(columns={'s_name': 'translated_s_name',
                                   'sg_name': 's_name'},
                          inplace=True)
            dvq_df['ID'] = dvq_df.index
            
            self.bh.reciprocal_best_hits(qvd_df, dvq_df).to_csv(self.output_fn,
                                                                index=False)

        td = {'name': 'reciprocal_best_last',
              'title': title_with_actions,
              'actions': [cmd],
              'file_dep': [self.translated_x_db_fn,
                           self.db_x_translated_fn],
              'targets': [self.output_fn],
              'clean': [clean_targets]}
        
        return dict_to_task(td)

    def rename_task(self):
        return rename_task(self.transcriptome_fn,
                           self.renamed_fn,
                           name_map_fn=self.name_map_fn)

    def translate_task(self):
        return transeq_task(self.renamed_fn,
                            self.translated_fn)

    def format_transcriptome_task(self):
        return lastdb_task(self.translated_fn,
                           prot=True)

    def format_database_task(self):
        return lastdb_task(self.database_fn,
                           prot=True)

    def align_transcriptome_task(self):
        return lastal_task(self.translated_fn,
                           self.database_fn + '.lastdb',
                           self.translated_x_db_fn,
                           translate=False, 
                           cutoff=self.cutoff,
                           n_threads=self.n_threads)

    def align_database_task(self):
        return lastal_task(self.database_fn,
                           self.translated_fn + '.lastdb',
                           self.db_x_translated_fn,
                           translate=False, 
                           cutoff=self.cutoff,
                           n_threads=self.n_threads)

    def tasks(self):
        yield self.rename_task()
        yield self.translate_task()
        yield self.format_transcriptome_task()
        yield self.format_database_task()
        yield self.align_database_task()
        yield self.align_transcriptome_task()
        yield self.reciprocal_best_last_task()


class ConditionalReciprocalBestLAST(ReciprocalBestLAST):

    def __init__(self, transcriptome_fn, database_fn, output_fn=None,
                 cutoff=.00001, n_threads=1):

        self.crbl_output_fn = output_fn
        if output_fn is None:
            self.crbl_output_fn = '{0}.crbl.{1}.csv'.format(self.transcriptome_fn,
                                                  self.database_fn)

        super(ConditionalReciprocalBestLAST, self).__init__(transcriptome_fn,
                                                            database_fn,
                                                            output_fn=None,
                                                            cutoff=cutoff,
                                                            n_threads=n_threads)

    def scale_evalue(self, df, name='E'):
        df['E_s'] = df[name]
        df.loc[df['E_s'] == 0.0, 'E_s'] = 1e-300
        df['E_s'] = -np.log10(df['E_s'])

    def generate_model_task(self):

        def model_fit(rbh_df):
            cols = {'s_aln_len_A': 'length', 'E_A': 'E'}
            data = rbh_df[cols.keys()].rename(columns=cols)
            data.sort_values('length', inplace=True)
            scale_evalue(data)
            
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
                hits = df[(df['length'] >= fit_row.left) & (df['length'] < fit_row.right)]
                return hits['E_s'].mean()
            fit['fit'] = fit.apply(bin_mean, args=(data,), axis=1)
            
            return fit.dropna(), data

        def filter_from_model(model_df, hits_df, rbh_df, scale_evalue=False):
        
            if scale_evalue:
                scale_evalue(hits_df) 
            
            comp_df = pd.merge(hits_df[hits_df['ID'].isin(rbh_df['ID_A']) == False], model_df, 
                               left_on='s_aln_len', right_on='center')
        
            return comp_df[comp_df['E_s'] >= comp_df['fit']]

        def cmd():
            rbh_df = pd.read_csv(self.output_fn)
            hits_df = MafParser(self.translated_x_db_fn).read()

            model, data = model_fit(rbh_df)
            crbl_df = filter_from_model(model, hits_df, scale_evalue=True)

            # get subset of columns in rbh_df matching the crbl_df,
            # rename them and concat results
            rename = {}
            for col in crbl_df:
                if col + '_A' in rbh_df:
                    rename[col +'_A'] = col
                else:
                    rename[col] = col
            reciprocals = pd.concat([rbh_df.rename(columns=rename)[rename.values()], crbl_df],
                                    axis=1)
            # save results
            reciprocals.to_csv(self.crbl_output_fn, index=False)


