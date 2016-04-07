#/usr/bin/env python3

from doit.tools import run_once, create_folder, title_with_actions
from doit.task import clean_targets, dict_to_task
import sys

from .hits import BestHits
from .last import lastdb_task, lastal_task, MafParser
from .transeq import transeq_task, rename_task

class ReciprocalBestLAST(object):

    def __init__(self, transcriptome_fn, database_fn, output_fn,
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
        self.crbl_fn = '{0}.crbl.{1}.csv'.format(self.transcriptome_fn,
                                                 self.database_fn)

        self.bh = BestHits(comparison_cols=['E', 'q_aln_len'])
        self.output_fn = output_fn

    def reciprocal_best_last_task(self):
        
        def cmd():
            qvd_df = MafParser(self.translated_x_db_fn).read()
            qvd_df[['qg_name', 'q_frame']] = qvd_df.q_name.str.partition('_')[[0,2]]
            qvd_df.rename(columns={'q_name': 'translated_q_name',
                                   'qg_name': 'q_name'},
                          inplace=True)

            dvq_df = MafParser(self.db_x_translated_fn).read()
            dvq_df[['sg_name', 'frame']] = dvq_df.s_name.str.partition('_')[[0,2]]
            dvq_df.rename(columns={'s_name': 'translated_s_name',
                                   'sg_name': 's_name'},
                          inplace=True)
            
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
