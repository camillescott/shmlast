from .hits import BestHits
from .last import lastdb_task, lastal_task, MafParser
from .transeq import transeq_task, rename_task
from .util import ShortenedPythonAction, title, Move


class ShmlastApp(TaskLoader):

    def __init__(self, directory, config=None):
        super(ShmlastApp, self).__init__()
        self.directory = directory
        try:
            mkdir(directory)
        except OSError:
            pass
        self.doit_config = {}
        if config is not None:
            self.doit_config.update(config)

    def tasks(self):
        raise NotImplementedError()

    def load_tasks(self, cmd, opt_values, pos_args):
        return list(self.tasks()), self.doit_config

    def run(self, doit_args=None, move=False):
        if doit_args is None:
            doit_args = ['run']
        runner = DoitMain(self)
        if move:
            with Move(self.directory):
                return runner.run(doit_args)
        else:
            return runner.run(doit_args)


class RBL(ShmlastApp):

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
        self.output_prefix = output_fn
        if self.output_fn is None:
            self.output_prefix = '{0}.x.{1}.rbl'.format(base(self.transcriptome_fn),
                                                        base(self.database_fn))
            self.output_fn = self.output_prefix + '.csv'
        self.unmapped_output_fn = self.output_prefix + '.unmapped.csv'
        

    def reciprocal_best_last_task(self):
       
        def do_reciprocals():
            rbh_df, qvd_df, dvq_df = RBL.get_reciprocals(self.translated_x_db_fn,
                                                     self.db_x_translated_fn,
                                                     self.bh)
            q_names = pd.read_csv(self.name_map_fn)
            d_names = pd.read_csv(self.database_name_map_fn)

            rbh_df.to_csv(self.unmapped_output_fn, index=False)
            qvd_df.to_csv(self.translated_x_db_fn + '.csv', index=False)
            dvq_df.to_csv(self.db_x_translated_fn + '.csv', index=False)

            rbh_df = self.backmap(rbh_df, q_names, d_names)
            qvd_df = self.backmap(qvd_df, q_names, d_names)
            dvq_df = self.backmap(dvq_df, q_names, d_names)

            rbh_df.to_csv(self.output_fn, index=False)
            qvd_df.to_csv(self.translated_x_db_fn + '.mapped.csv', index=False)
            dvq_df.to_csv(self.db_x_translated_fn + '.mapped.csv', index=False)

        td = {'name': 'reciprocal_best_last',
              'title': title,
              'actions': [ShortenedPythonAction(do_reciprocals)],
              'file_dep': [self.translated_x_db_fn,
                           self.db_x_translated_fn],
              'targets': [self.output_fn + '.csv',
                          self.translated_x_db_fn + '.csv',
                          self.db_x_translated_fn + '.csv',
                          self.output_fn + '.mapped.csv',
                          self.translated_x_db_fn + '.mapped.csv',
                          self.db_x_translated_fn + '.mapped.csv'],
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

