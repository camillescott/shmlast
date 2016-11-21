from os import path, getcwd, mkdir

from doit.tools import run_once, create_folder
from doit.task import clean_targets, dict_to_task
from doit.cmd_base import TaskLoader
from doit.doit_cmd import DoitMain
import pandas as pd

from .crbl import (get_reciprocal_best_last_translated, backmap_names,
                   scale_evalues, fit_crbh_model, filter_hits_from_model,
                   plot_crbh_fit)

from .profile import StartProfiler, profile_task
from .last import lastdb_task, lastal_task, MafParser
from .translate import translate_task, rename_task
from .util import ShortenedPythonAction, title, hidden_fn
from .util import create_doit_task as doit_task


class ShmlastApp(TaskLoader):

    def __init__(self, directory=None, config=None):
        super(ShmlastApp, self).__init__()

        if directory is None:
            directory = getcwd()
        self.directory = directory
        try:
            mkdir(self.directory)
        except OSError:
            pass
        self.doit_config = {'verbosity': 2}
        if config is not None:
            self.doit_config.update(config)

    def tasks(self):
        raise NotImplementedError()

    def load_tasks(self, cmd, opt_values, pos_args):
        return list(self.tasks()), self.doit_config

    def run(self, doit_args=None, move=False, profile_fn=None):
        if doit_args is None:
            doit_args = ['run']
        runner = DoitMain(self)

        print('\n--- Begin Task Execution ---')
        if profile_fn is not False and doit_args[0] == 'run':
            with StartProfiler(filename=profile_fn):
                return runner.run(doit_args)
        else:
            return runner.run(doit_args)


class RBL(ShmlastApp):

    def __init__(self, query_fn, database_fn, output_fn=None,
                 cutoff=.00001, n_threads=1, pbs=None, directory=None):
        '''Generate and manage the pydoit tasks for the RBL pipeline.

        Args:
            query_fn (str): The query filename.
            database_fn (str): The database filename.
            output_fn (str): Filename to store results in.
            cutoff (float): The score cutoff.
            n_threads (int): Number of threads to run on.
            pbs (bool): Whether to generate parallel command for PBS
                environment.
            directory (str): The directory to run tasks in.
        '''

        self.query_fn = query_fn
        self.renamed_query_fn = hidden_fn(path.basename(self.query_fn))
        self.query_name_map_fn = self.renamed_query_fn + '.names.csv'
        self.translated_query_fn = self.renamed_query_fn + '.pep'

        self.database_fn = database_fn
        self.renamed_database_fn = hidden_fn(path.basename(self.database_fn))
        self.database_name_map_fn = self.renamed_database_fn + '.names.csv'

        self.n_threads = n_threads
        self.pbs = pbs
        self.cutoff = cutoff

        self.db_x_query_fn = '{0}.x.{1}.maf'.format(self.renamed_database_fn,
                                                    self.translated_query_fn.strip('.'))

        self.query_x_db_fn = '{0}.x.{1}.maf'.format(self.translated_query_fn,
                                                    self.renamed_database_fn.strip('.'))
        
        self.output_fn = output_fn
        if self.output_fn is None:
            self.output_fn = '{q}.x.{d}.rbl.csv'.format(q=path.basename(self.query_fn),
                                                        d=path.basename(self.database_fn))
        self.unmapped_output_fn = hidden_fn(self.output_fn)
        
        dep_file = '.{0}.shmlast.doit'.format(path.basename(self.query_fn))
        super(RBL, self).__init__(directory=directory, 
                                  config={'dep_file': dep_file})

    @doit_task
    @profile_task
    def reciprocal_best_last_task(self):
       
        def do_reciprocals():
            rbh_df, qvd_df, dvq_df = get_reciprocal_best_last_translated(self.query_x_db_fn,
                                                                         self.db_x_query_fn)
            q_names = pd.read_csv(self.query_name_map_fn)
            d_names = pd.read_csv(self.database_name_map_fn)

            rbh_df.to_csv(self.unmapped_output_fn, index=False)
            rbh_df = backmap_names(rbh_df, q_names, d_names)
            rbh_df.to_csv(self.output_fn, index=False)

        td = {'name': 'reciprocal_best_last',
              'title': title,
              'actions': [ShortenedPythonAction(do_reciprocals)],
              'file_dep': [self.query_x_db_fn,
                           self.db_x_query_fn,
                           self.query_name_map_fn,
                           self.database_name_map_fn],
              'targets': [self.unmapped_output_fn,
                          self.output_fn],
              'clean': [clean_targets]}
        
        return td

    def rename_transcriptome_task(self):
        return rename_task(self.query_fn,
                           self.renamed_query_fn,
                           name_map_fn=self.query_name_map_fn)

    def rename_database_task(self):
        return rename_task(self.database_fn,
                           self.renamed_database_fn,
                           prefix='db',
                           name_map_fn=self.database_name_map_fn)

    def translate_task(self):
        return translate_task(self.renamed_query_fn,
                              self.translated_query_fn)

    def format_transcriptome_task(self):
        return lastdb_task(self.translated_query_fn,
                           prot=True)

    def format_database_task(self):
        return lastdb_task(self.renamed_database_fn,
                           prot=True)

    def align_transcriptome_task(self):
        return lastal_task(self.translated_query_fn,
                           self.renamed_database_fn,
                           self.query_x_db_fn,
                           translate=False, 
                           cutoff=self.cutoff,
                           n_threads=self.n_threads,
                           pbs=self.pbs)

    def align_database_task(self):
        return lastal_task(self.renamed_database_fn,
                           self.translated_query_fn,
                           self.db_x_query_fn,
                           translate=False, 
                           cutoff=self.cutoff,
                           n_threads=self.n_threads,
                           pbs=self.pbs)


    def tasks(self):
        '''Iterator over all tasks in pipeline.
        '''
        yield self.rename_transcriptome_task()
        yield self.rename_database_task()
        yield self.translate_task()
        yield self.format_transcriptome_task()
        yield self.format_database_task()
        yield self.align_database_task()
        yield self.align_transcriptome_task()
        yield self.reciprocal_best_last_task()


class CRBL(RBL):

    def __init__(self, query_fn, database_fn, output_fn=None,
                 model_fn=None, cutoff=.00001, n_threads=1, pbs=None):
        '''Generate and manage the pydoit tasks for the CRBL pipeline.

        Args:
            query_fn (str): The query filename.
            database_fn (str): The database filename.
            output_fn (str): Filename to store results in.
            model_fn (str): Filename to store the model in.
            cutoff (float): The score cutoff.
            n_threads (int): Number of threads to run on.
            pbs (bool): Whether to generate parallel command for PBS
                environment.
            directory (str): The directory to run tasks in.
        '''
        prefix = '{q}.x.{d}.crbl'.format(q=path.basename(query_fn),
                                         d=path.basename(database_fn))

        self.crbl_output_fn = output_fn
        if output_fn is None:
            self.crbl_output_fn = prefix + '.csv'
        self.unmapped_crbl_output_fn = hidden_fn(self.crbl_output_fn)

        self.model_fn = model_fn
        if model_fn is None:
            self.model_fn = prefix + '.model.csv'
            self.model_plot_fn = prefix + '.model.plot.pdf'
        else:
            self.model_plot_fn = self.model_fn + '.plot.pdf'

        super(CRBL, self).__init__(query_fn,
                                    database_fn,
                                    output_fn=None,
                                    cutoff=cutoff,
                                    n_threads=n_threads,
                                    pbs=pbs)

    @doit_task
    @profile_task
    def crbl_fit_and_filter_task(self):

        def do_crbl_fit_and_filter():
            rbh_df, hits_df, _ = get_reciprocal_best_last_translated(self.query_x_db_fn,
                                                                     self.db_x_query_fn)
            
            model_df = fit_crbh_model(rbh_df)
            model_df.to_csv(self.model_fn, index=False)
            model_df = pd.read_csv(self.model_fn)

            filtered_df = filter_hits_from_model(model_df, rbh_df, hits_df)
            results = pd.concat([rbh_df, filtered_df], axis=0)
            del results['translated_q_name']

            q_names = pd.read_csv(self.query_name_map_fn)
            d_names = pd.read_csv(self.database_name_map_fn)
            
            results = backmap_names(results, q_names, d_names)
            results.to_csv(self.crbl_output_fn, index=False)

            plot_crbh_fit(model_df, hits_df, self.model_plot_fn)

        td = {'name': 'fit_and_filter_crbl_hits',
              'title': title,
              'actions': [ShortenedPythonAction(do_crbl_fit_and_filter)],
              'file_dep': [self.query_x_db_fn,
                           self.query_name_map_fn,
                           self.database_name_map_fn],
              'targets': [self.crbl_output_fn, 
                          self.model_plot_fn,
                          self.model_fn],
              'clean': [clean_targets]}
        
        return td

    def tasks(self):
        '''Iterator over all pipeline tasks.
        '''
        for tsk in super(CRBL, self).tasks():
            if tsk.name != 'reciprocal_best_last':
                yield tsk
        yield self.crbl_fit_and_filter_task()
        
