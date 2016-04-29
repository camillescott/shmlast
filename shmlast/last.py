#!/usr/bin/env python

from doit.task import clean_targets
from doit.tools import LongRunning
import glob
import numpy as np
import os
import pandas as pd

from .util import create_doit_task as doit_task
from .util import which, parallel_fasta, title

LASTAL_CFG = { "params": "",
               "frameshift": 15 }
LASTDB_CFG = { "params": "-w3" }

def clean_lastdb(db_prefix):
    files = glob.glob('{0}.*'.format(db_prefix))
    for fn in files:
        try:
            os.remove(fn)
        except OSError as e:
            pass


@doit_task
def lastdb_task(db_fn, db_out_suffix=None, prot=True, cfg=LASTDB_CFG,
                use_existing=None):
    '''Create a pydoit task to run lastdb.

    Args:
        db_fn (str): The FASTA file to format.
        db_out_prefix (str): Prefix for the database files. Defaults
            to `<db_fn>.lastdb`.
        cfg (dict): Config for the command. Shoud contain an entry
            named "params" storing a str.
        prot (bool): True if a protein FASTA, False otherwise.
    Returns:
        dict: A pydoit task.
    '''

    exc = which('lastdb')
    params = cfg['params']
    if use_existing is not None:
        db_out_suffix = use_existing
    if db_out_suffix is None:
        db_out_suffix = '.lastdb'
    db_out_prefix = db_fn + db_out_suffix
    
    cmd = [exc]
    if prot:
        cmd.append('-p')
    cmd.extend([db_out_prefix, db_fn])
    cmd = ' '.join(cmd)

    name = 'lastdb:{0}'.format(os.path.basename(db_out_prefix))

    tskd = {'name': name,
            'title': title,
            'actions': [cmd],
            'targets': [db_out_prefix + '.prj'],
            'uptodate': [True],
            'clean': [clean_targets,
                      (clean_lastdb, [db_out_prefix])]}
    if not use_existing:
        tskd['file_dep'] = [db_fn]

    return tskd

@doit_task
def lastal_task(query, db, out_fn, cutoff=0.00001, n_threads=1,
                    translate=False, cfg=LASTAL_CFG, pbs=False):
    '''Create a pydoit task to run lastal

    Args:
        query (str): The file with the query sequences.
        db (str): The database file prefix.
        out_fn (str): Destination file for alignments.
        translate (bool): True if query is a nucleotide FASTA.
        n_threads (int): Number of threads to run with.
        cfg (dict): Config, must contain key params holding str.
    Returns:
        dict: A pydoit task.
    '''

    lastal_exc = which('lastal')
    parallel_exc = which('parallel-fasta')

    params = cfg['params']
    lastal_cmd = [lastal_exc]
    if translate:
        lastal_cmd.append('-F' + str(cfg['frameshift']))
    if cutoff is not None:
        cutoff = round(1.0 / cutoff)
        lastal_cmd.append('-D' + str(cutoff))
    lastal_cmd.append(db)
    lastal_cmd = ' '.join(lastal_cmd)

    cmd = parallel_fasta(query, out_fn, lastal_cmd, n_threads, pbs=pbs)

    name = 'lastal:{0}'.format(os.path.join(out_fn))

    return {'name': name,
            'title': title,
            'actions': [LongRunning(cmd)], 
            'targets': [out_fn],
            'file_dep': [query, db + '.prj'],
            'clean': [clean_targets]}


class MafParser(object):

    def __init__(self, filename, aln_strings=False, chunksize=10000, **kwargs):
        self.chunksize = chunksize
        self.filename = filename
        self.aln_strings = aln_strings
        self.LAMBDA = None
        self.K = None

    def read(self):
        '''Read the entire file at once and return as a single DataFrame.
        '''
        return pd.concat(iter(self), ignore_index=True)

    def __iter__(self):
        '''Iterator yielding DataFrames of length chunksize holding MAF alignments.

        An extra column is added for bitscore, using the equation described here:
            http://last.cbrc.jp/doc/last-evalues.html

        Args:
            fn (str): Path to the MAF alignment file.
            chunksize (int): Alignments to parse per iteration.
        Yields:
            DataFrame: Pandas DataFrame with the alignments.
        '''
        data = []
        with open(self.filename) as fp:
            while (True):
                try:
                    line = next(fp)
                except StopIteration:
                    break
                else:
                    line = line.strip()
                if not line:
                    continue
                if line.startswith('#'):
                    if 'lambda' in line:
                        meta = line.strip(' #').split()
                        meta = {k:v for k, _, v in map(lambda x: x.partition('='), meta)}
                        self.LAMBDA = float(meta['lambda'])
                        self.K = float(meta['K'])
                    else:
                        continue
                if line.startswith('a'):
                    cur_aln = {}

                    # Alignment info
                    tokens = line.split()
                    for token in tokens[1:]:
                        key, _, val = token.strip().partition('=')
                        cur_aln[key] = float(val)

                    # First sequence info
                    line = next(fp)
                    line = line.strip()
                    tokens = line.split()
                    cur_aln['s_name'] = tokens[1]
                    cur_aln['s_start'] = int(tokens[2])
                    cur_aln['s_aln_len'] = int(tokens[3])
                    cur_aln['s_strand'] = tokens[4]
                    cur_aln['s_len'] = int(tokens[5])
                    if self.aln_strings:
                        cur_aln['s_aln'] = tokens[6]

                    # First sequence info
                    line = next(fp)
                    line = line.strip()
                    tokens = line.split()
                    cur_aln['q_name'] = tokens[1]
                    cur_aln['q_start'] = int(tokens[2])
                    cur_aln['q_aln_len'] = int(tokens[3])
                    cur_aln['q_strand'] = tokens[4]
                    cur_aln['q_len'] = int(tokens[5])
                    if self.aln_strings:
                        cur_aln['q_aln'] = tokens[6]

                    data.append(cur_aln)
                    if len(data) >= self.chunksize:
                        if self.LAMBDA is None:
                            raise RuntimeError("old version of lastal; please update")
                        yield self._build_df(data)
                        data = []

        if data:
            yield self._build_df(data)

    def _build_df(self, data):

        def _fix_sname(name):
            new, _, _ = name.partition(',')
            return new

        df = pd.DataFrame(data)
        df['s_name'] = df['s_name'].apply(_fix_sname)
        setattr(df, 'LAMBDA', self.LAMBDA)
        setattr(df, 'K', self.K)
        df['bitscore'] = (self.LAMBDA * df['score'] - np.log(self.K)) / np.log(2)

        return df

