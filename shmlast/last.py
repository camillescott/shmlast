#!/usr/bin/env python

from doit.task import clean_targets
from doit.tools import LongRunning
import glob
from itertools import count
import numpy as np
import os
import pandas as pd

from .profile import profile_task
from .util import create_doit_task as doit_task
from .util import which, parallel_fasta, title

LASTAL_CFG = { "params": [],
               "frameshift": 15 }
LASTDB_CFG = { "params": ["-w3"] }

def clean_lastdb(db_prefix):
    files = glob.glob('{0}.*'.format(db_prefix))
    for fn in files:
        try:
            os.remove(fn)
        except OSError as e:
            pass


@doit_task
@profile_task
def lastdb_task(db_fn, db_out_prefix=None, prot=True,
                params=LASTDB_CFG['params'], task_dep=None):
    '''Create a pydoit task to run lastdb.

    WARNING: This does not define a file_dep, to make sure it doesn't
    get executed when the dependency and targets already exist. This means
    that if the db_fn is created by another task, it MUST be defined before the
    lastdb task. This is a kludge which will eventually be fixed...

    Args:
        db_fn (str): The FASTA file to format.
        db_out_prefix (str): Prefix for the database files. Same as db_fn
                             if None (default).
        prot (bool): True if a protein FASTA, False otherwise.
        params (list): A list of additional parameters.
    Returns:
        dict: A pydoit task.
    '''

    cmd = [which('lastdb')]
    if prot:
        cmd.append('-p')
    if params is not None:
        cmd.extend([str(p) for p in params])
    if db_out_prefix is None:
        db_out_prefix = db_fn
    cmd.extend([db_out_prefix, db_fn])
    cmd = ' '.join(cmd)

    name = 'lastdb:' + os.path.basename(db_out_prefix)

    task_d =  {'name': name,
              'title': title,
              'actions': [cmd],
              'targets': ['{0}.prj'.format(db_out_prefix)],
              'uptodate': [True],
              'clean': [clean_targets,
                        (clean_lastdb, [db_out_prefix])]}
    if task_dep is not None:
        task_d['task_dep'] = task_dep

    return task_d


@doit_task
@profile_task
def lastal_task(query, db, out_fn, translate=False,
                frameshift=LASTAL_CFG['frameshift'], cutoff=0.00001, 
                n_threads=1, pbs=None, params=None):
    '''Create a pydoit task to run lastal

    Args:
        query (str): The file with the query sequences.
        db (str): The database file prefix.
        out_fn (str): Destination file for alignments.
        translate (bool): True if query is a nucleotide FASTA.
        frameshift (int): Frameshift penalty for translated alignment.
        n_threads (int): Number of threads to run with.
    Returns:
        dict: A pydoit task.
    '''

    lastal_exc = which('lastal')
    name = 'lastal:{0}'.format(os.path.join(out_fn))

    cmd = [lastal_exc]
    if translate:
        cmd.append('-F' + str(frameshift))
    if cutoff is not None:
        cutoff = round(1.0 / cutoff, 2)
        cmd.append('-D' + str(cutoff))
    if params is not None:
        cmd.extend(params)
    cmd.append(db)

    cmd =  parallel_fasta(query, out_fn, cmd, n_threads, pbs=pbs)

    return {'name': name,
            'title': title,
            'actions': [cmd],
            'targets': [out_fn],
            'file_dep': [query, db + '.prj'],
            'clean': [clean_targets]}


def next_or_raise(fp):
    counter = count()
    def func(raise_exc=True):
        line = fp.readline()
        n = next(counter)
        if raise_exc is True and line == '':
            raise RuntimeError('Malformed MAF file (line {0})'.format(n))
        return line
    return func


class MafParser(object):

    def __init__(self, filename, aln_strings=False, chunksize=10000, **kwargs):
        '''Parser for LAST's MAF output.

        Args:
            filename (str): The MAF file.
            aln_strings (bool): If True, parse the alignment strings as well.
            chunksize (int): Size of chunks to parse.
        '''
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
            guarded_next = next_or_raise(fp)
            n = 0
            while (True):
                line = guarded_next(raise_exc=False)
                if line == '':
                    break
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
                    line = guarded_next()
                    tokens = line.split()
                    cur_aln['s_name'] = tokens[1]
                    cur_aln['s_start'] = int(tokens[2])
                    cur_aln['s_aln_len'] = int(tokens[3])
                    cur_aln['s_strand'] = tokens[4]
                    cur_aln['s_len'] = int(tokens[5])
                    if self.aln_strings:
                        cur_aln['s_aln'] = tokens[6]

                    # First sequence info
                    line = guarded_next() 
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

