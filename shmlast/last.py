#!/usr/bin/env python

from doit.task import clean_targets
from doit.tools import LongRunning
import glob
from itertools import count
import numpy as np
import os
import pandas as pd

from ope.io.maf import MafParser

from .profile import profile_task
from .util import create_doit_task as doit_task
from .util import which, title

float_info = np.finfo(float)

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
                n_threads=1, params=None):
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

    cmd = ['ope', 'parallel', '-j', n_threads, query,
           lastal_exc]
    if translate:
        cmd.append('-F' + str(frameshift))
    if cutoff is not None:
        cutoff = round(1.0 / cutoff, 2)
        cmd.append('-D' + str(cutoff))
    if params is not None:
        cmd.extend(params)
    cmd.extend([db, '>', out_fn])

    cmd = [str(token) for token in cmd]
    cmd = ' '.join(cmd)

    return {'name': name,
            'title': title,
            'actions': [cmd],
            'targets': [out_fn],
            'file_dep': [query, db + '.prj'],
            'clean': [clean_targets]}
