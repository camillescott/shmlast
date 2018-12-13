# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.
#
# Note: this code came from dammit (https://github.com/dib-lab/dammit)
# because Camille was too lazy to spin off a third package for the
# parsers only.
#

from itertools import count
from sys import stderr

import pandas as pd
import numpy as np


class EmptyFile(Exception):
    pass


def warn_empty(msg):
    '''Warn that a file is empty.'''
    print('\nWARNING: Empty file: {0}\n'.format(msg), file=stderr)


def next_or_raise(fp):
    '''Get the next line and raise an exception if its empty.
    '''
    counter = count()
    def func(raise_exc=True):
        line = fp.readline()
        n = next(counter)
        if raise_exc is True and line == '':
            raise RuntimeError('Malformed file (line {0})'.format(n))
        return line
    return func


def convert_dtypes(df, dtypes):
    '''Convert the columns of a DataFrame to the types specified
    in the given dictionary, inplace.

    Args:
        df (DataFrame): The DataFrame to convert.
        dtypes (dict): Dictionary mapping columns to types.
    '''

    for c in df.columns:
        try:
            df[c] = df[c].astype(dtypes[c])
        except KeyError:
            pass


class BaseParser(object):

    def __init__(self, filename):
        self.filename = filename

    def raise_empty(self):
        raise EmptyFile('Empty file: {0}'.format(self.filename))


class ChunkParser(BaseParser):

    def __init__(self, filename, chunksize=10000):
        '''
        Args:
            filename (str): Path to the file to parse.
            chunksize (int): Number of entries to yield per call.
        '''

        self.chunksize = chunksize
        super(ChunkParser, self).__init__(filename)

    def __iter__(self):
        raise NotImplementedError()
        yield

    def read(self):
        '''Read the entire file at once and return as a single DataFrame.
        '''
        try:
            return pd.concat(self, ignore_index=True)
        except (EmptyFile, ValueError) as e:
            # no objects, return an empty dataframe
            return self.empty()

    def empty(self):
        '''Get an empty DataFrame with the appropriate columns.
        '''

        df = pd.DataFrame(columns=[k for k, _ in self.columns])
        convert_dtypes(df, dict(self.columns))
        return df


class MafParser(ChunkParser):

    columns = [('E', float),
               ('EG2', float),
               ('q_aln_len', int),
               ('q_len', int),
               ('q_name', str),
               ('q_start', int),
               ('q_strand', str),
               ('s_aln_len', int),
               ('s_len', int),
               ('s_name', str),
               ('s_start', int),
               ('s_strand', str),
               ('score', float),
               ('bitscore', float)]

    def __init__(self, filename, aln_strings=False, chunksize=10000, **kwargs):
        self.aln_strings = aln_strings
        self.LAMBDA = None
        self.K = None
        super(MafParser, self).__init__(filename, chunksize=chunksize, **kwargs)

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
        n_entries = 0
        with open(self.filename) as fp:
            guarded_next = next_or_raise(fp)
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
                    n_entries += 1
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

        if n_entries == 0:
            self.raise_empty()
        if data:
            yield self._build_df(data)

    def _build_df(self, data):
        if not data:
            self.raise_empty()

        def _fix_sname(name):
            new, _, _ = name.partition(',')
            return new

        df = pd.DataFrame(data)
        df['s_name'] = df['s_name'].apply(_fix_sname)
        setattr(df, 'LAMBDA', self.LAMBDA)
        setattr(df, 'K', self.K)
        df['bitscore'] = (self.LAMBDA * df['score'] - np.log(self.K)) / np.log(2)

        return df
