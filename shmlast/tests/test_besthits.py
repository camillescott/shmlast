import pytest
import pandas as pd

from shmlast.tests.utils import datadir
from shmlast.hits import BestHits

def test_besthits_inplace(datadir):
    '''Test BestHits.best_hits with inplace=True
    '''
    input_df = pd.read_csv(datadir('besthits.input.csv'))
    expected_df = pd.read_csv(datadir('besthits.expected.csv'))
    expected_df.sort_values('E', inplace=True)
    expected_df.reset_index(inplace=True, drop=True)

    bh = BestHits()
    bh.best_hits(input_df)
    input_df.sort_values('E', inplace=True)
    input_df.reset_index(inplace=True, drop=True)
    
    assert expected_df.equals(input_df)

def test_besthits_non_inplace(datadir):
    '''Test BestHits.best_hits with inplace=False
    '''
    input_df = pd.read_csv(datadir('besthits.input.csv'))
    expected_df = pd.read_csv(datadir('besthits.expected.csv'))
    expected_df.sort_values('E', inplace=True)
    expected_df.reset_index(inplace=True, drop=True)

    bh = BestHits()
    results = bh.best_hits(input_df, inplace=False)
    results.sort_values('E', inplace=True)
    results.reset_index(inplace=True, drop=True)
    
    assert expected_df.equals(results)


