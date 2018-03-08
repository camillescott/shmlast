import pytest
import pandas as pd

from shmlast.tests.utils import datadir
from shmlast.hits import BestHits
from shmlast.crbl import scale_evalues


def check_df_equals(dfA, dfB, col='E'):
    dfA = dfA.sort_values(col, inplace=False).sort_index(axis=1)
    dfA.reset_index(inplace=True, drop=True)

    dfB = dfB.sort_values(col, inplace=False).sort_index(axis=1)
    dfB.reset_index(inplace=True, drop=True)

    print('Checking DataFrame Equality:\n', dfA)
    print(dfB)

    return dfA.equals(dfB)


def test_besthits_inplace(datadir):
    '''Test BestHits.best_hits with inplace=True
    '''
    input_df = pd.read_csv(datadir('query.maf.csv'))
    expected_df = pd.read_csv(datadir('besthits.expected.csv'))

    bh = BestHits()
    bh.best_hits(input_df, inplace=True)
    
    assert check_df_equals(expected_df, input_df)


def test_besthits_non_inplace(datadir):
    '''Test BestHits.best_hits with inplace=False
    '''
    input_df = pd.read_csv(datadir('query.maf.csv'))
    expected_df = pd.read_csv(datadir('besthits.expected.csv'))

    bh = BestHits()
    results_df = bh.best_hits(input_df, inplace=False)
    
    assert check_df_equals(expected_df, results_df)


def test_reciprocal_best_hits(datadir):
    query_df = pd.read_csv(datadir('query.maf.csv'))
    db_df = pd.read_csv(datadir('db.maf.csv'))
    expected_df = pd.read_csv(datadir('reciprocals.expected.csv'))

    bh = BestHits()
    results_df = bh.reciprocal_best_hits(query_df, db_df, inplace=False)
    
    assert check_df_equals(results_df, expected_df)

def test_scale_evalues():
    test_df = pd.DataFrame({'E': [.1, 0.01, 0.0]})
    expected_df = pd.DataFrame({'E_scaled': [1.0, 2.0, 307.652655568]})
    result_df, col_name = scale_evalues(test_df)

    assert list(result_df[col_name]) == \
        pytest.approx(list(expected_df['E_scaled']))
