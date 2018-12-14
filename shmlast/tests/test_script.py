import pytest
from os import path

from shmlast.tests.utils import runscript, datadir, N_THREADS

@pytest.mark.parametrize('n_threads', N_THREADS,
                         ids=lambda n: 'n_threads={0}'.format(n))
@pytest.mark.benchmark(group='rbl-script')
def test_rbl(tmpdir, datadir, n_threads, benchmark):
    query = datadir('sacPom.cdna.fa')
    database = datadir('sacPom.pep.fa')
    
    args = ['rbl', '--n_threads', str(n_threads),
            '-q', query, '-d', database]
    status, out, err = benchmark.pedantic(runscript,
                                          args=('shmlast', args),
                                          kwargs={'directory': str(tmpdir)},
                                          iterations=1,
                                          rounds=1)

    assert status == 0
    assert tmpdir.ensure('sacPom.cdna.fa.x.sacPom.pep.fa.rbl.csv')


@pytest.mark.parametrize('n_threads', N_THREADS,
                         ids=lambda n: 'n_threads={0}'.format(n))
@pytest.mark.benchmark(group='crbl-script')
def test_crbl(tmpdir, datadir, n_threads, benchmark):
    query = datadir('sacPom.cdna.fa')
    database = datadir('sacPom.pep.fa')
    out_fn = str(tmpdir.join('test.csv'))
    
    args = ['crbl', '--n_threads', str(n_threads),
            '-q', query, '-d', database, '-o', out_fn]
    status, out, err = benchmark.pedantic(runscript,
                                          args=('shmlast', args),
                                          kwargs={'directory': str(tmpdir)},
                                          iterations=1,
                                          rounds=1)
    assert status == 0
    assert tmpdir.ensure(out_fn)
