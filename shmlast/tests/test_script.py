
import pytest
from os import path

from shmlast.tests.utils import runscript, datadir

def test_rbl(tmpdir, datadir):
    query = datadir('sacPom.cdna.fa')
    database = datadir('sacPom.pep.fa')
    
    print(query, database, tmpdir)
    args = ['rbl', '-q', query, '-d', database]
    runscript('shmlast', args, directory=str(tmpdir))

    assert tmpdir.ensure('sacPom.cdna.fa.x.sacPom.pep.fa.rbl.csv')

def test_crbl(tmpdir, datadir):
    query = datadir('sacPom.cdna.fa')
    database = datadir('sacPom.pep.fa')
    out_fn = str(tmpdir.join('test.csv'))
    
    args = ['crbl', '-q', query, '-d', database, '-o', out_fn]
    status, out, err = runscript('shmlast', args, directory=str(tmpdir))
    assert status == 0
    assert tmpdir.ensure(out_fn)
