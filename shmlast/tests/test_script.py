
import pytest
from os import path

from shmlast.tests.utils import runscript, datadir

def test_rbl(tmpdir, datadir):
    query = datadir('sacPom.cdna.fa')
    database = datadir('sacPom.pep.fa')
    
    args = ['rbl', '-q', query, '-d', database]
    runscript('shmlast', args, directory=str(tmpdir))

    assert tmpdir.ensure('sacPom.cdna.fa.x.sacPom.pep.fa.rbl.csv')
