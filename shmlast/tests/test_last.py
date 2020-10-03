#!/usr/bin/env python
from __future__ import print_function

from ope.io.maf import MafParser
import pytest

import json
import glob
import os
import sys

from shmlast.tests.utils import datadir, run_task, run_tasks, check_status, touch, N_THREADS
from shmlast.last import lastal_task
from shmlast.last import lastdb_task

LASTDB_EXTENSIONS = ['.bck', '.des', '.prj', '.sds', '.ssp', '.suf', '.tis']


@pytest.fixture(scope='module')
def lastdb_dir(tmpdir_factory, datadir):
    d = tmpdir_factory.mktemp('sacpom_lastdb')
    with d.as_cwd():
        data = datadir('sacPom.pep.fa')
        task = lastdb_task(data, data, prot=True)
        result = run_tasks([task], ['run'])
        assert result == 0
    
        return d


def test_lastdb_task_nucl(tmpdir, datadir):
    with tmpdir.as_cwd():
        tf = datadir('test-transcript.fa')

        task = lastdb_task(tf, tf, prot=False)
        run_tasks([task], ['run'])
        status = check_status(task)
        print('PATH:', os.environ['PATH'], file=sys.stderr)

        for ext in LASTDB_EXTENSIONS:
            assert os.path.isfile(tf + ext)

        assert status.status == 'up-to-date'


def test_lastdb_task_prot(tmpdir, datadir):
    with tmpdir.as_cwd():
        tf = datadir('test-protein.fa')

        task = lastdb_task(tf, tf, prot=True)
        run_tasks([task], ['run'])
        status = check_status(task)
        
        for ext in LASTDB_EXTENSIONS:
            assert os.path.isfile(tf + ext)

        assert status.status == 'up-to-date'


def test_lastdb_task_existing(tmpdir, datadir):
    with tmpdir.as_cwd():
        tf = datadir('test-protein.fa')
        for ext in LASTDB_EXTENSIONS:
            touch(tf + ext)

        task = lastdb_task(tf, tf, prot=True)
        run_tasks([task], ['run'])
        print(task, file=sys.stderr)
        status = check_status(task)

        assert status.status == 'up-to-date'


def test_lastal_task_nucl_x_prot(tmpdir, datadir):
    with tmpdir.as_cwd():
        prot = datadir('test-protein.fa')
        tr = datadir('test-transcript.fa')
        out = tmpdir.join('test-out').strpath

        db_task = lastdb_task(prot, prot)
        aln_task = lastal_task(tr, prot, out,  
                                translate=True, 
                                cutoff=None)
        run_tasks([db_task, aln_task], ['run'])

        aln = ''.join(open(out).readlines())
        print(aln, file=sys.stderr)

        assert 'SPAC212_RecQ_type_DNA_helicase_PROTEIN' in aln
        assert 'SPAC212_RecQ_type_DNA_helicase_TRANSCRIPT' in aln
        assert 'lambda' in aln, 'lambda missing, wrong LAST version?'
        

def test_lastal_task_prot_x_prot(tmpdir, datadir):
    with tmpdir.as_cwd():
        prot = datadir('test-protein.fa')
        out = tmpdir.join('test-out').strpath
            
        db_task = lastdb_task(prot, prot)
        aln_task = lastal_task(prot, prot, out,
                                translate=False,
                                cutoff=None)
        run_tasks([db_task, aln_task], ['run'])

        aln = ''.join(open(out).readlines())
        print(aln, file=sys.stderr)

        assert aln.count('SPAC212_RecQ_type_DNA_helicase_PROTEIN') == 2
        assert 'EG2=0' in aln
        assert 'E=0' in aln
        assert 'lambda' in aln, 'lambda missing, wrong LAST version?'


def test_lastal_task_multithreaded(tmpdir, datadir):
    with tmpdir.as_cwd():
        for n_threads in (3,4,5):
            prot = datadir('test-protein.fa')
            tr = datadir('pom.50.fa')
            out_single = tmpdir.join('out-single').strpath
            out_multi = tmpdir.join('out-multi').strpath

            db_task = lastdb_task(prot, prot)
            aln_task_single = lastal_task(tr, prot, out_single, 
                                           translate=True, 
                                           cutoff=None)

            aln_task_multi = lastal_task(tr, prot, out_multi,
                                         translate=True, 
                                         cutoff=None,
                                         n_threads=n_threads)
            run_tasks([db_task, aln_task_multi, aln_task_single], 
                      ['run'])

            alns_single = MafParser(out_single).read()
            alns_multi = MafParser(out_multi).read()

            assert all(alns_single['E'].sort_values() == \
                       alns_multi['E'].sort_values())


@pytest.mark.parametrize('n_threads', N_THREADS,
                         ids=lambda n: 'n_threads={0}'.format(n))
@pytest.mark.benchmark(group='lastal')
def test_lastal_task_large(datadir, lastdb_dir, tmpdir_factory, benchmark, n_threads):
    with tmpdir_factory.mktemp('THREADS_{0}'.format(n_threads)).as_cwd():
        query    = datadir('sacPom.cdna.fa')
        database = str(lastdb_dir.join('sacPom.pep.fa'))
        output   = 'out'
        
        aln_task = lastal_task(query, database, output,
                               translate=True,
                               cutoff=None,
                               n_threads=n_threads)
        result = benchmark.pedantic(run_tasks,
                                    args=([aln_task], ['run']),
                                    iterations=1,
                                    rounds=1)
        assert result == 0
        

def test_lastal_task_uptodate(tmpdir, datadir):
    with tmpdir.as_cwd():
        prot = datadir('test-protein.fa')
        out = tmpdir.join('test-out').strpath

        db_task = lastdb_task(prot, prot)
        aln_task = lastal_task(prot, prot, out,
                                translate=False,
                                cutoff=None)
        # Run it once
        run_tasks([db_task, aln_task], ['run'])
        # Now run again and check the status
        #run_tasks(aln_tasks, ['run'])
        print(aln_task)
        status = check_status(aln_task, tasks=[aln_task, db_task])
        assert status.status == 'up-to-date'

