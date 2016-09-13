#!/usr/bin/env python3

from doit.task import clean_targets
import pandas as pd
import screed

from .profile import profile_task
from .util import create_doit_task as doit_task
from .util import ShortenedPythonAction, title, unwrap_fasta, which


@doit_task
@profile_task
def rename_task(input_fn, output_fn, name_map_fn='name_map.csv', prefix='tr'):
    
    def rename_input():
        name_map = []
        with open(output_fn, 'w') as output_fp:
            for n, record in enumerate(screed.open(input_fn)):
                new_name = '{0}{1}'.format(prefix, n)
                output_fp.write('>{0}\n{1}\n'.format(new_name,
                                                     record.sequence))
                name_map.append((record.name, new_name))

        pd.DataFrame(name_map,
                     columns=['old_name', 'new_name']).to_csv(name_map_fn,
                                                              index=False)

    return {'name': 'rename:{0}'.format(input_fn),
            'title': title,
            'actions': [ShortenedPythonAction(rename_input)],
            'targets': [output_fn, name_map_fn],
            'file_dep': [input_fn],
            'clean': [clean_targets]}


@doit_task
@profile_task
def transeq_task(input_fn, output_fn, clean=True, frame=6):

    exc = which('transeq')
    unwrap = unwrap_fasta()
    params = '-frame {frame}'.format(frame=frame)
    if clean:
        params = '{params} -clean'.format(params=params)
    
    cmd = [exc, '-sequence', input_fn, params, '-outseq',
           '/dev/stdout', '|', unwrap, '>', output_fn]
    cmd = ' '.join(cmd)

    return {'name': 'transeq:{0}'.format(input_fn),
            'title': title,
            'actions': [cmd],
            'targets': [output_fn],
            'file_dep': [input_fn],
            'clean': [clean_targets]}

