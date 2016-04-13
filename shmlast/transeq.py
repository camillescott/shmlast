#!/usr/bin/env python3

from doit.task import clean_targets
from doit.tools import title_with_actions
import pandas as pd
import screed

from .util import create_doit_task as doit_task
from .util import ShortenedPythonAction


@doit_task
def rename_task(input_fn, output_fn, name_map_fn='name_map.csv'):
    
    def rename_input():
        name_map = []
        with open(output_fn, 'w') as output_fp:
            for n, record in enumerate(screed.open(input_fn)):
                new_name = 'tr{0}'.format(n)
                output_fp.write('>{0}\n{1}\n'.format(new_name,
                                                     record.sequence))
                name_map.append((record.name, new_name))

        pd.DataFrame(name_map,
                     columns=['old_name', 'new_name']).to_csv(name_map_fn,
                                                              index=False)

    return {'name': 'rename_sequences',
            'title': title_with_actions,
            'actions': [ShortenedPythonAction(rename_input)],
            'targets': [output_fn, name_map_fn],
            'file_dep': [input_fn],
            'clean': [clean_targets]}


@doit_task
def transeq_task(input_fn, output_fn, clean=True, frame=6):

    params = '-frame {frame}'.format(frame=frame)
    if clean:
        params = '{params} -clean'.format(params=params)

    cmd = 'transeq -sequence {input_fn} {params}'\
          ' -outseq {output_fn}'.format(**locals())

    return {'name': 'transeq',
            'title': title_with_actions,
            'actions': [cmd],
            'targets': [output_fn],
            'file_dep': [input_fn],
            'clean': [clean_targets]}

