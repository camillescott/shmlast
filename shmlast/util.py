#!/usr/bin/env python

import hashlib
import os
from string import digits

from doit.tools import run_once, create_folder, title_with_actions, LongRunning
from doit.tools import PythonInteractiveAction
from doit.task import clean_targets, dict_to_task
from doit.cmd_base import TaskLoader
from doit.doit_cmd import DoitMain


def leftpad(s):
    return '\n'.join('    {0}'.format(i) for i in s.split('\n'))


def title(task):
    """return task name task actions"""
    if task.actions:
        title = "\n\t".join([leftpad(str(action)) for action in task.actions])
    # A task that contains no actions at all
    # is used as group task
    else:
        title = "Group: %s" % ", ".join(task.task_dep)
    return "%s => \n%s"% (task.name, title)


def create_doit_task(task_dict_func):
    '''Wrapper to decorate functions returning pydoit
    Task dictionaries and have them return pydoit Task
    objects
    '''
    def d_to_t(*args, **kwargs):
        ret_dict = task_dict_func(*args, **kwargs)
        return dict_to_task(ret_dict)
    return d_to_t


def task_str(task):
    return '{{ Task: {0}\n  actions: {1}\n  file_dep: {2}'\
           '\n  task_dep: {3}\n  targets: {4} }}'.format(task.name, task.actions,
                                                         task.file_dep, task.task_dep,
                                                         task.targets)


def prog_string(subcommand, version, action):
    s = ['shmlast {0} -- Camille Scott, 2016'.format(version)]
    s.append('\n')
    s.append(len(s[0]) * '-')
    s.append('\n')
    s.append('subcommand: {0}'.format(subcommand))
    s.append('\n')
    s.append('doit action: {0}'.format(action))
    return ''.join(s)


class ShortenedPythonAction(PythonInteractiveAction):

    def __str__(self):
        fullname = str(self.py_callable)[1:].split(' at ')[0]
        _, _, shortname = fullname.rpartition('.')
        return "Python: %s" % shortname


class DependencyError(RuntimeError):
    pass


def which(program, raise_err=True):
    '''Checks whether the given program (or program path) is valid and
    executable.

    NOTE: Sometimes copypasta is okay! This function came from stackoverflow:

        http://stackoverflow.com/a/377028/5109965

    Args:
        program (str): Either a program name or full path to a program.

    Returns:
        Return the path to the executable or None if not found
    '''
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    if raise_err:
        raise DependencyError('{0} not found; is it installed?'.format(program))
    else:
        return None


def parallel_fasta(input_filename, output_filename, command, n_jobs, pbs=None):

    exc = which('parallel')
    cmd = ['cat', input_filename, '|', exc, '--round-robin', '--pipe', '-L', 2,
           '-N', 10000, '--gnu']
    if pbs is not None:
        cmd.extend(['--sshloginfile', pbs, '--workdir $PWD'])
    else:
        cmd.extend(['-j', n_jobs])
    cmd.extend(['-a', input_filename])

    if isinstance(command, list):
        command = ' '.join(command)
    cmd.extend([command, '>', output_filename])
    return ' '.join(map(str, cmd))


def hidden_fn(fn):
    return '.{0}'.format(fn)
