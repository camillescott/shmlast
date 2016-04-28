#!/usr/bin/env python

import hashlib
import os
from string import digits

from doit.tools import run_once, create_folder, title_with_actions, LongRunning
from doit.tools import PythonInteractiveAction
from doit.task import clean_targets, dict_to_task
from doit.cmd_base import TaskLoader
from doit.doit_cmd import DoitMain


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
    s.append('\n\n')
    return ''.join(s)


def run_tasks(tasks, args, config={'verbosity': 2}):
   
    tasks = list(tasks)

    class Loader(TaskLoader):
        @staticmethod
        def load_tasks(cmd, opt_values, pos_args):
            return tasks, config
   
    return DoitMain(Loader()).run(args)


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

def filesize_commands(input_filename, n_jobs, n_nodes=None):
    if n_nodes is not None:
        n_jobs = n_jobs * n_nodes

    var=hashlib.sha224(input_filename.encode('utf-8')).hexdigest().translate({ord(k): None for k in digits})
    cmds = []
    cmds.append("export {var}_size=`du --apparent-size --block-size=1 pep.fa 2>/dev/null | awk {{'print $1'}}`".format(var=var))
    cmds.append("export {var}_block=`expr ${var}_size / {j}`".format(j=n_jobs, var=var))
    
    return cmds, "{var}_block".format(var=var)

def parallel_fasta(input_filename, n_jobs):
    cmds, block_var = filesize_commands(input_filename, n_jobs)

    exc = which('parallel')
    cmd = ['cat', input_filename, '|', exc, '--block', '$'+block_var,
           '--eta', '--pipe', '--recstart', '">"', '--gnu', '-j', str(n_jobs)]

    return cmds, ' '.join(cmd)

def multinode_parallel_fasta(input_filename, ppn, nodes):
    cmds, block_var = filesize_commands(input_filename, ppn, nodes)
    
    exc = which('parallel')
    cmd = ['cat', input_filename, '|', exc, '--block', '$'+block_var,
                 '--eta', '--pipe', '--recstart', '">"', '--gnu', '--jobs', str(ppn),
                 '--sshloginfile $PBS_NODEFILE', '--workdir $PWD']

    return cmds, ' '.join(cmd)

