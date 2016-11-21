import csv
import filelock
from functools import wraps
import sys
import time
import warnings
from contextlib import contextmanager
from os import path

from doit.action import PythonAction
from doit.task import Task as DoitTask
import six


def cleaned_actions(actions):
    '''Get a cleanup list of actions: Python actions have their <locals> portion
    stripped, which clutters up PythonActions that are closures.
    '''
    txt = ''
    for action in actions:
        txt_rep = six.text_type(action)
        if isinstance(action, PythonAction):
            # clean up inner fuctions in Python actions
            txt_rep = txt_rep.replace('<locals>.', '')
        else:
            txt_rep = txt_rep[:5] + '`' + txt_rep[5:] + '`'
        txt += "\n    * {0}".format(txt_rep)
    return txt


class Profiler(object):
    '''Thread-safe performance profiler.
    '''

    def __init__(self):
        self.running = False

    def start_profiler(self, filename=None, blockname='__main__'):
        '''Start the profiler, with results stored in the given filename.

        Args:
            filename (str): Path to store profiling results. If not given,
                uses a representation of the current time
            blockname (str): Name assigned to the main block.
        '''

        self.run_name = time.strftime("%a_%d_%b_%Y_%H%M%S", time.localtime())
        if filename is None:
            self.filename = '{0}.csv'.format(self.run_name)
        else:
            self.filename = filename
        self.run_name = time.ctime()
        self.start_time = time.time()
        self.blockname = blockname
        self.running = True
        self.lock = filelock.FileLock('{0}.lock'.format(self.filename))
        print('Profiling is ON:', self.filename, '\n', file=sys.stderr)

    def write_result(self, task_name, start_time, end_time, elapsed_time):
        '''Write results to the file, using the given task name as the
        name for the results block.

        Args:
            task_name (str): ID for the result row (the block profiled).
            start_time (float): Time of block start.
            end_time (float): Time of block end.
            elapsed_time (float): Total time.
        '''

        try:
            with self.lock.acquire(timeout=10):
                header = not path.isfile(self.filename)
                with open(self.filename, 'a') as fp:
                    writer = csv.writer(fp, delimiter=',')
                    if header:
                        writer.writerow(['run_id', 'block', 'start_t', 'end_t',
                                      'elapsed_t'])
                    row = [self.run_name, task_name, start_time, end_time, elapsed_time]
                    writer.writerow(row)
        except filelock.Timeout as e:
            warnings.warn(e, RuntimeWarning, stacklevel=1)

    def stop_profiler(self):
        '''Shut down the profiler and write the final elapsed time.
        '''
        self.end_time = time.time()
        elapsed = self.end_time - self.start_time
        self.write_result(self.blockname, self.start_time, self.end_time, elapsed)
        self.running = False
        return elapsed


class Timer(object):
    '''Simple timer class.
    '''

    def start(self):
        '''Start the timer.
        '''
        self.start_time = time.time()

    def stop(self):
        '''Stop the timer and return the elapsed time.
        '''
        self.end_time = time.time()
        return self.end_time - self.start_time


def title_without_profile_actions(task):
    """Generate title without profiling actions"""
    title = ''
    if task.actions:
        title = cleaned_actions(task.actions[1:-1])
    else:
        title = "Group: %s" % ", ".join(task.task_dep)
    return "%s: %s"% (task.name, title)


def setup_profiler():

    profiler = Profiler()
    @contextmanager
    def profiler_manager(filename=None, blockname='__main__'):
        profiler.start_profiler(filename=filename, blockname=blockname)
        yield
        profiler.stop_profiler()

    def add_profile_actions(task):
        timer = Timer()
        def start_profiling():
            if profiler.running:
                timer.start()

        def stop_profiling():
            if profiler.running:
                elapsed = timer.stop()
                profiler.write_result(task['name'], timer.start_time,
                                      timer.end_time, elapsed)
        
        if isinstance(task, DoitTask):
            actions = task._actions
            task.custom_title = title_without_profile_actions
        else:
            actions = task['actions']
            task['title'] = title_without_profile_actions

        actions.insert(0, start_profiling)
        actions.append(stop_profiling)

        return task
    
    def profile_decorator(task_func):

        @wraps(task_func)
        def func(*args, **kwargs):
            task = task_func(*args, **kwargs)
            return add_profile_actions(task)
        
        return func

    return profiler_manager, profile_decorator

StartProfiler, profile_task = setup_profiler()

