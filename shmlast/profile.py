import csv
import filelock
import sys
import time
import warnings
from contextlib import contextmanager
from os import path


class Profiler(object):

    def __init__(self):
        self.running = False

    def start_profiler(self, filename=None, blockname='__main__'):
        self.run_name = time.strftime("%a_%d_%b_%Y_%H%M%S", time.localtime())
        if filename is None:
            self.filename = '{0}.csv'.format(self.run_name)
        else:
            self.filename = filename
        self.run_name = time.ctime()
        self.start_time = time.time()
        self.blockname = blockname
        self.running = True
        self.lock = filelock.Filelock('{0}.lock'.format(self.filename))
        print('Profiling is ON:', self.filename, '\n', file=sys.stderr)

    def write_result(self, task_name, start_time, end_time, elapsed_time):
        try:
            with lock.acquire(timeout=10):
                header = not path.isfile(self.filename)
                with open(self.filename, 'a') as fp:
                    writer = csv.writer(fp, delimiter=',')
                    if header:
                        writer.write(['run_id', 'block', 'start_t', 'end_t',
                                      'elapsed_t'])
                    row = [self.run_name, task_name, start_time, end_time, elapsed_time]
                    writer.writerow(row)
        except filelock.Timeout as e:
            warnings.warn(e, RuntimeWarning, stacklevel=1)

    def stop_profiler(self):
        self.end_time = time.time()
        elapsed = self.end_time - self.start_time
        self.write_result(self.blockname, self.start_time, self.end_time, elapsed)
        self.running = False
        return elapsed


class Timer(object):

    def start(self):
        self.start_time = time.time()

    def stop(self):
        self.end_time = time.time()
        return self.end_time - self.start_time


def setup_profiler():

    profiler = Profiler()
    @contextmanager
    def profiler_manager(filename=None, blockname='__main__'):
        profiler.start_profiler(filename=filename, blockname=blockname)
        yield
        profiler.stop_profiler()

    def profile_decorator(task_dict_func):

        def func(*args, **kwargs):
            task = task_dict_func(*args, **kwargs)
            if profiler.running:
                timer = Timer()
                def start_profiling():
                    timer.start()

                def stop_profiling():
                    elapsed = timer.stop()
                    profiler.write_result(task['name'], timer.start_time,
                                         timer.end_time, elapsed)
                
                task['actions'].insert(0, start_profiling)
                task['actions'].append(stop_profiling)

            return task
        
        return func

    return profiler_manager, profile_decorator

StartProfiler, profile_task = setup_profiler()

