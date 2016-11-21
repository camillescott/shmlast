#!/usr/bin/env doit

DOIT_CONFIG = {'verbosity': 2}

def setupcmd(cmd):
    return ' '.join(['python', 'setup.py'] + cmd)

def task_install():
    
    return {'actions': [setupcmd(['install'])]}

def task_test():
    return {'actions': [setupcmd(['test'])]}

def task_publish():
    return {'actions': ['rm -rf shmlast/tests/__pycache__',
                        'rm -rf shmlast/__pycache__',
                        'pandoc README.md -o README.rst',
                        setupcmd(['sdist', 'upload'])],
            'file_dep': ['shmlast/VERSION'],
            'targets': ['README.rst']}
