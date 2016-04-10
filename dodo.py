#!/usr/bin/env doit

def setupcmd(cmd):
    return ' '.join(['python', 'setup.py'] + cmd)

def task_install():
    
    return {'actions': [setupcmd(['install'])]}

