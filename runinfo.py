#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Module to display various informations about files and local environment"""


import sys
import os.path as op
import subprocess
import inspect
import socket


def get_caller_module():
    return inspect.getframeinfo(inspect.getouterframes(inspect.currentframe())[1][0])[0]


def get_git_commit(module=None, modulepath=None, withfile=False, timeout=10):
    if modulepath is None:
        if module is None:
            module = sys.modules[__name__]

        moduledir, mfile = op.split(op.abspath(module.__file__))  #op.realpath
    else:
        moduledir, mfile = op.split(op.expandvars(op.expanduser(modulepath)))
    
    # '%x09' is a tab
    args = ['git', 'log', '-1',
            '--date=format:%Y-%m-%d %H:%M:%S',
            '--format=%h\t%ad\t%<(70,trunc)%s']
    if withfile:
        args.append(mfile)
    p = subprocess.Popen(args,
                         cwd=moduledir,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)

    try:
        out, err = p.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        p.kill()
        out, err = p.communicate()
        raise

    if err:
        print(err, file=sys.stderr)
    return  ['%s:%s' % (socket.gethostname(), moduledir)] + out.decode().rstrip().split('\t', maxsplit=2)


def print_git_commit(*args, sep='\n', **kwargs):
    print(sep.join(get_git_commit(*args, **kwargs)))
