#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Module to display various informations about files and local environment."""
# Alternative module name: silicotope


import sys
import os
import os.path as op
import subprocess
import inspect
import socket


def get_caller_module():
    return inspect.getframeinfo(inspect.getouterframes(inspect.currentframe())[1][0])[0]


def get_git_repo(module=None, modulepath=None):
    if modulepath is None:
        if module is None:
            module = sys.modules[__name__]

        moduledir, mfile = op.split(op.abspath(module.__file__))  #op.realpath
    else:
        moduledir, mfile = op.split(op.expandvars(op.expanduser(modulepath)))

    return moduledir, mfile


def run_git_command(args, moduledir, timeout=10):
    """Run a git command by moving to the appropriate directory first.
    By default the directory is this module directory, or a given loaded module.
    """
    p = subprocess.Popen(['git'] + args,
                         cwd=moduledir,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)

    try:
        out, err = p.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        p.kill()
        out, err = p.communicate()
        raise
    return out, err


def get_git_commit(moduledir, mfile=None, timeout=10):
    # '%x09' is a tab
    args = ['log', '-1',
            '--date=format:%Y-%m-%d %H:%M:%S',
            '--format=%h\t%ad\t%<(70,trunc)%s']
    if mfile:
        args.append(mfile)

    out, err = run_git_command(args, moduledir, timeout)

    if err:
        print(err, file=sys.stderr)
    return out.decode().rstrip().split('\t', maxsplit=2)


def print_git_commit(*args, sep='\n', **kwargs):
    print(sep.join(get_git_commit(*args, **kwargs)))


def get_unstaged_changed(moduledir, timeout=10):
    out, err = run_git_command(['diff', '--name-only'], moduledir, timeout)
    if err:
        print(err, file=sys.stderr)
    return out.decode().rstrip().split('\n')


def get_staged_changed(moduledir, timeout=10):
    out, err = run_git_command(['diff', '--name-only', '--staged'],
                               moduledir, timeout)
    if err:
        print(err, file=sys.stderr)
    return out.decode().rstrip().split('\n')


def print_git_state(module=None, modulepath=None, sep='\n', timeout=10):
    moduledir, mfile = get_git_repo(module, modulepath)
    state = ['%s:%s' % (socket.gethostname(), moduledir), '-'*50]
    state += get_git_commit(moduledir, timeout=timeout) + ['']
    state += ['# File %s' % mfile] + get_git_commit(moduledir, mfile, timeout=timeout) + ['']
    state += ['# Staged changes in:'] + get_staged_changed(moduledir, timeout) + ['']
    state += ['# Unstaged changes in:'] + get_unstaged_changed(moduledir, timeout) + ['']
    print(sep.join(state))



def redisplay():
    """Reset the correct DISPLAY environment variable, when using `tmux` over
    `ssh`."""
    correct_localhost = subprocess.check_output(['tmux', 'show-env', 'DISPLAY'])\
                            .decode()\
                            .replace('DISPLAY=', '')\
                            .rstrip()
    print('%s -> %s' % (os.environ['DISPLAY'], correct_localhost))
    os.environ['DISPLAY'] = correct_localhost
