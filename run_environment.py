#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Module to display various informations about files and local environment."""
# Alternative module name: silicotope, silicodrome


import sys
import os
import time
import os.path as op
import subprocess
import inspect
import socket

DATE_FORMAT = '%Y-%m-%d %H:%M:%S'

def get_caller_module():
    return inspect.getframeinfo(inspect.getouterframes(inspect.currentframe())[1][0])[0]


def get_git_repo(module=None, modulepath=None):
    if modulepath is None:
        if module is None:
            module = sys.modules[__name__]

        moduledir, mfile = op.split(op.abspath(module.__file__))  #op.realpath
    else:
        moduledir, mfile = op.split(op.expandvars(op.expanduser(modulepath)))
    repo, _ = run_git_command(['rev-parse', '--show-toplevel'], moduledir)
    repo = repo.strip()

    return repo, op.join(moduledir.replace(repo, ''), mfile)


def run_git_command(args, moduledir, timeout=10):
    """Run a git command by moving to the appropriate directory first.
    By default the directory is this module directory, or a given loaded module.
    """
    p = subprocess.Popen(['git', '--no-pager'] + args,
                         cwd=moduledir,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)

    try:
        out, err = p.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        p.kill()
        out, err = p.communicate()
        raise
    return out.decode(), err


def get_git_commit(moduledir, mfile=None, timeout=10):
    # '%x09' is a tab
    args = ['log', '-1',
            '--date=format:'+DATE_FORMAT,
            '--format=%h\t%ad\t%<(70,trunc)%s']
    if mfile:
        args.append(mfile)

    out, err = run_git_command(args, moduledir, timeout)

    if err:
        print(err, file=sys.stderr)
    return out.rstrip().split('\t', maxsplit=2)


def print_git_commit(*args, sep='\n', **kwargs):
    print(sep.join(get_git_commit(*args, **kwargs)))


def get_unstaged_changed(moduledir, timeout=10):
    # Interesting diff options: --name-status, --name-only
    out, err = run_git_command(['diff', '--name-status'], moduledir, timeout)
    if err:
        print(err, file=sys.stderr)
    if out:
        outlines = [(line, op.getmtime(op.join(moduledir, line.split('\t')[1])))
                    for line in out.rstrip().split('\n')]
        outlines.sort(key=lambda v: v[1])
    else:
        outlines = []
    return [line + '\t' + time.strftime(DATE_FORMAT, time.localtime(t))
            for line,t in outlines]


def get_staged_changed(moduledir, timeout=10):
    out, err = run_git_command(['diff', '--name-status', '--staged'],
                               moduledir, timeout)
    if err:
        print(err, file=sys.stderr)
    return out.rstrip().split('\n') if out else []


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


if __name__=='__main__':
    if len(sys.argv)==2 and sys.argv[1] == 'git_state':
        print_git_state(modulepath=(os.getcwd()+'/'))
    else:
        print('USAGE: ./run_environment.py git_state')
        sys.exit(2)

