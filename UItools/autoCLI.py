#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Module to automatically build CLI wrapper scripts for any module/class"""


from __future__ import print_function

from sys import version_info, exit, stderr, stdin
python3 = version_info[0] >= 3

import re
import logging
logger = logging.getLogger(__name__)
#logging.basicConfig(format='%(levelname)s:%(funcName)s:%(message)s')


# For automatic type conversion of command-line arguments
int_regex = re.compile(r'^-?[0-9]+$')
float_regex = re.compile(r'^-?[0-9](\.[0-9])?([eE]-?[0-9]+)')


def make_subparser_func(func, commanddest='commands'):
    """Transform a normal function so that it takes arguments from the Argparse args."""
    def subp_func(args):
        dictargs = vars(args)
        dictargs.pop(commanddest)
        dictargs.pop('func')
        return func(**dictargs)
    return subp_func


def build_cli_processor(obj, arg_converters, argstart=0, command=None, *args):
    """Typically use argstart=1 for class methods (removes `self`)"""
    help_opts = set(('-h', '--help', '-?'))

    if command is None or command in help_opts:
        validmethods = ['-h/--help/-?', 'pass'] + \
                       [m for m in dir(obj) if
                        callable(getattr(obj, m)) and not m.startswith("__")]
        print(__doc__ + '\nCommands:\n - ' + '\n - '.join(validmethods), file=stderr)
        exit()
    elif command == 'pass':
        processor = lambda self, *args: None
    else:
        processor = getattr(obj, command, None)

    if processor is None:
        raise ValueError('Method %r not found in %r' % (command, obj))

    proc_code = processor.__code__ if python3 else processor.func_code
    proc_args = proc_code.co_varnames[:proc_code.co_argcount]
    proc_defaults = [] if processor.__defaults__ is None else processor.__defaults__
    proc_doc = '' if processor.__doc__ is None else processor.__doc__

    n_mandatory = len(proc_args) - len(proc_defaults)
    proc_cli_args = [(a,) for a in proc_args[:n_mandatory]] + \
                    list(zip(proc_args[n_mandatory:], proc_defaults))

    #def cli_processor(args):
    help_text = "%s\n---\n%s\nArgs:\n- %s" % (command,
                                proc_doc,
                                "\n- ".join(("%s" % a) if len(a)==1
                                            else ("%s [%r]" % a)
                                            for a in proc_cli_args))
    
    not_enough_args = len(args) < len(proc_args) - len(proc_defaults)
    too_many_args = len(args) > len(proc_args)
    help_asked = help_opts.intersection(args)
    #if not args or args[0] in help_opts:
    if not_enough_args or too_many_args or help_asked:
        print(help_text)
        if not help_asked:
            exit_code = 1
            if not_enough_args:
                logger.error("Not enough arguments")
            elif too_many_args:
                logger.error("Too many arguments")
        else:
            exit_code = 0
        exit(exit_code)

    args = args[argstart:]
    proc_args = proc_args[argstart:]
    
    converted_args = []
    for expected, arg in zip(proc_args, args):
        
        # Default type conversion
        if arg == '-':
            arg = stdin
        elif arg.lower() == 'false':
            arg = False
        elif arg.lower() == 'true':
            arg = True
        elif arg.lower() == 'none':
            arg = None
        elif int_regex.match(arg):
            arg = int(arg)
        elif float_regex.match(arg):
            arg = float(arg)

        # Conversion based on the expected arg name.
        convert = arg_converters.get(expected.lower())
        if convert is not None:
            arg = convert(arg)

        converted_args.append(arg)

    return processor, converted_args
