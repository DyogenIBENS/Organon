#!/usr/bin/env python3


# Based on the failure of `parsoire.py`, I am trying to reimplement a parsing
# module. It should be able to do:
#
# - Provide a simple syntax to transform any file content into a structured object
#   with keys/values (Nested if wanted).
# - use capturing regexes
# - convert types
# - iteratively process lines from the file, unit of parsing by unit of parsing.
# - be constrained by the order of units (unless an unordered option is given?)
# - handle multiline patterns.
# - handle conditionals: if X: parse unit 0; else parse unit 1.
# - handle optional units.
# - handle repeatable units.
# - handle units that should be discarded.
#
# The returned object should ideally be dictionary-like, or a nested namespace,
# with path like wildcards keys allowed.

import re
from collections import OrderedDict

class ParseError(BaseException):
    pass


class ParseBloc(object):
    """Aim to structure multiple ParseUnits/ParseBlocs"""
    def __init__(self, name, units=None):#, struct=OrderedDict):
        self.name = name
        self.units = units or []
        #self.struct = struct

    def parse(self, line_iterable):
        lines = iter(line_iterable)

        parsed = OrderedDict()

        for unit in self.units:
            parsed[unit.name], lines = unit.parse(lines)
            if not lines: break

        return parsed


class ParseUnit(object):
    def __init__(self, name, pattern, keys=None, types=None,
                 optional=False, repeat=False, flags=0,
                 fixed=False):
        self.name = name
        self.regex = pattern if fixed else re.compile(pattern, flags=flags)
        self.keys = keys or []
        self.types = types or []

        assert (self.regex.groups == len(self.types) or len(self.types) == 0) and \
               (self.regex.groups == len(self.keys) or len(self.keys) == 0)

        self.optional = optional
        self.repeat = repeat
        self.fixed = fixed
        
    def parse(self, text):
        if self.fixed:
            raise NotImplementedError
        
        if self.repeat:
            matches = self.regex.finditer(text)
        else:
            m = self.regex.search(text)
            matches = [m] if m else []
        
        if not matches:
            if self.optional:
                return [], text
            else:
                raise ParseError

        allparsed = []

        for m in matches:
            if self.regex.groups:
                groups = m.groups()
            else:
                groups = [m.group()]

            if self.types:
                groups = [typ(g) for typ,g in zip(self.types, groups)]

            if self.keys:
                groups = OrderedDict((key, g) for key,g in zip(self.keys, groups))
            elif len(groups) == 1:
                groups, = groups
            
            allparsed.append(groups)

        assert allparsed

        if len(allparsed) == 1:
            allparsed, = allparsed

        return allparsed, text[(m.span[1]+1):]




        


