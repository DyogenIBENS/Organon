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
import logging
logger = logging.getLogger(__name__)
#logging.basicConfig(format="%(levelname)s:%(module)s:%(message)s")


class ParseUnitNotFoundError(LookupError):
    pass


def floatlist(line):
    return [float(x) for x in line.split()]


def dictfloatlist(line):
    key, *fields = line.split()
    return {key: [float(x) for x in fields]}


def strlist(line):
    return line.split()


class ParseBloc(object):
    """Aim to structure multiple ParseUnits/ParseBlocs"""
    def __init__(self, name, units=None):#, repeat=False, struct=OrderedDict):
        self.name = name
        self.units = units or []
        #self.struct = struct

    def parse(self, text):

        parsed = OrderedDict()

        for unit in self.units:
            parsed[unit.name], text = unit.parse(text)
            if not text: break

        return parsed, text


class ParseUnit(object):
    def __init__(self, name, pattern, keys=None, types=None,
                 optional=False, repeat=False, merge=False, flags=re.M,
                 fixed=False, convert=None):
        self.name = name
        if fixed:
            pattern = re.escape(pattern)
        self.regex = re.compile(pattern, flags=flags)
        self.keys = keys or []
        self.types = types or []
        self.convert = convert

        assert (self.regex.groups == len(self.types) or len(self.types) == 0) and \
               (self.regex.groups == len(self.keys) or len(self.keys) == 0)

        #assert merge is False or self.regex.groups > 1 and not self.keys, name

        self.optional = optional
        self.repeat = repeat
        self.merge = merge
        #self.fixed = fixed
        
    def parse(self, text):
        if self.repeat:
            matches = self.regex.finditer(text)
        else:
            m = self.regex.search(text)
            matches = [m] if m else []
        
        if not matches:
            if self.optional:
                return [], text
            else:
                lines = text.split('\n')
                nl = len(lines)
                lines = lines[:min(2, nl)]
                if nl > 4:
                    lines.append('[...]')
                lines.extend(lines[max(nl-2, nl//2+1):])
                raise ParseUnitNotFoundError("%s: '%s' in:\n%s" %(self.name,
                                                        self.regex.pattern,
                                                        '\n'.join(lines)))

        allparsed = []

        for m in matches:
            if self.regex.groups:
                groups = m.groups()
            else:
                groups = []
                # Would be more interesting not to store things when no capturing parentheses

            if self.types:
                groups = [typ(g) for typ,g in zip(self.types, groups)]

            if self.keys:
                groups = OrderedDict((key, g) for key,g in zip(self.keys, groups))
            
            if self.merge or len(groups) <= 1:
                if self.keys:
                    logger.warning('Useless key for single value. Discarding the key.')
                allparsed.extend(groups) #merge should make no difference if len(groups) was 1
            else:
                allparsed.append(groups)

        #assert allparsed

        if len(allparsed) == 1:
            allparsed, = allparsed

        if self.convert:
            allparsed = self.convert(allparsed)

        assert m is not None and m.span() is not None, self.name

        remainingtext = text[(m.span()[1]+1):]
        return allparsed, remainingtext




        


