#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Convenience classes such as namespaces"""


#import itertools as it
#from collections import OrderedDict
#from collections import abc
#from functools import singledispatch, wraps


#def dispatch_integerkey(methodname):
#    @wraps(method)
#    def newmethod(self, key, *a, **kw):
#        if isinstance(key, int):
#            return getattr(self.args, methodname)(key, *a, **kw)
#        else:
#            return getattr(self.kwargs, methodname)(key, *a, **kw)
#    return newmethod


class Args(object):
    """Object to hold unpacked arguments (*args, **kwargs).
    Iterate easily on all of its elements."""
    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

    def items(self):
        yield from enumerate(self.args)
        yield from self.kwargs.items()

    def values(self):
        yield from self.args
        yield from self.kwargs.values()

    def keys(self):
        return self.kwargs.keys()

    def __getitem__(self, key):
        try:
            return self.args[key]
        except TypeError:
            return self.kwargs[key]

    def __setitem__(self, key, value):
        if isinstance(key, int):
            self.args[key] = value
        else:
            self.kwargs[key] = value

    #__getitem__ = dispatch_integerkey('__getitem__')
    #__setitem__ = dispatch_integerkey('__setitem__')
    #pop = dispatch_integerkey('pop')
    def __contains__(self, value):
        return value in self.args or value in self.kwargs.values()

    def __iter__(self):  # Iterate from the list only (so that *unpack works as expected)
        yield from self.args

    @property
    def content(self):
        return (self.args, self.kwargs)

    @content.setter
    def content(self, content):
        self.args = [*content[0]]
        self.kwargs = dict(**content[1])

    @content.deleter
    def content(self):
        self.args = []
        self.kwargs = {}


    def __call__(self, func):
        """Call the given function with the stored arguments and kw-arguments"""
        return func(*self.args, **self.kwargs)


    def pop(self, *key_and_default):
        if len(key_and_default)>2:
            raise ValueError("pop() takes at most 2 arguments.")
        try:
            key = key_and_default[0]
        except IndexError:
            return self.args.pop()

        if insinstance(key, int):
            return self.args.pop(key)
        else:
            try:
                default = key_and_default[1]
            except IndexError:
                return self.kwargs.pop(key)
            return self.kwargs.pop(key, default)


    def append(self, value):
        self.args.append(value)


    def insert(self, index, value):
        self.args.insert(index, value)


    def update(self, *others, **others_kw):
        self.kwargs.update(*others, **others_kw)


    def __str__(self):
        return ('(' + ', '.join(
                    [str(a) for a in self.args]
                    + [('%s=%r' % ka) for ka in self.kwargs.items()])
                + ')')

    def __repr__(self):
        return '<%s%s at 0x%x>' %(self.__class__.__name__, self, id(self))
