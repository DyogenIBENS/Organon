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
        self.args = list(args)
        self.kwargs = kwargs

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

    def __nonzero__(self):
        return bool(self.args) or bool(self.kwargs)

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


    def __delitem__(self, key):
        if isinstance(key, (int, slice)):
            del self.args[key]
        else:
            del self.kwargs[key]

    @classmethod
    def fromitems(cls, *items):
        self = cls()
        for key, value in items:
            # or self.additem()
            if isinstance(key, int):
                self.args.append(value)
            else:
                self.kwargs[key] = value
        return self

    @classmethod
    def fromcontent(cls, content, kwargs=None):
        self = cls()
        if kwargs is not None:
            self.content = (list(content), kwargs)
        else:
            self.content = content
        return self

    #__getitem__ = dispatch_integerkey('__getitem__')
    #__setitem__ = dispatch_integerkey('__setitem__')
    #pop = dispatch_integerkey('pop')
    def __contains__(self, value):
        return value in self.args or value in self.kwargs.values()

    def __iter__(self):  # Iterate from the list only (so that *unpack works as expected)
        yield from self.args

    def __add__(self, other):
        result = Args.fromcontent(self.args+other.args, self.kwargs)
        result.kwargs.update(other.kwargs)  # Therefore not symmetrical
        return result

    def __iadd__(self, other):
        #if isinstance(other, (tuple, list)):
        try:
            self.args.extend(other.args)
            self.kwargs.update(other.kwargs)
        except AttributeError:
            self.args.extend(other[0])
            self.kwargs.update(other[1])

    def __len__(self):
        return len(self.args) + len(self.kwargs)


    def __call__(self, func, *args, **kwargs):
        """Call the given function with the stored arguments and kw-arguments"""
        return func(*args, *self.args, **kwargs, **self.kwargs)


    # Iterators
    def items(self):
        yield from enumerate(self.args)
        yield from self.kwargs.items()

    def values(self):
        yield from self.args
        yield from self.kwargs.values()

    def keys(self):
        return self.kwargs.keys()

    def pop(self, *key_and_default):
        if len(key_and_default)>2:
            raise ValueError("pop() takes at most 2 arguments.")
        try:
            key = key_and_default[0]
        except IndexError:
            return self.args.pop()

        if isinstance(key, int):
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
        """Only update the kwargs"""
        self.kwargs.update(*others, **others_kw)

    def additem(self, *item):
        """If key is an integer, simply ignore the position, and append."""
        #if len(item)==1:
        #    key, value = item[0]
        #elif len(item) == 2:
        #    key, value = item
        #else:
        #    raise ValueError('Arguments should be key, value of (key, value)')
        try:
            key, value = item
        except IndexError:
            try:
                key, value = item[0]
            except IndexError:
                raise ValueError('Arguments should be key, value of (key, value)')
        if isinstance(key, int):
            self.args.append(value)
        else:
            self.kwargs[key] = value

    def __str__(self):
        return ('(' + ', '.join(
                    [str(a) for a in self.args]
                    + [('%s=%r' % ka) for ka in self.kwargs.items()])
                + ')')

    def __repr__(self):
        return '<%s%s at 0x%x>' %(self.__class__.__name__, self, id(self))
