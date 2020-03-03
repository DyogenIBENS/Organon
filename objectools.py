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
    Iterate easily on all of its elements.

    Each "element" contained in Args is represented as an item (key, value):
    - if the element is in .args, then key is the index;
    - otherwise, it's the kwargs dictionary key.
    
    Therefore Args methods follow this logic:
    - items, __getitem__, __setitem__, pop.

    The Args.fromitems() methods requires "tuple-ized" items:
    >>> myargs = Args.fromitems((0, 'v0'), (1, 'v1'), ('a', 'va'), ('b', 'vb'))
    >>> print(myargs)
    (v0, v1, b='vb', a='va')
    
    But if key is an integer, the position is ignored, and item value appended.

    .args and .kwargs are unpacked separately:
    >>> [*myargs]  # Identical to [*myargs.args]
    ['1', '2']
    >>> {**myargs} # Identical to {**myargs.kwargs}
    {'a': 42, 'b': 55}

    CAUTION:
    The danger with this class is ambiguity in behavior when compared to python
    default behavior for list and dict:

    Example:
    - `x in myargs` checks whether x is in myargs.args or in myargs.kwargs.values()
    - but `myargs.remove(x)` removes the *value* x from args, but the *key* x from kwargs
      Use .pop for a consistent behavior.

    """
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
        """
        1st form: Args.fromcontent(args, kwargs)
        2d form: Args.fromcontent((args, kwargs))
        """
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
        #other = as_args(other)
        result = Args.fromcontent(self.args+other.args, self.kwargs)
        result.kwargs.update(other.kwargs)  # Therefore not symmetrical
        return result

    def __iadd__(self, other):
        #if isinstance(other, (tuple, list)):
        #other = as_args(other)
        self.args.extend(other.args)
        self.kwargs.update(other.kwargs)
        return self

    def __or__(self, other):  # Shouldn't it be __or__ ?
        #other = as_args(other)
        result = Args.fromcontent(self.args+[a for a in other if a not in self.args],
                                  self.kwargs)
        # This uses set logic on positional elements. I don't think this should be done in an
        # Args class
        result.kwargs.update(other.kwargs)
        return result

    def __ior__(self, other):
        #other = as_args(other)
        self.args.extend(other.args)
        self.kwargs.update(other.kwargs)
        return self

    #def __div__(self, other):
    #    self.

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
        #FIXME: should I also return range(len(self.args)) ?
        return self.kwargs.keys()

    def pop(self, *key_and_default):
        # I don't understand this implementation.
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
        """Convenience shortcut to Args.args.append"""
        self.args.append(value)


    def insert(self, index, value):
        """Convenience shortcut to Args.args.insert"""
        self.args.insert(index, value)

    # This failed attempt shows that you should do:
    #myargs += as_args(other_args_list_set_dict)
    #def update(self, *new_args, **new_kwargs):
    #    # Generate recursion error with self.update(self)
    #    self.args.extend(new_args)
    #    self.kwargs.update(**new_kwargs)
    #def extend(self, collection):


    #def push(self, *item):
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
        except ValueError:
            try:
                key, value = item[0]
            except ValueError:
                raise ValueError('Arguments should be key, value of (key, value)')
        if isinstance(key, int):
            self.args.append(value)
        else:
            self.kwargs[key] = value
    
    def popitem(self):
        try:
            return self.kwargs.popitem()
        except KeyError:
            return len(self.args)-1, self.args.pop()


    def remove_item(self, *item):
        """
        if key in kwargs, remove the *entry* (dictionary key)
        if value in args, remove the *value*.
        Raise error if not in any.
        """
        try:
            key, value = item
        except ValueError:
            try:
                key, value = item[0]
            except ValueError:
                raise ValueError('Arguments should be key, value of (key, value)')

        try:
            del self.kwargs[key]
        except KeyError:
            try:
                self.args.remove(value)
            except ValueError:
                raise ValueError('%r not in args and %s not in kwargs' % (value, key))


    def __str__(self):
        return (self.__class__.__name__ + '(' + ', '.join(
                    [repr(a) for a in self.args]
                    + [('%s=%r' % ka) for ka in self.kwargs.items()])
                + ')')

    def __repr__(self):
        return '<%s at 0x%x>' %(self, id(self))


def as_args(collection):
    """Convert list to Args, or dict to Args, and Args unchanged"""
    try:
        return Args.fromitems(*collection.items())
    except AttributeError:
        # collection is a list/tuple/set
        return Args(*collection)


def generic_update(current, new):
    """inplace, if element of new is already in current, do not add."""
    if isinstance(current, Args):
        current |= as_args(new)
        return
    #elif isinstance(current, (dict, set)):
    try:
        current.update(new)  # Don't try to mix sets and dicts here.
        return
    except AttributeError:
        pass
    #elif isinstance(new, (list, tuple)):
    try:
        current.extend((n for n in new if n not in current))
    except AttributeError:
        raise TypeError("Invalid type(current) = %s" % type(current))


def generic_extend(current, new):
    """inplace"""
    if isinstance(current, Args):
        current += as_args(new)
        return
    #elif isinstance(current, (dict, set)):
    try:
        current.update(new)  # Don't try to mix sets and dicts here.
        return
    except AttributeError:
        pass
    #elif isinstance(new, (list, tuple)):
    try:
        current.extend(new)
        return
    except AttributeError:
        raise TypeError("Invalid type(current) = %s" % type(current))


def generic_union(current, new):
    raise NotImplementedError
    if isinstance(current, dict):
        return {**current, **new}
    elif isinstance(current, set):
        return current | new
    elif isinstance(current, (list, tuple)):
        return current + type(new)((n for n in new if n not in current))
    elif isinstance(current, Args):
        return current + new  # Problem: adds element to .args (different from line above)


def generic_difference_update(current, new):  # Don't use.
    """"""
    raise NotImplementedError
    try:
        # If current is an Args
        for n in new:
            current.remove_item(new)
        return
    except AttributeError:
        pass

    try:
        # If current is a set
        current.difference_update(new)
        return
    except AttributeError:
        pass

    if isinstance(current, dict):
        for n in new:
            del current[n]
    elif isinstance(current, list):
        for n in new:
            try:
                current.remove(n)
            except ValueError:
                pass


def generic_remove_items(current, other):
    # This looks quite unsafe to me.
    #if isinstance(current, Args):
    try:
        for item in as_args(other):
            current.remove_item(item)
        return
    except AttributeError:
        pass
    try:
        for elem in other:
            current.remove(elem)
        return
    except AttributeError:
        pass
    #if isinstance(current, dict):
    for elem in other:
        del current[elem]


def generic_remove(current, element):
    """Remove item from Args, or element from set/list, or key from dict.
    Raises ValueError if not found."""
    try:
        # If current is an Args instance.
        current.remove_item(element)
        return
    except AttributeError:
        pass
    try:
        # If current is a set/list
        current.remove(element)
        return
    except AttributeError:
        pass

    # is a dict.
    del current[element]


#def generic_remove_item(current, *item):
