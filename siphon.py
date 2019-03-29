


from functools import partial, wraps

def dependency_func(firstget):
    dependency_name = firstget.__name__
    attr_name = '_' + dependency_name

    @wraps(firstget)
    def fget(self, *args, **kwargs):
        try:
            return getattr(self, attr_name)
        except AttributeError:
            #print('Setting attribute %r' % attr_name)
            value = firstget(self, *args, **kwargs)
            setattr(self, attr_name, value)
            return value
    return fget


def dependency(firstget):
    """Decorator to create special property, that gets defined only at the first access."""
    return property(dependency_func(firstget))


def auto_internalmethod(func, argkeys=None, kwargkeys=None):
    """Automatically call func as a method, by taking the arguments from
    the attributes of self."""
    varnames = func.__code__.co_varnames
    argnames = varnames[:func.__code__.co_argcount]
    ndefaults = len(func.__defaults__) if func.__defaults__ else 0
    mandatory_args = argnames[:-ndefaults] if ndefaults else argnames
    optional_args = argnames[-ndefaults:]
    
    @wraps(func)
    def internalmethod(self):
        internal_args = [getattr(self, v) for v in mandatory_args]
        for opt, default in zip(optional_args, func.__defaults__):
            internal_args.append(getattr(self, opt, default))
        if argkeys is not None:
            internal_args += argkeys
        internal_kw = {} if kwargkeys is None else {v: getattr(self, v)
                                                    for v in kwargkeys
                                                    if hasattr(self, v)}
        return func(*internal_args, **internal_kw)
    return internalmethod


# TODO: optional attr_name. Check that no two properties access the same attr.
#def tupledependency(firstget, n=1):
#    """firstget must return a tuple/iterable. Then, each item will be a property.
#    # Don't use lambdas! (they are all named "<lambda>")
#    """
#    # Could use n=func.__code__.__annotations__['return'] is exists.
#    dependency_name = firstget.__name__
#    attr_name = '_' + dependency_name
#
#    def item_fget(i, self, *args, **kwargs):
#        try:
#            return getattr(self, attr_name)[i]
#        except AttributeError:
#            print('Setting attribute %r', attr_name)
#            #print('Setting attribute %r' % attr_name)
#            value = firstget(self, *args, **kwargs)
#            setattr(self, attr_name, value)
#            return value[i]
#
#    #partialmethods?
#    return tuple(property(partial(item_fget, i)) for i in range(n))

