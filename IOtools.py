#!/usr/bin/env python3


from sys import stdin, stdout


class Stream(object):
    """Context Manager class (for use with `with` statement).
    
    Do the exact same as `open()`, but if filename is "-" or None:
        - open stdout for writing, or
        - open stdin for reading.
    """
    write_modes = ('w', 'x', 'a')

    def __init__(self, filename=None, *args, **kwargs):
        self.filename = filename
        if filename is None or filename == '-':
            mode = args[0] if args else kwargs.get('mode', 'r')
            if any(letter in mode for letter in self.write_modes):
                self.stream = stdout
            else:
                self.stream = stdin
        else:
            self.stream = open(filename, *args, **kwargs)

    def __enter__(self):
        return self.stream

    def __exit__(self, type, value, traceback):
        if self.filename != '-':
            self.stream.close()
