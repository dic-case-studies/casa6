from __future__ import absolute_import
"""Hook to allow user-specified customization code to run.

This code was imported from Python 2.7 for use with CASAtools...

However, some programs or sites may find it convenient to allow users
to have a standard customization file, which gets run when a program
requests it.  This module implements such a mechanism.  A program
that wishes to use the mechanism must execute the statement

    import ctuser

The ctuser module looks for a file .casa/toolrc.py in the user's home
directory and if it can be opened, execfile()s it in its own global
namespace.  Errors during this phase are not caught; that's up to the
program that imports the ctuser module, if it wishes.

The user's .casa/toolrc.py could conceivably test for sys.version if it
wishes to do different things depending on the Python version.

"""
import os

home = os.curdir                        # Default
if 'HOME' in os.environ:
    home = os.environ['HOME']
elif os.name == 'posix':
    home = os.path.expanduser("~/")
elif os.name == 'nt':                   # Contributed by Jeff Bauer
    if 'HOMEPATH' in os.environ:
        if 'HOMEDRIVE' in os.environ:
            home = os.environ['HOMEDRIVE'] + os.environ['HOMEPATH']
        else:
            home = os.environ['HOMEPATH']

configrc = os.path.join(home, ".casa/config.py")
try:
    from casatoolrc import *
except:
    try:
        f = open(configrc)
    except IOError:
        pass
    else:
        f.close()
        try:
            exec(open(configrc).read( ))
        except:
            import sys
            sys.stderr.write("error: evaluation of %s failed\n" % configrc)
