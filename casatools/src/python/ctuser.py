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
import sys
from contextlib import contextmanager

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

def _fileno(file_or_fd):
    fd = getattr(file_or_fd, 'fileno', lambda: file_or_fd)()
    if not isinstance(fd, int):
        raise ValueError("Expected a file (`.fileno()`) or a file descriptor")
    return fd

@contextmanager
def _stdout_redirected(to=os.devnull, stdout=None):
    if stdout is None:
       stdout = sys.stdout

    stdout_fd = _fileno(stdout)
    # copy stdout_fd before it is overwritten
    #NOTE: `copied` is inheritable on Windows when duplicating a standard stream
    with os.fdopen(os.dup(stdout_fd), 'wb') as copied: 
        stdout.flush()  # flush library buffers that dup2 knows nothing about
        try:
            os.dup2(_fileno(to), stdout_fd)  # $ exec >&to
        except ValueError:  # filename
            with open(to, 'wb') as to_file:
                os.dup2(to_file.fileno(), stdout_fd)  # $ exec > to
        try:
            yield stdout # allow code to be run with the redirected stdout
        finally:
            # restore stdout to its previous value
            #NOTE: dup2 makes stdout_fd inheritable unconditionally
            stdout.flush()
            os.dup2(copied.fileno(), stdout_fd)  # $ exec >&copied

@contextmanager
def _redirect_stdout(new_target):
    old_target, sys.stdout = sys.stdout, new_target # replace sys.stdout
    try:
        yield new_target # run some code with the replaced stdout
    finally:
        sys.stdout = old_target # restore to the previous value

def _merged_stderr_stdout():  # $ exec 2>&1
    return _stdout_redirected(to=sys.stdout, stdout=sys.stderr)


if len(sys.argv) > 0 and sys.argv[0] == '-m':
    ##  many packages use the idiom "python -m casatools --some-flag" to query for
    ##  casatools state. Output from user files causes problems for this.
    with _stdout_redirected(to=os.devnull), _merged_stderr_stdout():
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
                    sys.stderr.write("error: evaluation of %s failed\n" % configrc)
else:
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
                sys.stderr.write("error: evaluation of %s failed\n" % configrc)
