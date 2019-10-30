###
### This is the private module within casatasks which mediates access to
### user configuration options. It selects which user configuration file
### to load (only one should be loaded; if desired that file should
### load secondary files). config.py then pulls supported flags from from
### this file after validating them.
###
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

taskrc = os.path.join(home, ".casa/taskrc.py")
try:
    from casataskrc import *
except:
    try:
        f = open(taskrc)
    except IOError:
        pass
    else:
        f.close()
        try:
            exec(open(taskrc).read())
        except:
            import sys
            sys.stderr.write("error: evaluation of %s failed\n" % taskrc)
