###
### This is the casatasks configuration file. It loads user settings
### from userconfig.py. This file is available for use within casatasks
### as well as to users of the casatasks module.
###
import os as _os
import time as _time
from .private import userconfig as _uc

logfile = _os.path.join(_os.getcwd( ),'casa-'+_time.strftime("%Y%m%d-%H%M%S", _time.gmtime())+'.log')
if 'logfile' in dir(_uc) and type(_uc.logfile) is str:
    if _uc.logfile == "/dev/null":
        logfile="/dev/null"
    elif _os.path.exists(_uc.logfile):
        if _os.access(_uc.logfile,_os.W_OK):
            logfile = _uc.logfile
    elif _os.path.isdir(_os.path.dirname(_uc.logfile)) and \
         _os.access(_os.path.dirname(_uc.logfile),_os.W_OK):
        logfile = _os.path.realpath(_uc.logfile)
