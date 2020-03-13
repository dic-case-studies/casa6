###
### This is the casatasks configuration file. It loads user settings
### from userconfig.py. This file is available for use within casatasks
### as well as to users of the casatasks module.
###
import os as _os
import time as _time
from .private import userconfig as _uc
logfile = _os.path.realpath(_uc.logfile) if 'logfile' in dir(_uc) else _os.path.join(_os.getcwd( ),'casa-'+_time.strftime("%Y%m%d-%H%M%S", _time.gmtime())+'.log')
if 'telemetry_enabled' in dir(_uc):
    telemetry_enabled = _uc.telemetry_enabled
else:
    telemetry_enabled=True
