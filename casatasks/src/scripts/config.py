###
### This is the casatasks configuration file. It loads user settings
### from userconfig.py. This file is available for use within casatasks
### as well as to users of the casatasks module.
###
import os as _os
import time as _time
from .private import userconfig as _uc
logfile = _os.path.realpath(_uc.logfile) if 'logfile' in dir(_uc) else _os.path.join(_os.getcwd( ),'casa-'+_time.strftime("%Y%m%d-%H%M%S", _time.gmtime())+'.log')
rcdir = _os.path.realpath(_uc.rcdir) if 'rcdir' in dir(_uc) else None

# Telemetry and CrashReporter
telemetry_enabled = _uc.telemetry_enabled if 'telemetry_enabled' in dir(_uc) else True
crashreporter_enabled = _uc.crashreporter_enabled if 'crashreporter_enabled' in dir(_uc) else True
telemetry_log_directory = _uc.telemetry_log_directory if 'telemetry_log_directory' in dir(_uc) else None
telemetry_log_limit = _uc.telemetry_log_limit if 'telemetry_log_limit' in dir(_uc) else None
telemetry_log_size_interval = _uc.telemetry_log_size_interval if 'telemetry_log_size_interval' in dir(_uc) else None
telemetry_submit_interval = _uc.telemetry_submit_interval if 'telemetry_submit_interval' in dir(_uc) else None
