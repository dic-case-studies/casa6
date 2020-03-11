from casatasks import casalog as _clog
from datetime import datetime as _time
import casatasks

def start_log( tname, arguments ):
    spaces = ' '*(18-len(tname))
    try:
        if casatasks.telemetry_enabled:
            casatasks.telemetrylogger.logger.info(str(_time.now()) + ' Begin Task: ' + tname + spaces)
    except:
        pass
    _clog.origin(tname)
    _clog.post( '##########################################' )
    _clog.post( '##### Begin Task: ' + tname + spaces + ' #####' )
    _clog.post( '%s( %s )' % ( tname, ', '.join(arguments) ))
    return str(_time.now()),

def end_log( state, tname, result ):
    spaces = ' '*(18-len(tname))
    try:
        if casatasks.telemetry_enabled:
            casatasks.telemetrylogger.logger.info(str(_time.now()) + ' End Task: ' + tname + spaces)
    except:
        pass
    _clog.origin(tname)
    _clog.post( 'Result %s: %s' % (tname, repr(result)) )
    _clog.post( 'Task ' + tname + ' complete. Start time: ' + state[0] + ' End time: ' + str(_time.now()) )
    _clog.post( '##### End Task: ' + tname + '  ' + spaces + ' #####' )
    _clog.post( '##########################################' )
    return result
