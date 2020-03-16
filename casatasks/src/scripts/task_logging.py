from casatasks import casalog as _clog
from datetime import datetime as _time
import casatasks

def start_log( tname, arguments ):
    spaces = ' '*(18-len(tname))
    start_time = str(_time.now())
    try:
        if casatasks.config.telemetry_enabled:
            casatasks.telemetrylogger.logger.info(start_time + ' Begin Task: ' + tname + spaces)
    except:
        pass
    _clog.origin(tname)
    _clog.post( '##########################################' )
    _clog.post( '##### Begin Task: ' + tname + spaces + ' #####' )
    _clog.post( '%s( %s )' % ( tname, ', '.join(arguments) ))
    return start_time,

def end_log( state, tname, result ):
    spaces = ' '*(18-len(tname))
    end_time = str(_time.now())
    try:
        if casatasks.config.telemetry_enabled:
            casatasks.telemetrylogger.logger.info(end_time + ' End Task: ' + tname + spaces + 'Start time: ' + state[0])
    except:
        pass
    _clog.origin(tname)
    _clog.post( 'Result %s: %s' % (tname, repr(result)) )
    _clog.post( 'Task ' + tname + ' complete. Start time: ' + state[0] + ' End time: ' + end_time )
    _clog.post( '##### End Task: ' + tname + '  ' + spaces + ' #####' )
    _clog.post( '##########################################' )
    return result
