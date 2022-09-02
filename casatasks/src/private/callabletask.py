import functools
import inspect

from casatasks import casalog


def log_origin_setter(func):
    """
    This is a decorator function for a task calls other tasks.

    If the task is decorated it was called from a task and the caller set a log origin by casalog.origin(),
    then the decorator sets a caller's origin to the logger.
    Otherwise, it reads the function name and set an origin to the logger.
    So you don't need to set origin in the task any more.

    Usage:

    @log_origin_setter
    def sometask(..)
        pass

    def othertask(..)
        sometask(*args, **kwargs)  # logged "othertask::..."
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        __set_origin(inspect.stack(), casalog.getOrigin(), func.__name__)

        retval = func(*args, **kwargs)

        return retval

    return wrapper


def __set_origin(callstack, origin, funcname):
    for frame_info in callstack:
        if frame_info.function == origin:
            casalog.origin(origin)
            return
    casalog.origin(funcname)
