import functools
import inspect

from casatasks import casalog


def log_origin_setter(func):
    """
    This is a decorator function for a casatask calling other casatasks.

    This decorator is intended to be used on a "super" casatask calling "sub" casatasks.
    Its effect is to set the origin of messages logged by the "sub" tasks called by the "super" task to "super"
    For example:
    super_casatask::casa "Msg from super task"
    super_casatask::casa "Msg from sub-task1 called by super task"
    super_casatask::casa "Msg from sub-task1 called by super task"

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
