import functools
from casatasks import casalog


def log_origin_setter(func):
    """
    This is a decorator function for a task calls other tasks.

    if it get the parameter '__log_origin', read it and set origin to the casalog.
    otherwise it reads the function name and set origin to the logger.
    So you don't need to set origin in the task any more.

    Usage:

    @log_origin_setter
    def sometask(..)
        pass

    def othertask(..)
        kwargs['__log_origin'] = 'othertask'
        sometask(*args, **kwargs)  # logged "othertask::..."
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        caller = kwargs.pop('__log_origin', func.__name__)
        casalog.origin(caller)

        retval = func(*args, **kwargs)

        if caller != func.__name__:
            kwargs['__log_origin'] = caller

        return retval

    return wrapper
