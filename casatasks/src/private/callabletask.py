import functools
from casatasks import casalog


def callabletask_decorator(func):
    """
    This is a decorator function for a task calls other tasks.

    if it get the parameter '__taskcaller', read it and set origin to the casalog.
    otherwise it reads the function name and set origin to the logger.
    So you don't need to set origin in the task any more.

    Usage:

    @callabletask_decorator
    def sometask(..)
        pass

    def othertask(..)
        kwargs['__taskcaller'] = 'othertask'
        sometask(*args, **kwargs)  # logged "othertask::..."
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):

        caller = kwargs.pop('__taskcaller', func.__name__)
        casalog.origin(caller)

        retval = func(*args, **kwargs)

        if caller != func.__name__:
            kwargs['__taskcaller'] = caller

        return retval

    return wrapper