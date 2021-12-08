import functools
from casatasks import casalog


def callabletask_decorator(func):
    """
    This is a decorator function for a task calls other tasks.

    if it get the parameter '__taskcaller', read it and set origin to the casalog.
    otherwise it reads the function name and set origin to the logger.
    So you don't need to set origin in the task any more.
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):

        caller = func.__name__
        if '__taskcaller' in kwargs:
            caller = kwargs['__taskcaller']
            casalog.origin(caller)
            del(kwargs['__taskcaller'])
        else:
            casalog.origin(caller)

        retval = func(*args, **kwargs)

        if caller != func.__name__:
            kwargs['__taskcaller'] = caller

        return retval

    return wrapper