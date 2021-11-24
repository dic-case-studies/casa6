import functools
import re

from casatasks import casalog
import traceback


def callabletask_decorator(func):
    """
    This is a decorator function for tasks call from other tasks.
    Currently the decorator does:

       1) set origin to the logger
       2) handle exception

    So, you don't need to set origin in the task any more.
    Also, you don't need to write anything about error
    handling in the task. If you have something to do
    at the end of the task execution, those should be
    written in the destructor of worker class, not in
    the 'finally' block.
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

        retval = None
        # Any errors are handled outside the task.
        # however, the implementation below is effectively
        # equivalent to handling it inside the task.
        try:
            # execute task
            retval = func(*args, **kwargs)
        except Exception as e:
            traceback_info = __format_trace(traceback.format_exc())
            casalog.post(traceback_info, 'SEVERE')
            casalog.post(str(e), 'ERROR')
            raise
        finally:
            if caller != func.__name__:
                kwargs['__taskcaller'] = caller
        return retval
    return wrapper


def __format_trace(s):
    wexists = True
    regex = '.*callabletask\.py.*in wrapper.*'
    retval = s
    while wexists:
        ss = retval.split('\n')
        wexists = False
        for i in range(len(ss)):
            if re.match(regex, ss[i]):
                ss = ss[:i] + ss[i + 2:]
                wexists = True
                break
        retval = '\n'.join(ss)
    return retval
