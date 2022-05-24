import inspect
from types import CodeType

from casatasks import casalog
from casatools import quanta

from . import sdutil

qa = quanta()

"""
The following code is based on the mstransform code, with
task name and some task parameters modified.
To minimise code modification, the parameters used by
mstransform but not by sdpolaverage are kept and the
default values for mstransform are given to them.
(CAS-12083, 2019/1/22 WK)
"""


@sdutil.sdtask_decorator
def sdpolaverage(
        infile,
        datacolumn,
        antenna,
        field,
        spw,
        timerange,
        scan,
        intent,
        polaverage,
        outfile):

    # followings are parameters of mstransform but not used by sdpolaverage.
    # just putting default values
    do_timeaverage = False
    timebin = "0s"
    timespan = ""

    # Only parse timeaverage parameters when timebin > 0s
    if do_timeaverage:
        tb = qa.convert(qa.quantity(timebin), 's')['value']
        if not tb > 0:
            raise ValueError("Parameter timebin must be > '0s' to do time averaging")

    # extra parameter for do_mst
    ext_config = {"do_timeaverage": False,
                  "keepflags": True,
                  "do_check_tileshape": True,
                  "polaverage": polaverage,
                  "parse_chanaverage": False}

    sdpolaverage: CodeType = inspect.currentframe().f_code

    sdutil.do_mst(
        infile,
        datacolumn,
        field,
        spw,
        timerange,
        scan,
        antenna,
        timebin,
        timespan,
        outfile,
        intent,
        sdpolaverage,
        ext_config)

    sdutil.add_history(sdpolaverage, casalog, outfile)
