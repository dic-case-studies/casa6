import datetime
import inspect
import os
import shutil
from types import CodeType

from casatasks import casalog
from casatools import quanta

from . import sdutil

"""
The following code is based on the mstransform code, with
task name and some task parameters modified.
To minimise code modification, the parameters used by
mstransform but not by nrobeamaverage are kept as much as
possible, and the default values for mstransform are given
to them.
(CAS-12475, 2019/6/7 WK)
"""


@sdutil.sdtask_decorator
def nrobeamaverage(
        infile,
        datacolumn,
        field,
        spw,
        timerange,
        scan,
        beam,
        timebin,
        outfile):

    try:
        # set temporary data name
        tmpfile = 'tmp-nrobeamaverage-' + os.path.basename(infile.rstrip('/')) + '-' + \
                  "{0:%Y%m%d%H%M%S.%f}".format(datetime.datetime.now()) + '.ms'

        caller: CodeType = inspect.currentframe().f_code

        # data selection
        sdutil.do_mst(
            infile=infile,
            datacolumn=datacolumn,
            field=field,
            spw=spw,
            timerange=timerange,
            scan=scan,
            antenna='',
            timebin='0s',
            timespan='scan',
            outfile=tmpfile,
            intent='',
            caller=caller,
            ext_config={})

        # open tmpfile and rewrite antenna column of the ON spectra using beam
        idx_on = None
        with sdutil.table_manager(os.path.join(tmpfile, 'STATE')) as tb:
            ocol = tb.getcol('OBS_MODE')
            for i in range(len(ocol)):
                if ocol[i] == 'OBSERVE_TARGET#ON_SOURCE':
                    idx_on = i
                    break
        if idx_on is None:
            raise Exception('ON_SOURCE data not found.')

        with sdutil.table_manager(os.path.join(tmpfile, 'ANTENNA')) as tb:
            num_beams = len(tb.getcol('NAME'))
            _beam, min_beamid = get_beamid(beam, num_beams)

        with sdutil.table_manager(tmpfile, nomodify=False) as tb:
            acol = tb.getcol('ANTENNA1')
            scol = tb.getcol('STATE_ID')
            for i in range(len(acol)):
                if (acol[i] in _beam) and (scol[i] == idx_on):
                    acol[i] = min_beamid
            tb.putcol('ANTENNA1', acol)
            tb.putcol('ANTENNA2', acol)

        qa = quanta()
        tbin = qa.convert(qa.quantity(timebin), 's')['value']
        if tbin < 0:
            raise Exception("Parameter timebin must be > '0s' to do time averaging")
        do_timeaverage = (tbin > 0)

        ext_config = {'do_timeaverage': do_timeaverage}

        # time averaging
        sdutil.do_mst(
            infile=tmpfile,
            datacolumn=datacolumn,
            field='',
            spw='',
            timerange='',
            scan='',
            antenna='',
            timebin=timebin,
            timespan="scan",
            outfile=outfile,
            intent='',
            caller=caller,
            ext_config=ext_config)

        # History
        sdutil.add_history(caller, casalog, outfile)

    finally:
        # delete tmpfile
        if os.path.isdir(tmpfile):
            shutil.rmtree(tmpfile)


def get_beamid(beam, num_beams):
    _beam = beam
    try:
        if isinstance(_beam, str):
            _beam = _beam.strip().split(",")
        elif not isinstance(_beam, list):
            raise ValueError('the parameter beam must be list or string.')

        # the default case (beam='')
        if (len(_beam) == 1) and (_beam[0] == ''):
            _beam = []
            for i in range(num_beams):
                _beam.append(i)
        else:
            for i in range(len(_beam)):
                _beam[i] = int(_beam[i])
    except Exception as e:
        casalog.post("Error \'%s\' input beam ID is invalid" % (e))

    min_beamid = _beam[0]
    for i in range(len(_beam)):
        if _beam[i] < min_beamid:
            min_beamid = _beam[i]

    return _beam, min_beamid
