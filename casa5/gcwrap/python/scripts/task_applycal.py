from __future__ import absolute_import
import os
import time
import numpy as np

# get is_CASA6 and is_python3
from casatasks.private.casa_transition import *
if is_CASA6:
    from .. import casalog
    from .callibrary import *
    from . import flaghelper as fh
    from .parallel.parallel_data_helper import ParallelDataHelper
    from .parallel.parallel_task_helper import ParallelTaskHelper
    from .mstools import write_history
    from casatools import ms, calibrater
else:
    from taskinit import *
    from mstools import write_history
    from callibrary import *
    import flaghelper as fh
    from parallel.parallel_data_helper import ParallelDataHelper
    from parallel.parallel_task_helper import ParallelTaskHelper

    calibrater = cbtool
    ms = mstool

def applycal(
    vis=None,
    field=None,
    spw=None,
    intent=None,
    selectdata=None,
    timerange=None,
    uvrange=None,
    antenna=None,
    scan=None,
    observation=None,
    msselect=None,
    docallib=None,
    callib=None,
    gaintable=None,
    gainfield=None,
    interp=None,
    spwmap=None,
    calwt=None,
    parang=None,
    applymode=None,
    flagbackup=None,
    ):

    # Python script
    casalog.origin('applycal')

    # Take care of the trivial parallelization
    if ParallelDataHelper.isMMSAndNotServer(vis):
        
        # Back up the flags, if requested (and if necessary)
        if flagbackup and applymode != 'calonly' and applymode != 'trial':
            fh.backupFlags(aflocal=None, msfile=vis, prename='applycal')
            flagbackup = False
        
        # To be safe convert file names to absolute paths.
        gaintable = ParallelTaskHelper.findAbsPath(gaintable)
        helper = ParallelTaskHelper('applycal', locals())
        ret = helper.go()
        if ParallelTaskHelper.getAsyncMode():
            return ret
        else:
            return

    try:
        mycb = calibrater( )
        makecorr = True
        if (type(vis) == str) & os.path.exists(vis):
            # if applymode is flagonly don't add CORRECTED_DATA column
            # otherwise add CORRECTED_DATA column
            if applymode in ['flagonly', 'flagonlystrict', 'trial']:
                makecorr = False
            mycb.open(filename=vis, compress=False, addcorr=makecorr,
                      addmodel=False)
        else:
            raise ValueError( 'Visibility data set not found - please verify the name' )

        # enforce default if unspecified
        if applymode == '':
            applymode = 'calflag'

        # Back up the flags, if requested (and if necessary)
        if flagbackup and applymode != 'calonly' and applymode \
            != 'trial':
            fh.backupFlags(aflocal=None, msfile=vis, prename='applycal')

        # Do data selection according to selectdata
        if selectdata:
            # pass all data selection parameters in as specified
            mycb.selectvis(
                time=timerange,
                spw=spw,
                scan=scan,
                field=field,
                intent=intent,
                observation=str(observation),
                baseline=antenna,
                uvrange=uvrange,
                chanmode='none',
                msselect=msselect,
                )
        else:
            # selectdata=F, so time,scan,baseline,uvrange,msselect=''
            # using spw and field specifications only
            mycb.selectvis(
                time='',
                spw=spw,
                scan='',
                field=field,
                intent=intent,
                observation='',
                baseline='',
                uvrange='',
                chanmode='none',
                msselect='',
                )

        # Arrange applies....

        if docallib:
            # by cal library from file
            # parsing using c++ parser
            thiscallib=mycb.parsecallibfile(callib)
            mycb.setcallib(thiscallib)

        else:

            # by traditional parameters

            ngaintab = 0
            if gaintable != ['']:
                ngaintab = len(gaintable)

            ncalwt = len(calwt)
            if ncalwt == 1:
                calwt = [calwt[0] for i in range(ngaintab)]

            ngainfld = len(gainfield)
            nspwmap = len(spwmap)
            ninterp = len(interp)

            # handle list of list issues with spwmap
            if nspwmap > 0:
                if type(spwmap[0]) != list:
                    # first element not a list, only one spwmap specified
                    # make it a list of list
                    spwmap = [spwmap]
                    nspwmap = 1

            for igt in range(ngaintab):
                if gaintable[igt] != '':

                    # field selection is null unless specified
                    thisgainfield = ''
                    if igt < ngainfld:
                        thisgainfield = gainfield[igt]

                    # spwmap is null unless specifed
                    thisspwmap = [-1]
                    if igt < nspwmap:
                        thisspwmap = spwmap[igt]

                    # interp is 'linear' unless specified
                    thisinterp = 'linear'
                    if igt < ninterp:
                        if interp[igt] == '':
                            interp[igt] = thisinterp
                        thisinterp = interp[igt]

                    mycb.setapply(
                        t=0.0,
                        table=gaintable[igt],
                        field=thisgainfield,
                        calwt=calwt[igt],
                        spwmap=thisspwmap,
                        interp=thisinterp,
                        )

        # ...and now the specialized terms

        # Apply parallactic angle, if requested
        if parang:
            mycb.setapply(type='P')

        mycb.correct(applymode)

        # report what the flags did
        reportflags(mycb.activityrec())


        # write history
        try:
            param_names = \
                          applycal.__code__.co_varnames[:applycal.__code__.co_argcount]
            if is_python3:
                vars = locals( )
                param_vals = [vars[p] for p in param_names]
            else:
                param_vals = [eval(p) for p in param_names]
            write_history(
                ms(),
                vis,
                'applycal',
                param_names,
                param_vals,
                casalog,
                )
        except Exception as instance:
            casalog.post("*** Error \'%s\' updating HISTORY"
                         % instance, 'WARN')
    finally:
        mycb.close()

def reportflags(rec):
    try:
        if 'origin' in rec and rec['origin'] \
            == 'Calibrater::correct' and 'VisEquation' in rec:
            casalog.post('Calibration apply flagging statistics (among calibrateable spws):'
                         )
            VE = rec['VisEquation']
            nterm = len(VE)
            if nterm > 0:
                nVisTotal=VE['nVisTotal']

                casalog.post('  Total visibilities selected for correction (ncorr x nchan x nrow summed over spws) = '
                              + str(nVisTotal))
                casalog.post('  Flags:')
                partlog=False
                for iterm in range(nterm-1):    # one of the keys is nVisTotal; the rest are caltable indices
                    VEi = VE['*' + str(iterm + 1)]

                    nVisThis=VEi['ndata']
                    partial='  '
                    if nVisThis<nVisTotal:
                        partlog=True
                        partial='**'
                    flstr = '   ' + VEi['type']
                    flstr += ': '
                    flstr += 'In: ' + str(VEi['nflagIn'])
                    flstr += ' / '+str(nVisThis)+partial
                    flstr += ' (' + str(100. * VEi['nflagIn']/nVisThis) + '%) --> '
                    flstr += 'Out: ' + str(VEi['nflagOut'])
                    flstr += ' / '+str(nVisThis)+partial
                    flstr += ' (' + str(100. * VEi['nflagOut']/nVisThis) + '%)'

                    if 'table' in VEi:
                        flstr += ' (' + VEi['table'] + ')'
                    casalog.post(flstr)
                if partlog:
                    casalog.post('     ** = Denotes caltable that only corrected a subset of total selected visibilities')

    except Exception as instance:
        # complain mildly, but don't alarm
        casalog.post('Error formatting some or all of the applycal flagging log info: '
                      + str(instance), 'SEVERE')


