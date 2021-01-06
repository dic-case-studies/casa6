from __future__ import absolute_import
import os
import numpy as np

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from .callibrary import *
    from casatools import calibrater
    from casatasks import casalog
else:
    from callibrary import *
    from taskinit import *

    calibrater = cbtool

def fringefit(vis=None,caltable=None,
              field=None,spw=None,intent=None,
              selectdata=None,timerange=None,antenna=None,scan=None,
              observation=None, msselect=None,
              solint=None,combine=None,refant=None,
              minsnr=None,zerorates=None,globalsolve=None,niter=None,
              delaywindow=None,ratewindow=None,append=None,
              corrdepflags=None,
              docallib=None,callib=None,gaintable=None,gainfield=None,interp=None,spwmap=None,
              paramactive=None,
              parang=None):

    #Python script
    casalog.origin('fringefit')

    try: 
        mycb = calibrater()

        if ((type(vis)==str) & (os.path.exists(vis))):
            mycb.open(filename=vis,compress=False,addcorr=False,addmodel=False)
        else:
            raise ValueError('Visibility data set not found - please verify the name')

        # Do data selection according to selectdata
        if (selectdata):
            # pass all data selection parameters in as specified
            mycb.selectvis(time=timerange,spw=spw, scan=scan, field=field,
                           intent=intent, observation=str(observation),
                           baseline=antenna,chanmode='none',
                           msselect=msselect)
        else:
            # selectdata=F, so time,scan,baseline,msselect=''
            # using spw and field specifications only
            mycb.selectvis(time='',spw=spw,scan='',field=field,intent=intent,
                           observation='', baseline='',
                           chanmode='none', msselect='')

        # signal use of correlation-dependent flags, if requested
        if corrdepflags:
            mycb.setcorrdepflags(True)

                        
        # Arrange applies....
            
        if docallib:
            # by cal library from file
            mycallib=callibrary()
            mycallib.read(callib)
            mycb.setcallib(mycallib.cld)

        else:
            if paramactive is None or paramactive==[]:
                paramactive=[True, True, False]
            else:
                if len(paramactive)!=3:
                    casalog.post("paramactive: " + paramactive)
                    raise ValueError( 'Error: paramactive vector must have exactly three entries' )
            # Have to solve for peculiar phase!
            paramactive.insert(0, True)

            # by traditional parameters

            ngaintab = 0;
            if (gaintable!=['']):
                ngaintab=len(gaintable)

            ngainfld = len(gainfield)
            nspwmap = len(spwmap)
            ninterp = len(interp)

            # handle list of list issues with spwmap
            if (nspwmap>0):
                if (type(spwmap[0])!=list):
                    # first element not a list, only one spwmap specified
                    # make it a list of list
                    spwmap=[spwmap];
                    nspwmap=1;

            for igt in range(ngaintab):
                if (gaintable[igt]!=''):

                    # field selection is null unless specified
                    thisgainfield=''
                    if (igt<ngainfld):
                        thisgainfield=gainfield[igt]

                    # spwmap is null unless specifed
                    thisspwmap=[-1]
                    if (igt<nspwmap):
                        thisspwmap=spwmap[igt];

                    # interp is 'linear' unless specified
                    thisinterp='linear'
                    if (igt<ninterp):
                        if (interp[igt]==''):
                            interp[igt]=thisinterp;
                        thisinterp=interp[igt];

                    mycb.setapply(t=0.0,table=gaintable[igt],field=thisgainfield,
                                  calwt=True,spwmap=thisspwmap,interp=thisinterp)

            if len(delaywindow) != 2:
                delaywindow = [-1e6, 1e6]
            if len(ratewindow) != 2:
                ratewindow = [-1e6, 1e6]

        # ...and now the specialized terms
        # (BTW, interp irrelevant for these, since they are evaluated)
                
        # Apply parallactic angle, if requested
        if parang: mycb.setapply(type='P')

        # Set up for solving; only support one gaintype
        mycb.setsolve(type="FRINGE",t=solint,refant=refant,preavg=0.01,
                      minsnr=minsnr,combine=combine,
                      zerorates=zerorates,
                      globalsolve=globalsolve,
                      niter=niter,
                      delaywindow=delaywindow,
                      ratewindow=ratewindow,
                      paramactive=paramactive,
                      table=caltable,append=append)
        mycb.solve()

        reportsolvestats(mycb.activityrec());

    finally:
        mycb.close()

def reportsolvestats(rec):
    if (list(rec.keys()).count('origin')==1 and
        rec['origin']=='Calibrater::genericGatherAndSolve'):
        casalog.post("Calibration solve statistics per spw:  (expected/attempted/succeeded):")
        nexp=rec['nExpected']
        natt=rec['nAttempt']
        nsuc=rec['nSucceed']
        for ispw in range(len(nexp)):
            solstatstr="  Spw "+str(ispw)+": "
            solstatstr+=str(nexp[ispw])+"/"
            solstatstr+=str(natt[ispw])+"/"
            solstatstr+=str(nsuc[ispw])
            casalog.post(solstatstr)
