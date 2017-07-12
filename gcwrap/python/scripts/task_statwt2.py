from taskinit import mstool, tbtool, casalog, write_history
import flaghelper

def statwt2(
    vis, field, spw, intent, array, observation, combine,
    timebin, chanbin, minsamp, statalg, fence, center,
    lside, zscore, maxiter, excludechans, wtrange,
    flagbackup, preview, datacolumn
):
    casalog.origin('statwt2')
    try:
        if (flagbackup):
            casalog.post('Backup original flags before applying new flags')
            flaghelper.backupFlags(aflocal=None, msfile=vis, prename='statwt2')
        myms = mstool()
        myms.open(vis, nomodify=False)
        sel = {}
        selectdata = not (type(field) == str and len(field) == 0)
        if (selectdata):
            sel['spw'] = spw 
            #sel['time'] = timerange
            sel['field'] = field
            #sel['baseline'] = antenna
            #sel['scan'] = scan
            sel['scanintent'] = intent
            #sel['polarization'] = correlation
            #sel['uvdist'] = uvrange
            sel['observation'] = str(observation)
            sel['array'] = array
            #sel['feed'] = feed
            # Select the data. Only-parse is set to false.
            myms.msselect(sel, False)
        myms.statwt2(
            combine=combine, timebin=timebin, chanbin=chanbin,
            minsamp=minsamp, statalg=statalg, fence=fence,
            center=center, lside=lside, zscore=zscore,
            maxiter=maxiter, excludechans=excludechans,
            wtrange=wtrange, preview=preview, datacolumn=datacolumn
        ) 
        return True
    except Exception, instance:
        casalog.post( '*** Error ***'+str(instance), 'SEVERE' )
        raise
    finally:
        myms.done()
