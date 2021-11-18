from __future__ import absolute_import
import os

from casatasks.private.casa_transition import *
if is_CASA6:
    from casatasks import casalog
    from casatools import table as tbtool
    from casatools import ms as mstool
    from .mstools import write_history
else:
    from taskinit import *
    from mstools import write_history

def conjugatevis(vis,spwlist=[],outputvis="",overwrite=False):
    """:
    Change the sign of the phases in all visibility columns

    Keyword arguments:
    vis -- Name of input visibility file
        default: none; example='3C273XC1.ms'
    spwlist -- Select spectral window
        default: [] all spws will be conjugated; example: spw=[1,2]
    outputvis -- name of output visibility file
        default: 'conjugated_'+vis; example= 'conjugated.ms'
    overwrite -- Overwrite the outputvis if it exists
        default=False; example: overwrite=True

    """

    #Python script

    _tb = tbtool()
        
    try:
        casalog.origin('conjugatevis')
        myddlist = []
        _tb.open(vis+'/SPECTRAL_WINDOW')
        maxspw = _tb.nrows()-1
        _tb.close()
        if (type(spwlist)==type(1)):
            spwlist = [spwlist]
        elif(spwlist==None or spwlist==[] or spwlist==""):
            spwlist = []
            casalog.post("Will conjugate visibilities for all spws.", 'INFO')
        if not spwlist==[]:
            try:
                _tb.open(vis+'/DATA_DESCRIPTION')
                for k in spwlist:
                    if (k<-1 or k>maxspw):
                        raise RuntimeError("Error: max valid spw id is "+str(maxspw))
                    else:
                        for j in range(0,_tb.nrows()):
                            if(_tb.getcell("SPECTRAL_WINDOW_ID",j)==k and not (j in myddlist)):
                                myddlist = myddlist + [j]
                #end for k
                _tb.close()
                casalog.post('DATA_DESC_IDs to process: '+str(myddlist), 'INFO')
            except Exception as exc:
                raise RuntimeError('Error reading DATA_DESCRIPTION table: {}'.format(exc))
        #endif
        outname = 'conjugated_'+vis
        if not (outputvis==""):
            if((type(outputvis)!=str) or (len(outputvis.split()) < 1)):
                raise ValueError('parameter outputvis is invalid')
            outname = outputvis
        if not overwrite and os.path.exists(outname):
            raise RuntimeError('outputvis '+outname+' exists and you did not permit overwrite')
        os.system('rm -rf '+outname)
        os.system('cp -R '+vis+' '+outname)
        _tb.open(outname, nomodify=False)
        if _tb.iswritable():
            if(spwlist==[]):
                for colname in [ 'DATA', 'CORRECTED_DATA', 'FLOAT_DATA' ]:
                    if colname in _tb.colnames():
                        casalog.post('Conjugating '+str(colname), 'INFO')
                        for i in range(0,_tb.nrows()):
                            a = _tb.getcell(colname, i)
                            a = a.conjugate()
                            _tb.putcell(colname, i, a)
            else:
                for colname in [ 'DATA', 'CORRECTED_DATA', 'FLOAT_DATA' ]:
                    if colname in _tb.colnames():
                        casalog.post('Conjugating '+str(colname), 'INFO')
                        for i in range(0,_tb.nrows()):
                            if(_tb.getcell("DATA_DESC_ID",i) in myddlist):
                                a = _tb.getcell(colname, i)
                                a = a.conjugate()
                                _tb.putcell(colname, i, a)
            #endif
            _tb.flush()
            _tb.close()
            casalog.post('Created '+str(outname), 'INFO')
        else:
            _tb.close()
            casalog.post('Cannot write to output MS '+str(outname), 'WARN')

        # Write history to output MS 
        try:
            param_names = conjugatevis.__code__.co_varnames[:conjugatevis.__code__.co_argcount]
            if is_python3:
                vars = locals()
                param_vals = [vars[p] for p in param_names]
            else:
                param_vals = [eval(p) for p in param_names]
            write_history(mstool(), outname, 'conjugatevis', param_names,
                          param_vals, casalog)

        except Exception as instance:
            casalog.post("*** Error \'%s\' updating HISTORY" % (instance), 'WARN')

    finally:
        _tb.close()
