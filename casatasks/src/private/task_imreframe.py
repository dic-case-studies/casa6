from __future__ import absolute_import
import os
import shutil

# get is_CASA6 and is_python3
from casatasks.private.casa_transition import *
if is_CASA6:
    from casatools import image, measures, quanta
    from casatasks import casalog
    # this is a local tool
    qa = quanta()
else:
    from taskinit import *
    # uses qa, not a local tool

if is_python3:
    str_lower = str.lower
else:
    import string
    str_lower = string.lower
    
def imreframe(imagename=None, output=None, outframe=None, epoch=None, restfreq=None):
    try:
        casalog.origin('imreframe')
        if(output==imagename):
            output=''
        needregrid=False
        outframe=str_lower(outframe)
        if(((outframe == 'topo') or (outframe=='geo')) and (epoch != '')):
            needregrid=True
        if is_CASA6:
            myia = image()
            me = measures()
        else:
            myia,me=gentools(['ia', 'me'])
        myia.open(imagename)
        c=myia.coordsys()
        me.doframe(me.observatory(c.telescope()))
        me.doframe(c.referencevalue('m', 'direction')['measure']['direction'])
        me.doframe(c.epoch())
        myret = c.findcoordinate('spectral')
        hasspec = myret['return']
        pixax = myret['pixel']
        worldax = myret['world']
        if(not hasspec):
            raise RuntimeError('Could not find spectral axis')
        if(outframe != ''):
            c.setconversiontype(spectral=outframe)
        if(restfreq != ''):
            c.setrestfrequency(value=qa.quantity(restfreq, 'Hz'))
        if(epoch != ''):
            c.setepoch(me.epoch('utc', epoch))
            me.doframe(me.epoch('utc', epoch))
        if(not needregrid):
            if(output != ''):
                myia.fromimage(outfile=output, infile=imagename, overwrite=True)
                myia.close()
                myia.open(output)
            myia.setcoordsys(c.torecord())
            myia.done()
        else:
            c.setreferencecode(outframe, 'spectral', True)
            reffreq=c.referencevalue('m', 'spectral')['measure']['spectral']['frequency']
            newreffreq=me.measure(reffreq, outframe)
            c.setreferencevalue(qa.tos(newreffreq['m0']), 'spectral')
            outname='_temp_regrid_image' if(output=='') else output
            shp=myia.shape()             
            ib=myia.regrid(outfile=outname, shape=shp, csys=c.torecord(), 
                           axes=pixax, overwrite=True, asvelocity=False)
            ib.setcoordsys(c.torecord())
            if(output==''):
                myia.done()
                ib.rename(name=imagename, overwrite=True)
            myia.done()
            ib.done()
        
    finally:
        if('myia' in locals()):
            myia.close()
        if('ib' in locals()):
            ib.close()
