from __future__ import absolute_import
import os

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import componentlist, imager, measures, quanta, ms, table
    from casatasks import casalog
    from .setjy_helper import testerrs
    from . import solar_system_setjy as SSSetjy

    _qa = quanta( )
else:
    from taskinit import casalog, cltool, imtool, metool, qa
    from plotcomp import plotcomp
    from predictcomp_helper import *
    import pylab as pl

    componentlist = cltool
    imager = imtool
    measures = metool

    # not a local tool
    _qa = qa

def predictcomp(objname=None, standard=None, epoch=None,
                minfreq=None, maxfreq=None, nfreqs=None, prefix=None,
                antennalist=None, showplot=None, savefig=None, symb=None,
                include0amp=None, include0bl=None, blunit=None, bl0flux=None):
    """
    Writes a component list named clist to disk and returns a dict of
    {'clist': clist,
     'objname': objname,
     'angdiam': angular diameter in radians (if used in clist),
     'standard': standard,
     'epoch': epoch,
     'freqs': numpy.array of frequencies, in GHz,
     'uvrange': numpy.array of baseline lengths, in m,
     'amps':  numpy.array of predicted visibility amplitudes, in Jy,
     'savedfig': False or, if made, the filename of a plot.}
    or False on error.

    objname: An object supported by standard.
    standard: A standard for calculating flux densities, as in setjy.
              Default: 'Butler-JPL-Horizons 2010'
    epoch: The epoch to use for the calculations.   Irrelevant for
           extrasolar standards.
    minfreq: The minimum frequency to use.
             Example: '342.0GHz'
    maxfreq: The maximum frequency to use.
             Default: minfreq
             Example: '346.0GHz'
             Example: '', anything <= 0, or None: use minfreq.
    nfreqs:  The number of frequencies to use.
             Default: 1 if minfreq == maxfreq,
                      2 otherwise.
    prefix: The component list will be saved to
              prefix + '<objname>_spw0_<minfreq><epoch>.cl'
            Default: ''
    antennalist: An array configuration file as used by simdata.
                 If given, a plot of S vs. |u| will be made.
                 Default: '' (None, just make clist.)
    showplot: Whether or not to show the plot on screen.
              Subparameter of antennalist.
              Default: Necessarily False if antennalist is not specified.
                       True otherwise.
    savefig: Filename for saving a plot of S vs. |u|.
             Subparameter of antennalist.
             Default: False (necessarily if antennalist is not specified)
             Examples: True (save to prefix + '.png')
                       'myplot.png' (save to myplot.png) 
    symb: One of matplotlib's codes for plot symbols: .:,o^v<>s+xDd234hH|_
          default: '.'
    include0amp: Force the lower limit of the amplitude axis to 0.
                 Default: False
    include0bl: Force the lower limit of the baseline length axis to 0.
    blunit: Unit of the baseline length 
    bl0flux: show zero baseline flux
    """
    retval = False
    try:
        casalog.origin('predictcomp')
        # some parameter minimally required
        if objname=='':
          raise Exception("Error, objname is undefined")
        if minfreq=='':
          raise Exception("Error, minfreq is undefined")
        minfreqq = _qa.quantity(minfreq)
        minfreqHz = _qa.convert(minfreqq, 'Hz')['value']
        try:
            maxfreqq = _qa.quantity(maxfreq)
        except Exception(instance):
            maxfreqq = minfreqq
        frequnit = maxfreqq['unit']
        maxfreqHz = _qa.convert(maxfreqq, 'Hz')['value']
        if maxfreqHz <= 0.0:
            maxfreqq = minfreqq
            maxfreqHz = minfreqHz
        if minfreqHz != maxfreqHz:
            if nfreqs < 2:
                nfreqs = 2
        else:
            nfreqs = 1
        freqs = numpy.linspace(minfreqHz, maxfreqHz, nfreqs)
        
        myme = measures( )
        mepoch = myme.epoch('UTC', epoch)
        #if not prefix:
            ## meanfreq = {'value': 0.5 * (minfreqHz + maxfreqHz),
            ##             'unit': frequnit}
            ## prefix = "%s%s_%.7g" % (objname, epoch.replace('/', '-'),
            ##                         minfreqq['value'])
            ## if minfreqHz != maxfreqHz:
            ##     prefix += "to" + maxfreq
            ## else:
            ##     prefix += minfreqq['unit']
            ## prefix += "_"
        #    prefix = ''
        
        #
        if not prefix:
          if not os.access("./",os.W_OK):
            casalog.post("No write access in the current directory, trying to write cl to /tmp...","WARN")
            prefix="/tmp/"
            if not os.access(prefix, os.W_OK):
              casalog.post("No write access to /tmp to write cl file", "SEVERE")
              return False
        else:
          prefixdir=os.path.dirname(prefix)
          if prefixdir=='/' and len(prefix)>1: 
             prefix = prefix+'/'
             prefixdir = os.path.dirname(prefix)
          if not os.path.exists(prefixdir):
            prefixdirs = prefixdir.split('/')
            if prefixdirs[0]=="" and len(prefixdirs)>1:
              rootdir = "/" + prefixdirs[1]
            else:
              rootdir = "./"
            if os.access(rootdir,os.W_OK):
              if prefixdir!='':
                os.makedirs(prefixdir) 
            else:
              casalog.post("No write access to "+rootdir+" to write cl file", "SEVERE")
              return False
       
        # Get clist
        myim = imager( )
        if hasattr(myim, 'predictcomp'):
            casalog.post('local im instance created', 'DEBUG1')
        else:
            casalog.post('Error creating a local im instance.', 'SEVERE')
            return False
        #print "FREQS=",freqs
         # output CL file name is fixed : prefix+"spw0_"+minfreq+mepoch.cl
        minfreqGHz = _qa.convert(_qa.quantity(minfreq), 'GHz')['value']
        decimalfreq = minfreqGHz - int(minfreqGHz)
        decimalepoch =  mepoch['m0']['value'] - int(mepoch['m0']['value'])
        if decimalfreq == 0.0:
            minfreqGHzStr = str(int(minfreqGHz))+'GHz'
        else :
            minfreqGHzStr = str(minfreqGHz)+'GHz'

        if decimalepoch == 0.0:
            epochStr = str(int(mepoch['m0']['value']))+'d'
        else:
            epochStr=str(mepoch['m0']['value'])+'d'
        outfilename = "spw0_"+objname+"_"+minfreqGHzStr+epochStr+'.cl'
        outfilename = prefix+outfilename
        if (os.path.exists(outfilename) and os.path.isdir(outfilename)) :

          shutil.rmtree(outfilename)
          casalog.post("Removing the existing componentlist, "+outfilename)
 
        if standard=='Butler-JPL-Horizons 2012':
            clist = _predictSolarObjectCompList(objname, mepoch, freqs.tolist(), prefix)
        else:
            clist = myim.predictcomp(objname, standard, mepoch, freqs.tolist(), prefix)
        #print "created componentlist =",clist
        if os.path.isdir(clist):
            # The spw0 is useless here, but it is added by FluxStandard for the sake of setjy.
            casalog.post('The component list was saved to ' + clist)
            
            retval = {'clist': clist,
                      'objname': objname,
                      'standard': standard,
                      'epoch': mepoch,
                      'freqs (GHz)': 1.0e-9 * freqs,
                      'antennalist': antennalist}
            mycl = componentlist( )
            mycl.open(clist)
            comp = mycl.getcomponent(0)
            zeroblf=comp['flux']['value']
            if standard=='Butler-JPL-Horizons 2012':
              f0=comp['spectrum']['frequency']['m0']['value']
            else:
              f0=retval['freqs (GHz)'][0]
            casalog.post("Zero baseline flux %s @ %sGHz " % (zeroblf, f0),'INFO')
            mycl.close(False)               # False prevents the stupid warning.
            for k in ('shape', 'spectrum'):
                retval[k] = comp[k]
            if antennalist:
                retval['spectrum']['bl0flux']={}
                retval['spectrum']['bl0flux']['value']=zeroblf[0]
                retval['spectrum']['bl0flux']['unit']='Jy'
                retval['savedfig'] = None
                #retval['savedfig'] = savefig
                if not bl0flux:
                  zeroblf=[0.0]
                # no GUIs in this part of CASA
                #retval.update(plotcomp(retval, showplot, wantdict=True, symb=symb,
                #                       include0amp=include0amp, include0bl=include0bl, blunit=blunit, bl0flux=zeroblf[0]))
            else:
                retval['savedfig'] = None
        else:
            casalog.post("There was an error in making the component list.",
                         'SEVERE')

    except Exception as instance:
        casalog.post(str(instance), 'SEVERE')
        raise
    return retval

#import pylab as pl
#from taskinit import *
#import shutil
#from setjy_helper import *

def _predictSolarObjectCompList(objname, epoch, freqs, prefix):
    """
    predictcomp functionality using the new models
    set flux density of a solar system object using Bryan Butler's new
    python model calculation code.
    """
    retval = True
    cleanupcomps = False # leave genenerated cl files
    nfreqs=-1

    myms = ms( )
    mytb = table( )
    mycl = componentlist( )
    myme = measures( )

    #freqinc=freqs[0]*1e-6
    if len(freqs) == 1:
      freqinc=freqs[0]*1e-4
    else:
      freqinc=abs(freqs[1]-freqs[0])*1e-2
    freqlist=[]
    for freq in freqs:
      minf = freq-freqinc
      maxf = freq+freqinc
      freqlist.append([minf,maxf]) 
    #mepoch = myme.epoch('UTC', epoch)
    if epoch['m0']['value']==0.0:
        casalog.post('Invalid epoch, '+str(epoch['m0']['value'])+str(epoch['m0']['unit']),'SEVERE');
        raise Exception("Error")
    epochv = epoch['m0']['value'] 


    # turn user input epoch to mjd

    #print "sending objname=",objname, " epochv=",epochv, " freqlist=",freqlist
    observatory='ALMA'
    ss_setjy=SSSetjy.solar_system_setjy()
    (errcodes, fluxes, fluxerrs, sizes, dirs)=\
       ss_setjy.solar_system_fd(source_name=objname, MJDs=[epochv], frequencies=freqlist, observatory=observatory, casalog=casalog)
  
    #print "fluxes from ss_setjy=",fluxes
    #if errcodes[0][0] > 0:
    #    raise ValueError("cannot determined flux")
    
    reterr=testerrs(errcodes[0],objname) 
    if reterr == 2: 
        #raise Exception, "Error" 
        casalog.post("Flux densities cannot be determined","SEVERE")
        raise Exception

    dirstring = [dirs[0]['refer'],_qa.tos(dirs[0]['m0']),_qa.tos(dirs[0]['m1'])]
    # setup componentlists
    # need to set per dir
    # if scalebychan=F, len(freqs) corresponds to nspw selected

    #clpath='/tmp/'
    clpath='./'
    # use the first input frequency
    freqlabel = '%.3fGHz' % (freqs[0]/1.e9)
    tmlabel = '%.1fd' % epoch['m0']['value']
    clabel = objname+'_spw0_'+freqlabel+'_'+tmlabel
    clname = clabel+'.cl'
    if prefix: clname=prefix+clname

    if(os.path.exists(clname)):
        print("Removing previous cl file,", clname)
        try:
            shutil.rmtree(clname)
        except:
            print("shutil.rmtree failed")
    index= 2.0
    sptype = 'spectral index'
    #index= 0.0
    #sptype = 'constant'
    
    #print "fluxes=",fluxes
    #print "fluxerrs-=",fluxerrs
    #print "sizes=",sizes
    #print "dirs=",dirs 
    mycl.addcomponent(flux=fluxes[0][0],fluxunit='Jy', polarization="Stokes", dir=dirs[0],
         shape='disk', majoraxis=str(sizes[0][0])+'arcsec', minoraxis=str(sizes[0][1])+'arcsec',
    #     positionangle=str(sizes[0][2])+'arcsec', freq=['LSRK',str(freqs[0])+'Hz'],
         positionangle=str(sizes[0][2])+'arcsec', freq=['TOPO',str(freqs[0])+'Hz'],
         spectrumtype=sptype, index=index, label=clabel)
    # set tabular, use original input freqs
    if len(freqs)>1:
        mycl.setspectrum(which=0, type='tabular', tabularfreqs=freqs, tabularflux=fluxes[0],
                           tabularframe='TOPO')

    mycl.rename(clname)
    mycl.close(False)
    return (clname)
