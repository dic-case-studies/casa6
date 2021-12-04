from __future__ import absolute_import
import os

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import imager, ms, image, quanta
    from casatasks import casalog

    _im = imager( )
    _ms = ms( )
    _ia = image( )
else:
    from taskinit import *
    _im,_ms,_ia=gentools(['im','ms','ia'])

    quanta = qatool

def ft(vis=None,field=None,spw=None,model=None,nterms=None,reffreq=None,complist=None,incremental=None, usescratch=None):
       """ Insert a source model into the MODEL_DATA column of a visibility set:

       A source model (souce.model image) or components list is converted into a
       model visibility that is inserted into the MODEL_DATA column.  This is
       needed to use resolved source in gaincal and in fluxscale.  (Setjy will
       automatically make this ft step.)

       The sources currently available are 3C48, 3C138, 3C147, 3C286
       at 1.4, 5.0, 8.4, 15, 22, 43 GHz.  Their location is site
       dependent.  In Charlottesville and at the AOC, the models are
       in /usr/lib/casapy/data/nrao/VLA/CalModels.

       Keyword arguments:
       vis -- Name of input visibility file
               default: none; example: vis='ngc5921.ms'
       field -- Field name list
               default: '' ==> all
               NOTE: each source must be specified in a multi-source vis.
               field = '1328+307'  specifies source '1328+307'
               field = '4' specified field with index 4
       spw -- Spw selection
               default: spw = '' (all spw)
       model -- Name of input model image
               default: None;
               example: model='/usr/lib/casapy/data/nrao/VLA/CalModels/3C286_X.im'
               Note: The model visibilities are scaled from the model frequency
                     to the observed frequency of the data.
       nterms -- Number of terms used to model the sky frequency dependence
                 default: 1
                 example : nterms=3  represents a 2nd order Taylor-polynomial in frequency
                           and is to be used along with 3 model-image names. 
                           model=['xxx.image.tt0','xxx.image.tt1', 'xxx.image.tt2']
          reffreq -- Reference-frequency about which this Taylor-expansion is defined.
       complist -- Name of component list
               default: None; ; example: complist='test.cl'
               components tool not yet available
       incremental -- Add model visibility to the existing MODEL_DATA visibilties
               default: False; example: incremental=True

       """
       casalog.origin('ft')

       #Python script

       # Check if datafile exists and open it
       if ((type(vis)==str) & (os.path.exists(vis))):
               _im.open(vis, usescratch=usescratch)
       else:
               raise ValueError('Visibility data set not found - please verify the name')

       # Select data
       _im.selectvis(field=field,spw=spw)

       # Define image co-ordinates (all defaults)
       #_im.defineimage()

       # Check 'model'. The 'xml' allows a variant => do the checking here.
       if( (not type(model)==str) and (not (type(model)==list) ) ) :
               raise ValueError('The model image must be a string or a list of strings (or \'\' or [])')

       # If model is a single string, make it a list
       if( type(model)==str ):
               model = [model];

       # Check that either a model or a complist has been given.
       if( (model==[] or model==['']) and complist=='' ):
               raise ValueError('Please specify a model image or component list to ft')

       #model is a list now. Check that all elements are strings. If so, check file existence too.
       if( type(model)==list ):
               for onemodel in model:
                      if(not type(onemodel)==str):
                            raise ValueError('Model image names must be strings')
                      if( (not onemodel=='') and (not os.path.exists(onemodel)) ):
                            raise ValueError('Model image '+onemodel+' cannot be found')

       # Check complist : one string : name of complist file. Check existance on disk.
       if( (not complist=='') and (not os.path.exists(complist)) ):
               raise ValueError('Componentlist '+complist+' cannot be found')


       # If nterms>1, then check that len(model)=nterms [ no multifield for now ]
       # Call _im.settaylorterms()
       #
       if (nterms > 1) :
               if(type(model)==str or (not (type(model)==list and len(model)==nterms)) ):
                       raise ValueError('For nterms>1, please provide a list of nterms model-image names')
               # parse the reference-frequency field.
               qat=quanta();
               try:
                  rff=qat.canonical(reffreq);
               except Exception:
                  msg = '*** Error *** In conversion of reffreq=\'',reffreq,'\' to a numerical value'
                  raise RuntimeError(msg)

               reffreqVal=rff['value'];  # This is the frequency in Hz
               if(reffreqVal==0.0):   # if unspecified, set the default from the model image
                       _ia.open(model[0]);
                       icsys = _ia.coordsys();
                       _ia.close();
                       reffreqVal=icsys.referencevalue(type='spectral')['numeric'][0];
                       casalog.post('Using reference frequency from model image : '+str(reffreqVal)+' Hz');
               else:
                       casalog.post('Using reference frequency : '+str(reffreqVal)+' Hz');
               # set nterms and ref-freq
               _im.settaylorterms(ntaylorterms=nterms,reffreq=reffreqVal)

       # Just checking...
       if (nterms < 1) :
               raise ValueError('nterms must be greater than or equal to 1')


       # Do the forward transform and close.
       _im.ft(model=model,complist=complist,incremental=incremental)
       _im.close()


       #write history
       _ms.open(vis,nomodify=False)
       _ms.writehistory(message='taskname = ft',origin='ft')
       _ms.writehistory(message='vis         = "'+str(vis)+'"',origin='ft')
       _ms.writehistory(message='field       = "'+str(field)+'"',origin='ft')
       _ms.writehistory(message='spw         = "'+str(spw)+'"',origin='ft')
       _ms.writehistory(message='model       = "'+str(model)+'"',origin='ft')
       _ms.writehistory(message='complist    = "'+str(complist)+'"',origin='ft')
       _ms.writehistory(message='incremental = "'+str(incremental)+'"',origin='ft')
       _ms.close()
