from __future__ import absolute_import
from __future__ import print_function

from scipy import fftpack
import numpy as np
import shutil
import os
import time

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import image, quanta, table, imager
    from casatasks import casalog, imsubimage, feather
else:
    from taskinit import *
    from tasks import *
    image = iatool
    imager = imtool
    quanta = qatool
    casalog = casac.logsink()

_ia = image()
_qa = quanta()

class SDINT_helper:

#    def __init__(self):
#      print 'Init Helper'

################################################
    def getFreqList(self,imname=''):

      print("imname=",imname)
      _ia.open(imname)
      csys =_ia.coordsys()
      shp = _ia.shape()
      _ia.close()

      if(csys.axiscoordinatetypes()[3] == 'Spectral'):
           restfreq = csys.referencevalue()['numeric'][3]#/1.0e+09; # convert more generally..
           freqincrement = csys.increment()['numeric'][3]# /1.0e+09;
           freqlist = [];
           for chan in range(0,shp[3]):
                 freqlist.append(restfreq + chan * freqincrement);
      elif(csys.axiscoordinatetypes()[3] == 'Tabular'):
           freqlist = (csys.torecord()['tabular2']['worldvalues']) # /1.0e+09;
      else:
           casalog.post('Unknown frequency axis. Exiting.','SEVERE');
           return False;

      return freqlist

################################################

    def copy_restoringbeam(self,fromthis='',tothis=''):
        _ib = image()
#        ia.open(fromthis);
#        ib.open(tothis)
        freqlist = self.getFreqList(fromthis)
#        print freqlist
        for i in range(len(freqlist)):
            _ia.open(fromthis);
            beam = _ia.restoringbeam(channel = i);
            _ia.close()
            #print beam
            _ia.open(tothis)
            _ia.setrestoringbeam(beam = beam, channel = i, polarization = 0);
            _ia.close()
#        ib.close();
#        ia.close();	
        
################################################

    def feather_int_sd(self,sdcube='', intcube='', jointcube='',sdgain=1.0,dishdia=100.0, usedata='sdint'): 
#, pbcube='',applypb=False, pblimit=0.2):
        """
        Run the feather task to combine the SD and INT Cubes. 
        
        There's a bug in feather for cubes. Hence, do each channel separately.
        FIX feather and then change this....

        TODO : Add the effdishdia  usage to get freq-indep feathering.

        """
         
        ### Do the feathering.
        if usedata=='sdint':
            ## Feather runs in a loop on chans internally, but there are issues with open tablecache images
            ## Also, no way to set effective dish dia separately for each channel.
            #feather(imagename = jointcube, highres = intcube, lowres = sdcube, sdfactor = sdgain, effdishdiam=-1)

            freqlist = self.getFreqList(sdcube)
            
            os.system('rm -rf '+jointcube)
            os.system('cp -r ' + intcube + ' ' + jointcube)
            
            _ib = image()
            
            _ia.open(jointcube)
            
            for i in range(len(freqlist)):	

                freqdishdia = dishdia ## * (freqlist[0] / freqlist[i]) # * 0.5
                
                os.system('rm -rf tmp_*')
                imsubimage(imagename = sdcube, outfile = 'tmp_sdplane', chans = str(i));
                imsubimage(imagename = intcube, outfile = 'tmp_intplane', chans = str(i));

                #feather(imagename = 'tmp_jointplane', highres = 'tmp_intplane', lowres = 'tmp_sdplane', sdfactor = sdgain, effdishdiam=freqdishdia)
                # feathering via toolkit
                try: 
                    imFea=imager( )
                    imFea.setvp(dovp=True)
                    imFea.setsdoptions(scale=sdgain)
                    imFea.feather(image='tmp_jointplane',highres='tmp_intplane',lowres='tmp_sdplane', effdishdiam=freqdishdia)
                    imFea.done( )
                    del imFea
                except Exception as instance:
                    print('*** Error *** %s' % instance)
                    raise 

                _ib.open('tmp_jointplane')
                pixjoint = _ib.getchunk()
                _ib.close()
                _ia.putchunk(pixjoint, blc=[0,0,0,i])
            
            _ia.close()

        if usedata=='sd':
            ## Copy sdcube to joint.
            os.system('rm -rf '+jointcube)
            os.system('cp -r ' + sdcube + ' ' + jointcube)
        if usedata=='int':
            ## Copy intcube to joint
            os.system('rm -rf '+jointcube)
            os.system('cp -r ' + intcube + ' ' + jointcube)

################################################

    def modify_with_pb(self, inpcube='', pbcube='', action='mult',pblimit=0.2, freqdep=True):
        """
        Multiply or divide by the PB

        freqdep = True :  Channel by channel
        freqdep = False : Before/After deconvolution, use a freq-independent PB from the middle of the list
        """

        freqlist = self.getFreqList(inpcube)

        _ia.open(inpcube)
        shp=_ia.shape()
        _ia.close()

        midchan = int(len(freqlist)/2)
        _ia.open(pbcube)
        pbplane = _ia.getchunk(blc=[0,0,0,midchan],trc=[shp[0],shp[1],0,midchan])
        _ia.close()


        for i in range(len(freqlist)):

            if freqdep==True:
                _ia.open(pbcube)
                pbplane = _ia.getchunk(blc=[0,0,0,i],trc=[shp[0],shp[1],0,i])
                _ia.close()

            _ia.open(inpcube)
            implane = _ia.getchunk(blc=[0,0,0,i],trc=[shp[0],shp[1],0,i])

            outplane = pbplane.copy()
            outplane.fill(0.0)

            if action=='mult':
                pbplane[pbplane<pblimit]=0.0
                outplane = implane * pbplane
            else:
                implane[pbplane<pblimit]=0.0
                pbplane[pbplane<pblimit]=1.0
                outplane = implane / pbplane

            _ia.putchunk(outplane, blc=[0,0,0,i])
            _ia.close()

            self.addmask(inpcube,pbcube,pblimit)

################################################
    def addmask(self, inpimage='',pbimage='',pblimit=0.2):
        _ia.open(inpimage)
        _ia.calcmask(mask='"'+pbimage+'"'+'>'+str(pblimit));
        _ia.close()
        

################################################
    def cube_to_taylor_sum(self, cubename='', mtname='',reffreq='1.5GHz',nterms=2,dopsf=False):
        """
        Convert Cubes (output of major cycle) to Taylor weighted averages (inputs to the minor cycle)
        Input : Cube
        Output : Set of images with suffix : .tt0, .tt1, etc...
        """

        refnu = _qa.convert( _qa.quantity(reffreq) ,'Hz' )['value']

        pix=[]

        num_terms=nterms

        if dopsf==True:
            num_terms=2*nterms-1
 
        for tt in range(0,num_terms):
            _ia.open(mtname+'.tt'+str(tt))
            pix.append( _ia.getchunk() )
            _ia.close()
            pix[tt].fill(0.0)

        _ia.open(cubename)
        shp = _ia.shape()
        _ia.close()

        freqlist = self.getFreqList(cubename)
        for i in range(len(freqlist)):
            wt = (freqlist[i] - refnu)/refnu
            _ia.open(cubename)
            implane = _ia.getchunk(blc=[0,0,0,i],trc=[shp[0],shp[1],0,i])
            _ia.close()
            for tt in range(0,num_terms):
                pix[tt] = pix[tt] + (wt**tt) * implane

#        ia.close()

        for tt in range(0,num_terms):
            _ia.open(mtname+'.tt'+str(tt))
            _ia.putchunk(pix[tt]/len(freqlist))
            _ia.close()

################################################
    def taylor_model_to_cube(self, cubename='', mtname='',reffreq='1.5GHz',nterms=2):
        """
        Convert Taylor coefficients (output of minor cycle) to cube (input to major cycle)
        Input : Set of images with suffix : .tt0, .tt1, etc...
        Output : Cube
        """

        if not os.path.exists(cubename+'.model'):
            shutil.copytree(cubename+'.psf', cubename+'.model')
            _ia.open(cubename+'.model')
            _ia.set(0.0)
            _ia.setrestoringbeam(remove=True)
            _ia.setbrightnessunit('Jy/pixel')
            _ia.close()


        refnu = _qa.convert( _qa.quantity(reffreq) ,'Hz' )['value']

        pix=[]

        for tt in range(0,nterms):
            _ia.open(mtname+'.model.tt'+str(tt))
            pix.append( _ia.getchunk() )
            _ia.close()

        _ia.open(cubename+'.model')
        shp = _ia.shape()
        _ia.close()

        implane = pix[0].copy()

        freqlist = self.getFreqList(cubename+'.psf')
        for i in range(len(freqlist)):
            wt = (freqlist[i] - refnu)/refnu
            implane.fill(0.0)
            for tt in range(0,nterms):
                implane = implane + (wt**tt) * pix[tt]
            _ia.open(cubename+'.model')
            _ia.putchunk(implane, blc=[0,0,0,i])
            _ia.close()

##################################################


    def calc_sd_residual(self,origcube='', modelcube = '', residualcube = '', psfcube=''):
##, pbcube='', applypb=False, pblimit=0.2):
        """
        Residual = Original - ( PSF * Model )
        """

#        ia_orig = iatool()
#        ia_model = iatool()
#        ia_residual = iatool()
#        ia_psf = iatool()
#
#        ia_orig.open(origcube)
#        ia_model.open(modelcube)
#        ia_residual.open(residualcube)
#        ia_psf.open(psfcube)

        freqlist = self.getFreqList(origcube)
        
        _ia.open(origcube)
        shp = _ia.shape()
        _ia.close()

        for i in range(0,len(freqlist)):

            _ia.open(origcube)
            im_orig = _ia.getchunk(blc=[0,0,0,i],trc=[shp[0],shp[1],0,i])
            _ia.close()

            _ia.open(modelcube)
            im_model = _ia.getchunk(blc=[0,0,0,i],trc=[shp[0],shp[1],0,i])
            _ia.close()

            _ia.open(psfcube)
            im_psf = _ia.getchunk(blc=[0,0,0,i],trc=[shp[0],shp[1],0,i])
            _ia.close()

            smoothedim = self.myconvolve(im_model[:,:,0,0], im_psf[:,:,0,0])

            if( np.nansum(im_orig)==0.0):
                smoothedim.fill(0.0)

            im_residual=im_psf.copy() ## Just to init the shape of this thing
            im_residual[:,:,0,0] = im_orig[:,:,0,0] - smoothedim

            _ia.open(residualcube)
            _ia.putchunk(im_residual, blc=[0,0,0,i])
            _ia.close()

#        ia_orig.close()
#        ia_model.close()
#        ia_residual.close()
#        ia_psf.close()

##################################################

    def myconvolve(self,im1,im2):
        
        t3 = time.time()
        
        shp = im1.shape
        pshp = (shp[0]*2, shp[1]*2)
        pim1 = np.zeros(pshp,'float')
        pim2 = np.zeros(pshp,'float')
        
        pim1[shp[0]-shp[0]//2 : shp[0]+shp[0]//2, shp[1]-shp[1]//2 : shp[1]+shp[1]//2] = im1
        pim2[shp[0]-shp[0]//2 : shp[0]+shp[0]//2, shp[1]-shp[1]//2 : shp[1]+shp[1]//2] = im2
        
        fftim1 = fftpack.fftshift( fftpack.fft2( fftpack.ifftshift( pim1 ) ) )
        fftim2 = fftpack.fftshift( fftpack.fft2( fftpack.ifftshift( pim2 ) ) )
        fftprod = fftim1*fftim2
        psmoothedim = np.real(fftpack.fftshift( fftpack.ifft2( fftpack.ifftshift( fftprod ) ) ) )
        
        smoothedim = psmoothedim[shp[0]-shp[0]//2 : shp[0]+shp[0]//2, shp[1]-shp[1]//2 : shp[1]+shp[1]//2]
        
        return smoothedim

##########################################
