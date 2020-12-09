################################################
# Task to make masks.
#  reorganized after diesuccsion on Oct1, 2012
#  with Jeurgen and Urvashi
#
#  modified by TT
# based on the original code, 
# v1.0: 2012.03.20, U.Rau
#
################################################
# Notes (to self) - TT 
# 1. expanding one mask to another 
#    e.g.) expanding a continuum mask (both image mask/boolean mask)
#          channel mask 
# 2. part of copy mode func.: merging of different types of masks 
#    e.g.) inpimage and inpmask are lists of the mask to be merged
#          output mask is written to either outimage or outmask as embedded
#           T/F mask 
# 3. copying mask to another or create a new one
#    regrid if necessary (i.e. if the coords are different) 
# ----------------------------------------------
# basic rules:
# for mask image (1/0 mask): as it is
# for internal mask : parent_imagename:mask_name
#
# For input,
# inpimage is the casa image
# - mode='list': list internal masks of inpimage
# - other mode: used as a template for output
#               if region files are specified -> make mask specifeid with the regions on to inpimage
#               output is '' => modified inpimage unless overwrite=F else exception
#
# if inpmask='': use inpimage as input mask
# if inpmask='mask0' or other embedded mask name of the inpimage, 
# use that T/F mask
# 
# =expand=
# case1: on the same image (outimage=''), expand mask image from 
# prev. run etc. No regriding. Use nearest chan mask
# image to expand.  
# 1.a: inpimage is clean image mask (1s and 0s)
#      i) outimage != inpimage and outmask='' => new expanded mask image to outimage
#     ii) outimage != inpimage and outmask!='' => convert expanded mask image to T/F mask to store inside outimage
#    iii) outimage ==inpimage and outmask='' => update input mask image by expanding it 
#     iv) outimage ==inpimage and outmask!=''=> update input image with the expanded T/F mask
# 1.b: if inpmask!='', do T/F mask to 1/0 image mask conversion, then do as in 1.a 
 
# case2: outimage is in diffirent coords. (need to regrid)
#
#################
# tests 
# 1. for input: mask image 
# 2. for input: mask/or regular image with internal mask
# 3. for input: mask image; for output: mask image with different spectral grid
# 4. for input: mask/regular image with internal mask; for output: image with 
#      internal mask with different spectral grid
###################
from __future__ import absolute_import
import os
import shutil
import numpy as np

# get is_CASA6 and is_python3
from casatasks.private.casa_transition import *
if is_CASA6:
    from casatools import image, regionmanager, imager, table, quanta
    from casatasks import casalog
    from .imtools import pixelmask2cleanmask
else:
    from taskinit import *
    from recipes.pixelmask2cleanmask import pixelmask2cleanmask

    image = iatool
    regionmanager = rgtool
    quanta = qatool
    table = tbtool
    imager = imtool

if is_python3:
    import subprocess
    subprocess_getoutput = subprocess.getoutput
else:
    import commands
    subprocess_getoutput = commands.getoutput

_ia = image()
_rg = regionmanager()
_qa = quanta()

pid = str(os.getpid())
debug = False 
#debug = True 

def makemask(mode,inpimage, inpmask, output, overwrite, inpfreqs, outfreqs):
    """
    make /manipulate masks

    """
    _ia = image()
    _rg = regionmanager()
    _im = imager()
    casalog.origin('makemask')
    #casalog.post("params(mode,inpimage,inpmask,output,overwrite)=",mode,inpimage,inpmask,output,overwrite)

    try:
        # temp files
        tmp_maskimage='__tmp_makemaskimage_'+pid
        tmp_outmaskimage='__tmp_outmakemaskimage_'+pid
        tmp_regridim='__tmp_regridim_'+pid

        #intial cleanup to make sure nothing left from previous runs
        tmplist = [tmp_maskimage,tmp_outmaskimage,tmp_regridim]
        cleanuptempfiles(tmplist)
   
        # do parameter check first
        # check names of inpimage, inpmask check for existance
        # inpimage == output (exact match) then check overwrite
        #   => T overwrite inpimage
        #   => F exception

        # check inpimage 
        if (['list','copy','expand'].count(mode)==1):
            if inpimage=='': raise ValueError("inpimage is empty")
            if not os.path.isdir(inpimage):
                raise RuntimeError("inpimage=%s does not exist" % inpimage)
       
        # === list mode ===
        if mode == 'list':
           inpOK=checkinput(inpimage)
           if inpOK: 
              
               if _ia.isopen(): _ia.close()
               _ia.open(inpimage)
               inmasklist=_ia.maskhandler('get')
               # now ia.maskhandler returns ['T'] if no internal mask is there...
               if inmasklist.count('T')!=0:
                   inmasklist.remove('T')
               if len(inmasklist) ==0:
                   casalog.post('No internal (T/F) masks were found in %s' % (inpimage),'INFO')
               else:
                   defaultmaskname=_ia.maskhandler('default')[0]
                   printinmasks=''
                   for mname in inmasklist:
                       if mname==defaultmaskname:
                           printinmasks+='\''+mname+'\''+'(default)'
                       else:
                           printinmasks+='\''+mname+'\''
                       if mname != inmasklist[-1]:
                           printinmasks+=', '
                 
                   casalog.post('Internal (T/F) masks in %s: %s' % (inpimage, printinmasks),'INFO')
               _ia.close()
 
        # === setdefaultmask mode ===
        elif mode == 'setdefaultmask':
            inpOK=checkinput(inpmask)
            if inpOK:
                (parentimage,bmask)=extractmaskname(inpmask)
                if bmask=='':
                    raise RuntimeError("Missing an internal mask name")
                if _ia.isopen(): _ia.close()
                _ia.open(parentimage)
                defaultmaskname=_ia.maskhandler('default')[0]
                inmasklist=_ia.maskhandler('get')
                if defaultmaskname==bmask:
                    casalog.post('No change. %s is already a default internal mask' % bmask, 'INFO')
                else:
                    _ia.maskhandler('set',bmask)
                    casalog.post('Set %s as a default internal mask' % bmask, 'INFO')
                    if len(inmasklist)>1:
                        casalog.post('Current internal masks are %s' % str(inmasklist), 'INFO')
                _ia.close()

        # === delete mode ===
        elif mode == 'delete':
            inpOK=checkinput(inpmask)
            if inpOK:
                (parentimage,bmask)=extractmaskname(inpmask)
                if bmask=='':
                    raise RuntimeError("Missing an internal mask name")
                _ia.open(parentimage)
                casalog.post('Deleting the internal mask, %s ' % bmask, 'INFO')
                defaultmaskname=_ia.maskhandler('default')[0]
                _ia.maskhandler('delete',bmask)
                inmasklist=_ia.maskhandler('get')
                if inmasklist.count('T')!=0:
                    inmasklist.remove('T')
                if len(inmasklist) !=0 and defaultmaskname==bmask:
                    _ia.maskhandler('set',inmasklist[0])
                    casalog.post('Set %s as a default internal mask' % inmasklist[0], 'INFO')
                    if len(inmasklist)>1:
                        casalog.post('Current internal masks are %s' % str(inmasklist), 'INFO')
          
                _ia.close()

        else:
           # PREPROCESS STAGE for mode='copy' and 'expand'
           #DEBUG
           #casalog.post("mode=",mode)
           # copy can have multiple input masks, expand has only one.
           # check inpimage, inpmask, output, overwrite
           # 
           storeinmask = False # used to check if output to a new internal mask
           isNewOutfile = False

           inpOK=checkinput(inpimage)
           if inpOK:
               (immask,inmask)=extractmaskname(inpimage)
          
           # seperate text files(region files), images(with full ':'), and direct region 
           # input mask(s)
           if inpmask=='':
              raise ValueError("Input errror. The inpmask parameter is not specified.")
           if type(inpmask)!=list: 
              inpmask=[inpmask]
           
           # check if inpmask contains region file or region specification
           rgfiles=[]
           imgfiles=[]
           rglist=[]
           bmasks=[]
           for masklet in inpmask:
               # is text file?
               if type(masklet)==str: # text file or image
                   if os.path.exists(masklet):
                       # region file or image
                       if (subprocess_getoutput('file '+masklet).count('directory')):
                          if os.path.exists(masklet+'/table.f1'):
                              #casalog.post("%s is not in a recognized format for inpmask, ignored." % masklet, 'WARN') 
                              raise ValueError("%s is not in a recognized format for inpmask" % masklet)
                          else:
                          # probably image file (with no mask extension)
                              imgfiles.append(masklet)
                       # text file (CRTF)
                       elif (subprocess_getoutput('file '+masklet).count('text')):
                          rgfiles.append(masklet)
                       else:
                          #casalog.post("%s does not recognized format for inpmask, ignored." % masklet, 'WARN')
                          raise ValueError("%s is not in a recognized format for inpmask" % masklet)
                   else:
                       # direct string specification
                       if masklet.count('[') and masklet.count(']'): # rough check on region specification 
                           rglist.append(masklet) 
                       # extract internal mask from the input image
                       else:
                           (parentim, mask)=extractmaskname(masklet)
                           if mask!='':
                               bmasks.append(masklet)
                           else:
                               raise ValueError("%s is not an existing file/image or a region format" % masklet)
      
           # expand allows only a string for inpmask
           if mode=='expand':
               if type(inpmask)==list:
                   inpmask=inpmask[0] 
           # check for overwrite condition
           if output=='':
               if overwrite: 
                   output=inpimage
               else: 
                   raise ValueError("output is not specified. If you want to overwrite inpimage, please set overwrite=True")

           if inpimage==output:
               #if overwrite:
               #    tmp_outmaskimage=tmp_maskimage
               #else:
               if not overwrite:
                   raise ValueError("output=inpimage. If you want to overwrite inpimage, please set overwrite=True")

           outparentim=output
           outbmask=''
           if os.path.isdir(output): 
               if not overwrite:
                    raise RuntimeError("output=%s exists. If you want to overwrite it, please set overwrite=True" % output)
           else:
               (outparentim, outbmask)=extractmaskname(output)
               if outbmask!='':
                   (parentimexist,maskexist)=checkinmask(outparentim,outbmask)    
                   if parentimexist and maskexist: 
                       if not overwrite:
                           raise ValueError("output=%s exists. If you want to overwrite it, please set overwrite=True" % output)
                       else:
                        casalog.post("Will overwrite the existing internal mask, %s in %s" % (outbmask,outparentim))
                        storeinmask=True

                   #if parentimexist and not maskexist:
                   else:
                        storeinmask = True
                        if not parentimexist: isNewOutfile=True
               else:
                  outparentim=output
               
           #casalog.post("param checks before branching out for mode=========="))
           #casalog.post("storeinmask = ",storeinmask)
           #casalog.post("output=",output, " is exist?=",os.path.isdir(output))
           #casalog.post("outparentim=",outparentim, " is exist?=",os.path.isdir(outparentim))

    # MAIN PROCESS for 'copy' or 'expand' mode
    # the following code is somewhat duplicated among the modes but keep separated from each mode
    # for now.... copy now handle merge as well
    # === old copy mode === NOW combined to 'merge mode'
#        #if mode=='copy':
#           #casalog.post("Copy mode")
#           needregrid=True
#           #if outimage=='':
#               #overwrite
#           #    outimage=inpimage
#
#           if not os.path.isdir(outimage):
#               needregrid=False
#
#          if inpmask!='':
#          # need to extract the mask and put in tmp_maskimage
#              pixelmask2cleanmask(imagename=inpimage, maskname=inpmask, maskimage=tmp_maskimage, usemasked=True)    
#          else:
#              shutil.copytree(inpimage, tmp_maskimage)
#           if needregrid:
#               casalog.post("Regridding...",'DEBUG1')
#               regridmask(tmp_maskimage,outimage,tmp_outmaskimage)
#               # regrid may produce <1.0 pixels for mask so be sure to its all in 1.0
#               #_ia.open(tmp_outmaskimage) 
#               #_ia.calc('iif (%s>0.0 && %s<1.0,1,%s)'%(tmp_outmaskimage,tmp_outmaskimage,tmp_outmaskimage))
#               #_ia.close()
#               #casalog.post("Copying regrid output=",tmp_outmaskimage)
#           else:
#               shutil.copytree(tmp_maskimage,tmp_outmaskimage)
#          if outmask!='':
#          #convert the image mask to T/F mask
#               if not os.path.isdir(outimage):
#                   shutil.copytree(inpimage,outimage)
#               #
#               _ia.open(outimage)
#              casalog.post("convert the output image mask to T/F mask")
#              _ia.calcmask(mask='%s<0.5' % tmp_outmaskimage,name=outmask,asdefault=True)
#               _ia.done()
#           else:
#               # if regridded - tmp_outmaskimage is created by regridmask
#               # if not, tmp_outmaskimage=outimage
#               _ia.open(tmp_outmaskimage)
#               _ia.rename(outimage,overwrite=True)
#               _ia.done()

                
    # === expand mode ===
        if mode=='expand':
            _rg = regionmanager()
            needtoregrid=False
            bychanindx=False

            try:
                # These coordsys objects will need to be closed, if created
                incsys = None
                inmaskcsys = None
                ocsys = None

                #casalog.post("expand mode main processing blocks...")
                # do not allow list in this mode (for inpimage and inpmask) - maybe this is redundant now
                if type(inpmask)==list:
                    raise TypeError('A list for inpmask is not allowed for mode=expand')

                # input image info, actually this will be output coordinates
                _ia.open(inpimage)
                inshp = _ia.shape()
                incsys = _ia.coordsys()
                _ia.close() 
                #casalog.post("inpimage=",inpimage," is exist?=",os.path.isdir(inpimage))
                #casalog.post(" inshp for inpimage=",inshp)

                # prepare working input mask image (tmp_maskimage)
                if debug: casalog.post("prepare working input image (tmp_maskimage)...")
                if inpmask!='': # inpmask is either image mask or T/F mask now
                  # need to extract the mask and put in tmp_maskimage
                  # Note: changed usemasked=F, so that True (unmasked) part to be used. CAS- 
                  # ==> tmp_maskimage is an input mask image
                    (parentimage,bmask)=extractmaskname(inpmask)
                    if bmask!='':
                        pixelmask2cleanmask(imagename=parentimage, maskname=bmask, maskimage=tmp_maskimage, usemasked=False)    
                        #_ia.open(tmp_maskimage)
                    else:
                        if debug: 
                            casalog.post("parentimage=" + parentimage + " exist?=" + os.path.isdir(parentimage))
                            casalog.post("tmp_maskimage=" + tmp_maskimage + " exist?=" + os.path.isdir(tmp_maskimage))
                        # copy of inpimage in tmp_maskimage
                        _ia.fromimage(outfile=tmp_maskimage, infile=parentimage)
                else:
                    raise ValueError("inpmask must be specified")
                if _ia.isopen(): _ia.close() 
                #setting up the output image (copy from inpimage or template)
                if not os.path.isdir(outparentim):
                    #shutil.copytree(inpimage,tmp_outmaskimage)
                    _ia.fromshape(outfile=tmp_outmaskimage,shape=inshp, csys=incsys.torecord()) 
                    _ia.close() 
                    needtoregrid=False
                else:
                    # if inpimage == output, tmp_outmaskimage is already created... 
                    if not os.path.isdir(tmp_outmaskimage):
                       shutil.copytree(outparentim,tmp_outmaskimage)
                if debug: casalog.post("done setting up the out image=" + tmp_outmaskimage)
                # if inpfreq/outfreq are channel indices (int) then
                # regrid in x y coords only and extract specified channel mask
                # to specified output channels. (no regriding in spectral axis)
                # if inpfreqs/outfreqs are velocity or freqs, 
                # it assumes it is expressed in the range with minval~maxval
                # create subimage of the input mask with the range,
                # do regrid with the subimage to output.
          
                # decide to regrid or not
                # 1. the case all channels are selected for input and output, simply regrid
                # 2. if inpfreqs and outfreqs are integers (= channel indices), regrid only in
                #    first and second axes (e.g. ra,dec) and no regridding along spectral axis
                # 3. if inpfreqs and outfreqs are ranges in freq or vel, make subimage and regrid
                _ia.open(tmp_maskimage)
                inmaskshp = _ia.shape()
                inmaskcsys = _ia.coordsys()
                _ia.close()
                regridmethod = 'linear'
                if 'spectral2' in inmaskcsys.torecord():
                    inspecdelt = inmaskcsys.torecord()['spectral2']['wcs']['cdelt']
                    _ia.open(tmp_outmaskimage)
                    ocsys=_ia.coordsys()
                    oshp=_ia.shape() 
                    _ia.close()
                    outspecdelt = ocsys.torecord()['spectral2']['wcs']['cdelt']
                    if outspecdelt < inspecdelt:
                       regridmethod='nearest'

                if inmaskshp[3]!=1 and ((inpfreqs==[] and outfreqs==[]) \
                    or (inpfreqs=='' and outfreqs=='')):
                    # unless inpimage is continuum, skip chan selection part and regrid 
                    needtoregrid=True
                    # detach input(tmp) image and open output tmp image
                #    _ia.open(tmp_outmaskimage)
                else: 
                #    if _ia.isopen():
                #        if _ia.name(strippath=True)!=tmp_maskimage:
                #            _ia.close()
                #            _ia.open(tmp_maskimage)
                #    else:
                #        _ia.open(tmp_maskimage)

                    #if inshp[3]!=1: casalog.post("inpmask is continuum..","INFO")
                    if inmaskshp[3]==1: casalog.post("inpmask is continuum..","INFO")
                    # selection by channel indices (number) 
                    # if both inpfreqs and outfreqs are int skip regridding
                    # if outfreqs is vel or freq ranges, try regridding 
                    if inpfreqs==[[]] or inpfreqs==[]: 
                        # select all channels for input
                        inpfreqs=list(range(inmaskshp[3]))

                    # check inpfreqs and outfreqs types
                    # index based
                    selmode='bychan'
                    if type(inpfreqs)==list:
                        if type(inpfreqs[0])==int:
                            if type(outfreqs)==list and (len(outfreqs)==0 or type(outfreqs[0])==int):
                                selmode='bychan'
                            elif type(outfreqs)==str:
                                #if inpfreqs[0]==0: #contintuum -allow index-type specification
                                if len(inpfreqs)==1: #contintuum -allow index-type specification
                                    selmode='byvf'
                                else:
                                    raise TypeError("Mixed types in infreqs and outfreqs are not allowed") 
                        else:
                            raise TypeError("Non-integer in inpfreq is not supported") 
                    # by velocity or frequency
                    elif type(inpfreqs)==str:
                        if type(outfreqs)!=str:
                            raise TypeError("Mixed types in infreqs and outfreqs") 
                        selmode='byvf'
                    else:
                        raise TypeError("Wrong type for infreqs")

                    # inpfreqs and outfreqs are list of int
                    # match literally without regridding.
                    if selmode=='bychan': 
                        casalog.post("selection of input and output ranges by channel")
                        
                        if _ia.isopen():
                            _ia.close()
                        if outfreqs==[] or outfreqs==[[]]:
                            outchans=[]
                        else:
                            outchans=outfreqs
                        expandchanmask(tmp_maskimage,inpfreqs,tmp_outmaskimage,outchans)
                        _ia.open(tmp_outmaskimage)

                    elif selmode=='byvf': # outfreqs are quantities (freq or vel)
                        casalog.post("selection of input/output ranges by frequencies/velocities")
                        
                        # do it for input mask image (not the template )
                        inpfreqlist = translatefreqrange(inpfreqs,inmaskcsys)
                        #casalog.post("inpfreqlist=",inpfreqlist)
                        # close input image
                        if _ia.isopen():
                            _ia.close()
                        
                        #regrid to output image coordinates
                        if len(inpfreqlist)==1: # continuum
                            #do not regrid, use input image
                            shutil.copytree(tmp_maskimage,tmp_regridim)
                        else:
                            if debug: casalog.post("regrid the mask,tmp_maskimage=" + tmp_maskimage + " tmp_regridim=" + tmp_regridim)
                            #shutil.copytree(tmp_maskimage,'before_regrid_tmp_maskimage')
                            regridmask(tmp_maskimage,inpimage,tmp_regridim,chanrange=inpfreqlist,method=regridmethod)
                            #regridmask(tmp_maskimage,inpimage,tmp_regridim,chanrange=inpfreqlist)
                            # find edge masks (nonzero planes)
                            ##shutil.copytree(tmp_regridim,'saved_'+tmp_regridim)
                            if _ia.isopen(): _ia.close()
                            _ia.open(tmp_regridim)
                            sh=_ia.shape()
                            chanlist = list(range(sh[3]))
                            indlo=0
                            indhi=0
                            for i in chanlist:
                                sl1=[0,0,0,i]
                                sl2=[sh[0]-1,sh[1]-1,sh[2]-1,i]
                                psum = _ia.getchunk(sl1,sl2).sum()
                                pmsum = _ia.getchunk(sl1,sl2,getmask=True).sum()
                                if pmsum!=0 and psum>0.0: 
                                    indlo=i
                                    break
                            chanlist.reverse()
                            for i in chanlist:
                                sl1=[0,0,0,i]
                                sl2=[sh[0]-1,sh[1]-1,sh[2]-1,i]
                                psum = _ia.getchunk(sl1,sl2).sum()
                                if psum>0.0: 
                                    indhi=i
                                    break
                            if indhi < indlo:
                                raise RuntimeError("Incorrectly finding edges of input masks! Probably some logic error in the code!!!") 
                            else:
                                casalog.post("Determined non-zero channel range to be "+str(indlo)+"~"+str(indhi), 'DEBUG1')

                        # find channel indices for given outfreqs
                        #_ia.open(tmp_outmaskimage)
                        #ocsys=_ia.coordsys()
                        #oshp=_ia.shape() 
                        outfreqlist = translatefreqrange(outfreqs,ocsys)
                        rtn=ocsys.findcoordinate('spectral')
                        px=rtn['pixel'][0]
                        wc=rtn['world'][0]
                        world=ocsys.referencevalue()
                        # assume chanrange are in freqs
                        world['numeric'][wc]=_qa.convert(_qa.quantity(outfreqlist[0]),'Hz')['value']
                        p1 = ocsys.topixel(world)['numeric'][px]
                        world['numeric'][wc]=_qa.convert(_qa.quantity(outfreqlist[1]),'Hz')['value']
                        p2 = ocsys.topixel(world)['numeric'][px]
                        casalog.post("translated channel indices:"+_qa.tos(outfreqlist[0])+"->p1="+str(p1)+\
                        " "+_qa.tos(outfreqlist[0])+"->  p2="+str(p2))
                        if len(inpfreqs)==1:
                            inpfreqchans=inpfreqs
                        elif inpfreqs.find('~'):
                            inpfreqchans=list(range(indlo,indhi+1))
                        else:
                            inpfreqchans=[indlo,indhi]
                        outfreqchans=list(range(int(round(p1)),int(round(p2))+1))
                        #casalog.post("inpfreqchans=",inpfreqchans)
                        #casalog.post("outfreqchans=",outfreqchans)
                        expandchanmask(tmp_regridim,inpfreqchans,tmp_outmaskimage,outfreqchans)
                        #shutil.copytree(tmp_regridim,'my_tmp_regrid') 
                        #shutil.copytree(tmp_outmaskimage,'my_tmp_outmaskimage') 

#                       usechanims={}  # list of input mask to be use for each outpfreq
#                       for i in outfreqchans:
#                           nearestch = findnearest(inpfreqchans,i)
#                           usechanims[i]=nearestch 
#                        # put masks from inp image channel by channel
#                       for j in outfreqs:
#                           pix = refchanchunk[usechanims[j]-refchanst]
#                           #_ia.putchunk(pixels=pix,blc=[inshp[0]-1,inshp[1]-1,0,j])
#                           _ia.putchunk(pixels=pix.transpose(),blc=[0,0,0,j])
                        needtoregrid=False

                        if _ia.isopen(): _ia.close()
                        _ia.open(tmp_outmaskimage)
                # 
                
                if needtoregrid:
                    # closing current output image
                    if _ia.isopen():
                        _ia.close()
                    _ia.open(tmp_maskimage)
                    #os.system('cp -r %s beforeregrid.im' % tmp_maskimage)
                    if os.path.isdir(tmp_outmaskimage):
                        #casalog.post("Removing %s" % tmp_outmaskimage)
                        shutil.rmtree(tmp_outmaskimage)
                    #regridmask(tmp_maskimage,outparentim,tmp_outmaskimage)
                    regridmask(tmp_maskimage,inpimage,tmp_outmaskimage,method=regridmethod)
                    _ia.remove()
                    #casalog.post("closing after regrid")
                    _ia.open(tmp_outmaskimage) # reopen output tmp image
                
                # for debugging
                #os.system('cp -r '+outparentim+" expandmode-copy-"+outparentim)
                #os.system('cp -r '+tmp_outmaskimage+" expandmode-copy-"+tmp_outmaskimage)
                if outbmask!='':
                    #convert the image mask to T/F mask
                    casalog.post("Convert the image mask to T/F mask",'INFO')
                    # regions will be masked if == 0.0 for a new outfile, if outfile exists 
                    # the pixel values inside specified mask is preserved and the rest is masked
                    if os.path.isdir(outparentim):
                      _ia.calcmask(mask='%s==1.0' % tmp_outmaskimage,name=outbmask,asdefault=True)
                    else:
                      _ia.calcmask(mask='%s!=0.0' % tmp_outmaskimage,name=outbmask,asdefault=True)
                if storeinmask:
                    isNewFile=False
                    if not os.path.isdir(outparentim):
                      makeEmptyimage(inpimage,outparentim)
                      isNewFile=True
                    _ia.open(outparentim)
                    if isNewFile:
                      _ia.set(1) 
                      # if output image exist its image pixel values will not be normalized the region
                      # outside input mask will be masked.
                    #check 
                    curinmasks = _ia.maskhandler('get') 
                    if outbmask in curinmasks:
                       if  not overwrite:
                           raise RuntimeError("Internal mask,"+outbmask+" exists. Please set overwrite=True.")
                       else:
                           _ia.maskhandler('delete',outbmask)
                    
                    _ia.maskhandler('copy',[tmp_outmaskimage+':'+outbmask, outbmask])
                    _ia.maskhandler('set',outbmask)
                    _ia.done()
                    casalog.post("Output the mask to %s in %s" % (outbmask,outparentim) ,"INFO")
                else:
                    ow = False
                    if  inpimage==output:
                        casalog.post("Updating "+output+" with new mask","INFO")
                        ow=True
                    else:
                        if os.path.isdir(outparentim):
                            casalog.post(outparentim+" exists, overwriting","INFO")
                            ow=True
                        else:
                            casalog.post("Output the mask to "+outparentim ,"INFO")
                    _ia.rename(outparentim,ow)
                    _ia.done()

            except Exception as instance:
                casalog.post("*** Error (1) *** %s" % instance, 'ERROR')
                if _ia.isopen():
                    _ia.close()
                _ia.done()
                raise
            finally:
                if inmaskcsys:
                    inmaskcsys.done()
                if incsys:
                    incsys.done()
                if ocsys:
                    ocsys.done()
                if os.path.exists(tmp_maskimage):
                    shutil.rmtree(tmp_maskimage)
                if os.path.exists(tmp_regridim):
                    shutil.rmtree(tmp_regridim)
                if os.path.exists(tmp_outmaskimage):
                    shutil.rmtree(tmp_outmaskimage)

    # === main process for copy mode: also does merge of masks ===
        # copy is a just special case of merge mode
        # CHANGE:
        # all input masks should be specified in inpmask 
        # type of inpmask accepted: 1/0 mask, T/F mask, region file, and region expression in a string 
        # already stored internally in seperate lists
        #   rgfiles - region files (binary or CRTF-format text file)
        #   imgfiles - 1/0 image masks
        #   rglist - region expression in strings
        #   bmasks - T/F internal masks
        #
        # avaialble parameters: inpimage (string) , inpmask (list/string), output(string)
        # input inpimage as a template or image used for defining regions when it is given in inpmask 
        # inpmask as list of all the masks to be merged (image masks, T/F internal masks, regions)

        if mode=='copy':
            sum_tmp_outfile='__tmp_outputmask_'+pid
            tmp_inmask='__tmp_frominmask_'+pid
            tmp_allrgmaskim='__tmp_fromAllRgn_'+pid
            tmp_rgmaskim='__tmp_fromRgn_'+pid
            # making sure to remove pre-existing temp files
            cleanuptempfiles([sum_tmp_outfile, tmp_inmask, tmp_allrgmaskim, tmp_rgmaskim])
            usedimfiles=[]
            usedbmasks=[]
            usedrgfiles=[]
            usedrglist=[]
            # This coordsys object used in several places below will need to be closed, if created
            tcsys = None
            try:
                # check outparentim - image part of output and set as a template image
                if not (os.path.isdir(outparentim) or (outparentim==inpimage)):
                    # Output is either a new image or the same as inpimage

                    # figure out which input mask to be used as template
                    # if inpimage is defined use the first one else try the first one
                    # inpmask
                    #if output=='':
                    #    if type(inpimage)==list:
                    #         raise Exception, "inputimage must be a string"
                    #    elif type(inpimage)==str:
                    #         outimage=inpimage
                    #    casalog.post("No outimage is specified. Will overwrite input image: "+outimage,'INFO')

                    #if type(inpimage)==list and len(inpimage)!=0:
                    #    tmp_template=inpimage[0]
                    #elif inpimage!='' and inpimage!=[]:
                    #    tmp_template=inpimage # string
                    #tmp_template=inpimage # string
                    #else:
                    #     if type(inpmask)==list and len(inpmask)!=0:
                    #         fsep=inpmask[0].rfind(':')
                    #         if fsep ==-1:
                    #            raise IOError, "Cannot resolve inpmask name, check the format"
                    #        else:
                    #            tmp_template=inpmask[0][:inpmask[0].rfind(':')]
                    #    elif inpmask!='' and inpmask!=[]:
                    #        # this is equivalent to 'copy' the inpmask
                    #        tmp_template=inpmask #string
                    
                    # create an empty image with the coords from inpimage
                    makeEmptyimage(inpimage,sum_tmp_outfile)    
                    #casalog.post("making an empty image from inpimage to sum_tmp_outfile")
                else:
                    #use output image(either the same as the input image or other existing image) 
                    #  - does not do zeroeing out, so output image is only modified *for outbmask!=''*    
                    shutil.copytree(outparentim,sum_tmp_outfile)
                    # temporary clear out the internal masks from the working image
                    if _ia.isopen(): _ia.close()
                    _ia.open(sum_tmp_outfile)
                    if (len(imgfiles) or len(rglist) or len(rgfiles)):
                        # do zeroeing out working base image for masks 
                        _ia.set(0)
                    origmasks = _ia.maskhandler('get') 
                    _ia.maskhandler('delete',origmasks)
                    _ia.close()
                     
                if len(imgfiles)>0:
                    # summing all the images
                    casalog.post('Summing all mask images in inpmask and  normalized to 1 for mask','INFO')
                    for img in imgfiles:
                        #tmpregrid='__tmp_regrid.'+img
                        dirname = os.path.dirname(img)
                        basename = os.path.basename(img)
                        tmpregrid= dirname+'/'+'__tmp_regrid.'+basename if len(dirname) else '__tmp_regrid.'+basename
                        tmpregrid+='_'+pid
                        if os.path.exists(tmpregrid):
                            shutil.rmtree(tmpregrid)
                        # regrid to output image coords
                        try:
                            regridmask(img,sum_tmp_outfile,tmpregrid)
                            addimagemask(sum_tmp_outfile,tmpregrid)
                            usedimfiles.append(img)
                        finally:
                            shutil.rmtree(tmpregrid)
                    # get boolean masks
                    #  will work in image(1/0) masks
                
                    if debug: 
                        casalog.post(("imgfiles=" + imgfiles))
                        shutil.copytree(sum_tmp_outfile,sum_tmp_outfile+"_imagemaskCombined")
                       
                if len(bmasks)>0:
                    casalog.post('Summing all T/F mask in inpmask and normalized to 1 for mask','INFO')
                    for msk in bmasks:
                        (imname,mskname) = extractmaskname(msk)
                        #if msk.find(':')<0:
                        #    # assume default mask
                        #    msk=msk+':mask0'
                        #imname=msk[:msk.rfind(':')]
                        if _ia.isopen(): _ia.close()
                        _ia.open(imname)
                        inmasks=_ia.maskhandler('get')
                        _ia.close()
                        if not inmasks.count(mskname):
                            raise TypeError(mskname+" does not exist in "+imname+" -available masks:"+str(inmasks))
                        # move T/F mask to image mask
                        # changed to usemasked=False as of CAS-5443  

                        pixelmask2cleanmask(imname, mskname, tmp_inmask, False)    
                        cleanuptempfiles(['__tmp_fromTFmask_'+pid]) 
                        regridmask(tmp_inmask,sum_tmp_outfile,'__tmp_fromTFmask_'+pid)
                        # for T/F mask do AND operation
                        _ia.open(sum_tmp_outfile)
                        if _ia.statistics()['sum'][0] != 0:
                            multiplyimagemask(sum_tmp_outfile,'__tmp_fromTFmask_'+pid)
                        else:
                            addimagemask(sum_tmp_outfile,'__tmp_fromTFmask_'+pid)
                        _ia.close()
                        usedbmasks.append(msk)
                        # need this temp file for the process later
                        ###shutil.rmtree('__tmp_fromTFmask') 
                        if os.path.isdir(tmp_inmask):
                           shutil.rmtree(tmp_inmask) 
                        # if overwriting to inpimage and if not writing to in-mask, delete the boolean mask
                        if outparentim==inpimage and inpimage==imname:
                            if outbmask=="":
                                _ia.open(imname)
                                _ia.maskhandler('delete',[mskname])
                                _ia.close()
                        _ia.open(imname)
                        _ia.close()
                      
                if debug: casalog.post("check rgfiles and rglist")

                if len(rgfiles)>0 or len(rglist)>0:
                    # create an empty image with input image coords.
                    #casalog.post("Using %s as a template for regions" % inpimage )
                    if _ia.isopen(): _ia.close()
                    _ia.open(inpimage)
                    tshp=_ia.shape()
                    tcsys=_ia.coordsys() 
                    _ia.fromshape(tmp_allrgmaskim,shape=tshp, csys=tcsys.torecord(),overwrite=True)
                    _ia.done()    
                    if os.path.isdir(tmp_rgmaskim):
                       shutil.rmtree(tmp_rgmaskim)
                    shutil.copytree(tmp_allrgmaskim,tmp_rgmaskim)

                    if len(rgfiles)>0:
                        # Here only CRTF format file is expected 
                        nrgn=0
                        for rgn in rgfiles:
                            subnrgn=0
                            with open(rgn) as f:
                                # check if it has '#CRTF' in the first line
                                firstline=f.readline()
                                if firstline.count('#CRTF')==0:
                                    raise Exception("Input text file does not seems to be in a correct format \
                                                              (must contains #CRTF in the first line)")
                                else:     
                                    try:
                                        # reset temp mask image
                                        if _ia.isopen(): _ia.close()
                                        _ia.open(tmp_rgmaskim)
                                        _ia.set(pixels=0.0)
                                        _ia.close()
                                        #casalog.post... "tshp=",tshp
                                        #casalog.post... "tcsys.torecord=",tcsys.torecord()
                                        inrgn=_rg.fromtextfile(rgn, tshp, tcsys.torecord())
                                        #casalog.post... "inrgn=",inrgn
                                        _im.regiontoimagemask(tmp_rgmaskim,region=inrgn)
                                        addimagemask(tmp_allrgmaskim,tmp_rgmaskim)
                                        #shutil.rmtree(tmp_rgmaskim)
                                        subnrgn +=1
                                        nrgn +=1
                                    except:
                                        break
                            if subnrgn>0:
                                usedrgfiles.append(rgn)    
                        casalog.post("Converted %s regions from %s region files to image mask" % (nrgn,len(rgfiles)),"INFO")
                                                
                    if debug: casalog.post("processing rglist...")
                    if len(rglist)>0:
                        #casalog.post( "Has rglist...")
                        nrgn=0
                        for rgn in rglist:
                            # reset temp mask image
                            if _ia.isopen(): _ia.close()
                            _ia.open(tmp_rgmaskim)
                            _ia.set(pixels=0.0)
                            _ia.close()
                            inrgn=_rg.fromtext(rgn, tshp, tcsys.torecord())
                            _im.regiontoimagemask(tmp_rgmaskim,region=inrgn)
                            addimagemask(tmp_allrgmaskim,tmp_rgmaskim)
                            #shutil.rmtree(tmp_rgmaskim)
                            usedrglist.append(rgn)
                            nrgn+=1
                        casalog.post("Converted %s regions to image mask" % (nrgn),"INFO")
                
                        
                    if debug: casalog.post("Regirdding...")
                    # regrid if necessary
                    regridded_mask='__tmp_regrid_allrgnmask_'+pid
                    regridmask(tmp_allrgmaskim, sum_tmp_outfile,regridded_mask)
                    addimagemask(sum_tmp_outfile,regridded_mask)
                    #shutil.rmtree('__tmp_regridded_allrgnmask')
                    casalog.post("Added mask based on regions to output mask","INFO")
                    #cleanup
                    for tmpfile in [tmp_allrgmaskim,tmp_rgmaskim,regridded_mask]:
                        if os.path.isdir(tmpfile):
                            shutil.rmtree(tmpfile)
                                        
                # merge the bmasks with AND
                if os.path.exists('__tmp_fromTFmask_'+pid) and len(bmasks)>0:
                    if _ia.isopen(): _ia.close()
                    _ia.open(sum_tmp_outfile)
                    if _ia.statistics()['sum'][0]!=0:
                       multiplyimagemask(sum_tmp_outfile,'__tmp_fromTFmask_'+pid)
                    else:
                       addimagemask(sum_tmp_outfile,'__tmp_fromTFmask_'+pid)
                    _ia.done()
                    shutil.rmtree('__tmp_fromTFmask_'+pid) 
                if outbmask!='':
                    casalog.post('Putting mask in T/F','INFO')
                    if _ia.isopen(): _ia.close()
                    try:
                        _ia.open(sum_tmp_outfile)
                        _ia.calcmask(mask='%s==1.0' % sum_tmp_outfile,name=outbmask,asdefault=True)
                    # mask only pixel == 0.0 (for a new outfile), mask region !=1.0 and preserve
                    # the pixel values if outfile exists
                    #if os.path.isdir(outparentim):
                    #  _ia.calcmask(mask='%s==1.0' % sum_tmp_outfile,name=outbmask,asdefault=True)
                    #else:
                    #  _ia.calcmask(mask='%s!=0.0' % sum_tmp_outfile,name=outbmask,asdefault=True)
                    finally:
                        _ia.done()
                if debug: shutil.copytree(sum_tmp_outfile,sum_tmp_outfile+"_afterCoverttoTFmask")
                # if outfile exists initially outfile is copied to sum_tmp_outfile
                # if outfile does not exist initially sum_tmp_outfile is a copy of inpimage
                # so rename it with overwrite=T all the cases
                #casalog.post("open sum_tmp_outfile=",sum_tmp_outfile)
                if storeinmask:
                    if debug: casalog.post("Storeinmask......")
                    # by a request in CAS-6912 no setting of 1 for copying mask to the 'in-mask'
                    # (i.e. copy the values of inpimage as well for this mode)
                    # Replace
                    #isNewfile = False
                    #if not os.path.isdir(outparentim):
                    if isNewOutfile:
                      #makeEmptyimage(inpimage,outparentim)
                      #isNewfile=True

                      shutil.copytree(inpimage,outparentim)
                    if _ia.isopen(): _ia.close()
                    _ia.open(outparentim)
                    if isNewOutfile:
                      oldmasklist = _ia.maskhandler('get')
                      if oldmasklist[0]!='T':
                        # clean up existing internal masks for the copied image 
                        if outbmask in oldmasklist:
                          _ia.maskhandler('delete',outbmask)
                    if (maskexist and overwrite): 
                      if debug: casalog.post("outbmask=" + outbmask + " exists... deleting it")
                      _ia.maskhandler('delete',outbmask)    
                    _ia.maskhandler('copy',[sum_tmp_outfile+':'+outbmask, outbmask])    
                    _ia.maskhandler('set',outbmask)
                    _ia.done()
                    outputmsg="to create an output mask: %s in %s" % (outbmask,outparentim)
                else:
                    if debug: casalog.post("store as an image mask......")
                    if _ia.isopen(): _ia.close()
                    _ia.open(sum_tmp_outfile) 
                    _ia.rename(outparentim,overwrite=True) 
                    _ia.done()
                    outputmsg="to create an output mask: %s " % outparentim

                casalog.post("Merged masks from:","INFO")
                if len(usedimfiles)>0:
                    casalog.post("mask image(s): "+str(usedimfiles),"INFO")
                if len(usedbmasks)>0:
                    casalog.post("internal mask(s): "+str(usedbmasks),"INFO")
                if len(usedrgfiles)>0:
                    casalog.post("region txt file(s): "+str(usedrgfiles),"INFO")
                if len(usedrglist)>0:
                    casalog.post("region(s) from direct input: "+str(usedrglist),"INFO")
                casalog.post(outputmsg,"INFO")
                
            except Exception as instance:
                raise RuntimeError("*** Error (2), in mode copy: *** %s" % instance)
            finally:
                if os.path.exists(sum_tmp_outfile):
                    shutil.rmtree(sum_tmp_outfile)
                if os.path.exists(tmp_inmask):
                    shutil.rmtree(tmp_inmask)
                if os.path.exists(tmp_allrgmaskim):
                    shutil.rmtree(tmp_allrgmaskim)
                if os.path.exists(tmp_rgmaskim):
                    shutil.rmtree(tmp_rgmaskim)
                if tcsys:
                    tcsys.done()


                if type(inpimage)==list:
                   for im in inpimage:
                       basename = os.path.basename(im)
                       dirname = os.path.dirname(im)
                       tempregridname = dirname+'/__tmp_regrid.'+basename if len(dirname) else '__tmp_regrid.'+basename
                       tempregridname+='_'+pid
                       if os.path.isdir(tempregridname):
                            shutil.rmtree(tempregridname)
                       
             
                 
    # === draw mode ===
    # disabled - drawmaskinimage (working with viewer) a bit flaky
    # when run in succession.
    #    if mode=='draw':
    #        #call drawmaskinimage
    #        from recipes.drawmaskinimage import drawmaskinimage
    #        drawmaskinimage(inpimage,outmask)
            

    finally:
        # final clean up 
        if os.path.isdir(tmp_maskimage):
            shutil.rmtree(tmp_maskimage)
        if os.path.isdir(tmp_outmaskimage):
            shutil.rmtree(tmp_outmaskimage)
        if os.path.isdir(tmp_regridim):
            shutil.rmtree(tmp_regridim)
        if os.path.exists('__tmp_fromTFmask_'+pid):
            shutil.rmtree('__tmp_fromTFmask_'+pid)

def findnearest(arr, val):
    if type(arr)==list:
        arr = np.array(arr) 
    indx = np.abs(arr - val).argmin()
    return arr[indx] 

def regridmask(inputmask,template,outputmask,axes=[3,0,1],method='linear',chanrange=None):
    '''
    Regrid input mask (image) to output mask using a template.
    Currently the template must be a CASA image.
    The default interpolation method is set to 'linear' (since 'nearest'
    sometime fails).
    '''
    #casalog.post("Regrid..")
    #casalog.post("inputmask=",inputmask," template=",template," outputmask=",outputmask)
    if not os.path.isdir(template):
        raise OSError("template image %s does not exist" % template)
    
    _ia = image()
    try:
        _tb = table()
        inputmaskcopy = "_tmp_copy_"+os.path.basename(inputmask)
        cleanuptempfiles([inputmaskcopy])
        shutil.copytree(inputmask,inputmaskcopy)
        _ia.open(template)
        ocsys = _ia.coordsys()
        oshp = _ia.shape()
    finally:
        _ia.done()
    _tb.open(template)
    defTelescope = _tb.getkeywords()['coords']['telescope']
    _tb.close()
    _tb.open(inputmaskcopy, nomodify=False) 
    keys = _tb.getkeywords()  
    if keys['coords']['telescope']=="UNKNOWN":
        if defTelescope =="UNKNOWN":
            raise ValueError("UNKNOWN Telescope for %s " % inputmask)
        else:
            keys['coords']['telescope']=defTelescope
    _tb.putkeywords(keys)     
    _tb.close()

    if _ia.isopen(): _ia.close()
    _ia.open(inputmaskcopy)
    # check axis order, if necessary re-interprete input axes correctly 
    # assumed order of axes 
    reforder=['Right Ascension', 'Declination', 'Stokes', 'Frequency']
    axisorder=_ia.summary(list=False)['axisnames'].tolist()
    # check if all 4 axes exist
    errmsg = ""
    for axname in reforder:
      if axisorder.count(axname) == 0:
        errmsg += axname+" "
    if len(errmsg) != 0:
      errmsg = "There is no "+errmsg+" axes inpimage. ia.adddegaxis or importfits with defaultaxes=True can solve this problem"    
      raise Exception(errmsg)

    tmp_axes=[]
    for axi in axes:
        tmp_axes.append(axisorder.index(reforder[axi]))        
    axes=tmp_axes
    if type(chanrange)==list and len(chanrange)==2:
        incsys=_ia.coordsys()
        spaxis=incsys.findcoordinate('spectral')['world']
        # create subimage based on the inpfreqs range
        inblc=chanrange[0]
        intrc=chanrange[1]
        casalog.post("Regridmask: spaxis=%s, inblc=%s, intrc=%s" % (spaxis,inblc,intrc), 'DEBUG1')
        rgn = _rg.wbox(blc=inblc,trc=intrc,pixelaxes=spaxis.tolist(),csys=incsys.torecord())
        incsys.done()
    else:
        rgn={}     
    # for continuum case
    ir = None
    if oshp[tmp_axes[0]]==1:
       axes=[0,1]
    try:
        #check for an appropriate decimation factor
        min_axlen=min(oshp[:2])
        if min_axlen < 30:
            decfactor=min_axlen//3
            if decfactor==0: decfactor=1
        else:
            decfactor=10
        ir=_ia.regrid(outfile=outputmask,shape=oshp,csys=ocsys.torecord(),axes=axes,region=rgn,method=method,decimate=decfactor)       
        
    except:
        pass
    finally:
        ocsys.done()
        _ia.remove()
        _ia.done()
        # to ensure to create 1/0 mask image
        #ir.calc('iif (%s>0.0 && %s<1.0,1,%s)'%(outputmask,outputmask,outputmask))
        # treat everything not = 0.0 to be mask
        if (ir):
            try:
                ir.calc('iif (abs("%s")>0.0,1,"%s")'%(outputmask,outputmask),False)
            finally:
                ir.done()
        if os.path.isdir(inputmaskcopy):
            shutil.rmtree(inputmaskcopy)

def addimagemask(sumimage, imagetoadd, threshold=0.0):
    """
    add image masks (assumed the images are already in the same coordinates)
    """
    _ia = image()
    try:
        #casalog.post("addimagemask: sumimage=",sumimage," imagetoadd=",imagetoadd)
        _ia.open(sumimage)
        _ia.calc('iif ("'+imagetoadd+'">'+str(threshold)+',("'+sumimage+'"+"'+imagetoadd+'")/("'+sumimage+'"+"'+imagetoadd+'"),"'+sumimage+'")',False)
        # actually should be AND?
        #_ia.calc('iif ('+imagetoadd+'>'+str(threshold)+','+sumimage+'*'+imagetoadd+','+sumimage+')')
        #_ia.calc('iif ('+imagetoadd+'>'+str(threshold)+',('+sumimage+'*'+imagetoadd+')/('+sumimage+'*'+imagetoadd+'),'+sumimage+')')
    finally:
        _ia.close()
    
def multiplyimagemask(sumimage, imagetomerge):
    """
    multiple image masks (do AND operation, assumed the images are already in the same coordinates)
    to use for merging of two image masks originated from T/F masks or merging between mask image
    and a T/F mask originated mask image
    """
    _ia = image()
    try:
        _ia.open(sumimage)
        _ia.calc('iif ("'+imagetomerge+'"!=0.0,("'+sumimage+'"*"'+imagetomerge+'"),0.0 )',False)
        _ia.calc('iif ("'+sumimage+'"!=0.0,("'+sumimage+'")/("'+sumimage+'"),"'+sumimage+'")',False)
    finally:
        _ia.close()

def expandchanmask(inimage,inchans,outimage,outchans):
    """
    expand masks in channel direction,and insert then
    to output image with the same coordinates (post-regridded)
    only differ by channels
    """
    _ia = image()
    # input image
    _ia.open(inimage)
    inshp=_ia.shape()
    refchanst=inchans[0]
    refchanen=inchans[-1]
    #casalog.post("refchanst=",refchanst," refchanen=",refchanen," inshp=",inshp," inchans=",inchans)
    slst = [0,0,0,refchanst]
    slen = [inshp[0]-1,inshp[1]-1,0,refchanen]
    casalog.post("getting chunk at blc="+str(slst)+" trc="+str(slen),'DEBUG1')
    refchanchunk=_ia.getchunk(blc=slst,trc=slen)
    refchanchunk=refchanchunk.transpose()
    _ia.close()
    #casalog.post("refchanchunk:shape=",refchanchunk.shape)

    _ia.open(outimage)
    # need find nearest inchan
    # store by chan indices (no regrid)
    outshp=_ia.shape()
    if outchans==[]:
        #select all channels
        outchans=list(range(outshp[3]))
    usechanims={}  # list of input mask to be use for each outpfreq
    for i in outchans:
        nearestch = findnearest(inchans,i)
        usechanims[i]=nearestch
    #casalog.post("usechanims=",usechanims)
    casalog.post("Mapping of channels: usechanims="+str(usechanims),'DEBUG1')
    for j in outchans:
        pix = refchanchunk[usechanims[j]-refchanst]
        #casalog.post("pix=",pix)
        #casalog.post("pix.shape=",pix.shape)
        #casalog.post("inshp=",inshp, ' j=',j)
        #_ia.putchunk(pixels=pix,blc=[inshp[0]-1,inshp[1]-1,0,j])
        _ia.putchunk(pixels=pix.transpose(),blc=[0,0,0,j])
        #casalog.post("DONE putchunk for j=", j)
    _ia.done()

def translatefreqrange(freqrange,csys):
    """
    convert the range in list
    mainly for frequeny and velocity range determination
    """
    if type(freqrange)==list and type(freqrange[0])==int:
        #do nothing
        return freqrange
    elif type(freqrange)==str:
        freqlist=freqrange.split('~') 
        for i in list(range(len(freqlist))):
            if freqlist[i].find('m/s') > -1:
               fq = _qa.quantity(freqlist[i])
               vf=csys.velocitytofrequency(value=fq['value'],velunit=fq['unit'])
               freqlist[i]=str(vf[0])+'Hz'
        return freqlist
    else:
        raise TypeError("Cannot understand frequency range")

def checkinput(inpname):
    """
    do existance check on image and internal mask 
    """
    _ia = image()
    (parentimage,tfmaskname)=extractmaskname(inpname)
    (parentimexist,tfmaskexist)=checkinmask(parentimage,tfmaskname)
    if parentimexist:
        if tfmaskname=='':
            return True # only the image
        else:
            if not tfmaskexist: 
                _ia.open(parentimage)
                inmasklist=_ia.maskhandler('get')
                _ia.close()
                raise Exception("Cannot find the internal mask, %s. Candidate mask(s) are %s" % (tfmaskname, str(inmasklist)))
            else:
                return True # image mask and internal mask
    else:
        raise Exception("Cannot find the image=%s" % parentimage) 
   

def checkinmask(parentimage,tfmaskname):
    """
    check existance of the internal mask
    """
    _ia = image()
    if os.path.isdir(parentimage):
        if tfmaskname!='':
            _ia.open(parentimage)
            inmasks=_ia.maskhandler('get')
            _ia.done()
            if not any(tfmaskname in msk for msk in inmasks):
               return (True, False)
            else:
               return (True, True) # image mask and internal mask
        else:
            return (True, False)
    else:
       return (False,False) 

def extractmaskname(maskname):
    """
    split out imagename and maskname from a maskname string
    returns (parentimage, internalmask)
    """
    # the image file name may contains ':' some cases
    # take last one in split list as an internal mask name

    # Try to avoid issues with ':' included in paths when given absolute paths
    dirname = None
    if os.path.isabs(maskname):
        dirname = os.path.dirname(maskname)
        maskname = os.path.basename(maskname)

    indx = maskname.find(':') 
    for i in list(range(len(maskname))):
        if indx>-1:
            indx += maskname[indx+1:].find(':') 
            indx +=1
        else:
            break
    if indx != -1: 
        parentimage = maskname[:indx]
        maskn = maskname[indx+1:]
    else:
        parentimage = maskname
        maskn = ''

    if dirname:
        parentimage = os.path.join(dirname, parentimage)

    return (parentimage, maskn)

def makeEmptyimage(template,outimage):
    """
    make an empty image with the coords
    from template
    """

    _ia = image()
    _ia.open(template)
    inshp=_ia.shape()
    incsys=_ia.coordsys()
    _ia.fromshape(outimage,shape=inshp,csys=incsys.torecord())
    incsys.done()
    _ia.done()

def cleanuptempfiles(tempfilelist):
    """
    clean up tempfilelist
    """
    for fn in tempfilelist:
        if os.path.isdir(fn):
            shutil.rmtree(fn)
        elif os.path.isfile(fn):
           os.remove(fn)
 
