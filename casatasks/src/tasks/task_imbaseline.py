# image-based line finding and baseline subtraction.

import os
import sys
import shutil
from typing import List
import numpy as np
from numpy import array, uint64
from collections import Counter

from casatools import image, table, regionmanager, quanta, singledishms
from casatools import ms as mstool
from casatasks import casalog
from casatasks.private import sdutil
from casatasks.private.sdutil import table_manager
from casatasks.private.ialib import write_image_history

ia = image()
qa = quanta()
ms = mstool()
sdms = singledishms()

class ImBaselineVals:
    """ImBaselineVals is to handle input data set of imbaseline task.
    The task contains five parts, imsmooth, MS creating, sdsmooth, sdbaseline, image output.
    Each parts require different parameters, so the task has a very large number of parameters.
    In order to simplify the handling of parameters, all parameters are held in an instance of 
    this class, all methods get/set all parameters from it.
    """

    def __init__(self, imagename: str = '', linefile: str = '', output_cont: bool = False, bloutput: str = '', maskmode: str = 'list',
                chans: str = '', thresh: float = 5.0, avg_limit: int = 4, minwidth: int = 4, edge: List[int] = [0, 0],
                blfunc: str = 'poly', order: int = 5, npiece: int = 3, applyfft: bool = True, fftthresh: float = 3.0,
                addwn: List[int] = [0], rejwn: List[int] = [], blparam: str = '', clipniter: int = 0, clipthresh: float = 3.0,
                dirkernel: str = 'none', major: str = '', minor: str = '', pa: str = '', kimage: str = '',
                scale: float = -1.0, spkernel: str = 'none', kwidth: int = 5):
        self.imagename = imagename
        self.linefile = linefile
        self.output_cont = output_cont

        # sdbaseline
        self.sdbaseline_bloutput = bloutput
        self.sdbaseline_maskmode = maskmode.lower() # list(default)/auto
        self.sdbaseline_chans = chans           # maskmode = list
        self.sdbaseline_thresh = thresh         # maskmode = auto
        self.sdbaseline_avg_limit = avg_limit   # maskmode = auto
        self.sdbaseline_minwidth = minwidth     # maskmode = auto
        self.sdbaseline_edge = edge             # maskmode = auto
        self.sdbaseline_blfunc = blfunc.lower() # poly(default)/chebyshev/cspline/sinusoid/variable
        self.sdbaseline_order = order           # blfunc = poly/chebyshev
        self.sdbaseline_npiece = npiece         # blfunc = cspline
        self.sdbaseline_applyfft = applyfft     # blfunc = sinusoid
        self.sdbaseline_fftmethod = 'fft'       # blfunc = sinusoid
        self.sdbaseline_fftthresh = fftthresh   # blfunc = sinusoid
        self.sdbaseline_addwn = addwn           # blfunc = sinusoid
        self.sdbaseline_rejwn = rejwn           # blfunc = sinusoid
        self.sdbaseline_blparam = blparam       # blfunc = variable
        self.sdbaseline_clipniter = clipniter   # blfunc = variable
        self.sdbaseline_clipthresh = clipthresh # blfunc = variable
        self.sdbaseline_antenna = ''
        self.sdbaseline_field = ''
        self.sdbaseline_spw = ''
        self.sdbaseline_timerenge = ''
        self.sdbaseline_scan = ''
        self.sdbaseline_pol = ''
        self.sdbaseline_intent = ''
        self.sdbaseline_reindex = True
        self.sdbaseline_blmode = 'fit'
        self.sdbaseline_dosubtract = True
        self.sdbaseline_blformat = 'text'
        self.sdbaseline_bloutput = ''
        self.sdbaseline_updateweight = False
        self.sdbaseline_sigmavalue = 'stddev'   # maybe not use
        self.sdbaseline_showprogress = False    # not use
        self.sdbaseline_minnrow = 1000          # not use

        # imsmooth
        self.imsmooth_kernel = dirkernel        # none(default)/gaussian/boxcar/image
        self.imsmooth_major = major             # dirkernel = gaussian/boxcar
        self.imsmooth_minor = minor             # dirkernel = gaussian/boxcar
        self.imsmooth_pa = pa                   # dirkernel = gaussian/boxcar
        self.imsmooth_kimage = kimage           # dirkernel = image
        self.imsmooth_scale = scale             # dirkernel = image
        self.imsmooth_targetres = False
        self.imsmooth_mask = ''
        self.imsmooth_beam = {}
        self.imsmooth_region = ''
        self.imsmooth_box = ''
        self.imsmooth_chans = ''
        self.imsmooth_stokes = ''
        self.imsmooth_stretch = False

        # sdsmooth
        self.sdsmooth_kernel = spkernel         # none(default)/gaussian/boxcar
        self.sdsmooth_kwidth = kwidth           # gaussian/boxcar
        self.sdsmooth_spw = ''
        self.sdsmooth_field = ''
        self.sdsmooth_antenna = ''
        self.sdsmooth_timerange = ''
        self.sdsmooth_scan = ''
        self.sdsmooth_pol = ''
        self.sdsmooth_intent = ''
        self.sdsmooth_reindex = True

        # imbaseline local
        self.temporary_vis = 'temp.ms'
        self.imsmooth_output = 'imsmooth_output.image'
        self.sdsmooth_output = 'sdsmooth_output.ms'
        self.sdbaseline_output = 'sdbaseline_output.ms'
        self.output_cont_file = os.path.basename(self.imagename)+".cont"
        self.datacolumn = 'DATA'
        self.overwrite = True

        self.__check_args()

    def __check_args(self):
        self.__check_arg_imagename()
        self.__check_arg_linefile()
        self.__check_arg_maskmode()
        self.__check_arg_blfunc()
        self.__check_arg_dirkernel()
        self.__check_arg_spkernel()

    def __check_arg_imagename(self):
        if not os.path.exists(self.imagename):
            raise ValueError(f'Error: file {self.imagename} is not found.', 'SEVERE')

    def __check_arg_linefile(self):
        if self.linefile == '':
            self.linefile = os.path.basename(self.imagename).rstrip('/') + '_bs'
        if os.path.exists(self.linefile):
            raise ValueError(f'Error: file {self.linefile} already exists, please delete before continuing.', 'SEVERE')

    def __check_arg_maskmode(self):
        maskmode_list = self.sdbaseline_maskmode == 'list'
        maskmode_auto = self.sdbaseline_maskmode == 'auto'
        if not ( maskmode_list or maskmode_auto ):
            raise ValueError('Unsupported maskmode, ' + self.sdbaseline_maskmode)
        if maskmode_list:
            self.sdbaseline_thresh = 5.0
            self.sdbaseline_avg_limit = self.sdbaseline_minwidth = 4
            self.sdbaseline_edge = [0, 0]
        else:
            self.sdbaseline_chans = ''

    def __check_arg_blfunc(self):
        blfunc_poly = self.sdbaseline_blfunc == 'poly'
        blfunc_chebyshev = self.sdbaseline_blfunc == 'chebyshev'
        blfunc_cspline = self.sdbaseline_blfunc == 'cspline'
        blfunc_sinusoid = self.sdbaseline_blfunc == 'sinusoid'
        blfunc_variable = self.sdbaseline_blfunc == 'variable'
        if not (blfunc_poly or blfunc_chebyshev or blfunc_cspline or blfunc_sinusoid or blfunc_variable):
            raise ValueError('Unsupported blfunc, ' + self.sdbaseline_blfunc)
        
        if not (blfunc_poly or blfunc_chebyshev):
            self.sdbaseline_order = 5
        if not blfunc_cspline:
            self.sdbaseline_npiece = 3
        if not blfunc_sinusoid:
            self.sdbaseline_applyfft = True
            self.sdbaseline_fftthresh = 3.0
            self.sdbaseline_addwn = [0]
            self.sdbaseline_rejwn = []
        if not blfunc_variable:
            self.sdbaseline_blparam = ''
        else:
            self.sdbaseline_clipniter = 0
            self.sdbaseline_clipthresh = 3.0
            if not os.path.exists(self.sdbaseline_blparam):
                raise ValueError("input file '%s' does not exists" % self.blparam)

    def __check_arg_dirkernel(self):
        self.dir_none = self.imsmooth_kernel == 'none'
        self.dir_ikernel = self.imsmooth_kernel == 'image'
        self.dir_bkernel = self.imsmooth_kernel == 'boxcar'
        self.dir_gkernel = self.imsmooth_kernel == 'gaussian'
        if not ( self.dir_none or self.dir_gkernel or self.dir_bkernel or self.dir_ikernel ):
            raise ValueError('Unsupported direction smoothing kernel, ' + self.imsmooth_kernel)
        
        if self.dir_none:
            self.imsmooth_output = self.imagename
        elif self.dir_ikernel:
            self.imsmooth_major = self.imsmooth_minor = self.imsmooth_pa = ''
            if self.imsmooth_kimage != '' and not os.path.exists(self.imsmooth_kimage):
                raise ValueError(f'Error: file {self.imsmooth_kimage} is not found.', 'SEVERE')
        else:
            self.imsmooth_kimage = ''
            self.imsmooth_scale = -1.0

    def __check_arg_spkernel(self):
        self.sp_none = self.sdsmooth_kernel == 'none'
        self.sp_bkernel = self.sdsmooth_kernel == 'boxcar'
        self.sp_gkernel = self.sdsmooth_kernel == 'gaussian'
        if not ( self.sp_none or self.sp_bkernel or self.sp_gkernel ):
            raise ValueError('Unsupported spectral smoothing kernel, ' + self.sdsmooth_kernel)
        
        if self.sp_none:
            self.sdsmooth_output = self.temporary_vis
            self.sdsmooth_kwidth = 5
    
    def imsmooth_enable(self):
        return not self.dir_none
    
    def sdsmooth_enable(self):
        return not self.sp_none

    def convert_sdbaselining_dict(self):
        """convert properties to a dictionary for sdbaseline process"""
        return { 'blfunc' : self.sdbaseline_blfunc,
                 'datacolumn' : self.datacolumn.lower(),
                 'outfile' : self.sdbaseline_output,
                 'bloutput' : ','.join(self.sdbaseline_bloutput),
                 'dosubtract' : self.sdbaseline_dosubtract,
                 'spw' : self.sdbaseline_spw,
                 'pol' : self.sdbaseline_pol,
                 'linefinding' : (self.sdbaseline_maskmode=='auto'),
                 'threshold' : self.sdbaseline_thresh,
                 'avg_limit' : self.sdbaseline_avg_limit,
                 'minwidth' : self.sdbaseline_minwidth,
                 'edge' : self.sdbaseline_edge,
                 'order' : self.sdbaseline_order,
                 'npiece' : self.sdbaseline_npiece,
                 'applyfft' : self.sdbaseline_applyfft,
                 'fftmethod' : self.sdbaseline_fftmethod,
                 'fftthresh' : self.sdbaseline_fftthresh,
                 'addwn' : self.sdbaseline_addwn,
                 'rejwn' : self.sdbaseline_rejwn,
                 'clip_threshold_sigma' : self.sdbaseline_clipthresh,
                 'num_fitting_max' : self.sdbaseline_clipniter+1,
                 'blparam' : self.sdbaseline_blparam,
                 'verbose' : False,
                 'updateweight' : self.sdbaseline_updateweight,
                 'sigmavalue' : self.sdbaseline_sigmavalue }
            

@sdutil.sdtask_decorator
def imbaseline(imagename, linefile, output_cont, bloutput, maskmode, chans, thresh, avg_limit, minwidth, edge, blfunc, order, 
               npiece, applyfft, fftthresh, addwn, rejwn, blparam, clipniter, clipthresh, dirkernel, major, minor, pa, kimage, 
               scale, spkernel, kwidth):
    """
    THE MAIN METHOD OF IMBASELINE.
    """
    vals = ImBaselineVals(imagename, linefile, output_cont, bloutput, maskmode, chans, thresh, avg_limit, minwidth,
                         edge, blfunc, order, npiece, applyfft, fftthresh, addwn, rejwn, blparam, clipniter,
                         clipthresh, dirkernel, major, minor, pa, kimage, scale, spkernel, kwidth)

    prepare(vals)

    # imsmooth -> convert casaimage to MS -> sdsmooth
    # -> sdbaseline -> convert MS to image

    # imsmooth
    if vals.imsmooth_enable():
        print("start imsmooth")
        Imsmooth(vals).execute()
        print("end imsmooth")

    # casaimage -> MS
    Image2MSConverter(vals).convert()

    # sdsmooth
    if vals.sdsmooth_enable():
        print("start sdsmooth")
        Sdsmooth(vals).execute()
        print("end sdsmooth")

    # sdbaseline
    print("start sdbaseline")
    Sdbaseline(vals).execute()
    print("end sdbaseline")

    # MS -> casaimage
    MS2ImageConverter(vals).convert()

    # end

def prepare(vals: ImBaselineVals = None):
    ia.dohistory(False)


class Imsmooth():
    """
    imsmooth execution class
    This code is based on task_imsmooth.
    """

    def __init__(self, vals: ImBaselineVals = None):
        self.vals = vals

    def execute(self):

        try:
            ia.open(self.vals.imagename)
            mycsys = ia.coordsys()
            myrg = regionmanager()
            self.vals.reg = myrg.frombcs(csys=mycsys.torecord(), shape=ia.shape(), chans=self.vals.imsmooth_chans)
        finally:
            myrg.done()
            mycsys.done()
        
        outia = None
        try:
            if self.vals.dir_gkernel:
                outia = self.__imsmooth_gckernel()
            elif self.vals.dir_bkernel:
                outia = self.__imsmooth_bkernel()
            elif self.vals.dir_ikernel:
                outia = self.__imsmooth_ikernel()

        finally:
            ia.done()
            if outia: outia.done()

    def __imsmooth_ikernel(self):
        # image kernel for image smoothing
        return ia.convolve(
                        outfile=self.vals.imsmooth_output, kernel=self.vals.imsmooth_kimage, scale=self.vals.imsmooth_scale, 
                        region=self.vals.reg, mask=self.vals.imsmooth_mask, overwrite=self.vals.overwrite, stretch=self.vals.imsmooth_stretch 
                    )

    def __imsmooth_bkernel(self):
        # boxcar kernel for image smoothing
        if not self.vals.imsmooth_major or not self.vals.imsmooth_minor:
            raise ValueError("Both major and minor must be specified.")

        casalog.post( "ia.sepconvolve( axes=[0,1],"+\
                                "types=['boxcar','boxcar' ],"+\
                                "widths=[ "+str(self.vals.imsmooth_minor)+", "+str(self.vals.imsmooth_major)+" ],"+ \
                                "region="+str(self.vals.reg)+",outfile="+self.vals.imsmooth_output+" )",\
                                'DEBUG2' )
        return ia.sepconvolve(
                        axes=[0,1], types=['box','box'], widths=[ self.vals.imsmooth_minor, self.vals.imsmooth_major ],
                        region=self.vals.reg, outfile=self.vals.imsmooth_output, mask=self.vals.imsmooth_mask, 
                        overwrite=self.vals.overwrite, stretch=self.vals.imsmooth_stretch 
                    )

    def __imsmooth_gckernel(self):
        # gaussian kernel for image smoothing
        if not self.vals.imsmooth_major:
            raise ValueError("Major axis must be specified")
        if not self.vals.imsmooth_minor:
            raise ValueError("Minor axis must be specified")
        if not self.vals.imsmooth_pa:
            raise ValueError("Position angle must be specified")
            
        return ia.convolve2d(
                        axes=[0,1], region=self.vals.reg, major=self.vals.imsmooth_major,
                        minor=self.vals.imsmooth_minor, pa=self.vals.imsmooth_pa, outfile=self.vals.imsmooth_output,
                        mask=self.vals.imsmooth_mask, stretch=self.vals.imsmooth_stretch, targetres=self.vals.imsmooth_targetres,
                        beam=self.vals.imsmooth_beam, overwrite=self.vals.overwrite
                    )


class Image2MSConverter():
    """
    The class does convert a Casa image to a MeasurementSet.
    """

    def __init__(self, vals: ImBaselineVals = None):
        self.vals = vals

    def convert(self):
        self.__get_image_params_from_imsmooth_output()
        self.__create_empty_ms()
        self.__copy_image_array_to_ms()

    def __get_image_params_from_imsmooth_output(self):
        try:
            ia.open(self.vals.imsmooth_output)
            cs = ia.coordsys()
            self.vals.imshape = ia.shape()
            self.vals.diraxis = cs.findcoordinate('direction')['world']
            self.vals.spaxis = cs.findcoordinate('spectral')['world'] # 3 or 2
            self.vals.polaxis = cs.findcoordinate('stokes')['world']  # 2 or 3 or None
        finally:
            cs.done()
            ia.done()

        assert len(self.vals.diraxis) > 0

        self.vals.spaxisexist = len(self.vals.spaxis) == 1
        self.vals.polaxisexist = len(self.vals.polaxis) == 1

        self.vals.dirshape = self.vals.imshape[self.vals.diraxis]
        self.vals.imnrow = np.prod(self.vals.dirshape)

        if self.vals.spaxisexist:
            self.vals.imnchan = self.vals.imshape[self.vals.spaxis[0]]
        else:
            self.vals.imnchan = 1
        if self.vals.polaxisexist:
            self.vals.imnpol = self.vals.imshape[self.vals.polaxis[0]]
        else:
            self.vals.imnpol = 1
        print(f'image shape is {self.vals.imshape}, direciton {self.vals.dirshape} ({self.vals.imnrow} pixels), npol {self.vals.imnpol}, nchan {self.vals.imnchan}')

        # if imnchan is too few, say, <10, sdbaseline should abort
        if self.vals.imnchan < 10:
            print(f'nchan {self.vals.imnchan} is too few to perform baseline subtraction')
            return False
        return True

    def __create_empty_ms(self):

        self.__check_ms_path()
        self.__create_maintable()
        self.__create_antenna_table()
        self.__create_data_description_table()
        self.__create_feed_table()
        self.__create_field_table()
        self.__create_flag_cmd_table()
        self.__create_history_table()
        self.__create_observation_table()
        self.__create_pointing_table()
        self.__create_polarization_table()
        self.__create_processor_table()
        self.__create_source_table()
        self.__create_special_window_table()
        self.__create_state_table()

    def __create_maintable(self):

        tb = table()
        try:
            tb.create(self.vals.temporary_vis, EmptyMSBaseInformation.ms_desc, dminfo=EmptyMSBaseInformation.ms_dminfo)
            tb.putkeyword(keyword="MS_VERSION", value=2)
            nrow = tb.nrows()
            nrow_req = self.vals.imnrow * self.vals.imnpol
            nrow = tb.nrows()
            if nrow != nrow_req:
                tb.addrows(nrow_req)
            ddid = tb.getcol('DATA_DESC_ID')
            ddid[:] = 0
            tb.putcol('DATA_DESC_ID', ddid)
            dummy = np.zeros(nrow_req, dtype=int)
            tb.putcol('ANTENNA1', dummy)
            tb.putcol('ANTENNA2', dummy)
            tb.putcol('STATE_ID', dummy)
            time_list = 4304481539.999771 + np.arange(nrow_req) * 30.0
            tb.putcol('TIME', time_list)
            print(f'number of rows {nrow}, number of image pixels {self.vals.imnrow}, number of pols {self.vals.imnpol}, required rows {nrow_req}')
        finally:
            tb.close()

    def __create_antenna_table(self):
        self.__create_subtable("ANTENNA", EmptyMSBaseInformation.antenna_desc, EmptyMSBaseInformation.antenna_dminfo)
        with table_manager(os.path.join(self.vals.temporary_vis, 'ANTENNA'), nomodify=False) as tb:
            tb.addrows(1)

    def __create_data_description_table(self):
        self.__create_subtable("DATA_DESCRIPTION", EmptyMSBaseInformation.data_description_desc, EmptyMSBaseInformation.data_description_dminfo)
        with table_manager(os.path.join(self.vals.temporary_vis, 'DATA_DESCRIPTION'), nomodify=False) as tb:
            tb.addrows(1)
            tb.putcell('SPECTRAL_WINDOW_ID', 0, 0)
            tb.putcell('POLARIZATION_ID', 0, 0)

    def __create_feed_table(self):
        self.__create_subtable("FEED", EmptyMSBaseInformation.feed_desc, EmptyMSBaseInformation.feed_dminfo)
        with table_manager(os.path.join(self.vals.temporary_vis, 'FEED'), nomodify=False) as tb:
            tb.addrows(1)

    def __create_field_table(self):
        self.__create_subtable("FIELD", EmptyMSBaseInformation.field_desc, EmptyMSBaseInformation.field_dminfo)
        with table_manager(os.path.join(self.vals.temporary_vis, 'FIELD'), nomodify=False) as tb:
            tb.addrows(1)
            tb.putcell('DELAY_DIR', 0, np.zeros((2,1)))
            tb.putcell('PHASE_DIR', 0, np.zeros((2,1)))
            tb.putcell('REFERENCE_DIR', 0, np.zeros((2,1)))

    def __create_flag_cmd_table(self):
        self.__create_subtable("FLAG_CMD", EmptyMSBaseInformation.flag_cmd_desc, EmptyMSBaseInformation.flag_cmd_dminfo)

    def __create_history_table(self):
        self.__create_subtable("HISTORY", EmptyMSBaseInformation.history_desc, EmptyMSBaseInformation.history_dminfo)

    def __create_observation_table(self):
        self.__create_subtable("OBSERVATION", EmptyMSBaseInformation.observation_desc, EmptyMSBaseInformation.observation_dminfo)
        with table_manager(os.path.join(self.vals.temporary_vis, 'OBSERVATION'), nomodify=False) as tb:
            tb.addrows(1)

    def __create_pointing_table(self):
        self.__create_subtable("POINTING", EmptyMSBaseInformation.pointing_desc, EmptyMSBaseInformation.pointing_dminfo)

    def __create_polarization_table(self):
        self.__create_subtable("POLARIZATION", EmptyMSBaseInformation.polarization_desc, EmptyMSBaseInformation.polarization_dminfo)
        with table_manager(os.path.join(self.vals.temporary_vis, 'POLARIZATION'), nomodify=False) as tb:
            corr_type = np.ones(1, dtype=int)
            corr_product = np.ones(2, dtype=int).reshape((2, 1))
            if tb.nrows() == 0:
                tb.addrows(1)
            tb.putcell('NUM_CORR', 0, 1)
            tb.putcell('CORR_TYPE', 0, corr_type)
            tb.putcell('CORR_PRODUCT', 0, corr_product)

    def __create_processor_table(self):
        self.__create_subtable("PROCESSOR", EmptyMSBaseInformation.processor_desc, EmptyMSBaseInformation.processor_dminfo)

    def __create_source_table(self):
        self.__create_subtable("SOURCE", EmptyMSBaseInformation.source_desc, EmptyMSBaseInformation.source_dminfo)

    def __create_special_window_table(self):
        self.__create_subtable("SPECTRAL_WINDOW", EmptyMSBaseInformation.special_window_desc, EmptyMSBaseInformation.special_window_dminfo)
        with table_manager(os.path.join(self.vals.temporary_vis, 'SPECTRAL_WINDOW'), nomodify=False) as tb:
            cw = np.ones(self.vals.imnchan, dtype=float) * 1e6
            cf = 1e9 + np.arange(self.vals.imnchan, dtype=float) * 1e6

            if tb.nrows() == 0:
                tb.addrows(1)
            tb.putcell('NUM_CHAN', 0, self.vals.imnchan)
            tb.putcell('CHAN_FREQ', 0, cf)
            tb.putcell('CHAN_WIDTH', 0, cw)
            tb.putcell('RESOLUTION', 0, cw)
            tb.putcell('EFFECTIVE_BW', 0, cw)
            tb.putcell('TOTAL_BANDWIDTH', 0, cw.sum())

    def __create_state_table(self):
        self.__create_subtable("STATE", EmptyMSBaseInformation.state_desc, EmptyMSBaseInformation.state_dminfo)
        with table_manager(os.path.join(self.vals.temporary_vis, 'STATE'), nomodify=False) as tb:
            if tb.nrows() == 0:
                tb.addrows(1)
            tb.putcell('OBS_MODE', 0, 'OBSERVE_TARGET#ON_SOURCE_IMAGE_DOMAIN')

    def __check_ms_path(self, overWrite: bool = True):
        exists = os.path.exists(self.vals.temporary_vis)
        if overWrite and exists:
            shutil.rmtree(self.vals.temporary_vis)
        return exists

    def __create_subtable(self, subtable: str, desc: str, dminfo: str):
        tb = table()
        try:
            tb.create(f"{self.vals.temporary_vis}/{subtable}", desc, dminfo=dminfo)
            tb.open(self.vals.temporary_vis, nomodify=False)
            tb.putkeyword(subtable, f"Table: {self.vals.temporary_vis}/{subtable}")
        finally:
            tb.close()

    def __copy_image_array_to_ms(self):
        # get image array and mask from the image
        ia.open(self.vals.imsmooth_output)
        arr = ia.getchunk()
        msk = ia.getchunk(getmask=True)
        ia.done()

        # put image array slices to MS DATA column
        # axis indices for spatial, spectral and polarization axes
        xax, yax = self.vals.diraxis
        spax = self.vals.spaxis[0]
        if self.vals.polaxisexist:
            polax = self.vals.polaxis[0]
        else:
            arr = np.expand_dims(arr, axis=3)
            msk = np.expand_dims(msk, axis=3)
            polax = 3
        casalog.post(f'axis index: {xax} {yax} {polax} {spax}')

        # which data column to use
        with table_manager(self.vals.temporary_vis, nomodify=False) as tb:
            # also set FLAG, SIGMA, WEIGHT, and UVW
            # r2p is a mapping between MS row and image array slice
            xax, yax = self.vals.diraxis

            index_list = [0, 0, 0, 0]
            index_list[spax] = np.arange(self.vals.imnchan)
            r2p = []
            nx, ny = self.vals.dirshape
            irow = 0
            wgt = np.ones(1, dtype=float)
            uvw = np.zeros(3, dtype=float)
            for ix in range(nx):
                index_list[xax] = ix
                for iy in range(ny):
                    index_list[yax] = iy
                    for ip in range(self.vals.imnpol):
                        index_list[polax] = ip
                        slice = tuple(index_list)
                        subarr = arr[slice]
                        submsk = np.logical_not(msk[slice])
                        print(f'slice={slice}, subarr={subarr[200:-200]}')
                        r2p.append(slice)
                        tb.putcell(self.vals.datacolumn, irow, np.expand_dims(subarr, axis=0))
                        tb.putcell('FLAG', irow, np.expand_dims(submsk, axis=0))
                        tb.putcell('SIGMA', irow, wgt)
                        tb.putcell('WEIGHT', irow, wgt)
                        tb.putcell('UVW', irow, uvw)
                        irow += 1


class Sdsmooth():
    """
    sdsmooth execution class
    This code is based on task_sdsmooth.
    """
    def __init__(self, vals: ImBaselineVals = None):
        self.vals = vals

    def execute(self):
        sdms.open(self.vals.temporary_vis)
        sdms.set_selection(spw=self.vals.sdsmooth_spw, field=self.vals.sdsmooth_field, antenna=self.vals.sdsmooth_antenna,
                            timerange=self.vals.sdsmooth_timerange, scan=self.vals.sdsmooth_scan, polarization=self.vals.sdsmooth_pol,
                            intent=self.vals.sdsmooth_intent, reindex=self.vals.sdsmooth_reindex)
        sdms.smooth(type=self.vals.sdsmooth_kernel, width=self.vals.sdsmooth_kwidth, outfile=self.vals.sdsmooth_output,
                    datacolumn=self.vals.datacolumn.lower())
        sdms.close()


class Sdbaseline():
    """
    sdbaseline execution class
    This code is based on task_sdbaseline.
    """
    def __init__(self, vals: ImBaselineVals = None):
        self.vals = vals
    
    def execute(self):

        try:
            self.__prepare_sdbaseline()

            self.__output_bloutput_text_header()
            
            if self.vals.sdbaseline_blfunc == 'variable':
                sorttab_info = self.__remove_sorted_table_keyword(self.vals.sdsmooth_output)

            selected_spw = sdutil.get_spwids( ms.msseltoindex(vis=self.vals.sdsmooth_output, spw=self.vals.sdbaseline_spw, 
                                                            field=self.vals.sdbaseline_field, baseline='', time='',
                                                            scan=self.vals.sdbaseline_scan) )
            sdms.open(self.vals.sdsmooth_output)
            sdms.set_selection(spw=selected_spw, field=self.vals.sdbaseline_field, antenna=self.vals.sdbaseline_antenna,
                                timerange=self.vals.sdbaseline_timerenge, scan=self.vals.sdbaseline_scan,
                                polarization=self.vals.sdbaseline_pol, intent=self.vals.sdbaseline_intent,
                                reindex=self.vals.sdbaseline_reindex)
            func, params = self.__prepare_for_baselining(sdms, self.vals.convert_sdbaselining_dict())
            func(**params)
            sdms.close()
            
            if self.vals.sdbaseline_blfunc == 'variable':
                self.__restore_sorted_table_keyword(self.vals.sdsmooth_output, sorttab_info)

        except Exception:
            raise

    def __prepare_sdbaseline(self):
        blparam_file = self.vals.sdsmooth_output + '_blparam.txt'
        if os.path.exists(blparam_file):
            self.__remove_data(blparam_file)  # CAS-11781
        
        if self.vals.sdbaseline_spw == '': self.vals.sdbaseline_spw = '*'

        if self.vals.sdbaseline_blfunc == 'sinusoid':
            self.vals.sdbaseline_addwn = sdutil.parse_wavenumber_param(self.vals.sdbaseline_addwn)
            self.vals.sdbaseline_rejwn = sdutil.parse_wavenumber_param(self.vals.sdbaseline_rejwn)
        
        self.vals.sdbaseline_bloutput = self.__prepare_for_blformat_bloutput(self.vals.sdsmooth_output, self.vals.sdbaseline_bloutput, True)

    def __remove_data(self, filename):
        if os.path.exists(filename):
            if os.path.isdir(filename):
                shutil.rmtree(filename)
            elif os.path.isfile(filename):
                os.remove(filename)
            else:
                # could be a symlink
                os.remove(filename)

    def __remove_sorted_table_keyword(self, infile):
        res = {'is_sorttab': False, 'sorttab_keywd': '', 'sorttab_name': ''}
        with table_manager(infile, nomodify=False) as tb:
            try:
                sorttab_keywd = 'SORTED_TABLE'
                if sorttab_keywd in tb.keywordnames():
                    res['is_sorttab'] = True
                    res['sorttab_keywd'] = sorttab_keywd
                    res['sorttab_name'] = tb.getkeyword(sorttab_keywd)
                    tb.removekeyword(sorttab_keywd)
            except Exception:
                raise

        return res

    def __restore_sorted_table_keyword(self, infile, sorttab_info):
        if sorttab_info['is_sorttab'] and sorttab_info['sorttab_name'] != '':
            with table_manager(infile, nomodify=False) as tb:
                try:
                    tb.putkeyword(sorttab_info['sorttab_keywd'],
                                sorttab_info['sorttab_name'])
                except Exception:
                    raise

    def __prepare_for_baselining(self, sdms, keywords):
        params = {}
        funcname = 'subtract_baseline'

        blfunc = keywords['blfunc']
        keys = ['datacolumn', 'outfile', 'bloutput', 'dosubtract', 'spw', 
                'updateweight', 'sigmavalue']
        if blfunc in ['poly', 'chebyshev']:
            keys += ['blfunc', 'order']
        elif blfunc == 'cspline':
            keys += ['npiece']
            funcname += ('_' + blfunc)
        elif blfunc =='sinusoid':
            keys += ['applyfft', 'fftmethod', 'fftthresh', 'addwn', 'rejwn']
            funcname += ('_' + blfunc)
        elif blfunc == 'variable':
            keys += ['blparam', 'verbose']
            funcname += ('_' + blfunc)
        else:
            raise ValueError("Unsupported blfunc = %s" % blfunc)
        if blfunc != 'variable':
            keys += ['clip_threshold_sigma', 'num_fitting_max']
            keys += ['linefinding', 'threshold', 'avg_limit', 'minwidth', 'edge']
            
        for key in keys: params[key] = keywords[key]

        baseline_func = getattr(sdms, funcname)

        return baseline_func, params


    def __prepare_for_blformat_bloutput(self, infile, bloutput, overwrite):
        # force to string list
        blformat = ['text']
        bloutput = [bloutput]

        if self.__has_duplicate_nonnull_element_ex(bloutput, blformat):
            raise ValueError('duplicate elements in bloutput.')

        # fill bloutput items to be output, then rearrange them
        # in the order of blformat_item.
        return self.__normalise_bloutput(infile, blformat, bloutput, overwrite)

    def __has_duplicate_nonnull_element(self, in_list):
        #return True if in_list has duplicated elements other than ''
        duplicates = [key for key, val in Counter(in_list).items() if val > 1]
        len_duplicates = len(duplicates)
        
        if (len_duplicates >= 2):
            return True
        elif (len_duplicates == 1):
            return (duplicates[0] != '')
        else: #len_duplicates == 0
            return False


    def __has_duplicate_nonnull_element_ex(self, lst, base):
        # lst and base must have the same length.
        #
        # (1) extract elements from lst and make a new list
        #     if the element of base with the same index
        #     is not ''.
        # (2) check if the list made in (1) has duplicated
        #     elements other than ''.
        
        return self.__has_duplicate_nonnull_element([lst[i] for i in range(len(lst)) if base[i] != ''])

    def __normalise_bloutput(self, infile, blformat, bloutput, overwrite):
        normalised_bloutput = []
        for item in zip(['csv', 'text', 'table'], ['csv', 'txt', 'bltable']):
            normalised_bloutput.append(
                self.__get_normalised_name(infile, blformat, bloutput, item[0], item[1], overwrite))
        return normalised_bloutput

    def __get_normalised_name(self, infile, blformat, bloutput, name, ext, overwrite):
        fname = ''
        blformat_lower = [s.lower() for s in blformat]
        if (name in blformat_lower):
            fname = bloutput[blformat_lower.index(name)]
            if (fname == ''):
                fname = infile + '_blparam.' + ext
        if os.path.exists(fname):
            if overwrite:
                os.system('rm -rf %s' % fname)
            else:
                raise Exception(fname + ' exists.')
        return fname

    def __output_bloutput_text_header(self):
        if self.vals.sdbaseline_bloutput[1] == '': return
        
        with open(self.vals.sdbaseline_bloutput[1], 'w') as f:

            info = [['Source Table', self.vals.sdsmooth_output],
                    ['Output File', self.vals.sdbaseline_output if (self.vals.sdbaseline_output != '') else self.vals.sdsmooth_output],
                    ['Mask mode', self.vals.sdbaseline_maskmode]]
            separator = '#' * 60 + '\n'
            
            f.write(separator)
            for i in range(len(info)):
                f.write('%12s: %s\n' % tuple(info[i]))
            f.write(separator)
            f.write('\n')


class MS2ImageConverter():
    """
    The class does convert a MeasurementSet to a Casa image.
    """

    def __init__(self, vals: ImBaselineVals = None):
        self.vals = vals

    def convert(self):
        print("start imaging")
        self.__make_output_file() # mask data is copied in this method 

        self.__make_image_array()
        print("end arraying")
        if self.vals.output_cont:
            self.__output_cont_image()
        self.__output_image()
        print("end imaging")

    def __output_image(self):
        try:
            ia.open(self.vals.linefile)
            ia.putchunk(pixels=self.array, locking=True)

            try:
                param_names = imbaseline.__code__.co_varnames[:imbaseline.__code__.co_argcount]
                vars = locals( )
                param_vals = [vars[p] for p in param_names]
                write_image_history(ia, sys._getframe().f_code.co_name, param_names, param_vals, casalog)
            except Exception as instance:
                casalog.post("*** Error \'%s\' updating HISTORY" % (instance), 'WARN')
        finally:
            ia.done()

    def __output_cont_image(self):
        try:
            ia.open(self.vals.output_cont_file)
            ia.putchunk(pixels = ia.getchunk() - self.array, locking = True)
        finally:
            ia.done()

    def __make_image_array(self):
        nx, ny = self.vals.dirshape
        self.array = np.empty((nx,ny,self.vals.imnchan))
        if self.vals.polaxisexist:
            self.array = np.expand_dims(self.array, self.vals.polaxis[0])
        pos = 0
        with table_manager(self.vals.sdbaseline_output) as tb:
            ndim = tb.getcell(self.vals.datacolumn, 0).ndim
            for i in range(nx):
                for j in range(ny):
                    if ndim == 3:
                        self.array[i][j][0] =  tb.getcell(self.vals.datacolumn, pos)[0].real
                    else:
                        self.array[i][j] =  tb.getcell(self.vals.datacolumn, pos)[0].real
                    pos += 1

    def __make_output_file(self):
        shutil.copytree(self.vals.imagename, self.vals.linefile)
        if self.vals.output_cont:
            shutil.copytree(self.vals.imagename, self.vals.output_cont_file)



class EmptyMSBaseInformation:
    """
    The Parameters class for creating an empty MeasurementSet.

    This class contains dictionaries to create an empty MS using table.create(), and it has no method.
    Dictionaries have two types; desc(desctiption) and dminfo(data management infomation), 
    these are used as arguments of table.create(), and for a table creating, it needs a desc dict and a dminfo dict.
    so there are dicts of twice of table amount in a MeasurementSet.
    """

    ms_desc = {
        'ANTENNA1': {'comment': 'ID of first antenna in interferometer',
                     'dataManagerGroup': 'SSM',
                     'dataManagerType': 'StandardStMan',
                     'keywords': {},
                     'maxlen': 0,
                     'option': 0,
                     'valueType': 'int'},
        'ANTENNA2': {'comment': 'ID of second antenna in interferometer',
                     'dataManagerGroup': 'SSM',
                     'dataManagerType': 'StandardStMan',
                     'keywords': {},
                     'maxlen': 0,
                     'option': 0,
                     'valueType': 'int'},
        'ARRAY_ID': {'comment': 'ID of array or subarray',
                     'dataManagerGroup': 'ISMData',
                     'dataManagerType': 'IncrementalStMan',
                     'keywords': {},
                     'maxlen': 0,
                     'option': 0,
                     'valueType': 'int'},
        'DATA': {'comment': 'The data column',
                 'dataManagerGroup': 'TiledData',
                 'dataManagerType': 'TiledShapeStMan',
                 'keywords': {},
                 'maxlen': 0,
                 'ndim': 2,
                 'option': 0,
                 'valueType': 'complex'},
        'DATA_DESC_ID': {'comment': 'The data description table index',
                         'dataManagerGroup': 'SSM',
                         'dataManagerType': 'StandardStMan',
                         'keywords': {},
                         'maxlen': 0,
                         'option': 0,
                         'valueType': 'int'},
        'EXPOSURE': {'comment': 'The effective integration time',
                     'dataManagerGroup': 'ISMData',
                     'dataManagerType': 'IncrementalStMan',
                     'keywords': {'QuantumUnits': np.array(['s'], dtype='<U16')},
                     'maxlen': 0,
                     'option': 0,
                     'valueType': 'double'},
        'FEED1': {'comment': 'The feed index for ANTENNA1',
                  'dataManagerGroup': 'ISMData',
                  'dataManagerType': 'IncrementalStMan',
                  'keywords': {},
                  'maxlen': 0,
                  'option': 0,
                  'valueType': 'int'},
        'FEED2': {'comment': 'The feed index for ANTENNA2',
                  'dataManagerGroup': 'ISMData',
                  'dataManagerType': 'IncrementalStMan',
                  'keywords': {},
                  'maxlen': 0,
                  'option': 0,
                  'valueType': 'int'},
        'FIELD_ID': {'comment': 'Unique id for this pointing',
                     'dataManagerGroup': 'ISMData',
                     'dataManagerType': 'IncrementalStMan',
                     'keywords': {},
                     'maxlen': 0,
                     'option': 0,
                     'valueType': 'int'},
        'FLAG': {'comment': 'The data flags, array of bools with same shape as data',
                 'dataManagerGroup': 'TiledFlag',
                 'dataManagerType': 'TiledShapeStMan',
                 'keywords': {},
                 'maxlen': 0,
                 'ndim': 2,
                 'option': 0,
                 'valueType': 'boolean'},
        'FLAG_CATEGORY': {'comment': 'The flag category, NUM_CAT flags for each datum',
                          'dataManagerGroup': 'TiledFlagCategory',
                          'dataManagerType': 'TiledShapeStMan',
                          'keywords': {'CATEGORY': np.array(['FLAG_CMD', 'ORIGINAL', 'USER'], dtype='<U16')},
                          'maxlen': 0,
                          'ndim': 3,
                          'option': 0,
                          'valueType': 'boolean'},
        'FLAG_ROW': {'comment': 'Row flag - flag all data in this row if True',
                     'dataManagerGroup': 'ISMData',
                     'dataManagerType': 'IncrementalStMan',
                     'keywords': {},
                     'maxlen': 0,
                     'option': 0,
                     'valueType': 'boolean'},
        'INTERVAL': {'comment': 'The sampling interval',
                     'dataManagerGroup': 'ISMData',
                     'dataManagerType': 'IncrementalStMan',
                     'keywords': {'QuantumUnits': np.array(['s'], dtype='<U16')},
                     'maxlen': 0,
                     'option': 0,
                     'valueType': 'double'},
        'OBSERVATION_ID': {'comment': 'ID for this observation, index in OBSERVATION table',
                           'dataManagerGroup': 'ISMData',
                           'dataManagerType': 'IncrementalStMan',
                           'keywords': {},
                           'maxlen': 0,
                           'option': 0,
                           'valueType': 'int'},
        'PROCESSOR_ID': {'comment': 'Id for backend processor, index in PROCESSOR table',
                         'dataManagerGroup': 'ISMData',
                         'dataManagerType': 'IncrementalStMan',
                         'keywords': {},
                         'maxlen': 0,
                         'option': 0,
                         'valueType': 'int'},
        'SCAN_NUMBER': {'comment': 'Sequential scan number from on-line system',
                        'dataManagerGroup': 'ISMData',
                        'dataManagerType': 'IncrementalStMan',
                        'keywords': {},
                        'maxlen': 0,
                        'option': 0,
                        'valueType': 'int'},
        'SIGMA': {'comment': 'Estimated rms noise for channel with unity bandpass response',
                  'dataManagerGroup': 'TiledSigma',
                  'dataManagerType': 'TiledShapeStMan',
                  'keywords': {},
                  'maxlen': 0,
                  'ndim': 1,
                  'option': 0,
                  'valueType': 'float'},
        'STATE_ID': {'comment': 'ID for this observing state',
                     'dataManagerGroup': 'ISMData',
                     'dataManagerType': 'IncrementalStMan',
                     'keywords': {},
                     'maxlen': 0,
                     'option': 0,
                     'valueType': 'int'},
        'TIME': {'comment': 'Modified Julian Day',
                 'dataManagerGroup': 'ISMData',
                 'dataManagerType': 'IncrementalStMan',
                 'keywords': {'MEASINFO': {'Ref': 'TAI', 'type': 'epoch'},
                              'QuantumUnits': np.array(['s'], dtype='<U16')},
                 'maxlen': 0,
                 'option': 0,
                 'valueType': 'double'},
        'TIME_CENTROID': {'comment': 'Modified Julian Day',
                          'dataManagerGroup': 'ISMData',
                          'dataManagerType': 'IncrementalStMan',
                          'keywords': {'MEASINFO': {'Ref': 'TAI', 'type': 'epoch'},
                                       'QuantumUnits': np.array(['s'], dtype='<U16')},
                          'maxlen': 0,
                          'option': 0,
                          'valueType': 'double'},
        'UVW': {'comment': 'Vector with uvw coordinates (in meters)',
                'dataManagerGroup': 'TiledUVW',
                'dataManagerType': 'TiledColumnStMan',
                'keywords': {'MEASINFO': {'Ref': 'J2000', 'type': 'uvw'},
                             'QuantumUnits': np.array(['m', 'm', 'm'], dtype='<U16')},
                'maxlen': 0,
                'ndim': 1,
                'option': 5,
                'shape': np.array([3]),
                'valueType': 'double'},
        'WEIGHT': {'comment': 'Weight for each polarization spectrum',
                   'dataManagerGroup': 'TiledWgt',
                   'dataManagerType': 'TiledShapeStMan',
                   'keywords': {},
                   'maxlen': 0,
                   'ndim': 1,
                   'option': 0,
                   'valueType': 'float'},
        'WEIGHT_SPECTRUM': {'comment': 'Weight for each data point',
                            'dataManagerGroup': 'TiledWgtSpectrum',
                            'dataManagerType': 'TiledShapeStMan',
                            'keywords': {},
                            'maxlen': 0,
                            'ndim': 2,
                            'option': 0,
                            'valueType': 'float'},
    }

    ms_dminfo = {
        '*1': {'COLUMNS': array(['ARRAY_ID', 'EXPOSURE', 'FEED1', 'FEED2', 'FIELD_ID', 'FLAG_ROW',
                                'INTERVAL', 'OBSERVATION_ID', 'PROCESSOR_ID', 'SCAN_NUMBER',
                                 'STATE_ID', 'TIME', 'TIME_CENTROID'], dtype='<U16'),
               'NAME': 'ISMData',
               'SEQNR': 0,
               'SPEC': {'BUCKETSIZE': 62456, 'MaxCacheSize': 1, 'PERSCACHESIZE': 1},
               'TYPE': 'IncrementalStMan'},
        '*2': {'COLUMNS': array(['ANTENNA1', 'ANTENNA2', 'DATA_DESC_ID'], dtype='<U16'),
               'NAME': 'SSM',
               'SEQNR': 1,
               'SPEC': {'BUCKETSIZE': 32768,
                        'IndexLength': 198,
                        'MaxCacheSize': 2,
                        'PERSCACHESIZE': 2},
               'TYPE': 'StandardStMan'},
        '*3': {'COLUMNS': array(['DATA'], dtype='<U16'),
               'NAME': 'TiledData',
               'SEQNR': 2,
               'SPEC': {'DEFAULTTILESHAPE': array([2,   63, 1040]),
                        'HYPERCUBES': {'*1': {'BucketSize': 1048320,
                                              'CellShape': array([2, 63]),
                                              'CubeShape': array([2,    63, 22653]),
                                              'ID': {},
                                              'TileShape': array([2,   63, 1040])}},
                        'IndexSize': 1,
                        'MAXIMUMCACHESIZE': 0,
                        'MaxCacheSize': 0,
                        'SEQNR': 2},
               'TYPE': 'TiledShapeStMan'},
        '*4': {'COLUMNS': array(['FLAG'], dtype='<U16'),
               'NAME': 'TiledFlag',
               'SEQNR': 3,
               'SPEC': {'DEFAULTTILESHAPE': array([2,   63, 1040]),
                        'HYPERCUBES': {'*1': {'BucketSize': 16380,
                                              'CellShape': array([2, 63]),
                                              'CubeShape': array([2,    63, 22653]),
                                              'ID': {},
                                              'TileShape': array([2,   63, 1040])}},
                        'IndexSize': 1,
                        'MAXIMUMCACHESIZE': 0,
                        'MaxCacheSize': 0,
                        'SEQNR': 3},
               'TYPE': 'TiledShapeStMan'},
        '*5': {'COLUMNS': array(['FLAG_CATEGORY'], dtype='<U16'),
               'NAME': 'TiledFlagCategory',
               'SEQNR': 4,
               'SPEC': {'DEFAULTTILESHAPE': array([2,   63,    1, 1040]),
                        'HYPERCUBES': {},
                        'IndexSize': 0,
                        'MAXIMUMCACHESIZE': 0,
                        'MaxCacheSize': 0,
                        'SEQNR': 4},
               'TYPE': 'TiledShapeStMan'},
        '*6': {'COLUMNS': array(['WEIGHT_SPECTRUM'], dtype='<U16'),
               'NAME': 'TiledWgtSpectrum',
               'SEQNR': 5,
               'SPEC': {'DEFAULTTILESHAPE': array([2,   63, 1040]),
                        'HYPERCUBES': {'*1': {'BucketSize': 524160,
                                              'CellShape': array([2, 63]),
                                              'CubeShape': array([2,    63, 22653]),
                                              'ID': {},
                                              'TileShape': array([2,   63, 1040])}},
                        'IndexSize': 1,
                        'MAXIMUMCACHESIZE': 0,
                        'MaxCacheSize': 0,
                        'SEQNR': 5},
               'TYPE': 'TiledShapeStMan'},
        '*7': {'COLUMNS': array(['UVW'], dtype='<U16'),
               'NAME': 'TiledUVW',
               'SEQNR': 6,
               'SPEC': {'DEFAULTTILESHAPE': array([3, 1024]),
                        'HYPERCUBES': {'*1': {'BucketSize': 24576,
                                              'CellShape': array([3]),
                                              'CubeShape': array([3, 22653]),
                                              'ID': {},
                                              'TileShape': array([3, 1024])}},
                        'MAXIMUMCACHESIZE': 0,
                        'MaxCacheSize': 0,
                        'SEQNR': 6},
               'TYPE': 'TiledColumnStMan'},
        '*8': {'COLUMNS': array(['WEIGHT'], dtype='<U16'),
               'NAME': 'TiledWgt',
               'SEQNR': 7,
               'SPEC': {'DEFAULTTILESHAPE': array([2, 1040]),
                        'HYPERCUBES': {'*1': {'BucketSize': 8320,
                                              'CellShape': array([2]),
                                              'CubeShape': array([2, 22653]),
                                              'ID': {},
                                              'TileShape': array([2, 1040])}},
                        'IndexSize': 1,
                        'MAXIMUMCACHESIZE': 0,
                        'MaxCacheSize': 0,
                        'SEQNR': 7},
               'TYPE': 'TiledShapeStMan'},
        '*9': {'COLUMNS': array(['SIGMA'], dtype='<U16'),
               'NAME': 'TiledSigma',
               'SEQNR': 8,
               'SPEC': {'DEFAULTTILESHAPE': array([2, 1040]),
                        'HYPERCUBES': {'*1': {'BucketSize': 8320,
                                              'CellShape': array([2]),
                                              'CubeShape': array([2, 22653]),
                                              'ID': {},
                                              'TileShape': array([2, 1040])}},
                        'IndexSize': 1,
                        'MAXIMUMCACHESIZE': 0,
                        'MaxCacheSize': 0,
                        'SEQNR': 8},
               'TYPE': 'TiledShapeStMan'}
    }

    antenna_desc = {
        'DISH_DIAMETER': {'comment': 'Physical diameter of dish',
                          'dataManagerGroup': 'StandardStMan',
                          'dataManagerType': 'StandardStMan',
                          'keywords': {'QuantumUnits': array(['m'], dtype='<U16')},
                          'maxlen': 0,
                          'option': 0,
                          'valueType': 'double'},
        'FLAG_ROW': {'comment': 'Flag for this row',
                     'dataManagerGroup': 'StandardStMan',
                     'dataManagerType': 'StandardStMan',
                     'keywords': {},
                     'maxlen': 0,
                     'option': 0,
                     'valueType': 'boolean'},
        'MOUNT': {'comment': 'Mount type e.g. alt-az, equatorial, etc.',
                  'dataManagerGroup': 'StandardStMan',
                  'dataManagerType': 'StandardStMan',
                  'keywords': {},
                  'maxlen': 0,
                  'option': 0,
                  'valueType': 'string'},
        'NAME': {'comment': 'Antenna name, e.g. VLA22, CA03',
                 'dataManagerGroup': 'StandardStMan',
                 'dataManagerType': 'StandardStMan',
                 'keywords': {'ARRAY_NAME': 'VLA'},
                 'maxlen': 0,
                 'option': 0,
                 'valueType': 'string'},
        'OFFSET': {'comment': 'Axes offset of mount to FEED REFERENCE point',
                   'dataManagerGroup': 'StandardStMan',
                   'dataManagerType': 'StandardStMan',
                   'keywords': {'MEASINFO': {'Ref': 'ITRF', 'type': 'position'},
                                'QuantumUnits': array(['m', 'm', 'm'], dtype='<U16')},
                   'maxlen': 0,
                   'ndim': 1,
                   'option': 5,
                   'shape': array([3]),
                   'valueType': 'double'},
        'POSITION': {'comment': 'Antenna X,Y,Z phase reference position',
                     'dataManagerGroup': 'StandardStMan',
                     'dataManagerType': 'StandardStMan',
                     'keywords': {'ARRAY_POSITION': array([0., 0., 0.]),
                                  'MEASINFO': {'Ref': 'ITRF', 'type': 'position'},
                                  'QuantumUnits': array(['m', 'm', 'm'], dtype='<U16')},
                     'maxlen': 0,
                     'ndim': 1,
                     'option': 5,
                     'shape': array([3]),
                     'valueType': 'double'},
        'STATION': {'comment': 'Station (antenna pad) name',
                    'dataManagerGroup': 'StandardStMan',
                    'dataManagerType': 'StandardStMan',
                    'keywords': {},
                    'maxlen': 0,
                    'option': 0,
                    'valueType': 'string'},
        'TYPE': {'comment': 'Antenna type (e.g. SPACE-BASED)',
                 'dataManagerGroup': 'StandardStMan',
                 'dataManagerType': 'StandardStMan',
                 'keywords': {},
                 'maxlen': 0,
                 'option': 0,
                 'valueType': 'string'},
        '_define_hypercolumn_': {},
        '_keywords_': {'DEGPDY': 360.9856341442,
                       'GSTIA0': 3.5030897164680597,
                       'RDATE': 4304481539.999771,
                       'TIMSYS': 'TAI'},
        '_private_keywords_': {}
    }

    antenna_dminfo = {
        '*1': {'COLUMNS': array(['DISH_DIAMETER', 'FLAG_ROW', 'MOUNT', 'NAME', 'OFFSET', 'POSITION',
                                'STATION', 'TYPE'], dtype='<U16'),
               'NAME': 'StandardStMan',
               'SEQNR': 0,
               'SPEC': {'BUCKETSIZE': 3332,
                        'IndexLength': 126,
                        'MaxCacheSize': 2,
                        'PERSCACHESIZE': 2},
               'TYPE': 'StandardStMan'
               }
    }

    data_description_desc = {
        'FLAG_ROW': {'comment': 'Flag this row',
                     'dataManagerGroup': 'StandardStMan',
                     'dataManagerType': 'StandardStMan',
                     'keywords': {},
                     'maxlen': 0,
                     'option': 0,
                     'valueType': 'boolean'},
        'POLARIZATION_ID': {'comment': 'Pointer to polarization table',
                            'dataManagerGroup': 'StandardStMan',
                            'dataManagerType': 'StandardStMan',
                            'keywords': {},
                            'maxlen': 0,
                            'option': 0,
                            'valueType': 'int'},
        'SPECTRAL_WINDOW_ID': {'comment': 'Pointer to spectralwindow table',
                               'dataManagerGroup': 'StandardStMan',
                               'dataManagerType': 'StandardStMan',
                               'keywords': {},
                               'maxlen': 0,
                               'option': 0,
                               'valueType': 'int'},
        '_define_hypercolumn_': {},
        '_keywords_': {},
        '_private_keywords_': {}
    }

    data_description_dminfo = {
        '*1': {'COLUMNS': array(['FLAG_ROW', 'POLARIZATION_ID', 'SPECTRAL_WINDOW_ID'], dtype='<U19'),
               'NAME': 'StandardStMan',
               'SEQNR': 0,
               'SPEC': {'BUCKETSIZE': 260,
                        'IndexLength': 126,
                        'MaxCacheSize': 2,
                        'PERSCACHESIZE': 2},
               'TYPE': 'StandardStMan'
               }
    }

    feed_desc = {
        'ANTENNA_ID': {'comment': 'ID of antenna in this array',
                       'dataManagerGroup': 'StandardStMan',
                       'dataManagerType': 'StandardStMan',
                       'keywords': {},
                       'maxlen': 0,
                       'option': 0,
                       'valueType': 'int'},
        'BEAM_ID': {'comment': 'Id for BEAM model',
                    'dataManagerGroup': 'StandardStMan',
                    'dataManagerType': 'StandardStMan',
                    'keywords': {},
                    'maxlen': 0,
                    'option': 0,
                    'valueType': 'int'},
        'BEAM_OFFSET': {'comment': 'Beam position offset (on sky but in antennareference frame)',
                        'dataManagerGroup': 'StandardStMan',
                        'dataManagerType': 'StandardStMan',
                        'keywords': {'MEASINFO': {'Ref': 'J2000', 'type': 'direction'},
                                     'QuantumUnits': array(['rad', 'rad'], dtype='<U16')},
                        'maxlen': 0,
                        'ndim': 2,
                        'option': 0,
                        'valueType': 'double'},
        'FEED_ID': {'comment': 'Feed id',
                    'dataManagerGroup': 'StandardStMan',
                    'dataManagerType': 'StandardStMan',
                    'keywords': {},
                    'maxlen': 0,
                    'option': 0,
                    'valueType': 'int'},
        'INTERVAL': {'comment': 'Interval for which this set of parameters is accurate',
                     'dataManagerGroup': 'StandardStMan',
                     'dataManagerType': 'StandardStMan',
                     'keywords': {'QuantumUnits': array(['s'], dtype='<U16')},
                     'maxlen': 0,
                     'option': 0,
                     'valueType': 'double'},
        'NUM_RECEPTORS': {'comment': 'Number of receptors on this feed (probably 1 or 2)',
                          'dataManagerGroup': 'StandardStMan',
                          'dataManagerType': 'StandardStMan',
                          'keywords': {},
                          'maxlen': 0,
                          'option': 0,
                          'valueType': 'int'},
        'POLARIZATION_TYPE': {'comment': 'Type of polarization to which a given RECEPTOR responds',
                              'dataManagerGroup': 'StandardStMan',
                              'dataManagerType': 'StandardStMan',
                              'keywords': {},
                              'maxlen': 0,
                              'ndim': 1,
                              'option': 0,
                              'valueType': 'string'},
        'POL_RESPONSE': {'comment': 'D-matrix i.e. leakage between two receptors',
                         'dataManagerGroup': 'StandardStMan',
                         'dataManagerType': 'StandardStMan',
                         'keywords': {},
                         'maxlen': 0,
                         'ndim': 2,
                         'option': 0,
                         'valueType': 'complex'},
        'POSITION': {'comment': 'Position of feed relative to feed reference position',
                     'dataManagerGroup': 'StandardStMan',
                     'dataManagerType': 'StandardStMan',
                     'keywords': {'MEASINFO': {'Ref': 'ITRF', 'type': 'position'},
                                  'QuantumUnits': array(['m', 'm', 'm'], dtype='<U16')},
                     'maxlen': 0,
                     'ndim': 1,
                     'option': 5,
                     'shape': array([3]),
                     'valueType': 'double'},
        'RECEPTOR_ANGLE': {'comment': 'The reference angle for polarization',
                           'dataManagerGroup': 'StandardStMan',
                           'dataManagerType': 'StandardStMan',
                           'keywords': {'QuantumUnits': array(['rad'], dtype='<U16')},
                           'maxlen': 0,
                           'ndim': 1,
                           'option': 0,
                           'valueType': 'double'},
        'SPECTRAL_WINDOW_ID': {'comment': 'ID for this spectral window setup',
                               'dataManagerGroup': 'StandardStMan',
                               'dataManagerType': 'StandardStMan',
                               'keywords': {},
                               'maxlen': 0,
                               'option': 0,
                               'valueType': 'int'},
        'TIME': {'comment': 'Midpoint of time for which this set of parameters is accurate',
                 'dataManagerGroup': 'StandardStMan',
                 'dataManagerType': 'StandardStMan',
                 'keywords': {'MEASINFO': {'Ref': 'TAI', 'type': 'epoch'},
                              'QuantumUnits': array(['s'], dtype='<U16')},
                 'maxlen': 0,
                 'option': 0,
                 'valueType': 'double'},
        '_define_hypercolumn_': {},
        '_keywords_': {},
        '_private_keywords_': {}}

    feed_dminfo = {
        '*1': {'COLUMNS': array(['ANTENNA_ID', 'BEAM_ID', 'BEAM_OFFSET', 'FEED_ID', 'INTERVAL',
                                 'NUM_RECEPTORS', 'POLARIZATION_TYPE', 'POL_RESPONSE', 'POSITION',
                                 'RECEPTOR_ANGLE', 'SPECTRAL_WINDOW_ID', 'TIME'], dtype='<U19'),
               'NAME': 'StandardStMan',
               'SEQNR': 0,
               'SPEC': {'BUCKETSIZE': 3072,
                        'IndexLength': 126,
                        'MaxCacheSize': 2,
                        'PERSCACHESIZE': 2},
               'TYPE': 'StandardStMan'}}

    field_desc = {
        'CODE': {'comment': 'Special characteristics of field, e.g. Bandpass calibrator',
                 'dataManagerGroup': 'StandardStMan',
                 'dataManagerType': 'StandardStMan',
                 'keywords': {},
                 'maxlen': 0,
                 'option': 0,
                 'valueType': 'string'},
        'DELAY_DIR': {'comment': 'Direction of delay center (e.g. RA, DEC)as polynomial in time.',
                      'dataManagerGroup': 'StandardStMan',
                      'dataManagerType': 'StandardStMan',
                      'keywords': {'MEASINFO': {'Ref': 'J2000', 'type': 'direction'},
                                   'QuantumUnits': array(['rad', 'rad'], dtype='<U16')},
                      'maxlen': 0,
                      'ndim': 2,
                      'option': 0,
                      'valueType': 'double'},
        'FLAG_ROW': {'comment': 'Row Flag',
                     'dataManagerGroup': 'StandardStMan',
                     'dataManagerType': 'StandardStMan',
                     'keywords': {},
                     'maxlen': 0,
                     'option': 0,
                     'valueType': 'boolean'},
        'NAME': {'comment': 'Name of this field',
                 'dataManagerGroup': 'StandardStMan',
                 'dataManagerType': 'StandardStMan',
                 'keywords': {},
                 'maxlen': 0,
                 'option': 0,
                 'valueType': 'string'},
        'NUM_POLY': {'comment': 'Polynomial order of _DIR columns',
                     'dataManagerGroup': 'StandardStMan',
                     'dataManagerType': 'StandardStMan',
                     'keywords': {},
                     'maxlen': 0,
                     'option': 0,
                     'valueType': 'int'},
        'PHASE_DIR': {'comment': 'Direction of phase center (e.g. RA, DEC).',
                      'dataManagerGroup': 'StandardStMan',
                      'dataManagerType': 'StandardStMan',
                      'keywords': {'MEASINFO': {'Ref': 'J2000', 'type': 'direction'},
                                   'QuantumUnits': array(['rad', 'rad'], dtype='<U16')},
                      'maxlen': 0,
                      'ndim': 2,
                      'option': 0,
                      'valueType': 'double'},
        'REFERENCE_DIR': {'comment': 'Direction of REFERENCE center (e.g. RA, DEC).as polynomial in time.',
                          'dataManagerGroup': 'StandardStMan',
                          'dataManagerType': 'StandardStMan',
                          'keywords': {'MEASINFO': {'Ref': 'J2000', 'type': 'direction'},
                                       'QuantumUnits': array(['rad', 'rad'], dtype='<U16')},
                          'maxlen': 0,
                          'ndim': 2,
                          'option': 0,
                          'valueType': 'double'},
        'SOURCE_ID': {'comment': 'Source id',
                      'dataManagerGroup': 'StandardStMan',
                      'dataManagerType': 'StandardStMan',
                      'keywords': {},
                      'maxlen': 0,
                      'option': 0,
                      'valueType': 'int'},
        'TIME': {'comment': 'Time origin for direction and rate',
                 'dataManagerGroup': 'StandardStMan',
                 'dataManagerType': 'StandardStMan',
                 'keywords': {'MEASINFO': {'Ref': 'TAI', 'type': 'epoch'},
                              'QuantumUnits': array(['s'], dtype='<U16')},
                 'maxlen': 0,
                 'option': 0,
                 'valueType': 'double'},
        '_define_hypercolumn_': {},
        '_keywords_': {},
        '_private_keywords_': {}}

    field_dminfo = {
        '*1': {'COLUMNS': array(['CODE', 'DELAY_DIR', 'FLAG_ROW', 'NAME', 'NUM_POLY', 'PHASE_DIR',
                                 'REFERENCE_DIR', 'SOURCE_ID', 'TIME'], dtype='<U16'),
               'NAME': 'StandardStMan',
               'SEQNR': 0,
               'SPEC': {'BUCKETSIZE': 2052,
                        'IndexLength': 126,
                        'MaxCacheSize': 2,
                        'PERSCACHESIZE': 2},
               'TYPE': 'StandardStMan'}}

    flag_cmd_desc = {
        'APPLIED': {'comment': 'True if flag has been applied to main table',
                    'dataManagerGroup': 'StandardStMan',
                    'dataManagerType': 'StandardStMan',
                    'keywords': {},
                    'maxlen': 0,
                    'option': 0,
                    'valueType': 'boolean'},
        'COMMAND': {'comment': 'Flagging command',
                    'dataManagerGroup': 'StandardStMan',
                    'dataManagerType': 'StandardStMan',
                    'keywords': {},
                    'maxlen': 0,
                    'option': 0,
                    'valueType': 'string'},
        'INTERVAL': {'comment': 'Time interval for which this flag is valid',
                     'dataManagerGroup': 'StandardStMan',
                     'dataManagerType': 'StandardStMan',
                     'keywords': {'QuantumUnits': array(['s'], dtype='<U16')},
                     'maxlen': 0,
                     'option': 0,
                     'valueType': 'double'},
        'LEVEL': {'comment': 'Flag level - revision level ',
                  'dataManagerGroup': 'StandardStMan',
                  'dataManagerType': 'StandardStMan',
                  'keywords': {},
                  'maxlen': 0,
                  'option': 0,
                  'valueType': 'int'},
        'REASON': {'comment': 'Flag reason',
                   'dataManagerGroup': 'StandardStMan',
                   'dataManagerType': 'StandardStMan',
                   'keywords': {},
                   'maxlen': 0,
                   'option': 0,
                   'valueType': 'string'},
        'SEVERITY': {'comment': 'Severity code (0-10) ',
                     'dataManagerGroup': 'StandardStMan',
                     'dataManagerType': 'StandardStMan',
                     'keywords': {},
                     'maxlen': 0,
                     'option': 0,
                     'valueType': 'int'},
        'TIME': {'comment': 'Midpoint of interval for which this flag is valid',
                 'dataManagerGroup': 'StandardStMan',
                 'dataManagerType': 'StandardStMan',
                 'keywords': {'MEASINFO': {'Ref': 'TAI', 'type': 'epoch'},
                              'QuantumUnits': array(['s'], dtype='<U16')},
                 'maxlen': 0,
                 'option': 0,
                 'valueType': 'double'},
        'TYPE': {'comment': 'Type of flag (FLAG or UNFLAG)',
                 'dataManagerGroup': 'StandardStMan',
                 'dataManagerType': 'StandardStMan',
                 'keywords': {},
                 'maxlen': 0,
                 'option': 0,
                 'valueType': 'string'},
        '_define_hypercolumn_': {},
        '_keywords_': {},
        '_private_keywords_': {}}

    flag_cmd_dminfo = {
        '*1': {'COLUMNS': array(['APPLIED', 'COMMAND', 'INTERVAL', 'LEVEL', 'REASON', 'SEVERITY',
                                 'TIME', 'TYPE'], dtype='<U16'),
               'NAME': 'StandardStMan',
               'SEQNR': 0,
               'SPEC': {'BUCKETSIZE': 1924,
                        'IndexLength': 118,
                        'MaxCacheSize': 2,
                        'PERSCACHESIZE': 2},
               'TYPE': 'StandardStMan'}}

    history_desc = {
        'APPLICATION': {'comment': 'Application name',
                        'dataManagerGroup': 'StandardStMan',
                        'dataManagerType': 'StandardStMan',
                        'keywords': {},
                        'maxlen': 0,
                        'option': 0,
                        'valueType': 'string'},
        'APP_PARAMS': {'comment': 'Application parameters',
                       'dataManagerGroup': 'StandardStMan',
                       'dataManagerType': 'StandardStMan',
                       'keywords': {},
                       'maxlen': 0,
                       'ndim': 1,
                       'option': 0,
                       'valueType': 'string'},
        'CLI_COMMAND': {'comment': 'CLI command sequence',
                        'dataManagerGroup': 'StandardStMan',
                        'dataManagerType': 'StandardStMan',
                        'keywords': {},
                        'maxlen': 0,
                        'ndim': 1,
                        'option': 0,
                        'valueType': 'string'},
        'MESSAGE': {'comment': 'Log message',
                    'dataManagerGroup': 'StandardStMan',
                    'dataManagerType': 'StandardStMan',
                    'keywords': {},
                    'maxlen': 0,
                    'option': 0,
                    'valueType': 'string'},
        'OBJECT_ID': {'comment': 'Originating ObjectID',
                      'dataManagerGroup': 'StandardStMan',
                      'dataManagerType': 'StandardStMan',
                      'keywords': {},
                      'maxlen': 0,
                      'option': 0,
                      'valueType': 'int'},
        'OBSERVATION_ID': {'comment': 'Observation id (index in OBSERVATION table)',
                           'dataManagerGroup': 'StandardStMan',
                           'dataManagerType': 'StandardStMan',
                           'keywords': {},
                           'maxlen': 0,
                           'option': 0,
                           'valueType': 'int'},
        'ORIGIN': {'comment': '(Source code) origin from which message originated',
                   'dataManagerGroup': 'StandardStMan',
                   'dataManagerType': 'StandardStMan',
                   'keywords': {},
                   'maxlen': 0,
                   'option': 0,
                   'valueType': 'string'},
        'PRIORITY': {'comment': 'Message priority',
                     'dataManagerGroup': 'StandardStMan',
                     'dataManagerType': 'StandardStMan',
                     'keywords': {},
                     'maxlen': 0,
                     'option': 0,
                     'valueType': 'string'},
        'TIME': {'comment': 'Timestamp of message',
                 'dataManagerGroup': 'StandardStMan',
                 'dataManagerType': 'StandardStMan',
                 'keywords': {'MEASINFO': {'Ref': 'TAI', 'type': 'epoch'},
                              'QuantumUnits': array(['s'], dtype='<U16')},
                 'maxlen': 0,
                 'option': 0,
                 'valueType': 'double'},
        '_define_hypercolumn_': {},
        '_keywords_': {},
        '_private_keywords_': {}}

    history_dminfo = {
        '*1': {'COLUMNS': array(['APPLICATION', 'APP_PARAMS', 'CLI_COMMAND', 'MESSAGE', 'OBJECT_ID',
                                 'OBSERVATION_ID', 'ORIGIN', 'PRIORITY', 'TIME'], dtype='<U16'),
               'NAME': 'StandardStMan',
               'SEQNR': 0,
               'SPEC': {'BUCKETSIZE': 2816,
                        'IndexLength': 126,
                        'MaxCacheSize': 2,
                        'PERSCACHESIZE': 2},
               'TYPE': 'StandardStMan'}}

    observation_desc = {
        'FLAG_ROW': {'comment': 'Row flag',
                     'dataManagerGroup': 'StandardStMan',
                     'dataManagerType': 'StandardStMan',
                     'keywords': {},
                     'maxlen': 0,
                     'option': 0,
                     'valueType': 'boolean'},
        'LOG': {'comment': 'Observing log',
                'dataManagerGroup': 'StandardStMan',
                'dataManagerType': 'StandardStMan',
                'keywords': {},
                'maxlen': 0,
                'ndim': 1,
                'option': 0,
                'valueType': 'string'},
        'OBSERVER': {'comment': 'Name of observer(s)',
                     'dataManagerGroup': 'StandardStMan',
                     'dataManagerType': 'StandardStMan',
                     'keywords': {},
                     'maxlen': 0,
                     'option': 0,
                     'valueType': 'string'},
        'PROJECT': {'comment': 'Project identification string',
                    'dataManagerGroup': 'StandardStMan',
                    'dataManagerType': 'StandardStMan',
                    'keywords': {},
                    'maxlen': 0,
                    'option': 0,
                    'valueType': 'string'},
        'RELEASE_DATE': {'comment': 'Release date when data becomes public',
                         'dataManagerGroup': 'StandardStMan',
                         'dataManagerType': 'StandardStMan',
                         'keywords': {'MEASINFO': {'Ref': 'TAI', 'type': 'epoch'},
                                      'QuantumUnits': array(['s'], dtype='<U16')},
                         'maxlen': 0,
                         'option': 0,
                         'valueType': 'double'},
        'SCHEDULE': {'comment': 'Observing schedule',
                     'dataManagerGroup': 'StandardStMan',
                     'dataManagerType': 'StandardStMan',
                     'keywords': {},
                     'maxlen': 0,
                     'ndim': 1,
                     'option': 0,
                     'valueType': 'string'},
        'SCHEDULE_TYPE': {'comment': 'Observing schedule type',
                          'dataManagerGroup': 'StandardStMan',
                          'dataManagerType': 'StandardStMan',
                          'keywords': {},
                          'maxlen': 0,
                          'option': 0,
                          'valueType': 'string'},
        'TELESCOPE_NAME': {'comment': 'Telescope Name (e.g. WSRT, VLBA)',
                           'dataManagerGroup': 'StandardStMan',
                           'dataManagerType': 'StandardStMan',
                           'keywords': {},
                           'maxlen': 0,
                           'option': 0,
                           'valueType': 'string'},
        'TIME_RANGE': {'comment': 'Start and end of observation',
                       'dataManagerGroup': 'StandardStMan',
                       'dataManagerType': 'StandardStMan',
                       'keywords': {'MEASINFO': {'Ref': 'TAI', 'type': 'epoch'},
                                    'QuantumUnits': array(['s'], dtype='<U16')},
                       'maxlen': 0,
                       'ndim': 1,
                       'option': 5,
                       'shape': array([2]),
                       'valueType': 'double'},
        '_define_hypercolumn_': {},
        '_keywords_': {},
        '_private_keywords_': {}}

    observation_dminfo = {
        '*1': {'COLUMNS': array(['FLAG_ROW', 'LOG', 'OBSERVER', 'PROJECT', 'RELEASE_DATE',
                                 'SCHEDULE', 'SCHEDULE_TYPE', 'TELESCOPE_NAME', 'TIME_RANGE'],
                                dtype='<U16'),
               'NAME': 'StandardStMan',
               'SEQNR': 0,
               'SPEC': {'BUCKETSIZE': 3076,
                        'IndexLength': 126,
                        'MaxCacheSize': 2,
                        'PERSCACHESIZE': 2},
               'TYPE': 'StandardStMan'}}

    pointing_desc = {
        'ANTENNA_ID': {'comment': 'Antenna Id',
                       'dataManagerGroup': 'SSMPointing',
                       'dataManagerType': 'StandardStMan',
                       'keywords': {},
                       'maxlen': 0,
                       'option': 0,
                       'valueType': 'int'},
        'DIRECTION': {'comment': 'Antenna pointing direction as polynomial in time',
                      'dataManagerGroup': 'ISMPointing',
                      'dataManagerType': 'IncrementalStMan',
                      'keywords': {'MEASINFO': {'Ref': 'J2000', 'type': 'direction'},
                                   'QuantumUnits': array(['rad', 'rad'], dtype='<U16')},
                      'maxlen': 0,
                      'ndim': 2,
                      'option': 0,
                      'valueType': 'double'},
        'INTERVAL': {'comment': 'Time interval',
                     'dataManagerGroup': 'ISMPointing',
                     'dataManagerType': 'IncrementalStMan',
                     'keywords': {'QuantumUnits': array(['s'], dtype='<U16')},
                     'maxlen': 0,
                     'option': 0,
                     'valueType': 'double'},
        'NAME': {'comment': 'Pointing position name',
                 'dataManagerGroup': 'ISMPointing',
                 'dataManagerType': 'IncrementalStMan',
                 'keywords': {},
                 'maxlen': 0,
                 'option': 0,
                 'valueType': 'string'},
        'NUM_POLY': {'comment': 'Series order',
                     'dataManagerGroup': 'ISMPointing',
                     'dataManagerType': 'IncrementalStMan',
                     'keywords': {},
                     'maxlen': 0,
                     'option': 0,
                     'valueType': 'int'},
        'TARGET': {'comment': 'target direction as polynomial in time',
                   'dataManagerGroup': 'ISMPointing',
                   'dataManagerType': 'IncrementalStMan',
                   'keywords': {'MEASINFO': {'Ref': 'J2000', 'type': 'direction'},
                                'QuantumUnits': array(['rad', 'rad'], dtype='<U16')},
                   'maxlen': 0,
                   'ndim': -1,
                   'option': 0,
                   'valueType': 'double'},
        'TIME': {'comment': 'Time interval midpoint',
                 'dataManagerGroup': 'ISMPointing',
                 'dataManagerType': 'IncrementalStMan',
                 'keywords': {'MEASINFO': {'Ref': 'TAI', 'type': 'epoch'},
                              'QuantumUnits': array(['s'], dtype='<U16')},
                 'maxlen': 0,
                 'option': 0,
                 'valueType': 'double'},
        'TIME_ORIGIN': {'comment': 'Time origin for direction',
                        'dataManagerGroup': 'ISMPointing',
                        'dataManagerType': 'IncrementalStMan',
                        'keywords': {'MEASINFO': {'Ref': 'TAI', 'type': 'epoch'},
                                     'QuantumUnits': array(['s'], dtype='<U16')},
                        'maxlen': 0,
                        'option': 0,
                        'valueType': 'double'},
        'TRACKING': {'comment': 'Tracking flag - True if on position',
                     'dataManagerGroup': 'ISMPointing',
                     'dataManagerType': 'IncrementalStMan',
                     'keywords': {},
                     'maxlen': 0,
                     'option': 0,
                     'valueType': 'boolean'},
        '_define_hypercolumn_': {},
        '_keywords_': {},
        '_private_keywords_': {}}

    pointing_dminfo = {
        '*1': {'COLUMNS': array(['DIRECTION', 'INTERVAL', 'NAME', 'NUM_POLY', 'TARGET', 'TIME',
                                 'TIME_ORIGIN', 'TRACKING'], dtype='<U16'),
               'NAME': 'ISMPointing',
               'SEQNR': 0,
               'SPEC': {'BUCKETSIZE': 38436, 'MaxCacheSize': 1, 'PERSCACHESIZE': 1},
               'TYPE': 'IncrementalStMan'},
        '*2': {'COLUMNS': array(['ANTENNA_ID'], dtype='<U16'),
               'NAME': 'SSMPointing',
               'SEQNR': 1,
               'SPEC': {'BUCKETSIZE': 32768,
                        'IndexLength': 118,
                        'MaxCacheSize': 2,
                        'PERSCACHESIZE': 2},
               'TYPE': 'StandardStMan'}}

    polarization_desc = {
        'CORR_PRODUCT': {'comment': 'Indices describing receptors of feed going into correlation',
                         'dataManagerGroup': 'StandardStMan',
                         'dataManagerType': 'StandardStMan',
                         'keywords': {},
                         'maxlen': 0,
                         'ndim': 2,
                         'option': 0,
                         'valueType': 'int'},
        'CORR_TYPE': {'comment': 'The polarization type for each correlation product, as a Stokes enum.',
                      'dataManagerGroup': 'StandardStMan',
                      'dataManagerType': 'StandardStMan',
                      'keywords': {},
                      'maxlen': 0,
                      'ndim': 1,
                      'option': 0,
                      'valueType': 'int'},
        'FLAG_ROW': {'comment': 'Row flag',
                     'dataManagerGroup': 'StandardStMan',
                     'dataManagerType': 'StandardStMan',
                     'keywords': {},
                     'maxlen': 0,
                     'option': 0,
                     'valueType': 'boolean'},
        'NUM_CORR': {'comment': 'Number of correlation products',
                     'dataManagerGroup': 'StandardStMan',
                     'dataManagerType': 'StandardStMan',
                     'keywords': {},
                     'maxlen': 0,
                     'option': 0,
                     'valueType': 'int'},
        '_define_hypercolumn_': {},
        '_keywords_': {},
        '_private_keywords_': {}}

    polarization_dminfo = {
        '*1': {'COLUMNS': array(['CORR_PRODUCT', 'CORR_TYPE', 'FLAG_ROW', 'NUM_CORR'], dtype='<U16'),
               'NAME': 'StandardStMan',
               'SEQNR': 0,
               'SPEC': {'BUCKETSIZE': 644,
                        'IndexLength': 126,
                        'MaxCacheSize': 2,
                        'PERSCACHESIZE': 2},
               'TYPE': 'StandardStMan'}}

    processor_desc = {
        'FLAG_ROW': {'comment': 'Row flag',
                     'dataManagerGroup': 'StandardStMan',
                     'dataManagerType': 'StandardStMan',
                     'keywords': {},
                     'maxlen': 0,
                     'option': 0,
                     'valueType': 'boolean'},
        'MODE_ID': {'comment': 'Processor mode id',
                    'dataManagerGroup': 'StandardStMan',
                    'dataManagerType': 'StandardStMan',
                    'keywords': {},
                    'maxlen': 0,
                    'option': 0,
                    'valueType': 'int'},
        'SUB_TYPE': {'comment': 'Processor sub type',
                     'dataManagerGroup': 'StandardStMan',
                     'dataManagerType': 'StandardStMan',
                     'keywords': {},
                     'maxlen': 0,
                     'option': 0,
                     'valueType': 'string'},
        'TYPE': {'comment': 'Processor type',
                 'dataManagerGroup': 'StandardStMan',
                 'dataManagerType': 'StandardStMan',
                 'keywords': {},
                 'maxlen': 0,
                 'option': 0,
                 'valueType': 'string'},
        'TYPE_ID': {'comment': 'Processor type id',
                    'dataManagerGroup': 'StandardStMan',
                    'dataManagerType': 'StandardStMan',
                    'keywords': {},
                    'maxlen': 0,
                    'option': 0,
                    'valueType': 'int'},
        '_define_hypercolumn_': {},
        '_keywords_': {},
        '_private_keywords_': {}}

    processor_dminfo = {
        '*1': {'COLUMNS': array(['FLAG_ROW', 'MODE_ID', 'SUB_TYPE', 'TYPE', 'TYPE_ID'], dtype='<U16'),
               'NAME': 'StandardStMan',
               'SEQNR': 0,
               'SPEC': {'BUCKETSIZE': 1028,
                        'IndexLength': 118,
                        'MaxCacheSize': 2,
                        'PERSCACHESIZE': 2},
               'TYPE': 'StandardStMan'}}

    source_desc = {
        'CALIBRATION_GROUP': {'comment': 'Number of grouping for calibration purpose.',
                                         'dataManagerGroup': 'StandardStMan',
                                         'dataManagerType': 'StandardStMan',
                                         'keywords': {},
                                         'maxlen': 0,
                                         'option': 0,
                                         'valueType': 'int'},
        'CODE': {'comment': 'Special characteristics of source, e.g. Bandpass calibrator',
                            'dataManagerGroup': 'StandardStMan',
                            'dataManagerType': 'StandardStMan',
                            'keywords': {},
                            'maxlen': 0,
                            'option': 0,
                            'valueType': 'string'},
        'DIRECTION': {'comment': 'Direction (e.g. RA, DEC).',
                                 'dataManagerGroup': 'StandardStMan',
                                 'dataManagerType': 'StandardStMan',
                                 'keywords': {'MEASINFO': {'Ref': 'J2000', 'type': 'direction'},
                                              'QuantumUnits': array(['rad', 'rad'], dtype='<U16')},
                                 'maxlen': 0,
                                 'ndim': 1,
                                 'option': 5,
                                 'shape': array([2]),
                                 'valueType': 'double'},
        'INTERVAL': {'comment': 'Interval of time for which this set of parameters is accurate',
                                'dataManagerGroup': 'StandardStMan',
                                'dataManagerType': 'StandardStMan',
                                'keywords': {'QuantumUnits': array(['s'], dtype='<U16')},
                                'maxlen': 0,
                                'option': 0,
                                'valueType': 'double'},
        'NAME': {'comment': 'Name of source as given during observations',
                            'dataManagerGroup': 'StandardStMan',
                            'dataManagerType': 'StandardStMan',
                            'keywords': {},
                            'maxlen': 0,
                            'option': 0,
                            'valueType': 'string'},
        'NUM_LINES': {'comment': 'Number of spectral lines',
                                 'dataManagerGroup': 'StandardStMan',
                                 'dataManagerType': 'StandardStMan',
                                 'keywords': {},
                                 'maxlen': 0,
                                 'option': 0,
                                 'valueType': 'int'},
        'POSITION': {'comment': 'Position (e.g. for solar system objects',
                                'dataManagerGroup': 'StandardStMan',
                                'dataManagerType': 'StandardStMan',
                                'keywords': {'MEASINFO': {'Ref': 'ITRF', 'type': 'position'},
                                             'QuantumUnits': array(['m', 'm', 'm'], dtype='<U16')},
                                'maxlen': 0,
                                'ndim': -1,
                                'option': 0,
                                'valueType': 'double'},
        'PROPER_MOTION': {'comment': 'Proper motion',
                                     'dataManagerGroup': 'StandardStMan',
                                     'dataManagerType': 'StandardStMan',
                                     'keywords': {'QuantumUnits': array(['rad/s'], dtype='<U16')},
                                     'maxlen': 0,
                                     'ndim': 1,
                                     'option': 5,
                                     'shape': array([2]),
                                     'valueType': 'double'},
        'REST_FREQUENCY': {'comment': 'Line rest frequency',
                                      'dataManagerGroup': 'StandardStMan',
                                      'dataManagerType': 'StandardStMan',
                                      'keywords': {'MEASINFO': {'Ref': 'LSRK', 'type': 'frequency'},
                                                   'QuantumUnits': array(['Hz'], dtype='<U16')},
                                      'maxlen': 0,
                                      'ndim': -1,
                                      'option': 0,
                                      'valueType': 'double'},
        'SOURCE_ID': {'comment': 'Source id',
                                 'dataManagerGroup': 'StandardStMan',
                                 'dataManagerType': 'StandardStMan',
                                 'keywords': {},
                                 'maxlen': 0,
                                 'option': 0,
                                 'valueType': 'int'},
        'SOURCE_MODEL': {'comment': 'Component Source Model',
                                    'dataManagerGroup': 'StandardStMan',
                                    'dataManagerType': 'StandardStMan',
                                    'keywords': {},
                                    'maxlen': 0,
                                    'option': 0,
                                    'valueType': 'record'},
        'SPECTRAL_WINDOW_ID': {'comment': 'ID for this spectral window setup',
                                          'dataManagerGroup': 'StandardStMan',
                                          'dataManagerType': 'StandardStMan',
                                          'keywords': {},
                                          'maxlen': 0,
                                          'option': 0,
                                          'valueType': 'int'},
        'SYSVEL': {'comment': 'Systemic velocity at reference',
                              'dataManagerGroup': 'StandardStMan',
                              'dataManagerType': 'StandardStMan',
                              'keywords': {'MEASINFO': {'Ref': 'LSRK', 'type': 'radialvelocity'},
                                           'QuantumUnits': array(['m/s'], dtype='<U16')},
                              'maxlen': 0,
                              'ndim': -1,
                              'option': 0,
                              'valueType': 'double'},
        'TIME': {'comment': 'Midpoint of time for which this set of parameters is accurate.',
                            'dataManagerGroup': 'StandardStMan',
                            'dataManagerType': 'StandardStMan',
                            'keywords': {'MEASINFO': {'Ref': 'TAI', 'type': 'epoch'},
                                         'QuantumUnits': array(['s'], dtype='<U16')},
                            'maxlen': 0,
                            'option': 0,
                            'valueType': 'double'},
        'TRANSITION': {'comment': 'Line Transition name',
                                  'dataManagerGroup': 'StandardStMan',
                                  'dataManagerType': 'StandardStMan',
                                  'keywords': {},
                                  'maxlen': 0,
                                  'ndim': -1,
                                  'option': 0,
                                  'valueType': 'string'},
        '_define_hypercolumn_': {},
        '_keywords_': {},
        '_private_keywords_': {}}

    source_dminfo = {
        '*1': {'COLUMNS': array(['CALIBRATION_GROUP', 'CODE', 'DIRECTION', 'INTERVAL', 'NAME',
                                 'NUM_LINES', 'POSITION', 'PROPER_MOTION', 'REST_FREQUENCY',
                                              'SOURCE_ID', 'SOURCE_MODEL', 'SPECTRAL_WINDOW_ID', 'SYSVEL',
                                              'TIME', 'TRANSITION'], dtype='<U19'),
               'NAME': 'StandardStMan',
               'SEQNR': 0,
               'SPEC': {'BUCKETSIZE': 4224,
                        'IndexLength': 126,
                        'MaxCacheSize': 2,
                        'PERSCACHESIZE': 2},
               'TYPE': 'StandardStMan'}}

    special_window_desc = {
        'CHAN_FREQ': {'comment': 'Center frequencies for each channel in the data matrix',
                      'dataManagerGroup': 'StandardStMan',
                      'dataManagerType': 'StandardStMan',
                                         'keywords': {'MEASINFO': {'TabRefCodes': array([0,  1,  2,  3,  4,  5,  6,  7,  8, 64], dtype=uint64),
                                                                   'TabRefTypes': array(['REST', 'LSRK', 'LSRD', 'BARY', 'GEO', 'TOPO', 'GALACTO', 'LGROUP',
                                                                                        'CMB', 'Undefined'], dtype='<U16'),
                                                                   'VarRefCol': 'MEAS_FREQ_REF',
                                                                   'type': 'frequency'},
                                                      'QuantumUnits': array(['Hz'], dtype='<U16')},
                                         'maxlen': 0,
                                         'ndim': 1,
                                         'option': 0,
                                         'valueType': 'double'},
        'CHAN_WIDTH': {'comment': 'Channel width for each channel',
                       'dataManagerGroup': 'StandardStMan',
                       'dataManagerType': 'StandardStMan',
                                          'keywords': {'QuantumUnits': array(['Hz'], dtype='<U16')},
                                          'maxlen': 0,
                                          'ndim': 1,
                                          'option': 0,
                                          'valueType': 'double'},
        'EFFECTIVE_BW': {'comment': 'Effective noise bandwidth of each channel',
                         'dataManagerGroup': 'StandardStMan',
                         'dataManagerType': 'StandardStMan',
                                            'keywords': {'QuantumUnits': array(['Hz'], dtype='<U16')},
                                            'maxlen': 0,
                                            'ndim': 1,
                                            'option': 0,
                                            'valueType': 'double'},
        'FLAG_ROW': {'comment': 'Row flag',
                     'dataManagerGroup': 'StandardStMan',
                     'dataManagerType': 'StandardStMan',
                                        'keywords': {},
                                        'maxlen': 0,
                                        'option': 0,
                                        'valueType': 'boolean'},
        'FREQ_GROUP': {'comment': 'Frequency group',
                       'dataManagerGroup': 'StandardStMan',
                       'dataManagerType': 'StandardStMan',
                                          'keywords': {},
                                          'maxlen': 0,
                                          'option': 0,
                                          'valueType': 'int'},
        'FREQ_GROUP_NAME': {'comment': 'Frequency group name',
                            'dataManagerGroup': 'StandardStMan',
                            'dataManagerType': 'StandardStMan',
                                               'keywords': {},
                                               'maxlen': 0,
                                               'option': 0,
                                               'valueType': 'string'},
        'IF_CONV_CHAIN': {'comment': 'The IF conversion chain number',
                          'dataManagerGroup': 'StandardStMan',
                          'dataManagerType': 'StandardStMan',
                                             'keywords': {},
                                             'maxlen': 0,
                                             'option': 0,
                                             'valueType': 'int'},
        'MEAS_FREQ_REF': {'comment': 'Frequency Measure reference',
                          'dataManagerGroup': 'StandardStMan',
                          'dataManagerType': 'StandardStMan',
                                             'keywords': {},
                                             'maxlen': 0,
                                             'option': 0,
                                             'valueType': 'int'},
        'NAME': {'comment': 'Spectral window name',
                 'dataManagerGroup': 'StandardStMan',
                 'dataManagerType': 'StandardStMan',
                                    'keywords': {},
                                    'maxlen': 0,
                                    'option': 0,
                                    'valueType': 'string'},
        'NET_SIDEBAND': {'comment': 'Net sideband',
                         'dataManagerGroup': 'StandardStMan',
                         'dataManagerType': 'StandardStMan',
                                            'keywords': {},
                                            'maxlen': 0,
                                            'option': 0,
                                            'valueType': 'int'},
        'NUM_CHAN': {'comment': 'Number of spectral channels',
                     'dataManagerGroup': 'StandardStMan',
                     'dataManagerType': 'StandardStMan',
                                        'keywords': {},
                                        'maxlen': 0,
                                        'option': 0,
                                        'valueType': 'int'},
        'REF_FREQUENCY': {'comment': 'The reference frequency',
                          'dataManagerGroup': 'StandardStMan',
                          'dataManagerType': 'StandardStMan',
                                             'keywords': {'MEASINFO': {'TabRefCodes': array([0,  1,  2,  3,  4,  5,  6,  7,  8, 64], dtype=uint64),
                                                                       'TabRefTypes': array(['REST', 'LSRK', 'LSRD', 'BARY', 'GEO', 'TOPO', 'GALACTO', 'LGROUP',
                                                                                             'CMB', 'Undefined'], dtype='<U16'),
                                                                       'VarRefCol': 'MEAS_FREQ_REF',
                                                                       'type': 'frequency'},
                                                          'QuantumUnits': array(['Hz'], dtype='<U16')},
                                             'maxlen': 0,
                                             'option': 0,
                                             'valueType': 'double'},
        'RESOLUTION': {'comment': 'The effective noise bandwidth for each channel',
                       'dataManagerGroup': 'StandardStMan',
                       'dataManagerType': 'StandardStMan',
                                          'keywords': {'QuantumUnits': array(['Hz'], dtype='<U16')},
                                          'maxlen': 0,
                                          'ndim': 1,
                                          'option': 0,
                                          'valueType': 'double'},
        'TOTAL_BANDWIDTH': {'comment': 'The total bandwidth for this window',
                            'dataManagerGroup': 'StandardStMan',
                            'dataManagerType': 'StandardStMan',
                                               'keywords': {'QuantumUnits': array(['Hz'], dtype='<U16')},
                                               'maxlen': 0,
                                               'option': 0,
                                               'valueType': 'double'},
        '_define_hypercolumn_': {},
        '_keywords_': {},
        '_private_keywords_': {}}

    special_window_dminfo = {
        '*1': {'COLUMNS': array(['CHAN_FREQ', 'CHAN_WIDTH', 'EFFECTIVE_BW', 'FLAG_ROW',
                                 'FREQ_GROUP', 'FREQ_GROUP_NAME', 'IF_CONV_CHAIN', 'MEAS_FREQ_REF',
                                 'NAME', 'NET_SIDEBAND', 'NUM_CHAN', 'REF_FREQUENCY', 'RESOLUTION',
                                 'TOTAL_BANDWIDTH'], dtype='<U16'),
               'NAME': 'StandardStMan',
               'SEQNR': 0,
               'SPEC': {'BUCKETSIZE': 2948,
                        'IndexLength': 126,
                        'MaxCacheSize': 2,
                        'PERSCACHESIZE': 2},
               'TYPE': 'StandardStMan'}}

    state_desc = {
        'CAL': {'comment': 'Noise calibration temperature',
                'dataManagerGroup': 'StandardStMan',
                'dataManagerType': 'StandardStMan',
                'keywords': {'QuantumUnits': array(['K'], dtype='<U16')},
                'maxlen': 0,
                'option': 0,
                'valueType': 'double'},
        'FLAG_ROW': {'comment': 'Row flag',
                     'dataManagerGroup': 'StandardStMan',
                     'dataManagerType': 'StandardStMan',
                     'keywords': {},
                     'maxlen': 0,
                     'option': 0,
                     'valueType': 'boolean'},
        'LOAD': {'comment': 'Load temperature',
                 'dataManagerGroup': 'StandardStMan',
                 'dataManagerType': 'StandardStMan',
                 'keywords': {'QuantumUnits': array(['K'], dtype='<U16')},
                 'maxlen': 0,
                 'option': 0,
                 'valueType': 'double'},
        'OBS_MODE': {'comment': 'Observing mode, e.g., OFF_SPECTRUM',
                     'dataManagerGroup': 'StandardStMan',
                     'dataManagerType': 'StandardStMan',
                     'keywords': {},
                     'maxlen': 0,
                     'option': 0,
                     'valueType': 'string'},
        'REF': {'comment': 'True for a reference observation',
                'dataManagerGroup': 'StandardStMan',
                'dataManagerType': 'StandardStMan',
                'keywords': {},
                'maxlen': 0,
                'option': 0,
                'valueType': 'boolean'},
        'SIG': {'comment': 'True for a source observation',
                'dataManagerGroup': 'StandardStMan',
                'dataManagerType': 'StandardStMan',
                'keywords': {},
                'maxlen': 0,
                'option': 0,
                'valueType': 'boolean'},
        'SUB_SCAN': {'comment': 'Sub scan number, relative to scan number',
                     'dataManagerGroup': 'StandardStMan',
                     'dataManagerType': 'StandardStMan',
                     'keywords': {},
                     'maxlen': 0,
                     'option': 0,
                     'valueType': 'int'},
        '_define_hypercolumn_': {},
        '_keywords_': {},
        '_private_keywords_': {}}

    state_dminfo = {
        '*1': {'COLUMNS': array(['CAL', 'FLAG_ROW', 'LOAD', 'OBS_MODE', 'REF', 'SIG', 'SUB_SCAN'],
                                dtype='<U16'),
               'NAME': 'StandardStMan',
               'SEQNR': 0,
               'SPEC': {'BUCKETSIZE': 1036,
                        'IndexLength': 118,
                        'MaxCacheSize': 2,
                        'PERSCACHESIZE': 2},
               'TYPE': 'StandardStMan'}}
