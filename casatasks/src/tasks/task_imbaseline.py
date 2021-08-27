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
from casatasks.private.ialib import write_image_history

ia = image()
qa = quanta()
ms = mstool()
sdms = singledishms()

class ImBaselineVals:

    def __init__(self, imagename, linefile, output_cont, bloutput, maskmode, chans, thresh, avg_limit, minwidth, edge,
                 blfunc, order, npiece, applyfft, fftthresh, addwn, rejwn, blparam, clipniter, clipthresh, dirkernel,
                 major, minor, pa, kimage, scale, spkernel, kwidth):
        self.imagename = imagename
        self.linefile = linefile
        self.output_cont = output_cont

        # sdbaseline
        self.sdbaseline_bloutput = bloutput
        self.sdbaseline_maskmode = maskmode     # list(default)/auto
        self.sdbaseline_chans = chans           # maskmode = list
        self.sdbaseline_thresh = thresh         # maskmode = auto
        self.sdbaseline_avg_limit = avg_limit   # maskmode = auto
        self.sdbaseline_minwidth = minwidth     # maskmode = auto
        self.sdbaseline_edge = edge             # maskmode = auto
        self.sdbaseline_blfunc = blfunc         # poly(default)/chebyshev/cspline/sinusoid/variable
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
        self.datacolumn = 'DATA'
        self.overwrite = True

    def check_args(self):
        self._check_linefile()
        self._check_dirkernel()
        self._check_spkernel()
        self._check_blparam()

    def _check_linefile(self):
        if len(self.linefile) > 0 and os.path.exists(self.linefile):
            raise ValueError(f'Error: file {self.linefile} already exists, please delete before continuing.', 'SEVERE')
        else:
            casalog.post("The linefile parameter is empty, consequently the"
                         + " spectral line image will NOT be\nsaved on disk.", 'WARN')

    def _check_dirkernel(self): # default is 'none' but its not found in imsmooth
        self.dir_ikernel = self.imsmooth_kernel == 'image'
        self.dir_bkernel = self.imsmooth_kernel == 'boxcar'
        self.dir_gkernel = self.imsmooth_kernel == 'gaussian'
        if not ( self.dir_gkernel or self.dir_bkernel or self.dir_ikernel ):
            raise ValueError('Unsupported direction smoothing kernel, ' + self.imsmooth_kernel)

    def _check_spkernel(self): # default is 'none' but its not found in sdsmooth
        self.sp_bkernel = self.sdsmooth_kernel == 'boxcar'
        self.sp_gkernel = self.sdsmooth_kernel == 'gaussian'
        if not ( self.sp_bkernel or self.sp_gkernel ):
            raise ValueError('Unsupported spectral smoothing kernel, ' + self.spkernel)
    
    def _check_blparam(self):
        if self.sdbaseline_blfunc == 'variable' and not os.path.exists(self.sdbaseline_blparam):
            raise ValueError("input file '%s' does not exists" % self.blparam)


@sdutil.sdtask_decorator
def imbaseline(
        imagename: str = '', linefile: str = '', output_cont: bool = False, bloutput: str = '', maskmode: str = 'list',
        chans: str = '', thresh: float = 5.0, avg_limit: int = 4, minwidth: int = 4, edge: List[int] = [0, 0],
        blfunc: str = 'poly', order: int = 5, npiece: int = 3, applyfft: bool = True, fftthresh: float = 3.0,
        addwn: List[int] = [0], rejwn: List[int] = [], blparam: str = '', clipniter: int = 0, clipthresh: float = 3.0,
        dirkernel: str = 'none', major: str = '', minor: str = '', pa: str = '', kimage: str = '',
        scale: float = -1.0, spkernel: str = 'none', kwidth: int = 5):
    vals = ImBaselineVals(imagename, linefile, output_cont, bloutput, maskmode, chans, thresh, avg_limit, minwidth,
                         edge, blfunc, order, npiece, applyfft, fftthresh, addwn, rejwn, blparam, clipniter,
                         clipthresh, dirkernel, major, minor, pa, kimage, scale, spkernel, kwidth)
    vals.check_args()

    prepare(vals)

    imsmooth(vals)

    convert_image_to_ms(vals)

    sdsmooth(vals)

    sdbaseline(vals)

    imaging(vals)


def prepare(vals: ImBaselineVals = None):
    ia.dohistory(False)

def convert_image_to_ms(vals: ImBaselineVals = None):
    create_empty_ms(vals)
    copy_image_array_to_ms(vals)


### imsmooth ###

def imsmooth(vals: ImBaselineVals = None):

    ia.open(vals.imagename)
    mycsys = ia.coordsys()
    myrg = regionmanager()
    vals.reg = myrg.frombcs(csys=mycsys.torecord(), shape=ia.shape(), chans=vals.imsmooth_chans)
    myrg.done()
    mycsys.done()

    if not vals.dir_ikernel:
        if isinstance(vals.imsmooth_major, (int, float)):
            vals.imsmooth_major = str(vals.imsmooth_major)+'arcsec'
        if isinstance(vals.imsmooth_minor, (int, float)):
            vals.imsmooth_minor = str(vals.imsmooth_minor)+'arcsec'
        if isinstance(vals.imsmooth_pa, (int, float)):
            vals.imsmooth_pa = str(vals.imsmooth_pa)+'deg'
    
    outia = None
    try:
        if vals.dir_gkernel:
            outia = _imsmooth_gckernel(vals)
        elif vals.dir_bkernel:
            outia = _imsmooth_bkernel(vals)
        elif vals.dir_ikernel:
            outia = _imsmooth_ikernel(vals)

    finally:
        ia.done()
        if outia: outia.done()
    
    return _get_image_params_from_imsmooth_output(vals)

def _imsmooth_ikernel(vals):
    return ia.convolve(
                    outfile=vals.imsmooth_output, kernel=vals.imsmooth_kimage, scale=vals.imsmooth_scale, region=vals.reg,
                    mask=vals.imsmooth_mask, overwrite=vals.overwrite, stretch=vals.imsmooth_stretch 
                )

def _imsmooth_bkernel(vals):
    if not vals.imsmooth_major or not vals.imsmooth_minor:
        raise ValueError("Both major and minor must be specified.")
                # BOXCAR KERNEL
                #
                # Until convolve2d supports boxcar we will need to
                # use sepconvolve to do this.
                #
                # BIG NOTE!!!!!
                # According to Gaussian2D documentation the default position
                # angle aligns the major axis along the y-axis, which typically
                # be lat.  So this means that we need to use the major quantity
                # on the y axis (or 1) for sepconvolve.

    casalog.post( "ia.sepconvolve( axes=[0,1],"+\
                            "types=['boxcar','boxcar' ],"+\
                            "widths=[ "+str(vals.imsmooth_minor)+", "+str(vals.imsmooth_major)+" ],"+ \
                            "region="+str(vals.reg)+",outfile="+vals.imsmooth_output+" )",\
                            'DEBUG2' )
                #retValue = ia.sepconvolve( axes=[0,1], types=['box','box' ],\
                #                           widths=[ minor, major ], \
                #                           region=reg,outfile=outfile )
    return ia.sepconvolve(
                    axes=[0,1], types=['box','box'], widths=[ vals.imsmooth_minor, vals.imsmooth_major ],
                    region=vals.reg, outfile=vals.imsmooth_output,
                    mask=vals.imsmooth_mask, overwrite=vals.overwrite, stretch=vals.imsmooth_stretch 
                )

def _imsmooth_gckernel(vals):

    if not vals.imsmooth_major:
        raise ValueError("Major axis must be specified")
    if not vals.imsmooth_minor:
        raise ValueError("Minor axis must be specified")
    if not vals.imsmooth_pa:
        raise ValueError("Position angle must be specified")
        
    return ia.convolve2d(
                    axes=[0,1], region=vals.reg, major=vals.imsmooth_major,
                    minor=vals.imsmooth_minor, pa=vals.imsmooth_pa, outfile=vals.imsmooth_output,
                    mask=vals.imsmooth_mask, stretch=vals.imsmooth_stretch, targetres=vals.imsmooth_targetres,
                    beam=vals.imsmooth_beam, overwrite=vals.overwrite
                )

def _get_image_params_from_imsmooth_output(vals: ImBaselineVals = None):
    try:
        ia.open(vals.imsmooth_output)
        cs = ia.coordsys()
        vals.imshape = ia.shape()
        vals.diraxis = cs.findcoordinate('direction')['world']
        vals.spaxis = cs.findcoordinate('spectral')['world'] # 3 or 2
        vals.polaxis = cs.findcoordinate('stokes')['world']  # 2 or 3 or None
    finally:
        cs.done()
        ia.done()

    assert len(vals.diraxis) > 0

    vals.spaxisexist = len(vals.spaxis) == 1
    vals.polaxisexist = len(vals.polaxis) == 1

    vals.dirshape = vals.imshape[vals.diraxis]
    vals.imnrow = np.prod(vals.dirshape)

    if vals.spaxisexist:
        vals.imnchan = vals.imshape[vals.spaxis[0]]
    else:
        vals.imnchan = 1
    if vals.polaxisexist:
        vals.imnpol = vals.imshape[vals.polaxis[0]]
    else:
        vals.imnpol = 1
    print(f'image shape is {vals.imshape}, direciton {vals.dirshape} ({vals.imnrow} pixels), npol {vals.imnpol}, nchan {vals.imnchan}')

    # if imnchan is too few, say, <10, sdbaseline should abort
    if vals.imnchan < 10:
        print(f'nchan {vals.imnchan} is too few to perform baseline subtraction')
        return False
    return True


### creating empty MS ###

def create_empty_ms(vals: ImBaselineVals = None):

    _check_ms_path(vals.temporary_vis)
    _create_maintable(vals)
    _create_antenna_table(vals)
    _create_data_description_table(vals)
    _create_feed_table(vals)
    _create_field_table(vals)
    _create_flag_cmd_table(vals)
    _create_history_table(vals)
    _create_observation_table(vals)
    _create_pointing_table(vals)
    _create_polarization_table(vals)
    _create_processor_table(vals)
    _create_source_table(vals)
    _create_special_window_table(vals)
    _create_state_table(vals)

def _create_maintable(vals: ImBaselineVals = None, overWrite: bool = True):

    tb = table()
    try:
        tb.create(vals.temporary_vis, EmptyMSBaseInformation.ms_desc, dminfo=EmptyMSBaseInformation.ms_dminfo)
        tb.putkeyword(keyword="MS_VERSION", value=2)
        nrow = tb.nrows()
        nrow_req = vals.imnrow * vals.imnpol
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
        print(f'number of rows {nrow}, number of image pixels {vals.imnrow}, number of pols {vals.imnpol}, required rows {nrow_req}')
    finally:
        tb.close()

def _create_antenna_table(vals: ImBaselineVals = None):
    __create_subtable(vals.temporary_vis, "ANTENNA", EmptyMSBaseInformation.antenna_desc, EmptyMSBaseInformation.antenna_dminfo)
    with sdutil.tbmanager(os.path.join(vals.temporary_vis, 'ANTENNA'), nomodify=False) as tb:
        tb.addrows(1)

def _create_data_description_table(vals: ImBaselineVals = None):
    __create_subtable(vals.temporary_vis, "DATA_DESCRIPTION", EmptyMSBaseInformation.data_description_desc, EmptyMSBaseInformation.data_description_dminfo)
    with sdutil.tbmanager(os.path.join(vals.temporary_vis, 'DATA_DESCRIPTION'), nomodify=False) as tb:
        tb.addrows(1)
        tb.putcell('SPECTRAL_WINDOW_ID', 0, 0)
        tb.putcell('POLARIZATION_ID', 0, 0)

def _create_feed_table(vals: ImBaselineVals = None):
    __create_subtable(vals.temporary_vis, "FEED", EmptyMSBaseInformation.feed_desc, EmptyMSBaseInformation.feed_dminfo)
    with sdutil.tbmanager(os.path.join(vals.temporary_vis, 'FEED'), nomodify=False) as tb:
        tb.addrows(1)

def _create_field_table(vals: ImBaselineVals = None):
    __create_subtable(vals.temporary_vis, "FIELD", EmptyMSBaseInformation.field_desc, EmptyMSBaseInformation.field_dminfo)
    with sdutil.tbmanager(os.path.join(vals.temporary_vis, 'FIELD'), nomodify=False) as tb:
        tb.addrows(1)
        tb.putcell('DELAY_DIR', 0, np.zeros((2,1)))
        tb.putcell('PHASE_DIR', 0, np.zeros((2,1)))
        tb.putcell('REFERENCE_DIR', 0, np.zeros((2,1)))

def _create_flag_cmd_table(vals: ImBaselineVals = None):
    __create_subtable(vals.temporary_vis, "FLAG_CMD", EmptyMSBaseInformation.flag_cmd_desc, EmptyMSBaseInformation.flag_cmd_dminfo)

def _create_history_table(vals: ImBaselineVals = None):
    __create_subtable(vals.temporary_vis, "HISTORY", EmptyMSBaseInformation.history_desc, EmptyMSBaseInformation.history_dminfo)

def _create_observation_table(vals: ImBaselineVals = None):
    __create_subtable(vals.temporary_vis, "OBSERVATION", EmptyMSBaseInformation.observation_desc, EmptyMSBaseInformation.observation_dminfo)
    with sdutil.tbmanager(os.path.join(vals.temporary_vis, 'OBSERVATION'), nomodify=False) as tb:
        tb.addrows(1)

def _create_pointing_table(vals: ImBaselineVals = None):
    __create_subtable(vals.temporary_vis, "POINTING", EmptyMSBaseInformation.pointing_desc, EmptyMSBaseInformation.pointing_dminfo)

def _create_polarization_table(vals: ImBaselineVals = None):
    __create_subtable(vals.temporary_vis, "POLARIZATION", EmptyMSBaseInformation.polarization_desc, EmptyMSBaseInformation.polarization_dminfo)
    with sdutil.tbmanager(os.path.join(vals.temporary_vis, 'POLARIZATION'), nomodify=False) as tb:
        corr_type = np.ones(1, dtype=int)
        corr_product = np.ones(2, dtype=int).reshape((2, 1))
        if tb.nrows() == 0:
            tb.addrows(1)
        tb.putcell('NUM_CORR', 0, 1)
        tb.putcell('CORR_TYPE', 0, corr_type)
        tb.putcell('CORR_PRODUCT', 0, corr_product)

def _create_processor_table(vals: ImBaselineVals = None):
    __create_subtable(vals.temporary_vis, "PROCESSOR", EmptyMSBaseInformation.processor_desc, EmptyMSBaseInformation.processor_dminfo)

def _create_source_table(vals: ImBaselineVals = None):
    __create_subtable(vals.temporary_vis, "SOURCE", EmptyMSBaseInformation.source_desc, EmptyMSBaseInformation.source_dminfo)

def _create_special_window_table(vals: ImBaselineVals = None):
    __create_subtable(vals.temporary_vis, "SPECTRAL_WINDOW", EmptyMSBaseInformation.special_window_desc, EmptyMSBaseInformation.special_window_dminfo)
    with sdutil.tbmanager(os.path.join(vals.temporary_vis, 'SPECTRAL_WINDOW'), nomodify=False) as tb:
        cw = np.ones(vals.imnchan, dtype=float) * 1e6
        cf = 1e9 + np.arange(vals.imnchan, dtype=float) * 1e6

        if tb.nrows() == 0:
            tb.addrows(1)
        tb.putcell('NUM_CHAN', 0, vals.imnchan)
        tb.putcell('CHAN_FREQ', 0, cf)
        tb.putcell('CHAN_WIDTH', 0, cw)
        tb.putcell('RESOLUTION', 0, cw)
        tb.putcell('EFFECTIVE_BW', 0, cw)
        tb.putcell('TOTAL_BANDWIDTH', 0, cw.sum())

def _create_state_table(vals: ImBaselineVals = None):
    __create_subtable(vals.temporary_vis, "STATE", EmptyMSBaseInformation.state_desc, EmptyMSBaseInformation.state_dminfo)
    with sdutil.tbmanager(os.path.join(vals.temporary_vis, 'STATE'), nomodify=False) as tb:
        if tb.nrows() == 0:
            tb.addrows(1)
        tb.putcell('OBS_MODE', 0, 'OBSERVE_TARGET#ON_SOURCE_IMAGE_DOMAIN')

def _check_ms_path(ms: str = None, overWrite: bool = True):
    exists = os.path.exists(ms)
    if overWrite and exists:
        shutil.rmtree(ms)
    return exists

def __create_subtable(ms: str, subtable: str, desc: str, dminfo: str):
    tb = table()
    try:
        tb.create(f"{ms}/{subtable}", desc, dminfo=dminfo)
        tb.open(ms, nomodify=False)
        tb.putkeyword(subtable, f"Table: {ms}/{subtable}")
    finally:
        tb.close()


def copy_image_array_to_ms(vals: ImBaselineVals = None):
    # get image array and mask from the image
    ia.open(vals.imsmooth_output)
    arr = ia.getchunk()
    msk = ia.getchunk(getmask=True)
    ia.done()

    # put image array slices to MS DATA column
    # axis indices for spatial, spectral and polarization axes
    xax, yax = vals.diraxis
    spax = vals.spaxis[0]
    if vals.polaxisexist:
        polax = vals.polaxis[0]
    else:
        arr = np.expand_dims(arr, axis=3)
        msk = np.expand_dims(msk, axis=3)
        polax = 3
    casalog.post(f'axis index: {xax} {yax} {polax} {spax}')

    # which data column to use
    with sdutil.tbmanager(vals.temporary_vis, nomodify=False) as tb:
        # also set FLAG, SIGMA, WEIGHT, and UVW
        # r2p is a mapping between MS row and image array slice
        xax, yax = vals.diraxis

        index_list = [0, 0, 0, 0]
        index_list[spax] = np.arange(vals.imnchan)
        r2p = []
        nx, ny = vals.dirshape
        irow = 0
        wgt = np.ones(1, dtype=float)
        uvw = np.zeros(3, dtype=float)
        for ix in range(nx):
            index_list[xax] = ix
            for iy in range(ny):
                index_list[yax] = iy
                for ip in range(vals.imnpol):
                    index_list[polax] = ip
                    slice = tuple(index_list)
                    subarr = arr[slice]
                    submsk = np.logical_not(msk[slice])
                    print(f'slice={slice}, subarr={subarr[200:-200]}')
                    r2p.append(slice)
                    tb.putcell(vals.datacolumn, irow, np.expand_dims(subarr, axis=0))
                    tb.putcell('FLAG', irow, np.expand_dims(submsk, axis=0))
                    tb.putcell('SIGMA', irow, wgt)
                    tb.putcell('WEIGHT', irow, wgt)
                    tb.putcell('UVW', irow, uvw)
                    irow += 1


### sdsmooth ###

def sdsmooth(vals: ImBaselineVals = None):
    sdms.open(vals.temporary_vis)
    sdms.set_selection(spw=vals.sdsmooth_spw, field=vals.sdsmooth_field, 
                           antenna=vals.sdsmooth_antenna,
                           timerange=vals.sdsmooth_timerange, scan=vals.sdsmooth_scan,
                           polarization=vals.sdsmooth_pol, intent=vals.sdsmooth_intent,
                           reindex=vals.sdsmooth_reindex)
    sdms.smooth(type=vals.sdsmooth_kernel, width=vals.sdsmooth_kwidth, outfile=vals.sdsmooth_output, datacolumn=vals.datacolumn.lower())
    sdms.close()


### sdbaseline ###

def sdbaseline(vals: ImBaselineVals = None):
    print("start sdbaseline")
    vals.sdbaseline_blfunc = vals.sdbaseline_blfunc.lower()
    
    if isinstance(vals.sdbaseline_fftthresh, str):
        vals.sdbaseline_fftthresh = vals.sdbaseline_fftthresh.lower()

    try:

        blparam_file = vals.sdsmooth_output + '_blparam.txt'
        if os.path.exists(blparam_file):
            _remove_data(blparam_file)  # CAS-11781

        spw = vals.sdbaseline_spw
        if vals.sdbaseline_spw == '': spw = '*'
        # blmode = 'fit'

        if vals.sdbaseline_blfunc == 'sinusoid':
            vals.sdbaseline_addwn = sdutil.parse_wavenumber_param(vals.sdbaseline_addwn)
            vals.sdbaseline_rejwn = sdutil.parse_wavenumber_param(vals.sdbaseline_rejwn)
            _check_fftthresh(vals.sdbaseline_fftthresh)

        blformat, vals.sdbaseline_bloutput = _prepare_for_blformat_bloutput(vals.sdsmooth_output, 'text', vals.sdbaseline_bloutput, True)

        _output_bloutput_text_header(blformat, vals.sdbaseline_bloutput,
                                    vals.sdbaseline_blfunc, vals.sdbaseline_maskmode,
                                    vals.sdsmooth_output, vals.sdbaseline_output)
        
        if vals.sdbaseline_blfunc == 'variable':
            sorttab_info = _remove_sorted_table_keyword(vals.sdsmooth_output)

        selection = ms.msseltoindex(vis=vals.sdsmooth_output, spw=spw, field=vals.sdbaseline_field, 
                                    baseline='', time='', scan=vals.sdbaseline_scan)
        sdms.open(vals.sdsmooth_output)
        sdms.set_selection(spw=sdutil.get_spwids(selection),
                            field=vals.sdbaseline_field, antenna=vals.sdbaseline_antenna,
                            timerange=vals.sdbaseline_timerenge, scan=vals.sdbaseline_scan,
                            polarization=vals.sdbaseline_pol, intent=vals.sdbaseline_intent,
                            reindex=vals.sdbaseline_reindex)
        params, func = _prepare_for_baselining(sdms, blfunc=vals.sdbaseline_blfunc,
                                                datacolumn=vals.datacolumn.lower(),
                                                outfile=vals.sdbaseline_output,
                                                bloutput=','.join(vals.sdbaseline_bloutput),
                                                dosubtract=vals.sdbaseline_dosubtract,
                                                spw=spw,
                                                pol=vals.sdbaseline_pol,
                                                linefinding=(vals.sdbaseline_maskmode=='auto'),
                                                threshold=vals.sdbaseline_thresh,
                                                avg_limit=vals.sdbaseline_avg_limit,
                                                minwidth=vals.sdbaseline_minwidth,
                                                edge=vals.sdbaseline_edge,
                                                order=vals.sdbaseline_order,
                                                npiece=vals.sdbaseline_npiece,
                                                applyfft=vals.sdbaseline_applyfft,
                                                fftmethod=vals.sdbaseline_fftmethod,
                                                fftthresh=vals.sdbaseline_fftthresh,
                                                addwn=vals.sdbaseline_addwn,
                                                rejwn=vals.sdbaseline_rejwn,
                                                clip_threshold_sigma=vals.sdbaseline_clipthresh,
                                                num_fitting_max=vals.sdbaseline_clipniter+1,
                                                blparam=vals.sdbaseline_blparam,
                                                verbose=False,
                                                updateweight=vals.sdbaseline_updateweight,
                                                sigmavalue=vals.sdbaseline_sigmavalue
                                                )
        func(**params)
        sdms.close()
        
        if vals.sdbaseline_blfunc == 'variable':
            _restore_sorted_table_keyword(vals.sdsmooth_output, sorttab_info)


        # Write history to outfile
        # param_names = imbaseline.__code__.co_varnames[:imbaseline.__code__.co_argcount]
        # vars = locals()
        # param_vals = [vars[p] for p in param_names]
        # write_history(ms, outfile, 'imbaseline', param_names, param_vals, casalog)


    except Exception:
        raise

blformat_item = ['csv', 'text', 'table']
blformat_ext  = ['csv', 'txt',  'bltable']

def _remove_data(filename):
    if os.path.exists(filename):
        if os.path.isdir(filename):
            shutil.rmtree(filename)
        elif os.path.isfile(filename):
            os.remove(filename)
        else:
            # could be a symlink
            os.remove(filename)

def _check_fftthresh(fftthresh):
    has_valid_type = isinstance(fftthresh, float) or isinstance(fftthresh, int) or isinstance(fftthresh, str)
    if not has_valid_type:
        raise ValueError('fftthresh must be float or integer or string.')

    not_positive_mesg = 'threshold given to fftthresh must be positive.'
    
    if isinstance(fftthresh, str):
        try:
            val_not_positive = False
            if (3 < len(fftthresh)) and (fftthresh[:3] == 'top'):
                val_top = int(fftthresh[3:])
                if (val_top <= 0):
                    val_not_positive = True
            elif (5 < len(fftthresh)) and (fftthresh[-5:] == 'sigma'):
                val_sigma = float(fftthresh[:-5])
                if (val_sigma <= 0.0):
                    val_not_positive = True
            else:
                val_sigma = float(fftthresh)
                if (val_sigma <= 0.0):
                    val_not_positive = True
            
            if val_not_positive:
                raise ValueError(not_positive_mesg)
        except Exception as e:
            if (str(e) == not_positive_mesg):
                raise
            else:
                raise ValueError('fftthresh has a wrong format.')

    else:
        if (fftthresh <= 0.0):
            raise ValueError(not_positive_mesg)

def _remove_sorted_table_keyword(infile):
    res = {'is_sorttab': False, 'sorttab_keywd': '', 'sorttab_name': ''}
    with sdutil.tbmanager(infile, nomodify=False) as tb:
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

def _restore_sorted_table_keyword(infile, sorttab_info):
    if sorttab_info['is_sorttab'] and (sorttab_info['sorttab_name'] != ''):
        with sdutil.tbmanager(infile, nomodify=False) as tb:
            try:
                tb.putkeyword(sorttab_info['sorttab_keywd'],
                              sorttab_info['sorttab_name'])
            except Exception:
                raise

def _prepare_for_baselining(sdms, **keywords):
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
    if blfunc!= 'variable':
        keys += ['clip_threshold_sigma', 'num_fitting_max']
        keys += ['linefinding', 'threshold', 'avg_limit', 'minwidth', 'edge']
    for key in keys: params[key] = keywords[key]

    baseline_func = getattr(sdms, funcname)

    return params, baseline_func


def _prepare_for_blformat_bloutput(infile, blformat, bloutput, overwrite):
    # force to string list
    blformat = __force_to_string_list(blformat, 'blformat')
    bloutput = __force_to_string_list(bloutput, 'bloutput')

    # the default bloutput value '' is expanded to a list 
    # with length of blformat, and with '' throughout.
    if (bloutput == ['']): bloutput *= len(blformat)

    # check length
    if (len(blformat) != len(bloutput)):
        raise ValueError('blformat and bloutput must have the same length.')

    # check duplication
    if __has_duplicate_nonnull_element(blformat):
        raise ValueError('duplicate elements in blformat.')
    if __has_duplicate_nonnull_element_ex(bloutput, blformat):
        raise ValueError('duplicate elements in bloutput.')

    # fill bloutput items to be output, then rearrange them
    # in the order of blformat_item.
    bloutput = __normalise_bloutput(infile, blformat, bloutput, overwrite)

    return blformat, bloutput

def __force_to_string_list(s, name):
    mesg = '%s must be string or list of string.' % name
    if isinstance(s, str): s = [s]
    elif isinstance(s, list):
        for i in range(len(s)):
            if not isinstance(s[i], str):
                raise ValueError(mesg)
    else:
        raise ValueError(mesg)
    return s

def __has_duplicate_nonnull_element(in_list):
    #return True if in_list has duplicated elements other than ''
    duplicates = [key for key, val in Counter(in_list).items() if val > 1]
    len_duplicates = len(duplicates)
    
    if (len_duplicates >= 2):
        return True
    elif (len_duplicates == 1):
        return (duplicates[0] != '')
    else: #len_duplicates == 0
        return False


def __has_duplicate_nonnull_element_ex(lst, base):
    # lst and base must have the same length.
    #
    # (1) extract elements from lst and make a new list
    #     if the element of base with the same index
    #     is not ''.
    # (2) check if the list made in (1) has duplicated
    #     elements other than ''.
    
    return __has_duplicate_nonnull_element(
        [lst[i] for i in range(len(lst)) if base[i] != ''])

def __normalise_bloutput(infile, blformat, bloutput, overwrite):
    normalised_bloutput = []
    for item in zip(blformat_item, blformat_ext):
        normalised_bloutput.append(
            __get_normalised_name(infile, blformat, bloutput, item[0], item[1], overwrite))
    return normalised_bloutput

def __get_normalised_name(infile, blformat, bloutput, name, ext, overwrite):
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

def _output_bloutput_text_header(blformat, bloutput, blfunc, maskmode, infile, outfile):
    fname = bloutput[blformat_item.index('text')]
    if (fname == ''): return
    
    f = open(fname, 'w')

    blf = blfunc.lower()
    if (blf == 'poly'):
        ftitles = ['Fit order']
    elif (blf == 'chebyshev'):
        ftitles = ['Fit order']
    elif (blf == 'cspline'):
        ftitles = ['nPiece']
    elif (blf=='sinusoid'):
        ftitles = ['applyFFT', 'fftMethod', 'fftThresh', 'addWaveN', 'rejWaveN']
    elif (blf=='variable'):
        ftitles = []
    else:
        raise ValueError("Unsupported blfunc = %s" % blfunc)
        

    mm = maskmode.lower()
    if (mm == 'auto'):
        mtitles = ['Threshold', 'avg_limit', 'Edge']
    elif (mm == 'list'):
        mtitles = []
    else: # interact
        mtitles = []

    ctitles = ['clipThresh', 'clipNIter']

    info = [['Source Table', infile],
            ['Output File', outfile if (outfile != '') else infile],
            ['Mask mode', maskmode]]

    separator = '#' * 60 + '\n'
    
    f.write(separator)
    for i in range(len(info)):
        f.write('%12s: %s\n' % tuple(info[i]))
    f.write(separator)
    f.write('\n')
    f.close()


### imaging ###

def imaging(vals: ImBaselineVals = None):
    print("start imaging")
    vals.output_imagename = os.path.basename(vals.imagename) + "_out"
    _make_output_file(vals)

    nx, ny = vals.dirshape
    arr = np.empty((nx,ny,vals.imnchan))
    if vals.polaxisexist:
        arr = np.expand_dims(arr, vals.polaxis)
    pos = 0
    with sdutil.tbmanager(vals.sdbaseline_output) as tb:
        ndim = tb.getcell(vals.datacolumn, 0).ndim
        for i in range(nx):
            for j in range(ny):
                if ndim == 3:
                    arr[i][j][0] =  tb.getcell(vals.datacolumn, pos)[0].real
                else:
                    arr[i][j] =  tb.getcell(vals.datacolumn, pos)[0].real
                pos += 1

    print("end arraying")

    ia.open(vals.output_imagename)
    ia.putchunk(pixels=arr, locking=True)
    ia.done()
    print("end imaging")


def _make_output_file(vals: ImBaselineVals = None):
    shutil.copytree(vals.imagename, vals.output_imagename)


class EmptyMSBaseInformation:

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


# temporary main

workdir = 'working'
if os.path.exists(workdir):
    shutil.rmtree(workdir)
os.mkdir(workdir)
os.chdir(workdir)

try:
    sample_image = '/nfsstore/casa_share/casatestdata/image/ref_multipix.signalband'
    #sample_image = '/nfsstore/casa_share/casatestdata/image/expected.im'

    # pipeline_regression_dir = '/nfsstore/nakazato/pipeline/regression_ref'
    # sample_project = 'm100'
    # sample_fits_image = 'uid___A002_X345678_X90.M100_sci.virtspw23.cube.I.sd.fits'
    # sample_image = 'm100_spw23.image'

    # importfits(fitsimage=f'{pipeline_regression_dir}/{sample_project}/products/{sample_fits_image}', imagename=sample_image)

    imbaseline(imagename=sample_image, linefile="output_ms.image", dirkernel="boxcar", spkernel="boxcar", major = "5arcsec", minor = "2arcsec", pa = "20deg", 
    blfunc='poly', applyfft=True)

finally:
    os.chdir("..")
