# image-based line finding and baseline subtraction.

from abc import abstractmethod
import os
import sys
import shutil
from typing import Any, Dict, List, Tuple
import numpy as np
from numpy import array, uint64
import uuid

from casatools import image, table
from casatasks import casalog
from casatasks.private.sdutil import table_manager, tool_manager, sdtask_decorator
from casatasks.private.ialib import write_image_history
from casatasks.private.task_imsmooth import imsmooth
from casatasks.private.task_sdsmooth import sdsmooth
from casatasks.private.task_sdbaseline import sdbaseline

ia = image()

DATACOLUMN = 'DATA'
OVERWRITE = True
debug = False


@sdtask_decorator
def imbaseline(imagename=None, linefile=None, output_cont=None, bloutput=None, maskmode=None, chans=None, thresh=None,
               avg_limit=None, minwidth=None, edge=None, blfunc=None, order=None, npiece=None, applyfft=None, fftthresh=None,
               addwn=None, rejwn=None, blparam=None, clipniter=None, clipthresh=None, dirkernel=None, major=None, minor=None,
               pa=None, kimage=None, scale=None, spkernel=None, kwidth=None):
    """
    THE MAIN METHOD OF IMBASELINE.

    All specifications of arguments are defined in:
    https://open-jira.nrao.edu/browse/CAS-13520
    """
    __validate_imagename(imagename)
    linefile = __prepare_linefile(linefile, imagename)
    if not output_cont:
        output_cont = False
    output_cont_file = os.path.basename(imagename) + ".cont"
    erase_queue = EraseQueue()

    prepare()

    imsmooth_output = execute_imsmooth(imagename, dirkernel, major, minor, pa, kimage, scale, erase_queue)

    image2ms_output, image2ms_params = execute_image2ms(imsmooth_output, DATACOLUMN, erase_queue)

    sdsmooth_output = execute_sdsmooth(image2ms_output, DATACOLUMN.lower(), spkernel, kwidth, erase_queue)

    sdbaseline_output = execute_sdbaseline(sdsmooth_output, DATACOLUMN.lower(), bloutput, maskmode, chans, thresh, avg_limit,
                                           minwidth, edge, blfunc, order, npiece, applyfft, fftthresh, addwn, rejwn, blparam,
                                           clipniter, clipthresh, erase_queue)

    execute_ms2image(sdbaseline_output, linefile, imagename, DATACOLUMN, output_cont, output_cont_file, image2ms_params)

    erase_queue.clear(False)


def execute_imsmooth(imagename, dirkernel, major, minor, pa, kimage, scale, erase_queue):
    output_filepath = imagename

    if __validate_imsmooth_execution(dirkernel):
        casalog.post("execute image smoothing", "INFO")
        imsmooth_output = __generate_temporary_filename("dirsmooth-", "im")
        args, kargs = ImsmoothParams(imagename, imsmooth_output, dirkernel, major, minor, pa, kimage, scale)()
        imsmooth(*args, **kargs)
        output_filepath = imsmooth_output
        erase_queue.append(Erasable(imsmooth_output))
    else:
        casalog.post("omit image smoothing", "INFO")
        erase_queue.append(Unerasable(output_filepath))

    return output_filepath


def execute_image2ms(infile, datacolumn, erase_queue) -> Tuple[Any]:
    casalog.post("convert casaimage to MeasurementSet", "INFO")
    img2ms_output = __generate_temporary_filename("img2ms-", "ms")
    i2mp = image2ms(Image2MSParams(infile, img2ms_output, datacolumn))
    erase_queue.append(Erasable(img2ms_output))
    return img2ms_output, i2mp


def execute_sdsmooth(infile=None, datacolumn=None, spkernel=None, kwidth=None, erase_queue=None):
    sdsmooth_output = infile

    if __validate_sdsmooth_execution(spkernel):
        casalog.post("execute spectral smoothing", "INFO")
        sdsmooth_output = __generate_temporary_filename("spsmooth-", "ms")
        sdsmooth(**SdsmoothParams(infile, sdsmooth_output, datacolumn, spkernel, kwidth)())
        erase_queue.append(Erasable(sdsmooth_output))
    else:
        casalog.post("omit spectral smoothing", "INFO")

    return sdsmooth_output


def execute_sdbaseline(infile, datacolumn, bloutput, maskmode, chans, thresh, avg_limit, minwidth, edge, blfunc, order, npiece,
                       applyfft, fftthresh, addwn, rejwn, blparam, clipniter, clipthresh, erase_queue):
    casalog.post("execute spectral baselining", "INFO")
    sdbaseline_output = __generate_temporary_filename("baseline-", "ms")
    sdbaseline(**SdbaselineParams(infile, sdbaseline_output, datacolumn, bloutput, maskmode, chans, thresh, avg_limit, minwidth,
               edge, blfunc, order, npiece, applyfft, fftthresh, addwn, rejwn, blparam, clipniter, clipthresh)())
    erase_queue.append(Erasable(sdbaseline_output))

    return sdbaseline_output


def execute_ms2image(infile, linefile, imagefile, datacolumn, output_cont, output_cont_file, i2ms):
    casalog.post("convert MeasurementSet to casaimage", "INFO")
    ms2image(MS2ImageParams(infile, linefile, imagefile, datacolumn, output_cont, output_cont_file, i2ms))


def __validate_imagename(imagename):
    if not os.path.exists(imagename):
        raise ValueError(f'Error: file {imagename} is not found.', 'SEVERE')


def __prepare_linefile(linefile, imagename):
    if linefile == '':
        linefile = os.path.basename(imagename).rstrip('/') + '_bs'
    if not OVERWRITE and os.path.exists(linefile):
        raise ValueError(f'Error: file {linefile} already exists, please delete before continuing.', 'SEVERE')
    return linefile


def __generate_temporary_filename(prefix: str='', ext: str='') -> str:
    if ext != '':
        ext = '.' + ext
    while True:
        filename = prefix + str(uuid.uuid4()) + ext
        if not os.path.exists(filename):
            return filename


def __validate_imsmooth_execution(dirkernel):
    def is_valid_kernel(kernel):
        return kernel == 'none', kernel == 'image' or kernel == 'boxcar' or kernel == 'gaussian'

    none, valid = is_valid_kernel(dirkernel)
    if not none and not valid:
        raise ValueError(f'Unsupported direction smoothing kernel, {dirkernel}', 'SEVERE')
    if valid:
        return True
    return False


def __validate_sdsmooth_execution(spkernel):
    def is_valid_kernel(kernel):
        return kernel == 'none', kernel == 'boxcar' or kernel == 'gaussian'

    none, valid = is_valid_kernel(spkernel)
    if not none and not valid:
        raise ValueError(f'Unsupported spectral smoothing kernel, {spkernel}', 'SEVERE')
    if valid:
        return True
    return False


def prepare():
    ia.dohistory(False)


class Validable:

    @abstractmethod
    def validate(self):
        raise RuntimeError('Not implemented')


class AbstractEraseable:

    def __init__(self, filename):
        self.path = filename

    @abstractmethod
    def erase(self):
        raise RuntimeError('Not implemented')


class Erasable(AbstractEraseable):

    def erase(self, dry_run=True):
        if not dry_run:
            casalog.post(f"erase file:{self.path}", "DEBUG2")
            if os.path.exists(self.path):
                shutil.rmtree(self.path)


class Unerasable(AbstractEraseable):

    def erase(self, dry_run=True):
        casalog.post(f"un-erase file:{self.path}", "DEBUG2")


class EraseQueue:
    def __init__(self) -> None:
        self.queue = []

    def append(self, file: AbstractEraseable=None):
        if isinstance(file, AbstractEraseable):
            self.queue.append(file)
        else:
            raise ValueError('cannot append to erase queue')

    def clear(self, dry_run=True):
        for obj in self.queue:
            obj.erase(dry_run)
        self.queue.clear()


class ImsmoothParams(Validable):

    TARGETRES = False
    MASK = ''
    BEAM = {}
    REGION = ''
    BOX = ''
    CHANS = ''
    STOKES = ''
    STRETCH = False
    OVERWRITE = True

    def __init__(self, infile: str=None, outfile: str=None, dirkernel: str='none', major: str='', minor: str='', pa: str='',
                 kimage: str='', scale: int=-1.0) -> None:
        self.infile = infile
        self.outfile = outfile
        self.kernel = dirkernel if dirkernel is not None else 'none'       # none(default)/gaussian/boxcar/image
        self.major = major if major is not None else ''                    # dirkernel = gaussian/boxcar
        self.minor = minor if minor is not None else ''                    # dirkernel = gaussian/boxcar
        self.pa = pa if pa is not None else ''                             # dirkernel = gaussian/boxcar
        self.kimage = kimage if kimage is not None else ''                 # dirkernel = image
        self.scale = scale if scale is not None else -1.0                  # dirkernel = image
        self.validate()

    def validate(self):
        self.__validate_dirkernel()

    def __validate_dirkernel(self):
        if self.kernel == 'image':
            self.major = self.minor = self.pa = ''
            if self.kimage != '' and not os.path.exists(self.kimage):
                raise ValueError(f'Error: file {self.kimage} is not found.', 'SEVERE')
        else:  # bkernel/gkernel
            self.kimage = ''
            self.scale = -1.0

    def __call__(self) -> Tuple:
        return [self.infile, self.kernel, self.major, self.minor, self.pa, self.TARGETRES, self.kimage, self.scale,
                self.REGION, self.BOX, self.CHANS, self.STOKES, self.MASK, self.outfile, self.STRETCH, self.OVERWRITE,
                self.BEAM], dict(__taskcaller__="imbaseline")


class SdsmoothParams(Validable):

    SPW = ''
    FIELD = ''
    ANTENNA = ''
    TIMERANGE = ''
    SCAN = ''
    POL = ''
    INTENT = ''
    REINDEX = True
    OVERWRITE = True

    def __init__(self, infile: str=None, outfile: str=None, datacolumn: str=None, spkernel: str='none', kwidth: int=5) -> None:
        self.infile = infile
        self.outfile = outfile
        self.datacolumn = datacolumn
        self.kernel = spkernel if spkernel is not None else 'none'   # none(default)/gaussian/boxcar
        self.kwidth = kwidth if kwidth is not None else 5            # gaussian/boxcar
        self.validate()

    def validate(self):
        pass

    def __call__(self):
        return dict(infile=self.infile, datacolumn=self.datacolumn, antenna=self.ANTENNA, field=self.FIELD, spw=self.SPW,
                    timerange=self.TIMERANGE, scan=self.SCAN, pol=self.POL, intent=self.INTENT, reindex=self.REINDEX,
                    kernel=self.kernel, kwidth=self.kwidth, outfile=self.outfile, overwrite=self.OVERWRITE,
                    __taskcaller__="imbaseline")


class SdbaselineParams(Validable):

    ANTENNA = ''
    FIELD = ''
    SPW = ''
    TIMERANGE = ''
    SCAN = ''
    POL = ''
    INTENT = ''
    REINDEX = True
    BLMODE = 'fit'
    DOSUBTRACT = True
    BLFORMAT = 'text'
    BLTABLE = ''
    UPDATEWEIGHT = False
    SIGMAVALUE = 'stddev'
    SHOWPROGRESS = False
    MINNROW = 1000
    FFTMETHOD = 'fft'
    VERBOSE = False
    OVERWRITE = True

    def __init__(self, infile: str=None, outfile: str=None, datacolumn: str=None, bloutput: str='', maskmode: str='list',
                 chans: str='', thresh: float=5.0, avg_limit: int=4, minwidth: int=4, edge: List[int]=[0, 0], blfunc: str='poly',
                 order: int=5, npiece: int=3, applyfft: bool=True, fftthresh: float=3.0, addwn: List=[0], rejwn: List=[],
                 blparam: str='', clipniter: int=0, clipthresh: float=3.0) -> None:
        self.infile = infile
        self.outfile = outfile
        self.datacolumn = datacolumn
        self.bloutput = bloutput if bloutput is not None else ''
        self.maskmode = maskmode.lower() if maskmode is not None else 'list'  # list(default)/auto
        self.chans = chans if chans is not None else ''                       # maskmode = list
        self.thresh = thresh if thresh is not None else 5.0                  # maskmode = auto
        self.avg_limit = avg_limit if avg_limit is not None else 4            # maskmode = auto
        self.minwidth = minwidth if minwidth is not None else 4               # maskmode = auto
        self.edge = edge if edge is not None else [0, 0]                      # maskmode = auto
        self.blfunc = blfunc.lower() if blfunc is not None else 'poly'        # poly(default)/chebyshev/cspline/sinusoid/variable
        self.order = order if order is not None else 5                        # blfunc = poly/chebyshev
        self.npiece = npiece if npiece is not None else 3                     # blfunc = cspline
        self.applyfft = applyfft if applyfft is not None else True            # blfunc = sinusoid
        self.fftthresh = fftthresh if fftthresh is not None else 3.0          # blfunc = sinusoid
        self.addwn = addwn if addwn is not None else [0]                      # blfunc = sinusoid
        self.rejwn = rejwn if rejwn is not None else []                       # blfunc = sinusoid
        self.blparam = blparam if blparam is not None else ''                 # blfunc = variable
        self.clipniter = clipniter if clipniter is not None else 0
        self.clipthresh = clipthresh if clipthresh is not None else 3.0

    def validate(self):
        self.__validate_maskmode()
        self.__validate_blfunc()

    def __validate_maskmode(self):
        maskmode_list = self.maskmode == 'list'
        maskmode_auto = self.maskmode == 'auto'
        if not (maskmode_list or maskmode_auto):
            raise ValueError(f'Unsupported maskmode, {self.maskmode}', 'SEVERE')

    def __validate_blfunc(self):
        def is_valid_blfunc(self):
            return self.blfunc == 'poly' or self.blfunc == 'chebyshev' or self.blfunc == 'cspline' or self.blfunc == 'sinusoid' \
                or self.blfunc == 'variable'

        if not is_valid_blfunc():
            raise ValueError(f'Unsupported blfunc, {self.blfunc}', 'SEVERE')

        if self.blfunc == 'variable' and not os.path.exists(self.blparam):
            raise ValueError(f"input file '{self.blparam}' does not exists", 'SEVERE')

    def __call__(self) -> Dict[str, Any]:
        return dict(infile=self.infile, datacolumn=self.datacolumn, antenna=self.ANTENNA, field=self.FIELD, spw=self.SPW,
                    timerange=self.TIMERANGE, scan=self.SCAN, pol=self.POL, intent=self.INTENT, reindex=self.REINDEX,
                    maskmode=self.maskmode, thresh=self.thresh, avg_limit=self.avg_limit, minwidth=self.minwidth, edge=self.edge,
                    blmode=self.BLMODE, dosubtract=self.DOSUBTRACT, blformat=self.BLFORMAT, bloutput=self.bloutput, bltable=self.BLTABLE,
                    blfunc=self.blfunc, order=self.order, npiece=self.npiece, applyfft=self.applyfft, fftmethod=self.FFTMETHOD,
                    fftthresh=self.fftthresh, addwn=self.addwn, rejwn=self.rejwn, clipthresh=self.clipthresh, clipniter=self.clipniter,
                    blparam=self.blparam, verbose=self.VERBOSE, updateweight=self.UPDATEWEIGHT, sigmavalue=self.SIGMAVALUE,
                    showprogress=self.SHOWPROGRESS, minnrow=self.MINNROW, outfile=self.outfile, overwrite=self.OVERWRITE,
                    __taskcaller__="imbaseline")


class Image2MSParams(Validable):

    def __init__(self, infile: str=None, outfile: str=None, datacolumn: str='DATA') -> None:
        self.infile = infile
        self.outfile = outfile
        self.datacolumn = datacolumn
        self.validate()

    def validate(self):
        self.__validate_outfile()

    def __validate_outfile(self):
        if os.path.exists(self.outfile):
            raise ValueError(f"Folder exists:{self.outfile}")


def image2ms(params: Image2MSParams=None):
    __get_image_params_from_imsmooth_output(params)
    __create_empty_ms(params)
    __copy_image_array_to_ms(params)
    return params


def __get_axis_position(val: array=None):
    if val is not None and len(val) > 0:
        return val[0].item()
    return -1


def __get_image_params_from_imsmooth_output(params):
    with tool_manager(params.infile, image) as _ia:
        try:
            cs = _ia.coordsys()
            params.im_shape = _ia.shape()
            params.axis_dir = cs.findcoordinate('direction')['world']
            params.axis_sp = __get_axis_position(cs.findcoordinate('spectral')['world'])  # 3 or 2 or -1
            params.axis_pol = __get_axis_position(cs.findcoordinate('stokes')['world'])   # 2 or 3 or -1
        finally:
            cs.done()

    assert len(params.axis_dir) > 0

    params.dir_shape = params.im_shape[params.axis_dir]
    params.im_nrow = np.prod(params.dir_shape)
    params.im_nchan = params.im_shape[params.axis_sp] if params.axis_sp > 0 else 1
    params.im_npol = params.im_shape[params.axis_pol] if params.axis_pol > 0 else 1

    casalog.post(f'image shape is {params.im_shape}, direciton {params.dir_shape} ({params.im_nrow} pixels), '
                 f'npol {params.im_npol}, nchan {params.im_nchan}', 'DEBUG2')

    # if im_nchan is too few, say, <10, sdbaseline should abort
    if params.im_nchan < 10:
        casalog.post(f'nchan {params.im_nchan} is too few to perform baseline subtraction', 'DEBUG2')
        return False
    return True


def __create_empty_ms(params):
    __check_ms_path(params)
    __create_maintable(params)
    __create_antenna_table(params)
    __create_data_description_table(params)
    __create_feed_table(params)
    __create_field_table(params)
    __create_flag_cmd_table(params)
    __create_history_table(params)
    __create_observation_table(params)
    __create_pointing_table(params)
    __create_polarization_table(params)
    __create_processor_table(params)
    __create_source_table(params)
    __create_special_window_table(params)
    __create_state_table(params)


def __create_maintable(params):
    tb = table()
    try:
        tb.create(params.outfile, EmptyMSBaseInformation.ms_desc, dminfo=EmptyMSBaseInformation.ms_dminfo)
        tb.putkeyword(keyword="MS_VERSION", value=2)
        nrow = tb.nrows()
        nrow_req = params.im_nrow * params.im_npol
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
        casalog.post(f'number of rows {nrow}, number of image pixels {params.im_nrow}, number of pols {params.im_npol}, '
                     f'required rows {nrow_req}', 'DEBUG2')
    finally:
        tb.close()


def __create_antenna_table(params):
    __create_subtable(params.outfile, "ANTENNA", EmptyMSBaseInformation.antenna_desc, EmptyMSBaseInformation.antenna_dminfo)
    with table_manager(os.path.join(params.outfile, 'ANTENNA'), nomodify=False) as tb:
        tb.addrows(1)


def __create_data_description_table(params):
    __create_subtable(params.outfile, "DATA_DESCRIPTION", EmptyMSBaseInformation.data_description_desc,
                      EmptyMSBaseInformation.data_description_dminfo)
    with table_manager(os.path.join(params.outfile, 'DATA_DESCRIPTION'), nomodify=False) as tb:
        tb.addrows(1)
        tb.putcell('SPECTRAL_WINDOW_ID', 0, 0)
        tb.putcell('POLARIZATION_ID', 0, 0)


def __create_feed_table(params):
    __create_subtable(params.outfile, "FEED", EmptyMSBaseInformation.feed_desc, EmptyMSBaseInformation.feed_dminfo)
    with table_manager(os.path.join(params.outfile, 'FEED'), nomodify=False) as tb:
        tb.addrows(1)


def __create_field_table(params):
    __create_subtable(params.outfile, "FIELD", EmptyMSBaseInformation.field_desc, EmptyMSBaseInformation.field_dminfo)
    with table_manager(os.path.join(params.outfile, 'FIELD'), nomodify=False) as tb:
        tb.addrows(1)
        tb.putcell('DELAY_DIR', 0, np.zeros((2, 1)))
        tb.putcell('PHASE_DIR', 0, np.zeros((2, 1)))
        tb.putcell('REFERENCE_DIR', 0, np.zeros((2, 1)))


def __create_flag_cmd_table(params):
    __create_subtable(params.outfile, "FLAG_CMD", EmptyMSBaseInformation.flag_cmd_desc, EmptyMSBaseInformation.flag_cmd_dminfo)


def __create_history_table(params):
    __create_subtable(params.outfile, "HISTORY", EmptyMSBaseInformation.history_desc, EmptyMSBaseInformation.history_dminfo)


def __create_observation_table(params):
    __create_subtable(params.outfile, "OBSERVATION", EmptyMSBaseInformation.observation_desc, EmptyMSBaseInformation.observation_dminfo)
    with table_manager(os.path.join(params.outfile, 'OBSERVATION'), nomodify=False) as tb:
        tb.addrows(1)


def __create_pointing_table(params):
    __create_subtable(params.outfile, "POINTING", EmptyMSBaseInformation.pointing_desc, EmptyMSBaseInformation.pointing_dminfo)


def __create_polarization_table(params):
    __create_subtable(params.outfile, "POLARIZATION", EmptyMSBaseInformation.polarization_desc, EmptyMSBaseInformation.polarization_dminfo)
    with table_manager(os.path.join(params.outfile, 'POLARIZATION'), nomodify=False) as tb:
        corr_type = np.ones(1, dtype=int)
        corr_product = np.ones(2, dtype=int).reshape((2, 1))
        if tb.nrows() == 0:
            tb.addrows(1)
        tb.putcell('NUM_CORR', 0, 1)
        tb.putcell('CORR_TYPE', 0, corr_type)
        tb.putcell('CORR_PRODUCT', 0, corr_product)


def __create_processor_table(params):
    __create_subtable(params.outfile, "PROCESSOR", EmptyMSBaseInformation.processor_desc, EmptyMSBaseInformation.processor_dminfo)


def __create_source_table(params):
    __create_subtable(params.outfile, "SOURCE", EmptyMSBaseInformation.source_desc, EmptyMSBaseInformation.source_dminfo)


def __create_special_window_table(params):
    __create_subtable(params.outfile, "SPECTRAL_WINDOW", EmptyMSBaseInformation.special_window_desc, EmptyMSBaseInformation.special_window_dminfo)
    with table_manager(os.path.join(params.outfile, 'SPECTRAL_WINDOW'), nomodify=False) as tb:
        cw = np.ones(params.im_nchan, dtype=float) * 1e6
        cf = 1e9 + np.arange(params.im_nchan, dtype=float) * 1e6

        if tb.nrows() == 0:
            tb.addrows(1)
        tb.putcell('NUM_CHAN', 0, params.im_nchan)
        tb.putcell('CHAN_FREQ', 0, cf)
        tb.putcell('CHAN_WIDTH', 0, cw)
        tb.putcell('RESOLUTION', 0, cw)
        tb.putcell('EFFECTIVE_BW', 0, cw)
        tb.putcell('TOTAL_BANDWIDTH', 0, cw.sum())


def __create_state_table(params):
    __create_subtable(params.outfile, "STATE", EmptyMSBaseInformation.state_desc, EmptyMSBaseInformation.state_dminfo)
    with table_manager(os.path.join(params.outfile, 'STATE'), nomodify=False) as tb:
        if tb.nrows() == 0:
            tb.addrows(1)
        tb.putcell('OBS_MODE', 0, 'OBSERVE_TARGET#ON_SOURCE_IMAGE_DOMAIN')


def __check_ms_path(params, overWrite: bool = True):
    exists = os.path.exists(params.outfile)
    if overWrite and exists:
        shutil.rmtree(params.outfile)
    return exists


def __create_subtable(outfile, subtable: str, desc: str, dminfo: str):
    tb = table()
    try:
        tb.create(f"{outfile}/{subtable}", desc, dminfo=dminfo)
    finally:
        tb.close()
    with table_manager(outfile, nomodify=False) as tb:
        tb.putkeyword(subtable, f"Table: {outfile}/{subtable}")


def __copy_image_array_to_ms(params):
    __put_image_data_to_ms(params, *__get_image_value(params))


def __get_image_value(params):
    # get image array and mask from the image
    with tool_manager(params.infile, image) as _ia:
        arr = _ia.getchunk()
        msk = _ia.getchunk(getmask=True)

    # put image array slices to MS DATA column
    # axis indices for spatial, spectral and polarization axes
    xax, yax = params.axis_dir
    spax = params.axis_sp
    if params.axis_pol > 0:
        polax = params.axis_pol
    else:
        arr = np.expand_dims(arr, axis=3)
        msk = np.expand_dims(msk, axis=3)
        polax = 3
    casalog.post(f'axis index: {xax} {yax} {polax} {spax}', 'DEBUG2')
    return arr, msk, xax, yax, spax, polax


def __put_image_data_to_ms(params, image_array, mask_array, axis_x, axis_y, axis_sp, axis_pol):
    # which data column to use
    with table_manager(params.outfile, nomodify=False) as tb:
        # also set FLAG, SIGMA, WEIGHT, and UVW
        index_list = [0, 0, 0, 0]
        index_list[axis_sp] = np.arange(params.im_nchan)
        nx, ny = params.dir_shape
        irow = 0
        wgt = np.ones(1, dtype=float)
        uvw = np.zeros(3, dtype=float)
        for ix in range(nx):
            index_list[axis_x] = ix
            for iy in range(ny):
                index_list[axis_y] = iy
                for ip in range(params.im_npol):
                    index_list[axis_pol] = ip
                    slice = tuple(index_list)
                    cell = image_array[slice]
                    mask = np.logical_not(mask_array[slice])
                    casalog.post(f'slice={slice}, cell={cell}', 'DEBUG2')
                    tb.putcell(params.datacolumn, irow, np.expand_dims(cell, axis=0))
                    tb.putcell('FLAG', irow, np.expand_dims(mask, axis=0))
                    tb.putcell('SIGMA', irow, wgt)
                    tb.putcell('WEIGHT', irow, wgt)
                    tb.putcell('UVW', irow, uvw)
                    irow += 1


class MS2ImageParams(Validable):

    def __init__(self, infile: str=None, linefile: str='', imagefile: str='', datacolumn: str='DATA', output_cont: bool=False,
                 output_cont_file: str='', i2ms: Image2MSParams=None) -> None:
        self.infile = infile
        self.linefile = linefile
        self.imagename = imagefile
        self.datacolumn = datacolumn
        self.output_cont = output_cont
        self.output_cont_file = output_cont_file
        self.i2ms = i2ms
        self.validate()

    def validate(self):
        pass


def ms2image(params: MS2ImageParams):
    casalog.post("start imaging", "DEBUG2")
    try:
        __copy_image_file(params.imagename, params.linefile)  # mask data is also copied in this method
        if params.output_cont:
            __copy_image_file(infile=params.imagename, outfile=params.output_cont_file)
    except Exception as instance:
        casalog.post(f"*** Error '{instance}'", 'SEVERE')

    __make_image_array(params)
    casalog.post("end arraying", "DEBUG2")

    if params.output_cont:
        __output_cont_image(params)
    __output_image(params)
    casalog.post("end imaging", "DEBUG2")


def __copy_image_file(infile: str=None, outfile: str=None):
    if not os.path.exists(infile):
        raise Exception(f'Image files not found, infile:{infile}')
    try:
        _ia = image()
        ok = _ia.fromimage(infile=infile, outfile=outfile)
        if not ok:
            raise Exception(f'Some error occured, infile:{infile}, outfile:{outfile}')
    finally:
        _ia.done()


def __make_image_array(params):
    nx, ny = params.i2ms.dir_shape
    params.array = np.empty((nx, ny, params.i2ms.im_nchan))
    pos = 0
    with table_manager(params.infile) as tb:
        for i in range(nx):
            for j in range(ny):
                params.array[i][j] = tb.getcell(params.datacolumn, pos)[0].real
                pos += 1
    if params.i2ms.axis_pol > 0:
        params.array = np.expand_dims(params.array, params.i2ms.axis_pol)


def __output_image(params):
    with tool_manager(params.linefile, image) as _ia:
        _ia.putchunk(pixels=params.array, locking=True)
        try:
            param_names = imbaseline.__code__.co_varnames[:imbaseline.__code__.co_argcount]
            vars = locals()
            param_vals = [vars[p] for p in param_names]
            write_image_history(_ia, sys._getframe().f_code.co_name, param_names, param_vals, casalog)
        except Exception as instance:
            casalog.post(f"*** Error '{instance}' updating HISTORY", 'WARN')


def __output_cont_image(params):
    with tool_manager(params.output_cont_file, image) as _ia:
        _ia.putchunk(pixels=_ia.getchunk() - params.array, locking=True)


class EmptyMSBaseInformation:
    """The Parameters class for creating an empty MeasurementSet.

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
               'SPEC': {'DEFAULTTILESHAPE': array([2, 63, 1040]),
                        'HYPERCUBES': {'*1': {'BucketSize': 1048320,
                                              'CellShape': array([2, 63]),
                                              'CubeShape': array([2, 63, 22653]),
                                              'ID': {},
                                              'TileShape': array([2, 63, 1040])}},
                        'IndexSize': 1,
                        'MAXIMUMCACHESIZE': 0,
                        'MaxCacheSize': 0,
                        'SEQNR': 2},
               'TYPE': 'TiledShapeStMan'},
        '*4': {'COLUMNS': array(['FLAG'], dtype='<U16'),
               'NAME': 'TiledFlag',
               'SEQNR': 3,
               'SPEC': {'DEFAULTTILESHAPE': array([2, 63, 1040]),
                        'HYPERCUBES': {'*1': {'BucketSize': 16380,
                                              'CellShape': array([2, 63]),
                                              'CubeShape': array([2, 63, 22653]),
                                              'ID': {},
                                              'TileShape': array([2, 63, 1040])}},
                        'IndexSize': 1,
                        'MAXIMUMCACHESIZE': 0,
                        'MaxCacheSize': 0,
                        'SEQNR': 3},
               'TYPE': 'TiledShapeStMan'},
        '*5': {'COLUMNS': array(['FLAG_CATEGORY'], dtype='<U16'),
               'NAME': 'TiledFlagCategory',
               'SEQNR': 4,
               'SPEC': {'DEFAULTTILESHAPE': array([2, 63, 1, 1040]),
                        'HYPERCUBES': {},
                        'IndexSize': 0,
                        'MAXIMUMCACHESIZE': 0,
                        'MaxCacheSize': 0,
                        'SEQNR': 4},
               'TYPE': 'TiledShapeStMan'},
        '*6': {'COLUMNS': array(['WEIGHT_SPECTRUM'], dtype='<U16'),
               'NAME': 'TiledWgtSpectrum',
               'SEQNR': 5,
               'SPEC': {'DEFAULTTILESHAPE': array([2, 63, 1040]),
                        'HYPERCUBES': {'*1': {'BucketSize': 524160,
                                              'CellShape': array([2, 63]),
                                              'CubeShape': array([2, 63, 22653]),
                                              'ID': {},
                                              'TileShape': array([2, 63, 1040])}},
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
                                         'keywords': {'MEASINFO': {'TabRefCodes': array([0, 1, 2, 3, 4, 5, 6, 7, 8, 64], dtype=uint64),
                                                                   'TabRefTypes': array(['REST', 'LSRK', 'LSRD', 'BARY', 'GEO', 'TOPO', 'GALACTO',
                                                                                         'LGROUP', 'CMB', 'Undefined'], dtype='<U16'),
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
                                             'keywords': {'MEASINFO': {'TabRefCodes': array([0, 1, 2, 3, 4, 5, 6, 7, 8, 64], dtype=uint64),
                                                                       'TabRefTypes': array(['REST', 'LSRK', 'LSRD', 'BARY', 'GEO', 'TOPO', 'GALACTO',
                                                                                             'LGROUP', 'CMB', 'Undefined'], dtype='<U16'),
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
