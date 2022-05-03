# image-based line finding and baseline subtraction.

from abc import abstractmethod
import contextlib
import os
import shutil
import sys
from typing import Any, Dict, List, Tuple, Union
import uuid

import numpy as np
from numpy import array, uint64

from casatasks import casalog
from casatasks.private.ialib import write_image_history
from casatasks.private.sdutil import (sdtask_decorator, table_manager,
                                      tool_manager)
from casatasks.private.task_imsmooth import imsmooth
from casatasks.private.task_sdbaseline import sdbaseline
from casatasks.private.task_sdsmooth import sdsmooth
from casatasks.private.update_spw import sets_to_spwchan, spwchan_to_sets
from casatools import image, quanta, table

DATACOLUMN = 'DATA'
OVERWRITE = True

IMAGE_STACK_MAX_HEIGHT = 5
MS_STACK_MAX_HEIGHT = 3

qa = quanta()
do_not_erase_temporary_files = False
dump_tasks = False


class AbstractFolder:
    """Abstract class has Image/MeasurementSet file path.

    This class and child classes are wrapper of CasaImage/MeasurementSet file.
    The wrapped path could be decided to erase by which child classes are implemented.
    """

    has_file = False

    def __init__(self, file: str = None) -> None:
        """Initialize a Folder object."""
        if not os.path.exists(file):
            raise ValueError(f'file {file} is not found')
        self.path = file
        self.has_file = True

    @abstractmethod
    def erase(self) -> None:
        """Erase the file pointed path."""
        raise RuntimeError('Not implemented')


class _EraseableFolder(AbstractFolder):
    """Image/MeasurementSet file path class. The file path is permitted to erase."""

    def __init__(self, file: str = None) -> None:
        super().__init__(file)
        _eraseable_folder_register.register(self)

    def erase(self, dry_run: bool = True) -> None:
        if self.has_file:
            if dry_run:
                casalog.post(f'[DRY RUN] erase file: {self.path}', 'DEBUG2')
            else:
                casalog.post(f'erase file: {self.path}', 'DEBUG2')
                if os.path.exists(self.path):
                    shutil.rmtree(self.path)
                    if not os.path.exists(self.path):
                        self.has_file = False
        else:
            casalog.post(f'not found the file to erase: {self.path}', 'WARN')


class _UnerasableFolder(AbstractFolder):
    """Image/MeasurementSet file path class. The file path is NOT permitted to erase."""

    def erase(self, dry_run: bool = True) -> None:
        casalog.post(f'un-erase file: {self.path}', 'DEBUG2')


class _EraseableFolderRegister():
    """Class of the register of folders that need to be erased."""

    _register = dict()

    def register(self, folder: _EraseableFolder):
        if isinstance(folder, _EraseableFolder) and folder.path not in self._register.keys():
            self._register[folder.path] = folder
        else:
            raise ValueError('Irregal folder would be appended', 'SEVERE')

    def clear(self, dry_run: bool = True):
        if not dry_run:
            for path, folder in self._register.items():
                if not folder.has_file:
                    raise RuntimeError('Invalid code execution state', 'SEVERE')
                elif not os.path.exists(path):
                    raise RuntimeError(f'File not found: {path}', 'SEVERE')
                else:
                    folder.erase(dry_run)
        self._register.clear()
        casalog.post('cleaned up _EraseableFolderRegister', 'DEBUG2')

    def pop(self, key: str) -> _EraseableFolder:
        return self._register.pop(key)


_eraseable_folder_register = _EraseableFolderRegister()


class AbstractFileStack:
    """CasaImage/MeasurementSet file path stack to be processed by tasks in imbaseline.

    The paths of CasaImage or MeasurementSet are wrapped by AbstractFolder class.
    Implementation classes of AbstractFileStack are _EraseableFolder/Un_EraseableFolder, the
    _EraseableFolder class erases the path holden by a property 'path' when execute cleaning
    process, and the Un_EraseableFolder class doesn't erase it.
    If this class is used to stack a path of CasaImage, the bottom of it must be the input
    image(an argument 'imagename').
    """

    def __init__(self, top: AbstractFolder = None, max_height=None) -> None:
        """Initialize a FileStack object."""
        self.stack = []
        self.max_height = max_height
        if isinstance(top, AbstractFolder):
            self.push(top)

    def push(self, file: AbstractFolder = None) -> None:
        """Push a folder into the stack."""
        if not isinstance(file, AbstractFolder):
            raise ValueError(f'cannot append {file.path}')
        elif self.height() == self.max_height:
            raise RuntimeError('stack is full')
        else:
            casalog.post(f'push {file.path} into the stack', 'DEBUG2')
            self.stack.append(file)

    def pop(self) -> AbstractFolder:
        """Return and remove the top of the stack.

        The stack object should have input CasaImage or converted MeasurementSet at bottom,
        so if pop() is called when stack height is zero, then it raises RuntimeError.
        """
        if self.height() <= 1:
            raise RuntimeError('the stack cannot pop')
        return self.stack.pop()

    def peak(self) -> AbstractFolder:
        """Return a pointer of the top of the stack."""
        if len(self.stack) > 0:
            picked = self.stack[-1]
            casalog.post(f'pick from the stack: {picked.path}', 'DEBUG2')
            return self.stack[len(self.stack) - 1]
        else:
            raise RuntimeError('the stack is empty')

    def subpeak(self) -> AbstractFolder:
        """Return a pointer of a next of the top of the stack."""
        if len(self.stack) > 1:
            picked = self.stack[len(self.stack) - 2]
            casalog.post(f'pick from sub peak of the stack: {picked.path}', 'DEBUG2')
            return self.stack[len(self.stack) - 2]
        else:
            raise RuntimeError('the stack has only one stuff')

    def bottom(self) -> AbstractFolder:
        """Return a pointer of the bottom of the stack."""
        if len(self.stack) > 0:
            picked = self.stack[0]
            casalog.post(f'pick from bottom of the stack: {picked.path}', 'DEBUG2')
            return self.stack[0]
        else:
            raise RuntimeError('the stack has not have enough stuff')

    def clear(self, dry_run: bool = True) -> None:
        """Do erase method of all of the stack and clear the stack."""
        self.stack.clear()

    def height(self) -> int:
        """Return height of the stack."""
        return len(self.stack)


class _CasaImageStack(AbstractFileStack):
    """FileStack for CasaImage."""

    def __init__(self, top: AbstractFolder = None) -> None:
        super().__init__(top=top, max_height=IMAGE_STACK_MAX_HEIGHT)


class _MeasurementSetStack(AbstractFileStack):
    """FileStack for MeasurementSet."""

    def __init__(self) -> None:
        super().__init__(max_height=MS_STACK_MAX_HEIGHT)
        self.spsmoothed = False
        self.kwidth = 5


@contextlib.contextmanager
def _stack_manager(initial_image=None):
    image_stack = _CasaImageStack(top=_UnerasableFolder(initial_image))
    ms_stack = _MeasurementSetStack()
    try:
        yield image_stack, ms_stack
    finally:
        image_stack.clear()
        ms_stack.clear()


class AbstractValidatable:
    """Abstract class to enforce implementing validation method.

    A class that inherits this class is that require parameter validation of it.
    """

    @abstractmethod
    def validate(self) -> None:
        """Validate the object."""
        raise RuntimeError('Not implemented')


class _ImageShape(AbstractValidatable):
    """Shape parameters of input image.

    These parameters are been getting in Image2MS, using in MS2Image.
    """

    def __init__(self, im_shape: np.ndarray = None, axis_dir: np.ndarray = None,
                 axis_sp: int = None, axis_pol: int = None) -> None:
        self.im_shape = im_shape
        self.axis_dir = axis_dir
        self.axis_sp = axis_sp
        self.axis_pol = axis_pol
        self.dir_shape = self.im_shape[self.axis_dir]
        self.im_nrow = np.prod(self.dir_shape)
        self.im_nchan = self.im_shape[self.axis_sp] if self.axis_sp > 0 else 1
        self.im_npol = self.im_shape[self.axis_pol] if self.axis_pol > 0 else 1

    def validate(self) -> None:
        if not len(self.axis_dir):
            raise ValueError(f'invalid value: axis_dir {self.axis_dir}')

        if not (self.axis_sp in [-1, 2, 3] and self.axis_pol in [-1, 2, 3]
                and self.axis_sp != self.axis_pol):
            raise ValueError(f'invalid value: sp:{self.axis_sp} or pol:{self.axis_pol}')

        # if im_nchan is too few, say, <10, sdbaseline should abort
        if self.im_nchan < 10:
            raise ValueError(f'nchan {self.im_nchan} is too few to perform baseline subtraction')

        if len(self.im_shape) < 3:
            raise ValueError(f'invalid value: im_shape {self.im_shape}')

        if len(self.dir_shape) < 2:
            raise ValueError(f'invalid value: dir_shape {self.dir_shape}')


def _get_image_shape(imagepath: str) -> _ImageShape:
    if not os.path.exists(imagepath):
        raise ValueError(f"path '{imagepath}' is not found")

    shape = None
    with tool_manager(imagepath, image) as ia:
        try:
            cs = ia.coordsys()
            shape = _ImageShape(
                ia.shape(),
                cs.findcoordinate('direction')['world'],
                __get_axis_position(cs.findcoordinate('spectral')['world']),  # 3 or 2 or -1
                __get_axis_position(cs.findcoordinate('stokes')['world'])   # 2 or 3 or -1
            )
        finally:
            cs.done()

    if shape.im_shape.shape[0] < 3:
        raise ValueError(f"image '{imagepath}' is invalid")

    shape.validate()
    casalog.post(f'image shape is {shape.im_shape}, direciton {shape.dir_shape} '
                 f'({shape.im_nrow} pixels), npol {shape.im_npol}, nchan {shape.im_nchan}',
                 'DEBUG2')

    return shape


def __get_axis_position(val: array = None) -> int:
    if val is not None and len(val) > 0:
        return val[0].item()
    return -1


@sdtask_decorator
def imbaseline(imagename=None, linefile=None, output_cont=None, bloutput=None, maskmode=None,
               chans=None, thresh=None, avg_limit=None, minwidth=None, edge=None, blfunc=None,
               order=None, npiece=None, applyfft=None, fftthresh=None, addwn=None, rejwn=None,
               blparam=None, clipniter=None, clipthresh=None, dirkernel=None, major=None,
               minor=None, pa=None, kimage=None, scale=None, spkernel=None, kwidth=None) -> None:
    """Execute imbaseline.

    All specifications of arguments are defined in:
    https://open-jira.nrao.edu/browse/CAS-13520

    The task executes several processes as follows:
    (1) do direction plane smoothing of input casa image (execute imsmooth)
    (2) convert casa image into MeasurementSet (a pixel of image corresponds to a line of MS)
    (3) do spectral smoothing of MS (execute sdsmooth)
    (4) do baselining (execute sdbaseline)
    (5) convert MS into casa image, and subtract results
    """
    _validate_imagename(imagename)
    linefile = _prepare_linefile(linefile, imagename)

    with _stack_manager(imagename) as (image_stack, ms_stack):
        try:
            input_image_shape = _get_image_shape(image_stack.peak().path)

            # do direction plane smoothing
            _ImsmoothMethods.execute(dirkernel, major, minor, pa, kimage, scale, image_stack)

            # convert casaimage into MeasurementSet
            _Image2MSMethods.execute(DATACOLUMN, input_image_shape, image_stack, ms_stack)

            # do spectral smoothing
            _SdsmoothMethods.execute(DATACOLUMN, spkernel, kwidth, image_stack, ms_stack,
                                     input_image_shape)

            # do baselining
            _SdbaselineMethods.execute(DATACOLUMN, bloutput, maskmode, chans, thresh, avg_limit,
                                       minwidth, edge, blfunc, order, npiece, applyfft, fftthresh,
                                       addwn, rejwn, blparam, clipniter, clipthresh, image_stack,
                                       ms_stack, input_image_shape)

            # convert MeasurementSet into image and subtract results
            _ImageSubtractionMethods.execute(linefile, image_stack)

            if output_cont:
                _ImageSubtractionMethods.get_continuum_image(image_stack)
        finally:
            _do_post_processing(linefile)


def _validate_imagename(imagename: str = None) -> None:
    if not os.path.exists(imagename):
        raise ValueError(f'Error: file {imagename} is not found.', 'SEVERE')


def _prepare_linefile(linefile: str = None, imagename: str = None) -> str:
    if linefile == '' or linefile is None:
        linefile = os.path.basename(imagename).rstrip('/') + '_bs'
    if not OVERWRITE and os.path.exists(linefile):
        raise ValueError(f'Error: file {linefile} already exists, please delete before continuing.',
                         'SEVERE')
    return linefile


def _generate_temporary_filename(prefix: str = '', ext: str = '') -> str:
    if prefix and prefix[-1] != '-':
        prefix = prefix + '-'
    if ext != '':
        ext = '.' + ext
    while True:
        filename = prefix + str(uuid.uuid4()) + ext
        if not os.path.exists(filename):
            return filename


def _copy_image_file(infile: str = None, outfile: str = None) -> None:
    if not os.path.exists(infile):
        raise Exception(f'Image files not found, infile: {infile}')

    ia = image()
    try:
        ok = ia.fromimage(infile=infile, outfile=outfile)
        if not ok:
            raise Exception(f'Some error occured, infile: {infile}, outfile: {outfile}')
    finally:
        ia.done()


def _do_post_processing(outfile) -> None:
    """Execute some post-processes of imbaseline."""
    _eraseable_folder_register.clear(dry_run=do_not_erase_temporary_files)
    __write_image_history(outfile)


def __write_image_history(outfile) -> None:
    with tool_manager(outfile, image) as outia:
        try:
            param_names = imbaseline.__code__.co_varnames[:imbaseline.__code__.co_argcount]
            vars = locals()
            param_vals = [vars[p] for p in param_names]
            write_image_history(outia, sys._getframe().f_code.co_name, param_names,
                                param_vals, casalog)
        except Exception as instance:
            casalog.post(f'*** Error "{instance}" updating HISTORY', 'WARN')


class _ImageSubtractionMethods():

    @staticmethod
    def execute(linefile: str = None, image_stack: AbstractFileStack = None) -> None:
        """Execute image subtraction.

        The output image is computed as follows:
        - in case any smoothing (along direction and/or frequency) is executed
            output_image = input_image - (smoothed_image - smoothed_and_subtracted_image)
        - in case no smoothing is executed
            output_image = subtracted_image
        """
        if image_stack.height() <= 2:  # any smoothing were not executed
            output_image = image_stack.pop().path
            _eraseable_folder_register.pop(output_image)
            os.rename(output_image, linefile)
            image_stack.push(_UnerasableFolder(linefile))
        else:
            smoothed_image = image_stack.subpeak().path
            subtracted_image = image_stack.peak().path
            base_image = image_stack.bottom().path
            _copy_image_file(base_image, linefile)
            _ImageSubtractionMethods.__subtract_image(smoothed_image, subtracted_image)
            _ImageSubtractionMethods.__subtract_image(linefile, smoothed_image)
            image_stack.push(_UnerasableFolder(linefile))

    @staticmethod
    def get_continuum_image(image_stack: AbstractFileStack = None) -> None:
        """Compute 'input_image - output_image'."""
        base_image = image_stack.bottom().path
        output_image = os.path.basename(base_image) + '.cont'
        _copy_image_file(base_image, output_image)

        linefile = image_stack.peak().path
        _ImageSubtractionMethods.__subtract_image(output_image, linefile)

    @staticmethod
    def __subtract_image(operand_a: str = None, operand_b: str = None) -> None:
        """Subtract image chunk."""
        image_array = None
        with tool_manager(operand_b, image) as ia:
            image_array = ia.getchunk()

        with tool_manager(operand_a, image) as ia:
            ia.putchunk(pixels=ia.getchunk() - image_array, locking=True)


class _ImsmoothMethods():
    """Methoods for Imsmooth execution."""

    @staticmethod
    def execute(dirkernel: str = None, major: str = None, minor: str = None, pa: str = None,
                kimage: str = None, scale: float = None, stack: AbstractFileStack = None) -> None:
        """Call casatasks.imsmooth task if dirkernel is specified."""
        if not _ImsmoothMethods.require(dirkernel):
            casalog.post('omit image smoothing', 'INFO')
            return

        casalog.post('execute image smoothing', 'INFO')
        infile = stack.peak().path
        outfile = _generate_temporary_filename('dirsmooth', 'im')
        imsmooth(**_ImsmoothParams(infile, outfile, dirkernel, major, minor, pa, kimage, scale)())
        stack.push(_EraseableFolder(outfile))

    @staticmethod
    def require(dirkernel: str = 'none') -> None:
        if not dirkernel:
            dirkernel = 'none'
        if dirkernel == 'none':
            return False
        elif dirkernel in ['image', 'boxcar', 'gaussian']:
            return True
        else:
            raise ValueError(f'Unsupported direction smoothing kernel, {dirkernel}', 'SEVERE')


class _SdsmoothMethods():
    """Methoods for Sdsmooth execution."""

    @staticmethod
    def execute(datacolumn: str = None, spkernel: str = None, kwidth: int = None,
                image_stack: AbstractFileStack = None, ms_stack: AbstractFileStack = None,
                image_shape: _ImageShape = None) -> None:
        """Call casatasks.sdsmooth task if spkernel is specified."""
        if not _SdsmoothMethods.require(spkernel):
            casalog.post('omit spectral smoothing', 'INFO')
            return

        casalog.post('execute spectral smoothing', 'INFO')

        input_ms = ms_stack.peak().path
        output_ms = _generate_temporary_filename('spsmooth', 'ms')
        base_image = image_stack.bottom().path
        params = _SdsmoothParams(input_ms, output_ms, datacolumn.lower(), spkernel, kwidth)
        sdsmooth(**params())
        ms_stack.push(_EraseableFolder(output_ms))
        output_image = _MS2ImageMethods.convert(base_image, output_ms, image_shape, datacolumn)
        image_stack.push(_EraseableFolder(output_image))
        ms_stack.spsmoothed = True
        ms_stack.kwidth = params.kwidth

    @staticmethod
    def require(spkernel: str = 'none') -> None:
        if not spkernel:
            spkernel = 'none'
        if spkernel == 'none':
            return False
        elif spkernel in ['boxcar', 'gaussian']:
            return True
        else:
            raise ValueError(f'Unsupported spectral smoothing kernel, {spkernel}', 'SEVERE')


class _SdbaselineMethods():
    """Methoods for Sdbaseline execution."""

    @staticmethod
    def execute(datacolumn: str = None, bloutput: str = None, maskmode: str = None,
                chans: str = None, thresh: float = None, avg_limit: int = None,
                minwidth: int = None, edge: List[int] = None,
                blfunc: str = None, order: int = None, npiece: int = None, applyfft: bool = None,
                fftthresh: float = None, addwn: List[int] = None, rejwn: List[int] = None,
                blparam: str = None, clipniter: int = None, clipthresh: float = None,
                image_stack: AbstractFileStack = None, ms_stack: AbstractFileStack = None,
                image_shape: _ImageShape = None) -> None:
        """Call casatasks.sdbaseline task."""
        casalog.post('execute spectral baselining', 'INFO')
        input_ms = ms_stack.peak().path
        output_ms = _generate_temporary_filename('baseline', 'ms')
        base_image = image_stack.bottom().path
        sdbaseline(**_SdbaselineParams(input_ms, output_ms, datacolumn.lower(), bloutput, maskmode,
                                       chans, thresh, avg_limit, minwidth, edge, blfunc, order,
                                       npiece, applyfft, fftthresh, addwn, rejwn, blparam,
                                       clipniter, clipthresh, ms_stack.spsmoothed, image_shape,
                                       ms_stack.kwidth)())
        ms_stack.push(_EraseableFolder(output_ms))
        output_image = _MS2ImageMethods.convert(base_image, output_ms, image_shape, datacolumn)
        image_stack.push(_EraseableFolder(output_image))
        blparam_name = input_ms + '_blparam.' + _SdbaselineParams.FIXED_PARAM['blformat']
        if os.path.exists(blparam_name):
            _SdbaselineMethods.__rename_blparam_filename(blparam_name, base_image)

    @staticmethod
    def __rename_blparam_filename(filename: str = None, basename: str = None) -> str:
        if not os.path.exists(filename):
            return None
        newname = os.path.basename(basename) + '.ms_blparam.' + \
            _SdbaselineParams.FIXED_PARAM['blformat']
        if os.path.exists(newname):
            return filename
        try:
            os.rename(filename, newname)
        except Exception:
            casalog.post(f'rename failure:from {filename} to {newname}', 'SEVERE')
            return filename
        return newname


class _ImsmoothParams(AbstractValidatable):
    """Parameter manipulation class for execution of casatasks.imsmooth."""

    FIXED_PARAM = dict(
        targetres=False,
        mask='',
        beam={},
        region='',
        box='',
        chans='',
        stokes='',
        stretch=False,
        overwrite=True
    )

    def __init__(self, infile: str = None, outfile: str = None, dirkernel: str = 'none',
                 major: str = '', minor: str = '', pa: str = '', kimage: str = '',
                 scale: int = -1.0) -> None:
        self.infile = infile
        self.outfile = outfile

        # dirkernel options: none(default)/gaussian/boxcar/image
        self.kernel = dirkernel if dirkernel is not None else 'none'

        # subparameter for dirkernel = gaussian/boxcar
        self.major = major if major is not None else ''
        self.minor = minor if minor is not None else ''
        self.pa = pa if pa is not None else ''

        # subparameter for dirkernel = image
        self.kimage = kimage if kimage is not None else ''
        self.scale = scale if scale is not None else -1.0

        self.validate()

    def validate(self) -> None:
        self.__validate_dirkernel()

    def __validate_dirkernel(self) -> None:
        if self.kernel == 'image':
            self.major = self.minor = self.pa = ''
            if self.kimage != '' and not os.path.exists(self.kimage):
                raise ValueError(f'Error: file {self.kimage} is not found.', 'SEVERE')
        elif self.kernel == 'gaussian' or self.kernel == 'boxcar':
            self.kimage = ''
            self.scale = -1.0
        elif self.kernel != 'none':
            raise ValueError(f'Unsupported maskmode, {self.kernel}', 'SEVERE')

    def __call__(self) -> Union[List[Any], Dict[str, str]]:
        """Convert the class into arguments of imsmooth().

        __log_origin is for callabletask.log_origin_setter
        """
        retval = dict(self.FIXED_PARAM, imagename=self.infile, kernel=self.kernel, major=self.major,
                      minor=self.minor, pa=self.pa, kimage=self.kimage, scale=self.scale,
                      outfile=self.outfile, __log_origin='imbaseline')
        if dump_tasks:
            print(_dump_tasks('imsmooth', retval))
        return retval


class _SdsmoothParams(AbstractValidatable):
    """Parameter manipulation class for execution of casatasks.sdsmooth."""

    FIXED_PARAM = dict(
        spw='',
        field='',
        antenna='',
        timerange='',
        scan='',
        pol='',
        intent='',
        reindex=True,
        overwrite=True
    )

    def __init__(self, infile: str = None, outfile: str = None, datacolumn: str = None,
                 spkernel: str = 'none', kwidth: int = 5) -> None:
        self.infile = infile
        self.outfile = outfile
        self.datacolumn = datacolumn
        self.kernel = spkernel if spkernel is not None else 'none'   # none(default)/gaussian/boxcar
        self.kwidth = kwidth if kwidth is not None else 5            # gaussian/boxcar
        self.validate()

    def validate(self) -> None:
        self.__validate_spkernel()

    def __validate_spkernel(self) -> None:
        if self.kernel == 'none':
            self.kwidth = 5
        elif not (self.kernel == 'gaussian' or self.kernel == 'boxcar'):
            raise ValueError(f'Unsupported maskmode, {self.kernel}', 'SEVERE')

    def __call__(self) -> Dict[str, Any]:
        """Convert the class into arguments of sdsmooth().

        __log_origin is for sdutil.callabletask_decorator.
        """
        retval = dict(self.FIXED_PARAM, infile=self.infile, datacolumn=self.datacolumn,
                      kernel=self.kernel, kwidth=self.kwidth, outfile=self.outfile,
                      __log_origin='imbaseline')
        if dump_tasks:
            print(_dump_tasks('sdsmooth', retval))
        return retval


class _SdbaselineParams(AbstractValidatable):
    """Parameter manipulation class for execution of casatasks.sdbaseline."""

    FIXED_PARAM = dict(
        antenna='',
        field='',
        timerange='',
        scan='',
        pol='',
        intent='',
        reindex=True,
        blmode='fit',
        dosubtract=True,
        blformat='csv',
        bltable='',
        updateweight=False,
        sigmavalue='stddev',
        showprogress=False,
        minnrow=1000,
        fftmethod='fft',
        verbose=False,
        overwrite=True
    )

    def __init__(self, infile: str = None, outfile: str = None, datacolumn: str = None,
                 bloutput: str = '', maskmode: str = 'list', chans: str = '', thresh: float = 5.0,
                 avg_limit: int = 4, minwidth: int = 4, edge: List[int] = [0, 0],
                 blfunc: str = 'poly', order: int = 5, npiece: int = 3, applyfft: bool = True,
                 fftthresh: float = 3.0, addwn: List = [0],
                 rejwn: List = [], blparam: str = '', clipniter: int = 0, clipthresh: float = 3.0,
                 spsmoothed: bool = False, image_shape: _ImageShape = None,
                 kwidth: int = None) -> None:
        self.infile = infile
        self.outfile = outfile
        self.datacolumn = datacolumn
        self.bloutput = bloutput if bloutput is not None else ''

        # maskmode: list(default)/auto
        self.maskmode = maskmode.lower() if maskmode is not None else 'list'

        # subparam for maskmode = list
        self.spw = self.__chans2spw(chans, self.maskmode)

        # subparam for maskmode = auto
        self.thresh = thresh if thresh is not None else 5.0
        self.avg_limit = avg_limit if avg_limit is not None else 4
        self.minwidth = minwidth if minwidth is not None else 4
        self.edge = edge if edge is not None else [0, 0]

        # poly(default)/chebyshev/cspline/sinusoid/variable
        self.blfunc = blfunc.lower() if blfunc is not None else 'poly'

        # subparam for blfunc = poly/chebyshev
        self.order = order if order is not None else 5

        # subparam for blfunc = cspline
        self.npiece = npiece if npiece is not None else 3

        # subparam for blfunc = sinusoid
        self.applyfft = applyfft if applyfft is not None else True
        self.fftthresh = fftthresh if fftthresh is not None else 3.0
        self.addwn = addwn if addwn is not None else [0]
        self.rejwn = rejwn if rejwn is not None else []

        # subparam for blfunc = variable
        self.blparam = blparam if blparam is not None else ''
        self.clipniter = clipniter if clipniter is not None else 0
        self.clipthresh = clipthresh if clipthresh is not None else 3.0

        if spsmoothed:
            spw = '0'
            left_edge = kwidth // 2
            right_edge = image_shape.im_nchan - kwidth // 2
            if self.spw:
                coverage = set(list(range(left_edge, right_edge)))
                sets = dict()
                for i, s in spwchan_to_sets(self.infile, self.spw).items():
                    sets[i] = s & coverage
                spw = sets_to_spwchan(sets)
            else:
                spw = f'0:{left_edge}~{right_edge-1}'
            self.spw = spw
        if not self.spw:
            self.spw = '0'

    def __chans2spw(self, chans: str, maskmode) -> str:
        if not chans or maskmode != 'list':
            return ''
        return chans

    def validate(self) -> None:
        self.__validate_maskmode()
        self.__validate_blfunc()

    def __validate_maskmode(self) -> None:
        if not (self.maskmode == 'list' or self.maskmode == 'auto'):
            raise ValueError(f'Unsupported maskmode, {self.maskmode}', 'SEVERE')

    def __validate_blfunc(self) -> None:
        def is_valid_blfunc(self) -> None:
            return self.blfunc == 'poly' or self.blfunc == 'chebyshev' or self.blfunc == 'cspline' \
                or self.blfunc == 'sinusoid' or self.blfunc == 'variable'

        if not is_valid_blfunc(self):
            raise ValueError(f'Unsupported blfunc, {self.blfunc}', 'SEVERE')

        if self.blfunc == 'variable' and not os.path.exists(self.blparam):
            raise ValueError(f'input file {self.blparam} does not exists', 'SEVERE')

    def __call__(self) -> Dict[str, Any]:
        """Convert the class into arguments of sdbaseline().

        __log_origin is for sdutil.callabletask_decorator.
        """
        retval = dict(self.FIXED_PARAM, infile=self.infile, datacolumn=self.datacolumn,
                      maskmode=self.maskmode, thresh=self.thresh, avg_limit=self.avg_limit,
                      minwidth=self.minwidth, edge=self.edge, bloutput=self.bloutput,
                      blfunc=self.blfunc, order=self.order, npiece=self.npiece,
                      applyfft=self.applyfft, fftthresh=self.fftthresh, addwn=self.addwn,
                      rejwn=self.rejwn, clipthresh=self.clipthresh, clipniter=self.clipniter,
                      blparam=self.blparam, outfile=self.outfile, spw=self.spw,
                      __log_origin='imbaseline')
        if dump_tasks:
            print(_dump_tasks('sdbaseline', retval))

        return retval


class _Image2MSParams(AbstractValidatable):
    """Parameter manipulation class for executing image2ms()."""

    def __init__(self, infile: str = None, outfile: str = None, datacolumn: str = 'DATA',
                 input_image_shape: _ImageShape = None) -> None:
        self.infile = infile
        self.outfile = outfile
        self.datacolumn = datacolumn
        self.im_shape = input_image_shape.im_shape
        self.axis_dir = input_image_shape.axis_dir
        self.axis_sp = input_image_shape.axis_sp
        self.axis_pol = input_image_shape.axis_pol
        self.dir_shape = input_image_shape.dir_shape
        self.im_nrow = input_image_shape.im_nrow
        self.im_nchan = input_image_shape.im_nchan
        self.im_npol = input_image_shape.im_npol
        self.validate()

    def validate(self) -> None:
        self.__validate_outfile()

    def __validate_outfile(self) -> None:
        if os.path.exists(self.outfile):
            raise ValueError(f'Folder exists: {self.outfile}')


class _Image2MSMethods():
    """Methods for converting image to MeasurementSet."""

    @staticmethod
    def execute(datacolumn: str = None, input_image_shape: _ImageShape = None,
                image_stack: AbstractFileStack = None,
                ms_stack: AbstractFileStack = None) -> None:
        """Convert a casaimage to a MeasurementSet."""
        casalog.post('convert casaimage to MeasurementSet', 'INFO')
        infile = image_stack.peak().path
        outfile = _generate_temporary_filename('img2ms', 'ms')
        _Image2MSMethods.__image2ms(_Image2MSParams(infile, outfile, datacolumn, input_image_shape))
        ms_stack.push(_EraseableFolder(outfile))

    @staticmethod
    def __image2ms(params: _Image2MSParams = None) -> None:
        """Convert CasaImage into MeasurementSet."""
        _Image2MSMethods.__create_empty_ms(params)
        _Image2MSMethods.__put_parametes_from_image_to_ms(params)

    @staticmethod
    def __create_empty_ms(params: _Image2MSParams = None) -> None:
        _Image2MSMethods.__cleanup_ms_path(params)
        _Image2MSMethods.__create_maintable(params)
        _Image2MSMethods.__create_antenna_table(params)
        _Image2MSMethods.__create_data_description_table(params)
        _Image2MSMethods.__create_feed_table(params)
        _Image2MSMethods.__create_field_table(params)
        _Image2MSMethods.__create_flag_cmd_table(params)
        _Image2MSMethods.__create_history_table(params)
        _Image2MSMethods.__create_observation_table(params)
        _Image2MSMethods.__create_pointing_table(params)
        _Image2MSMethods.__create_polarization_table(params)
        _Image2MSMethods.__create_processor_table(params)
        _Image2MSMethods.__create_source_table(params)
        _Image2MSMethods.__create_special_window_table(params)
        _Image2MSMethods.__create_state_table(params)

    @staticmethod
    def __generate_time_list(nrow_req: int) -> np.ndarray:
        mjd_sec = qa.convert(qa.quantity('1995-04-13T09:19:00'), 's')['value']  # dummy timestamp
        interval = 30.0
        return mjd_sec + np.arange(nrow_req) * interval

    @staticmethod
    def __create_maintable(params: _Image2MSParams = None) -> None:
        tb = table()
        try:
            tb.create(params.outfile, _EmptyMSBaseInformation.ms_desc,
                      dminfo=_EmptyMSBaseInformation.ms_dminfo)
            tb.putkeyword(keyword='MS_VERSION', value=2)
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
            time_list = _Image2MSMethods.__generate_time_list(nrow_req)
            tb.putcol('TIME', time_list)
            casalog.post(f'number of rows {nrow}, number of image pixels {params.im_nrow}, '
                         f'number of pols {params.im_npol}, '
                         f'required rows {nrow_req}', 'DEBUG2')
        finally:
            tb.close()

    @staticmethod
    def __create_antenna_table(params: _Image2MSParams = None) -> None:
        _Image2MSMethods.__create_subtable(params.outfile,
                                           'ANTENNA',
                                           _EmptyMSBaseInformation.antenna_desc,
                                           _EmptyMSBaseInformation.antenna_dminfo)
        with table_manager(os.path.join(params.outfile, 'ANTENNA'), nomodify=False) as tb:
            tb.addrows(1)

    @staticmethod
    def __create_data_description_table(params: _Image2MSParams = None) -> None:
        _Image2MSMethods.__create_subtable(params.outfile,
                                           'DATA_DESCRIPTION',
                                           _EmptyMSBaseInformation.data_description_desc,
                                           _EmptyMSBaseInformation.data_description_dminfo)
        with table_manager(os.path.join(params.outfile, 'DATA_DESCRIPTION'), nomodify=False) as tb:
            tb.addrows(1)
            tb.putcell('SPECTRAL_WINDOW_ID', 0, 0)
            tb.putcell('POLARIZATION_ID', 0, 0)

    @staticmethod
    def __create_feed_table(params: _Image2MSParams = None) -> None:
        _Image2MSMethods.__create_subtable(params.outfile,
                                           'FEED',
                                           _EmptyMSBaseInformation.feed_desc,
                                           _EmptyMSBaseInformation.feed_dminfo)
        with table_manager(os.path.join(params.outfile, 'FEED'), nomodify=False) as tb:
            tb.addrows(1)

    @staticmethod
    def __create_field_table(params: _Image2MSParams = None) -> None:
        _Image2MSMethods.__create_subtable(params.outfile,
                                           'FIELD',
                                           _EmptyMSBaseInformation.field_desc,
                                           _EmptyMSBaseInformation.field_dminfo)
        with table_manager(os.path.join(params.outfile, 'FIELD'), nomodify=False) as tb:
            tb.addrows(1)
            tb.putcell('DELAY_DIR', 0, np.zeros((2, 1)))
            tb.putcell('PHASE_DIR', 0, np.zeros((2, 1)))
            tb.putcell('REFERENCE_DIR', 0, np.zeros((2, 1)))

    @staticmethod
    def __create_flag_cmd_table(params: _Image2MSParams = None) -> None:
        _Image2MSMethods.__create_subtable(params.outfile,
                                           'FLAG_CMD',
                                           _EmptyMSBaseInformation.flag_cmd_desc,
                                           _EmptyMSBaseInformation.flag_cmd_dminfo)

    @staticmethod
    def __create_history_table(params: _Image2MSParams = None) -> None:
        _Image2MSMethods.__create_subtable(params.outfile,
                                           'HISTORY',
                                           _EmptyMSBaseInformation.history_desc,
                                           _EmptyMSBaseInformation.history_dminfo)

    @staticmethod
    def __create_observation_table(params: _Image2MSParams = None) -> None:
        _Image2MSMethods.__create_subtable(params.outfile,
                                           'OBSERVATION',
                                           _EmptyMSBaseInformation.observation_desc,
                                           _EmptyMSBaseInformation.observation_dminfo)
        with table_manager(os.path.join(params.outfile, 'OBSERVATION'), nomodify=False) as tb:
            tb.addrows(1)

    @staticmethod
    def __create_pointing_table(params: _Image2MSParams = None) -> None:
        _Image2MSMethods.__create_subtable(params.outfile,
                                           'POINTING',
                                           _EmptyMSBaseInformation.pointing_desc,
                                           _EmptyMSBaseInformation.pointing_dminfo)

    @staticmethod
    def __create_polarization_table(params: _Image2MSParams = None) -> None:
        _Image2MSMethods.__create_subtable(params.outfile,
                                           'POLARIZATION',
                                           _EmptyMSBaseInformation.polarization_desc,
                                           _EmptyMSBaseInformation.polarization_dminfo)
        with table_manager(os.path.join(params.outfile, 'POLARIZATION'), nomodify=False) as tb:
            corr_type = np.ones(1, dtype=int)
            corr_product = np.ones(2, dtype=int).reshape((2, 1))
            if tb.nrows() == 0:
                tb.addrows(1)
            tb.putcell('NUM_CORR', 0, 1)
            tb.putcell('CORR_TYPE', 0, corr_type)
            tb.putcell('CORR_PRODUCT', 0, corr_product)

    @staticmethod
    def __create_processor_table(params: _Image2MSParams = None) -> None:
        _Image2MSMethods.__create_subtable(params.outfile,
                                           'PROCESSOR',
                                           _EmptyMSBaseInformation.processor_desc,
                                           _EmptyMSBaseInformation.processor_dminfo)

    @staticmethod
    def __create_source_table(params: _Image2MSParams = None) -> None:
        _Image2MSMethods.__create_subtable(params.outfile,
                                           'SOURCE',
                                           _EmptyMSBaseInformation.source_desc,
                                           _EmptyMSBaseInformation.source_dminfo)

    @staticmethod
    def __create_special_window_table(params: _Image2MSParams = None) -> None:
        _Image2MSMethods.__create_subtable(params.outfile,
                                           'SPECTRAL_WINDOW',
                                           _EmptyMSBaseInformation.special_window_desc,
                                           _EmptyMSBaseInformation.special_window_dminfo)
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

    @staticmethod
    def __create_state_table(params: _Image2MSParams = None) -> None:
        _Image2MSMethods.__create_subtable(params.outfile,
                                           'STATE',
                                           _EmptyMSBaseInformation.state_desc,
                                           _EmptyMSBaseInformation.state_dminfo)
        with table_manager(os.path.join(params.outfile, 'STATE'), nomodify=False) as tb:
            if tb.nrows() == 0:
                tb.addrows(1)
            tb.putcell('OBS_MODE', 0, 'OBSERVE_TARGET#ON_SOURCE_IMAGE_DOMAIN')

    @staticmethod
    def __cleanup_ms_path(params, overwrite: bool = True) -> None:
        exists = os.path.exists(params.outfile)
        if overwrite and exists:
            shutil.rmtree(params.outfile)
        return exists

    @staticmethod
    def __create_subtable(outfile: str = None, subtable: str = None,
                          desc: str = None, dminfo: str = None) -> None:
        tb = table()
        try:
            tb.create(f'{outfile}/{subtable}', desc, dminfo=dminfo)
        finally:
            tb.close()
        with table_manager(outfile, nomodify=False) as tb:
            tb.putkeyword(subtable, f'Table: {outfile}/{subtable}')

    @staticmethod
    def __put_parametes_from_image_to_ms(params: _Image2MSParams = None) -> None:
        _Image2MSMethods.__put_image_parameters_into_ms(
            params, *_Image2MSMethods.__get_image_parameters(params))

    @staticmethod
    def __get_image_parameters(params: _Image2MSParams = None) -> Tuple[np.array, int]:
        # get image array and mask from the image
        with tool_manager(params.infile, image) as ia:
            arr = ia.getchunk()
            msk = ia.getchunk(getmask=True)

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

    @staticmethod
    def __put_image_parameters_into_ms(params: _Image2MSParams, image_array: np.array,
                                       mask_array: np.array, axis_x: int,
                                       axis_y: int, axis_sp: int, axis_pol: int) -> None:
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


class _MS2ImageMethods():

    @staticmethod
    def convert(base_image: str = None, input_ms: str = None, input_image_shape: _ImageShape = None,
                datacolumn: str = None) -> None:
        output_image = _MS2ImageMethods.__change_file_extension(input_ms, 'im')
        _copy_image_file(base_image, output_image)  # mask data is also copied in this method

        image_array = _MS2ImageMethods.__make_image_array(input_image_shape, input_ms, datacolumn)
        _MS2ImageMethods.__output_image(output_image, image_array)

        return output_image

    @staticmethod
    def __change_file_extension(path: str = None, ext: str = None) -> str:
        if not os.path.exists(path):
            RuntimeError(f'cannot find path: {path}')
        base, oldext = os.path.splitext(path)
        new_path = base_path = f'{base}.{ext}'
        i = 0
        while True:
            if not os.path.exists(new_path):
                break
            new_path = base_path + f'.{i}'
            i += 1
            if i > 10:
                raise RuntimeError('some file system error occured')
        return new_path

    @staticmethod
    def __make_image_array(input_image_shape: _ImageShape = None, infile: str = None,
                           datacolumn: str = None) -> np.array:
        nx, ny = input_image_shape.dir_shape
        image_array = np.empty((nx, ny, input_image_shape.im_nchan))
        pos = 0
        with table_manager(infile) as tb:
            for i in range(nx):
                for j in range(ny):
                    image_array[i][j] = tb.getcell(datacolumn, pos)[0].real
                    pos += 1
        if input_image_shape.axis_pol > 0:
            image_array = np.expand_dims(image_array, input_image_shape.axis_pol)
        return image_array

    @staticmethod
    def __output_image(outfile: str = None, image_array: np.array = None) -> None:
        with tool_manager(outfile, image) as ia:
            ia.putchunk(pixels=image_array, locking=True)


def _dump_tasks(taskname: str, vals: dict):
    cmd = f'{taskname}('
    arr = []
    for key, val in vals.items():
        if key != '__log_origin':
            quote = '' if type(val) in (bool, int, float) else '\''
            arr.append(f'{key}={quote}{val}{quote}')
    cmd += ', '.join(arr)
    cmd += ')'
    return cmd


class _EmptyMSBaseInformation:
    """The Parameters class for creating an empty MeasurementSet.

    This class contains dictionaries to create an empty MS using table.create(),
    and it has no method. Dictionaries have two types; desc(desctiption) and
    dminfo(data management infomation), these are used as arguments of table.create(),
    and for a table creating, it needs a desc dict and a dminfo dict.
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
                          'keywords': {'CATEGORY':
                                       np.array(['FLAG_CMD', 'ORIGINAL', 'USER'], dtype='<U16')},
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
        '*1':
        {'COLUMNS': array(['FLAG_ROW', 'POLARIZATION_ID', 'SPECTRAL_WINDOW_ID'], dtype='<U19'),
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
        'REFERENCE_DIR': {
            'comment': 'Direction of REFERENCE center (e.g. RA, DEC).as polynomial in time.',
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
        'CORR_TYPE': {'comment':
                      'The polarization type for each correlation product, as a Stokes enum.',
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
        '*1': {'COLUMNS': array(['CORR_PRODUCT', 'CORR_TYPE',
                                 'FLAG_ROW', 'NUM_CORR'], dtype='<U16'),
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
        '*1': {'COLUMNS': array(['FLAG_ROW', 'MODE_ID',
                                 'SUB_TYPE', 'TYPE', 'TYPE_ID'], dtype='<U16'),
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
                      'keywords':
                          {'MEASINFO':
                              {'TabRefCodes': array([0, 1, 2, 3, 4, 5, 6, 7, 8, 64], dtype=uint64),
                               'TabRefTypes': array(['REST', 'LSRK', 'LSRD', 'BARY',
                                                     'GEO', 'TOPO', 'GALACTO',
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
        'REF_FREQUENCY': {
            'comment': 'The reference frequency',
            'dataManagerGroup': 'StandardStMan',
            'dataManagerType': 'StandardStMan',
            'keywords':
                {'MEASINFO':
                    {'TabRefCodes': array([0, 1, 2, 3, 4, 5, 6, 7, 8, 64], dtype=uint64),
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
        'TOTAL_BANDWIDTH': {
            'comment': 'The total bandwidth for this window',
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
