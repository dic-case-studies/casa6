from __future__ import absolute_import
import os
import numpy
import traceback
import functools
import re
import abc
import contextlib

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import table, calibrater, imager
    from casatools import ms as mstool
    from casatools.platform import bytes2str
    from casatasks import casalog
else:
    from taskinit import casalog, gentools
    # make CASA5 tools constructors look like CASA6 tools
    from taskinit import tbtool as table
    from taskinit import cbtool as calibrater
    from taskinit import imtool as imager
    from taskinit import mstool


@contextlib.contextmanager
def toolmanager(vis, ctor, *args, **kwargs):
    if is_CASA6:
        # this is the only syntax allowed in CASA6, code in CASA6 should be converted to
        # call this method with a tool constructor directly
        tool = ctor()
    else:
        # CASA5 code can invoke this with a tool name, shared CASA5 and CASA6 source
        # uses the CASA6 syntax - use callable to tell the difference
        if callable(ctor):
            tool = ctor()
        else:
            # assume the argument is string and use it to get the appropriate tool constructor
            # the original argument name here was 'tooltype'
            tool = gentools([ctor])[0]
    if vis:
        tool.open(vis, *args, **kwargs)
    try:
        yield tool
    finally:
        if "close" in dir(tool):
            tool.close()
        elif "done" in dir(tool):
            tool.done()


def tbmanager(vis, *args, **kwargs):
    return toolmanager(vis, table, *args, **kwargs)


def cbmanager(vis, *args, **kwargs):
    return toolmanager(vis, calibrater, *args, **kwargs)


def measuresmanager(*args, **kwargs):
    return toolmanager(None, measures, *args, **kwargs)


def is_ms(filename):
    if (os.path.isdir(filename) and os.path.exists(filename+'/table.info') and os.path.exists(filename+'/table.dat')):
        f = open(filename + '/table.info')
        if is_CASA6:
            l = bytes2str(f.readline())
        else:
            l = f.readline()
        f.close()
        if (l.find('Measurement Set') != -1):
            return True
        else:
            return False
    else:
        return False


@contextlib.contextmanager
def table_selector(table, taql, *args, **kwargs):
    with tbmanager(table, *args, **kwargs) as tb:
        tsel = tb.query(taql)
        try:
            yield tsel
        finally:
            tsel.close()


def sdtask_decorator(func):
    """
    This is a decorator function for sd tasks.
    Currently the decorator does:

       1) set origin to the logger
       2) handle exception

    So, you don't need to set origin in the task any more.
    Also, you don't need to write anything about error
    handling in the task. If you have something to do
    at the end of the task execution, those should be
    written in the destructor of worker class, not in
    the 'finally' block.
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # set origin
        casalog.origin(func.__name__)

        retval = None
        # Any errors are handled outside the task.
        # however, the implementation below is effectively
        # equivalent to handling it inside the task.
        try:
            # execute task
            retval = func(*args, **kwargs)
        except Exception as e:
            traceback_info = format_trace(traceback.format_exc())
            casalog.post(traceback_info, 'SEVERE')
            casalog.post(str(e), 'ERROR')
            raise
        return retval
    return wrapper


def format_trace(s):
    wexists = True
    regex = '.*sdutil\.py.*in wrapper.*'
    retval = s
    while wexists:
        ss = retval.split('\n')
        wexists = False
        for i in range(len(ss)):
            if re.match(regex, ss[i]):
                ss = ss[:i] + ss[i+2:]
                wexists = True
                break
        retval = '\n'.join(ss)
    return retval


class sdtask_manager(object):
    def __init__(self, cls, args):
        self.cls = cls
        self.args = args

    def __enter__(self):
        self.obj = self.cls(**self.args)
        return self.obj

    def __exit__(self, exc_type, exc_value, traceback):
        # explicitly call destructure to make sure it is called here
        self.obj.__del__()
        del self.obj
        if exc_type:
            return False
        else:
            return True


class sdtask_interface(object):
    """
    The sdtask_interface defines a common interface for sdtask_worker
    class. All worker classes can be used as follows:

       worker = sdtask_worker(**locals())
       worker.initialize()
       worker.execute()
       worker.finalize()
       del worker

    Derived classes must implement the above three methods: initialize(),
    execute(), and finalize().
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, **kwargs):
        for (k, v) in kwargs.items():
            setattr(self, k, v)
        # special treatment for selection parameters
        select_params = ['scan', 'pol', 'beam']
        for param in select_params:
            if hasattr(self, param):
                setattr(self, param+'no', getattr(self, param))
                # casalog.post("renaming self.%s -> self.%sno='%s'" % (param, param, getattr(self, param)))
                delattr(self, param)

    def __del__(self):
        pass

    @abc.abstractmethod
    def initialize(self):
        raise NotImplementedError('initialize is abstract method')

    @abc.abstractmethod
    def execute(self):
        raise NotImplementedError('execute is abstract method')

    @abc.abstractmethod
    def finalize(self):
        raise NotImplementedError('finalize is abstract method')


class sdtask_template_imaging(sdtask_interface):
    """
    The sdtask_template_imaging is a template class for worker
    class of imaging related sdtasks. It partially implement initialize()
    and finalize() using internal methods that must be implemented
    in the derived classes. For initialize(), derived classes
    must implement compile(), which sets up imaging parameters.
    You can implement paramter_check() to do any task specific parameter
    check in initialize().
    For finalize(), derived classes can implement cleanup().
    """
    def __init__(self, **kwargs):
        super(sdtask_template_imaging, self).__init__(**kwargs)
        self.is_table_opened = False
        self.is_imager_opened = False
        self.table = table()
        self.imager = imager()
        # workaround for sdtpimaging
        if not hasattr(self, 'infiles') and hasattr(self, 'infile'):
            self.infiles = [self.infile]

        self.__set_infiles()
        self.__set_subtable_name()

    def __del__(self, base=sdtask_interface):
        # table and imager must be closed when the instance
        # is deleted
        self.close_table()
        self.close_imager()
        self.cleanup()
        super(sdtask_template_imaging, self).__del__()

    def open_table(self, name, nomodify=True):
        if self.is_table_opened:
            casalog.post('Close table before re-open', priority='WARN')
            return
        self.table.open(name, nomodify=nomodify)
        self.is_table_opened = True

    def close_table(self):
        if self.is_table_opened:
            self.table.close()
        self.is_table_opened = False

    def open_imager(self, name=''):
        if self.is_imager_opened:
            casalog.post('Close imager before re-open', priority='WARN')
            return
        self.imager.open(name)
        self.is_imager_opened = True

    def close_imager(self):
        if self.is_imager_opened:
            self.imager.close()
        self.is_imager_opened = False

    def initialize(self):
        # infiles must be MS
        for idx in range(len(self.infiles)):
            if not is_ms(self.infiles[idx]):
                msg='input data sets must be in MS format'
                raise Exception(msg)

        self.parameter_check()
        self.compile()

    def finalize(self):
        pass

    def parameter_check(self):
        pass

    def compile(self):
        pass

    def cleanup(self):
        pass

    def __set_subtable_name(self):
        self.open_table(self.infiles[0])
        keys = self.table.getkeywords()
        self.close_table()
        self.field_table = get_subtable_name(keys['FIELD'])
        self.spw_table = get_subtable_name(keys['SPECTRAL_WINDOW'])
        self.source_table = get_subtable_name(keys['SOURCE'])
        self.antenna_table = get_subtable_name(keys['ANTENNA'])
        self.polarization_table = get_subtable_name(keys['POLARIZATION'])
        self.observation_table = get_subtable_name(keys['OBSERVATION'])
        self.pointing_table = get_subtable_name(keys['POINTING'])
        self.data_desc_table = get_subtable_name(keys['DATA_DESCRIPTION'])
        self.pointing_table = get_subtable_name(keys['POINTING'])

    def __set_infiles(self):
        if type(self.infiles) == str:
            self.infiles = [self.infiles]


def get_abspath(filename):
    return os.path.abspath(expand_path(filename))


def expand_path(filename):
    return os.path.expanduser(os.path.expandvars(filename))


def assert_outfile_canoverwrite_or_nonexistent(outfile=None, outform=None, overwrite=None):
    if not overwrite and (outform.upper() != "ASCII"):
        filename = get_abspath(outfile)
        if os.path.exists(filename):
            mesg = "Output file '%s' exists." % (outfile)
            raise Exception(mesg)


def convert_antenna_spec_autocorr(antenna):
    """Convert antenna (baseline) specification(s) to include autocorr data.

    Args:
        antenna (str): antenna specification

    Returns:
        str: tweaked antenna specification
    """
    if len(antenna) == 0:
        return antenna
    elif antenna.find(';') >= 0:
        # antenna selection is semi-colon separated list of baseline
        # specifications: 'SEL1;SEL2...'
        return ';'.join(map(convert_antenna_spec_autocorr, antenna.split(';')))
    elif antenna.find('&') < 0:
        # no '&' in the selection string
        #  -> 'ANT&&&'
        return antenna + '&&&'
    elif antenna.endswith('&&'):
        # 'ANT&&' or 'ANT&&&'
        #  -> as is
        return antenna
    elif antenna.endswith('&'):
        # 'ANT&'
        #  -> 'ANT&&&'
        return antenna.strip('&') + '&&&'
    else:
        # 'ANT1&ANT2' or 'ANT1&&ANT2'
        #  -> 'ANT1&&&;ANT2&&&'
        specs = [a for a in antenna.split('&') if len(a) > 0]
        return ';'.join(map(convert_antenna_spec_autocorr, specs))


def get_antenna_selection_include_autocorr(msname, antenna):
    """Get antenna selection string that includes autocorr data.

    Args:
        msname (str): name of MS
        antenna (str): antenna selection string

    Raises:
        RuntimeError: failed to handle antenna selection string

    Returns:
        str: antenna selection string including autocorr data
    """
    if len(antenna) == 0:
        # if no selection is specified, do nothing
        return antenna

    # test if given antenna selection is valid and if contains any autocorr data
    ms = mstool()
    sel = ms.msseltoindex(msname, baseline=antenna)
    if any([b[0] == b[1] for b in sel['baselines']]):
        antenna_autocorr = antenna
    else:
        antenna_autocorr = convert_antenna_spec_autocorr(antenna)
        casalog.post(
            'Tweaked antenna selection to include autocorr data: original "{}" tweaked "{}"'.format(
                antenna, antenna_autocorr
            )
        )
        # test if tweaked selection is valid
        sel = ms.msseltoindex(msname, baseline=antenna_autocorr)
        if all([b[0] != b[1] for b in sel['baselines']]):
            raise RuntimeError('Cannot handle antenna selection properly. Abort.')
    return antenna_autocorr


def get_nx_ny(n):
    nl = _to_list(n, int)
    if not nl:  # check for numpy int types
        nl = _to_list(n, numpy.integer)
    if len(nl) == 1:
        nx = ny = nl[0]
    else:
        nx = nl[0]
        ny = nl[1]
    return (nx, ny)


def get_cellx_celly(c,unit='arcsec'):
    if isinstance(c, str):
        cellx = celly = c
    elif type(c) in (list, tuple, numpy.ndarray):
        if len(c) == 1:
            cellx = celly = __to_quantity_string(c[0], unit)
        elif len(c) > 1:
            cellx = __to_quantity_string(c[0], unit)
            celly = __to_quantity_string(c[1], unit)
        else:
            cellx = celly = ''
    else:
        cellx = celly = __to_quantity_string(c, unit)
    return (cellx, celly)


def __to_quantity_string(v, unit='arcsec'):
    if isinstance(v, str):
        return v
    else:
        return '%s%s' % (v, unit)


def get_subtable_name(v):
    return v.replace('Table:', '').strip()


def get_spwids(selection, infile=None):
    # return a comma-separated string of spw IDs.
    # selection should be an output of ms.msseltoindex()

    spw_list = selection['spw']
    if len(spw_list) == 0:
        if infile is None:
            raise Exception("infile is needed when selection['spw'] is empty.")
        with tbmanager(os.path.join(infile, 'DATA_DESCRIPTION')) as tb:
            spw_list = tb.getcol('SPECTRAL_WINDOW_ID')

    l = []
    for item in spw_list:
        l.append(str(item))
    return ','.join(l)


def parse_wavenumber_param(wn):
    if isinstance(wn, list):
        _check_positive_or_zero(wn)
        wn.sort()
        return ','.join(_get_strlist(wn))
    elif isinstance(wn, tuple):
        _check_positive_or_zero(wn)
        wn_list = list(wn)
        wn_list.sort()
        return ','.join(_get_strlist(wn_list))
    elif isinstance(wn, int):
        _check_positive_or_zero(wn)
        return str(wn)
    elif isinstance(wn, str):
        if ',' in wn:                            # cases 'a,b,c,...'
            val0 = wn.split(',')
            _check_positive_or_zero(val0)
            val = []
            for v in val0: val.append(int(v))
            val.sort()
            res = list(set(val))  # uniq
        elif '-' in wn:                          # case 'a-b' : return [a,a+1,...,b-1,b]
            val = wn.split('-')
            _check_positive_or_zero(val)
            val = [int(val[0]), int(val[1])]
            val.sort()
            res = [i for i in range(val[0], val[1]+1)]
        elif '~' in wn:                          # case 'a~b' : return [a,a+1,...,b-1,b]
            val = wn.split('~')
            _check_positive_or_zero(val)
            val = [int(val[0]), int(val[1])]
            val.sort()
            res = [i for i in range(val[0], val[1]+1)]
        elif wn[:2] == '<=' or wn[:2] == '=<':   # cases '<=a','=<a' : return [0,1,...,a-1,a]
            val = wn[2:]
            _check_positive_or_zero(val)
            res = [i for i in range(int(val)+1)]
        elif wn[-2:] == '>=' or wn[-2:] == '=>': # cases 'a>=','a=>' : return [0,1,...,a-1,a]
            val = wn[:-2]
            _check_positive_or_zero(val)
            res = [i for i in range(int(val)+1)]
        elif wn[0] == '<':                       # case '<a' :         return [0,1,...,a-2,a-1]
            val = wn[1:]
            _check_positive_or_zero(val, False)
            res = [i for i in range(int(val))]
        elif wn[-1] == '>':                      # case 'a>' :         return [0,1,...,a-2,a-1]
            val = wn[:-1]
            _check_positive_or_zero(val, False)
            res = [i for i in range(int(val))]
        elif wn[:2] == '>=' or wn[:2] == '=>':   # cases '>=a','=>a' : return [a,-999], which is
                                                 #                     then interpreted in C++
                                                 #                     side as [a,a+1,...,a_nyq]
                                                 #                     (CAS-3759)
            val = wn[2:]
            _check_positive_or_zero(val)
            res = [int(val), -999]
        elif wn[-2:] == '<=' or wn[-2:] == '=<': # cases 'a<=','a=<' : return [a,-999], which is
                                                 #                     then interpreted in C++
                                                 #                     side as [a,a+1,...,a_nyq]
                                                 #                     (CAS-3759)
            val = wn[:-2]
            _check_positive_or_zero(val)
            res = [int(val), -999]
        elif wn[0] == '>':                       # case '>a' :         return [a+1,-999], which is
                                                 #                     then interpreted in C++
                                                 #                     side as [a+1,a+2,...,a_nyq]
                                                 #                     (CAS-3759)
            val0 = wn[1:]
            val = int(val0)+1
            _check_positive_or_zero(val)
            res = [val, -999]
        elif wn[-1] == '<':                      # case 'a<' :         return [a+1,-999], which is
                                                 #                     then interpreted in C++
                                                 #                     side as [a+1,a+2,...,a_nyq]
                                                 #                     (CAS-3759)
            val0 = wn[:-1]
            val = int(val0)+1
            _check_positive_or_zero(val)
            res = [val, -999]
        else:
            _check_positive_or_zero(wn)
            res = [int(wn)]

        # return res
        return ','.join(_get_strlist(res))
    else:
        msg = 'wrong value given for addwn/rejwn'
        raise RuntimeError(msg)


def _check_positive_or_zero(param, allowzero=True):
    msg = 'wrong value given for addwn/rejwn'
    try:
        if isinstance(param, list) or isinstance(param, tuple):
            for i in range(len(param)):
                _do_check_positive_or_zero(int(param[i]), allowzero)
        elif isinstance(param, int):
            _do_check_positive_or_zero(param, allowzero)
        elif isinstance(param, str):
            _do_check_positive_or_zero(int(param), allowzero)
        else:
            raise RuntimeError(msg)
    except:
        raise RuntimeError(msg)


def _get_strlist(param):
    res = []
    for i in range(len(param)):
        res.append(str(param[i]))
    return res


def _do_check_positive_or_zero(param, allowzero):
    msg = 'wrong value given for addwn/rejwn'
    if (param < 0) or ((param == 0) and not allowzero):
        raise RuntimeError(msg)


def _is_sequence_or_number(param, ptype=int):
    """
    Returns true if input is an array type or a number with a give data type.
    Arguments
        param : an array or number to test
        ptype : the data type that param should be.
    """
    if hasattr(param, '__iter__'):
        out = True
        for p in param:
            out &= isinstance(p, ptype)
        return out
    else:
        return isinstance(param, ptype)


def _to_list(param, ptype=int, convert=False):
    """
    Convert a number, an array type or a string to a list.
    The function returns None if input values are not ptype and convert=False.
    When convert is True, force converting input values to a list of ptype.
    """
    if isinstance(param, ptype):  # a string or a number
        if ptype is str:
            return param.split()
        elif convert:
            return [ptype(param)]
        else:
            return [param]
    if _is_sequence_or_number(param, ptype):
        return list(param)
    elif convert:
        return [ptype(p) for p in param]
    return None
