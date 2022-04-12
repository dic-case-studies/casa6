from collections import Counter
import datetime
import os
import shutil

from casatasks import casalog
from casatools import ms as mstool
from casatools import singledishms

from . import sdutil
from .mstools import write_history

ms = mstool()


@sdutil.callable_sdtask_decorator
def sdbaseline(infile=None, datacolumn=None, antenna=None, field=None,
               spw=None, timerange=None, scan=None, pol=None, intent=None,
               reindex=None, maskmode=None, thresh=None, avg_limit=None,
               minwidth=None, edge=None, blmode=None, dosubtract=None,
               blformat=None, bloutput=None, bltable=None, blfunc=None,
               order=None, npiece=None, applyfft=None, fftmethod=None,
               fftthresh=None, addwn=None, rejwn=None, clipthresh=None,
               clipniter=None, blparam=None, verbose=None,
               updateweight=None, sigmavalue=None,
               showprogress=None, minnrow=None,
               outfile=None, overwrite=None):

    temp_outfile = ''

    try:
        # CAS-12985 requests the following params be given case insensitively,
        # so they need to be converted to lowercase here (2021/1/28 WK)
        blfunc = blfunc.lower()
        blmode = blmode.lower()
        fftmethod = fftmethod.lower()
        if isinstance(fftthresh, str):
            fftthresh = fftthresh.lower()

        if not os.path.exists(infile):
            raise ValueError("infile='" + str(infile) + "' does not exist.")
        if (outfile == '') or not isinstance(outfile, str):
            outfile = infile.rstrip('/') + '_bs'
            casalog.post("outfile is empty or non-string. set to '" + outfile + "'")
        if (not overwrite) and dosubtract and os.path.exists(outfile):
            raise ValueError("outfile='%s' exists, and cannot overwrite it." % (outfile))

        if (maskmode == 'interact'):
            raise ValueError("maskmode='%s' is not supported yet" % maskmode)
        if (blfunc == 'variable') and not os.path.exists(blparam):
            raise ValueError("input file '%s' does not exists" % blparam)

        if (spw == ''):
            spw = '*'

        if (blmode == 'fit'):
            temp_outfile = _do_fit(infile, datacolumn, antenna, field, spw, timerange, scan,
                                   pol, intent, reindex, maskmode, thresh, avg_limit, minwidth,
                                   edge, dosubtract, blformat, bloutput, blfunc, order, npiece,
                                   applyfft, fftmethod, fftthresh, addwn, rejwn, clipthresh,
                                   clipniter, blparam, verbose, updateweight, sigmavalue,
                                   outfile, overwrite)
        elif (blmode == 'apply'):
            _do_apply(infile, datacolumn, antenna, field, spw, timerange, scan, pol, intent,
                      reindex, bltable, updateweight, sigmavalue, outfile, overwrite)

        # Remove {WEIGHT|SIGMA}_SPECTRUM columns if updateweight=True (CAS-13161)
        if updateweight:
            with sdutil.table_manager(outfile, nomodify=False) as mytb:
                cols_spectrum = ['WEIGHT_SPECTRUM', 'SIGMA_SPECTRUM']
                cols_remove = [col for col in cols_spectrum if col in mytb.colnames()]
                if len(cols_remove) > 0:
                    mytb.removecols(' '.join(cols_remove))

        # Write history to outfile
        if dosubtract:
            param_names = sdbaseline.__code__.co_varnames[:sdbaseline.__code__.co_argcount]
            var_local = locals()
            param_vals = [var_local[p] for p in param_names]
            write_history(ms, outfile, 'sdbaseline', param_names,
                          param_vals, casalog)

    finally:
        if (not dosubtract):
            # Remove (skeleton) outfile
            if temp_outfile != '':
                outfile = temp_outfile
            remove_data(outfile)


blformat_item = ['csv', 'text', 'table']
blformat_ext = ['csv', 'txt', 'bltable']

mesg_invalid_wavenumber = 'wrong value given for addwn/rejwn'


def remove_data(filename):
    if not os.path.exists(filename):
        return

    if os.path.isdir(filename):
        shutil.rmtree(filename)
    elif os.path.isfile(filename):
        os.remove(filename)
    else:
        # could be a symlink
        os.remove(filename)


def is_empty(blformat):
    """Check if blformat is empty.

    returns True if blformat is None, '', [] and
    a string list containing only '' (i.e., ['', '', ..., ''])
    """
    if isinstance(blformat, list):
        return all(map(is_empty, blformat))

    return not blformat


def prepare_for_blformat_bloutput(infile, blformat, bloutput, overwrite):
    # force to string list
    blformat = force_to_string_list(blformat, 'blformat')
    bloutput = force_to_string_list(bloutput, 'bloutput')

    # the default bloutput value '' is expanded to a list
    # with length of blformat, and with '' throughout.
    if (bloutput == ['']):
        bloutput *= len(blformat)

    # check length
    if (len(blformat) != len(bloutput)):
        raise ValueError('blformat and bloutput must have the same length.')

    # check duplication
    if has_duplicate_nonnull_element(blformat):
        raise ValueError('duplicate elements in blformat.')
    if has_duplicate_nonnull_element_ex(bloutput, blformat):
        raise ValueError('duplicate elements in bloutput.')

    # fill bloutput items to be output, then rearrange them
    # in the order of blformat_item.
    bloutput = normalise_bloutput(infile, blformat, bloutput, overwrite)

    return blformat, bloutput


def force_to_string_list(s, name):
    mesg = '%s must be string or list of string.' % name
    if isinstance(s, str):
        s = [s]
    elif isinstance(s, list):
        for i in range(len(s)):
            if not isinstance(s[i], str):
                raise ValueError(mesg)
    else:
        raise ValueError(mesg)
    return s


def has_duplicate_nonnull_element(in_list):
    # return True if in_list has duplicated elements other than ''
    duplicates = [key for key, val in Counter(in_list).items() if val > 1]
    len_duplicates = len(duplicates)

    if (len_duplicates >= 2):
        return True
    elif (len_duplicates == 1):
        return (duplicates[0] != '')
    else:  # len_duplicates == 0
        return False


def has_duplicate_nonnull_element_ex(lst, base):
    # lst and base must have the same length.
    #
    # (1) extract elements from lst and make a new list
    #     if the element of base with the same index
    #     is not ''.
    # (2) check if the list made in (1) has duplicated
    #     elements other than ''.

    return has_duplicate_nonnull_element(
        [lst[i] for i in range(len(lst)) if base[i] != ''])


def normalise_bloutput(infile, blformat, bloutput, overwrite):
    return [get_normalised_name(infile, blformat, bloutput, item[0], item[1], overwrite)
            for item in zip(blformat_item, blformat_ext)]


def get_normalised_name(infile, blformat, bloutput, name, ext, overwrite):
    fname = ''
    blformat_lower = [s.lower() for s in blformat]
    if (name in blformat_lower):
        fname = bloutput[blformat_lower.index(name)]
        if (fname == ''):
            fname = infile + '_blparam.' + ext
    if os.path.exists(fname):
        if overwrite:
            remove_data(fname)
        else:
            raise ValueError(fname + ' exists.')
    return fname


def output_bloutput_text_header(blformat, bloutput, blfunc, maskmode, infile, outfile):
    fname = bloutput[blformat_item.index('text')]
    if (fname == ''):
        return

    with open(fname, 'w') as f:
        info = [['Source Table', infile],
                ['Output File', outfile if (outfile != '') else infile],
                ['Mask mode', maskmode]]

        separator = '#' * 60 + '\n'

        f.write(separator)
        for i in range(len(info)):
            f.write('%12s: %s\n' % tuple(info[i]))
        f.write(separator)
        f.write('\n')


def get_temporary_file_name(basename):
    name = basename + '_sdbaseline_pid' + str(os.getpid()) + '_' \
        + datetime.datetime.now().strftime('%Y%m%d%H%M%S%f')
    return name


def parse_wavenumber_param(wn):
    if isinstance(wn, bool):
        raise ValueError(mesg_invalid_wavenumber)
    elif isinstance(wn, list):
        __check_positive_or_zero(wn)
        wn_uniq = list(set(wn))
        wn_uniq.sort()
        return ','.join(__get_strlist(wn_uniq))
    elif isinstance(wn, tuple):
        __check_positive_or_zero(wn)
        wn_uniq = list(set(wn))
        wn_uniq.sort()
        return ','.join(__get_strlist(wn_uniq))
    elif isinstance(wn, int):
        __check_positive_or_zero(wn)
        return str(wn)
    elif isinstance(wn, str):
        if '.' in wn:
            # case of float value as string
            raise ValueError(mesg_invalid_wavenumber)
        elif ',' in wn:
            # cases 'a,b,c,...'
            val0 = wn.split(',')
            __check_positive_or_zero(val0)
            val = []
            for v in val0:
                val.append(int(v))
            val.sort()
            res = list(set(val))  # uniq
        elif '-' in wn:
            # case 'a-b' : return [a,a+1,...,b-1,b]
            val = wn.split('-')
            __check_positive_or_zero(val)
            val = [int(val[0]), int(val[1])]
            val.sort()
            res = [i for i in range(val[0], val[1] + 1)]
        elif '~' in wn:
            # case 'a~b' : return [a,a+1,...,b-1,b]
            val = wn.split('~')
            __check_positive_or_zero(val)
            val = [int(val[0]), int(val[1])]
            val.sort()
            res = [i for i in range(val[0], val[1] + 1)]
        elif wn[:2] == '<=' or wn[:2] == '=<':
            # cases '<=a','=<a' : return [0,1,...,a-1,a]
            val = wn[2:]
            __check_positive_or_zero(val)
            res = [i for i in range(int(val) + 1)]
        elif wn[-2:] == '>=' or wn[-2:] == '=>':
            # cases 'a>=','a=>' : return [0,1,...,a-1,a]
            val = wn[:-2]
            __check_positive_or_zero(val)
            res = [i for i in range(int(val) + 1)]
        elif wn[0] == '<':
            # case '<a' :         return [0,1,...,a-2,a-1]
            val = wn[1:]
            __check_positive_or_zero(val, False)
            res = [i for i in range(int(val))]
        elif wn[-1] == '>':
            # case 'a>' :         return [0,1,...,a-2,a-1]
            val = wn[:-1]
            __check_positive_or_zero(val, False)
            res = [i for i in range(int(val))]
        elif wn[:2] == '>=' or wn[:2] == '=>':
            # cases '>=a','=>a' : return [a,-999], which is
            #                     then interpreted in C++
            #                     side as [a,a+1,...,a_nyq]
            #                     (CAS-3759)
            val = wn[2:]
            __check_positive_or_zero(val)
            res = [int(val), -999]
        elif wn[-2:] == '<=' or wn[-2:] == '=<':
            # cases 'a<=','a=<' : return [a,-999], which is
            #                     then interpreted in C++
            #                     side as [a,a+1,...,a_nyq]
            #                     (CAS-3759)
            val = wn[:-2]
            __check_positive_or_zero(val)
            res = [int(val), -999]
        elif wn[0] == '>':
            # case '>a' :         return [a+1,-999], which is
            #                     then interpreted in C++
            #                     side as [a+1,a+2,...,a_nyq]
            #                     (CAS-3759)
            val0 = wn[1:]
            val = int(val0) + 1
            __check_positive_or_zero(val)
            res = [val, -999]
        elif wn[-1] == '<':
            # case 'a<' :         return [a+1,-999], which is
            #                     then interpreted in C++
            #                     side as [a+1,a+2,...,a_nyq]
            #                     (CAS-3759)
            val0 = wn[:-1]
            val = int(val0) + 1
            __check_positive_or_zero(val)
            res = [val, -999]
        else:
            # case 'a'
            __check_positive_or_zero(wn)
            res = [int(wn)]

        # return res
        return ','.join(__get_strlist(res))
    else:
        raise ValueError(mesg_invalid_wavenumber)


def __get_strlist(param):
    return [str(p) for p in param]


def check_fftthresh(fftthresh):
    """Validate fftthresh value.

    The fftthresh must be one of the following:
    (1) positive value (float, integer or string)
    (2) 'top' + positive integer value
    (3) positive float value + 'sigma'
    """
    has_invalid_type = False
    val_not_positive = False

    if isinstance(fftthresh, bool):
        # Checking for bool must precede checking for integer
        has_invalid_type = True
    elif isinstance(fftthresh, int) or isinstance(fftthresh, float):
        if (fftthresh <= 0.0):
            val_not_positive = True
    elif isinstance(fftthresh, str):
        try:
            if (3 < len(fftthresh)) and (fftthresh[:3] == 'top'):
                if (int(fftthresh[3:]) <= 0):
                    val_not_positive = True
            elif (5 < len(fftthresh)) and (fftthresh[-5:] == 'sigma'):
                if (float(fftthresh[:-5]) <= 0.0):
                    val_not_positive = True
            else:
                if (float(fftthresh) <= 0.0):
                    val_not_positive = True
        except Exception:
            raise ValueError('fftthresh has a wrong format.')
    else:
        has_invalid_type = True

    if has_invalid_type:
        raise ValueError('fftthresh must be float or integer or string.')
    if val_not_positive:
        raise ValueError('threshold given to fftthresh must be positive.')


def __check_positive_or_zero(param, allowzero=True):
    if isinstance(param, list) or isinstance(param, tuple):
        for i in range(len(param)):
            __do_check_positive_or_zero(int(param[i]), allowzero)
    elif isinstance(param, int):
        __do_check_positive_or_zero(param, allowzero)
    elif isinstance(param, str):
        __do_check_positive_or_zero(int(param), allowzero)
    else:
        raise ValueError(mesg_invalid_wavenumber)


def __do_check_positive_or_zero(param, allowzero):
    if (param < 0) or ((param == 0) and not allowzero):
        raise ValueError(mesg_invalid_wavenumber)


def prepare_for_baselining(**keywords):
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
    elif blfunc == 'sinusoid':
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
    for key in keys:
        params[key] = keywords[key]

    baseline_func = getattr(keywords['sdms'], funcname)

    return params, baseline_func


def remove_sorted_table_keyword(infile):
    res = {'is_sorttab': False, 'sorttab_keywd': '', 'sorttab_name': ''}

    with sdutil.table_manager(infile, nomodify=False) as tb:
        sorttab_keywd = 'SORTED_TABLE'
        if sorttab_keywd in tb.keywordnames():
            res['is_sorttab'] = True
            res['sorttab_keywd'] = sorttab_keywd
            res['sorttab_name'] = tb.getkeyword(sorttab_keywd)
            tb.removekeyword(sorttab_keywd)

    return res


def restore_sorted_table_keyword(infile, sorttab_info):
    if sorttab_info['is_sorttab'] and (sorttab_info['sorttab_name'] != ''):
        with sdutil.table_manager(infile, nomodify=False) as tb:
            tb.putkeyword(sorttab_info['sorttab_keywd'], sorttab_info['sorttab_name'])


def _do_apply(infile, datacolumn, antenna, field, spw, timerange, scan, pol, intent,
              reindex, bltable, updateweight, sigmavalue, outfile, overwrite):
    if not os.path.exists(bltable):
        raise ValueError("file specified in bltable '%s' does not exist." % bltable)

    # Note: the condition "infile != outfile" in the following line is for safety
    # to prevent from accidentally removing infile by setting outfile=infile.
    # Don't remove it.
    if overwrite and (infile != outfile) and os.path.exists(outfile):
        remove_data(outfile)

    sorttab_info = remove_sorted_table_keyword(infile)

    with sdutil.tool_manager(infile, singledishms) as mysdms:
        selection = ms.msseltoindex(vis=infile, spw=spw, field=field,
                                    baseline=antenna, time=timerange,
                                    scan=scan)
        mysdms.set_selection(spw=sdutil.get_spwids(selection), field=field,
                             antenna=antenna, timerange=timerange,
                             scan=scan, polarization=pol, intent=intent,
                             reindex=reindex)
        mysdms.apply_baseline_table(bltable=bltable,
                                    datacolumn=datacolumn,
                                    spw=spw,
                                    updateweight=updateweight,
                                    sigmavalue=sigmavalue,
                                    outfile=outfile)

    restore_sorted_table_keyword(infile, sorttab_info)


def _do_fit(infile, datacolumn, antenna, field, spw, timerange, scan, pol, intent,
            reindex, maskmode, thresh, avg_limit, minwidth, edge, dosubtract, blformat,
            bloutput, blfunc, order, npiece, applyfft, fftmethod, fftthresh, addwn,
            rejwn, clipthresh, clipniter, blparam, verbose, updateweight, sigmavalue,
            outfile, overwrite):

    temp_outfile = ''

    if (not dosubtract) and is_empty(blformat):
        raise ValueError("blformat must be specified when dosubtract is False")

    blformat, bloutput = prepare_for_blformat_bloutput(infile, blformat, bloutput, overwrite)

    output_bloutput_text_header(blformat, bloutput, blfunc, maskmode, infile, outfile)

    # Set temporary name for output MS if dosubtract is False and outfile exists
    # for not removing/overwriting outfile that already exists
    if os.path.exists(outfile):
        # Note: the condition "infile != outfile" in the following line is for safety
        # to prevent from accidentally removing infile by setting outfile=infile
        # Don't remove it.
        if dosubtract and overwrite and (infile != outfile):
            remove_data(outfile)
        elif (not dosubtract):
            outfile = get_temporary_file_name(infile)
            temp_outfile = outfile

    if (blfunc == 'variable'):
        sorttab_info = remove_sorted_table_keyword(infile)
    elif (blfunc == 'sinusoid'):
        addwn = parse_wavenumber_param(addwn)
        rejwn = parse_wavenumber_param(rejwn)
        check_fftthresh(fftthresh)

    with sdutil.tool_manager(infile, singledishms) as mysdms:
        selection = ms.msseltoindex(vis=infile, spw=spw, field=field, baseline=antenna,
                                    time=timerange, scan=scan)
        mysdms.set_selection(spw=sdutil.get_spwids(selection), field=field, antenna=antenna,
                             timerange=timerange, scan=scan, polarization=pol, intent=intent,
                             reindex=reindex)
        params, func = prepare_for_baselining(sdms=mysdms,
                                              blfunc=blfunc,
                                              datacolumn=datacolumn,
                                              outfile=outfile,
                                              bloutput=','.join(bloutput),
                                              dosubtract=dosubtract,
                                              spw=spw,
                                              pol=pol,
                                              linefinding=(maskmode == 'auto'),
                                              threshold=thresh,
                                              avg_limit=avg_limit,
                                              minwidth=minwidth,
                                              edge=edge,
                                              order=order,
                                              npiece=npiece,
                                              applyfft=applyfft,
                                              fftmethod=fftmethod,
                                              fftthresh=fftthresh,
                                              addwn=addwn,
                                              rejwn=rejwn,
                                              clip_threshold_sigma=clipthresh,
                                              num_fitting_max=clipniter + 1,
                                              blparam=blparam,
                                              verbose=verbose,
                                              updateweight=updateweight,
                                              sigmavalue=sigmavalue)
        func(**params)

    if (blfunc == 'variable'):
        restore_sorted_table_keyword(infile, sorttab_info)

    return temp_outfile
