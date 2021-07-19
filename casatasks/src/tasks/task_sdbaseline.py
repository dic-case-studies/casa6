import os
import shutil
from collections import Counter

# get is_CASA6 and is_python3
from casatasks.private.casa_transition import *

if is_CASA6:
    from casatasks import casalog
    from casatools import ms as mstool
    from casatools import msmetadata, singledishms, table

    from . import sdutil
    from .mstools import write_history

    ms = mstool()
    sdms = singledishms()
    tb = table()
    msmd = msmetadata()
else:
    import sdutil
    from mstools import write_history
    from taskinit import casalog, gentools
    ms, sdms, tb, msmd = gentools(['ms', 'sdms', 'tb', 'msmd'])

@sdutil.sdtask_decorator
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

    try:
        # CAS-12985 requests the following params be given case insensitively,
        # so they need to be converted to lowercase here (2021/1/28 WK)
        blfunc = blfunc.lower()
        blmode = blmode.lower()
        fftmethod = fftmethod.lower()
        if isinstance(fftthresh, str):
            fftthresh = fftthresh.lower()

        if not os.path.exists(infile):
            raise Exception("infile='" + str(infile) + "' does not exist.")
        if (outfile == '') or not isinstance(outfile, str):
            #casalog.post("type=%s, value=%s" % (type(outfile), str(outfile)))
            #raise ValueError, "outfile name is empty."
            outfile = infile.rstrip('/') + '_bs'
            casalog.post("outfile is empty or non-string. set to '" + outfile + "'")
        if os.path.exists(outfile) and not overwrite:
            raise Exception("outfile='%s' exists, and cannot overwrite it." % (outfile))
        if (maskmode == 'interact'):
            raise ValueError("maskmode='%s' is not supported yet" % maskmode)
        if (blfunc == 'variable' and not os.path.exists(blparam)):
            raise ValueError("input file '%s' does not exists" % blparam)
        blparam_file = infile + '_blparam.txt'
        if os.path.exists(blparam_file):
            remove_data(blparam_file)  # CAS-11781
        
        if (spw == ''): spw = '*'

        if (blmode == 'apply'):
            if not os.path.exists(bltable):
                raise ValueError("file specified in bltable '%s' does not exist." % bltable)

            sorttab_info = remove_sorted_table_keyword(infile)

            if overwrite and os.path.exists(outfile) and (infile != outfile):
                os.system('rm -rf %s' % outfile)

            selection = ms.msseltoindex(vis=infile, spw=spw, field=field, 
                                        baseline=antenna, time=timerange, 
                                        scan=scan)
            sdms.open(infile)
            sdms.set_selection(spw=sdutil.get_spwids(selection), field=field, 
                               antenna=antenna, timerange=timerange, 
                               scan=scan, polarization=pol, intent=intent,
                               reindex=reindex)
            sdms.apply_baseline_table(bltable=bltable,
                                      datacolumn=datacolumn,
                                      spw=spw,
                                      updateweight=updateweight,
                                      sigmavalue=sigmavalue,
                                      outfile=outfile)
            sdms.close()
            
            restore_sorted_table_keyword(infile, sorttab_info)
            
        elif (blmode == 'fit'):

            if(blfunc == 'sinusoid'):
                addwn = sdutil.parse_wavenumber_param(addwn)
                rejwn = sdutil.parse_wavenumber_param(rejwn)
                check_fftthresh(fftthresh)

            blformat, bloutput = prepare_for_blformat_bloutput(infile, blformat, bloutput, overwrite)

            output_bloutput_text_header(blformat, bloutput,
                                        blfunc, maskmode,
                                        infile, outfile)
            
            if (blfunc == 'variable'):
                sorttab_info = remove_sorted_table_keyword(infile)
        
            if overwrite and os.path.exists(outfile) and (infile != outfile):
                os.system('rm -rf %s' % outfile)

            selection = ms.msseltoindex(vis=infile, spw=spw, field=field, 
                                        baseline=antenna, time=timerange, 
                                        scan=scan)
            sdms.open(infile)
            sdms.set_selection(spw=sdutil.get_spwids(selection),
                               field=field, antenna=antenna,
                               timerange=timerange, scan=scan,
                               polarization=pol, intent=intent,
                               reindex=reindex)
            params, func = prepare_for_baselining(blfunc=blfunc,
                                                  datacolumn=datacolumn,
                                                  outfile=outfile,
                                                  bloutput=','.join(bloutput),
                                                  dosubtract=dosubtract,
                                                  spw=spw,
                                                  pol=pol,
                                                  linefinding=(maskmode=='auto'),
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
                                                  num_fitting_max=clipniter+1,
                                                  blparam=blparam,
                                                  verbose=verbose,
                                                  updateweight=updateweight,
                                                  sigmavalue=sigmavalue)
            func(**params)
            sdms.close()
            
            if (blfunc == 'variable'):
                restore_sorted_table_keyword(infile, sorttab_info)

        # Remove {WEIGHT|SIGMA}_SPECTRUM columns if updateweight=True (CAS-13161)
        if updateweight:
            with sdutil.table_manager(outfile, nomodify=False) as mytb:
                cols_remove = []
                for col in ['WEIGHT_SPECTRUM', 'SIGMA_SPECTRUM']:
                    if col in mytb.colnames():
                        cols_remove.append(col)
                if len(cols_remove) > 0:
                    mytb.removecols(' '.join(cols_remove))

        # Write history to outfile
        param_names = sdbaseline.__code__.co_varnames[:sdbaseline.__code__.co_argcount]
        if is_python3:
            vars = locals()
            param_vals = [vars[p] for p in param_names]
        else:
            param_vals = [eval(p) for p in param_names]
        write_history(ms, outfile, 'sdbaseline', param_names,
                      param_vals, casalog)


    except Exception:
        raise


blformat_item = ['csv', 'text', 'table']
blformat_ext  = ['csv', 'txt',  'bltable']


def remove_data(filename):
    if os.path.exists(filename):
        if os.path.isdir(filename):
            shutil.rmtree(filename)
        elif os.path.isfile(filename):
            os.remove(filename)
        else:
            # could be a symlink
            os.remove(filename)

def check_fftthresh(fftthresh):
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

def prepare_for_blformat_bloutput(infile, blformat, bloutput, overwrite):
    # force to string list
    blformat = force_to_string_list(blformat, 'blformat')
    bloutput = force_to_string_list(bloutput, 'bloutput')

    # the default bloutput value '' is expanded to a list 
    # with length of blformat, and with '' throughout.
    if (bloutput == ['']): bloutput *= len(blformat)

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
    if isinstance(s, str): s = [s]
    elif isinstance(s, list):
        for i in range(len(s)):
            if not isinstance(s[i], str):
                raise ValueError(mesg)
    else:
        raise ValueError(mesg)
    return s

def has_duplicate_nonnull_element(in_list):
    #return True if in_list has duplicated elements other than ''
    duplicates = [key for key, val in Counter(in_list).items() if val > 1]
    len_duplicates = len(duplicates)
    
    if (len_duplicates >= 2):
        return True
    elif (len_duplicates == 1):
        return (duplicates[0] != '')
    else: #len_duplicates == 0
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
    normalised_bloutput = []
    for item in zip(blformat_item, blformat_ext):
        normalised_bloutput.append(
            get_normalised_name(infile, blformat, bloutput, item[0], item[1], overwrite))
    return normalised_bloutput

def get_normalised_name(infile, blformat, bloutput, name, ext, overwrite):
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

def output_bloutput_text_header(blformat, bloutput, blfunc, maskmode, infile, outfile):
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
    
    
def remove_sorted_table_keyword(infile):
    res = {'is_sorttab': False, 'sorttab_keywd': '', 'sorttab_name': ''}
    with sdutil.table_manager(infile, nomodify=False) as tb:
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

def restore_sorted_table_keyword(infile, sorttab_info):
    if sorttab_info['is_sorttab'] and (sorttab_info['sorttab_name'] != ''):
        with sdutil.table_manager(infile, nomodify=False) as tb:
            try:
                tb.putkeyword(sorttab_info['sorttab_keywd'],
                              sorttab_info['sorttab_name'])
            except Exception:
                raise
