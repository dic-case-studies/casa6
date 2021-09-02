import os

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatasks import casalog
    from . import sdutil
    from .sdutil import calibrater_manager
else:
    from taskinit import casalog
    from . import sdutil
    from sdutil import cbmanager as calibrater_manager

DEFAULT_VALUE = {'interp': 'linear',
                 'spwmap': [-1]}

def parse_interp_item(interp):
    assert isinstance(interp, str)
    if len(interp) == 0:
        return DEFAULT_VALUE['interp']
    else:
        return interp

def parse_interp(interp, index):
    assert index >= 0
    if isinstance(interp, str):
        # interp is a string that is valid to all applytables
        return parse_interp_item(interp)
    elif hasattr(interp, '__iter__'):
        # interp is a list of strings
        if index >= len(interp):
            # wrong index or empty list
            return DEFAULT_VALUE['interp']
        else:
            # interp is a list of strings
            return parse_interp_item(interp[index])
    assert False

def parse_spwmap_item(spwmap):
    assert hasattr(spwmap, '__iter__')
    if len(spwmap) == 0:
        return DEFAULT_VALUE['spwmap']
    else:
        return spwmap    

def parse_spwmap(spwmap, index):
    assert hasattr(spwmap, '__iter__')
    assert index >= 0
    if len(spwmap) == 0:
        # empty list
        return DEFAULT_VALUE['spwmap']
    elif all(map(lambda x: hasattr(x, '__iter__'), spwmap)):
        # spwmap is list-of-list
        if index >= len(spwmap):
            # maybe wrong index
            return DEFAULT_VALUE['spwmap']
        else:
            return parse_spwmap_item(spwmap[index])
    elif all(map(lambda x: isinstance(x, int), spwmap)):
        # maybe single spwmap that is valid to all applytables
        return spwmap
    assert False

@sdutil.sdtask_decorator
def sdgaincal(infile=None, calmode=None, radius=None, smooth=None, 
              antenna=None, field=None, spw=None, scan=None, intent=None, 
              applytable=None, interp=None, spwmap=None, outfile='', overwrite=False):

    # outfile must be specified
    if (outfile == '') or not isinstance(outfile, str):
        raise ValueError("outfile is empty.")

    # overwrite check
    if os.path.exists(outfile) and not overwrite:
        raise RuntimeError(outfile + ' exists.')

    if infile is None or not isinstance(infile, str) or not os.path.exists(infile):
        raise RuntimeError('infile not found - please verify the name')

    # Calibrater tool
    with calibrater_manager(infile) as mycb:

        # select data
        if isinstance(antenna, str) and len(antenna) > 0:
            baseline = '{ant}&&&'.format(ant=antenna)
        else:
            baseline = ''
        mycb.selectvis(spw=spw, scan=scan, field=field, intent=intent, baseline=baseline)
        
        # set apply
        casalog.post('interp="{0}" spwmap={1}'.format(interp, spwmap))
        if isinstance(applytable, str):
            if len(applytable) > 0:
                thisinterp = parse_interp(interp, 0)
                thisspwmap = parse_spwmap(spwmap, 0)
                casalog.post('thisinterp="{0}" thisspwmap={1}'.format(thisinterp, thisspwmap))
                mycb.setapply(table=applytable, interp=thisinterp, spwmap=thisspwmap)
        elif hasattr(applytable, '__iter__'):
            # list type 
            for i in range(len(applytable)):
                table = applytable[i]
                if isinstance(table, str) and len(table) > 0:
                    thisinterp = parse_interp(interp, i)
                    thisspwmap = parse_spwmap(spwmap, i)
                    casalog.post('thisinterp="{0}" thisspwmap={1}'.format(thisinterp, thisspwmap))
                    mycb.setapply(table=table, interp=thisinterp, spwmap=thisspwmap)
                else:
                    raise RuntimeError('wrong type of applytable item ({0}). it should be string'.format(type(table)))
        else:
            raise RuntimeError('wrong type of applytable ({0}). it should be string or list'.format(type(applytable)))
        
        # set solve
        if calmode == 'doublecircle':
            if radius is None:
                raise RuntimeError('radius must be specified.')
            elif not isinstance(radius, str):
                rcenter = '%sarcsec'%(radius)
            else:
                try:
                    # if radius is a string only consists of numeric value without unit, 
                    # it will succeed.
                    rcenter = '%sarcsec'%(float(radius))
                except:
                    # if the above fails, it may indicate that the string contains unit
                    rcenter = radius
            mycb.setsolve(type='SDGAIN_OTFD', table=outfile, radius=rcenter, smooth=smooth)
        else:
            raise RuntimeError('Unknown calibration mode: \'{mode}\''.format(mode=calmode))

        # solve
        mycb.solve()


