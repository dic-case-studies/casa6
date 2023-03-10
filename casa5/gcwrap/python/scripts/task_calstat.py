from __future__ import absolute_import

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatasks import casalog
    from casatools import table
    _tb = table( )
else:
    from taskinit import *
    _tb = tb

def calstat(caltable=None,
            axis=None,
            datacolumn=None):

    casalog.origin('calstat')

    _tb.close()
    _tb.open(caltable)

    if axis in ['amp', 'amplitude', 'phase', 'imag', 'imaginary', 'real']:
        complex_type = axis
        col = datacolumn
    else:
        complex_type = ''
        col = axis

    if col == 'corrected':
        col = 'CORRECTED_DATA'
    if col == 'model':
        col = 'MODEL_DATA'
    s = _tb.statistics(column=col.upper(),
                      complex_value=complex_type)
        
    _tb.close()

    for stats in s.keys():
        casalog.post(stats + " values --- ", "NORMAL")
        
        if s[stats]['npts'] > 0:
            casalog.post("         -- number of points [npts]:           " + str(int(round(s[stats]['npts']))), "NORMAL")
            casalog.post("         -- minimum value [min]:               " + str(s[stats]['min'  ]), "NORMAL")
            casalog.post("         -- maximum value [max]:               " + str(s[stats]['max'  ]), "NORMAL")
            casalog.post("         -- Sum of values [sum]:               " + str(s[stats]['sum'  ]), "NORMAL")
            casalog.post("         -- Sum of squared values [sumsq]:     " + str(s[stats]['sumsq']), "NORMAL")

        casalog.post(stats + " statistics --- ", "NORMAL")
        if s[stats]['npts'] > 0:
                casalog.post("        -- Mean of the values [mean]:                 " + str(s[stats]['mean']), "NORMAL")
                casalog.post("        -- Variance of the values [var]:              " + str(s[stats]['var']), "NORMAL")
                casalog.post("        -- Standard deviation of the values [stddev]: " + str(s[stats]['stddev']), "NORMAL")
                casalog.post("        -- Root mean square [rms]:                    " + str(s[stats]['rms']), "NORMAL")
                casalog.post("        -- Median of the pixel values [median]:       " + str(s[stats]['median']), "NORMAL")
                casalog.post("        -- Median of the deviations [medabsdevmed]:   " + str(s[stats]['medabsdevmed']), "NORMAL")
                casalog.post("        -- Quartile [quartile]:                       " + str(s[stats]['quartile']), "NORMAL")
        else:
            casalog.post(stats + " -- No valid points found", "WARN")

    return s
