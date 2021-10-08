########################################################################
# Task to print a summary of the content of an SDM dataset
#
# v1.0: 2012.02.13, M. Caillat
#
import os
from casatools import sdm
from casatasks import casalog

def asdmsummary(asdm=None):
    """Prints a description of the content of an SDM dataset to the CASA logger.

    Keyword argument:

    asdm -- Name of the input SDM directory.

    """
    casalog.origin('asdmsummary')
    dataset = sdm(asdm)
    summary = dataset.summarystr( )
    if len(summary) == 1:
        casalog.post('error generating summary: ' + summary[0], 'SEVERE')
    else:
        casalog.post(summary)
