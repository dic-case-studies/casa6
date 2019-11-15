from __future__ import absolute_import
import os
import glob

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
       from casatools import ctsys
       from casatasks import casalog
else:
       from taskinit import *
       from tasksinfo import *

def rmtables(tablenames=None):
       """ Removes tables cleanly 
           Arguments may contain *, ?. Ranges [] also supported but not ~ expansion.
       """

       casalog.origin('rmtables')
       tablelist = []
       for tables in tablenames :
          for table in glob.glob(tables) :
             tablelist.append(table)
       for table in tablelist :
          casalog.post('Removing '+table)
       if is_CASA6:
              ctsys.removetable(tablelist)
       else:
              cu.removetable(tablelist)

