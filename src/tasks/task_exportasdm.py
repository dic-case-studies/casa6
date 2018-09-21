import os
from casatasks import casalog
from casatools import sdm, table, quanta, ms

def exportasdm(vis=None, asdm=None, datacolumn=None, archiveid=None, rangeid=None,
               subscanduration=None, sbduration=None, apcorrected=None,
               verbose=None, showversion=None, useversion=None):
        """ Convert a CASA visibility file (MS) into an ALMA or EVLA Science Data Model.
                                          
        Keyword arguments:
        vis       -- MS name,
             default: none

        asdm -- Name of output ASDM file (directory),
             default: none; example: asdm='ExecBlock3'

        datacolumn -- specifies which of the MS data columns (DATA,
                  CORRECTED_DATA, or MODEL_DATA) should be used as the
                  visibilities in the ASDM, default: DATA

        archiveid -- the X0 in uid://X0/X1/X<running>
                  default: "S0"

        rangeid -- the X1 in uid://X0/X1/X<running>
                  default: "X1"

        subscanduration -- maximum duration of a subscan in the output ASDM
                  default: "24h"

        sbduration -- maximum duration of a scheduling block in the output ASDM
                  default: "2700s"

        apcorrected -- If true, the data in column datacolumn should be regarded
                  as having atmospheric phase correction, default: False

        verbose     -- produce log output, default: True

        showversion -- report the version of the ASDM class set, 
                 default: True

        useversion -- Selects the version of MS2asdm to be used (presently only \'v3\' is available)
              default: v3

        """
        #Python script

        # make ms and tb tool local 
        _ms = ms( )
        _tb = table( )
        _qa = quanta( )

        try:
                casalog.origin('exportasdm')

                if not (type(vis)==str) or not (os.path.exists(vis)):
                        raise Exception('Visibility data set not found - please verify the name')
                
                if (asdm == ""):
                        raise Exception("Must provide output data set name in parameter asdm.")
                
                if os.path.exists(asdm):
                        raise Exception("Output ASDM %s already exists - will not overwrite." % asdm)

                # determine sb and subscan duration
                ssdur_secs = 24.*3600 # default is one day, i.e. there will be only one subscan per scan
                if not(subscanduration==""):
                        if (_qa.canonical(subscanduration)['unit'].find('s') < 0):
                                raise TypeError("subscanduration is not a valid time quantity: %s" % subscanduration)
                        else:
                                ssdur_secs = _qa.canonical(subscanduration)['value']

                sbdur_secs = 2700. # default is 45 minutes
                if not(sbduration==""):
                        if (_qa.canonical(sbduration)['unit'].find('s') < 0):
                                raise TypeError("sbduration is not a valid time quantity: %s" % sbduration)
                        else:
                                sbdur_secs = _qa.canonical(sbduration)['value']

                # create timesorted copy of the input ms
                os.system('rm -rf '+vis+'-tsorted')
                _ms.open(vis)
                _ms.timesort(vis+'-tsorted')
                _ms.close()

                tsortvis = vis+'-tsorted'

                # Prepare for actual exportasdm
                casalog.post("Checking timesorted MS for potential problems ... ")
                _tb.open(tsortvis)
                a = _tb.getcol('PROCESSOR_ID')
                a0 = a[0]
                candoit = True
                for i in range(0,len(a)-1):
                        if(a[i]!=a[i+1]):
                                candoit = False
                                break
                _tb.close()

                if candoit:
                        casalog.post("   Checking if PROCESSOR and MAIN need modifications ...")
                        _tb.open(tsortvis+'/PROCESSOR')
                        nprocrows = _tb.nrows()
                        _tb.close()
                        if ((nprocrows>0) and (a0>-1)):
                                _tb.open(tsortvis+'/PROCESSOR')
                                therow = _tb.nrows()-1 
                                mode0 = _tb.getcell('MODE_ID',a0)
                                _tb.close()
                                offset = 1
                                if nprocrows>1:
                                        casalog.post("   Modifying PROCESSOR subtable ...")
                                        while (nprocrows>1 and therow>0):
                                                _tb.open(tsortvis+'/PROCESSOR', nomodify=False)
                                                therow = _tb.nrows()-offset
                                                if(_tb.getcell('MODE_ID',therow)!=mode0):
                                                        _tb.removerows([therow])
                                                else:
                                                        offset += 1
                                                nprocrows = _tb.nrows()
                                        casalog.post("... done.")

                                        casalog.post("   Modifying processor ids in main table ...")
                                        a = a - a # set all precessor ids to zero
                                        _tb.open(tsortvis, nomodify=False)
                                        _tb.putcol('PROCESSOR_ID', a)
                                        _tb.close()
                                        casalog.post(" ... done.")
                                else:
                                        casalog.post("   No modifications to proc id in PROCESSOR and MAIN necessary.")
                        casalog.post("   Checking if SYSCAL needs modifications ...")
                        if(os.path.exists(tsortvis+'/SYSCAL')):
                                for cname in ['TANT_SPECTRUM',
                                              'TSYS_SPECTRUM',
                                              'TANT_TSYS_SPECTRUM',
                                              'TCAL_SPECTRUM',
                                              'TRX_SPECTRUM',
                                              'TSKY_SPECTRUM',
                                              'PHASE_DIFF_SPECTRUM']:
                                        _tb.open(tsortvis+'/SYSCAL', nomodify=False)
                                        if(cname in _tb.colnames()):
                                                cdesc = _tb.getcoldesc(cname)
                                                if 'ndim' in cdesc and (cdesc['ndim']==-1):
                                                        _tb.removecols([cname])
                                                        casalog.post('   Removed empty array column '+cname+' from table SYSCAL.')
                                        _tb.close()

                        casalog.post("   Checking if OBSERVATION needs modifications ...")
                        _tb.open(tsortvis+'/OBSERVATION')
                        nobsrows = _tb.nrows()
                        _tb.close()
                        if(nobsrows>0):
                                _tb.open(tsortvis+'/OBSERVATION', nomodify=False)
                                cdesc = _tb.getcoldesc('LOG')
                                if 'ndim' in cdesc and (cdesc['ndim']>0):
                                        b = _tb.getvarcol('LOG')
                                        if not (type(b['r1'])==bool):
                                                kys = b.keys()
                                                modified = False
                                                b2 = []
                                                for i in range(0,len(kys)):
                                                        k = 'r'+str(i+1)
                                                        if (b[k][0] == [''])[0]:
                                                                b[k][0] = ["-"]
                                                                modified = True
                                                        b2.append([b[k][0][0]])
                                                if modified:
                                                        _tb.putcol('LOG',b2)
                                                        casalog.post("   Modified log column in OBSERVATION table.")
                                _tb.close()
                        casalog.post("Done.")
                else:
                        raise Exception("More than one processor id in use in the main table. Cannot proceed.")

                _sdm = sdm(asdm)
                rval = _sdm.fromms(vis, datacolumn, archiveid, rangeid, ssdur_secs, sbdur_secs, apcorrected, verbose)
                os.system('rm -rf '+vis+'-tsorted')
                return rval
        
        except Exception as instance:
                casalog.post("Error: %s" % instance, 'SEVERE')
                raise
