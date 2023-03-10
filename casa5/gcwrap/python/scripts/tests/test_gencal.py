from __future__ import absolute_import
from __future__ import print_function
import os
import sys
import shutil
import numpy as np
import numpy.ma as ma
import unittest

from casatasks.private.casa_transition import is_CASA6
from casatestutils import testhelper as th

if is_CASA6:
    from casatools import ctsys, table
    from casatasks import gencal, rmtables
    from casatasks.private import tec_maps    

    _tb= table()

    datapath=ctsys.resolve('/unittest/gencal/')
else:
    from __main__ import default
    from tasks import gencal, rmtables
    from taskinit import *
    from recipes import tec_maps
    
    _tb=tbtool()
    
    datapath=os.environ.get('CASAPATH').split()[0]+'/casatestdata/unittest/gencal/'

# input data
evndata = 'n08c1.ms'
vlbadata = 'ba123a.ms'
vlbacal = os.path.join(datapath,'ba123a.gc')
evncal = os.path.join(datapath,'n08c1.tsys')

caltab = 'cal.A'
evncopy = 'evn_copy.ms'
vlbacopy = 'vlba_copy.ms'

'''
Unit tests for gencal 
'''
#
# ToDo:
# add more tests
# once more independent tests (e.g. comparison
# the AIPS REWAY results) add reference mses
# and do tests against them
# 

# Pick up alternative data directory to run tests on MMSs
testmms = False
if 'TEST_DATADIR' in os.environ:   
    DATADIR = str(os.environ.get('TEST_DATADIR'))+'/gencal/'
    if os.path.isdir(DATADIR):
        testmms = True
        datapath = DATADIR
    else:
        print('WARN: directory '+DATADIR+' does not exist')

print('gencal tests will use data from '+datapath)         


class gencal_antpostest(unittest.TestCase):

    # Input and output names
    msfile = 'tdem0003gencal.ms'
    # used for test_antpos_auto_evla_CAS13057
    msfile2 = 'auto_antposcorr_evla_gencal.ms'
#    if testmms:
#        msfile = 'tdem0003gencal.mms'
    caltable = 'anpos.cal'
    reffile1 = os.path.join(datapath+'evla_reference/','anpos.manual.cal')
    reffile2 = os.path.join(datapath+'evla_reference/','anpos.auto.cal')
    reffile3 = os.path.join(datapath+'evla_reference/','anpos.autoCAS13057.cal')
    res = False

    def setUp(self):
        if (os.path.exists(self.msfile)):
            shutil.rmtree(self.msfile)
        if (os.path.exists(self.msfile2)):
            shutil.rmtree(self.msfile2)

        shutil.copytree(os.path.join(datapath,self.msfile), self.msfile, symlinks=True)
        shutil.copytree(os.path.join(datapath,self.msfile2), self.msfile2, symlinks=True)

    def tearDown(self):
        if (os.path.exists(self.msfile)):
            shutil.rmtree(self.msfile)
        if (os.path.exists(self.msfile2)):
            shutil.rmtree(self.msfile2)

        shutil.rmtree(self.caltable,ignore_errors=True)

    def test_antpos_manual(self):
        """
        gencal: test manual antenna position correction 
        """
        gencal(vis=self.msfile,
               caltable=self.caltable, 
               caltype='antpos',
               antenna='ea12,ea22',
               parameter=[-0.0072,0.0045,-0.0017, -0.0220,0.0040,-0.0190])

        self.assertTrue(os.path.exists(self.caltable))

        # ToDo:check generated caltable. Wait for new caltable
        
        # Compare with reference file from the repository
        reference = self.reffile1
        self.assertTrue(th.compTables(self.caltable, reference, ['WEIGHT','OBSERVATION_ID']))

    def test_antpos_auto_evla(self):
        """
        gencal: test automated antenna position correction
        """
        # check if the URL is reachable
        if is_CASA6:
            from urllib.request import urlopen
            from urllib.error import URLError
        else:
            from urllib2 import urlopen, URLError
    
        # current EVLA baseline correction URL
        evlabslncorrURL="http://www.vla.nrao.edu/cgi-bin/evlais_blines.cgi?Year="
        try: 
          urlaccess=urlopen(evlabslncorrURL+"2010", timeout=60.0) 
          gencal(vis=self.msfile,
                 caltable=self.caltable,
                 caltype='antpos',
                 antenna='')

          self.assertTrue(os.path.exists(self.caltable))
          
          # ToDo: check for generated caltable
          
          # Compare with reference file from the repository
          reference = self.reffile2
          self.assertTrue(th.compTables(self.caltable, reference, ['WEIGHT','OBSERVATION_ID']))

        except URLError as err:
          print("Cannot access %s , skip this test" % evlabslncorrURL)
          self.res=True

    def test_antpos_auto_evla_CAS13057(self):
        """
        gencal: test a bugfix of CAS-13057 for automated antenna position correction
        """
        # check if the URL is reachable
        if is_CASA6:
            from urllib.request import urlopen
            from urllib.error import URLError
        else:
            from urllib2 import urlopen, URLError
    
        # current EVLA baseline correction URL
        evlabslncorrURL="http://www.vla.nrao.edu/cgi-bin/evlais_blines.cgi?Year="
        try: 
          urlaccess=urlopen(evlabslncorrURL+"2019", timeout=60.0) 
          gencal(vis=self.msfile2,
                 caltable=self.caltable,
                 caltype='antpos',
                 antenna='')

          self.assertTrue(os.path.exists(self.caltable))
          # Compare with reference file from the repository
          reference = self.reffile3
          self.assertTrue(th.compTables(self.caltable, reference, ['WEIGHT','OBSERVATION_ID']))

        except URLError as err:
          print("Cannot access %s , skip this test" % evlabslncorrURL)
          self.res=True



class test_gencal_antpos_alma(unittest.TestCase):
    """
    Tests the automatic generation of antenna position corrections for ALMA

    New REST web service:
    https://bitbucket.sco.alma.cl/projects/ALMA/repos/almasw/browse/CONTROL-SERVICES/PositionsService


    Old SOAP web service:
    http://asa.alma.cl/axis2/services/TMCDBAntennaPadService?wsdl
    Example minimalistic use of a client to query the service:
      from suds.client import Client
      srv_wsdl_url = 'http://asa.alma.cl/axis2/services/TMCDBAntennaPadService?wsdl'
      ws_cli = Client(srv_wsdl_url)
      resp = ws_cli.service.getAntennaPositions("CURRENT.AOS", "DA49",
                                                "2017-01-30T01:53:54")
    """

    # setup of the ALMA TMC DB AntennaPadService
    ALMA_SRV_WSDL_URL = 'http://asa.alma.cl/axis2/services/TMCDBAntennaPadService?wsdl'

    # For this MS, there is position information for 25 out of the 29 antennas
    # (at 2013-11-15T10:26:19)
    ALMA_MS = 'uid___A002_X72c4aa_X8f5_scan21_spw18_field2_corrXX.ms'
    CAL_TYPE = 'antpos'
    REF_CALTABLE_MANUAL = os.path.join(datapath,'alma_reference/A002_X72c4aa_ref_ant_pos.manual.cal')
    REF_CALTABLE_AUTO = os.path.join(datapath,'alma_reference/A002_X72c4aa_ref_ant_pos.auto.cal')
    IGNORE_COLS = ['WEIGHT','OBSERVATION_ID']

    def setUp(self):
        if (os.path.exists(self.ALMA_MS)):
            shutil.rmtree(self.ALMA_MS)

        shutil.copytree(os.path.join(datapath, self.ALMA_MS),
                        self.ALMA_MS, symlinks=True)

    def tearDown(self):
        if (os.path.exists(self.ALMA_MS)):
            shutil.rmtree(self.ALMA_MS)

    def remove_caltable(self, ct_name):
        """ Removes a cal table. ct_name: path to the caltable """
        import shutil
        shutil.rmtree(ct_name)

    def test_antpos_alma_manual(self):
        """
        gencal: manual antenna position correction on ALMA table
        """

        out_caltable = 'ant_pos_man.cal'
        gencal(vis=self.ALMA_MS,
               caltable=out_caltable,
               caltype=self.CAL_TYPE,
               antenna='DV07,DV10,DV11',
               parameter=[-0.0072,0.0045,-0.0017, -0.0220,0.0040,-0.0190])

        self.assertTrue(os.path.exists(out_caltable),
                        "The output cal table should have been created")

        # Compare against ref file
        self.assertTrue(th.compTables(out_caltable,
                                      self.REF_CALTABLE_MANUAL,
                                      self.IGNORE_COLS))

        self.remove_caltable(out_caltable)

    @unittest.skip('SOAP AntennaPad Positions SOAP service needs to be removed once the '
                   'TMCDB based auto correction in gencal is validated.')
    def tmp_disabled_test_antpos_alma_server_SOAP_methods(self):
        """
        gencal: connection to alma TCM DB AntennaPadService for ALMA
        """
        try:
            # these imports don't work in CASA6 - test is being skipped so not important
            import urllib2
            from suds.client import Client
            ws_cli = Client(self.ALMA_SRV_WSDL_URL)

            # Basic check that the schema has the minimum requirement
            method_name = 'getAntennaPositions'
            self.assertTrue(callable(getattr(ws_cli.service, method_name)),
                            'The client service should have this method: {}, and '
                            'it should be callable.'.format(method_name))
        except ImportError as exc:
            print('Cannot import required dependencies to query the ALMA TCM DB web service')
            raise
        except urllib2.URLError as exc:
            print('Connection/network error while querying the ALMA TCM DB web service')
            raise

    @unittest.skip('SOAP AntennaPad Positions SOAP service needs to be removed once the '
                   'TMCDB based auto correction in gencal is validated.')
    def tmp_disabled_test_antpos_auto_alma_SOAP_empty_query(self):
        """
        gencal: empty query (empty antennas list) to the (old) SOAP TCMDB AntennaPadService
        web service (ALMA)
        """
        try:
            import correct_ant_posns_alma as almacor

            resp = almacor.query_tmcdb_antennas_rest([], '2017-01-01T16:53:54.000')
            if resp:
                raise RuntimeError('Unexpected response for an empty query: {0}'.
                                   format(resp))
        except ImportError:
            print('Cannot import required dependencies to query the ALMA TCM DB web service')
            raise
        except urllib2.URLError as exc:
            print('Connection/network error while querying the ALMA TCM DB web service')
            raise

    @unittest.skip('SOAP AntennaPad Positions SOAP service needs to be removed once the '
                   'TMCDB based auto correction in gencal is validated.')
    def tmp_disabled_test_antpos_auto_web_srv_SOAP_alma(self):
        """
        gencal: auto gencal using data from TCM DB AntennaPadService (ALMA)
        """

        import urllib2

        out_caltable = 'ant_pos_web_srv.cal'
        try:
            # This will import the required libraries, urllib2, suds, etc.
            # Coul also use additional parameters: antenna='', parameter=''
            gencal(vis=self.ALMA_MS, caltable=out_caltable, caltype=self.CAL_TYPE)
        except ImportError:
            print('Cannot import required dependencies to query the ALMA TCM DB web service')
            raise
        except urllib2.URLError:
            print('Connection/network error while querying the ALMA TCM DB web service')
            raise

        self.assertTrue(os.path.exists(out_caltable),
                        "The output cal table should have been created: {0}".
                        format(out_caltable))

        # Compare against ref file
        self.assertTrue(th.compTables(out_caltable,
                                      self.REF_CALTABLE_AUTO,
                                      self.IGNORE_COLS))
        self.remove_caltable(out_caltable)


    @unittest.skip('REST Position service needs validation and final deployment')
    def tmp_disabled_test_antpos_auto_alma_REST_empty_query(self):
        """
        gencal: empty query (empty antennas list) to the (new) REST TCMDB Positions
        web service (ALMA)
        """
        import urllib2

        TEST_HOSTNAME = 'https://2018may.asa-test.alma.cl'

        hostname = TEST_HOSTNAME
        port = 443
        api = 'antenna-position/position/antenna'
        try:
            import requests
            import correct_ant_posns_alma as almacor

            tstamp = '2017-01-01T16:53:54.000'
            # query via correct_ant_posns function
            resp = almacor.query_tmcdb_antennas_rest([], tstamp)
            if resp:
                raise RuntimeError('Unexpected response for an empty query: {0}'.
                                   format(resp))

            # query directly via requests
            url = '{}:{}/{}?antenna={}&timestamp={}'.format(hostname, port, api, '',
                                                            '2017-01-01T16:53:54.000')
            resp = requests.get(url)
            if resp:
                raise RuntimeError('Unexpected response for an empty query: {0}'.
                                   format(resp))
        except ImportError:
            print('Cannot import required dependencies to query the ALMA TCM DB web service')
            raise
        except urllib2.URLError as exc:
            print('Connection/network error while querying the ALMA TCM DB web service')
            raise

    @unittest.skip('REST Position service needs validation and final deployment')
    def tmp_disabled_test_antpos_auto_web_srv_REST_alma(self):
        """
        gencal: auto gencal using data from TCMDB Positions service (ALMA)
        """

        import urllib2

        out_caltable = 'ant_pos_web_srv.cal'
        try:
            # This will import the required libraries, urllib2, suds, etc.
            # Coul also use additional parameters: antenna='', parameter=''
            gencal(vis=self.ALMA_MS, caltable=out_caltable, caltype=self.CAL_TYPE)
        except urllib2.URLError:
            print('Connection/network error while querying the ALMA TCMDB Positions web service')
            raise

        self.assertTrue(os.path.exists(out_caltable),
                        "The output cal table should have been created: {0}".
                        format(out_caltable))

        # Compare against ref file
        self.assertTrue(th.compTables(out_caltable,
                                      self.REF_CALTABLE_AUTO,
                                      self.IGNORE_COLS))
        self.remove_caltable(out_caltable)

class gencal_test_tec_vla(unittest.TestCase):

    # Input and output names
    msfile = 'tdem0003gencal.ms'
    igsfile= 'igsg1160.10i'
    tecfile= msfile+'.IGS_TEC.im'
    rmstecfile= msfile+'.IGS_RMS_TEC.im'
    caltable= msfile+'_tec.cal'

    # NEAL: Please check that these setUp and tearDown functions are ok

    def setUp(self):
        self.tearDown()
        shutil.copytree(os.path.join(datapath,self.msfile), self.msfile, symlinks=True)

    def tearDown(self):
        if os.path.exists(self.msfile):
            shutil.rmtree(self.msfile)

        if os.path.exists(self.igsfile):
                os.remove(self.igsfile)

        shutil.rmtree(self.tecfile,ignore_errors=True)
        shutil.rmtree(self.rmstecfile,ignore_errors=True)
        shutil.rmtree(self.caltable,ignore_errors=True)

    def test_tec_maps(self):
        """
        gencal: very basic test of tec_maps and gencal(caltype='tecim')
        """

        try:
            tec_maps.create0(self.msfile)
            gencal(vis=self.msfile, caltable=self.caltable, caltype='tecim',infile=self.msfile+'.IGS_TEC.im')

            self.assertTrue(os.path.exists(self.caltable))

            _tb.open(self.caltable)
            nrows=_tb.nrows()
            dtecu=abs(13.752-np.mean(_tb.getcol('FPARAM'))/1e16)
            _tb.close()
            
            #print(str(nrows)+' '+str(dtecu))

            self.assertTrue(nrows==1577)
            self.assertTrue(dtecu<1e-3)
            
            
        except:
            # should catch case of internet access failure?
            raise

class gencal_gaincurve_test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        shutil.copytree(os.path.join(datapath,evndata), evncopy)
        shutil.copytree(os.path.join(datapath,vlbadata), vlbacopy)

    def setUp(self):
        if not is_CASA6:
            default(gencal)

    def tearDown(self):
        rmtables(caltab)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(evncopy)
        shutil.rmtree(vlbacopy)

    def test_gainCurve(self):
        ''' Test calibration table produced when gencal is run on an MS with a GAIN_CURVE table '''

        gencal(vis=vlbacopy, caltable=caltab, caltype='gc')

        self.assertTrue(os.path.exists(caltab))
        self.assertTrue(th.compTables(caltab, vlbacal, ['WEIGHT']))

    def test_noGainCurve(self):
        ''' Test that when gencal is run on an MS with no GAIN_CURVE table it creates no calibration table '''

        try:
            gencal(vis=evncopy, caltable=caltab, caltype='gc')
        except:
            pass

        self.assertFalse(os.path.exists(caltab))

class gencal_tsys_test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        shutil.copytree(os.path.join(datapath,evndata), evncopy)
        shutil.copytree(os.path.join(datapath,vlbadata), vlbacopy)

    def setUp(self):
        if not is_CASA6:
            default(gencal)

    def tearDown(self):
        rmtables(caltab)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(evncopy)
        shutil.rmtree(vlbacopy)

    def test_tsys(self):
        ''' Test calibration table produced when gencal is run on an MS with a SYSCAL table'''

        gencal(vis=evncopy, caltable=caltab, caltype='tsys', uniform=False)

        self.assertTrue(os.path.exists(caltab))
        self.assertTrue(th.compTables(caltab, evncal, ['WEIGHT']))

    def test_tsys_nan(self):
        ''' Test calibration table produced when gencal is run on an MS with a SYSCAL table that contains NaNs'''

        # Change negative values in SYSCAL to NaNs.
        # This should result in the same calibration table entries
        # being flagged.
        _tb.open(evncopy + '/SYSCAL', nomodify=False)
        tsys = _tb.getcol('TSYS')
        tsys = np.where(tsys < 0, float('nan'), tsys)
        _tb.putcol('TSYS', tsys)
        _tb.close()

        gencal(vis=evncopy, caltable=caltab, caltype='tsys', uniform=False)

        self.assertTrue(os.path.exists(caltab))
        self.assertTrue(th.compTables(caltab, evncal, ['FPARAM', 'WEIGHT']))


def suite():
    return [gencal_antpostest,
            test_gencal_antpos_alma,
            gencal_test_tec_vla,
            gencal_gaincurve_test,
            gencal_tsys_test]

if is_CASA6:
    if __name__ == '__main__':
        unittest.main()
