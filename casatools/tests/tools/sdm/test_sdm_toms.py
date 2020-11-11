import os
import sys
import shutil
import numpy
from casatools import ctsys, sdm, ms, table, quanta, measures, calibrater
### for testhelper import
sys.path.append(os.path.dirname(os.path.abspath(os.path.dirname(__file__))))
import testhelper as th
from xmlhelper import readXML
import unittest

myname = 'test_sdm_toms'

datapath=ctsys.resolve('unittest/sdmtool/')

# default ASDM dataset name
myasdm_dataset_name = 'uid___X5f_X18951_X1'
myms_dataset_name = 'M51.ms'

# name of the resulting MS
msname = myasdm_dataset_name+'.ms'

# name of the exported ASDM
asdmname = myms_dataset_name+'.asdm'

# name of the reimported MS
reimp_msname = 'reimported-'+myms_dataset_name

# make local copies of the tools
tblocal = table()
mslocal = ms()

# imported from casatasks convertephem
# only used in test7_lazy3 which is currently disabled because that option from asdm2MS
# is not (yet?) available in sdm.toms
def converttopoephem2geo(tablename='', outtablename='', overwrite=True):
    """
    converttopoephem2geo

    Convert the given topo ephemeris table to the geocentric ref frame.
    Converted are the RA, DEC, and RadVel columns only

    tablename -- name of the TOPO frame ephemeris table (in canonical CASA format)

    outtablename -- name of the output GEO frame ephemeris table

    """

    if type(tablename)!=str or tablename=='':
        casalog.post('Invalid parameter tablename', 'WARN')
        return False

    if type(outtablename)!=str or outtablename=='':
        casalog.post('Invalid parameter outtablename', 'WARN')
        return False

    #convert RA, DEC, and RadVel
    tbt = table( )
    met = measures( )
    qat = quanta( )

    tbt.open(tablename)
    ra = tbt.getcol('RA')
    dec = tbt.getcol('DEC')
    radvel = tbt.getcol('RadVel')
    radvelunit = 'km/s'
    tmpkw = tbt.getcolkeywords('RadVel')
    if tmpkw.has_key('UNIT'):
        radvelunit = tmpkw['UNIT']
    elif tmpkw.has_key('QuantumUnits'):
        radvelunit = tmpkw['QuantumUnits'][0]
    else:
        casalog.post('Cannot determine units of radial velocity column. Assuming km/s.', 'WARN')
    mjd = tbt.getcol('MJD')
    kw = tbt.getkeywords()
    tbt.close()

    geodist = kw['GeoDist'] # (km above reference ellipsoid)
    geolat = kw['GeoLat'] # (deg)
    geolong = kw['GeoLong'] # (deg)

    if kw.has_key('obsloc'):
        obsloc = kw['obsloc']
    else:
        casalog.post('Ephemeris does not have the obsloc keyword.', 'INFO')
        if (geodist==geolat==geolong==0.):
            casalog.post('   Assuming obsloc == GEOCENTRIC since lat, long, and dist are zero.', 'INFO')
            obsloc='GEOCENTRIC'
        else:
            obsloc='unknown'

    oldref = 'J2000'
    newref = 'ICRS'

    if kw.has_key('posrefsys'):
        posref = kw['posrefsys']
    else:
        casalog.post('Ephemeris does not have the posrefsys keyword. Assuming ICRF/J2000.0', 'WARN')
        posref = 'ICRF/J2000.0'

    if not ('J2000' in posref):
        if 'ICRS' in posref:
            oldref = 'ICRS'
        else:
            casalog.post('Position reference is '+posref+' is not supported, yet.', 'WARN')
            return False

    newposref = 'ICRF/'+newref

    if oldref!='ICRS':
        casalog.post('Position reference is '+oldref+'. Will convert positions to '+newref, 'INFO')

    if obsloc=='GEOCENTRIC':
        casalog.post('Obsloc is already GEOCENTRIC. Nothing to be done.', 'INFO')
        return True

    mepos = {'m0': {'value': geolong, 'unit': 'deg'},
             'm1': {'value': geolat, 'unit': 'deg'}, # latitude
             'm2': {'value': geodist, 'unit': 'km'}, # alt above ref ellipsoid
             'refer': 'WGS84',
             'type': 'position'}
    met.doframe(mepos)

    newra=[]
    newdec=[]
    newrho=[]
    newradvel=[]

    for i in range(0, len(mjd)):
        memjd={'m0': {'value': mjd[i], 'unit': 'd'},
               'refer': 'UTC',
               'type': 'epoch'}
        met.doframe(memjd)

        olddir={'m0': {'value': ra[i], 'unit': 'deg'},
                'm1': {'value': dec[i], 'unit': 'deg'},
                'refer': oldref,
                'type': 'direction'}
        met.doframe(olddir)
        newdir=met.measure(olddir, newref)

        tmpnewra = qat.convert(newdir['m0'],'deg')['value']
        if tmpnewra<0:
            tmpnewra += 360.
        newra.append(tmpnewra)
        newdec.append(qat.convert(newdir['m1'],'deg')['value'])

        oldradvel={'m0': {'value': radvel[i], 'unit': radvelunit},
                   'refer': 'TOPO',
                   'type': 'radialvelocity'}

        newradvelme = met.measure(oldradvel, 'GEO')
        newradvel.append(qat.convert(newradvelme['m0'], radvelunit)['value'])

    # create the converted table
    safetycopyname=tablename+'.orig'
    if overwrite:
        if outtablename==tablename:
            os.system('cp -R '+tablename+' '+safetycopyname)
        else:
            os.system('rm -rf '+outtablename)
            os.system('cp -R '+tablename+' '+outtablename)
    else:
        if os.path.exists(outtablename):
            casalog.post('Output table '+outtablename+' already exists.', 'WARN')
            return False
        os.system('cp -R '+tablename+' '+outtablename)

    try:
        tbt.open(outtablename, nomodify=False)
        tbt.putcol('RA', newra)
        tbt.putcol('DEC', newdec)
        tbt.putcol('RadVel', newradvel)
        tbt.putkeyword('GeoDist', 0.)
        tbt.putkeyword('GeoLat', 0.)
        tbt.putkeyword('GeoLong', 0.)
        tbt.putkeyword('obsloc', 'GEOCENTRIC')
        tbt.putkeyword('posrefsys', newposref)
        tbt.close()
    except Exception as instance:
        casalog.post("*** Error \'%s\' " % (instance), 'SEVERE')
        if overwrite and outtablename==tablename:
            casalog.post('Conversion in situ was not possible. Restoring original ephemeris ...', 'INFO')
            os.system('rm -rf '+tablename)
            os.system('mv '+safetycopyname+' '+tablename)
        return False

    if overwrite and outtablename==tablename:
        os.system('rm -rf '+safetycopyname)

    return True

def checktable(thename, theexpectation):
    global msname, myname
    tblocal.open(msname+"/"+thename)
    if thename == "":
        thename = "MAIN"
    for mycell in theexpectation:
        print("%s: comparing %s" % (myname, mycell))
        value = tblocal.getcell(mycell[0], mycell[1])
        in_agreement = False
        
        # bool and str cases
        if type(value) is bool or type(value) is str:
            # must be equal
            in_agreement = (value == mycell[2])
        else:
            # see if value is array
            try:
                isarray = value.__len__
            except:
                # it's not an array
                # zero tolerance? also works for bool and string scalar values
                if mycell[3] == 0:
                    in_agreement = (value == mycell[2])
                else:
                    in_agreement = ( abs(value - mycell[2]) < mycell[3])
            else:
                # it's an array
                # zero tolerance?
                if mycell[3] == 0:
                    in_agreement =  (value == mycell[2]).all()
                else:
                    in_agreement = (abs(value - mycell[2]) < mycell[3]).all()
        if not in_agreement:
            print("%s:  Error in MS subtable %s:" % (myname, thename))
            print("     column %s row %s contains %s" % (mycell[0], mycell[1], value))
            print("     expected value is %s" % mycell[2])
            tblocal.close()
            return False
    tblocal.close()
    print("%s: table %s as expected." % (myname, thename))
    return True

#########################

def verify_asdm(asdmname, withPointing):
    print("Verifying asdm %s" % asdmname)
    if(not os.path.exists(asdmname)):
        print("asdm %s doesn't exist." % asdmname)
        raise Exception
    # test for the existence of all obligatory tables
    allTables = [ "Antenna.xml",
                  "ASDM.xml",
                 # "CalData.xml",
                 # "CalDelay.xml",
                 # "CalReduction.xml",
                  "ConfigDescription.xml",
                  "CorrelatorMode.xml",
                  "DataDescription.xml",
                  "ExecBlock.xml",
                  "Feed.xml",
                  "Field.xml",
                 #"FocusModel.xml",
                 #"Focus.xml",
                  "Main.xml",
                  "PointingModel.xml",
                  "Polarization.xml",
                  "Processor.xml",
                  "Receiver.xml",
                  "SBSummary.xml",
                  "Scan.xml",
                  "Source.xml",
                  "SpectralWindow.xml",
                  "State.xml",
                  "Station.xml",
                  "Subscan.xml",
                  "SwitchCycle.xml"
                  ]
    isOK = True
    for fileName in allTables:
        filePath = asdmname+'/'+fileName
        if(not os.path.exists(filePath)):
            print("ASDM table file %s doesn't exist." % filePath)
            isOK = False
        else:
            # test if well formed
            rval = os.system('xmllint --noout '+filePath)
            if(rval !=0):
                print("Table %s is not a well formed XML document." % filePath)
                isOK = False

    print("Note: xml validation not possible since ASDM DTDs (schemas) not yet online.")

    if(not os.path.exists(asdmname+"/ASDMBinary")):
        print("ASDM binary directory %s/ASDMBinary doesn't exist." % asdmname)
        isOK = False

    if(withPointing and not os.path.exists(asdmname+"/Pointing.bin")):
        print("ASDM binary file %s/Pointing.bin doesn't exist." % asdmname)
        isOK = False

    if (not isOK):
        raise Exception

# Setup for different data importing
class test_base(unittest.TestCase):
    def setUp_m51(self):
        res = None
        if(os.path.exists(myasdm_dataset_name)):
            shutil.rmtree(myasdm_dataset_name)

        shutil.copytree(os.path.join(datapath,myasdm_dataset_name), myasdm_dataset_name)
        shutil.copytree(os.path.join(datapath,myms_dataset_name), myms_dataset_name)

    def setUp_xosro(self):
        self.asdm = 'X_osro_013.55979.93803716435'
        if(not os.path.lexists(self.asdm)):
            os.system('ln -s '+datapath+self.asdm)

    def setUp_autocorr(self):
        self.asdm = 'AutocorrASDM'
        if(not os.path.lexists(self.asdm)):
            os.system('ln -s '+datapath+self.asdm +' '+self.asdm)

    def setUp_acaex(self):
        res = None
        myasdmname = 'uid___A002_X72bc38_X000' # ACA example ASDM with mixed pol/channelisation
        os.system('ln -sf '+datapath+myasdmname)

    def setUp_12mex(self):
        res = None
        myasdmname = 'uid___A002_X71e4ae_X317_short' # 12m example ASDM with mixed pol/channelisation
        os.system('ln -sf '+datapath+myasdmname)

    def setUp_eph(self):
        res = None
        myasdmname = 'uid___A002_X997a62_X8c-short' # 12m example ASDM with ephemerides
        os.system('ln -sf '+datapath+myasdmname)

    def setUp_flags(self):
        res = None
        myasdmname = 'test_uid___A002_X997a62_X8c-short' # Flag.xml is modified
        os.system('ln -sf '+datapath+myasdmname)

    def setUp_SD(self):
        res = None
        myasdmname = 'uid___A002_X6218fb_X264' # Single-dish ASDM
        os.system('ln -sf '+datapath+myasdmname)

    def setUp_numbin(self):
        res = None
        # need full copies as this test involves renaming some xml files
        for this_asdm_name in ['alma_numbin_mixed','evla_numbin_2','evla_numbin_4']:
            if (os.path.exists(this_asdm_name)):
                shutil.rmtree(this_asdm_name)
            shutil.copytree(os.path.join(datapath,this_asdm_name), this_asdm_name)

###########################
# beginning of actual test
class asdm_import1(test_base):

    def setUp(self):
        self.setUp_m51()

    def tearDown(self):
        os.system('rm -rf '+myasdm_dataset_name)
        shutil.rmtree(myms_dataset_name)
        shutil.rmtree(msname,ignore_errors=True)
        shutil.rmtree(msname+'.flagversions',ignore_errors=True)
        for thisdir in ['reimported-M51.ms','reimported-M51.ms.flagversions','M51.ms.asdm','myinput.ms']:
            shutil.rmtree(thisdir,ignore_errors=True)
    
    def test1(self):
        '''Asdm-import: Test good v1.2 input with filler v3 and inverse filler v3 '''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }

        mysdm = sdm(myasdm_dataset_name)
        self.res = mysdm.toms(useversion='v3')
        self.assertTrue(self.res)
        print( "%s: Success! Now checking output ..." % myname)
        mscomponents = set(["table.dat",
#                            "table.f0",
                            "table.f1",
                            "table.f2",
                            "table.f3",
                            "table.f4",
                            "table.f5",
                            "table.f6",
                            "table.f7",
                            "table.f8",
                            "ANTENNA/table.dat",
                            "DATA_DESCRIPTION/table.dat",
                            "FEED/table.dat",
                            "FIELD/table.dat",
                            "FLAG_CMD/table.dat",
                            "HISTORY/table.dat",
                            "OBSERVATION/table.dat",
                            "POINTING/table.dat",
                            "POLARIZATION/table.dat",
                            "PROCESSOR/table.dat",
                            "SOURCE/table.dat",
                            "SPECTRAL_WINDOW/table.dat",
                            "STATE/table.dat",
                            "SYSCAL/table.dat",
                            "ANTENNA/table.f0",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.f0",
                            "FIELD/table.f0",
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.f0",
                            "POINTING/table.f0",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.f0",
                            "SOURCE/table.f0",
                            "SPECTRAL_WINDOW/table.f0",
                            "STATE/table.f0",
                            "SYSCAL/table.f0"
                            ])
        for name in mscomponents:
            if not os.access(msname+"/"+name, os.F_OK):
                print("%s: Error  %s/%s doesn't exist..." % (myname,msname,name))
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+msname+'/'+name+' does not exist'
            else:
                print("%s: %s present." % (myname,name))
        print("%s: MS exists. All tables present. Try opening as MS ..." % myname)
        try:
            mslocal.open(msname)
        except:
            print("%s: Error Cannot open MS table %s" % (myname, msname))
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+msname
        else:
            mslocal.close()
            print("%s: OK. Checking tables in detail ..." % myname)

            # check main table first
            name = ""
            #             col name, row number, expected value, tolerance
            expected = [
                         ['UVW',       42, [ 0., 0., 0. ], 1E-8],
                         ['EXPOSURE',  42, 1.008, 0],
                         ['DATA',      42, [ [10.5526886+0.0j] ], 1E-7]
                         ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table MAIN failed'

            expected = [
    # old values using TAI     ['UVW',       638, [-65.07623467,   1.05534109, -33.65801386], 1E-8],
                         ['UVW',       638, [-65.14758508, 1.13423277, -33.51712451], 1E-7],
                         ['EXPOSURE',  638, 1.008, 0],
                         ['DATA',      638, [ [0.00362284+0.00340279j] ], 1E-8]
                         ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table MAIN failed'

            name = "ANTENNA"
            expected = [ ['OFFSET',       1, [ 0.,  0.,  0.], 0],
                         ['POSITION',     1, [2202242.5520, -5445215.1570, -2485305.0920], 0.0001],
                         ['DISH_DIAMETER',1, 12.0, 0]
                         ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table ANTENNA failed'

            name = "POINTING"
            expected = [ ['DIRECTION',       10, [[ 1.94681283],[ 1.19702955]], 1E-8],
                         ['INTERVAL',        10, 0.048, 0],
                         ['TARGET',          10, [[ 1.94681283], [ 1.19702955]], 1E-8],
                         ['TIME',            10, 4758823736.016000, 1E-6],
                         ['TIME_ORIGIN',     10, 0., 0],
                         ['POINTING_OFFSET', 10, [[ 0.],[ 0.]], 0],
                         ['ENCODER',         10, [ 1.94851533,  1.19867576], 1E-8 ]
                         ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table POINTING failed'

            name = "PROCESSOR"
            expected = [ ['FLAG_ROW',            0, False, 0],
                         ['MODE_ID',             0, 0, 0],
                         ['TYPE',                0, 'CORRELATOR', 0],
                         ['TYPE_ID',             0, -1, 0],
                         ['SUB_TYPE',            0, 'ALMA_BASELINE', 0]
                       ]
            results = checktable(name, expected)
            if not results:
                retValue['success'] = False
                retValue['error_msgs'] = retValue['error_msgs']+'Check of table PROCESSOR failed'
                
        self.assertTrue(retValue['success'],retValue['error_msgs'])

        myvis = myms_dataset_name
        os.system('rm -rf exportasdm-output.asdm myinput.ms')
        os.system('cp -R ' + myvis + ' myinput.ms')
        try:
            print("\n>>>> Test of exportasdm v3: input MS is %s" % myvis)
            print("(a simulated input MS with pointing table)")
            tosdm = sdm('exportasdm-output.asdm')
            rval = tosdm.fromms( 'myinput.ms', archiveid="S002", apcorrected=False, useversion='v3' )
            if not rval:
                raise Exception
            os.system('rm -rf '+asdmname+'; mv exportasdm-output.asdm '+asdmname)
            verify_asdm(asdmname, True)
        except:
            print('%s: *** Unexpected error exporting MS to ASDM, regression failed ***' % myname)
            raise

        try:
            print("Reimporting the created ASDM (v3)....")
            fromsdm = sdm(asdmname)
            fromsdm.toms(vis=reimp_msname, wvr_corrected_data='no', useversion='v3')
            print("Testing existence of reimported MS ....")
            if(not os.path.exists(reimp_msname)):
                print("MS %s doesn't exist." % reimp_msname)
                raise Exception
            print("Testing equivalence of the original and the reimported MS.")
            tblocal.open(myms_dataset_name)
            nrowsorig = tblocal.nrows()
            print("Original MS contains %s integrations." % nrowsorig)
            tblocal.close()
            tblocal.open(reimp_msname)
            nrowsreimp = tblocal.nrows()
            tblocal.close()
            print("Reimported MS contains %s integrations." % nrowsreimp)
            if(not nrowsreimp==nrowsorig):
                print("Numbers of integrations disagree.")
                raise Exception
        except:
            print('%s: *** Unexpected error reimporting the exported ASDM, regression failed ***' % myname)
            raise

class asdm_import2(test_base):

    def setUp(self):
        self.setUp_m51()

    def tearDown(self):
        shutil.rmtree(myasdm_dataset_name)
        shutil.rmtree(myms_dataset_name)
        shutil.rmtree(msname,ignore_errors=True)
        shutil.rmtree(msname+'.flagversions',ignore_errors=True)
        for thisdir in ['reimported-M51.ms','reimported-M51.ms.flagversions','M51.ms.asdm','myinput.ms']:
            shutil.rmtree(thisdir,ignore_errors=True)

    def test_import2(self):
        '''Asdm-import: Test good v1.2 input with filler v3 and inverse filler v3 '''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }

        mysdm = sdm(myasdm_dataset_name)
        self.res = mysdm.toms(useversion='v3')
        self.assertTrue(self.res)
        print("%s: Success! Now checking output ..." % myname)
        mscomponents = set(["table.dat",
#                            "table.f0",
                            "table.f1",
                            "table.f2",
                            "table.f3",
                            "table.f4",
                            "table.f5",
                            "table.f6",
                            "table.f7",
                            "table.f8",
                            "ANTENNA/table.dat",
                            "DATA_DESCRIPTION/table.dat",
                            "FEED/table.dat",
                            "FIELD/table.dat",
                            "FLAG_CMD/table.dat",
                            "HISTORY/table.dat",
                            "OBSERVATION/table.dat",
                            "POINTING/table.dat",
                            "POLARIZATION/table.dat",
                            "PROCESSOR/table.dat",
                            "SOURCE/table.dat",
                            "SPECTRAL_WINDOW/table.dat",
                            "STATE/table.dat",
                            "SYSCAL/table.dat",
                            "ANTENNA/table.f0",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.f0",
                            "FIELD/table.f0",
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.f0",
                            "POINTING/table.f0",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.f0",
                            "SOURCE/table.f0",
                            "SPECTRAL_WINDOW/table.f0",
                            "STATE/table.f0",
                            "SYSCAL/table.f0"
                            ])
        for name in mscomponents:
            if not os.access(msname+"/"+name, os.F_OK):
                print("%s: Error  %s/%s doesn't exist ..." % (myname,msname,name))
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+msname+'/'+name+' does not exist'
            else:
                print("%s: %s present." % (myname,name))
        print("%s: MS exists. All tables present. Try opening as MS ..." % myname)
        try:
            mslocal.open(msname)
        except:
            print("%s: Error Cannot open MS table %s" % (myname,msname))
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+msname
        else:
            mslocal.close()
            print("%s: OK. Checking tables in detail ..." % myname)

            # check main table first
            name = ""
            #             col name, row number, expected value, tolerance
            expected = [
                         ['UVW',       42, [ 0., 0., 0. ], 1E-8],
                         ['EXPOSURE',  42, 1.008, 0],
                         ['DATA',      42, [ [10.5526886+0.0j] ], 1E-7]
                         ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table MAIN failed'

            expected = [
    # old values using TAI     ['UVW',       638, [-65.07623467,   1.05534109, -33.65801386], 1E-8],
                         ['UVW',       638, [-65.14758508, 1.13423277, -33.51712451], 1E-7],
                         ['EXPOSURE',  638, 1.008, 0],
                         ['DATA',      638, [ [0.00362284+0.00340279j] ], 1E-8]
                         ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table MAIN failed'

            name = "ANTENNA"
            expected = [ ['OFFSET',       1, [ 0.,  0.,  0.], 0],
                         ['POSITION',     1, [2202242.5520, -5445215.1570, -2485305.0920], 0.0001],
                         ['DISH_DIAMETER',1, 12.0, 0]
                         ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table ANTENNA failed'

            name = "POINTING"
            expected = [ ['DIRECTION',       10, [[ 1.94681283],[ 1.19702955]], 1E-8],
                         ['INTERVAL',        10, 0.048, 0],
                         ['TARGET',          10, [[ 1.94681283], [ 1.19702955]], 1E-8],
                         ['TIME',            10, 4758823736.016000, 1E-6],
                         ['TIME_ORIGIN',     10, 0., 0],
                         ['POINTING_OFFSET', 10, [[ 0.],[ 0.]], 0],
                         ['ENCODER',         10, [ 1.94851533,  1.19867576], 1E-8 ]
                         ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table POINTING failed'

            name = "PROCESSOR"
            expected = [ ['FLAG_ROW',            0, False, 0],
                         ['MODE_ID',             0, 0, 0],
                         ['TYPE',                0, 'CORRELATOR', 0],
                         ['TYPE_ID',             0, -1, 0],
                         ['SUB_TYPE',            0, 'ALMA_BASELINE', 0]
                       ]
            results = checktable(name, expected)
            if not results:
                retValue['success'] = False
                retValue['error_msgs'] = retValue['error_msgs']+'Check of table PROCESSOR failed'
                
        self.assertTrue(retValue['success'],retValue['error_msgs'])

        myvis = myms_dataset_name
        os.system('rm -rf exportasdm-output.asdm myinput.ms')
        os.system('cp -R ' + myvis + ' myinput.ms')
        try:
            print("\n>>>> Test of exportasdm v3: input MS  is %s" % myvis)
            print("(a simulated input MS with pointing table)")
            tosdm = sdm('exportasdm-output.asdm')
            rval = tosdm.fromms( 'myinput.ms', archiveid="S002", apcorrected=False, useversion='v3' )
            if not rval:
                raise Exception
            os.system('rm -rf '+asdmname+'; mv exportasdm-output.asdm '+asdmname)
            verify_asdm(asdmname, True)
        except:
            print('%s: *** Unexpected error exporting MS to ASDM, regression failed ***' % myname)
            raise

        try:
            print("Reimporting the created ASDM (v3)....")
            fromsdm = sdm(asdmname)
            fromsdm.toms(vis=reimp_msname, wvr_corrected_data='no', useversion='v3')
            print("Testing existence of reimported MS ....")
            if(not os.path.exists(reimp_msname)):
                print("MS %s doesn't exist." % reimp_msname)
                raise Exception
            print("Testing equivalence of the original and the reimported MS.")
            tblocal.open(myms_dataset_name)
            nrowsorig = tblocal.nrows()
            print("Original MS contains %s integrations." % nrowsorig)
            tblocal.close()
            tblocal.open(reimp_msname)
            nrowsreimp = tblocal.nrows()
            tblocal.close()
            print("Reimported MS contains %s integrations." % nrowsreimp)
            if(not nrowsreimp==nrowsorig):
                print("Numbers of integrations disagree.")
                raise Exception
        except:
            print('%s: *** Unexpected error reimporting the exported ASDM, regression failed ***' % myname)
            raise

class asdm_import3(test_base):

    def setUp(self):
        self.setUp_xosro()

    def tearDown(self):
        os.system('rm -rf '+self.asdm)

    def test_CAS4532(self):
        '''importasdm CAS-4532: white spaces on Antenna.xml'''
        # The X_osro_scan1/Antenna.xml and SpectralWindow.xml
        # contain white spaces between some of the contents and
        # the tags. This should not cause any error in the XML
        # parser from readXML

        flagcmddict = readXML(self.asdm, 0.0)
        self.assertTrue(flagcmddict, 'Some XML file may contain white spaces not handled by readXML')

        self.assertEqual(flagcmddict.keys().__len__(),214)


class asdm_import5(test_base):

    def setUp(self):
        self.setUp_m51()

    def tearDown(self):
        shutil.rmtree(myasdm_dataset_name)
        shutil.rmtree(myms_dataset_name)
        shutil.rmtree(msname,ignore_errors=True)
        shutil.rmtree(msname+'.flagversions',ignore_errors=True)
        os.system('rm -rf reference.ms* reimported-M51.ms*')


    def test1_lazy1(self):
        '''Asdm-import: Test good v1.2 input with default filler in lazy mode, with_pointing_correction=True'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }

        mysdm = sdm(myasdm_dataset_name)
        self.res = mysdm.toms(lazy=True, with_pointing_correction=True)
        self.assertTrue(self.res)
        print("%s: Success! Now checking output ..." % myname)
        mscomponents = set(["table.f17asdmindex",
                            "ANTENNA/table.dat",
                            "DATA_DESCRIPTION/table.dat",
                            "FEED/table.dat",
                            "FIELD/table.dat",
                            "FLAG_CMD/table.dat",
                            "HISTORY/table.dat",
                            "OBSERVATION/table.dat",
                            "POINTING/table.dat",
                            "POLARIZATION/table.dat",
                            "PROCESSOR/table.dat",
                            "SOURCE/table.dat",
                            "SPECTRAL_WINDOW/table.dat",
                            "STATE/table.dat",
                            "SYSCAL/table.dat",
                            "ANTENNA/table.f0",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.f0",
                            "FIELD/table.f0",
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.f0",
                            "POINTING/table.f0",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.f0",
                            "SOURCE/table.f0",
                            "SPECTRAL_WINDOW/table.f0",
                            "STATE/table.f0",
                            "SYSCAL/table.f0"
                            ])
        for name in mscomponents:
            if not os.access(msname+"/"+name, os.F_OK):
                print('%s: Error  %s/%s doesn\'t exist ...' % (myname,msname,name))
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+msname+'/'+name+' does not exist'
            else:
                print("%s: %s present." % (myname,name))
        print("%s: MS exists. All tables present. Try opening as MS ..." % myname)
        try:
            mslocal.open(msname)
        except:
            print("%s: Error  Cannot open MS table %s" % (myname,msname))
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+msname
        else:
            mslocal.close()
            print("%s: OK. Checking tables in detail ..." % myname)

            fromsdm = sdm(myasdm_dataset_name)
            fromsdm.toms(vis='reference.ms', lazy=False, overwrite=True)

            if(os.path.exists('reference.ms')):
                retValue['success'] = th.compTables('reference.ms', msname, ['FLAG', 'FLAG_CATEGORY', 'DATA','WEIGHT_SPECTRUM'],
                                                    0.001)
                retValue['success'] = th.compVarColTables('reference.ms', msname, 'DATA',1E-5)
                retValue['success'] = th.compVarColTables('reference.ms', msname, 'FLAG')

                for subtname in ["ANTENNA",
                                 "DATA_DESCRIPTION",
                                 "FEED",
                                 "FIELD",
                                 "FLAG_CMD",
                                 "OBSERVATION",
                                 #"POINTING", # expect difference since with_pointing_correction=True
                                 "POLARIZATION",
                                 "PROCESSOR",
                                 "SOURCE",
                                 "SPECTRAL_WINDOW",
                                 "STATE",
                                 "SYSCAL"]:

                    print("\n*** Subtable %s" % subtname)
                    excllist = []
                    if subtname=='SOURCE':
                        excllist=['POSITION', 'TRANSITION', 'REST_FREQUENCY', 'SYSVEL']
                    if subtname=='SPECTRAL_WINDOW':
                        excllist=['CHAN_FREQ', 'CHAN_WIDTH', 'EFFECTIVE_BW', 'RESOLUTION']
                        for colname in excllist:
                            retValue['success'] = th.compVarColTables('reference.ms/SPECTRAL_WINDOW',
                                                                      msname+'/SPECTRAL_WINDOW', colname, 0.01) and retValue['success']

                    try:
                        retValue['success'] = th.compTables('reference.ms/'+subtname,
                                                            msname+'/'+subtname, excllist,
                                                            0.01) and retValue['success']
                    except:
                        retValue['success'] = False
                        print("ERROR for table %s" % subtname)

                print("\n*** Subtable POINTING")
                try:
                    retValue['success'] = retValue['success'] and not (th.compTables('reference.ms/POINTING', # expect difference since with_pointing_correction=True
                                                                                     msname+'/POINTING', [], 0.0))
                except:
                    retValue['success'] = False
                    print("ERROR: POINTING tables should differ in this test.")




        self.assertTrue(retValue['success'],retValue['error_msgs'])


class asdm_import6(test_base):

    def setUp(self):
        self.setUp_acaex()

    def tearDown(self):
        myasdmname = 'uid___A002_X72bc38_X000'
        os.system('rm '+myasdmname) # a link
        shutil.rmtree(msname,ignore_errors=True)
        shutil.rmtree(msname+'.flagversions',ignore_errors=True)
        shutil.rmtree(myasdmname+'.ms',ignore_errors=True)
        shutil.rmtree(myasdmname+'.ms.flagversions',ignore_errors=True)
        os.system('rm -rf reference.ms*')

    def test6_lazy1(self):
        '''Asdm-import: Test good ACA ASDM with mixed pol/channelisation input with default filler in lazy mode'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }

        myasdmname = 'uid___A002_X72bc38_X000'
        themsname = myasdmname + ".ms"

        mysdm = sdm(myasdmname)
        self.res = mysdm.toms(vis=themsname, lazy=True, scans='0:1~3') # only the first 3 scans to save time
        self.assertTrue(self.res)

        #test that scratch columns can be created from lazy import
        cblocal = calibrater( )
        cblocal.open(themsname)
        cblocal.close( )
        tblocal = table( )
        tblocal.open(themsname)
        self.assertIn('CORRECTED_DATA', tblocal.getdesc())
        self.assertIn('MODEL_DATA', tblocal.getdesc())
        tblocal.close()

        print("%s: Success! Now checking output ..." % myname)
        mscomponents = set(["ANTENNA/table.dat",
                            "DATA_DESCRIPTION/table.dat",
                            "FEED/table.dat",
                            "FIELD/table.dat",
                            "FLAG_CMD/table.dat",
                            "HISTORY/table.dat",
                            "OBSERVATION/table.dat",
                            "POINTING/table.dat",
                            "POLARIZATION/table.dat",
                            "PROCESSOR/table.dat",
                            "SOURCE/table.dat",
                            "SPECTRAL_WINDOW/table.dat",
                            "STATE/table.dat",
                            "SYSCAL/table.dat",
                            "ANTENNA/table.f0",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.f0",
                            "FIELD/table.f0",
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.f0",
                            "POINTING/table.f0",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.f0",
                            "SOURCE/table.f0",
                            "SPECTRAL_WINDOW/table.f0",
                            "STATE/table.f0",
                            "SYSCAL/table.f0"
                            ])
        for name in mscomponents:
            if not os.access(themsname+"/"+name, os.F_OK):
                print("%s: Error  %s/%s doesn't exist ..." % (myname,themsname,name))
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+themsname+'/'+name+' does not exist'
            else:
                print("%s: %s present." % (myname,name))
        print("%s: MS exists. All tables present. Try opening as MS ..." % myname)
        try:
            mslocal.open(themsname)
        except:
            print("%s: Error  Cannot open MS table %s" % (myname,themsname))
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+themsname
        else:
            mslocal.close()
            print("%s: OK. Checking tables in detail ..." % myname)

            fromsdm = sdm(myasdmname)
            fromsdm.toms(vis='reference.ms', lazy=False, overwrite=True, scans='0:1~3')

            if(os.path.exists('reference.ms')):
                retValue['success'] = th.checkwithtaql("select from [select from reference.ms orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "
                                                    +themsname+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(near(t1.DATA,t2.DATA, 1.e-06)))") == 0
                if not retValue['success']:
                    print("ERROR: DATA does not agree with reference.")
                else:
                    print("DATA columns agree.")
                retValueTmp = th.checkwithtaql("select from [select from reference.ms orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "
                                            +themsname+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(t1.FLAG==t2.FLAG)) ") == 0
                if not retValueTmp:
                    print("ERROR: FLAG does not agree with reference.")
                else:
                    print("FLAG columns agree.")

                retValue['success'] = retValue['success'] and retValueTmp

                for subtname in ["ANTENNA",
                                 "DATA_DESCRIPTION",
                                 "FEED",
                                 "FIELD",
                                 "FLAG_CMD",
                                 "OBSERVATION",
                                 "POINTING",
                                 "POLARIZATION",
                                 "PROCESSOR",
                                 "SOURCE",
                                 "SPECTRAL_WINDOW",
                                 "STATE",
                                 "SYSCAL"]:

                    print("\n*** Subtable %s" % subtname)
                    excllist = []
                    if subtname=='SOURCE':
                        excllist=['POSITION', 'TRANSITION', 'REST_FREQUENCY', 'SYSVEL']
                    if subtname=='SYSCAL':
                        excllist=['TANT_SPECTRUM', 'TANT_TSYS_SPECTRUM']
                    if subtname=='SPECTRAL_WINDOW':
                        excllist=['CHAN_FREQ', 'CHAN_WIDTH', 'EFFECTIVE_BW', 'RESOLUTION']
                        for colname in excllist:
                            retValue['success'] = th.compVarColTables('reference.ms/SPECTRAL_WINDOW',
                                                                      themsname+'/SPECTRAL_WINDOW', colname, 0.01) and retValue['success']

                    try:
                        retValue['success'] = th.compTables('reference.ms/'+subtname,
                                                            themsname+'/'+subtname, excllist,
                                                            0.01) and retValue['success']
                    except:
                        retValue['success'] = False
                        print("ERROR for table %s" % subtname)


        self.assertTrue(retValue['success'],retValue['error_msgs'])


class asdm_import7(test_base):

    def setUp(self):
        self.setUp_12mex()
        self.setUp_eph()
        self.setUp_SD()

    def tearDown(self):
        for myasdmname in ['uid___A002_X71e4ae_X317_short', 'uid___A002_X997a62_X8c-short', 'uid___A002_X6218fb_X264']:
            os.system('rm -f '+myasdmname) # a link
            shutil.rmtree(myasdmname+".ms",ignore_errors=True)
            shutil.rmtree(myasdmname+'.ms.flagversions',ignore_errors=True)
        shutil.rmtree("reference.ms",ignore_errors=True)
        shutil.rmtree("reference.ms.flagversions",ignore_errors=True)

    def test7_lazy1(self):
        '''Asdm-import: Test good 12 m ASDM with mixed pol/channelisation input with default filler in lazy mode'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }

        myasdmname = 'uid___A002_X71e4ae_X317_short'
        themsname = myasdmname+".ms"

        mysdm = sdm(myasdmname)
        self.res = mysdm.toms(vis=themsname, lazy=True, scans='0:1~4') # only the first 4 scans to save time
        self.assertTrue(self.res)
        print("%s: Success! Now checking output ..." % myname)
        mscomponents = set(["ANTENNA/table.dat",
                            "CALDEVICE/table.dat",
                            "DATA_DESCRIPTION/table.dat",
                            "FEED/table.dat",
                            "FIELD/table.dat",
                            "FLAG_CMD/table.dat",
                            "HISTORY/table.dat",
                            "OBSERVATION/table.dat",
                            "POINTING/table.dat",
                            "POLARIZATION/table.dat",
                            "PROCESSOR/table.dat",
                            "SOURCE/table.dat",
                            "SPECTRAL_WINDOW/table.dat",
                            "STATE/table.dat",
                            "SYSCAL/table.dat",
                            "SYSPOWER/table.dat",
                            "WEATHER/table.dat",
                            "ANTENNA/table.f0",
                            "CALDEVICE/table.f0",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.f0",
                            "FIELD/table.f0",
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.f0",
                            "POINTING/table.f0",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.f0",
                            "SOURCE/table.f0",
                            "SPECTRAL_WINDOW/table.f0",
                            "STATE/table.f0",
                            "SYSCAL/table.f0",
                            "SYSPOWER/table.f0",
                            "WEATHER/table.f0"
                            ])
        for name in mscomponents:
            if not os.access(themsname+"/"+name, os.F_OK):
                print("%s: Error  %s/%s doesn't exist ..." % (myname,themsname,name))
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+themsname+'/'+name+' does not exist'
            else:
                print("%s: %s present." % (myname,name))
        print("%s: MS exists. All tables present. Try opening as MS ..." % myname)
        try:
            mslocal.open(themsname)
            mslocal.close()
            print("%s: MS can be opened. Now testing the changing of the asdmref ..." % myname)
            mslocal.open(themsname)
            mslocal.asdmref("./moved_"+myasdmname)
            mslocal.close()
            os.system("mv "+myasdmname+" moved_"+myasdmname)

            mslocal.open(themsname)

        except:
            print("%s: Error  Cannot open MS table %s" % (myname,themsname))
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+themsname
        else:
            mslocal.close()
            print("%s: OK. Checking tables in detail ..." % myname)
            fromsdm = sdm("moved_"+myasdmname)
            fromsdm.toms(vis='reference.ms', lazy=False, overwrite=True, scans='0:1~4')
            if(os.path.exists('reference.ms')):
                retValue['success'] = th.checkwithtaql("select from [select from reference.ms orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "
                                                    +themsname+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(near(t1.DATA,t2.DATA, 1.e-06)))") == 0
                if not retValue['success']:
                    print("ERROR: DATA does not agree with reference.")
                else:
                    print("DATA columns agree.")

                retValueTmp = th.checkwithtaql("select from [select from reference.ms orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "
                                                    +themsname+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(near(t1.WEIGHT,t2.WEIGHT, 1.e-06)))") == 0
                if not retValueTmp:
                    print("ERROR: WEIGHT does not agree with reference.")
                else:
                    print("WEIGHT columns agree.")

                retValueTmp2 = th.checkwithtaql("select from [select from reference.ms orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "
                                            +themsname+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(t1.FLAG==t2.FLAG)) ") == 0
                if not retValueTmp2:
                    print("ERROR: FLAG does not agree with reference.")
                else:
                    print("FLAG columns agree.")

                retValue['success'] = retValue['success'] and retValueTmp and retValueTmp2

                for subtname in ["ANTENNA",
                                 "CALDEVICE",
                                 "DATA_DESCRIPTION",
                                 "FEED",
                                 "FIELD",
                                 "FLAG_CMD",
                                 "OBSERVATION",
                                 "POLARIZATION",
                                 "PROCESSOR",
                                 "SOURCE",
                                 "SPECTRAL_WINDOW",
                                 "STATE",
                                 "SYSCAL",
                                 "SYSPOWER",
                                 "WEATHER"]:

                    print("\n*** Subtable %s" % subtname)
                    excllist = []
                    if subtname=='CALDEVICE':
                        excllist=['NOISE_CAL','CAL_EFF']
                    if subtname=='SOURCE':
                        excllist=['POSITION', 'TRANSITION', 'REST_FREQUENCY', 'SYSVEL']
                    if subtname=='SYSCAL':
                        excllist=['TANT_SPECTRUM', 'TANT_TSYS_SPECTRUM']
                    if subtname=='SPECTRAL_WINDOW':
                        excllist=['CHAN_FREQ', 'CHAN_WIDTH', 'EFFECTIVE_BW', 'RESOLUTION', 'ASSOC_SPW_ID', 'ASSOC_NATURE']
                        for colname in excllist:
                            if colname!='ASSOC_NATURE':
                                retValue['success'] = th.compVarColTables('reference.ms/SPECTRAL_WINDOW',
                                                                          themsname+'/SPECTRAL_WINDOW', colname, 0.01) and retValue['success']
                    if subtname=='POLARIZATION':
                        excllist=['CORR_TYPE', 'CORR_PRODUCT']
                        for colname in excllist:
                            retValue['success'] = th.compVarColTables('reference.ms/POLARIZATION',
                                                                      themsname+'/POLARIZATION', colname, 0.01) and retValue['success']
                    try:
                        retValue['success'] = th.compTables('reference.ms/'+subtname,
                                                            themsname+'/'+subtname, excllist,
                                                            0.01) and retValue['success']
                    except:
                        retValue['success'] = False
                        print("ERROR for table %s" % subtname)

                try:
                    # test that the PRESSURE column has the expected units
                    wxcalOK = tblocal.open(themsname+'/WEATHER')
                    if wxcalOK:
                        wxcalOK = tblocal.getcolkeyword("PRESSURE","QuantumUnits")[0] == 'hPa'
                        tblocal.close()
                    retValue['success'] = wxcalOK and retValue['success']
                    if not wxcalOK:
                        print("PRESSURE column in WEATHER table is missing or has incorrect units")
                except:
                    retValue['success'] = False
                    print("ERROR getting units of PRESSURE column in WEATHER table.")

                try:
                    # test that the SDM_WINDOW_FUNCTION column exists and has the exepcted values
                    winFuncOK = tblocal.open(themsname+'/SPECTRAL_WINDOW')
                    if winFuncOK:
                        winFuncCol = tblocal.getcol('SDM_WINDOW_FUNCTION')
                        tblocal.close()
                        # test values here
                        # expect 55 rows, rows 1:24 are HANNING, the rest are UNIFORM
                        indx = numpy.arange(len(winFuncCol))
                        winFuncOK = winFuncOK and (numpy.array_equal(indx[winFuncCol=="HANNING"],(numpy.arange(24)+1)))
                        winFuncOK = winFuncOK and (len(indx[winFuncCol=="UNIFORM"])==31)

                    retValue['success'] = winFuncOK and retValue['success']
                    if not winFuncOK:
                        print("SDM_WINDOW_FUNCTION column in the SPECTRAL_WINDOW table is missing or has incorrect values")
                except:
                    retValue['success'] = False
                    print("ERROR checking the value of the SDM_WINDOW_FUNCTION column in the SPECTRAL_WINDOW table.")

        os.system("mv moved_"+myasdmname+" "+myasdmname)

        self.assertTrue(retValue['success'],retValue['error_msgs'])

    def test7_lazy2(self):
        '''Asdm-import: Test good 12 m ASDM with mixed pol/channelisation input with default filler in lazy mode with reading the BDF flags'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }

        myasdmname = 'uid___A002_X71e4ae_X317_short'
        themsname = myasdmname+".ms"
        mysdm = sdm(myasdmname)
        self.res = mysdm.toms(vis=themsname, lazy=True, bdfflags=True)
        self.assertTrue(self.res)
        print("%s: Success! Now checking output ..." % myname)
        mscomponents = set(["ANTENNA/table.dat",
                            "DATA_DESCRIPTION/table.dat",
                            "FEED/table.dat",
                            "FIELD/table.dat",
                            "FLAG_CMD/table.dat",
                            "HISTORY/table.dat",
                            "OBSERVATION/table.dat",
                            "POINTING/table.dat",
                            "POLARIZATION/table.dat",
                            "PROCESSOR/table.dat",
                            "SOURCE/table.dat",
                            "SPECTRAL_WINDOW/table.dat",
                            "STATE/table.dat",
                            "SYSCAL/table.dat",
                            "ANTENNA/table.f0",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.f0",
                            "FIELD/table.f0",
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.f0",
                            "POINTING/table.f0",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.f0",
                            "SOURCE/table.f0",
                            "SPECTRAL_WINDOW/table.f0",
                            "STATE/table.f0",
                            "SYSCAL/table.f0"
                            ])
        for name in mscomponents:
            if not os.access(themsname+"/"+name, os.F_OK):
                print("%s: Error  %s/%s doesn't exist ..." % (myname,themsname,name))
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+themsname+'/'+name+' does not exist'
            else:
                print("%s: %s present." % (myname,name))
        print("%s: MS exists. All tables present. Try opening as MS ..." % myname)
        try:
            mslocal.open(themsname)
        except:
            print("%s: Error  Cannot open MS table %s" % (myname,themsname))
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+themsname
        else:
            mslocal.close()
            print("%s: OK. Checking tables in detail ..." % myname)
            fromsdm = sdm(myasdmname)
            fromsdm.toms(vis='reference.ms', lazy=True, overwrite=True, bdfflags=False)
            if(os.path.exists('reference.ms')):
                retValue['success'] = th.checkwithtaql("select from [select from reference.ms orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "
                                                    +themsname+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(near(t1.DATA,t2.DATA, 1.e-06)))") == 0
                if not retValue['success']:
                    print("ERROR: DATA does not agree with reference.")
                else:
                    print("DATA columns agree.")
                retValueTmp = th.checkwithtaql("select from [select from reference.ms orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "
                                            +themsname+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(t1.FLAG==t2.FLAG)) ") != 0
                if not retValueTmp:
                    print("ERROR: FLAG columns do agree with reference but they shouldn't.")
                else:
                    print("FLAG columns do not agree as expected.")

                retValue['success'] = retValue['success'] and retValueTmp

                for subtname in ["ANTENNA",
                                 "DATA_DESCRIPTION",
                                 "FEED",
                                 "FIELD",
                                 "FLAG_CMD",
                                 "OBSERVATION",
                                 "POLARIZATION",
                                 "PROCESSOR",
                                 "SOURCE",
                                 "SPECTRAL_WINDOW",
                                 "STATE",
                                 "SYSCAL"]:

                    print("\n*** Subtable %s" % subtname)
                    excllist = []
                    if subtname=='SOURCE':
                        excllist=['POSITION', 'TRANSITION', 'REST_FREQUENCY', 'SYSVEL']
                    if subtname=='SYSCAL':
                        excllist=['TANT_SPECTRUM', 'TANT_TSYS_SPECTRUM']
                    if subtname=='SPECTRAL_WINDOW':
                        excllist=['CHAN_FREQ', 'CHAN_WIDTH', 'EFFECTIVE_BW', 'RESOLUTION', 'ASSOC_SPW_ID', 'ASSOC_NATURE']
                        for colname in excllist:
                            if colname!='ASSOC_NATURE':
                                retValue['success'] = th.compVarColTables('reference.ms/SPECTRAL_WINDOW',
                                                                          themsname+'/SPECTRAL_WINDOW', colname, 0.01) and retValue['success']
                    if subtname=='POLARIZATION':
                        excllist=['CORR_TYPE', 'CORR_PRODUCT']
                        for colname in excllist:
                            retValue['success'] = th.compVarColTables('reference.ms/POLARIZATION',
                                                                      themsname+'/POLARIZATION', colname, 0.01) and retValue['success']
                    try:
                        retValue['success'] = th.compTables('reference.ms/'+subtname,
                                                            themsname+'/'+subtname, excllist,
                                                            0.01) and retValue['success']
                    except:
                        retValue['success'] = False
                        print("ERROR for table %s" % subtname)


        self.assertTrue(retValue['success'],retValue['error_msgs'])

    @unittest.skip("uses --inerpolate-ephemeris option in asdm2MS not (yet?) available in sdm.toms")
    def test7_lazy3(self):
        '''Asdm-import: Test good 12 m ASDM with Ephemeris table in lazy mode'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }

        myasdmname = 'uid___A002_X997a62_X8c-short'
        themsname = myasdmname+".ms"
        mysdm = sdm(myasdmname)
        self.res = mysdm.toms( vis=themsname, lazy=True, convert_ephem2geo=True,
                               process_pointing=False, flagbackup=False )
        self.assertEqual(self.res, None)
        print("%s: Success! Now checking output ..." % myname)
        mscomponents = set(["FIELD/table.dat",
                            "FIELD/EPHEM0_Mars_57034.9.tab",
                            "FIELD/EPHEM1_Titania_57034.9.tab"
                            ])
        for name in mscomponents:
            if not os.access(themsname+"/"+name, os.F_OK):
                print("%s: Error  %s/%s doesn't exist ..." % (myname,themsname,name))
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+themsname+'/'+name+' does not exist'
            else:
                print("%s: %s present." % (myname,name))
        print("%s: MS exists. All relevant tables present. Try opening as MS ..." % myname)
        try:
            mslocal.open(themsname)
        except:
            print("%s: Error  Cannot open MS table %" % (myname,themsname))
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+themsname


        print("%s:  testing FIELD values in ms.summary()" % myname)
        try:
            mssum = mslocal.summary()
            # only Mars appears here because this short SDM only contains a single scan and that uses Mars
            self.assertEqual(mssum['scan_1']['0']['FieldName'],'Mars')
            self.assertAlmostEqual(mssum['field_0']['direction']['m0']['value'],-0.4770797859505159,15)
            self.assertAlmostEqual(mssum['field_0']['direction']['m1']['value'],-0.2154815444753364,15)
        except:
            print("%s: Error ms summary has an unexpect source or direction value" % myname)
            retValue['success']=False
            retValue['error_msg']=retValue['err_msg']+'Unexpected source or direction value in ms summary '+thismsname + '\n'


        mslocal.close()

        ephems = []
        # values from indivual rows were picked for no particular reason except verify they've not changed
        ephems.append({'name':"FIELD/EPHEM0_Mars_57034.9.tab",
                       'nrows':27,
                       'rows':[{'row':10,
                                'values':{'MJD':57035.041666666664,
                                          'RA':332.7140437500001,
                                          'DEC':-12.327346944444447,
                                          'Rho':2.024609480125507,
                                          'RadVel':723729.77502873}},
                               {'row':22,
                                'values':{'MJD':57035.208333333336,
                                          'RA':332.8387870833333,
                                          'DEC':-12.2793975,
                                          'Rho':2.0254053468626436,
                                          'RadVel':705588.202526264}}
                               ]
                       }
                      )

        ephems.append({'name':"FIELD/EPHEM1_Titania_57034.9.tab",
                       'nrows':45,
                       'rows':[{'row':17,
                                'values':{'MJD':57035.055555555555,
                                          'RA':11.813166666666666,
                                          'DEC':4.365749999999999,
                                          'Rho':20.150883673698488,
                                          'RadVel':2730048.0839084117}},
                               {'row':40,
                                'values':{'MJD':57035.21527777778,
                                          'RA':11.816041666666667,
                                          'DEC':4.3661111111111115,
                                          'Rho':20.153736461701364,
                                          'RadVel':2711142.1699538543}}
                               ]
                       }
                      )

        for ephem in ephems:
            print("%s: Testing various things in ephemeris %s ..." % (myname,ephem['name']))

            tblocal.open(themsname+"/"+ephem['name'])
            kw = tblocal.getkeywords()
            nrows = tblocal.nrows()
            if not nrows==ephem['nrows']:
                print("%s: Error. unexpected number of rows in ephemeris : %s" % (myname,ephem['name']))
                retValue['success']=False
                retValue['error_msg']=retValue['error_msgs']+' Unexpected number of rows in ephemeris table :'+ ephem['name'] + '\n'

            for row in ephem['rows']:
                thisRow = row['row']
                for colname in row['values']:
                    thisVal = tblocal.getcell(colname,thisRow)
                    self.assertAlmostEqual(thisVal,row['values'][colname],10)

            # unfilled columns
            self.assertEqual((tblocal.getcol('diskLong') != 0.0).sum(),0)
            self.assertEqual((tblocal.getcol('diskLat') != 0.0).sum(),0)

            tblocal.close()
            geodist = kw['GeoDist'] # (km above reference ellipsoid)
            geolat = kw['GeoLat'] # (deg)
            geolong = kw['GeoLong'] # (deg)
            if not (geodist==geolat==geolong==0.):
                print("%s: ERROR." % myname)
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+' Ephemeris was not converted to GEO for '+themsname+'\n'
            prsys = kw['posrefsys']
            if not (prsys=="ICRF/ICRS"):
                print("%s: ERROR." % myname)
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+' posrefsys keyword is not ICRF/ICRS '+themsname+'\n'

        # fill and request an interpolated table.  Tests asdm2MS directly as this option isn't
        # available in importasdm

        print("%s filling an interpolated version of the same ephemeris" % myname)
        themsname_interp = myasdmname+".interp.ms"
        execute_string = "asdm2MS --no-pointing --interpolate-ephemeris 'yes' " + myasdmname + ' ' + themsname_interp
        print('%s executing : %s' % (myname,execute_string))
        exitcode = os.system(execute_string)
        self.assertEqual(exitcode,0)
        ce.convert2geo(themsname_interp, '*') # convert the ephemeris to GEO
        # note that the recalculation of UVW and the adjustment of the SOURCE table are not
        # done here the way they would be done if filled via importasdm
        print("%s: Success! Now checking output ..." % myname)
        for name in mscomponents:
            if not os.access(themsname_interp+"/"+name, os.F_OK):
                print("%s: Error  %s/%s doesn't exist ..." % (myname,themsname_interp,name))
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+themsname_interp+'/'+name+' does not exist'
            else:
                print("%s: %s present." % (myname,name))
        print("%s: MS exists. All relevant tables present. Try opening as MS ..." % myname)
        try:
            mslocal.open(themsname_interp)
        except:
            print("%s: Error  Cannot open MS table %s" % (myname,themsname_interp))
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+themsname_interp
        print("%s:  testing FIELD values in ms.summary()" % myname)
        try:
            mssum = mslocal.summary()
            # only Mars appears here because this short SDM only contains a single scan and that uses Mars
            self.assertEqual(mssum['scan_1']['0']['FieldName'],'Mars')
            # difference here is < 0".0004 of the above, non-interpolated value
            self.assertAlmostEqual(mssum['field_0']['direction']['m0']['value'],-0.4770797877079177,15)
            # difference here is < 0".00005 of the above, non-interpolated value
            self.assertAlmostEqual(mssum['field_0']['direction']['m1']['value'],-0.2154815442529733,15)
        except:
            print("%s: Error ms summary has an unexpect source or direction value" % myname)
            retValue['success']=False
            retValue['error_msg']=retValue['err_msg']+'Unexpected source or direction value in ms summary '+thismsname + '\n'

        mslocal.close()
        ephems = []
        # values from indivual rows were picked for no particular reason except verify they've not changed
        # these rows
        ephems.append({'name':"FIELD/EPHEM0_Mars_57034.9.tab",
                       'nrows':361,
                       'rows':[{'row':105,
                                'values':{'MJD':57035.008000000001630,
                                          'RA':332.688983703339886,
                                          'DEC':-12.337033046664128,
                                          'Rho':2.024447666954067,
                                          'RadVel':722270.482337458524853}},
                               {'row':320,
                                'values':{'MJD':57035.222999999998137,
                                          'RA':332.849807533339913,
                                          'DEC':-12.275182948886352,
                                          'Rho':2.025474163348287,
                                          'RadVel':702390.653405150165781}}
                               ]
                       }
                      )

        ephems.append({'name':"FIELD/EPHEM1_Titania_57034.9.tab",
                       'nrows':306,
                       'rows':[{'row':95,
                                'values':{'MJD':57035.033000000003085,
                                          'RA':11.812802333333494,
                                          'DEC':4.365715333333369,
                                          'Rho':20.150479182212013,
                                          'RadVel':2725243.279430569149554}},
                               {'row':250,
                                'values':{'MJD':57035.188000000001921,
                                          'RA':11.815509000000159,
                                          'DEC':4.366057555555590,
                                          'Rho':20.153251583732832,
                                          'RadVel':2721431.250284913461655}}
                               ]
                       }
                      )

        for ephem in ephems:
            print("%s: Testing various things in ephemeris %s ..." % (myname,ephem['name']))

            tblocal.open(themsname_interp+"/"+ephem['name'])
            kw = tblocal.getkeywords()
            nrows = tblocal.nrows()
            if not nrows==ephem['nrows']:
                print("%s: Error. unexpected number of rows in ephemeris : %s" % (myname,ephem['name']))
                retValue['success']=False
                retValue['error_msg']=retValue['error_msgs']+' Unexpected number of rows in ephemeris table :'+ ephem['name'] + '\n'

            for row in ephem['rows']:
                thisRow = row['row']
                for colname in row['values']:
                    thisVal = tblocal.getcell(colname,thisRow)
                    self.assertAlmostEqual(thisVal,row['values'][colname],10)

            # unfilled columns
            self.assertEqual((tblocal.getcol('diskLong') != 0.0).sum(),0)
            self.assertEqual((tblocal.getcol('diskLat') != 0.0).sum(),0)

            tblocal.close()
            geodist = kw['GeoDist'] # (km above reference ellipsoid)
            geolat = kw['GeoLat'] # (deg)
            geolong = kw['GeoLong'] # (deg)
            if not (geodist==geolat==geolong==0.):
                print("%s: ERROR." % myname)
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+' Ephemeris was not converted to GEO for '+themsname_interp+'\n'
            prsys = kw['posrefsys']
            if not (prsys=="ICRF/ICRS"):
                print("%s: ERROR." % myname)
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+' posrefsys keyword is not ICRF/ICRS '+themsname_interp+'\n'

        self.assertTrue(retValue['success'],retValue['error_msgs'])
        print("%s: OK." % myname)


    def test7_lazy4(self):
        '''Asdm-import: Test good 12 m ASDM with mixed pol/channelisation input with default filler in lazy mode selecting only AUTO data, writing to FLOAT_DATA'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }

        myasdmname = 'uid___A002_X71e4ae_X317_short'
        themsname = myasdmname+".ms"
        mysdm = sdm(myasdmname)
        self.res = mysdm.toms(vis=themsname, ocorr_mode="ao", lazy=True, scans='0:1~4') # only the first 4 scans to save time
        self.assertTrue(self.res)
        print("%s: Success! Now checking output ..." % myname)
        mscomponents = set(["ANTENNA/table.dat",
                            "CALDEVICE/table.dat",
                            "DATA_DESCRIPTION/table.dat",
                            "FEED/table.dat",
                            "FIELD/table.dat",
                            "FLAG_CMD/table.dat",
                            "HISTORY/table.dat",
                            "OBSERVATION/table.dat",
                            "POINTING/table.dat",
                            "POLARIZATION/table.dat",
                            "PROCESSOR/table.dat",
                            "SOURCE/table.dat",
                            "SPECTRAL_WINDOW/table.dat",
                            "STATE/table.dat",
                            "SYSCAL/table.dat",
                            "SYSPOWER/table.dat",
                            "ANTENNA/table.f0",
                            "CALDEVICE/table.f0",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.f0",
                            "FIELD/table.f0",
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.f0",
                            "POINTING/table.f0",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.f0",
                            "SOURCE/table.f0",
                            "SPECTRAL_WINDOW/table.f0",
                            "STATE/table.f0",
                            "SYSCAL/table.f0",
                            "SYSPOWER/table.f0"
                            ])
        for name in mscomponents:
            if not os.access(themsname+"/"+name, os.F_OK):
                print("%s: Error  %s/%s doesn't exist ..." % (myname,themsname,name))
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+themsname+'/'+name+' does not exist'
            else:
                print("%s: %s present." % (myname,name))
        print("%s: MS exists. All tables present. Try opening as MS ..." % myname)
        try:
            mslocal.open(themsname)
            mslocal.close()
            print("%s: MS can be opened. Now testing the changing of the asdmref ..." % myname)
            mslocal.open(themsname)
            mslocal.asdmref("./moved_"+myasdmname)
            mslocal.close()
            os.system("mv "+myasdmname+" moved_"+myasdmname)

            mslocal.open(themsname)

        except:
            print("%s: Error  Cannot open MS table %s" % (myname,themsname))
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+themsname
        else:
            mslocal.close()
            print("%s: OK. Checking tables in detail ..." % myname)
            fromsdm = sdm("moved_"+myasdmname)
            fromsdm.toms(vis='reference.ms', ocorr_mode="ao", lazy=False, overwrite=True, scans='0:1~3')

            if(os.path.exists('reference.ms')):
                retValue['success'] = th.checkwithtaql("select from [select from reference.ms orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "
                                                    +themsname+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(near(t1.FLOAT_DATA,t2.FLOAT_DATA, 1.e-06)))") == 0
                if not retValue['success']:
                    print("ERROR: DATA does not agree with reference.")
                else:
                    print("DATA columns agree.")

                retValueTmp = th.checkwithtaql("select from [select from reference.ms orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "
                                                    +themsname+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(near(t1.WEIGHT,t2.WEIGHT, 1.e-06)))") == 0
                if not retValueTmp:
                    print("ERROR: WEIGHT does not agree with reference.")
                else:
                    print("WEIGHT columns agree.")

                retValueTmp2 = th.checkwithtaql("select from [select from reference.ms orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "
                                            +themsname+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(t1.FLAG==t2.FLAG)) ") == 0
                if not retValueTmp2:
                    print("ERROR: FLAG does not agree with reference.")
                else:
                    print("FLAG columns agree.")

                retValue['success'] = retValue['success'] and retValueTmp and retValueTmp2

                for subtname in ["ANTENNA",
                                 "DATA_DESCRIPTION",
                                 "FEED",
                                 "FIELD",
                                 "FLAG_CMD",
                                 "OBSERVATION",
                                 "POLARIZATION",
                                 "PROCESSOR",
                                 "SOURCE",
                                 "SPECTRAL_WINDOW",
                                 "STATE",
                                 "SYSCAL"]:

                    print("\n*** Subtable %s" % subtname)
                    excllist = []
                    if subtname=='SOURCE':
                        excllist=['POSITION', 'TRANSITION', 'REST_FREQUENCY', 'SYSVEL']
                    if subtname=='SYSCAL':
                        excllist=['TANT_SPECTRUM', 'TANT_TSYS_SPECTRUM']
                    if subtname=='SPECTRAL_WINDOW':
                        excllist=['CHAN_FREQ', 'CHAN_WIDTH', 'EFFECTIVE_BW', 'RESOLUTION', 'ASSOC_SPW_ID', 'ASSOC_NATURE']
                        for colname in excllist:
                            if colname!='ASSOC_NATURE':
                                retValue['success'] = th.compVarColTables('reference.ms/SPECTRAL_WINDOW',
                                                                          themsname+'/SPECTRAL_WINDOW', colname, 0.01) and retValue['success']
                    if subtname=='POLARIZATION':
                        excllist=['CORR_TYPE', 'CORR_PRODUCT']
                        for colname in excllist:
                            retValue['success'] = th.compVarColTables('reference.ms/POLARIZATION',
                                                                      themsname+'/POLARIZATION', colname, 0.01) and retValue['success']
                    try:
                        retValue['success'] = th.compTables('reference.ms/'+subtname,
                                                            themsname+'/'+subtname, excllist,
                                                            0.01) and retValue['success']
                    except:
                        retValue['success'] = False
                        print("ERROR for table %s" % subtname)

        os.system("mv moved_"+myasdmname+" "+myasdmname)

        self.assertTrue(retValue['success'],retValue['error_msgs'])

    def test7_lazy5(self):
        '''Asdm-import: Test TP asdm with default filler in lazy mode selecting only AUTO data, writing to FLOAT_DATA'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }

        myasdmname = 'uid___A002_X6218fb_X264'
        themsname = myasdmname+".ms"
        mysdm = sdm(myasdmname)
        self.res = mysdm.toms(vis=themsname, ocorr_mode="ao", bdfflags=True, applyflags=True, lazy=True)
        self.assertTrue(self.res)
        print("%s: Success! Now checking output ..." % myname)
        mscomponents = set(["ANTENNA/table.dat",
                            "CALDEVICE/table.dat",
                            "DATA_DESCRIPTION/table.dat",
                            "FEED/table.dat",
                            "FIELD/table.dat",
                            "FLAG_CMD/table.dat",
                            "HISTORY/table.dat",
                            "OBSERVATION/table.dat",
                            "POINTING/table.dat",
                            "POLARIZATION/table.dat",
                            "PROCESSOR/table.dat",
                            "SOURCE/table.dat",
                            "SPECTRAL_WINDOW/table.dat",
                            "STATE/table.dat",
                            "SYSCAL/table.dat",
                            "SYSPOWER/table.dat",
                            "WEATHER/table.dat",
                            "ANTENNA/table.f0",
                            "CALDEVICE/table.f0",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.f0",
                            "FIELD/table.f0",
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.f0",
                            "POINTING/table.f0",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.f0",
                            "SOURCE/table.f0",
                            "SPECTRAL_WINDOW/table.f0",
                            "STATE/table.f0",
                            "SYSCAL/table.f0",
                            "SYSPOWER/table.f0",
                            "WEATHER/table.f0"
                            ])
        for name in mscomponents:
            if not os.access(themsname+"/"+name, os.F_OK):
                print("%s: Error  %s/%s doesn't exist ..." % (myname,themsname,name))
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+themsname+'/'+name+' does not exist'
            else:
                print("%s: %s present." % (myname,name))
        print("%s: MS exists. All tables present. Try opening as MS ..." % myname)
        try:
            mslocal.open(themsname)
            mslocal.close()
            print("%s: MS can be opened. Now testing the changing of the asdmref ..." %  myname)
            mslocal.open(themsname)
            mslocal.asdmref("./moved_"+myasdmname)
            mslocal.close()
            os.system("mv "+myasdmname+" moved_"+myasdmname)

            mslocal.open(themsname)

        except:
            print("%s: Error  Cannot open MS table %s" % (myname,themsname))
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+themsname
        else:
            mslocal.close()
            print("%s: OK. Checking tables in detail ..." % myname)
            fromsdm = sdm("moved_"+myasdmname)
            fromsdm.toms(vis='reference.ms', ocorr_mode="ao", lazy=False, bdfflags=True, applyflags=True, overwrite=True)

            if(os.path.exists('reference.ms')):
                retValue['success'] = th.checkwithtaql("select from [select from reference.ms orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "
                                                    +themsname+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(near(t1.FLOAT_DATA,t2.FLOAT_DATA, 1.e-06)))") == 0
                if not retValue['success']:
                    print("ERROR: DATA does not agree with reference.")
                else:
                    print("DATA columns agree.")

                retValueTmp = th.checkwithtaql("select from [select from reference.ms orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "
                                                    +themsname+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(near(t1.WEIGHT,t2.WEIGHT, 1.e-06)))") == 0
                if not retValueTmp:
                    print("ERROR: WEIGHT does not agree with reference.")
                else:
                    print("WEIGHT columns agree.")

                retValueTmp2 = th.checkwithtaql("select from [select from reference.ms orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "
                                            +themsname+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(t1.FLAG==t2.FLAG)) ") == 0
                if not retValueTmp2:
                    print("ERROR: FLAG does not agree with reference.")
                else:
                    print("FLAG columns agree.")

                retValue['success'] = retValue['success'] and retValueTmp and retValueTmp2

                for subtname in ["ANTENNA",
                                 "DATA_DESCRIPTION",
                                 "FEED",
                                 "FIELD",
                                 "FLAG_CMD",
                                 "OBSERVATION",
                                 "POLARIZATION",
                                 "PROCESSOR",
                                 "SOURCE",
                                 "SPECTRAL_WINDOW",
                                 "STATE",
                                 "SYSCAL",
                                 "WEATHER"]:

                    print("\n*** Subtable %s" % subtname)
                    excllist = []
                    if subtname=='SOURCE':
                        excllist=['POSITION', 'TRANSITION', 'REST_FREQUENCY', 'SYSVEL']
                    if subtname=='SYSCAL':
                        excllist=['TANT_SPECTRUM', 'TANT_TSYS_SPECTRUM']
                    if subtname=='SPECTRAL_WINDOW':
                        excllist=['CHAN_FREQ', 'CHAN_WIDTH', 'EFFECTIVE_BW', 'RESOLUTION', 'ASSOC_SPW_ID', 'ASSOC_NATURE']
                        for colname in excllist:
                            if colname!='ASSOC_NATURE':
                                retValue['success'] = th.compVarColTables('reference.ms/SPECTRAL_WINDOW',
                                                                          themsname+'/SPECTRAL_WINDOW', colname, 0.01) and retValue['success']
                    if subtname=='POLARIZATION':
                        excllist=['CORR_TYPE', 'CORR_PRODUCT']
                        for colname in excllist:
                            retValue['success'] = th.compVarColTables('reference.ms/POLARIZATION',
                                                                      themsname+'/POLARIZATION', colname, 0.01) and retValue['success']
                    try:
                        retValue['success'] = th.compTables('reference.ms/'+subtname,
                                                            themsname+'/'+subtname, excllist,
                                                            0.01) and retValue['success']
                    except:
                        retValue['success'] = False
                        print("ERROR for table %s" % subtname)

        os.system("mv moved_"+myasdmname+" "+myasdmname)

        self.assertTrue(retValue['success'],retValue['error_msgs'])

    @unittest.skip("uses --checkdupints option in asdm2MS not available in sdm.toms")
    def test7_skiprows1(self):
        '''Asdm-import: Test TP asdm, comparing output when duplicate DATA rows are skipped versus not-skipped, lazy and regular, with bdflagging on'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }

        myasdmname = 'uid___A002_X6218fb_X264'
        themsname = myasdmname+".ms"

        # the same tests are done for lazy being False and True
        for lazy in [False, True]:
            print("%s Test loop with lazy = %s" % (myname,lazy))
            # always start with a clean slate
            shutil.rmtree(themsname,True)
            shutil.rmtree('referemce.ms',True)
            # use importasdm, which always looks for and skips duplicate DATA rows.
            lsdm = sdm(myasdmname)
            self.res = lsdm.toms(vis=themsname, ocorr_mode="ao", bdfflags=True, lazy=lazy, overwrite=True)
            self.assertTrue(self.res)
            print("%s: Success! Now checking output ..." % myname)
            mscomponents = set(["ANTENNA/table.dat",
                                "CALDEVICE/table.dat",
                                "DATA_DESCRIPTION/table.dat",
                                "FEED/table.dat",
                                "FIELD/table.dat",
                                "FLAG_CMD/table.dat",
                                "HISTORY/table.dat",
                                "OBSERVATION/table.dat",
                                "POINTING/table.dat",
                                "POLARIZATION/table.dat",
                                "PROCESSOR/table.dat",
                                "SOURCE/table.dat",
                                "SPECTRAL_WINDOW/table.dat",
                                "STATE/table.dat",
                                "SYSCAL/table.dat",
                                "SYSPOWER/table.dat",
                                "WEATHER/table.dat",
                                "ANTENNA/table.f0",
                                "CALDEVICE/table.f0",
                                "DATA_DESCRIPTION/table.f0",
                                "FEED/table.f0",
                                "FIELD/table.f0",
                                "FLAG_CMD/table.f0",
                                "HISTORY/table.f0",
                                "OBSERVATION/table.f0",
                                "POINTING/table.f0",
                                "POLARIZATION/table.f0",
                                "PROCESSOR/table.f0",
                                "SOURCE/table.f0",
                                "SPECTRAL_WINDOW/table.f0",
                                "STATE/table.f0",
                                "SYSCAL/table.f0",
                                "SYSPOWER/table.f0",
                                "WEATHER/table.f0"
                                ])
            for name in mscomponents:
                if not os.access(themsname+"/"+name, os.F_OK):
                    print("%s: Error  %s/%s doesn't exist ..." % (myname,themsname,name))
                    retValue['success']=False
                    retValue['error_msgs']=retValue['error_msgs']+themsname+'/'+name+' does not exist'
                else:
                    print("%s: %spresent." % (myname,name))
            print("%s: MS exists. All tables present. Try opening as MS ..." % myname)
            try:
                mslocal.open(themsname)
                print("%s: MS can be opened" % myname)
                mslocal.close()

            except:
                print("%s: Error  Cannot open MS table %s" % (myname,themsname))
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+themsname
            else:
                print("%s: OK. Generating a reference MS with first integration checking turned off" % myname)

                # this must be done using asdm2MS and bdflags2MS directly
                asdm2MScmd = 'asdm2MS --ocm "ao" --checkdupints false'
                if lazy:
                    asdm2MScmd = asdm2MScmd + " -lazy"
                asdm2MScmd = asdm2MScmd + " " + myasdmname + " reference.ms"
                print('%s Running asdm2MS standalone invoked as:' % myname)
                print(asdm2MScmd)
                exitcode = os.system(asdm2MScmd)
                if exitcode != 0:
                    print("%s asdm2MS terminated with exit code %s" % (myname,exitcode))
                    retValue['success'] = False
                    retValue['error_msgs']=retValue['error_msgs']+' standalone execution of asdm2MS failed'
                    # this should break out of the main loop over lazy values
                    break

                bdflags2MScmd = 'bdflags2MS -f ALL --ocm "ao" --checkdupints false'
                if lazy:
                    bdflags2MScmd = bdflags2MScmd + " -lazy=true"
                bdflags2MScmd = bdflags2MScmd + " " + myasdmname + " reference.ms"
                print('%s Running bdflags2MS standalone invoked as:' % myname)
                print(bdflags2MScmd)
                exitcode = os.system(bdflags2MScmd)
                if exitcode != 0:
                    print("%s bdflags2MS terminated with exit code %s" % (myname,exitcode))
                    retValue['success'] = False
                    retValue['error_msgs']=retValue['error_msgs']+' standalone execution of bdflags2MS failed'
                    # this should break out of the main loop over lazy values
                    break

                # at this point, reference.ms should exist, with all of the auto rows (ao), including duplicates, with BDF flags applied
                if(os.path.exists('reference.ms')):
                    # and the test here is to make sure that the two MSs are identical in the Main table. The same values in
                    # the same order except for gaps where rows from the larger table are skipped.

                    # expected size (rows)
                    msSize = 31972
                    refSize = 32040

                    # expected gaps start at these rows and are always 4 rows long, row numbers in reference.ms
                    gaps = [2280,3176,5328,7120,9276,10172,11432,12328,14480,15376,18752,19648,21800,26636,27896,28792,30944]

                    mstb = table( )
                    mstb.open(themsname)
                    if mstb.nrows() != msSize:
                        print('%s MS size is not of the expected number of rows : %s != %s' % (myname,mstb.nrows(),msSize))
                        retValue['success'] = False
                        retValue['error_msgs'] = retValue['error_msgs'] + 'bad size for MS'

                    reftb = table( )
                    reftb.open('reference.ms')
                    if reftb.nrows() != refSize:
                        print('%s Reference MS size is not of the expected number of rows : %s != %s' % (myname,reftb.nrows(),refSize))
                        retValue['success'] = False
                        retValue['error_msgs'] = retValue['error_msgs'] + 'bad size for reference MS'

                    if retValue['success']:
                        refrow = 0
                        msrow = 0
                        gapStart = -1
                        gapIndex = -1
                        # assumes set of column names is the same
                        cols = mstb.colnames()

                        while((msrow < mstb.nrows()) and (refrow < reftb.nrows())):
                            matchFound = True
                            for col in cols:
                                if mstb.iscelldefined(col,msrow) != reftb.iscelldefined(col,refrow):
                                    # defined in one cell, not defined in the other
                                    matchFound = False
                                else:
                                    if (mstb.iscelldefined(col,msrow)):
                                        msVal = mstb.getcell(col,msrow)
                                        refVal = reftb.getcell(col,refrow)
                                        # assumes the type is the same for the same column in both tables
                                        if isinstance(msVal,numpy.ndarray):
                                            matchFound = numpy.array_equal(msVal,refVal)
                                        else:
                                            matchFound = msVal == refVal
                                if not matchFound:
                                    # no point in checking the other columns
                                    break
                            if matchFound:
                                if gapStart >= 0:
                                    # a gap has ended, verify that it was exactly 4 rows long
                                    if (refrow-gapStart) != 4:
                                        print('%s Unexpected gap length not equal to 4 rows. Gap length = %s starting at row %s' % (myname,refrow-gapStart,refRow))
                                        retValue['success'] = False
                                        retValue['error_msg'] = 'Unexpected gap length not equal to 4 rows'
                                    gapStart = -1
                                msrow = msrow + 1
                            else:
                                # do not increment msrow in this case, keep looking for a match
                                if gapStart < 0:
                                    # new gap, increment index and verify it's at the expected place
                                    gapIndex = gapIndex+1
                                    gapStart = refrow
                                    if gapIndex > len(gaps):
                                        print('%s Unexpected gap seen past end of known gaps. Starting at row %s' % (myname,refrow))
                                        retValue['success'] = False
                                        retValue['error_msg'] = 'Unexpected gap after end of known gaps'
                                    else:
                                        if gapStart != gaps[gapIndex]:
                                            print('%s Unexpected gap start at row %s expected at row %s' % (myname,gapStart,gaps[gapIndex]))
                                            retValue['success'] = False
                                            retValue['error_msg'] = 'Unexpected gap start row'
                            # refrow is always incremented
                            refrow = refrow + 1

                            # bail out on failure
                            if not retValue['success']:
                                break

                    reftb.close()
                    mstb.close()
            # break out of the lazy loop on failure
            if not retValue['success']:
                break

        self.assertTrue(retValue['success'],retValue['error_msgs'])

    def test7_bdflags1(self):
        '''Asdm-import: Test good 12 m ASDM with mixed pol/channelisation input with default filler selecting "co" on output and using the BDF flags'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }

        myasdmname = 'uid___A002_X71e4ae_X317_short'
        themsname = myasdmname+".ms"
        mysdm = sdm(myasdmname)
        self.res = mysdm.toms(vis=themsname, ocorr_mode="co", bdfflags=True)
        self.assertTrue(self.res)
        print("%s: Success! Now checking output ..." % myname)
        mscomponents = set(["ANTENNA/table.dat",
                            "DATA_DESCRIPTION/table.dat",
                            "FEED/table.dat",
                            "FIELD/table.dat",
                            "FLAG_CMD/table.dat",
                            "HISTORY/table.dat",
                            "OBSERVATION/table.dat",
                            "POINTING/table.dat",
                            "POLARIZATION/table.dat",
                            "PROCESSOR/table.dat",
                            "SOURCE/table.dat",
                            "SPECTRAL_WINDOW/table.dat",
                            "STATE/table.dat",
                            "SYSCAL/table.dat",
                            "ANTENNA/table.f0",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.f0",
                            "FIELD/table.f0",
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.f0",
                            "POINTING/table.f0",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.f0",
                            "SOURCE/table.f0",
                            "SPECTRAL_WINDOW/table.f0",
                            "STATE/table.f0",
                            "SYSCAL/table.f0"
                            ])
        for name in mscomponents:
            if not os.access(themsname+"/"+name, os.F_OK):
                print("%s: Error  %s/%s doesn't exist ..." % (myname,themsname,name))
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+themsname+'/'+name+' does not exist'
            else:
                print("%s: %s present." % (myname,name))
        print("%s: MS exists. All tables present. Try opening as MS ..." % myname)
        try:
            mslocal.open(themsname)
        except:
            print("%s: Error  Cannot open MS table %s" % (myname,themsname))
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+themsname
        else:
            mslocal.close()
            print("%s: OK. Checking tables in detail ..." % myname)
            fromsdm = sdm(myasdmname)
            fromsdm.toms(vis='reference.ms', overwrite=True, ocorr_mode="co", bdfflags=False)

            if(os.path.exists('reference.ms')):
                retValue['success'] = th.checkwithtaql("select from [select from reference.ms orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "
                                                    +themsname+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(near(t1.DATA,t2.DATA, 1.e-06)))") == 0
                if not retValue['success']:
                    print("ERROR: DATA does not agree with reference.")
                else:
                    print("DATA columns agree.")

                retValueTmp = th.checkwithtaql("select from [select from reference.ms orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "
                                            +themsname+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(t1.FLAG==t2.FLAG)) ") != 0
                if not retValueTmp:
                    print("ERROR: FLAG columns do agree with reference but they shouldn't.")
                else:
                    print("FLAG columns do not agree as expected.")

                retValue['success'] = retValue['success'] and retValueTmp

                for subtname in ["ANTENNA",
                                 "DATA_DESCRIPTION",
                                 "FEED",
                                 "FIELD",
                                 "FLAG_CMD",
                                 "OBSERVATION",
                                 "POLARIZATION",
                                 "PROCESSOR",
                                 "SOURCE",
                                 "SPECTRAL_WINDOW",
                                 "STATE",
                                 "SYSCAL"]:

                    print("\n*** Subtable %s" % subtname)
                    excllist = []
                    if subtname=='SOURCE':
                        excllist=['POSITION', 'TRANSITION', 'REST_FREQUENCY', 'SYSVEL']
                    if subtname=='SYSCAL':
                        excllist=['TANT_SPECTRUM', 'TANT_TSYS_SPECTRUM']
                    if subtname=='SPECTRAL_WINDOW':
                        excllist=['CHAN_FREQ', 'CHAN_WIDTH', 'EFFECTIVE_BW', 'RESOLUTION', 'ASSOC_SPW_ID', 'ASSOC_NATURE']
                        for colname in excllist:
                            if colname!='ASSOC_NATURE':
                                retValue['success'] = th.compVarColTables('reference.ms/SPECTRAL_WINDOW',
                                                                          themsname+'/SPECTRAL_WINDOW', colname, 0.01) and retValue['success']
                    if subtname=='POLARIZATION':
                        excllist=['CORR_TYPE', 'CORR_PRODUCT']
                        for colname in excllist:
                            retValue['success'] = th.compVarColTables('reference.ms/POLARIZATION',
                                                                      themsname+'/POLARIZATION', colname, 0.01) and retValue['success']
                    try:
                        retValue['success'] = th.compTables('reference.ms/'+subtname,
                                                            themsname+'/'+subtname, excllist,
                                                            0.01) and retValue['success']
                    except:
                        retValue['success'] = False
                        print("ERROR for table %s" % subtname)
        self.assertTrue(retValue['success'],retValue['error_msgs'])


class asdm_import8(test_base):
    # these are more like unit tests, difficult to test without invoking all of importasdm
    # currently this is just tests on SDM_NUM_BIN
    
    def setUp(self):
        self.setUp_numbin()

    def tearDown(self):
        for this_asdm_name in ['alma_numbin_mixed','evla_numbin_2','evla_numbin_4']:
            os.system('rm -rf '+this_asdm_name+"*")

    def doNumTest(self, testName, asdm_name, ms_name, spWin_name, execBlock_name, expWinFunCol, expNumBinCol, expResCol):
        retValue = {'success': True, 'error_msgs': '' } 
        print(testName,": testing SDM columns in",asdm_name,' writing to ', ms_name)

        originalSpWin = None
        originalExecBlock = None

        if spWin_name is not None:
            print(testName,": using",spWin_name,"for SpectralWindow.xml")
            spWin_path = os.path.join(asdm_name,spWin_name)
            originalSpWin = os.path.join(asdm_name,"SpectralWindow.xml.original")
            if not os.path.exists(spWin_path):
                msg = spWin_path+" does not exist"
                retValue['success'] = False
                retValue['error_msgs'] = msg
                return retValue
            if os.path.exists(originalSpWin):
                msg = originalSpWin+" already exists, will not overwrite"
                retValue['success'] = False
                retValue['error_msgs'] = msg
                return retValue
            shutil.move(os.path.join(asdm_name,'SpectralWindow.xml'),originalSpWin)
            shutil.copyfile(spWin_path,os.path.join(asdm_name,'SpectralWindow.xml'))

        if execBlock_name is not None:
            print(testName,": using",execBlock_name,"for ExecBlock.xml")
            execBlock_path = os.path.join(asdm_name,execBlock_name)
            originalExecBlock = os.path.join(asdm_name,"ExecBlock.xml.original")
            if not os.path.exists(execBlock_path):
                msg = execBlock_path+" does not exist"
                retValue['success'] = False
                retValue['error_msgs'] = msg
                return retValue
            if os.path.exists(originalExecBlock):
                msg = originalExecBlock+" already exists, will not overwrite"
                retValue['success'] = False
                retValue['error_msgs'] = msg
                return retValue
            shutil.move(os.path.join(asdm_name,'ExecBlock.xml'),originalExecBlock)
            shutil.copyfile(execBlock_path,os.path.join(asdm_name,'ExecBlock.xml'))
 
        mysdm = sdm(asdm_name)
        self.res = mysdm.toms(vis=ms_name,lazy=True,process_syspower=False,process_caldevice=False,process_pointing=False,process_flags=False)

        # the only table this test cares about is SPECTRAL_WINDOW
        spwName = os.path.join(ms_name,"SPECTRAL_WINDOW")
        if not os.access(spwName,os.F_OK):
            print(testName,": Error ", spwName, "doesn't exist ...")
            retValue['success'] = False
            retValue['error_msgs']=spwName+' does not exist'
        else:
            ok = tblocal.open(spwName)
            if (ok):
                try:
                    winFunCol = tblocal.getcol('SDM_WINDOW_FUNCTION')
                    if not numpy.all(winFunCol==expWinFunCol):
                        retValue['success'] = False
                        msg = "ERROR Unexpected SDM_WINDOW_FUNCTION values when filling "+asdm_name
                        retValue['error_msgs']=msg
                        print(testName,":", msg)
                except:
                     retValue['success'] = False
                     msg = "ERROR getting/testing SDM_WINDOW_FUNCTION column in "+spwName
                     retValue['error_msgs']=msg
                     print(testName,":",msg)

                try:
                    numBinCol = tblocal.getcol('SDM_NUM_BIN')
                    if not numpy.all(numBinCol==expNumBinCol):
                        retValue['success'] = False
                        msg = "ERROR Unexpected SDM_NUM_BIN values when filling "+asdm_name
                        # there may already be messages in error_msgs
                        if len(retValue['error_msgs']>0):
                            retValue['error_msgs']=retValue['error_msgs']+'\n'
                        retValue['error_msgs']=retValue['error_msgs']+msg
                        print(testName,":",msg)            
                except:
                    retValue['success'] = False
                    msg = "ERROR getting/testing SDM_NUM_BIN column in "+spwName
                    # there may already be messages in error_msgs
                    if len(retValue['error_msgs'])>0:
                        retValue['error_msgs']=retValue['error_msgs']+'\n'
                    retValue['error_msgs']=retValue['error_msgs']+msg
                    print(testName,":",msg)

                # only test RESOLUTION values if expResCol is not None
                if expResCol is not None:
                    try:
                        resCol = tblocal.getcol('RESOLUTION')
                        # only test first value in each row, assumes all values in a row are equal
                        resCol = resCol[0,:]
                        if not numpy.all(resCol==expResCol):
                            retValue['success'] = False
                            msg = "ERROR Unexpected RESOLUTION values when filling "+asdm_name
                            # there may already be messages in error_msgs
                            if len(retValue['error_msgs']>0):
                                retValue['error_msgs']=retValue['error_msgs']+'\n'
                            retValue['error_msgs']=retValue['error_msgs']+msg
                            print(testName,":",msg)         
                    except:
                        retValue['success'] = False
                        msg = "ERROR getting/testing RESOLUTION column in "+spwName
                        # there may already be messages in error_msgs
                        if len(retValue['error_msgs'])>0:
                            retValue['error_msgs']=retValue['error_msgs']+'\n'
                        retValue['error_msgs']=retValue['error_msgs']+msg
                        print(testName,":",msg)
                tblocal.close()

            else:
                msg = "ERROR opening",spwName
                retValue['success'] = False
                retValue['error_msgs'] = msg
                print(testName,":",msg)

        if originalSpWin is not None:
            os.remove(os.path.join(asdm_name,'SpectralWindow.xml'))
            shutil.move(originalSpWin,os.path.join(asdm_name,'SpectralWindow.xml'))
            print(testName,": restored original SpectralWindow.xml")

        if originalExecBlock is not None:
            os.remove(asdm_name+'/ExecBlock.xml')
            shutil.move(originalExecBlock,os.path.join(asdm_name,'ExecBlock.xml'))
            print(testName,": restored original ExecBlock.xml")

        return retValue 

    def test_alma_numbin(self):
        retValue = {'success': True, 'error_msgs': '' } 

        # original SpectralWindow.xml and inferred numBin  values
        asdm_name = 'alma_numbin_mixed'
        ms_name = asdm_name+".ms"
        expWinFunCol = numpy.array(['UNIFORM']*5 + ['HANNING']*8 + ['UNIFORM']*4 + ['HANNING']*18 + ['UNIFORM']*42)
        # expected values, 8 @ 24, 2 @ 27,29,31,33, rest are 1
        expNumBinCol = numpy.ones(77,dtype=numpy.int32)
        expNumBinCol[25] = 8
        for indx in [27,29,31,33]:
            expNumBinCol[indx] = 2
        res = self.doNumTest(myname,asdm_name,ms_name,None,None,expWinFunCol,expNumBinCol,None)
        retValue['success'] = res['success']
        retValue['error_msgs'] = res['error_msgs']

        # SpectralWindow.xml with appropriate numBin values, should yield same column values
        ms_name = asdm_name+".numbin.ms"
        res = self.doNumTest(myname,asdm_name,ms_name,'SpectralWindow.xml.numBin',None,expWinFunCol,expNumBinCol,None)
        retValue['success'] = retValue['success'] and res['success']
        retValue['error_msgs'] = retValue['error_msgs'] + res['error_msgs']

        # SpectralWindow.xml with faked resolution and expectedBw values but no numBin, tests other inferred numBin values
        expNumBinCol[5] = 4
        expNumBinCol[7] = 16
        ms_name = asdm_name+".faked.ms"
        res = self.doNumTest(myname,asdm_name,ms_name,'SpectralWindow.xml.faked',None,expWinFunCol,expNumBinCol,None)
        retValue['success'] = retValue['success'] and res['success']
        retValue['error_msgs'] = retValue['error_msgs'] + res['error_msgs']

        # SpectralWindow.xml with faked resolution and expectedBw values and added numBin values, same expected values as previous test
        ms_name = asdm_name+".faked.numBin.ms"
        res = self.doNumTest(myname,asdm_name,ms_name,'SpectralWindow.xml.faked.numBin',None,expWinFunCol,expNumBinCol,None)
        retValue['success'] = retValue['success'] and res['success']
        retValue['error_msgs'] = retValue['error_msgs'] + res['error_msgs']

        self.assertTrue(retValue['success'],retValue['error_msgs'])

    def test_evla_numbin(self):
        retValue = {'success':True, 'error_msgs':''}

        # numbin=2 related tests
        sdm_name = 'evla_numbin_2'

        # original SpectralWindow.xml and inferred numBin values, all equal to 2
        expWinFunCol = numpy.array(['UNIFORM']*16)
        expNumBinCol = numpy.empty(16,dtype=numpy.int32)
        expNumBinCol.fill(2)
        # also should alter resolution to these expected values
        expResCol = numpy.empty(16)
        expResCol.fill(4000000.)
        ms_name = sdm_name+".ms"
        res = self.doNumTest(myname,sdm_name,ms_name,None,None,expWinFunCol,expNumBinCol,expResCol)
        retValue['success'] = res['success']
        retValue['error_msgs'] = res['error_msgs']

        # SpectralWindow.xml with numBin field and appropriately modified resolution, same expected values
        ms_name = sdm_name+".numBin.ms"
        res = self.doNumTest(myname,sdm_name,ms_name,'SpectralWindow.xml.numBin',None,expWinFunCol,expNumBinCol,expResCol)
        retValue['success'] = res['success']
        retValue['error_msgs'] = res['error_msgs']

        # SpectralWindow.xml with mostly numBin and alterned resolution, but one row has the original values, same expected values
        ms_name = sdm_name+".mixed.ms"
        res = self.doNumTest(myname,sdm_name,ms_name,'SpectralWindow.xml.mixed',None,expWinFunCol,expNumBinCol,expResCol)
        retValue['success'] = res['success']
        retValue['error_msgs'] = res['error_msgs']

        # original SpectralWindow.xml but with one non-physical (bad) value of resolution leading to the algorithm giving up and numBin=1 and resolution the original bad value
        ms_name = sdm_name+".bad.ms"
        expNumBinCol[0] = 1
        expResCol[0] = 9000000.
        res = self.doNumTest(myname,sdm_name,ms_name,'SpectralWindow.xml.bad',None,expWinFunCol,expNumBinCol,expResCol)
        retValue['success'] = res['success']
        retValue['error_msgs'] = res['error_msgs']

        # numbin=4 related tests
        sdm_name = 'evla_numbin_4'
        
        # original SpectralWindow.xml and inferred numBin values, all equal to 4
        expWinFunCol = numpy.array(['UNIFORM']*16)
        expNumBinCol = numpy.empty(16,dtype=numpy.int32)
        expNumBinCol.fill(4)
        # also should alter resolution to these expected values
        expResCol = numpy.empty(16)
        expResCol.fill(8000000.)
        ms_name = sdm_name+".ms"
        res = self.doNumTest(myname,sdm_name,ms_name,None,None,expWinFunCol,expNumBinCol,expResCol)
        retValue['success'] = res['success']
        retValue['error_msgs'] = res['error_msgs']

        # original SpectralWindow.xml with numBin field and altered resolution, same expected values
        ms_name = sdm_name+".numBin.ms"
        res = self.doNumTest(myname,sdm_name,ms_name,'SpectralWindow.xml.numBin',None,expWinFunCol,expNumBinCol,expResCol)
        retValue['success'] = res['success']
        retValue['error_msgs'] = res['error_msgs']

        # original SpectralWindow.xml with numBin field and original resolution
        # expected numBin is the same, expected resolution is now the original values
        expResCol /= 4.0
        ms_name = sdm_name+".onlyNumBin.ms"
        res = self.doNumTest(myname,sdm_name,ms_name,'SpectralWindow.xml.onlyNumBin',None,expWinFunCol,expNumBinCol,expResCol)
        retValue['success'] = res['success']
        retValue['error_msgs'] = res['error_msgs']

        # original SpectralWindow.xml and altered ExecBlock so that the telescope is UNKNOWN
        # numBin is all 1 and expected resolution is the original resolution
        expNumBinCol.fill(1)
        ms_name = sdm_name+".unknownTel.ms"
        res = self.doNumTest(myname,sdm_name,ms_name,None,'ExecBlock.xml.unknownTel',expWinFunCol,expNumBinCol,expResCol)
        retValue['success'] = res['success']
        retValue['error_msgs'] = res['error_msgs']

        self.assertTrue(retValue['success'],retValue['error_msgs'])

def suite():
    ### asdm_import4 exists in the importasdm test, it involves flagging - which is largely
    ### only implemented in the task, not this tool
    return [asdm_import1, asdm_import2, asdm_import3, asdm_import5, asdm_import6, asdm_import7, asdm_import8]

if __name__ == '__main__':
    unittest.main()
