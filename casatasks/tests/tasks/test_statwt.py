import os
import sys
import shutil
import unittest
import numpy as np
import numpy.ma as ma

from casatools import ctsys, table, ms
from casatasks import statwt

datadir = ctsys.resolve('regression/unittest/statwt')
src = os.path.join(datadir,'ngc5921.split_2.ms')

def get_weights(data, flags, exposures, combine_corr, target_exposure):  
    shape = data.shape
    ncorr_groups = 1 if combine_corr else shape[0]
    ncorr = shape[0]
    weights = np.zeros([shape[0], shape[1]])
    nrows = data.shape[2]
    for corr in range(ncorr_groups):
        end = corr + 1 if ncorr_groups > 1 else ncorr + 1
        var = variance(data[corr:end, :, :], flags[corr:end, :, :], exposures)
        if var == 0:
            weights[corr:end, :] = 0 
        else:
            weights[corr:end, :] = target_exposure/var
    return weights

def variance(data, flags, exposures):
    if flags.all():
        return 0
    expo = ma.masked_array(np.resize(exposures, data.shape), mask=flags)
    d = ma.array(data, mask=flags)
    myreal = np.real(d)
    myimag = np.imag(d)
    mean_r = np.sum(expo*myreal)/np.sum(expo)
    mean_i = np.sum(expo*myimag)/np.sum(expo)
    var_r = np.sum(expo * (myreal - mean_r)*(myreal - mean_r))/d.count()
    var_i = np.sum(expo * (myimag - mean_i)*(myimag - mean_i))/d.count()
    return (var_r + var_i)/2

def _get_dst_cols(dst, other="", dodata=True):
    mytb = table()
    mytb.open(dst)
    wt = mytb.getcol("WEIGHT")
    wtsp = mytb.getcol("WEIGHT_SPECTRUM")
    flag = mytb.getcol("FLAG")
    frow = mytb.getcol("FLAG_ROW")
    if dodata:
        data = mytb.getcol("CORRECTED_DATA")
    if len(other) > 0:
        if type(other) == type([]):
            othercol = []
            for x in other:
                othercol.append(mytb.getcol(x))
        else:
            othercol = mytb.getcol(other)
    mytb.close()
    cols = [wt, wtsp, flag, frow]
    if dodata:
        cols.append(data)
    if len(other) > 0:
        if type(other) == type([]):
            for x in othercol:
                cols.append(x)
        else:
            cols.append(othercol)
    return cols

def _get_table_cols(mytb):
    times = mytb.getcol("TIME")
    wt = mytb.getcol("WEIGHT")
    wtsp = mytb.getcol("WEIGHT_SPECTRUM")
    flag = mytb.getcol("FLAG")
    frow = mytb.getcol("FLAG_ROW")
    data = mytb.getcol("CORRECTED_DATA")
    return [times, wt, wtsp, flag, frow, data]

# combine correlations
def _variance(dr, di, flag, row):
    fr = numpy.extract(numpy.logical_not(flag[:,:,row]), dr[:,:,row])
    fi = numpy.extract(numpy.logical_not(flag[:,:,row]), di[:,:,row])
    if len(fr) <= 1:
        return 0
    else:
        vr = numpy.var(fr, ddof=1)
        vi = numpy.var(fi, ddof=1)
        return 2/(vr + vi)

# per correlation
def _variance2(dr, di, flag, corr, row):
    fr = numpy.extract(numpy.logical_not(flag[corr,:,row]), dr[corr,:,row])
    fi = numpy.extract(numpy.logical_not(flag[corr,:,row]), di[corr,:,row])
    if len(fr) <= 1:
        return 0
    else:
        vr = numpy.var(fr, ddof=1)
        vi = numpy.var(fi, ddof=1)
        return 2/(vr + vi)


class statwt_test(unittest.TestCase):
    def _check_weights(
        self, msname, mode, data_column, chan_flags, combine_corr
    ):
        if data_column.startswith('c'):
            colname = 'CORRECTED_DATA'
        elif data_column.startswith('d'):
            colname = 'DATA'
        else:
            raise Exception("Unhandled column spec " + data_column)
        if not mode.startswith('one'):
            raise Exception("Unhandled mode")
        tb = table()
        for ant1 in range(10):
            for ant2 in range((ant1 + 1), 10):
                query_str = 'ANTENNA1=' + str(ant1) + ' AND ANTENNA2=' + str(ant2)
                tb.open(msname)
                subt = tb.query(query_str)
                data = subt.getcol(colname)
                flags = subt.getcol('FLAG')
                exposures = subt.getcol('EXPOSURE')
                wt = subt.getcol('WEIGHT')
                wtsp = subt.getcol('WEIGHT_SPECTRUM')
                wt = subt.getcol('WEIGHT')
                subt.done()
                tb.done()
                if type(chan_flags) != type(None):
                    t_flags = np.expand_dims(np.expand_dims(chan_flags, 0), 2)
                    flags = np.logical_or(flags, t_flags)
                nrows = data.shape[2]
                for row in range(nrows):
                    if mode.startswith('one'):
                        start = row
                        end = row+1
                    if mode.startswith('one'):
                        start = row
                        end = row+1
                    # print 'baseline ' + str([ant1, ant2]) + ' row ' + str(row)
                    weights = get_weights(
                        data[:,:,start:end], flags[:, :, start:end],
                        exposures[start: end], combine_corr, exposures[row]
                    )
                    self.assertTrue(
                        np.allclose(weights, wtsp[:, :, row]), 'Failed wtsp, got '
                        + str(wtsp[:, :, row]) + '\nexpected ' + str(weights)
                        + '\nbaseline ' + str([ant1, ant2])
                        + '\nrow ' + str(row)
                    )
                    self.assertTrue(
                        np.allclose(np.median(weights, 1), wt[:, row]),
                        'Failed weight, got ' + str(wt[:, row])
                        + '\nexpected ' + str(np.median(weights, 1))
                        + '\nbaseline ' + str([ant1, ant2]) + '\nrow '
                        + str(row)
                    )
    def test_algorithm(self):
        """ Test the algorithm, includes excludechans tests"""
        mytb = table()
        mytb.open(src)
        expflag = mytb.getcol("FLAG")
        expfrow = mytb.getcol("FLAG_ROW")
        mytb.done()
        dst = "ngc5921.split.ms"
        rtol = 1e-7
        cflags = np.array(63 * [False])
        cflags[10:21] = True
        for combine in ["", "corr"]:
            c = 0
            for fitspw in ["0:0~9;21~62", "", "0:10~20"]:
                shutil.copytree(src, dst) 
                excludechans = c == 2
                statwt(
                    vis=dst, combine=combine, fitspw=fitspw,
                    excludechans=excludechans
                )
                chan_flags = cflags if fitspw else None
                self._check_weights(
                    dst, mode='one_to_one', data_column='c',
                    chan_flags=chan_flags, combine_corr=bool(combine)
                )
                shutil.rmtree(dst)
                c += 1 
    """
    def test_algorithm(self):
        "" Test the algorithm, includes fitspw, excludechans tests""
        mytb = table()
        mytb.open(src)
        expflag = mytb.getcol("FLAG")
        expfrow = mytb.getcol("FLAG_ROW")
        mytb.done()
        dst = "ngc5921.split.ms"
        rtol = 1e-7
        for combine in ["", "corr"]:
            c = 0
            for fitspw in ["0:0~9;21~63", "", "0:10~20"]:
                shutil.copytree(src, dst)
                excludechans = c == 2
                statwt(dst, combine=combine, fitspw=fitspw, excludechans=excludechans)
                [wt, wtsp, flag, frow, data] = _get_dst_cols(dst)
                actflag = flag.copy()
                if fitspw != "":
                    actflag[:, 10:21, :] = True
                dr = numpy.real(data)
                di = numpy.imag(data)
                myshape = wtsp.shape
                ncorr = myshape[0]
                nrow = myshape[2]
                if (combine == "corr"):
                    for row in range(nrow):
                        expec = _variance(dr, di, actflag, row)
                        self.assertTrue(
                            numpy.all(numpy.isclose(wt[:, row], expec, rtol=rtol)),
                            "WEIGHT fail at row" + str(row) + ". got: " + str(wt[:, row])
                            + " expec " + str(expec)
                        )
                        self.assertTrue(
                            len(numpy.unique(wtsp[:,:,row])) == 1,
                            "Weight values are not the same"
                        )
                        self.assertTrue(
                            numpy.all(numpy.isclose(wtsp[:,:,row], expec, rtol)),
                            "Incorrect weights"
                        )
                        if expec == 0:
                            self.assertTrue(numpy.all(flag[:,:,row]), "Not all flags are true")
                            self.assertTrue(frow[row], "FLAG_ROW is not true")
                        else:
                            self.assertTrue(
                                numpy.all(flag[:,:,row] == expflag[:,:,row]),
                                "FLAGs don't match"
                            )
                            self.assertTrue(
                                frow[row] == expfrow[row], "FLAG_ROW doesn't match"
                            )
                else:
                    for row in range(nrow):
                        for corr in range(ncorr):
                            expec = _variance2(dr, di, actflag, corr, row)
                            self.assertTrue(
                                numpy.isclose(wt[corr, row], expec, rtol=rtol),
                                "WEIGHT fail at row" + str(row) + ". got: " + str(wt[corr, row])
                                + " expec " + str(expec)
                            )
                            self.assertTrue(
                                len(numpy.unique(wtsp[corr,:,row])) == 1,
                                "Weight values are not the same"
                            )
                            self.assertTrue(
                                numpy.all(
                                    numpy.isclose(wtsp[corr,:,row], expec, rtol)), "Incorrect weights"
                            )
                            if expec == 0:
                                self.assertTrue(numpy.all(flag[corr,:,row]), "Not all flags are true")
                            else:
                                self.assertTrue(
                                    numpy.all(flag[corr,:,row] == expflag[corr,:,row]),
                                    "FLAGs don't match at row " + str(row) + ".\n"
                                    "Expected: " + str(expflag[corr,:,row]) + "\n"
                                    "Got     : " + str(flag[corr,:,row])
                                )
                        if (numpy.all(flag[:,:,row])):
                            self.assertTrue(frow[row], "FLAG_ROW is not true")
                        else:
                            self.assertFalse(frow[row], "FLAG_ROW is not false")
                shutil.rmtree(dst) 
                c += 1
    """
              
    def test_timebin(self):
        """ Test time binning"""
        dst = "ngc5921.split.timebin.ms"
        ref = os.path.join(datadir,"ngc5921.timebin300s_2.ms.ref")
        [refwt, refwtsp, refflag, reffrow, refdata] = _get_dst_cols(ref)
        rtol = 1e-7
        combine = "corr"
        for timebin in ["300s", 10]:
            shutil.copytree(src, dst) 
            statwt(dst, timebin=timebin, combine=combine)
            [tstwt, tstwtsp, tstflag, tstfrow, tstdata] = _get_dst_cols(dst)
            self.assertTrue(numpy.all(tstflag == refflag), "FLAGs don't match")
            self.assertTrue(numpy.all(tstfrow == reffrow), "FLAG_ROWs don't match")
            self.assertTrue(numpy.all(numpy.isclose(tstwt, refwt, rtol)), "WEIGHTs don't match")
            self.assertTrue(numpy.all(numpy.isclose(tstwtsp, refwtsp, rtol)), "WEIGHT_SPECTRUMs don't match")
            shutil.rmtree(dst)

    def test_chanbin(self):
        """Test channel binning"""
        dst = "ngc5921.split.chanbin_0.ms"
        rtol = 1e-7
        for combine in ["", "corr"]:
            if combine == "":
                ref = os.path.join(datadir,"ngc5921.chanbin_sepcorr_2.ms.ref")
            else:
                ref = os.path.join(datadir,"ngc5921.chanbin_combcorr_2.ms.ref")
            [refwt, refwtsp, refflag, reffrow, refdata] = _get_dst_cols(ref)
            for i in [0, 1]:
                for chanbin in ["195.312kHz", 8]:
                    if i == 1 and (combine == "corr" or chanbin == 8):
                        continue
                    shutil.copytree(src, dst)
                    if i == 0:
                        statwt(dst, chanbin=chanbin, combine=combine)
                    elif i == 1:
                        # check WEIGHT_SPECTRUM is created, only check once,
                        # this test is long as it is
                        # shutil.copytree(src, dst)
                        mytb = table()
                        mytb.open(dst, nomodify=False)
                        x = mytb.ncols()
                        self.assertTrue(mytb.removecols("WEIGHT_SPECTRUM"), "column not removed")
                        y = mytb.ncols()
                        self.assertTrue(y == x-1, "wrong number of columns")
                        mytb.done()
                        statwt(dst, chanbin=chanbin, combine=combine)
                    [tstwt, tstwtsp, tstflag, tstfrow, tstdata] = _get_dst_cols(dst)
                    self.assertTrue(numpy.all(tstflag == refflag), "FLAGs don't match")
                    self.assertTrue(numpy.all(tstfrow == reffrow), "FLAG_ROWs don't match")
                    self.assertTrue(numpy.all(numpy.isclose(tstwt, refwt, rtol)), "WEIGHTs don't match")
                    self.assertTrue(numpy.all(numpy.isclose(tstwtsp, refwtsp, rtol)), "WEIGHT_SPECTRUMs don't match")
                    shutil.rmtree(dst)

    def test_minsamp(self):
        """Test minimum number of points"""
        dst = "ngc5921.split.minsamp.ms"
        combine = "corr"
        trow = 12
        for minsamp in [60, 80]:
            shutil.copytree(src, dst)
            statwt(dst, minsamp=minsamp, combine=combine)
            [wt, wtsp, flag, frow, data] = _get_dst_cols(dst)
            if minsamp == 60:
                self.assertTrue((wt[:, trow] > 0).all(), "Incorrect weight row " + str(trow))
                self.assertTrue((wtsp[:, :, trow] > 0).all(), "Incorrect weight spectrum row " + str(trow))
                self.assertFalse(flag[:,:,trow].all(), "Incorrect flag row " + str(trow))
                self.assertFalse(frow[trow], "Incorrect flagrow row " + str(trow))
            else:
                self.assertTrue((wt[:, trow] == 0).all(), "Incorrect weight row " + str(trow))
                self.assertTrue((wtsp[:, :, trow] == 0).all(), "Incorrect weight spectrum row " + str(trow))
                self.assertTrue(flag[:,:,trow].all(), "Incorrect flag row " + str(trow))
                self.assertTrue(frow[trow], "Incorrect flagrow row " + str(trow))
            shutil.rmtree(dst)
            
    def test_fieldsel(self):
        """Test field selection"""
        dst = "ngc5921.split.fieldsel.ms"
        combine = "corr"
        [origwt, origwtsp, origflag, origfrow, origdata] = _get_dst_cols(src)
        rtol = 1e-7
        for field in ["2", "N5921_2"]:
            shutil.copytree(src, dst)
            statwt(dst, field=field, combine=combine)
            [wt, wtsp, flag, frow, data, field_id] = _get_dst_cols(dst, "FIELD_ID")
            nrow = len(frow)
            dr = numpy.real(data)
            di = numpy.imag(data)
            for row in range(nrow):
                if field_id[row] == 2:
                    expec = _variance(dr, di, flag, row)
                    self.assertTrue(
                        numpy.all(numpy.isclose(wt[:, row], expec, rtol=rtol)),
                        "WEIGHT fail at row" + str(row) + ". got: "
                        + str(wt[:, row]) + " expec " + str(expec)
                    )
                    self.assertTrue(
                        numpy.all(numpy.isclose(wtsp[:,:,row], expec, rtol)),
                        "Incorrect weight spectrum"   
                    )
                else:
                    self.assertTrue(
                        numpy.all(numpy.isclose(wt[:, row], origwt[:, row], rtol=rtol)),
                        "WEIGHT fail at row" + str(row) + ". got: " + str(wt[:, row])
                        + " expec " + str(origwt[:, row])
                    )
                    self.assertTrue(
                        numpy.all(numpy.isclose(wtsp[:,:,row], origwtsp[:,:,row], rtol)),
                        "Incorrect weight spectrum"   
                    )
            shutil.rmtree(dst)
            
    def test_spwsel(self):
        """Test spw selection"""
        dst = "ngc5921.split.spwsel.ms"
        combine = "corr"
        [origwt, origwtsp, origflag, origfrow, origdata] = _get_dst_cols(src)
        rtol = 1e-7
        spw="0"
        # data set only has one spw
        shutil.copytree(src, dst)
        statwt(dst, spw=spw, combine=combine)
        [wt, wtsp, flag, frow, data] = _get_dst_cols(dst)
        nrow = len(frow)
        dr = numpy.real(data)
        di = numpy.imag(data)
        for row in range(nrow):
            expec = _variance(dr, di, flag, row)
            self.assertTrue(
                numpy.all(numpy.isclose(wt[:, row], expec, rtol=rtol)),
                "WEIGHT fail at row" + str(row) + ". got: "
                + str(wt[:, row]) + " expec " + str(expec)
            )
            self.assertTrue(
                numpy.all(numpy.isclose(wtsp[:,:,row], expec, rtol)),
                "Incorrect weight spectrum"   
            )
        shutil.rmtree(dst)

    def test_scansel(self):
        """CAS-11858 Test scan selection"""
        dst = "ngc5921.split.scansel.ms"
        combine = "corr"
        [origwt, origwtsp, origflag, origfrow, origdata] = _get_dst_cols(src)
        rtol = 1e-7
        scan = "5"
        shutil.copytree(src, dst)
        statwt(dst, scan=scan, combine=combine)
        [wt, wtsp, flag, frow, data, scan_id] = _get_dst_cols(dst, "SCAN_NUMBER")
        nrow = len(frow)
        dr = numpy.real(data)
        di = numpy.imag(data)
        for row in range(nrow):
            if str(scan_id[row]) == scan:
                expec = _variance(dr, di, flag, row)
                self.assertTrue(
                    numpy.all(numpy.isclose(wt[:, row], expec, rtol=rtol)),
                    "WEIGHT fail at row" + str(row) + ". got: "
                    + str(wt[:, row]) + " expec " + str(expec)
                )
                self.assertTrue(
                    numpy.all(numpy.isclose(wtsp[:,:,row], expec, rtol)),
                    "Incorrect weight spectrum"   
                )
            else:
                self.assertTrue(
                    numpy.all(numpy.isclose(wt[:, row], origwt[:, row], rtol=rtol)),
                    "WEIGHT fail at row" + str(row) + ". got: " + str(wt[:, row])
                    + " expec " + str(origwt[:, row])
                )
                self.assertTrue(
                    numpy.all(numpy.isclose(wtsp[:,:,row], origwtsp[:,:,row], rtol)),
                    "Incorrect weight spectrum"   
                )
        shutil.rmtree(dst)
        
    def test_default_boundaries(self):
        """Test default scan, field, etc boundaries"""
        dst = "ngc5921.split.normalbounds.ms"
        timebin = "6000s"
        ref = os.path.join(datadir,"ngc5921.normal_bounds_2.ms.ref")
        rtol = 1e-7
        [expwt, expwtsp, expflag, expfrow, expdata] = _get_dst_cols(ref)
        # there are three field_ids, and there is a change in field_id when
        # there is a change in scan number, so specifying combine="field" in the
        # absence of "scan" will give the same result as combine=""
        for combine in ["corr", "corr,field"]:
            shutil.copytree(src, dst)
            statwt(dst, timebin=timebin, combine=combine)
            [gotwt, gotwtsp, gotflag, gotfrow, gotdata] = _get_dst_cols(dst)
            self.assertTrue(numpy.all(numpy.isclose(gotwt, expwt, rtol)))
            self.assertTrue(numpy.all(numpy.isclose(gotwtsp, expwtsp, rtol)))
            self.assertTrue(numpy.all(gotflag == expflag))
            self.assertTrue(numpy.all(gotfrow == expfrow))
            shutil.rmtree(dst)
        
    def test_no_scan_boundaries(self):
        """Test no scan boundaries"""
        dst = "ngc5921.no_scan_bounds.ms"
        timebin = "6000s"
        ref = os.path.join(datadir,"ngc5921.no_scan_bounds_2.ms.ref")
        rtol = 1e-7
        [expwt, expwtsp, expflag, expfrow, expdata] = _get_dst_cols(ref)
        combine = "corr, scan"
        shutil.copytree(src, dst)
        statwt(dst, timebin=timebin, combine=combine)
        [gotwt, gotwtsp, gotflag, gotfrow, gotdata] = _get_dst_cols(dst)
        self.assertTrue(numpy.all(numpy.isclose(gotwt, expwt, rtol)))
        self.assertTrue(numpy.all(numpy.isclose(gotwtsp, expwtsp, rtol)))
        self.assertTrue(numpy.all(gotflag == expflag))
        self.assertTrue(numpy.all(gotfrow == expfrow))
        shutil.rmtree(dst)
    
    def test_no_scan_nor_field_boundaries(self):
        """Test no scan nor field boundaries"""
        dst = "ngc5921.no_scan_nor_field_bounds.ms"
        timebin = "6000s"
        ref = os.path.join(datadir,"ngc5921.no_scan_nor_field_bounds_2.ms.ref")
        rtol = 1e-7
        [expwt, expwtsp, expflag, expfrow, expdata] = _get_dst_cols(ref)
        for combine in ["corr,scan,field", "corr,field,scan"]:
            shutil.copytree(src, dst)
            statwt(dst, timebin=timebin, combine=combine)
            [gotwt, gotwtsp, gotflag, gotfrow, gotdata] = _get_dst_cols(dst)
            self.assertTrue(numpy.all(numpy.isclose(gotwt, expwt, rtol)))
            self.assertTrue(numpy.all(numpy.isclose(gotwtsp, expwtsp, rtol)))
            self.assertTrue(numpy.all(gotflag == expflag))
            self.assertTrue(numpy.all(gotfrow == expfrow))
            shutil.rmtree(dst)
                
    def test_statalg(self):
        """Test statalg"""
        # just testing inputs
        dst = "ngc5921.split.statalg.ms"
        for statalg in ["cl", "ch", "h", "f", "bogus"]:
            shutil.copytree(src, dst)
            if statalg == "cl":
                statwt(vis=dst, statalg=statalg)
            elif statalg == "ch":
                self.assertTrue(statwt(vis=dst, statalg=statalg, zscore=5, maxiter=3))
            elif statalg == "h":
                self.assertTrue(statwt(vis=dst, statalg=statalg, fence=0.2))
            elif statalg == "f":
                self.assertTrue(statwt(vis=dst, statalg=statalg, center="median", lside=False))
            elif statalg == "bogus":
                self.assertRaises(Exception, statwt, vis=dst, statalg=statalg)
            shutil.rmtree(dst)
                
    def test_wtrange(self):
        """ Test weight range"""
        dst = "ngc5921.split.timebin.ms"
        ref = os.path.join(datadir,"ngc5921.timebin300s_2.ms.ref")
        [refwt, refwtsp, refflag, reffrow, refdata] = _get_dst_cols(ref)
        rtol = 1e-7
        combine = "corr"
        timebin = "300s"
        wtrange = [1, 2]
        shutil.copytree(src, dst) 
        statwt(dst, timebin=timebin, combine=combine, wtrange=wtrange)
        [tstwt, tstwtsp, tstflag, tstfrow, tstdata] = _get_dst_cols(dst)
        self.assertTrue(
            numpy.all(
                tstflag == numpy.logical_or(
                    refflag, numpy.logical_not(
                        numpy.logical_and(tstwtsp >= 1, tstwtsp <= 2)
                    )
                )
            ),
            "FLAGs don't match"
        )
        self.assertTrue(
            numpy.all(numpy.all(tstflag, axis=(0,1)) == tstfrow),
            "FLAG_ROWs don't match"
        )
        nrows = tstwtsp.shape[2]
        for row in range(nrows):
            rowwtsp = tstwtsp[:,:,row][numpy.logical_not(tstflag[:,:,row])]
            if (len(rowwtsp) == 0):
                expec = 0
            else:
                expec = rowwtsp[0]
            self.assertTrue(
                numpy.all(numpy.isclose(tstwt[:, row], expec, rtol)),
                "WEIGHTs don't match"
            )
        self.assertTrue(
            numpy.all(numpy.isclose(tstwtsp, refwtsp, rtol)),
            "WEIGHT_SPECTRUMs don't match"
        )
        shutil.rmtree(dst)

    def test_preview(self):
        """ Test preview mode"""
        dst = "ngc5921.split.preview.ms"
        [refwt, refwtsp, refflag, reffrow, refdata] = _get_dst_cols(src)
        rtol = 1e-7
        combine = "corr"
        timebin = "300s"
        wtrange = [1, 2]
        preview = True
        shutil.copytree(src, dst)
        statwt(
            dst, timebin=timebin, combine=combine, wtrange=wtrange, preview=preview
        )
        [tstwt, tstwtsp, tstflag, tstfrow, tstdata] = _get_dst_cols(dst)
        self.assertTrue(numpy.all(tstflag == refflag), "FLAGs don't match")
        self.assertTrue(numpy.all(tstfrow == reffrow), "FLAG_ROWs don't match")
        self.assertTrue(
            numpy.all(numpy.isclose(tstwt, refwt, rtol)), "WEIGHTs don't match"
        )
        self.assertTrue(
            numpy.all(numpy.isclose(tstwtsp, refwtsp, rtol)),
            "WEIGHT_SPECTRUMs don't match"
        )
        shutil.rmtree(dst)

    def test_data(self):
        """ Test using data column"""
        dst = "ngc5921.split.data.ms"
        ref = os.path.join(datadir,"ngc5921.data_col.ms.ref")
        [refwt, refwtsp, refflag, reffrow, refsig, refsigsp] = _get_dst_cols(
            ref, ["SIGMA", "SIGMA_SPECTRUM"], dodata=False
        )
        rtol = 1e-7
        combine = "corr"
        timebin = 10
        data = "data"
        mytb = table()
        shutil.copytree(src, dst)
        self.assertTrue(mytb.open(dst, nomodify=False))
        self.assertTrue(mytb.removecols("DATA"))
        self.assertTrue(mytb.renamecol("CORRECTED_DATA", "DATA"))
        mytb.done()
        statwt(dst, timebin=timebin, combine=combine, datacolumn=data)
        [tstwt, tstwtsp, tstflag, tstfrow, tstsig, tstsigsp] = _get_dst_cols(
            dst, ["SIGMA", "SIGMA_SPECTRUM"], False
        )
        self.assertTrue(numpy.all(tstflag == refflag), "FLAGs don't match")
        self.assertTrue(numpy.all(tstfrow == reffrow), "FLAG_ROWs don't match")
        self.assertTrue(
            numpy.all(numpy.isclose(tstwt, refwt, rtol)),
            "WEIGHTs don't match"
        )
        self.assertTrue(
            numpy.all(numpy.isclose(tstwtsp, refwtsp, rtol)),
            "WEIGHT_SPECTRUMs don't match"
        )
        self.assertTrue(
            numpy.all(numpy.isclose(tstsig, refsig)),
            "SIGMA is incorrect"
        )
        self.assertTrue(
            numpy.all(numpy.isclose(tstsigsp, refsigsp)),
            "SIGMA_SPECTRUM is incorrect"
        )
        shutil.rmtree(dst)

    def test_slding_time_window(self):
        """ Test sliding time window"""
        dst = "ngc5921.split.sliding_time_window.ms"
        ref = os.path.join(datadir,"ngc5921.slidingtimebin300s_2.ms.ref")
        [refwt, refwtsp, refflag, reffrow] = _get_dst_cols(ref, "", dodata=False)
        rtol = 1e-7
        timebin = "300s"
        shutil.copytree(src, dst)
        statwt(dst, timebin=timebin, slidetimebin=True)
        [tstwt, tstwtsp, tstflag, tstfrow] = _get_dst_cols(dst, "", False)
        self.assertTrue(numpy.all(tstflag == refflag), "FLAGs don't match")
        self.assertTrue(numpy.all(tstfrow == reffrow), "FLAG_ROWs don't match")
        self.assertTrue(
            numpy.all(numpy.isclose(tstwt, refwt, rtol)),
            "WEIGHTs don't match"
        )
        self.assertTrue(
            numpy.all(numpy.isclose(tstwtsp, refwtsp, rtol)),
            "WEIGHT_SPECTRUMs don't match"
        )
        shutil.rmtree(dst)

    def test_residual(self):
        """ Test using corrected_data - model_data column"""
        dst = "ngc5921.split.residualwmodel.ms"
        ref = os.path.join(datadir,"ngc5921.resid_with_model.ms.ref")
        [refwt, refwtsp, refflag, reffrow] = _get_dst_cols(ref, "", dodata=False)
        rtol = 1e-7
        data = "residual"
        mytb = table()
        shutil.copytree(src, dst)
        statwt(dst, datacolumn=data)
        [tstwt, tstwtsp, tstflag, tstfrow] = _get_dst_cols(dst, "", False)
        self.assertTrue(numpy.all(tstflag == refflag), "FLAGs don't match")
        self.assertTrue(numpy.all(tstfrow == reffrow), "FLAG_ROWs don't match")
        self.assertTrue(
            numpy.all(numpy.isclose(tstwt, refwt, rtol)),
            "WEIGHTs don't match"
        )
        self.assertTrue(
            numpy.all(numpy.isclose(tstwtsp, refwtsp, rtol)),
            "WEIGHT_SPECTRUMs don't match"
        )
        shutil.rmtree(dst)
            
    def test_residual_no_model(self):
        """ Test using corrected_data - model_data column"""
        dst = "ngc5921.split.residualwoutmodel.ms"
        ref = os.path.join(datadir,"ngc5921.resid_without_model.ms.ref")
        [refwt, refwtsp, refflag, reffrow] = _get_dst_cols(ref, "", dodata=False)
        rtol = 1e-6
        data = "residual"
        mytb = table()
        shutil.copytree(src, dst)
        self.assertTrue(mytb.open(dst, nomodify=False))
        self.assertTrue(mytb.removecols("MODEL_DATA"))
        mytb.done()
        statwt(dst, datacolumn=data)
        [tstwt, tstwtsp, tstflag, tstfrow] = _get_dst_cols(dst, "", False)
        self.assertTrue(numpy.all(tstflag == refflag), "FLAGs don't match")
        self.assertTrue(numpy.all(tstfrow == reffrow), "FLAG_ROWs don't match")
        self.assertTrue(
            numpy.all(numpy.isclose(tstwt, refwt, rtol)),
            "WEIGHTs don't match"
        )
        self.assertTrue(
            numpy.all(numpy.isclose(tstwtsp, refwtsp, rtol)),
            "WEIGHT_SPECTRUMs don't match"
        )
        shutil.rmtree(dst)

    def test_residual_data(self):
        """ Test using _data - model_data column"""
        dst = "ngc5921.split.residualdatawmodel.ms"
        ref = os.path.join(datadir,"ngc5921.residdata_with_model_2.ms.ref")
        [refwt, refwtsp, refflag, reffrow, refsig, refsigsp] = _get_dst_cols(
            ref, ["SIGMA", "SIGMA_SPECTRUM"], dodata=False
        )
        rtol = 1e-7
        data = "residual_data"
        mytb = table()
        shutil.copytree(src, dst)
        statwt(dst, datacolumn=data)
        [tstwt, tstwtsp, tstflag, tstfrow, tstsig, tstsigsp] = _get_dst_cols(
            dst, ["SIGMA", "SIGMA_SPECTRUM"], False
        )
        self.assertTrue(numpy.all(tstflag == refflag), "FLAGs don't match")
        self.assertTrue(numpy.all(tstfrow == reffrow), "FLAG_ROWs don't match")
        shutil.rmtree(dst)
        self.assertTrue(
            numpy.all(numpy.isclose(tstwt, refwt, rtol)),
            "WEIGHTs don't match"
        )
        self.assertTrue(
            numpy.all(numpy.isclose(tstwtsp, refwtsp, rtol)),
            "WEIGHT_SPECTRUMs don't match"
        )
        self.assertTrue(
            numpy.all(numpy.isclose(tstsig, refsig, rtol)),
            "SIGMAs don't match"
        )
        self.assertTrue(
            numpy.all(numpy.isclose(tstsigsp, refsigsp, rtol)),
            "SIGMA_SPECTRUMs don't match"
        )

    def test_residual_data_no_model(self):
        """ Test using data - default model """
        dst = "ngc5921.split.residualdatawoutmodel.ms"
        ref = os.path.join(datadir,"ngc5921.residdata_without_model_2.ms.ref")
        [refwt, refwtsp, refflag, reffrow, refsig, refsigsp] = _get_dst_cols(
            ref, ["SIGMA", "SIGMA_SPECTRUM"], dodata=False
        )
        rtol = 1e-7
        data = "residual_data"
        mytb = table()
        shutil.copytree(src, dst)
        self.assertTrue(mytb.open(dst, nomodify=False))
        self.assertTrue(mytb.removecols("MODEL_DATA"))
        mytb.done()
        statwt(dst, datacolumn=data)
        [tstwt, tstwtsp, tstflag, tstfrow, tstsigma, tstsigsp] = _get_dst_cols(
            dst, ["SIGMA", "SIGMA_SPECTRUM"], False
        )
        self.assertTrue(numpy.all(tstflag == refflag), "FLAGs don't match")
        self.assertTrue(numpy.all(tstfrow == reffrow), "FLAG_ROWs don't match")
        refsigma = 1/numpy.sqrt(refwt);
        numpy.place(refsigma, refwt == 0, -1)
        self.assertTrue(
            numpy.all(numpy.isclose(tstwt, refwt, rtol)),
            "WEIGHTs don't match"
        )
        self.assertTrue(
            numpy.all(numpy.isclose(tstwtsp, refwtsp, rtol)),
            "WEIGHT_SPECTRUMs don't match"
        )
        self.assertTrue(
            numpy.all(numpy.isclose(tstsigma, refsig, rtol)),
            "SIGMAs don't match"
        )
        self.assertTrue(
            numpy.all(numpy.isclose(tstsigsp, refsigsp, rtol)),
            "SIGMA_SPECTRUMs don't match"
        )
        shutil.rmtree(dst)

    def test_returned_stats(self):
        """ Test returned stats, CAS-10881"""
        dst = "ngc5921.split.statstest.ms"
        ref = os.path.join(datadir,"ngc5921.residdata_without_model_2.ms.ref")
        rtol = 1e-7
        shutil.copytree(src, dst)
        res = statwt(dst)
        self.assertTrue(
            numpy.isclose(res['mean'], 3.6326332, rtol),
            "mean is incorrect"
        )
        self.assertTrue(
            numpy.isclose(res['variance'], 6.6448922, rtol),
            "variance is incorrect"
        )
        shutil.rmtree(dst)

def suite():
    return [statwt_test]

if __name__ == '__main__':
    unittest.main()
