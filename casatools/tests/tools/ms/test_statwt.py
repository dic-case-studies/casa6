import os
import sys
import shutil
import unittest
import math
import numpy as np
import numpy.ma as ma
import numbers

from casatools import ctsys, table, ms

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
        othercol = mytb.getcol(other)
    mytb.close()
    cols = [wt, wtsp, flag, frow]
    if dodata:
        cols.append(data)
    if len(other) > 0:
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
        self, msname, row_to_rows, data_column, chan_flags, combine_corr
    ):
        if data_column.startswith('c'):
            colname = 'CORRECTED_DATA'
        elif data_column.startswith('d'):
            colname = 'DATA'
        else:
            raise Exception("Unhandled column spec " + data_column)
        # if not mode.startswith('one'):
        #    raise Exception("Unhandled mode")
        for ant1 in range(10):
            for ant2 in range((ant1 + 1), 10):
                query_str = 'ANTENNA1=' + str(ant1) + ' AND ANTENNA2=' + str(ant2)
                tb = table()
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
                    #if mode.startswith('one'):
                    #    start = row
                    #    end = row+1
                    #if mode.startswith('one'):
                    #    start = row
                    #    end = row+1
                    # print 'baseline ' + str([ant1, ant2]) + ' row ' + str(row)
                    start = row_to_rows[row][0]
                    end = row_to_rows[row][1]
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

    
    """
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
        """
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
        myms = ms()
        row_to_rows = []
        for row in range(60):
            row_to_rows.append((row, row+1))
        for combine in ["", "corr"]:
            c = 0
            for fitspw in ["0:0~9;21~62", "", "0:10~20"]:
                shutil.copytree(src, dst) 
                excludechans = c == 2
                # tool method
                myms.open(dst, nomodify=False)
                myms.statwt(
                    combine=combine, fitspw=fitspw,
                    excludechans=excludechans
                )
                myms.done()
                chan_flags = cflags if fitspw else None
                self._check_weights(
                    dst, row_to_rows=row_to_rows, data_column='c',
                    chan_flags=chan_flags, combine_corr=bool(combine)
                )
                """
                self._check_weights(
                    dst, mode='one_to_one', data_column='c',
                    chan_flags=chan_flags, combine_corr=bool(combine)
                )
                """
                shutil.rmtree(dst)
                c += 1               


    """
    def test_algorithm(self):
        "" Test the algorithm, includes fitspw, excludechan tests""
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
                myms = ms( )
                myms.open(dst, nomodify=False)
                myms.statwt(combine=combine, fitspw=fitspw, excludechans=excludechans)
                myms.done()
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
        # ref = datadir + "ngc5921.timebin300s_2.ms.ref"
        # [refwt, refwtsp, refflag, reffrow, refdata] = _get_dst_cols(ref)
        # rtol = 1e-7
        combine = "corr"
        r2r_300 = []
        for i in range(10):
            r2r_300.append((0, 10))
        for i in range(10, 12):
            r2r_300.append((10, 12))
        for i in range(12, 17):
            r2r_300.append((12, 17))
        for i in range(17, 22):
            r2r_300.append((17, 22))
        for i in range(22, 27):
            r2r_300.append((22, 27))
        for i in range(27, 32):
            r2r_300.append((27, 32))
        r2r_300.append((32, 33))
        for i in range(33, 35):
            r2r_300.append((33, 35))
        for i in range(35, 38):
            r2r_300.append((35, 38))
        for i in range(38, 43):
            r2r_300.append((38, 43))
        for i in range(43, 48):
            r2r_300.append((43, 48))
        for i in range(48, 53):
            r2r_300.append((48, 53))
        for i in range(53,56):
            r2r_300.append((53, 56))
        for i in range(56, 60):
            r2r_300.append((56, 60))
            
        r2r_10 = []
        for i in range(10):
            r2r_10.append((0, 10))
        for i in range(10, 12):
            r2r_10.append((2, 12))
        for i in range(12, 17):
            r2r_10.append((12, 17))
        for i in range(17, 27):
            r2r_10.append((17, 27))
        for i in range(27, 33):
            r2r_10.append((23, 33))
        for i in range(33, 35):
            r2r_10.append((33, 35))
        for i in range(35, 38):
            r2r_10.append((35, 38))
        for i in range(38, 48):
            r2r_10.append((38, 48))
        for i in range(48, 56):
            r2r_10.append((46, 56))
        for i in range(56, 60):
            r2r_10.append((56, 60))
        for timebin in ["300s", 10]:
            shutil.copytree(src, dst) 
            myms = ms()
            myms.open(dst, nomodify=False)
            myms.statwt(timebin=timebin, combine=combine)
            myms.done()
            if timebin == "300s":
                row_to_rows = r2r_300
            else:
                row_to_rows = r2r_10
            self._check_weights(
                dst, row_to_rows=row_to_rows, data_column='c',
                chan_flags=None, combine_corr=True
            )
            
            
            # [tstwt, tstwtsp, tstflag, tstfrow, tstdata] = _get_dst_cols(dst)
            # self.assertTrue(np.all(tstflag == refflag), "FLAGs don't match")
            # self.assertTrue(np.all(tstfrow == reffrow), "FLAG_ROWs don't match")
            # self.assertTrue(np.all(np.isclose(tstwt, refwt, rtol)), "WEIGHTs don't match")
            # self.assertTrue(np.all(np.isclose(tstwtsp, refwtsp, rtol)), "WEIGHT_SPECTRUMs don't match")
            
            shutil.rmtree(dst)

    
    """
    def test_timebin(self):
        "" Test time binning""
        dst = "ngc5921.split.timebin.ms"
        ref = os.path.join(datadir,"ngc5921.timebin300s_2.ms.ref")
        [refwt, refwtsp, refflag, reffrow, refdata] = _get_dst_cols(ref)
        rtol = 1e-7
        combine = "corr"
        for timebin in ["300s", 10]:
            shutil.copytree(ctsys.resolve(src), dst) 
            myms = ms( )
            myms.open(dst, nomodify=False)
            myms.statwt(timebin=timebin, combine=combine)
            myms.done()
            [tstwt, tstwtsp, tstflag, tstfrow, tstdata] = _get_dst_cols(dst)
            self.assertTrue(numpy.all(tstflag == refflag), "FLAGs don't match")
            self.assertTrue(numpy.all(tstfrow == reffrow), "FLAG_ROWs don't match")
            self.assertTrue(numpy.all(numpy.isclose(tstwt, refwt, rtol)), "WEIGHTs don't match")
            self.assertTrue(numpy.all(numpy.isclose(tstwtsp, refwtsp, rtol)), "WEIGHT_SPECTRUMs don't match")
            shutil.rmtree(dst)
    """
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
            for i in [0, 2]:
                for chanbin in ["195.312kHz", 8]:
                    if i == 2 and (combine == "corr" or chanbin == 8):
                        continue
                    shutil.copytree(ctsys.resolve(src), dst)
                    if i == 0:
                        myms = ms( )
                        myms.open(dst, nomodify=False)
                        myms.statwt(chanbin=chanbin, combine=combine)
                        myms.done()
                    elif i == 2:
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
                        myms = ms( )
                        myms.open(dst, nomodify=False)
                        myms.statwt(chanbin=chanbin, combine=combine)
                        myms.done()
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
            shutil.copytree(ctsys.resolve(src), dst)
            myms = ms( )
            myms.open(dst, nomodify=False)
            myms.statwt(minsamp=minsamp, combine=combine)
            myms.done()
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
            shutil.copytree(ctsys.resolve(src), dst)
            myms = ms( )
            myms.open(dst, nomodify=False)
            myms.statwt(timebin=timebin, combine=combine)
            myms.done()
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
        shutil.copytree(ctsys.resolve(src), dst)
        myms = ms( )
        myms.open(dst, nomodify=False)
        myms.statwt(timebin=timebin, combine=combine)
        myms.done()
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
            shutil.copytree(ctsys.resolve(src), dst)
            myms = ms( )
            myms.open(dst, nomodify=False)
            myms.statwt(timebin=timebin, combine=combine)
            myms.done()
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
            shutil.copytree(ctsys.resolve(src), dst)
            myms = ms( )
            myms.open(dst, nomodify=False)
            if statalg == "cl":
                self.assertTrue(myms.statwt(statalg=statalg))
            elif statalg == "ch":
                self.assertTrue(myms.statwt(statalg=statalg, zscore=5, maxiter=3))
            elif statalg == "h":
                self.assertTrue(myms.statwt(statalg=statalg, fence=0.2))
            elif statalg == "f":
                self.assertTrue(myms.statwt(statalg=statalg, center="median", lside=False))
            elif statalg == "bogus":
                self.assertRaises(Exception, myms.statwt, statalg=statalg)
            myms.done()
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
        shutil.copytree(ctsys.resolve(src), dst) 
        myms = ms( )
        myms.open(dst, nomodify=False)
        myms.statwt(timebin=timebin, combine=combine, wtrange=wtrange)
        myms.done()
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
        shutil.copytree(ctsys.resolve(src), dst)
        myms = ms( )
        myms.open(dst, nomodify=False)
        myms.statwt(
            timebin=timebin, combine=combine, wtrange=wtrange, preview=preview
        )
        myms.done()
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
        ref = os.path.join(datadir,"ngc5921.timebin300s_2.ms.ref")
        [refwt, refwtsp, refflag, reffrow] = _get_dst_cols(ref, "", dodata=False)
        rtol = 1e-7
        combine = "corr"
        timebin = 10
        data = "data"
        mytb = table()
        myms = ms( )
        shutil.copytree(ctsys.resolve(src), dst)
        self.assertTrue(mytb.open(dst, nomodify=False))
        self.assertTrue(mytb.removecols("DATA"))
        self.assertTrue(mytb.renamecol("CORRECTED_DATA", "DATA"))
        mytb.done()
        myms.open(dst, nomodify=False)
        myms.statwt(timebin=timebin, combine=combine, datacolumn=data)
        myms.done()
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

    def test_slding_time_window(self):
        """ Test sliding time window"""
        dst = "ngc5921.split.sliding_time_window.ms"
        ref = os.path.join(datadir,"ngc5921.slidingtimebin300s_2.ms.ref")
        [refwt, refwtsp, refflag, reffrow] = _get_dst_cols(ref, "", dodata=False)
        rtol = 1e-7
        timebin = "300s"
        myms = ms( )
        shutil.copytree(ctsys.resolve(src), dst)
        myms.open(dst, nomodify=False)
        myms.statwt(timebin=timebin, slidetimebin=True)
        myms.done()
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

def suite():
    return [statwt_test]

if __name__ == '__main__':
    unittest.main()
