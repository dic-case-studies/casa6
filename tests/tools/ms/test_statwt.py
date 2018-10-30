import os
import sys
import shutil
import unittest
import math
import numpy
import numbers

from casatools import ms as mstool
from casatools import table as tbtool
from casatools import ctsys

datadir = 'regression/unittest/statwt/'
src = datadir + 'ngc5921.split_2.ms'

def _get_dst_cols(dst, other="", dodata=True):
    mytb = tbtool()
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


    def test_algorithm(self):
        """ Test the algorithm, includes fitspw tests"""
        mytb = table()
        mytb.open(src)
        expflag = mytb.getcol("FLAG")
        expfrow = mytb.getcol("FLAG_ROW")
        mytb.done()
        dst = "ngc5921.split.ms"
        rtol = 1e-7
        for combine in ["", "corr"]:
            for fitspw in ["0:0~9;21~62", ""]:
                for i in [0]:
                    shutil.copytree(src, dst) 
                    myms = ms( )
                    if i == 0:
                        myms.open(dst, nomodify=False)
                        myms.statwt(combine=combine, fitspw=fitspw)
                        myms.done()
                    else:
                        statwt(dst, combine=combine, fitspw=fitspw)
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
               
    def test_timebin(self):
        """ Test time binning"""
        dst = "ngc5921.split.timebin.ms"
        ref = datadir + "ngc5921.timebin300s_2.ms.ref"
        [refwt, refwtsp, refflag, reffrow, refdata] = _get_dst_cols(ref)
        rtol = 1e-7
        combine = "corr"
        for timebin in ["300s", 10]:
            shutil.copytree(ctsys.resolve(src), dst) 
            myms = mstool()
            myms.open(dst, nomodify=False)
            myms.statwt(timebin=timebin, combine=combine)
            myms.done()
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
                ref = datadir + "ngc5921.chanbin_sepcorr_2.ms.ref"
            else:
                ref = datadir + "ngc5921.chanbin_combcorr_2.ms.ref"
            [refwt, refwtsp, refflag, reffrow, refdata] = _get_dst_cols(ref)
            for i in [0, 2]:
                for chanbin in ["195.312kHz", 8]:
                    if i == 2 and (combine == "corr" or chanbin == 8):
                        continue
                    shutil.copytree(ctsys.resolve(src), dst)
                    if i == 0:
                        myms = mstool()
                        myms.open(dst, nomodify=False)
                        myms.statwt(chanbin=chanbin, combine=combine)
                        myms.done()
                    elif i == 2:
                        # check WEIGHT_SPECTRUM is created, only check once,
                        # this test is long as it is
                        # shutil.copytree(src, dst)
                        mytb = tbtool()
                        mytb.open(dst, nomodify=False)
                        x = mytb.ncols()
                        self.assertTrue(mytb.removecols("WEIGHT_SPECTRUM"), "column not removed")
                        y = mytb.ncols()
                        self.assertTrue(y == x-1, "wrong number of columns")
                        mytb.done()
                        myms = mstool()
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
            myms = mstool()
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

    @unittest.skip("uses a task")
    def test_fieldsel(self):
        """Test field selection"""
        dst = "ngc5921.split.fieldsel.ms"
        combine = "corr"
        [origwt, origwtsp, origflag, origfrow, origdata] = _get_dst_cols(src)
        rtol = 1e-7
        for field in ["2", "N5921_2"]:
            shutil.copytree(ctsys.resolve(src), dst)
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
            
    @unittest.skip("uses a task")
    def test_spwsel(self):
        """Test spw selection"""
        dst = "ngc5921.split.spwsel.ms"
        combine = "corr"
        [origwt, origwtsp, origflag, origfrow, origdata] = _get_dst_cols(src)
        rtol = 1e-7
        spw="0"
        # data set only has one spw
        shutil.copytree(ctsys.resolve(src), dst)
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
        
    def test_default_boundaries(self):
        """Test default scan, field, etc boundaries"""
        dst = "ngc5921.split.normalbounds.ms"
        timebin = "6000s"
        ref = datadir + "ngc5921.normal_bounds_2.ms.ref"
        rtol = 1e-7
        [expwt, expwtsp, expflag, expfrow, expdata] = _get_dst_cols(ref)
        # there are three field_ids, and there is a change in field_id when
        # there is a change in scan number, so specifying combine="field" in the
        # absence of "scan" will give the same result as combine=""
        for combine in ["corr", "corr,field"]:
            shutil.copytree(ctsys.resolve(src), dst)
            myms = mstool()
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
        ref = datadir + "ngc5921.no_scan_bounds_2.ms.ref"
        rtol = 1e-7
        [expwt, expwtsp, expflag, expfrow, expdata] = _get_dst_cols(ref)
        combine = "corr, scan"
        shutil.copytree(ctsys.resolve(src), dst)
        myms = mstool()
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
        ref = datadir + "ngc5921.no_scan_nor_field_bounds_2.ms.ref"
        rtol = 1e-7
        [expwt, expwtsp, expflag, expfrow, expdata] = _get_dst_cols(ref)
        for combine in ["corr,scan,field", "corr,field,scan"]:
            shutil.copytree(ctsys.resolve(src), dst)
            myms = mstool()
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
            myms = mstool()
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
        ref = datadir + "ngc5921.timebin300s_2.ms.ref"
        [refwt, refwtsp, refflag, reffrow, refdata] = _get_dst_cols(ref)
        rtol = 1e-7
        combine = "corr"
        timebin = "300s"
        wtrange = [1, 2]
        shutil.copytree(ctsys.resolve(src), dst) 
        myms = mstool()
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
        myms = mstool()
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
        ref = datadir + "ngc5921.timebin300s_2.ms.ref"
        [refwt, refwtsp, refflag, reffrow] = _get_dst_cols(ref, "", dodata=False)
        rtol = 1e-7
        combine = "corr"
        timebin = 10
        data = "data"
        mytb = tbtool()
        myms = mstool()
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
        ref = datadir + "ngc5921.slidingtimebin300s_2.ms.ref"
        [refwt, refwtsp, refflag, reffrow] = _get_dst_cols(ref, "", dodata=False)
        rtol = 1e-7
        timebin = "300s"
        myms = mstool()
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
