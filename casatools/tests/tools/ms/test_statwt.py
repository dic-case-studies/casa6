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

# rows and target_row are the row numbers from the subtable formed
# by the baseline query
# In the chan_flags, a value of False means the channel is good (not flagged)
# so should be used. It follows the convention of the FLAGS column in the MS.

# EVEN IF THIS IS NO LONGER USED BY THE TESTS, IT SHOULDN'T BE DELETED BECAUSE
# IT IS USEFUL IN SANTIFY CHECKING NEW TESTS
def get_weights(
    data, flags, chan_flags, exposures, combine_corr, target_exposure, chanbins,
    target_flags, wtrange
):  
    shape = data.shape
    ncorr_groups = 1 if combine_corr else shape[0]
    nchanbins = 1 if chanbins is None else len(chanbins)
    tchanbins = chanbins
    if nchanbins == 1:
        tchanbins = [[0, shape[1]]]
    ncorr = shape[0]
    weights = np.zeros([shape[0], shape[1]])
    wt = np.zeros(shape[0])
    nrows = data.shape[2]
    median_axis = 1 if ncorr_groups > 1 else None
    mod_flags = target_flags[:]
    if type(chan_flags) == type(None):
        myflags = flags[:]
    else:
        t_flags = np.expand_dims(np.expand_dims(chan_flags, 0), 2)
        myflags = np.logical_or(flags, t_flags)
    for corr in range(ncorr_groups):
        end_corr = corr + 1 if ncorr_groups > 1 else ncorr + 1
        for cb in tchanbins:
            var = variance(
                data[corr:end_corr, cb[0]:cb[1], :],
                myflags[corr:end_corr, cb[0]:cb[1], :], exposures
            )
            if flags[corr:end_corr, cb[0]:cb[1]].all():
                weights[corr:end_corr, cb[0]:cb[1]] = 0
                mod_flags[corr:end_corr, cb[0]:cb[1]] = True
            if var == 0:
                weights[corr:end_corr, cb[0]:cb[1]] = 0 
                mod_flags[corr:end_corr, cb[0]:cb[1]] = True
            else:
                weights[corr:end_corr, cb[0]:cb[1]] = target_exposure/var
            if type(wtrange) != type(None):
                condition = np.logical_or(
                    np.less(
                        weights[corr:end_corr, cb[0]:cb[1]], wtrange[0]
                    ),
                    np.greater(
                        weights[corr:end_corr, cb[0]:cb[1]], wtrange[1]
                    )
                )
                exp_condition = np.expand_dims(condition, 2)
                weights[corr:end_corr, cb[0]:cb[1]] = np.where(
                    condition, 0, weights[corr:end_corr, cb[0]:cb[1]]
                )
                mod_flags[corr:end_corr, cb[0]:cb[1]] = np.where(
                    exp_condition, True, mod_flags[corr:end_corr, cb[0]:cb[1]]
                )
            mweights = ma.array(
                weights[corr:end_corr, :],
                mask=mod_flags[corr:end_corr, :]
            )
            wt[corr:end_corr] = np.median(mweights, median_axis)
    mod_flags = np.where(np.expand_dims(weights, 2) == 0, True, mod_flags)
    return (weights, wt, mod_flags)

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

class statwt_test(unittest.TestCase):
    
    def _check_weights(
        self, msname, row_to_rows, data_column, chan_flags, combine_corr,
        chanbins, wtrange
    ):
        if data_column.startswith('c'):
            colname = 'CORRECTED_DATA'
        elif data_column.startswith('d'):
            colname = 'DATA'
        else:
            raise Exception("Unhandled column spec " + data_column)
        for ant1 in range(10):
            for ant2 in range((ant1 + 1), 10):
                query_str = 'ANTENNA1=' + str(ant1) + ' AND ANTENNA2=' \
                     + str(ant2)
                tb.open(msname)
                subt = tb.query(query_str)
                data = subt.getcol(colname)
                flags = subt.getcol('FLAG')
                exposures = subt.getcol('EXPOSURE')
                wt = subt.getcol('WEIGHT')
                wtsp = subt.getcol('WEIGHT_SPECTRUM')
                flag_row = subt.getcol('FLAG_ROW')
                subt.done()
                tb.done()
                nrows = data.shape[2]
                for row in range(nrows):
                    start = row_to_rows[row][0]
                    end = row_to_rows[row][1]
                    (weights, ewt, mod_flags) = get_weights(
                        data[:,:,start:end], flags[:, :, start:end], chan_flags,
                        exposures[start: end], combine_corr, exposures[row],
                        chanbins, flags[:, :, row:row+1], wtrange
                    )
                    self.assertTrue(
                        np.allclose(weights, wtsp[:, :, row]),
                        'Failed wtsp, got ' + str(wtsp[:, :, row])
                        + '\nexpected ' + str(weights) + '\nbaseline '
                        + str([ant1, ant2]) + '\nrow ' + str(row)
                    )
                    self.assertTrue(
                        np.allclose(ewt, wt[:, row]),
                        'Failed weight, got ' + str(wt[:, row])
                        + '\nexpected ' + str(np.median(weights, 1))
                        + '\nbaseline ' + str([ant1, ant2]) + '\nrow '
                        + str(row)
                    )
                    self.assertTrue(
                        (mod_flags == np.expand_dims(flags[:, :, row], 2)).all(),
                        'Failed flag, got ' + str(flags[:, :, row])
                        + '\nexpected ' + str(mod_flags) + '\nbaseline '
                        + str([ant1, ant2]) + '\nrow ' + str(row)
                    )
                    eflag_row = mod_flags.all()
                    self.assertTrue(
                        (eflag_row == flag_row[row]).all(),
                        'Failed flag_row, got ' + str(flag_row[row])
                        + '\nexpected ' + str(eflag_row) + '\nbaseline '
                        + str([ant1, ant2]) + '\nrow ' + str(row)
                    )
                    # all flags must be True where wtsp = 0
                    self.assertTrue(np.extract(weights == 0, mod_flags).all())

    def compare(self, dst, ref):
        mytb = table()
        self.assertTrue(mytb.open(dst), "Table open failed for " + dst)
        [gtimes, gwt, gwtsp, gflag, gfrow, gdata] = _get_table_cols(mytb)
        mytb.done()
        ref = os.path.join(datadir, ref)
        print("ref in compare", ref)
        self.assertTrue(mytb.open(ref), "Table open failed for " + ref)
        [etimes, ewt, ewtsp, eflag, efrow, edata] = _get_table_cols(mytb)
        mytb.done()
        self.assertTrue(np.allclose(gwt, ewt), 'WEIGHT comparison failed')
        self.assertTrue(
            np.allclose(gwtsp, ewtsp), 'WEIGHT_SPECTRUM comparison failed'
        )
        self.assertTrue((gflag == eflag).all(), 'FLAG comparison failed')
        self.assertTrue((gfrow == efrow).all(), 'FLAG_ROW comparison failed')
        # all flags must be True where wtsp = 0
        self.assertTrue(np.extract(gwtsp == 0, gflag).all())


    def test_algorithm(self):
        """ Test the algorithm, includes excludechans tests"""
        mytb = table()
        dst = "ngc5921.split.ms"
        cflags = np.array(63 * [False])
        cflags[10:21] = True
        myms = ms()
        row_to_rows = []
        for row in range(60):
            row_to_rows.append((row, row+1))
        for combine in ["", "corr"]:
            c = 0
            for fitspw in ["0:0~9;21~62", "", "0:10~20"]:
                shutil.copytree(ctsys.resolve(src), dst)
                excludechans = c == 2
                myms.open(dst, nomodify=False)
                myms.statwt(
                    combine=combine, fitspw=fitspw,
                    excludechans=excludechans
                )
                myms.done()
                chan_flags = cflags if fitspw else None
                if combine == '':
                    if fitspw == '':
                        ref = 'ref_test_algorithm_sep_corr_no_fitspw.ms'
                    else: 
                        ref = 'ref_test_algorithm_sep_corr_fitspw.ms'
                else:
                    if fitspw == '':
                        ref = 'ref_test_algorithm_combine_corr_no_fitspw.ms'
                    else:
                        ref = 'ref_test_algorithm_combine_corr_has_fitspw.ms'
                self.compare(dst, ref)
                shutil.rmtree(dst)
                c += 1               

    def test_timebin(self):
        """ Test time binning"""
        dst = "ngc5921.split.timebin.ms"
        combine = "corr"
        for timebin in ["300s", 10]:
            shutil.copytree(src, dst) 
            myms = ms()
            myms.open(dst, nomodify=False)
            myms.statwt(timebin=timebin, combine=combine)
            myms.done()
            ref = 'ref_test_timebin_' + str(timebin) + '.ms'
            self.compare(dst, ref)
            shutil.rmtree(dst)
            
    def test_chanbin(self):
        """Test channel binning"""
        dst = "ngc5921.split.chanbin_0.ms"
        row_to_rows = []
        for i in range(60):
            row_to_rows.append([i, i+1])
        bins = [
            [0, 8], [8, 16], [16, 24], [24, 32], [32, 40], [40, 48],
            [48, 56], [56,63]
        ]
        myms = ms()
        for combine in ["", "corr"]:
            for i in [0, 2]:
                for chanbin in ["195.312kHz", 8]:
                    if i == 2 and combine != '' and chanbin != 8:
                        # only run the check for i == 2 once
                        continue
                    shutil.copytree(src, dst)
                    if i == 0:
                        myms.open(dst, nomodify=False)
                        myms.statwt(chanbin=chanbin, combine=combine)
                        myms.done()
                    elif i == 2:
                        # check WEIGHT_SPECTRUM is created, only check once,
                        # this test is long as it is
                        mytb = table()
                        mytb.open(dst, nomodify=False)
                        x = mytb.ncols()
                        self.assertTrue(
                            mytb.removecols("WEIGHT_SPECTRUM"),
                            "column not removed"
                        )
                        y = mytb.ncols()
                        self.assertTrue(y == x-1, "wrong number of columns")
                        mytb.done()
                        myms.open(dst, nomodify=False)
                        myms.statwt(chanbin=chanbin, combine=combine)
                        myms.done()
                    if combine == '':
                        ref = datadir + 'ref_test_chanbin_sep_corr.ms'
                    else:
                        ref = datadir + 'ref_test_chanbin_combine_corr.ms'
                    shutil.rmtree(dst)

    def test_minsamp(self):
        """Test minimum number of points"""
        dst = "ngc5921.split.minsamp.ms"
        combine = "corr"
        trow = 12
        for minsamp in [60, 80]:
            shutil.copytree(ctsys.resolve(src), dst)
            myms = ms()
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
        ref = 'ref_test_default_boundaries.ms'
        timebin = "6000s"
        # there are three field_ids, and there is a change in field_id when
        # there is a change in scan number, so specifying combine="field" in the
        # absence of "scan" will give the same result as combine=""
        row_to_rows = []
        for i in range(12):
            row_to_rows.append([0, 12])
        for i in range(12, 17):
            row_to_rows.append([12, 17])
        for i in range(17, 33):
            row_to_rows.append([17, 33])
        for i in range(33, 35):
            row_to_rows.append([33, 35])
        for i in range(35, 38):
            row_to_rows.append([35, 38])
        for i in range(38, 56):
            row_to_rows.append([38, 56])
        for i in range(56, 60):
            row_to_rows.append([56, 60])
        myms = ms()
        for combine in ["corr", "corr,field"]:
            shutil.copytree(src, dst)
            myms.open(dst, nomodify=False)
            myms.statwt(timebin=timebin, combine=combine)
            myms.done()
            self.compare(dst, ref)
            shutil.rmtree(dst)
        
    def test_no_scan_boundaries(self):
        """Test no scan boundaries"""
        dst = "ngc5921.no_scan_bounds.ms"
        timebin = "6000s"
        # ref = os.path.join(datadir, 'ref_test_no_scan_bounds.ms')
        ref = 'ref_test_no_scan_bounds.ms'
        combine = "corr, scan"
        shutil.copytree(ctsys.resolve(src), dst)
        myms = ms()
        myms.open(dst, nomodify=False)
        myms.statwt(timebin=timebin, combine=combine)
        myms.done()
        self.compare(dst, ref)
        shutil.rmtree(dst)
        
    def test_no_scan_nor_field_boundaries(self):
        """Test no scan nor field boundaries"""
        dst = "ngc5921.no_scan_nor_field_bounds.ms"
        timebin = "6000s"
        # ref = os.path.join(datadir, 'ref_test_no_scan_nor_field_bounds.ms')
        ref = 'ref_test_no_scan_nor_field_bounds.ms'
        for combine in ["corr,scan,field", "corr,field,scan"]:
            shutil.copytree(src, dst)
            myms = ms()
            myms.open(dst, nomodify=False)
            myms.statwt(timebin=timebin, combine=combine)
            myms.done()
            self.compare(dst, ref)
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
        """Test weight range"""
        dst = "ngc5921.split.timebin.ms"
        ref = "ref_test_wtrange_300s.ms"
        combine = "corr"
        timebin = "300s"
        wtrange = [1, 2]
        """
        row_to_rows = []
        for i in range(10):
            row_to_rows.append([0, 10])
        for i in range(2):
            row_to_rows.append([10, 12])
        for i in range(5):
            row_to_rows.append([12, 17])
        for i in range(5):
            row_to_rows.append([17, 22])
        for i in range(5):
            row_to_rows.append([22, 27])
        for i in range(5):
            row_to_rows.append([27, 32])
        for i in range(1):
            row_to_rows.append([32, 33])
        for i in range(2):
            row_to_rows.append([33, 35])
        for i in range(3):
            row_to_rows.append([35, 38])
        for i in range(5):
            row_to_rows.append([38, 43])
        for i in range(5):
            row_to_rows.append([43, 48])
        for i in range(5):
            row_to_rows.append([48, 53])
        for i in range(3):
            row_to_rows.append([53, 56])
        for i in range(4):
            row_to_rows.append([56, 60])
        """
        myms = ms()
        for i in [0, 1]:
            shutil.copytree(src, dst) 
            myms.open(dst, nomodify=False)
            myms.statwt(timebin=timebin, combine=combine, wtrange=wtrange)
            myms.done()
            self.compare(dst, ref)
            # self._check_weights(
            #    dst, row_to_rows, 'c', None, True, None, wtrange
            # )
            shutil.rmtree(dst)

    def test_preview(self):
        """ Test preview mode"""
        dst = "ngc5921.split.preview.ms"
        [refwt, refwtsp, refflag, reffrow, refdata] = _get_dst_cols(src)
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
        self.assertTrue(np.all(tstflag == refflag), "FLAGs don't match")
        self.assertTrue(np.all(tstfrow == reffrow), "FLAG_ROWs don't match")
        self.assertTrue(
            np.all(np.isclose(tstwt, refwt)), "WEIGHTs don't match"
        )
        self.assertTrue(
            np.all(np.isclose(tstwtsp, refwtsp)), "WEIGHT_SPECTRUMs don't match"
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
