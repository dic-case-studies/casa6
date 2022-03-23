#########################################################################
# test_regression_timeaverage_and_rflag.py
# Copyright (C) 2018
# Associated Universities, Inc. Washington DC, USA.
#
# This script is free software; you can redistribute it and/or modify it
# under the terms of the GNU Library General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
# License for more details.
#
# Based on examples from CAS-11910
#
# Rationale for Inclusion:
#    test flagdata autoflagging (rflag mode) with on-the-fly time average
#    and backpropagation of flags, when the MS has pre-existing flags
#    and other tricky traits.
#
# MS datasets used in the tests
#   Four_ants_3C286.ms:
#           simulated EVLA MS with 4 spws
#   3ctst.ms:
#           from the CAS-11910 original description. GMRT MS that has un-ordered (in time) rows.
#  That makes it particularly handy to test and expose issues in time average / AveragingTVI and
#  related functionality.
#
##########################################################################

import time
import shutil
import unittest

from casatools import ctsys
from casatasks import flagdata, mstransform

datapath = ctsys.resolve("regression/time_average_and_rflag/")

# Input MS
ms_four_ants = 'Four_ants_3C286.ms'
ms1 = "Four_ants_copy1.ms"
ms2 = "Four_ants_copy2.ms"

ms3ctst = "3ctst_gmrt_unordered_time_rows.ms"
ms3ctst_copy = "3ctst_copy.ms"

totaltime = 0
inittime = time.time()
ttime = inittime
steptime = []


def timing(stepname=''):
    global totaltime
    global inittime
    global ttime
    global steptime
    global step_title
    global mystep
    global thesteps
    thetime = time.time()
    dtime = thetime - ttime
    steptime.append(dtime)
    totaltime += dtime
    ttime = thetime
    casalog.origin('TIMING')
    casalog.post( 'Step '+stepname, 'WARN')
    casalog.post( 'Time now: '+str(ttime), 'WARN')
    casalog.post( 'Time used this step: '+str(dtime), 'WARN')
    casalog.post( 'Total time used so far: ' + str(totaltime), 'WARN')
    casalog.post( 'Step  Time used (s)     Fraction of total time (percent) [description]', 'WARN')
    for i in range(0, len(steptime)):
        casalog.post('  ' + stepname + '   ' + str(steptime[i]) + '  ' + str(steptime[i] / totaltime * 100.)
                     + ' [' + stepname + ']', 'WARN')


def flag_step(scale=7.0, timebin='1min', step_name='', field='', spw='' , ms=ms1):
    """Run flagdata with mode='rflag' and timeavg=True on selected data.
       Run and save the summary of flags before and after flagging."""

    # Using spw selection also in summary mode
    summary_orig = flagdata(vis=ms, spw=spw, mode='summary')
    casalog.post(' * Flagging {0}'.format(ms))
    flagdata(vis=ms, mode='rflag', field=field, spw=spw, timeavg=True,
             timebin=timebin, freqdevscale=scale, timedevscale=scale, datacolumn='data',
             extendflags=False, correlation='',combinescans=True, display='none',
             action='apply', flagbackup=False)

    # Using spw selection also in summary mode
    summary_after = flagdata(vis=ms, spw=spw, mode='summary')
    casalog.post(' * Step {0}. Total flagged before: {1} ({2:.3f}%), after: {3} ({4:.3f}%). '
                 'Used timebin: {5}, scale: {6}'.format(step_name,
                        summary_orig['flagged'], 100.0 * summary_orig['flagged'] / summary_orig['total'],
                        summary_after['flagged'], 100.0 * summary_after['flagged'] / summary_after['total'],
                        timebin, scale))

    return summary_orig['flagged'], summary_after['flagged']


def mstransform_flag_step(scale=7.0, timebin='1min', step_name=''):
    """
    Flags with time average, following two different approaches
    and return the flag summary of both approaches

    Approach 1:
        step 1: run mstransform on original MS with timeaverage=True and save ms_step1
        step 2: run flagdagta on ms_step1 output with mode=rflag (no timeavg)
                save summary of flags of ms_step1.

    Approach 2:
        step 1: run flagdata on original MS with mode=rflag and on-the-fly timeavg=True.
                This step will back-propagate and save the flags (it will not save any average data)
        step 2: run mstransform on flagged MS from previous step with timeaverage=True. This step
                will match the average dimension of approach 1 step 1 and save in a new ms_step2.
                save the summary of flags of ms_step2.

    We expect summary of output1 to be the same of summary of output2.
    """

    casalog.post(' * Flagging {0}'.format(ms3ctst_copy))

    # Approach 1
    # step 1: Run mstransform with timeaverage=True and save output ms_step1
    ms_step1 = 'ms_approach1.ms'
    mstransform(vis=ms3ctst_copy, outputvis=ms_step1, datacolumn='data', timeaverage=True,
                timebin=timebin)

    # step 2: Run flagdata mode=rflag on averaged data of ms_step1
    flagdata(vis=ms_step1, mode='rflag', field='', spw='', freqdevscale=scale,
             timedevscale=scale, datacolumn='data', extendflags=False, correlation='',
             combinescans=True, display='none', action='apply', flagbackup=False)
    # Save summary of average flags of ms_step1
    res1 = flagdata(vis=ms_step1, mode='summary')

    # Approach 2
    # step 1: Run flagdata mode=rflag with on-the-fly time average. Average flags are not saved, only
    #            back-ported flags are saved.
    flagdata(vis=ms3ctst_copy, mode='rflag', field='', spw='', timeavg=True,
             timebin=timebin, freqdevscale=scale, timedevscale=scale, datacolumn='data',
             extendflags=False, correlation='', combinescans=True, display='none',
             action='apply', flagbackup=False)

    # step 2: Run mstransform with timeaverage=True on flagged MS and save output ms_step2
    ms_step2 = 'ms_approach2.ms'
    mstransform(vis=ms3ctst_copy, outputvis=ms_step2, datacolumn='data', timeaverage=True, timebin=timebin)
    # Save summary of average flags of ms_step2
    res2 = flagdata(vis=ms_step2, mode='summary')

    casalog.post('**** Flags of Approach 1: {}'.format(res1))
    casalog.post('**** Flags of Approach 2: {}'.format(res2))

    return res1['flagged'], res2['flagged']

class Timeaverage_and_Rflag(unittest.TestCase):

    def setUp(self):
        shutil.copytree(datapath+ms_four_ants, ms1)
        shutil.copytree(datapath+ms_four_ants, ms2)
        shutil.copytree(datapath+ms3ctst, ms3ctst_copy)

    def tearDown(self):
        shutil.rmtree(ms1)
        shutil.rmtree(ms2)
        shutil.rmtree(ms3ctst_copy)

    def test_rflag_and_timeavg_OTF(self):
        """Test that flagdata with on-the-fly time averaging accumulates flags
            when using mode=rflag for different timebins and freqdevscale/timedevscale values"""

        # Using Four_ants_3C286 MS copy1 with
        # timebin=30s, timedevscale=freqdevscale=7, accumulate flags 3 times
        flag_30s_7_1_before, flag_30s_7_1_after = flag_step(7, '30s', step_name='1', ms=ms1)
        flag_30s_7_2_before, flag_30s_7_2_after = flag_step(7, '30s', step_name='2', ms=ms1)
        flag_30s_7_3_before, flag_30s_7_3_after = flag_step(7, '30s', step_name='3', ms=ms1)

        timing('0')

        # Using Four_ants_3C286 MS copy2 with
        # timebin=30s, timedevscale=freqdevscale=3, accumulate flags 7 times
        flag_30s_3_1_before, flag_30s_3_1_after = flag_step(3, '30s', step_name='1', ms=ms2)
        flag_30s_3_2_before, flag_30s_3_2_after = flag_step(3, '30s', step_name='2', ms=ms2)
        flag_30s_3_3_before, flag_30s_3_3_after = flag_step(3, '30s', step_name='3', ms=ms2)
        flag_30s_3_4_before, flag_30s_3_4_after = flag_step(3, '30s', step_name='4', ms=ms2)
        flag_30s_3_5_before, flag_30s_3_5_after = flag_step(3, '30s', step_name='5', ms=ms2)
        flag_30s_3_6_before, flag_30s_3_6_after = flag_step(3, '30s', step_name='6', ms=ms2)
        flag_30s_3_7_before, flag_30s_3_7_after = flag_step(3, '30s', step_name='7', ms=ms2)

        timing('1')

        # Using 3ctst MS with
        # timebin=3min, timedevscale=freqdevscale=7, accumulate 7 times
        flag_30s_3_1_before_3ct, flag_30s_3_1_after_3ct = flag_step(7, '3min', step_name='1', field='0', spw='0:20~350', ms=ms3ctst_copy)
        flag_30s_3_2_before_3ct, flag_30s_3_2_after_3ct = flag_step(7, '3min', step_name='2', field='0', spw='0:20~350', ms=ms3ctst_copy)
        flag_30s_3_3_before_3ct, flag_30s_3_3_after_3ct = flag_step(7, '3min', step_name='3', field='0', spw='0:20~350', ms=ms3ctst_copy)
        flag_30s_3_4_before_3ct, flag_30s_3_4_after_3ct = flag_step(7, '3min', step_name='4', field='0', spw='0:20~350', ms=ms3ctst_copy)
        flag_30s_3_5_before_3ct, flag_30s_3_5_after_3ct = flag_step(7, '3min', step_name='5', field='0', spw='0:20~350', ms=ms3ctst_copy)
        flag_30s_3_6_before_3ct, flag_30s_3_6_after_3ct = flag_step(7, '3min', step_name='6', field='0', spw='0:20~350', ms=ms3ctst_copy)
        flag_30s_3_7_before_3ct, flag_30s_3_7_after_3ct = flag_step(7, '3min', step_name='7', field='0', spw='0:20~350', ms=ms3ctst_copy)

        timing('2')

        # Expected summary flags for Four_ants_3C286 copy1 runs
        expected_30s_7_1_before = 211894.0
        expected_30s_7_2_before = 226824.0
        expected_30s_7_3_before = 227004.0

        expected_30s_7_1_after = expected_30s_7_2_before
        expected_30s_7_2_after = expected_30s_7_3_before
        expected_30s_7_3_after = 227004.0

        # Expected summary flags for Four_ants_3C286 copy2 runs
        expected_30s_3_1_before = 211894.0
        expected_30s_3_2_before = 471802.0
        expected_30s_3_3_before = 602582.0
        expected_30s_3_4_before = 687092.0
        expected_30s_3_5_before = 736226.0
        expected_30s_3_6_before = 776867.0
        expected_30s_3_7_before = 817331.0

        expected_30s_3_1_after = expected_30s_3_2_before
        expected_30s_3_2_after = expected_30s_3_3_before
        expected_30s_3_3_after = expected_30s_3_4_before
        expected_30s_3_4_after = expected_30s_3_5_before
        expected_30s_3_5_after = expected_30s_3_6_before
        expected_30s_3_6_after = expected_30s_3_7_before
        expected_30s_3_7_after = 843051.0

        # Expected summary flags of 3ctst copy runs
        expected_30s_3_1_before_3ct = 9593008.0
        expected_30s_3_2_before_3ct = 9611700.0
        expected_30s_3_3_before_3ct = 9614850.0
        expected_30s_3_4_before_3ct = 9615270.0
        expected_30s_3_5_before_3ct = 9615343.0
        expected_30s_3_6_before_3ct = 9615350.0
        expected_30s_3_7_before_3ct = 9615354.0

        expected_30s_3_1_after_3ct = expected_30s_3_2_before_3ct
        expected_30s_3_2_after_3ct = expected_30s_3_3_before_3ct
        expected_30s_3_3_after_3ct = expected_30s_3_4_before_3ct
        expected_30s_3_4_after_3ct = expected_30s_3_5_before_3ct
        expected_30s_3_5_after_3ct = expected_30s_3_6_before_3ct
        expected_30s_3_6_after_3ct = expected_30s_3_7_before_3ct
        expected_30s_3_7_after_3ct = 9615354.0

        observed_list = [flag_30s_7_1_before, flag_30s_7_2_before, flag_30s_7_3_before,
                         flag_30s_7_1_after, flag_30s_7_2_after, flag_30s_7_3_after,
                         flag_30s_3_1_before, flag_30s_3_2_before, flag_30s_3_3_before,
                         flag_30s_3_4_before, flag_30s_3_5_before, flag_30s_3_6_before, flag_30s_3_7_before,
                         flag_30s_3_1_after, flag_30s_3_2_after, flag_30s_3_3_after,
                         flag_30s_3_4_after, flag_30s_3_5_after, flag_30s_3_6_after, flag_30s_3_7_after,
                         flag_30s_3_1_before_3ct, flag_30s_3_2_before_3ct, flag_30s_3_3_before_3ct,
                         flag_30s_3_4_before_3ct, flag_30s_3_5_before_3ct, flag_30s_3_6_before_3ct, flag_30s_3_7_before_3ct,
                         flag_30s_3_1_after_3ct, flag_30s_3_2_after_3ct, flag_30s_3_3_after_3ct,
                         flag_30s_3_4_after_3ct, flag_30s_3_5_after_3ct, flag_30s_3_6_after_3ct, flag_30s_3_7_after_3ct]

        expected_list = [expected_30s_7_1_before, expected_30s_7_2_before, expected_30s_7_3_before,
                         expected_30s_7_1_after, expected_30s_7_2_after, expected_30s_7_3_after,
                         expected_30s_3_1_before, expected_30s_3_2_before, expected_30s_3_3_before,
                         expected_30s_3_4_before, expected_30s_3_5_before, expected_30s_3_6_before, expected_30s_3_7_before,
                         expected_30s_3_1_after, expected_30s_3_2_after, expected_30s_3_3_after,
                         expected_30s_3_4_after, expected_30s_3_5_after, expected_30s_3_6_after, expected_30s_3_7_after,
                         expected_30s_3_1_before_3ct, expected_30s_3_2_before_3ct, expected_30s_3_3_before_3ct,
                         expected_30s_3_4_before_3ct, expected_30s_3_5_before_3ct, expected_30s_3_6_before_3ct, expected_30s_3_7_before_3ct,
                         expected_30s_3_1_after_3ct, expected_30s_3_2_after_3ct, expected_30s_3_3_after_3ct,
                         expected_30s_3_4_after_3ct, expected_30s_3_5_after_3ct, expected_30s_3_6_after_3ct, expected_30s_3_7_after_3ct]

        # Check that expected and observed lists are the same or fail if they are not
        for l_obs, l_exp in zip(observed_list,expected_list):
            self.assertEqual(l_obs, l_exp, "Observed value {} does not equal expected value {}".format(l_obs, l_exp))

        timing('3')

        shutil.rmtree(ms3ctst_copy)
        shutil.copytree(datapath+ms3ctst, ms3ctst_copy)


class Mstransform_and_Flagdata_TimeAverage(unittest.TestCase):

    def setUp(self):
        shutil.copytree(datapath+ms3ctst, ms3ctst_copy)

    def tearDown(self):
        shutil.rmtree(ms3ctst_copy)
        shutil.rmtree("ms_approach1.ms")
        shutil.rmtree("ms_approach2.ms")

    def test_timeaverage_and_flags(self):
        """Run mstransform with timeaverage and compare with flagdata OTF time average"""

        # Compare flags of approach 1 and approach 2
        flags1, flags2 = mstransform_flag_step(step_name='1', timebin='2min', scale=7)
        self.assertEqual(flags1, flags2)

        shutil.rmtree("ms_approach1.ms")
        shutil.rmtree("ms_approach2.ms")

        flags1, flags2 = mstransform_flag_step(step_name='2', timebin='6min', scale=7)
        self.assertEqual(flags1, flags2)


if __name__ == "__main__":
    unittest.main()
