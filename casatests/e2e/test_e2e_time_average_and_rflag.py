#############################################################################
# $Id:$
# Test Name:                                                                #
# test_e2e_time_average_and_rflag.py                                        #
#                                                                           #
# Test examples from CAS-11910                                              #
#                                                                           #
# Rationale for Inclusion:                                                  #
#    test autoflagging (rflag) with on-the-fly time average                 #
#    and backpropagation of flags, when the MS has pre-existing flags       #
#    and other tricky traits.                                               #
#                                                                           #
# Input data:                                                               #
#     two MS                                                                #
#                                                                           #
#                                                                           #
#############################################################################

import time
import shutil
import unittest
CASA6 = False

try:
    import casatools
    from casatasks import flagdata, casalog, mstransform
    CASA6 = True
except:
    from tasks import *
    from taskinit import *


# TODO: Data, Need to ask where this data should go
# this is a placeholder location
if CASA6:
    casapath = ''
    datapath_four_ants = casatools.ctsys.resolve("regression/flagdata/Four_ants_3C286.ms/")
    datapath_3ctst = casatools.ctsys.resolve("regression/flagdata/3ctst.ms")
else:
    casapath = os.environ['CASAPATH'].split()[0]
    datapath = ''

if 'datasets' not in (locals()):
    myname = 'time_average_and_rflag :'
    mydict = { 1: 'Four_ants_3C286',
               2: '3ctst'}
else:
    myname = ' '
    mydict = (locals())['datasets']


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


my_dataset_name = "Four_ants_3C286.ms"

ms3ctst = "3ctst.ms"
ms3ctst_copy = "3ctst_copy.ms"

ms1 = "Four_ants_copy1.ms"
ms2 = "Four_ants_copy2.ms"


def flag_step(scale=7.0, timebin='1min', step_name='', field='', spw='' , ms=ms1):
    # Note using spw_sel also in the summary
    summary_orig = flagdata(vis=ms, spw=spw, mode='summary')
    casalog.post(' * Flagging {0}'.format(ms))
    # Note display='none', would be 'both' but that requires user interaction
    flagdata(vis=ms, mode='rflag', field=field, spw=spw, timeavg=True,
             timebin=timebin, freqdevscale=scale, timedevscale=scale, datacolumn='data',
             extendflags=False, correlation='',combinescans=True, display='none',
             action='apply', flagbackup=False)
    summary_after = flagdata(vis=ms, spw=spw, mode='summary')
    casalog.post(' * Step {0}. Total flagged before: {1} ({2:.3f}%), after: {3} ({4:.3f}%). '
                 'Used timebin: {5}, scale: {6}'.format(step_name,
                        summary_orig['flagged'], 100.0 * summary_orig['flagged'] / summary_orig['total'],
                        summary_after['flagged'], 100.0 * summary_after['flagged'] / summary_after['total'],
                        timebin, scale))

    return summary_orig['flagged'], summary_after['flagged']


def mstransform_flag_step(scale=7.0, timebin='1min', step_name=''):
    """
    Flags with time averge, following two alternative approaches. And checks (crude assert)
    that both produce identical results.

    Approach "step 1": command 1: mstransform with timeavg =>
       command 2: flagdata applied on the averaged MS (without timeavg).

    Approach "step 2": command 1: flagdata with timeavg on the fly =>
       command 2: mstransform with timeavg to match the dimensions of A).
    """
    summary_orig = flagdata(vis=ms3ctst_copy, spw='', mode='summary')
    casalog.post(' * Flagging {0}'.format(ms3ctst_copy))

    ms_step1 = 'ms_CAS-11910_mstransform_step_name{}_methodA.ms'.format(step_name)
    mstransform(vis=ms3ctst_copy, outputvis=ms_step1, datacolumn='data', timeaverage=True,
                timebin=timebin)
    flagdata(vis=ms_step1, mode='rflag', field='', spw='', freqdevscale=scale,
             timedevscale=scale, datacolumn='data', extendflags=False, correlation='',
             combinescans=True, display='none', action='apply', flagbackup=False)
    res1 = flagdata(vis=ms_step1, mode='summary')

    flagdata(vis=ms3ctst_copy, mode='rflag', field='', spw='', timeavg=True,
             timebin=timebin, freqdevscale=scale, timedevscale=scale, datacolumn='data',
             extendflags=False, correlation='', combinescans=True, display='none',
             action='apply', flagbackup=False)
    summary_after = flagdata(vis=ms3ctst_copy, spw='', mode='summary')
    ms_step2 = 'ms_CAS-11910_mstransform_step_name{}_methodB.ms'.format(step_name)
    mstransform(vis=ms3ctst_copy, outputvis=ms_step2, datacolumn='data', timeaverage=True, timebin=timebin)
    res2 = flagdata(vis=ms_step2, mode='summary')

    print('* Step {}. ResultA: {}'.format(step_name, res1))
    print('* Step {}. ResultB: {}'.format(step_name, res2))
    # Check reslts from flagdata-timeavg and mstransform-timeavg are identical
    assert res1 == res2
    summary_after = flagdata(vis=ms3ctst_copy, spw='', mode='summary')
    casalog.post(' * Step {0}. Total flagged before: {1} ({2:.3f}%), after: {3} ({4:.3f}%). '
                 'Used timebin: {5}, scale: {6}'.
                 format(step_name,
                        summary_orig['flagged'], 100.0 * summary_orig['flagged'] / summary_orig['total'],
                        summary_after['flagged'], 100.0 * summary_after['flagged'] / summary_after['total'],
                        timebin, scale))



class regression_time_average_rflag(unittest.TestCase):

    def setUp(self):

        if not CASA6:
            default(flagdata)

        shutil.copytree(my_dataset_name, ms1)
        shutil.copytree(my_dataset_name, ms2)
        shutil.copytree(ms3ctst, ms3ctst_copy)

    def tearDown(self):

        shutil.rmtree(ms1)
        shutil.rmtree(ms2)
        shutil.rmtree(ms3ctst_copy)
        shutil.rmtree("ms_CAS-11910_mstransform_step_name1_methodA.ms/")
        shutil.rmtree("ms_CAS-11910_mstransform_step_name1_methodB.ms/")

        shutil.rmtree("ms_CAS-11910_mstransform_step_name2_methodA.ms/")
        shutil.rmtree("ms_CAS-11910_mstransform_step_name2_methodB.ms/")

    def test_regression(self):

        flag_30s_7_1_before, flag_30s_7_1_after = flag_step(7, '30s', step_name='1', ms=ms1)
        flag_30s_7_2_before, flag_30s_7_2_after = flag_step(7, '30s', step_name='2', ms=ms1)
        flag_30s_7_3_before, flag_30s_7_3_after = flag_step(7, '30s', step_name='3', ms=ms1)

        timing('0')

        flag_30s_3_1_before, flag_30s_3_1_after = flag_step(3, '30s', step_name='1', ms=ms2)
        flag_30s_3_2_before, flag_30s_3_2_after = flag_step(3, '30s', step_name='2', ms=ms2)
        flag_30s_3_3_before, flag_30s_3_3_after = flag_step(3, '30s', step_name='3', ms=ms2)
        flag_30s_3_4_before, flag_30s_3_4_after = flag_step(3, '30s', step_name='4', ms=ms2)
        flag_30s_3_5_before, flag_30s_3_5_after = flag_step(3, '30s', step_name='5', ms=ms2)
        flag_30s_3_6_before, flag_30s_3_6_after = flag_step(3, '30s', step_name='6', ms=ms2)
        flag_30s_3_7_before, flag_30s_3_7_after = flag_step(3, '30s', step_name='7', ms=ms2)

        timing('1')

        flag_30s_3_1_before_3ct, flag_30s_3_1_after_3ct = flag_step(7, '3min', step_name='1', field='0', spw='0:20~350', ms=ms3ctst_copy)
        flag_30s_3_2_before_3ct, flag_30s_3_2_after_3ct = flag_step(7, '3min', step_name='2', field='0', spw='0:20~350', ms=ms3ctst_copy)
        flag_30s_3_3_before_3ct, flag_30s_3_3_after_3ct = flag_step(7, '3min', step_name='3', field='0', spw='0:20~350', ms=ms3ctst_copy)
        flag_30s_3_4_before_3ct, flag_30s_3_4_after_3ct = flag_step(7, '3min', step_name='4', field='0', spw='0:20~350', ms=ms3ctst_copy)
        flag_30s_3_5_before_3ct, flag_30s_3_5_after_3ct = flag_step(7, '3min', step_name='5', field='0', spw='0:20~350', ms=ms3ctst_copy)
        flag_30s_3_6_before_3ct, flag_30s_3_6_after_3ct = flag_step(7, '3min', step_name='6', field='0', spw='0:20~350', ms=ms3ctst_copy)
        flag_30s_3_7_before_3ct, flag_30s_3_7_after_3ct = flag_step(7, '3min', step_name='7', field='0', spw='0:20~350', ms=ms3ctst_copy)

        timing('2')

        lists_equal = True

        # Four ants outputs
        expected_30s_7_1_before = 211894.0
        expected_30s_7_2_before = 226824.0
        expected_30s_7_3_before = 227004.0

        expected_30s_7_1_after = 226824.0
        expected_30s_7_2_after = 227004.0
        expected_30s_7_3_after = 227004.0

        expected_30s_3_1_before = 211894.0
        expected_30s_3_2_before = 471802.0
        expected_30s_3_3_before = 602582.0
        expected_30s_3_4_before = 687092.0
        expected_30s_3_5_before = 736226.0
        expected_30s_3_6_before = 776867.0
        expected_30s_3_7_before = 817331.0

        expected_30s_3_1_after = 471802.0
        expected_30s_3_2_after = 602582.0
        expected_30s_3_3_after = 687092.0
        expected_30s_3_4_after = 736226.0
        expected_30s_3_5_after = 776867.0
        expected_30s_3_6_after = 817331.0
        expected_30s_3_7_after = 843051.0

        # 3ctst outputs
        expected_30s_3_1_before_3ct = 9593008.0
        expected_30s_3_2_before_3ct = 9611700.0
        expected_30s_3_3_before_3ct = 9614850.0
        expected_30s_3_4_before_3ct = 9615270.0
        expected_30s_3_5_before_3ct = 9615343.0
        expected_30s_3_6_before_3ct = 9615350.0
        expected_30s_3_7_before_3ct = 9615354.0

        expected_30s_3_1_after_3ct = 9611700.0
        expected_30s_3_2_after_3ct = 9614850.0
        expected_30s_3_3_after_3ct = 9615270.0
        expected_30s_3_4_after_3ct = 9615343.0
        expected_30s_3_5_after_3ct = 9615350.0
        expected_30s_3_6_after_3ct = 9615354.0
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

        ### mstransform test case

        #hutil.rmtree(ms3ctst_copy)
        #shutil.copytree(ms3ctst, ms3ctst_copy)

        #mstransform_flag_step(step_name='1', timebin='2min', scale=7)
        #mstransform_flag_step(step_name='2', timebin='6min', scale=7)


        for i in range(len(observed_list)):
            if observed_list[i] != expected_list[i]:
                lists_equal = False
                casalog.post("Observed value {} does not equal expected value {}".format(observed_list[i], expected_list[i]))
                print("Observed value {} does not equal expected value {}".format(observed_list[i], expected_list[i]))

        timing('3')

        if lists_equal:
            casalog.post("Regression PASSED")
            print("Regression PASSED")
            self.assertTrue(True)
        else:
            casalog.post("Regression FAILED")
            print("Regression FAILED")
            self.assertTrue(False)

        shutil.rmtree(ms3ctst_copy)
        shutil.copytree(ms3ctst, ms3ctst_copy)

        mstransform_flag_step(step_name='1', timebin='2min', scale=7)

        shutil.rmtree(ms3ctst_copy)
        shutil.copytree(ms3ctst, ms3ctst_copy)

        mstransform_flag_step(step_name='2', timebin='6min', scale=7)

        casalog.post("Done")

def suite():
    return[regression_time_average_rflag]

if __name__ == "__main__":
    unittest.main()
