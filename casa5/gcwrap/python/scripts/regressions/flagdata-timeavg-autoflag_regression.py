#############################################################################
# $Id:$
# Test Name:                                                                #
# flagdata-timeavg-autoplag_regression.py                                   #
#                                                                           #
# Test examples from CAS-11910                                              #
#                                                                           #
# Rationale for Inclusion:                                                  #
#    Need test of complete ALMA analysis chain                              #
#                                                                           #
# Input data:                                                               #
#     two ASDMs                                                             #
#     the clean masks                                                       #
#                                                                           #
#############################################################################

import shutil
import time

thesteps = []

step_title = {0: 'flag with scale 7 for 4 flags',
            1: 'flag with scale 3 for 4 flags',
            2: 'flag with scale 7 for ct',
            3: 'Checks'
            }

try:
    print("List of steps to be executed ... ", mysteps)
    thesteps = mysteps
except:
    print("global variable mysteps not set.")

if (thesteps==[]):
    thesteps = range(0, len(step_title))
    print("mysteps empty. Executing all steps: ", thesteps)

totaltime = 0
inittime = time.time()
ttime = inittime
steptime = []

def timing():
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
    casalog.post( 'Step '+str(mystep)+': '+step_title[mystep], 'WARN')
    casalog.post( 'Time now: '+str(ttime), 'WARN')
    casalog.post( 'Time used this step: '+str(dtime), 'WARN')
    casalog.post( 'Total time used so far: ' + str(totaltime), 'WARN')
    casalog.post( 'Step  Time used (s)     Fraction of total time (percent) [description]', 'WARN')
    for i in range(0, len(steptime)):
        casalog.post('  ' + str(thesteps[i]) + '   ' + str(steptime[i]) + '  ' + str(steptime[i] / totaltime * 100.)
                     + ' [' + step_title[thesteps[i]] + ']', 'WARN')


my_dataset_name = "Four_ants_3C286.ms"
ms3ctst = "3ctst.ms"

ms1 = "Four_ants_3C286.ms"
ms2 = "Four_ants_copy.ms"
shutil.copytree(ms1, ms2)


def flag_step(scale=7.0, timebin='1min', step_name='', field='', spw='' , ms=ms):
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
    print(' * Step {0}. Total flagged before: {1} ({2:.3f}%), after: {3} ({4:.3f}%). '
                 'Used timebin: {5}, scale: {6}'.format(step_name,
                                                        summary_orig['flagged'],
                                                        100.0 * summary_orig['flagged'] / summary_orig['total'],
                                                        summary_after['flagged'],
                                                        100.0 * summary_after['flagged'] / summary_after['total'],
                                                        timebin, scale))

    return summary_orig['flagged'], summary_after['flagged']

def unflag_then_flag(scale=7.0, timebin='1min', field='', spw='', ms=ms):
    flagdata(vis=ms, mode='unflag')
    summary_unflagged = flagdata(vis=ms, spw=spw, mode='summary')
    # Note display='none', would be 'both' but that requires user interaction
    flagdata(vis=ms, mode='rflag', field=field, spw=spw, timeavg=True,
             timebin=timebin, freqdevscale=scale, timedevscale=scale, datacolumn='data',
             extendflags=False, correlation='',combinescans=True, display='none',
             action='apply', flagbackup=False)
    summary_reflag = flagdata(vis=ms, spw=spw, mode='summary')
    casalog.post(' * Total flagged after unflagging: {0} ({1:.3f}%), and then after flagging from scratch: {2} ({3:.3f}%). '
                 'Used timebin: {4}, scale: {5}'.
                 format(summary_unflagged['flagged'], 100.0 * summary_unflagged['flagged']/summary_unflagged['total'],
                        summary_reflag['flagged'], 100.0 *summary_reflag['flagged'] / summary_reflag['total'],
                        timebin, scale))


# Step 0
mystep = 0
if (mystep in thesteps):
    print('Step ', mystep, step_title[mystep])

    flag_30s_7_1_before, flag_30s_7_1_after = flag_step(7, '30s', step_name='1', ms=ms1)
    flag_30s_7_2_before, flag_30s_7_2_after = flag_step(7, '30s', step_name='2', ms=ms1)
    flag_30s_7_3_before, flag_30s_7_3_after = flag_step(7, '30s', step_name='3', ms=ms1)

    timing()

# Step 1
mystep = 1
if (mystep in thesteps):
    print('Step ', mystep, step_title[mystep])

    flag_30s_3_1_before, flag_30s_3_1_after = flag_step(3, '30s', step_name='1', ms=ms2)
    flag_30s_3_2_before, flag_30s_3_2_after = flag_step(3, '30s', step_name='2', ms=ms2)
    flag_30s_3_3_before, flag_30s_3_3_after = flag_step(3, '30s', step_name='3', ms=ms2)
    flag_30s_3_4_before, flag_30s_3_4_after = flag_step(3, '30s', step_name='4', ms=ms2)
    flag_30s_3_5_before, flag_30s_3_5_after = flag_step(3, '30s', step_name='5', ms=ms2)
    flag_30s_3_6_before, flag_30s_3_6_after = flag_step(3, '30s', step_name='6', ms=ms2)
    flag_30s_3_7_before, flag_30s_3_7_after = flag_step(3, '30s', step_name='7', ms=ms2)


    timing()

# Step 2
mystep = 2
if (mystep in thesteps):
    print('Step ', mystep, step_title[mystep])

    flag_30s_3_1_before_3ct, flag_30s_3_1_after_3ct = flag_step(7, '3min', step_name='1', field='0', spw='0:20~350', ms=ms3ctst)
    flag_30s_3_2_before_3ct, flag_30s_3_2_after_3ct = flag_step(7, '3min', step_name='2', field='0', spw='0:20~350', ms=ms3ctst)
    flag_30s_3_3_before_3ct, flag_30s_3_3_after_3ct = flag_step(7, '3min', step_name='3', field='0', spw='0:20~350', ms=ms3ctst)
    flag_30s_3_4_before_3ct, flag_30s_3_4_after_3ct = flag_step(7, '3min', step_name='4', field='0', spw='0:20~350', ms=ms3ctst)
    flag_30s_3_5_before_3ct, flag_30s_3_5_after_3ct = flag_step(7, '3min', step_name='5', field='0', spw='0:20~350', ms=ms3ctst)
    flag_30s_3_6_before_3ct, flag_30s_3_6_after_3ct = flag_step(7, '3min', step_name='6', field='0', spw='0:20~350', ms=ms3ctst)
    flag_30s_3_7_before_3ct, flag_30s_3_7_after_3ct = flag_step(7, '3min', step_name='7', field='0', spw='0:20~350', ms=ms3ctst)


    timing()

regstate = False
# Step 3
mystep = 3
if (mystep in thesteps):
    print('Step ', mystep, step_title[mystep])

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

    for i in range(len(observed_list)):
        if observed_list[i] != expected_list[i]:
            lists_equal = False
            print("Observed value {} does not equal expected value {}".format(observed_list[i], expected_list[i]))


    shutil.rmtree(ms2)
    timing()

if lists_equal:
    print("Regression PASSED")
else:
    print("Regression FAILED")
print("Done")
