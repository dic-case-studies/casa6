"""almastktestutils - utility module to be used for the maintenance of 
   test_stk_alma_pipeline_imaging.py
   
   extract_expdict: extract the fiducial metrics stored with tests in the old format
                    of test_stk_alma_pipeline_imaging.py

"""
import ast
import copy
import json
import os
import sys


def extract_expdict(testlist=None, testsrcpath=None):
    """ Read old test_alma_stk_pipeline_imaging.py and extract
    exp dictionaries of the specific test and save to a json file
    """
    if testlist is None:
        testlist = []
    if not os.path.exists(testsrcpath):
        raise IOError("{} not found".format(testsrcpath))

    if not testsrcpath.endswith('/'):
        testsrcpath += '/'
    oldstktestfile = testsrcpath + "test_stk_alma_pipeline_imaging.py"
    print("testcode =", oldstktestfile)
    if type(testlist) == list() and len(testlist) == 0:
        print("Scan all tests in {}".format(oldstktestfile))
        sys.path.append(testsrcpath)
        from src.test_stk_alma_pipeline_imaging import Test_standard
        testlist.extend([testname for testname in dir(Test_standard) if testname.startswith('test_') is True])
    for testname in testlist:
        print("Processing {}.... ".format(testname))

        outfile = testname + '_exp_dicts.json'
        # srcdir='/export/home/murasame/casa/casa6/'
        # testdir='casatests/stakeholder/'
        # abspath = srcdir+testdir

        with open(oldstktestfile) as f:
            readDict = False
            exp_im_stats_str = ""
            exp_im_stats_str_end = False
            exp_mask_stats_str = ""
            exp_mask_stats_str_end = False
            exp_pb_stats_str = ""
            exp_pb_stats_str_end = False
            exp_psf_stats_str = ""
            exp_psf_stats_str_end = False
            exp_model_stats_str = ""
            exp_model_stats_str_end = False
            exp_resid_stats_str = ""
            exp_resid_stats_str_end = False
            exp_sumwt_stats_str = ""
            exp_sumwt_stats_str_end = False
            exp_wt_stats_str = ""
            exp_wt_stats_str_end = False
            exp_bmin_dict_str = ""
            exp_bmaj_dict_str = ""
            exp_pa_dict_str = ""
            exp_im1_stats_str = ""
            exp_im1_stats_str_end = False
            exp_model1_stats_str = ""
            exp_model1_stats_str_end = False
            exp_resid1_stats_str = ""
            exp_resid1_stats_str_end = False
            exp_sumwt1_stats_str = ""
            exp_sumwt1_stats_str_end = False

            for ln in f:
                # skip the comment lines
                if ln.lstrip().startswith("#"):
                    pass
                else:
                    if "def " + testname + "(" in ln:
                        print("test {} found!".format(testname))
                        readDict = True
                    elif 'def test_' in ln and readDict:
                        print("Finish reading {}".format(testname))
                        break

                    # read each exp_ dictionary here (stop reading when detect "}"
                    if readDict:
                        if not len(exp_im_stats_str) and "exp_im_stats = {" in ln:
                            exp_im_stats_str = ln.lstrip("exp_im_stats = ").rstrip("\n")
                        elif len(exp_im_stats_str) and not exp_im_stats_str_end:
                            exp_im_stats_str += ln.rstrip("\n").rstrip(" ").replace("\\", '').lstrip()
                            if ln.rstrip("\n").rstrip().endswith("}"):
                                exp_im_stats_str_end = True
                                print("Stop reading exp_im_stats")
                        elif not len(exp_mask_stats_str) and "exp_mask_stats = {" in ln:
                            exp_mask_stats_str = ln.lstrip("exp_mask_stats = ").rstrip("\n")
                        elif len(exp_mask_stats_str) and not exp_mask_stats_str_end:
                            exp_mask_stats_str += ln.rstrip("\n").rstrip(" ").replace("\\", '')
                            if ln.rstrip("\n").rstrip().endswith("}"):
                                exp_mask_stats_str_end = True
                                print("Stop reading exp_mask_stats")
                        elif not len(exp_pb_stats_str) and "exp_pb_stats = {" in ln:
                            exp_pb_stats_str = ln.lstrip("exp_pb_stats = ").rstrip("\n")
                        elif len(exp_pb_stats_str) and not exp_pb_stats_str_end:
                            exp_pb_stats_str += ln.rstrip("\n").rstrip(" ").replace("\\", '')
                            if ln.rstrip("\n").rstrip().endswith("}"):
                                exp_pb_stats_str_end = True
                                print("Stop reading exp_pb_stats")
                        elif not len(exp_psf_stats_str) and "exp_psf_stats = {" in ln:
                            exp_psf_stats_str = ln.lstrip("exp_psf_stats = ").rstrip("\n")
                        elif len(exp_psf_stats_str) and not exp_psf_stats_str_end:
                            exp_psf_stats_str += ln.rstrip("\n").rstrip(" ").replace("\\", '')
                            if ln.rstrip("\n").rstrip().endswith("}"):
                                exp_psf_stats_str_end = True
                                print("Stop reading exp_psf_stats")
                        elif not len(exp_model_stats_str) and "exp_model_stats = {" in ln:
                            exp_model_stats_str = ln.lstrip("exp_model_stats = ").rstrip("\n")
                        elif len(exp_model_stats_str) and not exp_model_stats_str_end:
                            exp_model_stats_str += ln.rstrip("\n").rstrip(" ").replace("\\", '')
                            if ln.rstrip("\n").rstrip().endswith("}"):
                                exp_model_stats_str_end = True
                                print("Stop reading exp_model_stats")
                        elif not len(exp_resid_stats_str) and "exp_resid_stats = {" in ln:
                            exp_resid_stats_str = ln.lstrip("exp_resid_stats = ").rstrip("\n")
                        elif len(exp_resid_stats_str) and not exp_resid_stats_str_end:
                            exp_resid_stats_str += ln.rstrip("\n").rstrip(" ").replace("\\", '')
                            if ln.rstrip("\n").rstrip().endswith("}"):
                                exp_resid_stats_str_end = True
                                print("Stop reading exp_resid_stats")
                        elif not len(exp_sumwt_stats_str) and "exp_sumwt_stats = {" in ln:
                            exp_sumwt_stats_str = ln.lstrip("exp_sumwt_stats = ").rstrip("\n")
                        elif len(exp_sumwt_stats_str) and not exp_sumwt_stats_str_end:
                            exp_sumwt_stats_str += ln.rstrip("\n").rstrip(" ").replace("\\", '')
                            if ln.rstrip("\n").rstrip().endswith("}"):
                                exp_sumwt_stats_str_end = True
                                print("Stop reading exp_sumwt_stats")
                        elif not len(exp_wt_stats_str) and "exp_wt_stats = {" in ln:
                            exp_wt_stats_str = ln.lstrip("exp_wt_stats = ").rstrip("\n")
                        elif len(exp_wt_stats_str) and not exp_wt_stats_str_end:
                            exp_wt_stats_str += ln.rstrip("\n").rstrip(" ").replace("\\", '')
                            if ln.rstrip("\n").rstrip().endswith("}"):
                                exp_wt_stats_str_end = True
                                print("Stop reading exp_wt_stats")
                        elif not len(exp_bmin_dict_str) and "exp_bmin_dict = {" in ln:
                            exp_bmin_dict_str = ln.lstrip("exp_bmin_dict = ").rstrip("\n")
                            # elif len(exp_bmin_dict_str) and not exp_bmin_dict_str_end:
                            #    exp_bmin_dict_str+=ln.rstrip("\n").rstrip(" ").replace("\\",'')
                            if ln.rstrip("\n").rstrip().endswith("}"):
                                print("Stop reading exp_bmin_dict")
                        elif not len(exp_bmaj_dict_str) and "exp_bmaj_dict = {" in ln:
                            exp_bmaj_dict_str = ln.lstrip("exp_bmaj_dict = ").rstrip("\n")
                            # elif len(exp_bmaj_dict_str) and not exp_bmaj_dict_str_end:
                            #    exp_bmaj_dict_str+=ln.rstrip("\n").rstrip(" ").replace("\\",'')
                            if ln.rstrip("\n").rstrip().endswith("}"):
                                print("Stop reading exp_bmaj_dict")
                        elif not len(exp_pa_dict_str) and "exp_pa_dict = {" in ln:
                            exp_pa_dict_str = ln.lstrip("exp_pa_dict = ").rstrip("\n")
                            # elif len(exp_pa_dict_str) and not exp_pa_dict_str_end:
                            #    exp_pa_dict_str+=ln.rstrip("\n").rstrip(" ").replace("\\",'')
                            if ln.rstrip("\n").rstrip().endswith("}"):
                                print("Stop reading exp_pa_dict")
                        elif not len(exp_im1_stats_str) and "exp_im1_stats = {" in ln:
                            exp_im1_stats_str = ln.lstrip("exp_im1_stats = ").rstrip("\n")
                        elif len(exp_im1_stats_str) and not exp_im1_stats_str_end:
                            exp_im1_stats_str += ln.rstrip("\n").rstrip(" ").replace("\\", '')
                            if ln.rstrip("\n").rstrip().endswith("}"):
                                exp_im1_stats_str_end = True
                                print("Stop reading exp_im1_stats")
                        elif not len(exp_model1_stats_str) and "exp_model1_stats = {" in ln:
                            exp_model1_stats_str = ln.lstrip("exp_model1_stats = ").rstrip("\n")
                        elif len(exp_model1_stats_str) and not exp_model1_stats_str_end:
                            exp_model1_stats_str += ln.rstrip("\n").rstrip(" ").replace("\\", '')
                            if ln.rstrip("\n").rstrip().endswith("}"):
                                exp_model1_stats_str_end = True
                                print("Stop reading exp_model1_stats")
                        elif not len(exp_resid1_stats_str) and "exp_resid1_stats = {" in ln:
                            exp_resid1_stats_str = ln.lstrip("exp_resid1_stats = ").rstrip("\n")
                        elif len(exp_resid1_stats_str) and not exp_resid1_stats_str_end:
                            exp_resid1_stats_str += ln.rstrip("\n").rstrip(" ").replace("\\", '')
                            if ln.rstrip("\n").rstrip().endswith("}"):
                                exp_resid1_stats_str_end = True
                                print("Stop reading exp_resid1_stats")
                        elif not len(exp_sumwt1_stats_str) and "exp_sumwt1_stats = {" in ln:
                            exp_sumwt1_stats_str = ln.lstrip("exp_sumwt1_stats = ").rstrip("\n")
                        elif len(exp_sumwt1_stats_str) and not exp_sumwt1_stats_str_end:
                            exp_sumwt1_stats_str += ln.rstrip("\n").rstrip(" ").replace("\\", '')
                            if ln.rstrip("\n").rstrip().endswith("}"):
                                exp_sumwt1_stats_str_end = True
                                print("Stop reading exp_sumwt1_stats")

            # if exp_im_stats_str_end and exp_mask_stats_str_end and exp_pb_stats_str_end:
            #    print("finish reading the test")
            #    readDict=False
            #    break

        outdict = {'exp_im_stats': ast.literal_eval(exp_im_stats_str),
                   'exp_mask_stats': ast.literal_eval(exp_mask_stats_str),
                   'exp_pb_stats': ast.literal_eval(exp_pb_stats_str),
                   'exp_psf_stats': ast.literal_eval(exp_psf_stats_str),
                   'exp_model_stats': ast.literal_eval(exp_model_stats_str),
                   'exp_resid_stats': ast.literal_eval(exp_resid_stats_str),
                   'exp_sumwt_stats': ast.literal_eval(exp_sumwt_stats_str)}
        if 'mosaic' in testname:
            outdict['exp_wt_stats'] = ast.literal_eval(exp_wt_stats_str)
        if 'cube' in testname and not ('eph' in testname):
            outdict['exp_bmin_dict'] = ast.literal_eval(exp_bmin_dict_str)
            outdict['exp_bmaj_dict'] = ast.literal_eval(exp_bmaj_dict_str)
            outdict['exp_pa_dict'] = ast.literal_eval(exp_pa_dict_str)
        if 'mtmfs' in testname:
            outdict['exp_im1_stats'] = ast.literal_eval(exp_im1_stats_str)
            outdict['exp_model1_stats'] = ast.literal_eval(exp_model1_stats_str)
            outdict['exp_resid1_stats'] = ast.literal_eval(exp_resid1_stats_str)
            outdict['exp_sumwt1_stats'] = ast.literal_eval(exp_sumwt1_stats_str)

        outdictwithTestname = {testname: outdict}
        with open(outfile, 'w') as outf:
            json.dump(outdictwithTestname, outf, indent=4)
        print("Extracted the fudicial value dictionaries is saved in a file, ", outfile)


def create_expdict_jsonfile(inmetricsfile, templatemetrics, outmetricsfile):
    """
    create the fiducial metric dictionaries
    from the corresponding metric dictionaries of the current
    run saved as a json file
    """
    infiles = [inmetricsfile, templatemetrics]
    for f in infiles:
        if not os.path.exists(f):
            raise Exception("%s is not found" % f)
    outDict = {}
    isbeaminfo = False
    with open(templatemetrics, 'r') as tmplf, open(inmetricsfile, 'r') as curf, open(outmetricsfile, 'w') as outf:
        tmplFidDict = json.load(tmplf)
        curDict = json.load(curf)
        testname = list(curDict.keys())[0]
        if 'test_' in testname:
            if testname in curDict:
                outDict[testname] = {}
                subOutDict = outDict[testname]
                for expkey in tmplFidDict[testname]:
                    print("Processing expkey=", expkey)
                    curkey = expkey[4:]  # name of the dict in the current metrics
                    print("curkey before mode=", curkey)
                    if not curkey.endswith("_dict"):
                        curkey += "_dict"
                    else:
                        isbeaminfo = True
                    print("Processing curkey=", curkey)
                    print("isbeaminfo=", isbeaminfo)
                    if curkey in curDict[testname]:
                        subOutDict[expkey] = {}
                        if isbeaminfo:
                            subOutDict[expkey] = curDict[testname][curkey]
                        else:
                            # loop through each metric inside the particular exp_ dict
                            for metrickey in tmplFidDict[testname][expkey]:
                                # check to see if the metric exist in the input (current) metric dict
                                if metrickey in curDict[testname][curkey]:
                                    print("expkey={}, metrickey={}, metricbool={}".format(expkey, metrickey,
                                                                                          tmplFidDict[testname][expkey][
                                                                                              metrickey][0]))
                                    subOutDict[expkey][metrickey] = [tmplFidDict[testname][expkey][metrickey][0],
                                                                     curDict[testname][curkey][metrickey]]
                                else:  # metric does not exist in input cur metric dicts
                                    raise Exception(
                                        "Missing the metric key={} in {}. Check the input file".format(metrickey,
                                                                                                       inmetricsfile))
                    else:
                        if curkey == 'bmin_dict' or curkey == 'bmaj_dict' or curkey == 'pa_dict':
                            print("Missing key={} in {}. The input json is probably made from serial run"
                                  .format(curkey, inmetricsfile))
                        else:
                            raise Exception("Missing key={} in {}. Check the input file".format(curkey, inmetricsfile))

                json.dump(outDict, outf)
            else:
                raise Exception(
                    "{} does not contain test name {} as a top level key. Please modify the input file".format(
                        inmetricsfile, testname))
        else:
            raise Exception(
                "{} does not contain a test name in the top key. Please modify the input file".format(tmplFidDict))


def create_combined_expdict_jsonfile(jsonlist, outjson, casaversion):
    """
    make a single json to be used in alma stakeholder tests
    """
    import copy
    outdict = {'casa_version': casaversion}
    with open(outjson, 'w') as outf:
        for jsonfile in jsonlist:
            try:
                with open(jsonfile, 'r') as inf:
                    testdict = json.load(inf)
                    if type(testdict) is dict:
                        testname = list(testdict.keys())[0]
                        if 'test_' in testname:
                            # check the version info of the test results
                            if 'casa_version' in testdict[testname]:
                                casaversionused = testdict['casa_version']
                                if casaversionused != casaversion:
                                    print('casa_version for {} is different from expected casaversion: {}'
                                          .format(casaversionused, casaversion))
                            else:
                                print("No casa_version info in the input dictionary. Skip the check.")
                            outdict[testname] = copy.deepcopy(testdict[testname])
            except RuntimeError:
                print("Errors in reading the input json file. Check the input")
        if outdict != {}:
            json.dump(outdict, outf)
        else:
            print("Error occured. No outfile is wriiten.")


def read_expdict_jsonfile(jsonfilename=None):
    """read the json file containing the fiducial metrics parameter values for all tests"""
    try:
        with open(jsonfilename, 'r') as fexp:
            return json.load(fexp)
    except RuntimeError:
        print("Error occurred in reading the json file for fiducial values")


def read_testcase_expdicts(jsonfilename, testcasename, version):
    try:
        with open(jsonfilename, 'r') as fexp:
            alltestdicts = json.load(fexp)
            # check a CASA version that exp_dicts based on
            if version != '':
                if 'casa_version' in alltestdicts.keys():
                    if version != alltestdicts['casa_version']:
                        raise SystemError('Mismatch in the fiducial data file version. The testcase expects fiducial '
                                          'values based on the CASA {} '.format(version))
            if testcasename in alltestdicts:
                return alltestdicts[testcasename]
            else:
                raise Exception("key {} is not found in {} exp_dicts".format(testcasename, jsonfilename))
    except RuntimeError:
        print("Error occurred in reading the json file for fiducial values")


def update_expdict_jsonfile(newexpdictlist, jsonfilename):
    """Convert current metrics parameter values stored in json per test to the exp_dict json.
       The output will be an updated json file named jsonfilename+"_update"

       newexpdictlist: a list of name of json files, which are produced by running
       save_to_dict() inside each stakeholder testcase or can be produced. The list can
       consist of only the testcases that need to be updated and the exp_dicts for other
       testcases not in the list won't be modified and copy to the new json as is.

       jsonfilename: current json file conta:w!ins all the fiducial metrics values

    """
    import copy
    tmplFidDict = read_expdict_jsonfile(jsonfilename)
    newjsonfile = jsonfilename.split(".json")[0] + "_update.json"
    with open(newjsonfile, 'w') as outf:
        outDict = copy.deepcopy(tmplFidDict)
        isbeaminfo = False
        for inmetricsfile in newexpdictlist:  # read current values for each testscase
            with open(inmetricsfile, 'r') as curf:
                curDict = json.load(curf)
                testname = list(curDict.keys())[0]
                if 'test_' in testname:
                    if testname in tmplFidDict:
                        outDict[testname] = {}
                        subOutDict = outDict[testname]
                        for expkey in tmplFidDict[testname]:
                            print("Processing expkey=", expkey)
                            curkey = expkey[4:]  # name of the dict in the current metrics
                            print("curkey before mode=", curkey)
                            if not curkey.endswith("_dict"):
                                curkey += "_dict"
                            else:
                                isbeaminfo = True
                            print("Processing curkey=", curkey)
                            print("isbeaminfo=", isbeaminfo)
                            if curkey in curDict[testname]:
                                subOutDict[expkey] = {}
                                if isbeaminfo:
                                    subOutDict[expkey] = curDict[testname][curkey]
                                else:
                                    # loop through each metric inside the particular exp_ dict
                                    for metrickey in tmplFidDict[testname][expkey]:
                                        # check to see if the metric exist in the input (current) metric dict
                                        if metrickey in curDict[testname][curkey]:
                                            print("expkey={}, metrickey={}, metricbool={}".format(expkey, metrickey,
                                                                                                  tmplFidDict[testname][
                                                                                                      expkey][
                                                                                                      metrickey][0]))
                                            subOutDict[expkey][metrickey] = [
                                                tmplFidDict[testname][expkey][metrickey][0],
                                                curDict[testname][curkey][metrickey]]
                                        else:  # metric does not exist in input cur metric dicts
                                            raise Exception("Missing the metric key={} in {}. Check the input file".
                                                            format(metrickey, inmetricsfile))
                            else:
                                if curkey == 'bmin_dict' or curkey == 'bmaj_dict' or curkey == 'pa_dict':
                                    print("Missing key={} in {}. The input json is probably made from serial run"
                                          .format(curkey, inmetricsfile))
                                else:
                                    raise Exception(
                                        "Missing key={} in {}. Check the input file".format(curkey, inmetricsfile))

                        json.dump(outDict, outf)
                    else:
                        raise Exception("{} does not contain test name {} as a top level key." +
                                        "Please modify the input file".format(inmetricsfile, testname))
                else:
                    raise Exception("{} does not contain a test name in the top key." +
                                    "Please modify the input file".format(tmplFidDict))


def compare_expdictjson(newjson, oldjson):
    """
    compare the two exp dicts - used to check updating of exp_dicts json file
    is done properly...
    """
    with open(oldjson, 'r') as fold, open(newjson, 'r') as fnew:
        newdict = json.load(fnew)
        olddict = json.load(fold)

        if newdict == olddict:
            return "The two json files are indentical"
        else:
            # level 0 (testcase level)
            newkey0list = list(newdict.keys())
            oldkey0list = list(olddict.keys())
            newonlykey0 = set(newkey0list).difference(oldkey0list)
            oldonlykey0 = set(oldkey0list).difference(newkey0list)
            commonkey0 = set(newkey0list).intersection(oldkey0list)

            finaldiffdict = {}
            for key0 in commonkey0:
                # do for each testcase, extract set of metrics  for each image type
                if type(newdict[key0]) == dict and type(olddict[key0]) == dict:
                    newkey1list = list(newdict[key0].keys())
                    oldkey1list = list(olddict[key0].keys())
                    newonlykey1 = set(newkey1list).difference(oldkey1list)
                    oldonlykey1 = set(oldkey1list).difference(newkey1list)
                    commonkey1 = set(newkey1list).intersection(oldkey1list)

                    for key1 in commonkey1:
                        newkey2list = list(newdict[key0][key1].keys())
                        oldkey2list = list(olddict[key0][key1].keys())
                        newonlykey2 = set(newkey2list).difference(oldkey2list)
                        oldonlykey2 = set(oldkey2list).difference(newkey2list)
                        commonkey2 = set(newkey2list).intersection(oldkey2list)
                        # Comparison of the values and test threshold type
                        diffkey2dict = {}
                        for key2 in commonkey2:  # a list containing [thres. type, []]
                            if type(olddict[key0][key1][key2]) == list:
                                if olddict[key0][key1][key2][1] != newdict[key0][key1][key2][1]:
                                    if key2 not in diffkey2dict:
                                        diffkey2dict[key2] = {}
                                    diffkey2dict[key2]['msg'] = 'diff in value(s)'
                                    diffkey2dict[key2]['json1'] = newdict[key0][key1][key2]
                                    diffkey2dict[key2]['json2'] = olddict[key0][key1][key2]
                                if olddict[key0][key1][key2][0] != newdict[key0][key1][key2][0]:
                                    if 'msg' in diffkey2dict[key2]:
                                        diffkey2dict[key2]['msg'] += ' and threshold type'
                                    else:
                                        diffkey2dict[key2]['msg'] = 'diff in threshold type'
                                        diffkey2dict[key2]['json1'] = newdict[key0][key1][key2][0]
                                        diffkey2dict[key2]['json2'] = olddict[key0][key1][key2][0]
                            else:
                                # non metric values (possibly comments or version info)
                                if olddict[key0][key1][key2] != newdict[key0][key1][key2]:
                                    diffkey2dict[key2] = \
                                        'json1: {}, json2: {}'.format(newdict[key0][key1][key2],
                                                                      olddict[key0][key1][key2])

                        if diffkey2dict != dict():
                            if key0 not in finaldiffdict:
                                finaldiffdict[key0] = {}
                            if key1 not in finaldiffdict:
                                finaldiffdict[key0][key1] = {}
                            finaldiffdict[key0][key1] = copy.deepcopy(diffkey2dict)
                        if newonlykey2 != set():
                            finaldiffdict[key0][key1]['metric keys only in json1'] = newonlykey2
                        if oldonlykey2 != set():
                            finaldiffdict[key0][key1]['metric keys only in json2'] = oldonlykey2
                    if newonlykey1 != set():
                        finaldiffdict[key0]['matric dict only in json1'] = newonlykey1
                    if oldonlykey1 != set():
                        finaldiffdict[key0]['matric dict only in json2'] = oldonlykey1
                else:
                    if newdict[key0] != olddict[key0]:
                        finaldiffdict[key0] = 'diff info json1: {}, json2: {}'.format(newdict[key0], olddict[key0])
            if newonlykey0 != set():
                finaldiffdict['testcase(s)/other info only in json1'] = newonlykey0
            if oldonlykey0 != set():
                finaldiffdict['testcase(s)/other info  only in json2'] = oldonlykey0
        return finaldiffdict
