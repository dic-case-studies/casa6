import shutil
import os
import inspect
import unittest
from collections import OrderedDict

from casatools import quanta, ctsys
from casatasks import flagdata
from casatasks.private import flaghelper as fh

qa = quanta()

# Path to data
datapath = ctsys.resolve("unittest/flagdata/")
#
# Test of flaghelper.py
#
"""
Test the following functions:
readFile
readFiles
readAndParse
parseDictionary
applyTimeBuffer
writeFlagCommands
parseAgents
"""


def _test_eq(result, total, flagged):

    print(
        "%s of %s data was flagged, expected %s of %s"
        % (result["flagged"], result["total"], flagged, total)
    )
    assert result["total"] == total, "%s data in total; %s expected" % (
        result["total"],
        total,
    )
    assert result["flagged"] == flagged, "%s flags set; %s expected" % (
        result["flagged"],
        flagged,
    )


def create_input(str_text, filename):
    """Save the string in a text file"""

    inp = filename
    cmd = str_text

    # remove file first
    if os.path.exists(inp):
        os.system("rm -f " + inp)

    # save to a file
    with open(inp, "w") as f:
        f.write(cmd)

    f.close()

    return


# Base class which defines setUp functions
# for importing different data sets
class test_base(unittest.TestCase):
    def setUp_onlineFlags(self):
        """Large file with online flags"""
        self.inpfile = "BigOnlineFlags.txt"

        if not os.path.exists(self.inpfile):
            os.system("cp -RH " + datapath + self.inpfile + " " + self.inpfile)

    def setUp_mockasdm(self):
        """Mock ASDM containing only XML files"""
        self.inpfile = "uid___MockASDM"

        if not os.path.exists(self.inpfile):
            os.system("cp -RH " + datapath + self.inpfile + " " + self.inpfile)


class test_flaghelper(test_base):
    @classmethod
    def tearDownClass(cls) -> None:
        os.system("rm -rf flaghelper*.txt")
        shutil.rmtree("uid___MockASDM", ignore_errors=True)

    def test_readFile1(self):
        """flaghelper: read a file from disk"""
        # creat input file
        myinput = (
            "scan='1~3' mode='manual'\n" + "scan='5' mode='manualflag'\n" "#scan='4'"
        )
        filename = "flaghelper1.txt"
        create_input(myinput, filename)

        alist = fh.readFile(filename)
        self.assertEqual(len(alist), 2)

    def test_readFile2(self):
        """flaghelper: read a file from disk with spaces in command"""
        # creat input file
        myinput = (
            "scan='1~3' mode='manual'\n" + "scan=' 5'      mode='manualflag'\n"
            "antenna='DV04 &&*'"
        )
        filename = "flaghelper2.txt"
        create_input(myinput, filename)

        alist = fh.readFile(filename)
        self.assertEqual(len(alist), 3)

    def test_readFiles(self):
        """flaghelper: read several files from disk"""
        myinput = "scan='1'\n" "scan='2'\n" "# a comment line\n" "empty_line"
        filename1 = "flaghelper3a.txt"
        create_input(myinput, filename1)

        # Create second input file
        myinput = "scan='5'\n" " \n" "scan='6'\n" "scan='7'"
        filename2 = "flaghelper3b.txt"
        create_input(myinput, filename2)

        # Create third input file
        myinput = "scan='4' mode='clip' clipminmax=[0,4]"
        filename3 = "flaghelper3c.txt"
        create_input(myinput, filename3)

        alist = fh.readFiles([filename1, filename2, filename3])
        self.assertEqual(len(alist), 7)

        astring = alist[2]
        self.assertEqual(astring, "empty_line")

    def test_readAndParse(self):
        """flaghelper: compare the read and parse from a file and from a list of strings"""
        # <startTime>4891227930515540000 <endTime>4891227932453838000
        # <startTime>4891228473545856000 <endTime>4891228473731891000
        # <startTime>4891226924455911000 <endTime>4891226927502314000
        # <startTime>4891228838164987000 <endTime>4891228838418996000
        # <startTime>4891228609440808000 <endTime>4891228612489617000

        online = [
            "antenna='DV03&&*' timerange='2013/11/15/10:25:30.516~2013/11/15/10:25:32.454'",
            "antenna='DA44&&*' timerange='2013/11/15/10:34:33.546~2013/11/15/10:34:33.732'",
            "antenna='DA46&&*' timerange='2013/11/15/10:08:44.456~2013/11/15/10:08:47.502'",
            "antenna='DV09&&*' timerange='2013/11/15/10:18:11.798~2013/11/15/10:18:13.837'",
            "antenna='DV05&&*' timerange='2013/11/15/10:40:38.165~2013/11/15/10:40:38.419'",
        ]

        myinput = (
            "antenna='DV03&&*' timerange='2013/11/15/10:25:30.516~2013/11/15/10:25:32.454'\n"
            "antenna='DA44&&*' timerange='2013/11/15/10:34:33.546~2013/11/15/10:34:33.732'\n"
            "antenna='DA46&&*' timerange='2013/11/15/10:08:44.456~2013/11/15/10:08:47.502'\n"
            "antenna='DV09&&*' timerange='2013/11/15/10:18:11.798~2013/11/15/10:18:13.837'\n"
            "antenna='DV05&&*' timerange='2013/11/15/10:40:38.165~2013/11/15/10:40:38.419'"
        )

        filename1 = "flaghelperonline1.txt"
        create_input(myinput, filename1)

        dlist1 = fh.readAndParse([filename1])
        self.assertEqual(len(dlist1), 5)

        # Use the list instead of the file
        dlist2 = fh.readAndParse(online)

        self.assertListEqual(dlist1, dlist2)

        # Compare with the original Flag.xml, second row
        orig_time_start = float(4891228473545856000) * 1.0e-9
        orig_time_end = float(4891228473731891000) * 1.0e-9

        proc_time = dlist2[1]["timerange"]
        t0, t1 = proc_time.split("~", 1)
        startTime = qa.totime(t0)["value"]
        startTimeSec = float(startTime * 24 * 3600)
        endTime = qa.totime(t1)["value"]
        endTimeSec = float(endTime * 24 * 3600)

        self.assertAlmostEqual(orig_time_start, startTimeSec, places=3)
        self.assertAlmostEqual(orig_time_end, endTimeSec, places=3)

    def test_readAndParseTbuff(self):
        """flaghelper: compare the read and parse and apply tbuff"""
        # MJD in seconds of timeranges are these
        # <startTime>4891227930515540000 <endTime>4891227932453838000
        # <startTime>4891228473545856000 <endTime>4891228473731891000
        # <startTime>4891226924455911000 <endTime>4891226927502314000
        # <startTime>4891228838164987000 <endTime>4891228838418996000
        # <startTime>4891228609440808000 <endTime>4891228612489617000

        online = [
            "antenna='DV03&&*' timerange='2013/11/15/10:25:30.516~2013/11/15/10:25:32.454'",
            "antenna='DA44&&*' timerange='2013/11/15/10:34:33.546~2013/11/15/10:34:33.732'",
            "antenna='DA46&&*' timerange='2013/11/15/10:08:44.456~2013/11/15/10:08:47.502'",
            "antenna='DV09&&*' timerange='2013/11/15/10:18:11.798~2013/11/15/10:18:13.837'",
            "antenna='DV05&&*' timerange='2013/11/15/10:40:38.165~2013/11/15/10:40:38.419'",
        ]

        myinput = (
            "antenna='DV03&&*' timerange='2013/11/15/10:25:30.516~2013/11/15/10:25:32.454'\n"
            "antenna='DA44&&*' timerange='2013/11/15/10:34:33.546~2013/11/15/10:34:33.732'\n"
            "antenna='DA46&&*' timerange='2013/11/15/10:08:44.456~2013/11/15/10:08:47.502'\n"
            "antenna='DV09&&*' timerange='2013/11/15/10:18:11.798~2013/11/15/10:18:13.837'\n"
            "antenna='DV05&&*' timerange='2013/11/15/10:40:38.165~2013/11/15/10:40:38.419'"
        )

        filename1 = "flaghelperonline2.txt"
        create_input(myinput, filename1)

        # First timerange from online before padding
        origt = timerange = "2013/11/15/10:25:30.516~2013/11/15/10:25:32.454"

        # Apply tbuff to timeranges
        timebuffer = 1.1
        dlist1 = fh.readAndParse([filename1], tbuff=timebuffer)
        self.assertEqual(len(dlist1), 5)

        # Get the first padded timerange from output
        padt = dlist1[0]["timerange"]

        # Revert the tbuff application manually
        t0, t1 = padt.split("~", 1)
        startTime = qa.totime(t0)["value"]
        startTimeSec = float((startTime * 24 * 3600) + timebuffer)
        startTimeSec = qa.quantity(startTimeSec, "s")
        paddedT0 = qa.time(startTimeSec, form="ymd", prec=9)[0]
        # end time
        endTime = qa.totime(t1)["value"]
        endTimeSec = float((endTime * 24 * 3600) - timebuffer)
        endTimeSec = qa.quantity(endTimeSec, "s")
        paddedT1 = qa.time(endTimeSec, form="ymd", prec=9)[0]

        newtimerange = paddedT0 + "~" + paddedT1

        # Compare with the original
        self.assertEqual(origt, newtimerange)

        # Compare with original values from Flag.xml
        xmlt0 = float(4891227930515540000) * 1.0e-9
        xmlt1 = float(4891227932453838000) * 1.0e-9

        self.assertAlmostEqual(xmlt0, startTimeSec["value"], places=3)
        self.assertAlmostEqual(xmlt1, endTimeSec["value"], places=3)

    def test_readAndParseIrregularTbuff(self):
        """flaghelper: compare the read and parse and apply of irregular tbuff"""
        # MJD in seconds of timeranges are these
        # <startTime>4891227930515540000 <endTime>4891227932453838000
        # <startTime>4891228473545856000 <endTime>4891228473731891000
        # <startTime>4891226924455911000 <endTime>4891226927502314000
        # <startTime>4891228838164987000 <endTime>4891228838418996000
        # <startTime>4891228609440808000 <endTime>4891228612489617000

        online = [
            "antenna='DV03&&*' timerange='2013/11/15/10:25:30.516~2013/11/15/10:25:32.454'",
            "antenna='DA44&&*' timerange='2013/11/15/10:34:33.546~2013/11/15/10:34:33.732'",
            "antenna='DA46&&*' timerange='2013/11/15/10:08:44.456~2013/11/15/10:08:47.502'",
            "antenna='DV09&&*' timerange='2013/11/15/10:18:11.798~2013/11/15/10:18:13.837'",
            "antenna='DV05&&*' timerange='2013/11/15/10:40:38.165~2013/11/15/10:40:38.419'",
        ]

        myinput = (
            "antenna='DV03&&*' timerange='2013/11/15/10:25:30.516~2013/11/15/10:25:32.454'\n"
            "antenna='DA44&&*' timerange='2013/11/15/10:34:33.546~2013/11/15/10:34:33.732'\n"
            "antenna='DA46&&*' timerange='2013/11/15/10:08:44.456~2013/11/15/10:08:47.502'\n"
            "antenna='DV09&&*' timerange='2013/11/15/10:18:11.798~2013/11/15/10:18:13.837'\n"
            "antenna='DV05&&*' timerange='2013/11/15/10:40:38.165~2013/11/15/10:40:38.419'"
        )

        filename1 = "flaghelperonline2.txt"
        create_input(myinput, filename1)

        # timeranges from online before padding, for comparison later
        timeranges = []
        for cmd in online:
            a, b = cmd.split(" ")
            b = b.lstrip("timerange=")
            timeranges.append(b.strip("'"))

        # Apply 2 values of tbuff to timeranges
        timebuffer = [0.4, 0.7]
        dlist1 = fh.readAndParse([filename1], tbuff=timebuffer)
        self.assertEqual(len(dlist1), 5)

        for n, cmd in enumerate(dlist1):
            padt = cmd["timerange"]
            #        padt = dlist1[0]['timerange']

            # Revert the tbuff application manually
            t0, t1 = padt.split("~", 1)
            startTime = qa.totime(t0)["value"]
            startTimeSec = float((startTime * 24 * 3600) + timebuffer[0])
            startTimeSec = qa.quantity(startTimeSec, "s")
            paddedT0 = qa.time(startTimeSec, form="ymd", prec=9)[0]
            # end time
            endTime = qa.totime(t1)["value"]
            endTimeSec = float((endTime * 24 * 3600) - timebuffer[1])
            endTimeSec = qa.quantity(endTimeSec, "s")
            paddedT1 = qa.time(endTimeSec, form="ymd", prec=9)[0]

            newtimerange = paddedT0 + "~" + paddedT1

            # Compare with the original
            self.assertEqual(timeranges[n], newtimerange)

    def test_parseDictionary(self):
        """flaghelper: read a file and parse to a dictionary"""
        # creat input file
        myinput = (
            "scan='1~3' mode='manual'\n"
            "scan=' 5'      mode='manualflag'\n"
            "antenna='DV04 &&*'\n"
            "#mode=shadow\n"
            "antenna='DV01&&*' timerange='2013/11/15/10:25:30.516~2013/11/15/10:25:32.454' reason='ACS_not in place'"
        )
        filename = "flaghelper4.txt"
        create_input(myinput, filename)

        alist = fh.readFile(filename)
        self.assertEqual(len(alist), 4)

        adict = fh.parseDictionary(alist)
        self.assertEqual(
            adict[3]["command"]["timerange"],
            "2013/11/15/10:25:30.516~2013/11/15/10:25:32.454",
        )
        self.assertEqual(adict[2]["command"]["antenna"], "DV04 &&*")

    def test_parseNoEval1(self):
        """flaghelper: CAS-6553 parse and evaluate a string with double whitespaces"""

        # Dividers for string
        first = " "
        second = "="
        reference = OrderedDict(
            [("mode", "extend"), ("antenna", "ea24"), ("flagnearfreq", True)]
        )

        # cmd with whitespace between pairs
        cmd = "mode='extend' antenna='ea24'  flagnearfreq=True"
        myparser = fh.Parser(first, second)

        res = myparser.parseNoEval(cmd)

        # evaluate parameters to fix single quote left dangling
        resdict = fh.evaluateParameters(res)

        self.assertDictEqual(
            reference,
            resdict,
            "Failed to parserNoEval() with whitespaces between pairs",
        )

    def test_parseNoEval2(self):
        """flaghelper: CAS-6553 parse and evaluate a string with whitespaces"""

        # Dividers for string
        first = " "
        second = "="
        reference = OrderedDict(
            [("mode", "extend"), ("antenna", "ea24"), ("flagnearfreq", True)]
        )

        # cmd with whitespace between pairs and at the end
        cmd = "mode='extend' antenna='ea24 '  flagnearfreq=True"
        myparser = fh.Parser(first, second)

        # test parseNoEval()
        res = myparser.parseNoEval(cmd)

        # evaluate parameters to fix single quote
        resdict = fh.evaluateParameters(res)
        self.assertDictEqual(
            reference, resdict, "Failed to parserNoEval() with whitespaces in value"
        )

    def test_parseNoEval3(self):
        """flaghelper: CAS-6553 parse and evaluate a string with whitespaces"""

        # Dividers for string
        first = " "
        second = "="
        reference = OrderedDict([("mode", "tfcrop"), ("antenna", "ea24")])

        # cmd with whitespace in the begin and end
        cmd = " mode='tfcrop' antenna='ea24' "
        myparser = fh.Parser(first, second)

        res = myparser.parseNoEval(cmd)

        # evaluate parameters to fix single quote left dangling
        resdict = fh.evaluateParameters(res)
        self.assertDictEqual(
            reference,
            resdict,
            "Failed to parserNoEval() with whitespaces at begin and end",
        )

    def test_evaluateParameters1(self):
        """flaghelper: parse and evaluate a string with a many extra whitespaces"""

        # Dividers for string
        first = " "
        second = "="
        reference = OrderedDict(
            [
                ("mode", "manual"),
                ("antenna", "ea24"),
                ("spw", "0"),
                ("reason", "MY WHITESPACES"),
            ]
        )

        # cmd with single quote inside a string
        cmd = " mode='manual'   antenna='ea24'      spw='0'   reason='MY WHITESPACES'"
        myparser = fh.Parser(first, second)

        res = myparser.parseNoEval(cmd)

        # evaluate parameters to fix single quote
        resdict = fh.evaluateParameters(res)
        self.assertDictEqual(
            reference, resdict, "Failed to evaluateParameters with many whitespaces"
        )

    @unittest.skip("CAS-6553 breaks this use-case.")
    def test_evaluateParameters2(self):
        """flaghelper: parse and evaluate a string with a single quote inside"""

        # Dividers for string
        first = " "
        second = "="
        reference = OrderedDict([("mode", "manual"), ("field", "It A'int a Field")])

        # cmd with single quote inside a string
        cmd = "mode='manual' field='It A'int a Field'"
        myparser = fh.Parser(first, second)

        res = myparser.parseNoEval(cmd)

        # evaluate parameters to fix single quote
        resdict = fh.evaluateParameters(res)
        self.assertDictEqual(
            reference,
            resdict,
            "Failed to evaluateParameters with single quote in value",
        )

    def test_evaluateFlagParameters1(self):
        """flaghelper: evaluate non-existing flagdata parameters"""
        cmd = [
            "antenna='''BK07"
            "'asdf timesrange='2013/01/31/08:09:55.248~2013/01/31/08:10:01.296' reason='quack'"
        ]

        # Get the defaults of each parameter
        fparams = dict(inspect.signature(flagdata).parameters.items())
        mydict = fh.parseDictionary(cmd)
        try:
            res = fh.evaluateFlagParameters(mydict, fparams)
        except ValueError as instance:
            print("Expected error: %s" % instance)
        else:
            self.assertNotEqual(
                res, True, "Some parameter in cmd should not evaluate True"
            )

    def test_evaluateFlagParameters3(self):
        """flaghelper: evaluate correct parameter types"""
        fparams = dict(inspect.signature(flagdata).parameters.items())

        cmd = ["antenna='DA31'", "mode='clip' clipzeros=True"]
        mydict = fh.parseDictionary(cmd)
        self.assertTrue(fh.evaluateFlagParameters(mydict, fparams))

    def test_evaluateFlagParameters4(self):
        """flaghelper: evaluate new quackinterval parameter type"""
        cmd = ["mode='quack' quackinterval=2", "mode='quack' quackinterval=3.1"]

        fparams = dict(inspect.signature(flagdata).parameters.items())

        mydict = fh.parseDictionary(cmd)

        self.assertTrue(fh.evaluateFlagParameters(mydict, fparams))

    def test_evaluateFlagParameters5(self):
        """flaghelper: evaluate wrong rflag parameter type"""
        cmd = [
            "antenna='''BK07"
            "'asdf timerange='2013/01/31/08:09:55.248~2013/01/31/08:10:01.296' reason='quack'",
            "mode='rflag' timedev={'threshold':3.0} freqdev=5.0",
        ]

        # Need to cast to dict because odict was giving erros from copy.deepcopy inside evaluateFlagParameters
        fparams = dict(inspect.signature(flagdata).parameters.items())
        mydict = fh.parseDictionary(cmd)

        try:
            res = fh.evaluateFlagParameters(mydict, fparams)
        except Exception as instance:
            print("Expected error: %s" % instance)
        else:
            self.assertNotEqual(
                res, True, "Some parameter in cmd should not evaluate True"
            )

    def test_evaluateDataSelectionParameters(self):
        """flaghelper: evaluate wrong type of data selection parameters"""
        fparams = dict(inspect.signature(flagdata).parameters.items())

        cmd = ["field=2 spw=1 observation=0", "mode='manual' spw='0'"]
        mydict = fh.parseDictionary(cmd)
        try:
            res = fh.evaluateFlagParameters(mydict, fparams)
        except Exception as instance:
            print("Expected error: %s" % instance)
        else:
            self.assertNotEqual(
                res, True, "Some parameter in cmd should not evaluate True"
            )

    def test_veryLongOnlineFlags(self):
        """flaghelper: evaluate a very long list of parameters"""
        # Get a big input file
        self.setUp_onlineFlags()
        fparams = dict(inspect.signature(flagdata).parameters.items())

        alist = fh.readFile(self.inpfile)
        adict = fh.parseDictionary(alist)
        os.system("rm -rf " + self.inpfile)
        self.assertTrue(fh.evaluateFlagParameters(adict, fparams))

    def test_parseXML1(self):
        """flaghelper: test parsing XML for online flags"""
        self.setUp_mockasdm()
        fdict = fh.parseXML(self.inpfile, 0.0)
        spw0 = fdict[0]["spw"]
        cmd0 = fdict[0]["command"]
        if "spw" in cmd0:
            is_spw = True
        self.assertTrue(spw0.__contains__("WVR#NOMINAL"))
        self.assertFalse(spw0.__contains__("WVR#Antenna"))
        self.assertEqual(len(fdict.keys()), 12572)
        self.assertTrue(is_spw)


if __name__ == "__main__":
    unittest.main()
