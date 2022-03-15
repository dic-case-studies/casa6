########################################################################
# test_task_slsearch.py
#
# Copyright (C) 2018
# Associated Universities, Inc. Washington DC, USA
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
# [Add the link to the JIRA ticket here once it exists]
#
# Based on the requirements listed in plone found here:
# https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.information.slsearch.html
#
#
##########################################################################
import os
import shutil
import unittest

from casatools import ctsys, spectralline, table
from casatasks import slsearch, casalog

good_table = "biglist.tbl"

_tb = table( )

def run_search(
    tab, outfile, freqrange, species, reconly,
    chemnames, qns, intensity, smu2, loga, el,
    eu, rrlinclude, rrlonly, verbose, logfile,
    append
):
    mysl = spectralline()
    restool = None
    if (not mysl.open(tab)):
        raise Exception
    try:
        return mysl.search(
            outfile=outfile, freqrange=freqrange,
            species=species, reconly=reconly,
            chemnames=chemnames, qns=qns,
            intensity=intensity, smu2=smu2, loga=loga,
            el=el, eu=eu, rrlinclude=rrlinclude,
            rrlonly=rrlonly, verbose=verbose, logfile=logfile,
            append=append
        )
    except:
        raise
    finally:
        mysl.done()

def run_slsearch(
    tab, outfile, freqrange, species, reconly,
    chemnames, qns, intensity, smu2, loga, el,
    eu, rrlinclude, rrlonly, verbose, logfile,
    append
):

    return slsearch(
        tablename=tab, outfile=outfile, freqrange=freqrange,
        species=species, reconly=reconly,
        chemnames=chemnames, qns=qns,
        intensity=intensity, smu2=smu2, loga=loga,
        el=el, eu=eu, rrlinclude=rrlinclude,
        rrlonly=rrlonly, verbose=verbose, logfile=logfile,
        append=append
    )


_mycount = 0

class slsearch_test(unittest.TestCase):
    
    def _testit(
        self, tab, outfile, freqrange, species, reconly,
        chemnames, qns, intensity, smu2, loga, el,
        eu, rrlinclude, rrlonly, verbose, logfile,
        append, nrows
    ):
        global _mycount
        mysl = spectralline()
        mytb = table()
        for i in [0, 1]:
            if (i==0):
                mysl = run_search(tab, outfile,
                    freqrange, species, reconly, chemnames,
                    qns, intensity, smu2, loga, el, eu,
                    rrlinclude, rrlonly, verbose, logfile, 
                    append
                )
            else:
                if (not outfile):
                    outfile = "count" + str(_mycount) + ".tbl"
                    _mycount = _mycount + 1
                self.assertTrue(
                    run_slsearch(
                        tab, outfile, freqrange, species,
                        reconly, chemnames, qns, intensity,
                        smu2, loga, el, eu, rrlinclude, rrlonly,
                        verbose, logfile, append
                    )
                )
                mysl.open(outfile)
            self.assertEqual(nrows, mysl.nrows())
            mysl.done()

            if (outfile):
                mytb.open(outfile)
                self.assertEqual(nrows, mytb.nrows())
                shutil.rmtree(outfile)
            mytb.done()

    
    def setUp(self):
        datapath=ctsys.resolve('unittest/slsearch/')
        shutil.copytree(os.path.join(datapath,good_table), good_table)

    def tearDown(self):
        shutil.rmtree(good_table)
        self.assertTrue(len(_tb.showcache()) == 0)

    def test_exceptions(self):
        """slsearch: Test various exception cases"""
        
        # check_search is used when run_search is expected to fail
        # includes closing the returned spectralline tool in
        # case run_search does not throw an exception
        def check_search(
                tab, outfile, freqrange, species, reconly,
                chemnames, qns, intensity, smu2, loga, el,
                eu, rrlinclude, rrlonly, verbose, logfile,
                append
        ):
            mysl = run_search(tab, outfile, freqrange, species, reconly,
                              chemnames, qns, intensity, smu2, loga, el,
                              eu, rrlinclude, rrlonly, verbose, logfile,
                              append)
            mysl.done()

            
        def testit(
            tab, outfile, freqrange, species, reconly,
            chemnames, qns, intensity, smu2, loga, el,
            eu, rrlinclude, rrlonly, verbose, logfile, 
            append
        ):
            for i in [0,1]:
                if (i==0):
                    self.assertRaises(
                        Exception, check_search, tab, outfile,
                        freqrange, species, reconly, chemnames,
                        qns, intensity, smu2, loga, el, eu,
                        rrlinclude, rrlonly, verbose, logfile, 
                        append
                    )
                    self.assertTrue(len(_tb.showcache()) == 0)
                else:
                    # CASA6 slsearch raises an exception, CASA5 returns None
                    self.assertRaises(
                            Exception, run_slsearch,
                            tab, outfile, freqrange, species,
                            reconly, chemnames, qns, intensity,
                            smu2, loga, el, eu, rrlinclude, rrlonly,
                            verbose, logfile, append
                    )

                    # either way, no tables should be open
                    self.assertTrue(len(_tb.showcache()) == 0)

        # bogus input table name
        # the version of testit used here throws an exception if the
        # expected exceptions or return values did not happen
        try:
            testit(
                tab="fred.tbl", outfile="x", freqrange=[0, 100], species=[],
                reconly=True, chemnames=[], qns=[], intensity=[-1], smu2=[-1],
                loga=[-1], el=[-1], eu=[-1], rrlinclude=True, rrlonly=True,
                verbose=True, logfile="", append=True
            )
        except:
            casalog.post("Failure in test_exceptions testing bogus input table name", 'SEVERE')
            raise

        # bad output name
        try:
            testit(
                tab=good_table, outfile="foo/bar/bad", freqrange=[0, 100], species=[],
                reconly=True, chemnames=[], qns=[], intensity=[-1], smu2=[-1],
                loga=[-1], el=[-1], eu=[-1], rrlinclude=True, rrlonly=True,
                verbose=True, logfile="", append=True
            )
        except:
            casalog.post("Failure in test_exceptions testing bad output name", 'SEVERE')
            raise
            
    def test_table(self):
        """ test various settings of the table parameter"""

        # no table name works because it defaults to the system spectral line table
        self._testit(
            tab="", outfile="", freqrange=[0, 100], species=[],
            reconly=True, chemnames=[], qns=[], intensity=[-1], smu2=[-1],
            loga=[-1], el=[-1], eu=[-1], rrlinclude=True, rrlonly=False,
            verbose=False, logfile="", append=True, nrows=59998
        )
        # test user specified table search
        self._testit(
            tab=good_table, outfile="", freqrange=[0, 100], species=[],
            reconly=True, chemnames=[], qns=[], intensity=[-1], smu2=[-1],
            loga=[-1], el=[-1], eu=[-1], rrlinclude=True, rrlonly=False,
            verbose=False, logfile="", append=True, nrows=67858
        )

    def test_outfile(self):
        """ test various settings of the outfile parameter"""

        outfile = "blah.tbl"
        self._testit(
            tab=good_table, outfile=outfile, freqrange=[0, 100], species=[],
            reconly=True, chemnames=[], qns=[], intensity=[-1], smu2=[-1],
            loga=[-1], el=[-1], eu=[-1], rrlinclude=True, rrlonly=False,
            verbose=False, logfile="", append=True, nrows=67858
        )
        
    def test_freqrange(self):
        """ test various settings of the freqrange parameter"""

        self._testit(
            tab=good_table, outfile="", freqrange=[0, 100], species=[],
            reconly=True, chemnames=[], qns=[], intensity=[-1], smu2=[-1],
            loga=[-1], el=[-1], eu=[-1], rrlinclude=True, rrlonly=False,
            verbose=False, logfile="", append=True, nrows=67858
        )
        self._testit(
            tab=good_table, outfile="", freqrange=[0, 10], species=[],
            reconly=True, chemnames=[], qns=[], intensity=[-1], smu2=[-1],
            loga=[-1], el=[-1], eu=[-1], rrlinclude=True, rrlonly=False,
            verbose=False, logfile="", append=True, nrows=15292
        )

    def test_species(self):
        """ test various settings of the species parameter"""

        self._testit(
            tab=good_table, outfile="", freqrange=[0, 100], species=['S18O'],
            reconly=True, chemnames=[], qns=[], intensity=[-1], smu2=[-1],
            loga=[-1], el=[-1], eu=[-1], rrlinclude=True, rrlonly=False,
            verbose=False, logfile="", append=True, nrows=9
        )
        self._testit(
            tab=good_table, outfile="", freqrange=[0, 100], species=["S18O","HC5Nv11=1"],
            reconly=True, chemnames=[], qns=[], intensity=[-1], smu2=[-1],
            loga=[-1], el=[-1], eu=[-1], rrlinclude=True, rrlonly=False,
            verbose=False, logfile="", append=True, nrows=81
        )
        
    def test_chemnames(self):
        """ test various settings of the chemnames parameter"""

        self._testit(
            tab=good_table, outfile="", freqrange=[0, 100], species=[],
            reconly=True, chemnames=["Silicon Monocarbide"], qns=[],
            intensity=[-1], smu2=[-1],
            loga=[-1], el=[-1], eu=[-1], rrlinclude=True, rrlonly=False,
            verbose=False, logfile="", append=True, nrows=6
        )
        self._testit(
            tab=good_table, outfile="", freqrange=[0, 100], species=[],
            reconly=True, chemnames=["Silicon Monocarbide", "Potassium chloride"],
            qns=[], intensity=[-1], smu2=[-1],
            loga=[-1], el=[-1], eu=[-1], rrlinclude=True, rrlonly=False,
            verbose=False, logfile="", append=True, nrows=32
        )

    def test_qns(self):
        """ test various settings of the qns parameter"""

        self._testit(
            tab=good_table, outfile="", freqrange=[0, 100], species=[],
            reconly=True, chemnames=[], qns=["11-10"],
            intensity=[-1], smu2=[-1],
            loga=[-1], el=[-1], eu=[-1], rrlinclude=True, rrlonly=False,
            verbose=False, logfile="", append=True, nrows=17
        )

    def test_intensity(self):
        """ test various settings of the intensity parameter"""

        self._testit(
            tab=good_table, outfile="", freqrange=[0, 100], species=[],
            reconly=True, chemnames=[], qns=[],
            intensity=[-10,-8], smu2=[-1],
            loga=[-1], el=[-1], eu=[-1], rrlinclude=True, rrlonly=False,
            verbose=False, logfile="", append=True, nrows=13447
        )
        self._testit(
            tab=good_table, outfile="", freqrange=[0, 100], species=[],
            reconly=True, chemnames=[], qns=[],
            intensity=[-10,-8], smu2=[-1],
            loga=[-1], el=[-1], eu=[-1], rrlinclude=False, rrlonly=False,
            verbose=False, logfile="", append=True, nrows=7626
        )

    def test_smu2(self):
        """ test various settings of the smu2 parameter"""

        self._testit(
            tab=good_table, outfile="", freqrange=[0, 100], species=[],
            reconly=True, chemnames=[], qns=[],
            intensity=[], smu2=[5, 10],
            loga=[-1], el=[-1], eu=[-1], rrlinclude=True, rrlonly=False,
            verbose=False, logfile="", append=True, nrows=12227
        )
        self._testit(
            tab=good_table, outfile="", freqrange=[0, 100], species=[],
            reconly=True, chemnames=[], qns=[],
            intensity=[-1], smu2=[5, 10],
            loga=[-1], el=[-1], eu=[-1], rrlinclude=False, rrlonly=False,
            verbose=False, logfile="", append=True, nrows=6406
        )

    def test_loga(self):
        """ test various settings of the loga parameter"""

        self._testit(
            tab=good_table, outfile="", freqrange=[0, 100], species=[],
            reconly=True, chemnames=[], qns=[],
            intensity=[], smu2=[],
            loga=[-6, -4], el=[-1], eu=[-1], rrlinclude=True, rrlonly=False,
            verbose=False, logfile="", append=True, nrows=22781
        )
        self._testit(
            tab=good_table, outfile="", freqrange=[0, 100], species=[],
            reconly=True, chemnames=[], qns=[],
            intensity=[], smu2=[],
            loga=[-6, -4], el=[-1], eu=[-1], rrlinclude=False, rrlonly=False,
            verbose=False, logfile="", append=True, nrows=16960
        )

    def test_eu(self):
        """ test various settings of the eu parameter"""

        self._testit(
            tab=good_table, outfile="", freqrange=[0, 100], species=[],
            reconly=True, chemnames=[], qns=[],
            intensity=[], smu2=[],
            loga=[], el=[-1], eu=[150,200], rrlinclude=True, rrlonly=False,
            verbose=False, logfile="", append=True, nrows=9079
        )
        self._testit(
            tab=good_table, outfile="", freqrange=[0, 100], species=[],
            reconly=True, chemnames=[], qns=[],
            intensity=[], smu2=[],
            loga=[], el=[-1], eu=[150,200], rrlinclude=False, rrlonly=False,
            verbose=False, logfile="", append=True, nrows=3258
        )
        self._testit(
            tab=good_table, outfile="", freqrange=[1, 1.1], species=[],
            reconly=True, chemnames=[], qns=[],
            intensity=[], smu2=[],
            loga=[], el=[], eu=[1581.52, 1581.53], rrlinclude=False, rrlonly=False,
            verbose=False, logfile="", append=True, nrows=0
        )
        self._testit(
            tab=good_table, outfile="", freqrange=[1, 1.1], species=[],
            reconly=True, chemnames=[], qns=[],
            intensity=[], smu2=[],
            loga=[], el=[], eu=[1581.57, 1581.58] , rrlinclude=False, rrlonly=False,
            verbose=False, logfile="", append=True, nrows=1
        )


    def test_el(self):
        """ test various settings of the el parameter"""

        self._testit(
            tab=good_table, outfile="", freqrange=[0, 100], species=[],
            reconly=True, chemnames=[], qns=[],
            intensity=[], smu2=[],
            loga=[], el=[150,200], eu=[], rrlinclude=True, rrlonly=False,
            verbose=False, logfile="", append=True, nrows=9023
        )
        self._testit(
            tab=good_table, outfile="", freqrange=[0, 100], species=[],
            reconly=True, chemnames=[], qns=[],
            intensity=[], smu2=[],
            loga=[], el=[150,200], eu=[], rrlinclude=False, rrlonly=False,
            verbose=False, logfile="", append=True, nrows=3202
        )
        self._testit(
            tab=good_table, outfile="", freqrange=[1, 1.1], species=[],
            reconly=True, chemnames=[], qns=[],
            intensity=[], smu2=[],
            loga=[], el=[1581.52, 1581.53], eu=[], rrlinclude=False, rrlonly=False,
            verbose=False, logfile="", append=True, nrows=1
        )
        self._testit(
            tab=good_table, outfile="", freqrange=[1, 1.1], species=[],
            reconly=True, chemnames=[], qns=[],
            intensity=[], smu2=[],
            loga=[], el=[1581.57, 1581.58] , eu=[], rrlinclude=False, rrlonly=False,
            verbose=False, logfile="", append=True, nrows=0
        )

    def test_rrlonly(self):
        """ test various settings of the rrlonly parameter"""

        self._testit(
            tab=good_table, outfile="", freqrange=[0, 100], species=[],
            reconly=True, chemnames=[], qns=[],
            intensity=[], smu2=[],
            loga=[], el=[], eu=[], rrlinclude=True, rrlonly=False,
            verbose=False, logfile="", append=True, nrows=67858
        )
        self._testit(
            tab=good_table, outfile="", freqrange=[0, 100], species=[],
            reconly=True, chemnames=[], qns=[],
            intensity=[], smu2=[],
            loga=[], el=[], eu=[], rrlinclude=True, rrlonly=True,
            verbose=False, logfile="", append=True, nrows=5821
        )

    def test_logfile(self):
        """ test various settings of the logfile and append parameters"""

        def count_lines(txtfile):
            count = 0
            with open(txtfile,'r')  as f:
                for count,l in enumerate(f,1):
                    pass
            return count
        
        logfile = "xx.log"
        # verbose = False so no logfile should be written
        self._testit(
            tab=good_table, outfile="", freqrange=[0, 100], species=[],
            reconly=True, chemnames=[], qns=[],
            intensity=[], smu2=[],
            loga=[], el=[], eu=[], rrlinclude=True, rrlonly=True,
            verbose=False, logfile=logfile, append=False, nrows=5821
        )
        self.assertFalse(os.path.exists(logfile))
        # verbose and overwrite
        self._testit(
            tab=good_table, outfile="", freqrange=[0, 100], species=[],
            reconly=True, chemnames=[], qns=[],
            intensity=[], smu2=[],
            loga=[], el=[], eu=[], rrlinclude=True, rrlonly=True,
            verbose=True, logfile=logfile, append=False, nrows=5821
        )
        self.assertTrue(os.path.exists(logfile))

        num_lines = count_lines(logfile)
        self.assertEqual(num_lines, 5822)
        # append (twice)
        self._testit(
            tab=good_table, outfile="", freqrange=[0, 100], species=[],
            reconly=True, chemnames=[], qns=[],
            intensity=[], smu2=[],
            loga=[], el=[], eu=[], rrlinclude=True, rrlonly=True,
            verbose=True, logfile=logfile, append=True, nrows=5821
        )
        self.assertTrue(os.path.exists(logfile))
        num_lines = count_lines(logfile)
        self.assertEqual(num_lines, 3*5822)
        os.remove(logfile)

if __name__ == '__main__':
    unittest.main()
