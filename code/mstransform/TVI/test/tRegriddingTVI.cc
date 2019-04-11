//# tRegriddingTVI: This file contains the unit tests of the RegriddingTVI class.
//#
//#  CASA - Common Astronomy Software Applications (http://casa.nrao.edu/)
//#  Copyright (C) Associated Universities, Inc. Washington DC, USA 2011, All rights reserved.
//#  Copyright (C) European Southern Observatory, 2011, All rights reserved.
//#
//#  This library is free software; you can redistribute it and/or
//#  modify it under the terms of the GNU Lesser General Public
//#  License as published by the Free software Foundation; either
//#  version 2.1 of the License, or (at your option) any later version.
//#
//#  This library is distributed in the hope that it will be useful,
//#  but WITHOUT ANY WARRANTY, without even the implied warranty of
//#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//#  Lesser General Public License for more details.
//#
//#  You should have received a copy of the GNU Lesser General Public
//#  License along with this library; if not, write to the Free Software
//#  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//#  MA 02111-1307  USA
//# $Id: $


#include <mstransform/TVI/test/tRegriddingTVI.h>

using namespace std;
using namespace casa;
using namespace casacore;
using namespace casa::vi;


//////////////////////////////////////////////////////////////////////////
// RegriddingTVITest class
//////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void RegriddingTVITest::generateTestFile()
{
    String path("");
    if (autoMode_p) path = String("/data/regression/unittest/mstransform/");
    copyTestFile(path,inpFile_p,testFile_p);

    return;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void RegriddingTVITest::generateReferenceFile()
{
    String path("");
    if (autoMode_p) path = String("/data/regression/unittest/mstransform/");
    copyTestFile(path,inpFile_p,referenceFile_p);

    return;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void RegriddingTVITest::initTestConfiguration(Record &configuration)
{
    testConfiguration_p = configuration;
    testConfiguration_p.define ("inputms", testFile_p);

    return;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void RegriddingTVITest::initReferenceConfiguration(Record &configuration)
{
    refConfiguration_p = configuration;
    refConfiguration_p.define ("inputms", referenceFile_p);
    refConfiguration_p.define ("regridms", true);

    return;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
RegriddingTVITest::RegriddingTVITest(): FreqAxisTVITest ()
{
    inpFile_p = String("CAS-7382.ms");
    testFile_p = String("CAS-7382.ms.test");
    referenceFile_p = String("CAS-7382.ms.ref");

    Record configuration;
    configuration.define ("spw", "1");
    configuration.define ("field", "Titan");
    configuration.define ("mode", "velocity");
    configuration.define ("width", "0.2km/s");
    configuration.define ("interpolation", "linear");
    configuration.define ("restfreq", "349.45370GHz");
    configuration.define ("outframe", "SOURCE");
    configuration.define ("datacolumn", String("ALL"));
    configuration.define ("reindex", false);

    init(configuration);
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
RegriddingTVITest::RegriddingTVITest(Record configuration): FreqAxisTVITest(configuration)
{
    init(configuration);
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void RegriddingTVITest::TestBody()
{
    SetUp();
    testCompareTransformedData();
    TearDown();

    return;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void RegriddingTVITest::testCompareTransformedData()
{
    // Declare working variables
    Float tolerance = 1E-5; // FLT_EPSILON is 1.19209290e-7F

    // Create MSTransformIterator pointing to reference file
    refConfiguration_p.define("factory",False);
    MSTransformIteratorFactory refFactory(refConfiguration_p);
    VisibilityIterator2 refTVI(refFactory);

    // Use MSTransformFactory to create a plain input VII
    testConfiguration_p.define("factory",True);
    MSTransformIteratorFactory plainVIFactory(testConfiguration_p);
    ViImplementation2 *inputVI = plainVIFactory.getInputVI()->getImpl();

    // Generate TVI to test
    RegriddingTVIFactory testFactory(testConfiguration_p,inputVI);
    VisibilityIterator2 testTVI(testFactory);

    // Determine columns to check
    VisBufferComponents2 columns;
    columns += VisBufferComponent2::NRows;
    columns += VisBufferComponent2::NChannels;
    columns += VisBufferComponent2::NCorrelations;
    columns += VisBufferComponent2::Frequencies;
    columns += VisBufferComponent2::FlagRow;
    columns += VisBufferComponent2::FlagCube;
    columns += VisBufferComponent2::VisibilityCubeObserved;
    columns += VisBufferComponent2::WeightSpectrum;
    columns += VisBufferComponent2::SigmaSpectrum;
    columns += VisBufferComponent2::Weight;
    columns += VisBufferComponent2::Sigma;

    // Compare
    SCOPED_TRACE("Comparing transformed data");
    compareVisibilityIterators(testTVI,refTVI,columns,tolerance);
}

//////////////////////////////////////////////////////////////////////////
// Googletest macros
//////////////////////////////////////////////////////////////////////////
TEST_F(RegriddingTVITest, testCompareTransformedData)
{
    testCompareTransformedData();
}

//////////////////////////////////////////////////////////////////////////
// main
//////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
    int ret;
    string parameter,value;
    Record configuration;
    Bool autoMode = true;

    for (unsigned short i=0;i<argc-1;i++)
    {
        parameter = string(argv[i]);
        value = string(argv[i+1]);

        if (parameter == string("-vis"))
        {
            configuration.define ("inputms", value);
            autoMode = false;
        }
        else if (parameter == string("-spw"))
        {
            configuration.define ("spw", value);
        }
    }

    if (autoMode)
    {
        ::testing::InitGoogleTest(&argc, argv);
        ret = RUN_ALL_TESTS();
    }
    else
    {
        RegriddingTVITest test(configuration);
        test.TestBody();
        if (test.getTestResult()) ret = 0;
        else ret = 1;
    }

    return ret;
}
