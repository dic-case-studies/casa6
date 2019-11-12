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


#include <mstransform/TVI/test/tStatWtVarianceAndWeightCalculator.h>

#include <casacore/scimath/StatsFramework/ClassicalStatistics.h>

using namespace std;
using namespace casa;
using namespace casacore;

StatWtVarianceAndWeightCalculatorTest::StatWtVarianceAndWeightCalculatorTest(): ::testing::Test () {}

StatWtVarianceAndWeightCalculatorTest::~StatWtVarianceAndWeightCalculatorTest() {}

void StatWtVarianceAndWeightCalculatorTest::TestBody() {
    SetUp();
    TearDown();
}

//////////////////////////////////////////////////////////////////////////
// Googletest macros
//////////////////////////////////////////////////////////////////////////
/*
TEST_F(StatWtTVITest, testCompareTransformedData)
{
	testCompareTransformedData();
}
*/

TEST (StatWtVarianceAndWeightCalculatorTest, testComputeVariance) {
    shared_ptr<
        ClassicalStatistics<
            Double, Array<Float>::const_iterator,
            Array<Bool>::const_iterator, Array<Double>::const_iterator
        >
    > statAlg(
        new ClassicalStatistics<
            Double, Array<Float>::const_iterator,
            Array<Bool>::const_iterator, Array<Double>::const_iterator
        >()
    );
    shared_ptr<map<uInt, pair<uInt, uInt>>> samples(
        new map<uInt, pair<uInt, uInt>>()
    );
    uInt spw = 5;
    (*samples)[spw] = make_pair(0, 0);
    StatWtVarianceAndWeightCalculator calc(statAlg, samples);
    uInt npol = 2;
    uInt nchan = 3;
    uInt nrow = 4;
    Cube<Complex> data(npol, nchan, nrow);
    Cube<Bool> flags(npol, nchan, nrow, False);
    Vector<Double> exposures(nrow);
    for (uInt i=0; i<npol; ++i) {
        for (uInt j=0; j<nchan; ++j) {
            for (uInt k=0; k<nrow; ++k) {
                data(i, j, k) = Complex((i + 0.5)*j*k, i*j*(k+0.5));
            }
        }
    }
    flags(0, 1, 1) = True;
    exposures[0] = 4.5;
    exposures[1] = 6;
    exposures[2] = 5;
    exposures[3] = 5.5;
    auto var = calc.computeVariance(data, flags, exposures, spw);
    cout << "var " << var << endl;
    cout << "samples " << *samples << endl;
    ASSERT_FLOAT_EQ(23.100544223707672, var);
}

//////////////////////////////////////////////////////////////////////////
// main
//////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
    /*
	int ret;
	string parameter, value;
	Record configuration;
	Bool autoMode = true;

	for (unsigned short i=0; i<argc-1; ++i) {
		parameter = string(argv[i]);
		value = string(argv[i+1]);
		if (parameter == string("-vis")) {
			configuration.define ("inputms", value);
			autoMode = false;
		}
		else if (parameter == string("-spw")) {
			configuration.define ("spw", value);
		}
	}
	if (autoMode) {
		::testing::InitGoogleTest(&argc, argv);
		ret = RUN_ALL_TESTS();
	}
	else {
        StatWtTVITest test(configuration);
        test.TestBody();
        ret = test.getTestResult() ? 0: 1;
	}
	return ret;
	*/
}
