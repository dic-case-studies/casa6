//# tMSTransformIterator: This file contains the unit tests of the MSTransformIterator class.
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

#ifndef MSTransformIteratorTest_H_
#define MSTransformIteratorTest_H_

#include "TestUtilsTVI.h"
#include <mstransform/MSTransform/MSTransformIteratorFactory.h>

using namespace std;
using namespace casa;
using namespace casa::vi;


/////////////////////////////////////////////////////////////////////////
// MSTransformIteratorTest class
/////////////////////////////////////////////////////////////////////////
class MSTransformIteratorTest: public FreqAxisTVITest {

public:

	MSTransformIteratorTest();
	MSTransformIteratorTest(casacore::Record configuration);
	~MSTransformIteratorTest();

    void TestBody();
    void testCompareTransformedData();

protected:

    void generateTestFile();
    void generateReferenceFile();
    void initTestConfiguration(casacore::Record &configuration);
    void initReferenceConfiguration(casacore::Record &configuration);
};

#endif /* MSTransformIteratorTest_H_ */
