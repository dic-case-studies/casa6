//# TestUtilsTVI.h This file contains the interface definition of the TestUtilsTVI class.
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

#ifndef TestUtilsTVI_H_
#define TestUtilsTVI_H_

// Google test
#include <gtest/gtest.h>

// casacore containers
#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Containers/Record.h>

// VI/VB framework
#include <msvis/MSVis/VisBuffer2.h>
#include <msvis/MSVis/VisibilityIterator2.h>
#include <msvis/MSVis/test/MsFactory.h>


namespace casa { //# NAMESPACE CASA - BEGIN

namespace vi { //# NAMESPACE VI - BEGIN

//////////////////////////////////////////////////////////////////////////
// FreqAxisTVITest class
//////////////////////////////////////////////////////////////////////////
class FreqAxisTVITest: public ::testing::Test {

public:

    FreqAxisTVITest();
    FreqAxisTVITest(casacore::Record configuration);
    virtual ~FreqAxisTVITest();

    void SetUp();
    void TearDown();
    casacore::Bool getTestResult() {return testResult_p;}

protected:

    void init(casacore::Record &configuration);
    virtual void generateTestFile() = 0;
    virtual void generateReferenceFile() = 0;
    virtual void initTestConfiguration(casacore::Record &configuration) = 0;
    virtual void initReferenceConfiguration(casacore::Record &configuration) = 0;

    casacore::Bool autoMode_p;
    casacore::Bool testResult_p;
    casacore::String inpFile_p;
    casacore::String testFile_p;
    casacore::String referenceFile_p;
    casacore::Record refConfiguration_p;
    casacore::Record testConfiguration_p;
};

class MsFactoryTVITester : public ::testing::Test
{
public:
  
    /*
     * Constructor: create the temporary dir and the MsFactory used later on
     * to create the MS.
     * @param testSubdir The subdirectory to create under the temporary dir
     * @param msName     The name of the MS to create (under testSubdir)
     */
    MsFactoryTVITester(const std::string& testSubdir ="MsFactoryTVITester",
                       const std::string& msName = "MsFactory");

    //Google test fixture constructor
    void SetUp();

    //Google test fixture destructor
    void TearDown();

    /*
     * Create the synthetic MS  
     */
    void createMS();

    /*
     * Create the TVI stack to access the MS. 
     */
    void instantiateVI(std::vector<ViiLayerFactory*>& factories);

    /*
     * Iterate the whole MS calling a user provided function.
     * The only useful case is having a lambda as a visitor function.
     * It can access the visibility buffer from variable vb_p
     */
    void visitIterator(std::function<void(void)> visitor);

    /* 
     * Return the MsFactory for full customization if needed 
     */
    casa::vi::test::MsFactory& getMsFactory();
  
    //Destructor
    ~MsFactoryTVITester();

protected:
    
    //The temporary dir where the synthetic MS is created  
    char tmpdir_p[_POSIX_PATH_MAX];

    //The helper class to create synthetic MS
    std::unique_ptr<casa::vi::test::MsFactory> msf_p;

    //The synthetic MS.
    std::unique_ptr<casacore::MeasurementSet> ms_p;

    //The VisibilityIterator2 used to iterate trough the data
    std::unique_ptr<VisibilityIterator2> vi_p;

    //The name of the ms created by MsFactory
    std::string msName_p;

    //The attached VisBuffer
    VisBuffer2 * vb_p;
};


//////////////////////////////////////////////////////////////////////////
// Convenience methods
//////////////////////////////////////////////////////////////////////////
template <class T> void compareVector(const casacore::Char* column,
                                      const casacore::Vector<T> &inp,
                                      const casacore::Vector<T> &ref,
                                      casacore::Float tolerance = FLT_EPSILON);

template <class T> void compareMatrix(const casacore::Char* column,
                                      const casacore::Matrix<T> &inp,
                                      const casacore::Matrix<T> &ref,
                                      casacore::Float tolerance = FLT_EPSILON);

template <class T> void compareCube(const casacore::Char* column,
                                    const casacore::Cube<T> &inp,
                                    const casacore::Cube<T> &ref,
                                    casacore::Float tolerance = FLT_EPSILON);

template <class T> void compareCubesVector(const casacore::Char* column,
                                           const casacore::Vector<casacore::Cube<T>> &inp,
                                           const casacore::Vector<casacore::Cube<T>> &ref,
                                           casacore::Float tolerance = FLT_EPSILON);

/*
 * Compare two VisBuffers
 * The components specified in columns will be compared using
 * a float point tolerance of tolerance.
 * If different data columns should be compared together (for instance
 * if DATA in refVb should be compared to CORRECTED in testVb) then
 * the datacolmap map can be used to define the map. The map keyword
 * is the column in testVb and the map value is the column in refVb.
 */
void compareVisBuffers(VisBuffer2 &testVb,
                       VisBuffer2 &refVb,
                       VisBufferComponents2 &columns,
                       casacore::Float tolerance = FLT_EPSILON,
                       std::map<casacore::MS::PredefinedColumns,casacore::MS::PredefinedColumns> *datacolmap = NULL);

/*
 * Compare two Vis iterators
 * This will run the full iteration and for each iteration a comparison
 * of the VisBuffers is done.
 * The VisBuffer components specified in columns will be compared using
 * a float point tolerance of tolerance.
 * If different data columns should be compared together (for instance
 * if DATA in refVb should be compared to CORRECTED in testVb) then
 * the datacolmap map can be used to define the map. The map keyword
 * is the column in testVb and the map value is the column in refVb.
 */
void compareVisibilityIterators(VisibilityIterator2 &testTVI,
                                VisibilityIterator2 &refTVI,
                                VisBufferComponents2 &columns,
                                casacore::Float tolerance = FLT_EPSILON,
                                std::map<casacore::MS::PredefinedColumns,casacore::MS::PredefinedColumns> *datacolmap = NULL);

/*
 * Compare two Vis iterators
 * This will run the full iteration and for each iteration the visitor
 * function will be run.
 * The visitor function can access variables testVb and refVb of type
 * VisBuffer2 *, which will give full access to the contents of both
 * VisBuffers for each iteration
 */
void compareVisibilityIterators(VisibilityIterator2 &testTVI,
                                VisibilityIterator2 &refTVI,
                                std::function<void(VisBuffer2* testVb, VisBuffer2* refVb)> visitor);

void copyTestFile(casacore::String &path,casacore::String &filename,casacore::String &outfilename);

const casacore::Cube<casacore::Complex> & getViscube(VisBuffer2 &vb,
                                                     casacore::MS::PredefinedColumns datacol,
                                                     std::map<casacore::MS::PredefinedColumns,casacore::MS::PredefinedColumns> *datacolmap);

const casacore::Vector<casacore::Cube<casacore::Complex>> & getViscubes(VisBuffer2 &vb,
                                                                        casacore::MS::PredefinedColumns datacol,
                                                                        std::map<casacore::MS::PredefinedColumns,casacore::MS::PredefinedColumns> *datacolmap);

void flagEachOtherChannel(VisibilityIterator2 &vi, bool undoChanbin, int chanbin = 1);

} //# NAMESPACE VI - END

} //# NAMESPACE CASA - END

#endif /* TestUtilsTVI_H_ */
