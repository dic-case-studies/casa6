//# tViiLayerFactory.cc: Tests Recursive factory for layered VI2s
//# Copyright (C) 1995,1999,2000,2001,2016
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This program is free software; you can redistribute it and/or modify it
//# under the terms of the GNU General Public License as published by the Free
//# Software Foundation; either version 2 of the License, or (at your option)
//# any later version.
//#
//# This program is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
//# more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with this program; if not, write to the Free Software Foundation, Inc.,
//# 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id$


#include <casa/aips.h>
#include <casa/Exceptions/Error.h>
#include <casacore/casa/OS/EnvVar.h>
#include <casacore/casa/OS/Path.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casa/iostream.h>
#include <msvis/MSVis/IteratingParameters.h>
#include <msvis/MSVis/ViiLayerFactory.h>
#include <msvis/MSVis/LayeredVi2Factory.h>
#include <msvis/MSVis/TransformingVi2.h>
#include <msvis/MSVis/SimpleSimVi2.h>
#include <msvis/MSVis/VisBuffer2.h>
#include <msvis/MSVis/test/TestUtilsTVI.h>
#include <casa/iomanip.h>
#include <gtest/gtest.h>

using namespace std;
using namespace casa;
using namespace casacore;
using namespace casa::vi;
using namespace casa::vi::test;


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

TEST( ViiLayerFactoryTest , ViiLayerFactoryBasicTest ) {
 
  // A very rudimentary test of a single layer

  SimpleSimVi2Parameters s0;
  SimpleSimVi2LayerFactory fac(s0);

  Vector<ViiLayerFactory*> facts(1);
  facts[0]=&fac;

  std::unique_ptr<VisibilityIterator2> vi(new VisibilityIterator2(facts));
  VisBuffer2 *vb = vi->getImpl()->getVisBuffer();

  // Has a viable VI2 been generated?
  Int niter(0);
  for (vi->originChunks();vi->moreChunks();vi->nextChunk()) {
    for (vi->origin();vi->more();vi->next()) {
      ASSERT_EQ(4,vb->nAntennas());
      ++niter;
    }
  }
  ASSERT_EQ(1,niter);
}


TEST( ViiLayerFactoryTest , ViiLayerFactoryRealDataBasicTest ) {
 
  // A very rudimentary test of a single layer, using a real MS


  String *casapath = new String[2];
  split(EnvironmentVariable::get("CASAPATH"), casapath, 2, String(" "));
  // Use of Path().absoluteName() absorbs relative stuff in casapath
  String mspath(Path(casapath[0]+"/data/regression/unittest/flagdata/Four_ants_3C286.ms").absoluteName());

  MeasurementSet ms(mspath);

  Double interval(60000.0);
  IteratingParameters ipar(interval);   // seems to include SCAN_NUMBER automatically?
  VisIterImpl2LayerFactory fac(&ms,ipar,False); 

  Vector<ViiLayerFactory*> facts(1);
  facts[0]=&fac;

  std::unique_ptr<VisibilityIterator2> vi(new VisibilityIterator2(facts));
  VisBuffer2 *vb = vi->getImpl()->getVisBuffer();

  // Has a viable VI2 been generated?
  Int chunk(0),niter(0);
  for (vi->originChunks();vi->moreChunks();vi->nextChunk(),++chunk) {
    vi->origin();
    /*
    cout << "ch="<< chunk
	 << " scan="<<vb->scan()(0)
	 << " field="<< vb->fieldId()(0)
	 << " spw="<< vb->spectralWindows()(0)
	 << endl;
    */

    for (vi->origin();vi->more();vi->next(),++niter) {

      /*
      cout << "*************************************" << endl;
      cout << "chunk="<< chunk << ";  niter=" << niter << endl;

      cout << " scan=" << vb->scan()(0) << endl;
      cout << " fieldId=" << vb->fieldId()(0) << endl;

      cout << "nAntennas=" << vb->nAntennas() << endl;
      cout << "nRows=" << vb->nRows() << endl;
      cout << "nChannels=" << vb->nChannels() << endl;
      cout << "nCorrelations=" << vb->nCorrelations() << endl;
      cout << "spectralWindows=" << vb->spectralWindows() << endl;
      */

      ASSERT_EQ(4,vb->nAntennas());
      ASSERT_EQ(6U,vb->nRows());
      ASSERT_EQ(64,vb->nChannels());  // all spws
      ASSERT_EQ(4,vb->nCorrelations());

    }
  }

  //cout << "chunk=" << chunk << endl;
  //  cout << "niter=" << niter << endl;

  ASSERT_EQ(32,chunk);
  ASSERT_EQ(2864,niter);
  delete[] casapath;
}

/* 
 * This class is a very simple TVI that passes all the requests for the
 * underlying VI (delegating to TransformingVi2) except that it swaps
 * the access to corrected data and data, i.e., when accessing the 
 * data corrected it returns the data column of the underlying VI.
 */
class DataSwappingTVI : public TransformingVi2
{
public:
  
  //Constructor
  DataSwappingTVI(ViImplementation2 * inputVii) :
    TransformingVi2 (inputVii)
  {
    setVisBuffer(createAttachedVisBuffer (VbRekeyable));
  }

  //Access to DATA CORRECTED returns underlying DATA column
  virtual void visibilityCorrected(casacore::Cube<casacore::Complex>& vis) const
  {
    VisBuffer2* vb = getVii()->getVisBuffer();
    vis = vb->visCube();
  }
  
  //Access to DATA returns underlying DATA CORRECTED column
  virtual void visibilityObserved(casacore::Cube<casacore::Complex>& vis) const
  {
    VisBuffer2* vb = getVii()->getVisBuffer();
    vis = vb->visCubeCorrected();
  }
  
  void origin()
  {
    // Drive underlying ViImplementation2
    getVii()->origin();

    // Synchronize own VisBuffer
    configureNewSubchunk();
  }

  void next()
  {
    // Drive underlying ViImplementation2
    getVii()->next();

    // Synchronize own VisBuffer
    configureNewSubchunk();
  }

};

/*
 * Factory that allows the creation of DataSwappingTVI classes.
 * This factory doesn't have any parameter to configure
 */
class DataSwappingTVILayerFactory : public ViiLayerFactory
{

public:

  DataSwappingTVILayerFactory()
  {
  }

  virtual ~DataSwappingTVILayerFactory() {};

protected:

  virtual ViImplementation2 * createInstance(ViImplementation2* vii0) const
  {
    ViImplementation2 *vii = new DataSwappingTVI(vii0);
    return vii;
  }
};

/*
 * Gtest fixture used to test the access to the different data columns.
 * This class will create a synthetic MS in a temporary directory with or 
 * without DATA CORRECTED columns. 
 * It will also create a stack of TVIs with two TVIs: 
 * a disk access layer and a swapping data TVI (DataSwappingTVI). This TVI
 * is optional.
 */
class DataAccessTest : public MsFactoryTVITester
{
public:
  
  /*
   * Constructor: create the temporary dir and the MsFactory used later on
   * to create the MS.
   */
  DataAccessTest() :
    MsFactoryTVITester("test_tViiLayerFactory","DataAccessTest")
  {
  }

  /*
   * Do not create DATA CORRECTED column. 
   * This function must be called before createTVIs()
   */
  void removeCorrectedData()
  {
    msf_p->removeColumn(MS::CORRECTED_DATA);
  }

  /*
   * Do not create DATA column. 
   * This function must be called before createTVIs()
   */
  void removeData()
  {
    msf_p->removeColumn(MS::DATA);
  }

  /*
   * Switch the creation of the DataSwappingTVI on top the disk access layer.
   * This function must be called before createTVIs(). 
   */
  void addSwappingDataTVI()
  {
    withSwappingDataTVI_p = true;
  }
  
  /*
   * Create the synthetic MS and the TVI stack to access it. 
   */
  void createTVIs()
  {
    createMS();
    
    //Create a disk layer type VI Factory
    IteratingParameters ipar;
    VisIterImpl2LayerFactory diskItFac(ms_p.get(),ipar,false);

    //Create a SwappingDataTVI Factory if requested
    std::unique_ptr<DataSwappingTVILayerFactory> swapFac;
    if(withSwappingDataTVI_p)
      swapFac.reset(new DataSwappingTVILayerFactory());

    //Create a layered factory with all the layers of factories
    size_t nFac = 1;
    if(withSwappingDataTVI_p) 
      nFac++;
    std::vector<ViiLayerFactory*> factories(nFac);
    factories[0]=&diskItFac;
    if(withSwappingDataTVI_p)
      factories[1]= swapFac.get();

    instantiateVI(factories);
  }

  //Destructor
  ~DataAccessTest()
  {
  }

  //Wether to use the DataSwappingTVI
  bool withSwappingDataTVI_p = false;
};
 

/*
 * This test will simply access the corrected data column of a
 * synthetic created MS
 */
TEST_F(DataAccessTest, AccessCorrectedData)
{
  createTVIs();

  //Traverse the iterator accessing the corrected data cube
  visitIterator([&]() -> void {vb_p->visCubeCorrected().shape();});
}

/*
 * This test will check that an exception is thrown if the
 * CORRECTED DATA column is missing in the MS.
 */
TEST_F(DataAccessTest, AccessCorrectedDataWhenMissing)
{
  removeCorrectedData();

  createTVIs();

  //Traverse the iterator accessing the corrected data cube. This should
  //throw, since it has been removed from the MS.
  ASSERT_THROW(visitIterator([&]() -> void {vb_p->visCubeCorrected().shape();}),
               AipsError);
}

/*
 * This test will access the corrected data column of a
 * synthetic MS that doesn't contain that column. The upper TVI, however
 * will swap the access to DATA column instead so the test should succeed
 */
TEST_F(DataAccessTest, AccessCorrectedDataInSwappingDataTVI)
{
  removeCorrectedData();

  addSwappingDataTVI();
  
  createTVIs();

  //Traverse the iterator accessing the corrected data cube. This works
  //despite removing the corrected data from disk because the upper TVI layer
  //(DataSwappingTVI) delivers DATA when CORRECTED DATA is requested.
  visitIterator([&]() -> void {vb_p->visCubeCorrected().shape();});
}

/*
 * This test will access the corrected data column of a
 * synthetic that doesn't contain column DATA. The upper TVI, however,
 * will deliver the underlying DATA column when asking for CORRECTED DATA,
 * and therefore the test will fail.
 */
TEST_F(DataAccessTest, AccessCorrectedDataInSwappingDataTVIWhenMissingData)
{
  removeData();

  addSwappingDataTVI();
  
  createTVIs();

  //Traverse the iterator accessing the corrected data cube. 
  //The swapping TVI will access the underlying DATA column when accessing
  //the CORRECTED DATA. Since DATA has been removed from the MS this should
  //throw
  ASSERT_THROW(visitIterator([&]() -> void {vb_p->visCubeCorrected().shape();}),
               AipsError);
}

/* 
 * This test will simply access the data column of a
 * synthetic created MS
 */
TEST_F(DataAccessTest, AccessData)
{
  createTVIs();

  //Traverse the iterator accessing the corrected data cube
  visitIterator([&]() -> void {vb_p->visCube().shape();});
}

/*
 * This test will check that an exception is thrown if the
 * CORRECTED DATA column is missing in the MS.
 */
TEST_F(DataAccessTest, AccessDataWhenMissing)
{
  removeData();

  createTVIs();

  //Traverse the iterator accessing the corrected data cube. This should
  //throw, since it has been removed from the MS.
  ASSERT_THROW(visitIterator([&]() -> void {vb_p->visCube().shape();}),
               AipsError);
}

/*
 * This test will access the corrected data column of a
 * synthetic that doesn't contain that column. The upper TVI, however
 * will swap the access to DATA column instead so the test should succeed
 */
TEST_F(DataAccessTest, AccessDataInSwappingDataTVI)
{
  removeData();

  addSwappingDataTVI();

  createTVIs();

  //Traverse the iterator accessing the corrected data cube. This works
  //despite removing the corrected data from disk because the upper TVI layer
  //(DataSwappingTVI) delivers DATA when CORRECTED DATA is requested.
  visitIterator([&]() -> void {vb_p->visCube().shape();});
}

TEST_F(DataAccessTest, AccessDataInSwappingDataTVIWhenMissingCorrectedData)
{
  removeCorrectedData();

  addSwappingDataTVI();

  createTVIs();

  //Traverse the iterator accessing the corrected data cube.
  //The swapping TVI will access the underlying DATA column when accessing
  //the CORRECTED DATA. Since DATA has been removed from the MS this should
  //throw
  ASSERT_THROW(visitIterator([&]() -> void {vb_p->visCube().shape();}),
               AipsError);
}

/*
 * Gtest fixture used to test the access to the different data columns
 * when the VI/VB2 is configured with rowBlocking().
 * This class will create a synthetic MS in a temporary directory.
 */
class DataAccessRowBlockingTest : public MsFactoryTVITester
{
public:

  /*
   * Constructor: create the temporary dir and the MsFactory used later on
   * to create the MS.
   */
  DataAccessRowBlockingTest() :
    MsFactoryTVITester("test_tViiLayerFactory","DataAccessRowBlockingTest"),
    nRowBlocking_p(-1),
    lastSPWEndTime_p(-1)
  {
  }

  /*
   * Create the synthetic MS and a VisibilityIteratorImpl2 to access it.
   */
  void createVI()
  {

    // If so specified create a last SPW which finishes before the other ones.
    // Without lastSPWEndTime_p enabled, it will create 4 SPWs (the default)
    // that are observed from time 0 to time 14.
    // With lastSPWEndTime_p enable, it will create 4 SPWs (now explicitely)
    // observed in time range [0, 14] and a last SPW observed in range [0, 9]
    // The default adds 4 antennas, i.e., 6 baselines. This gives 15*6=90 rows
    // per SPW and 10*6=60 rows in the case a 5th SPW is added.
    if(lastSPWEndTime_p != -1)
    {
      int nChannels = 100;
      double frequency = 1e11;
      double frequencyDelta = 1e9;
      std::string stokes("XX YY");
      double endTimeThisSPW = 10;
      msf_p->addSpectralWindows(4);
      msf_p->addSpectralWindow("SPW5", nChannels,
                               frequency, frequencyDelta, stokes, endTimeThisSPW);
    }

    // Create the MS
    createMS();

    // Create a disk layer type VI Factory
    IteratingParameters ipar;
    VisIterImpl2LayerFactory diskItFac(ms_p.get(),ipar,false);

    // Create a layered factory with all the layers of factories
    std::vector<ViiLayerFactory*> factories(1);
    factories[0]=&diskItFac;

    instantiateVI(factories, nRowBlocking_p);
  }

  // Destructor
  ~DataAccessRowBlockingTest()
  {
  }

  rownr_t nRowBlocking_p;

  double lastSPWEndTime_p;

};

/*
 * This test will access data when row blocking is configured
 */
TEST_F(DataAccessRowBlockingTest, AccessDataRowBlocking)
{
  // Set rowBlocking
  nRowBlocking_p = 7;

  // This creates a MS with a total of 360 rows in four SPWs.
  // The VI created to iterate that MS
  createVI();

  // Traverse the iterator accessing the proper number of rows.
  // There are 4 chunks (one per SPW), each with 90 rows. With rowBlocking
  // the number of rows in each subchunk is actually nRowBlocking+1
  size_t nRowsSubchunk = 90;
  size_t nSubChunk = nRowsSubchunk / (nRowBlocking_p + 1) + 1;
  size_t iSubChunk = 0;
  visitIterator([&]() -> void {
    if(iSubChunk < nSubChunk - 1)
      ASSERT_EQ(vb_p->nRows() , nRowBlocking_p + 1);
    else // Last subchunk has less rows than nRowBlocking_p + 1
      ASSERT_EQ(vb_p->nRows() , nRowsSubchunk - (nSubChunk - 1 ) * (nRowBlocking_p + 1));
    iSubChunk++;
    if(iSubChunk == nSubChunk)
      iSubChunk = 0;
    vb_p->visCube().shape();
  });

}

/*
 * This test will access data when row blocking is configured
 * and the chunk size of the last chunk is smaller than the rest.
 */
TEST_F(DataAccessRowBlockingTest, AccessDataRowBlockingDecreasingChunkSize)
{
  // Set rowBlocking
  nRowBlocking_p = 7;

  // Create a SPW with fewer timestamps
  lastSPWEndTime_p = 10;

  // This creates a total of 420 rows in five SPWs.
  createVI();

  // Traverse the iterator accessing the proper number of rows.
  // There are 5 chunks (one per SPW), 4 with 90 rows and a last one with 60.
  // With rowBlocking the number of rows in each subchunk is
  // actually nRowBlocking+1
  size_t iSubChunk = 0, iChunk = 0;
  size_t nChunk = 5;
  visitIterator([&]() -> void {
    size_t nRowsSubchunk = iChunk < nChunk - 1 ? 90 : 60;
    size_t nSubChunk = nRowsSubchunk / (nRowBlocking_p + 1) + 1;
    if(iSubChunk < nSubChunk - 1)
      ASSERT_EQ(vb_p->nRows() , nRowBlocking_p + 1);
    else // Last subchunk has less rows than nRowBlocking_p + 1
      ASSERT_EQ(vb_p->nRows() , nRowsSubchunk - (nSubChunk - 1 ) * (nRowBlocking_p + 1));
    iSubChunk++;
    if(iSubChunk == nSubChunk)
    {
      iSubChunk = 0;
      iChunk++;
    }
    vb_p->visCube().shape();
  });

}

/*
 * This class is a simplistic TVI that modifies the subtables
 * antenna, spw and dd
 */
class SubtableChangerTVI : public TransformingVi2
{
public:

  //Constructor
  SubtableChangerTVI(ViImplementation2 * inputVii) :
    TransformingVi2 (inputVii)
  {
    setVisBuffer(createAttachedVisBuffer (VbRekeyable));
    resetSubtables();
  }

  void origin()
  {
    // Drive underlying ViImplementation2
    getVii()->origin();

    // Synchronize own VisBuffer
    configureNewSubchunk();
  }

  void next()
  {
    // Drive underlying ViImplementation2
    getVii()->next();

    // Synchronize own VisBuffer
    configureNewSubchunk();
  }

  void
  originChunks(Bool forceRewind) override
  {
    // Drive underlying ViImplementation2
    getVii()->originChunks(forceRewind);

    // Potentially the new chunk can be from a different MS
    resetSubtables();
  }

  void
  nextChunk() override
  {
    // Drive underlying ViImplementation2
    getVii()->nextChunk();

    // Potentially the new chunk can be from a different MS
    resetSubtables();
  }

  void resetSubtables()
  {
    // Note that the creation of a new subtables is done using a
    // copy of the original subtables. However, to access these we
    // need to use the method table() of a given column (e. g. name() )
    // It would be better if the MSAntennaColumns object had
    // an getter to the MSAntenna object.
    // The same applies to the other subtables

    // Create antenna subtable
    auto& underlyingAntennaSubtablecols = getVii()->antennaSubtablecols();
    auto underlyingAntennaSubtable = underlyingAntennaSubtablecols.name().table();
    newAntennaSubtable_p = underlyingAntennaSubtable.copyToMemoryTable("SubtableChangerAntennaSubtable");
    newAntennaSubtablecols_p.reset(new MSAntennaColumns(newAntennaSubtable_p));
    // Add one antenna
    newAntennaSubtable_p.addRow();

    // Create DD subtable
    auto& underlyingDDSubtablecols = getVii()->dataDescriptionSubtablecols();
    auto underlyingDDSubtable = underlyingDDSubtablecols.spectralWindowId().table();
    newDDSubtable_p = underlyingDDSubtable.copyToMemoryTable("SubtableChangerDDSubtable");
    newDDSubtablecols_p.reset(new MSDataDescColumns(newDDSubtable_p));
    // Double the rows
    auto nrowDD = newDDSubtable_p.nrow();
    for(size_t irow = 0 ; irow < nrowDD; irow++)
      newDDSubtable_p.addRow();

    // Create spw subtable
    auto& underlyingSPWSubtablecols = getVii()->spectralWindowSubtablecols();
    auto underlyingSPWSubtable = underlyingSPWSubtablecols.name().table();
    newSPWSubtable_p = underlyingSPWSubtable.copyToMemoryTable("SubtableChangerSPWSubtable");
    newSPWSubtablecols_p.reset(new MSSpWindowColumns(newSPWSubtable_p));
    // Double the rows
    auto nrowSPW = newSPWSubtable_p.nrow();
    for(size_t irow = 0 ; irow < nrowSPW; irow++)
      newSPWSubtable_p.addRow();
  }


  const casacore::MSAntennaColumns& antennaSubtablecols() const override
  {
    return *newAntennaSubtablecols_p;
  }

  const casacore::MSDataDescColumns& dataDescriptionSubtablecols() const override
  {
    return *newDDSubtablecols_p;
  }

  const casacore::MSSpWindowColumns& spectralWindowSubtablecols() const override
  {
    return *newSPWSubtablecols_p;
  }

private:

  casacore::MSAntenna newAntennaSubtable_p;
  std::unique_ptr<casacore::MSAntennaColumns> newAntennaSubtablecols_p;

  casacore::MSSpectralWindow newSPWSubtable_p;
  std::unique_ptr<casacore::MSSpWindowColumns> newSPWSubtablecols_p;

  casacore::MSDataDescription newDDSubtable_p;
  std::unique_ptr<casacore::MSDataDescColumns> newDDSubtablecols_p;

};

/*
 * Factory that allows the creation of SubtableChangerTVI classes.
 * This factory doesn't have any parameter to configure
 */
class SubtableChangerTVILayerFactory : public ViiLayerFactory
{

public:

  SubtableChangerTVILayerFactory()
  {
  }

  virtual ~SubtableChangerTVILayerFactory() {};

protected:

  virtual ViImplementation2 * createInstance(ViImplementation2* vii0) const
  {
    ViImplementation2 *vii = new SubtableChangerTVI(vii0);
    return vii;
  }
};

/*
 * Gtest fixture used to test the creation of subtables by a TVI.
 * This class will create a synthetic MS in a temporary directory.
 * It will also create a stack of TVIs with two TVIs:
 * a disk access layer and a TVI that modifies subtables (SubtableChangerTVI).
 * The function checkSubtables can actually check that the subtables
 * have been changed as expected.
 */
class SubtableChangerTest : public MsFactoryTVITester
{
public:

  /*
   * Constructor: create the temporary dir and the MsFactory used later on
   * to create the MS.
   */
  SubtableChangerTest() : 
    MsFactoryTVITester("tViiLayerFactory","SubtableChangerTest")
  {
  }

  /*
   * Create the synthetic MS and the TVI stack to access it.
   */
  void createTVIs()
  {
    // Set the number of antennas and SPWs for the MS generated on disk
    nAntennas = 10;
    nSPWs = 10;

    msf_p->addAntennas(nAntennas);
    msf_p->addSpectralWindows(nSPWs);

    // Create synthethic MS using the msf_p factory
    createMS();

    // Create a disk layer type VI Factory
    IteratingParameters ipar;
    VisIterImpl2LayerFactory diskItFac(ms_p.get(),ipar,false);

    // Create a SubtableChangerTVI Factory
    SubtableChangerTVILayerFactory subtableChangerFac;

    // Create a layered factory with all the layers of factories
    size_t nFac = 2;
    std::vector<ViiLayerFactory*> facts(nFac);
    facts[0]=&diskItFac;
    facts[1]= &subtableChangerFac;

    // Finally create the top VI
    instantiateVI(facts);
  }

  void checkSubtables()
  {
    // Check the antenna tables size (the TVI has added one antenna)
    EXPECT_EQ(nAntennas + 1, vi_p->antennaSubtablecols().nrow());

    // Check the SPW tables size (the TVI has doubled the number)
    EXPECT_EQ(nSPWs * 2, vi_p->spectralWindowSubtablecols().nrow());

    // Check the DD tables size (the TVI has doubled the number)
    EXPECT_EQ(nSPWs * 2, vi_p->dataDescriptionSubtablecols().nrow());
  }

  // The number of antennas originally created
  size_t nAntennas;
  // The number of SPWs originally created
  size_t nSPWs;
};


/*
 * This test will simply access the corrected data column of a
 * synthetic created MS
 */
TEST_F(SubtableChangerTest, CheckSubtables)
{

  createTVIs();

  checkSubtables();
}

/*
 * Gtest fixture used to test defining sorting in a full generic way
 * for both the outer (chunk) and inner (subchunk) loops.
 * This class will create a synthetic MS in a temporary directory.
 */
class FullSortingDefinitionTest : public MsFactoryTVITester
{
public:

  /*
   * Constructor: create the temporary dir and the MsFactory used later on
   * to create the MS.
   */
  FullSortingDefinitionTest() :
    MsFactoryTVITester("tViiLayerFactory","FullSortingDefinitionTest")
  {
  }

  /*
   * Set the sortign functions for both the outer loop (chunk)
   * and inner loop (subchunk)
   */
  void setSortingDefinition(SortColumns& sortColumnsChunk,
                            SortColumns& sortColumnsSubchunk)
  {
    sortColumnsChunk_p = sortColumnsChunk;
    sortColumnsSubchunk_p = sortColumnsSubchunk;
  }

  /*
   * Create the synthetic MS and the TVI stack, which only contains
   * the disk layer VI,  to access it.
   */
  void createTVIs()
  {
    // Set the number of antennas and SPWs for the MS generated on disk
    size_t nAntennas = 10;
    nSPWs_p = 10;

    msf_p->addAntennas(nAntennas);
    msf_p->addSpectralWindows(nSPWs_p);

    // Create synthethic MS using the msf_p factory
    createMS();

    // Create a disk layer type VI Factory
    std::auto_ptr<VisIterImpl2LayerFactory> diskItFac;
    diskItFac.reset(new VisIterImpl2LayerFactory (ms_p.get(),sortColumnsChunk_p,
                                                  sortColumnsSubchunk_p,false));

    // Create a layered factory with all the layers of factories
    size_t nFac = 1;
    std::vector<ViiLayerFactory*> facts(nFac);
    facts[0]=diskItFac.get();

    // Finally create the top VI
    instantiateVI(facts);
  }

  // Sorting definitions for both inner and outer loop
  SortColumns sortColumnsChunk_p;
  SortColumns sortColumnsSubchunk_p;

  // The number of SPWs originally created
  size_t nSPWs_p;
};

/*
 * Comparison function that groups DDIds in "bins"
 */
class DDIdGroupCompare : public BaseCompare
{
public:
  explicit DDIdGroupCompare(size_t groupBin, bool ascending = true) : groupBin_p(groupBin), orderFactor_p(2*ascending-1) {}
  virtual ~DDIdGroupCompare() {}
  virtual int comp(const void * obj1, const void * obj2) const;
private:
  size_t groupBin_p;
  int orderFactor_p;
};

int DDIdGroupCompare::comp(const void * obj1, const void * obj2) const
{
  Int v1 = *(const Int*)obj1;
  Int v2 = *(const Int*)obj2;

  // The DDIds are binned in bins with a width of groupBin_p.
  // Simple integer division gives you that.
  Int bin1 = v1 / groupBin_p;
  Int bin2 = v2 / groupBin_p;

  // orderFactor_p is 1 or -1 depending on whether it has been requested
  // ascending or descending order.
  return (bin1==bin2 ? 0 : (bin1<bin2 ? -1*orderFactor_p : 1*orderFactor_p));
}

TEST_F(FullSortingDefinitionTest, GroupSPWs)
{
  // Create sorting functions for the chunk (grouping DDiD/SPWs)
  // and the subchunk (nothing).
  SortColumns sortColumnsChunk(false);
  SortColumns sortColumnsSubchunk(false);
  size_t spwBin = 2;
  CountedPtr<DDIdGroupCompare> cmpFunc(new DDIdGroupCompare(spwBin));
  sortColumnsChunk.addSortingColumn(MS::DATA_DESC_ID, cmpFunc);

  setSortingDefinition(sortColumnsChunk, sortColumnsSubchunk);

  createTVIs();

  Int spwGroup = 0;
  visitIterator([&]() -> void {
    // Get all SPWs for all rows of the VisBuffer and retain only the unique ones
    std::list<Int> spws (vb_p->spectralWindows().begin(), vb_p->spectralWindows().end());
    spws.sort();
    spws.unique();
    // The unique SPWs is equal the "bin" size. Note that the total
    // SPWs in this case is proportional to the bin size.
    ASSERT_EQ(spws.size(), spwBin);
    int i = 0;
    for(auto spw : spws)
    {
      // Check that for each "bin" of SPWs the id is as expected.
      ASSERT_EQ(Int(spwGroup * spwBin + i), spw);
      ++i;
    }
    spwGroup++;
    });

}

TEST_F(FullSortingDefinitionTest, GroupSPWsDescending)
{
  // Create sorting functions for the chunk (grouping DDiD/SPWs 
  // in descending order) and the subchunk (nothing).
  SortColumns sortColumnsChunk(false);
  SortColumns sortColumnsSubchunk(false);
  size_t spwBin = 2;
  CountedPtr<DDIdGroupCompare> cmpFunc(new DDIdGroupCompare(spwBin, false));
  sortColumnsChunk.addSortingColumn(MS::DATA_DESC_ID, cmpFunc);

  setSortingDefinition(sortColumnsChunk, sortColumnsSubchunk);

  createTVIs();

  // The first subchunk will have the largest SPWs Ids.
  Int spwGroup = nSPWs_p / spwBin - 1;
  visitIterator([&]() -> void {
    std::list<Int> spws (vb_p->spectralWindows().begin(), vb_p->spectralWindows().end());
    spws.sort();
    spws.unique();
    ASSERT_EQ(spws.size(), spwBin);
    int i = 0;
    for(auto spw : spws)
    {
      ASSERT_EQ(Int(spwGroup * spwBin + i), spw);
      ++i;
    }
    spwGroup--;
    });

}

TEST_F(FullSortingDefinitionTest, DDIdSortingBothInnerOuterLoop)
{
  // This test will create a sorting function for DDId and an empty
  // one. With loop=0 the DDId sorting function is applied to the chunks,
  // while with loop=1 it is applied to the subchunks. Both cases should
  // give the same results
  for (int loop = 0; loop < 1; loop++)
  {
    SortColumns sortColumns1(false);
    SortColumns sortColumns2(false);
    CountedPtr<ObjCompare<Int>> cmpFunc(new ObjCompare<Int>());
    sortColumns1.addSortingColumn(MS::DATA_DESC_ID, cmpFunc);

    if (loop==0)
      setSortingDefinition(sortColumns1, sortColumns2);
    else 
      setSortingDefinition(sortColumns2, sortColumns1);

    createTVIs();

    visitIterator([&]() -> void {
      std::list<Int> spws (vb_p->spectralWindows().begin(), vb_p->spectralWindows().end());
      spws.sort();
      spws.unique();
      ASSERT_EQ(spws.size(), (size_t)1);
      });
  }
}

/*
 * Comparison function that groups DDIds in "bins"
 */
class TimeDescendingCompare : public BaseCompare
{
public:
  explicit TimeDescendingCompare() {}
  virtual ~TimeDescendingCompare() {}
  virtual int comp(const void * obj1, const void * obj2) const;
};

int TimeDescendingCompare::comp(const void * obj1, const void * obj2) const
{
  return (*(const Double*)obj1  > *(const Double*)obj2  ?  -1 :
    (*(const Double*)obj1 == *(const Double*)obj2  ?  0 : 1));
}

TEST_F(FullSortingDefinitionTest, InnerLoopGroupingTimeChunksDescending)
{
  // This tests how to get subchunks in descending order.
  // The test indirectly demonstrates that time sorting in subchunk is
  // guaranteed (contrary to implementation before to CAS-12879)
  SortColumns sortColumnsChunk(false);
  SortColumns sortColumnsSubchunk(false);
  CountedPtr<TimeDescendingCompare> cmpFunc(new TimeDescendingCompare());
  sortColumnsSubchunk.addSortingColumn(MS::TIME, cmpFunc);

  setSortingDefinition(sortColumnsChunk, sortColumnsSubchunk);

  createTVIs();

  Double prevTime = -1;
  visitIterator([&]() -> void {
    std::list<Int> times (vb_p->time().begin(), vb_p->time().end());
    times.sort();
    times.unique();
    // A single timestamp per subchunk
    ASSERT_EQ(times.size(), (size_t)1);
    // Timestamp is lower than previous one.
    if(prevTime != -1.)
      ASSERT_LT(times.front(), prevTime);
    prevTime = times.front();
    });

}

// Comparison function that allows resetting the start of interval
class CompareTimeInterval : public BaseCompare
{
public:
  // Construct from the given interval values.
  CompareTimeInterval(Double interval, Double start) :
    itsInterval(interval), itsStart(start)
  {}
  
  virtual ~CompareTimeInterval() {};
  
  // Compare the interval the left and right value belong to.
  virtual int comp(const void * obj1, const void * obj2) const
  {
    Double v1 = *static_cast<const Double*>(obj1);
    Double v2 = *static_cast<const Double*>(obj2);
    // Shortcut if values are equal.
    if (v1 == v2) return 0;
    // The times are binned in bins with a width of interval_p.
    Double t1 = std::floor((v1 - itsStart) / itsInterval);
    Double t2 = std::floor((v2 - itsStart) / itsInterval);
    return (t1==t2  ?  0 : (t1<t2 ? -1 : 1));
  }

  void setStart(Double newStart) { 
  itsStart = newStart;}

private: 
  Double itsInterval;
  Double itsStart;
};

TEST_F(FullSortingDefinitionTest, ResetTimeStartInnerLoop)
{
  // This test demonstrates how to modify the sorting function 
  // before the subchunk loop.
  // A use case for that is time chunking resetting the bin 
  // at the beginning of each chunk.
  // This test uses a timeInterval of 1.1 in the inner loop and sets
  // the starting of interval to 0.55. That way a single timestamp 
  // will be included in each subchunk (time bin is 1 second) 
  // Without further support, the 5th interval from 4.95s to 6.05s 
  // will contain 2 timestamps. However, in this test, after the second chunk
  // (in 5 seconds interval), the starting of interval is reset so that
  // the fith interval is now 4.45 - 5.55 and contains again one single timestamp
  
  SortColumns sortColumnsChunk(false);
  SortColumns sortColumnsSubchunk(false);
  Double timeIntervalSubchunk = 1.1, timeStartSubchunk = -timeIntervalSubchunk/2;
  Double timeIntervalChunk = 5, timeStartChunk = 0;
  CountedPtr<CompareTimeInterval> cmpFuncSubchunk(new CompareTimeInterval(timeIntervalSubchunk, timeStartSubchunk));
  CountedPtr<CompareTimeInterval> cmpFuncChunk(new CompareTimeInterval(timeIntervalChunk, timeStartChunk));
  sortColumnsSubchunk.addSortingColumn(MS::TIME, cmpFuncSubchunk);
  sortColumnsChunk.addSortingColumn(MS::TIME, cmpFuncChunk);
  setSortingDefinition(sortColumnsChunk, sortColumnsSubchunk);

  createTVIs();

  for (vi_p->originChunks (); vi_p->moreChunks(); vi_p->nextChunk())
  {
    // Here is where the interval start is reset before subchunk loop starts.
    cmpFuncSubchunk->setStart(timeStartSubchunk);
    timeStartSubchunk += timeIntervalChunk;
    for (vi_p->origin(); vi_p->more (); vi_p->next())
    {
      std::list<Int> times (vb_p->time().begin(), vb_p->time().end());
      times.sort();
      times.unique();
      // A single timestamp per subchunk
      ASSERT_EQ(times.size(), (size_t)1);
    }
  }
}

TEST_F(FullSortingDefinitionTest, NoSortingFuncAtAll)
{
  // If no sorting functions are defined all rows are 
  // collected in a single subchunk
  SortColumns sortColumnsChunk(false);
  SortColumns sortColumnsSubchunk(false);

  setSortingDefinition(sortColumnsChunk, sortColumnsSubchunk);

  createTVIs();

  size_t totalSubchunks = 0;
  visitIterator([&]() -> void {totalSubchunks++;});
  ASSERT_EQ(totalSubchunks, (size_t)1);
}

class SeveralShapesInSubchunkTest : public MsFactoryTVITester {

public:

    SeveralShapesInSubchunkTest() :
        MsFactoryTVITester("SeveralShapesInSubchunkTest","MSFactoryCreated")
    {
    }

    void createTVIs()
    {
        //Setting the parameters to generate a synthetic MS 
        double timeInterval = 10;
        double timeSpan = 100;
        msf_p->setTimeInfo (0, timeSpan, timeInterval);
        msf_p->addAntennas(6);
        //Adding two spectral windows with different number of channels
        int nChannels = 100;
        double frequency = 1e11;
        double frequencyDelta = 1e9;
        std::string stokes("XX YY");
        msf_p->addSpectralWindow("SPW0", nChannels,
                                 frequency, frequencyDelta, stokes);
        nChannels = 50;
        frequency = 2e11;
        msf_p->addSpectralWindow("SPW1", nChannels,
                                 frequency, frequencyDelta, stokes);
        msf_p->addFeeds (10); //Need antenna and spw to be set up first

        //Creates the synthetic MS 
        createMS();

        //Create a VI Factory that access directly the MS (or the selected MS)
        //In order to have several SPWs in the same subchunk the DDID cannot be
        //in the sorting columns
        SortColumns sortColumnsChunk(false);
        SortColumns sortColumnsSubchunk(false);
        sortColumnsChunk.addSortingColumn(MS::ARRAY_ID, nullptr);
        sortColumnsChunk.addSortingColumn(MS::STATE_ID, nullptr);
        sortColumnsChunk.addSortingColumn(MS::FIELD_ID, nullptr);
        sortColumnsSubchunk.addSortingColumn(MS::TIME, nullptr);
        sortColumnsSubchunk.addSortingColumn(MS::ANTENNA1, nullptr);
        sortColumnsSubchunk.addSortingColumn(MS::ANTENNA2, nullptr);
        std::unique_ptr<VisIterImpl2LayerFactory> diskItFac;
        diskItFac.reset(new VisIterImpl2LayerFactory(ms_p.get(),sortColumnsChunk,
                                                     sortColumnsSubchunk, false));

        //Create a layered factory with all the layers of factories
        std::vector<ViiLayerFactory*> factories;
        factories.push_back(diskItFac.get());
        //Finally create the VI resulting from all the layered TVIs
        instantiateVI(factories);
    }

};

TEST_F(SeveralShapesInSubchunkTest, CheckVbColumns)
{
    createTVIs();
    size_t nRows = 0;
    //Check that after averaging with chanbin=5 the first spectral window has 20
    //channels and the second spectral window has 10 channel
    visitIterator([&]() -> void {//auto shape = vb_p->visCube().shape();
                                 //ASSERT_EQ(allEQ(vb_p->spectralWindows(), 2), true);
                                 ASSERT_EQ(vb_p->nChannelsPerShape(), IPosition(2, 100, 50));
                                 ASSERT_EQ(vb_p->nShapes(), 2u);
                                 ASSERT_EQ(vi_p->nSpectralWindows(), 2);
                                 ASSERT_EQ(vb_p->visCubes().size(),2u);
                                 ASSERT_EQ(vb_p->visCubes()[0].shape(),IPosition(3, 2, 100, 1));
                                 ASSERT_EQ(vb_p->visCubes()[1].shape(),IPosition(3, 2, 50, 1));
                                 ASSERT_EQ(vb_p->visCubesCorrected().size(),2u);
                                 ASSERT_EQ(vb_p->visCubesCorrected()[0].shape(),IPosition(3, 2, 100, 1));
                                 ASSERT_EQ(vb_p->visCubesModel().size(),2u);
                                 ASSERT_EQ(vb_p->visCubesModel()[0].shape(),IPosition(3, 2, 100, 1));
                                 ASSERT_EQ(vb_p->visCubesModel()[1].shape(),IPosition(3, 2, 50, 1));
                                 ASSERT_EQ(vb_p->sigmas().size(),2u);
                                 ASSERT_EQ(vb_p->sigmas()[0].shape(),IPosition(2, 2, 1));
                                 ASSERT_EQ(vb_p->sigmas()[1].shape(),IPosition(2, 2, 1));
                                 ASSERT_EQ(vb_p->sigmaSpectra().size(),2u);
                                 ASSERT_EQ(vb_p->sigmaSpectra()[0].shape(),IPosition(3, 2, 100, 1));
                                 ASSERT_EQ(vb_p->sigmaSpectra()[1].shape(),IPosition(3, 2, 50, 1));
                                 ASSERT_EQ(vb_p->weights().size(),2u);
                                 ASSERT_EQ(vb_p->weights()[0].shape(),IPosition(2, 2, 1));
                                 ASSERT_EQ(vb_p->weights()[1].shape(),IPosition(2, 2, 1));
                                 ASSERT_EQ(vb_p->weightSpectra().size(),2u);
                                 ASSERT_EQ(vb_p->weightSpectra()[0].shape(),IPosition(3, 2, 100, 1));
                                 ASSERT_EQ(vb_p->weightSpectra()[1].shape(),IPosition(3, 2, 50, 1));
                                 ASSERT_EQ(vb_p->flagCubes().size(),2u);
                                 ASSERT_EQ(vb_p->flagCubes()[0].shape(),IPosition(3, 2, 100, 1));
                                 ASSERT_EQ(vb_p->flagCubes()[1].shape(),IPosition(3, 2, 50, 1));
                                 nRows+=vb_p->visCubes()[0].shape()[2] + vb_p->visCubes()[1].shape()[2];
                                 });

    //All the original rows divided by 2
    size_t expectedRows = 300;

    //Check that we get the number of expected rows for each spw
    ASSERT_EQ(nRows, expectedRows);
}

