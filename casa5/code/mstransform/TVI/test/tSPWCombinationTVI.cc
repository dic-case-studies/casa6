//# tChannelAverageTVI: This file contains the unit tests of the ChannelAverageTVI class.
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


#include <limits>
#include <msvis/MSVis/SimpleSimVi2.h>
#include <msvis/MSVis/PassThroughTVI.h>
#include <mstransform/TVI/SPWCombinationTVI.h>
#include <msvis/MSVis/test/TestUtilsTVI.h>
#include <msvis/MSVis/test/MsFactory.h>

using namespace std;
using namespace casa;
using namespace casacore;
using namespace casa::vi;


TEST(SPWCombinationTVIExecuteSimulatedTest, UniformMS)
{
    //for(int nAnt = 4; nAnt <= 8; nAnt*=2)
    for(int nAnt = 4; nAnt <= 4; nAnt*=2)
    {
        //for(int nTimePerField = 1; nTimePerField <= 10; nTimePerField*=10)
        for(int nTimePerField = 1; nTimePerField <= 1; nTimePerField*=10)
        {
            //for(int nChannSpw0 = 64; nChannSpw0 <= 512; nChannSpw0*=2)
            for(int nChannSpw0 = 64; nChannSpw0 <= 64; nChannSpw0*=2)
            {
                //for(int nSpw = 1; nSpw <= 4; nSpw++)
                for(int nSpw = 1; nSpw <= 2; nSpw++)
                {
                    for(int nChannExponent = -2; nChannExponent <= 2; nChannExponent++)
                    {
                        Vector<Int> channelsPerSPW(nSpw);
                        iota(channelsPerSPW.begin(), channelsPerSPW.end(), 0);
                         
                        int nScan = 1;
                        int nCorr = 4;
                        int nField = 1;
                        double channWidth=1.0e6;
                        double deltaTime=1.0;
                        Complex visValue(1.0, 2.0);
                        Vector<Double> refFreqs(nSpw);
                        refFreqs[0] = 100.0e9;
                        for(Int iFreq = 1; iFreq < nSpw ; ++iFreq)
                            refFreqs[iFreq] = refFreqs[iFreq-1] + channelsPerSPW[iFreq-1] * channWidth;
                        SCOPED_TRACE(string("SPW combination with nSpw=") + 
                                            to_string(nSpw) + " nScan=" + to_string(nScan) +
                                            " nAnt " + to_string(nAnt) + " nCorr=" + to_string(nCorr) +
                                            " nTimePerField " + to_string(nTimePerField) + 
                                            " nChan=" + to_string(channelsPerSPW));
                        //Generating a uniform simulated Vi2 
                        SimpleSimVi2Parameters simParam(nField, nScan, 
                                                        nSpw, nAnt, nCorr,
                                                        Vector<Int>(nField, nTimePerField),
                                                        channelsPerSPW, "2016/01/06/00:00:00.",
                                                        deltaTime,
                                                        refFreqs,
                                                        Vector<Double>(nSpw, channWidth),
                                                        Matrix<Float>(1,1,1.0),
                                                        false,
                                                        Matrix<Float>(1,1,1.0),
                                                        Matrix<Float>(1,1,1.0),
                                                        false,
                                                        "circ",
                                                        false,
                                                        visValue,
                                                        false,
                                                        RowScope,
                                                        SubchunkScope);

                        SimpleSimVi2Factory simFactory(simParam);
                        VisibilityIterator2 *simVi = new VisibilityIterator2(simFactory);
                                
                        //Chaining a SPWCombinationTVI 
                        //after the simulated Vi2
                        SPWCombinationTVIFactory testFactory(simVi->getImpl());
                        VisibilityIterator2 testTVI(testFactory);

                        //Generating a simulated Vi2 with the  
                        //expected result (which is also uniform)
                        int nCombinedChannels = accumulate(channelsPerSPW.begin(), channelsPerSPW.end(), 0);
                        Vector<Double> refFreqCombined(1);
                        refFreqs[0] = 100.0e9;
                        SimpleSimVi2Parameters simResultParam(nField, nScan, 
                                                              1, nAnt, nCorr, 
                                                              Vector<Int>(nField, nTimePerField),
                                                              Vector<Int>(1, nCombinedChannels), "2016/01/06/00:00:00.",
                                                              deltaTime,
                                                              refFreqCombined,
                                                              Vector<Double>(1, channWidth),
                                                              Matrix<Float>(1,1,1.0),
                                                              false,
                                                              Matrix<Float>(1,1,1.0),
                                                              Matrix<Float>(1,1,1.0),
                                                              false,
                                                              "circ",
                                                              false,
                                                              visValue,
                                                              false,
                                                              RowScope,
                                                              SubchunkScope);

                        SimpleSimVi2Factory simResultFactory(simResultParam);
                        VisibilityIterator2 *simResultVi = new VisibilityIterator2(simResultFactory);

                        // Determine columns to check
                        VisBufferComponents2 columns;
                        columns += VisBufferComponent2::NRows;
                        columns += VisBufferComponent2::NChannels;
                        columns += VisBufferComponent2::NCorrelations;
                        //columns += VisBufferComponent2::Frequencies;
                        //columns += VisBufferComponent2::Scan;
                        //columns += VisBufferComponent2::Correlations;
                        //columns += VisBufferComponent2::Exposure;
                        //columns += VisBufferComponent2::SpectralWindows;
                        columns += VisBufferComponent2::VisibilityCubeObserved;
                        columns += VisBufferComponent2::VisibilityCubesObserved;
                        //columns += VisBufferComponent2::VisibilityCubeCorrected;
                        //columns += VisBufferComponent2::VisibilityCubeModel;

                        // Compare the combined spws 
                        SCOPED_TRACE("Comparing transformed data for simulated uniform ms");
                        double tolerance = std::numeric_limits<double>::epsilon();
                        compareVisibilityIterators(testTVI,*simResultVi, columns, tolerance);
                        compareVisibilityIterators(testTVI,*simResultVi, [&](VisBuffer2* testVb, VisBuffer2* refVb) -> void {
                                 ASSERT_EQ(testVb->spectralWindows().tovector(),
                                           (refVb->spectralWindows()+nSpw).tovector());});
                        compareVisibilityIterators(testTVI,*simResultVi, [&](VisBuffer2* testVb, VisBuffer2* refVb) -> void {
                                 ASSERT_EQ(testVb->getChannelNumbers(0).tovector(),
                                           refVb->getChannelNumbers(0).tovector());});
                    }
                }
            }
        }
    }
}


//////////////////////////////////////////////////////////////////////////
// main
//////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  int ret;

  ::testing::InitGoogleTest(&argc, argv);
  ret = RUN_ALL_TESTS();

  return ret;
}
