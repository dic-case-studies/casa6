//# tJJones: test generic corrections
//# Copyright (C) 2013,2017,2021
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

#include <casa/Arrays/Array.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Exceptions/Error.h>
#include <casa/iostream.h>
#include <casa/BasicMath/Math.h>

#include <synthesis/MeasurementComponents/StandardVisCal.h>
#include <synthesis/MeasurementComponents/SolveDataBuffer.h>
#include <synthesis/MeasurementComponents/MSMetaInfoForCal.h>
#include <msvis/MSVis/SimpleSimVi2.h>
#include <msvis/MSVis/VisBuffer2.h>

#include <gtest/gtest.h>

#include "VisCalTestBase_GT.h"

using namespace std;
using namespace casa;
using namespace casacore;
using namespace casa::vi;

// Control verbosity
#define JJONES_TEST_VERBOSE false

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

class JJonesTest : public VisCalTestBase {

public:

  Cube<Complex> x;

  JJonesTest() :
    VisCalTestBase(1,1,1,4,4,8,1,true,"lin"),
    x(4,1,nAnt,Complex(0.0))
  {
    // Create a JJones that swaps the polarizations for antenna 1
    for (Int iant=0; iant<nAnt; iant++) {
      if (iant == 1) {
	x(1,0,iant) = Complex(1.0);
	x(2,0,iant) = Complex(1.0);
      } else {
	x(0,0,iant) = Complex(1.0);
	x(3,0,iant) = Complex(1.0);
      }
    }

    if (JJONES_TEST_VERBOSE)
      summary("JJonesTest");
  }
};


TEST_F(JJonesTest, ApplyTest) {

  JJones Japp(msmc);
  Japp.setApply();

  for (Int ispw=0;ispw<nSpw;++ispw) {
    Japp.setMeta(0,0,0.0,ispw,ssvp.freqs(ispw),nChan);
    Japp.sizeApplyParCurrSpw(nChan);
    Japp.setApplyParCurrSpw(x,true,false);
  }

  for (vi2.originChunks();vi2.moreChunks();vi2.nextChunk()) {
    for (vi2.origin();vi2.more();vi2.next()) {

      Int ispw=vb2->spectralWindows()(0);
      Int obsid(vb2->observationId()(0));
      Int scan(vb2->scan()(0));
      Double timestamp(vb2->time()(0));
      Int fldid(vb2->fieldId()(0));
      Vector<Double> freqs(vb2->getFrequencies(0));
      Vector<Int> a1(vb2->antenna1());
      Vector<Int> a2(vb2->antenna2());

      vb2->resetWeightsUsingSigma();
      vb2->setVisCubeCorrected(vb2->visCube());
      vb2->setFlagCube(vb2->flagCube());

      Japp.setMeta(obsid,scan,timestamp,
		   ispw,freqs,
		   fldid);
      Japp.correct2(*vb2,false,false,false);  // (trial?,doWtSp?,dosync?)

      Cube<Complex> vC(vb2->visCube());
      Cube<Complex> vCC(vb2->visCubeCorrected());

      for (rownr_t irow=0;irow<vb2->nRows();++irow) {
	for (Int icor=0;icor<nCorr;icor++) {
	  for (Int ich=0;ich<nChan;ich++) {

	    IPosition pos1, pos2;
	    pos1 = IPosition(3,icor,ich,irow);
	    if (a1(irow) == 1)
	      pos2 = IPosition(3,(icor^2),ich,irow);
	    else if (a2(irow) == 1)
	      pos2 = IPosition(3,(icor^1),ich,irow);
	    else
	      pos2 = IPosition(3,icor,ich,irow);
	    ASSERT_NEAR(real(vC(pos1)),real(vCC(pos2)),1e-6);
	    ASSERT_NEAR(imag(vC(pos1)),imag(vCC(pos2)),1e-6);
	  }
	}
      }
    }
  }
}


class JfJonesTest : public VisCalTestBase {

public:

  Cube<Complex> x;

  JfJonesTest() :
    VisCalTestBase(1,1,1,4,4,8,1,true,"lin"),
    x(4,nChan,nAnt,Complex(0.0))
  {
    // Create a JfJones that swaps the polarizations for channels > 4
    // of antenna 1
    for (Int iant=0; iant<nAnt; iant++) {
      for (Int ich=0; ich<nChan; ich++) {
	if (iant == 1 && ich > 4) {
	  x(1,ich,iant) = Complex(1.0);
	  x(2,ich,iant) = Complex(1.0);
	} else {
	  x(0,ich,iant) = Complex(1.0);
	  x(3,ich,iant) = Complex(1.0);
	}
      }
    }

    if (JJONES_TEST_VERBOSE)
      summary("JfJonesTest");
  }
};


TEST_F(JfJonesTest, ApplyTest) {

  JfJones Japp(msmc);
  Japp.setApply();

  for (Int ispw=0;ispw<nSpw;++ispw) {
    Japp.setMeta(0,0,0.0,ispw,ssvp.freqs(ispw),nChan);
    Japp.sizeApplyParCurrSpw(nChan);
    Japp.setApplyParCurrSpw(x,true,false);
  }

  for (vi2.originChunks();vi2.moreChunks();vi2.nextChunk()) {
    for (vi2.origin();vi2.more();vi2.next()) {

      Int ispw=vb2->spectralWindows()(0);
      Int obsid(vb2->observationId()(0));
      Int scan(vb2->scan()(0));
      Double timestamp(vb2->time()(0));
      Int fldid(vb2->fieldId()(0));
      Vector<Double> freqs(vb2->getFrequencies(0));
      Vector<Int> a1(vb2->antenna1());
      Vector<Int> a2(vb2->antenna2());

      vb2->resetWeightsUsingSigma();
      vb2->setVisCubeCorrected(vb2->visCube());
      vb2->setFlagCube(vb2->flagCube());

      Japp.setMeta(obsid,scan,timestamp,
		   ispw,freqs,
		   fldid);
      Japp.correct2(*vb2,false,false,false);  // (trial?,doWtSp?,dosync?)

      Cube<Complex> vC(vb2->visCube());
      Cube<Complex> vCC(vb2->visCubeCorrected());

      for (rownr_t irow=0;irow<vb2->nRows();++irow) {
	for (Int icor=0;icor<nCorr;icor++) {
	  for (Int ich=0;ich<nChan;ich++) {

	    IPosition pos1, pos2;
	    pos1 = IPosition(3,icor,ich,irow);
	    if (a1(irow) == 1 && ich > 4)
	      pos2 = IPosition(3,(icor^2),ich,irow);
	    else if (a2(irow) == 1 && ich > 4)
	      pos2 = IPosition(3,(icor^1),ich,irow);
	    else
	      pos2 = IPosition(3,icor,ich,irow);
	    ASSERT_NEAR(real(vC(pos1)),real(vCC(pos2)),1e-6);
	    ASSERT_NEAR(imag(vC(pos1)),imag(vCC(pos2)),1e-6);
	  }
	}
      }
    }
  }
}
