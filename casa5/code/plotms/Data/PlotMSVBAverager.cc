//# PlotMSVBAverager.cc: Implementation of PlotMSVBAverager.h
//# Copyright (C) 1996,1997,1998,1999,2000,2002,2003
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU Library General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
//# License for more details.
//#
//# You should have received a copy of the GNU Library General Public License
//# along with this library; if not, write to the Free Software Foundation,
//# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//----------------------------------------------------------------------------

#include <plotms/Data/PlotMSVBAverager.h>
#include <casa/Exceptions/Error.h>
#include <casa/Arrays/ArrayLogical.h>
#include <casa/Arrays/ArrayMath.h>
#include <measures/Measures/Stokes.h>
#include <msvis/MSVis/VisBufferComponents2.h>
#include <cfloat>

#define PRTLEV_PMSVBA -1

using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN

//----------------------------------------------------------------------------

PlotMSVBAverager::PlotMSVBAverager(Int nAnt)
  : nAnt_p(nAnt),
    nCorr_p(0),
    nChan_p(0),
    nBlnMax_p(0),
    blnOK_p(),
    blnAve_p(false),
    antAve_p(false),
    inCoh_p(false),
    timeRef_p(0.0),
    minTime_p(0.0),
    maxTime_p(0.0),
    aveTime_p(0.0),
    aveInterval_p(0.0),
    blnWtSum_p(),
    vbWtSum_p(0.0),
    aveScan_p(0),
    doVC_p(false),
    doMVC_p(false),
    doCVC_p(false),
    doFC_p(false),
    doUVW_p(false),
    doWC_p(false),
    initialized_p(false),
    jcor_p(4,0),
    prtlev_(PRTLEV_PMSVBA)
{
// Construct from the number of antennas and the averaging interval
// Input:
// Output to private data:
//    nAnt_p               Int              No. of antennas
//
  // Note: this returns a VisBufferImpl2
  avBuf_p = vi::VisBuffer2::factory(vi::VbRekeyable);

  if (prtlev()>2) cout << "PMSVBA::PMSVBA()" << endl;

  // Ensure we do weightCube if any data is done
  doWC_p = doVC_p || doMVC_p || doCVC_p || doFC_p;

  // Nominal correlation order
  indgen(jcor_p);


};

//----------------------------------------------------------------------------

PlotMSVBAverager::~PlotMSVBAverager()
{
// Null default destructor
//
  if (prtlev() > 2) cout << "PMSVBA::~PMSVBA()" << endl;
};

//----------------------------------------------------------------------------

void PlotMSVBAverager::finalizeAverage()
{
// Normalize the current accumulation interval
// Output to private data:
//    avBuf_p         vi::VisBuffer2&       Averaged buffer
//  

  if (prtlev()>2) cout  << "  PMSVBA::finalizeAverage()" << endl;

  Int nBln = ntrue(blnOK_p);
  if (nBln>0 && vbWtSum_p>0.0) {

    Int obln = 0;

    // Divide by the weights
    for (Int ibln=0; ibln<nBlnMax_p; ++ibln) {
      if (blnOK_p(ibln)) {
	for (Int ichn=0; ichn<nChan_p; ++ichn) {
	  for (Int icor=0; icor<nCorr_p; ++icor) {
	    Float& thiswt(avgWeight_(icor,ichn,ibln));
	    // normalize ibln data at obln:
	    if (doWC_p && (thiswt > 0.0)) {
	      if (doVC_p) 
		avgVisCube_(icor,ichn,obln) =
		  avgVisCube_(icor,ichn,ibln) / thiswt;
	      if (doMVC_p) 
		avgModelCube_(icor,ichn,obln) =
		  avgModelCube_(icor,ichn,ibln) / thiswt;
	      if (doCVC_p) 
		avgCorrectedCube_(icor,ichn,obln) =
		  avgCorrectedCube_(icor,ichn,ibln) / thiswt;
	      if (doFC_p) 
		avgFloatCube_(icor,ichn,obln) =
		  avgFloatCube_(icor,ichn,ibln) / thiswt;
	    } // !flagCube & wt>0
	    else {
	      // make sure obln is zero
	      if (obln<ibln) {
		if (doVC_p)  avgVisCube_(icor,ichn,obln) = 0.0;
		if (doMVC_p) avgModelCube_(icor,ichn,obln) = 0.0;
		if (doCVC_p) avgCorrectedCube_(icor,ichn,obln) = 0.0;
		if (doFC_p)  avgFloatCube_(icor,ichn,obln) = 0.0;
	      }
	    }
	    // copy flags, weights to obln
	    if (obln < ibln) {
	      avgFlagCube_(icor,ichn,obln) = avgFlagCube_(icor,ichn,ibln);
	      if (doWC_p)
		avgWeight_(icor,ichn,obln) = thiswt;
	    }
	  } // icor
	} // ichn
	 
	if (doUVW_p && (blnWtSum_p(ibln) > 0.0))
	  avgUvw_.column(obln) =
	    avgUvw_.column(ibln) / blnWtSum_p(ibln);

	if (obln<ibln) {
	  avgFlagRow_(obln) = avgFlagRow_(ibln);
	  avgAntenna1_(obln) = avgAntenna1_(ibln);
	  avgAntenna2_(obln) = avgAntenna2_(ibln);
	}

	// increment the output row
	++obln;
      } // blnOK
    } // ibln

    // Contract the VisBuffer data, if necessary
    if (nBln < nBlnMax_p) {
      avBuf_p->setShape(nCorr_p, nChan_p, nBln, true);
      avgTime_.resize(nBln);
      avgTimeInterval_.resize(nBln);
      avgScan_.resize(nBln);
      avgAntenna1_.resize(nBln, true);
      avgAntenna2_.resize(nBln, true);
      if (doUVW_p)
	avgUvw_.resize(3,nBln,true);
      if (doVC_p)
	avgVisCube_.resize(nCorr_p,nChan_p,nBln,true);
      if (doMVC_p)
	avgModelCube_.resize(nCorr_p,nChan_p,nBln,true);
      if (doCVC_p)
	avgCorrectedCube_.resize(nCorr_p,nChan_p,nBln,true);
      if (doFC_p)
	avgFloatCube_.resize(nCorr_p,nChan_p,nBln,true);
      if (doWC_p)
	avgWeight_.resize(nCorr_p,nChan_p,nBln,true);
      avgFlagRow_.resize(nBln,true);
      avgFlagCube_.resize(nCorr_p,nChan_p,nBln,true);
      fieldid_.resize(nBln, true);
      spw_.resize(nBln, true);
      obsid_.resize(nBln, true);
      stateid_.resize(nBln, true);
    }

    // Time:
    //    aveTime_p/=vbWtSum_p;
    //    avBuf_p->time()=aveTime_p;
    avgTime_.set(timeRef_p + (maxTime_p + minTime_p)/2.0);  // must be center of the interval!
    avgTimeInterval_.set(maxTime_p - minTime_p);
    avgScan_.set(aveScan_p);

    // Assign averaged data
    if (doVC_p)  avBuf_p->setVisCube(avgVisCube_);
    if (doMVC_p) avBuf_p->setVisCubeModel(avgModelCube_);
    if (doCVC_p) avBuf_p->setVisCubeCorrected(avgCorrectedCube_);
    if (doFC_p) avBuf_p->setVisCubeFloat(avgFloatCube_);
    if (doUVW_p) avBuf_p->setUvw(avgUvw_);
    if (doWC_p)  avBuf_p->setWeightSpectrum(avgWeight_);
    avBuf_p->setFlagCube(avgFlagCube_);
    avBuf_p->setFlagRow(avgFlagRow_);
    avBuf_p->setAntenna1(avgAntenna1_);
    avBuf_p->setAntenna2(avgAntenna2_);
    avBuf_p->setTime(avgTime_);
    avBuf_p->setTimeInterval(avgTimeInterval_);
    avBuf_p->setScan(avgScan_);
    avBuf_p->setFieldId(fieldid_);
    avBuf_p->setSpectralWindows(spw_);
    avBuf_p->setObservationId(obsid_);
    avBuf_p->setStateId(stateid_);
    //    cout << "final amp = " << amplitude(avBuf_p->visCube()) << endl;
    //    cout << "final flagR = " << boolalpha << avBuf_p->flagRow() << endl;
    //    cout << "final flagC = " << boolalpha << avBuf_p->flagCube() << endl;

    // We need to be reinitialized to do more accumulating
    initialized_p = false;

  }
  else {

    // No rows!
    //    avBuf_p->nRow()=0;

    //    cout << avBuf_p->nRow() << endl;
    //    cout << "nBln = " << nBln << " (" << nBlnMax_p << ")" << endl;
    //    cout << "vbWtSum_p = " << vbWtSum_p << endl;
    //    cout << "avBuf_p->flagRow() = " << avBuf_p->flagRow() << endl;
    //    cout << "avBuf_p->antenna1() = " << avBuf_p->antenna1() << endl;
    //    cout << "avBuf_p->antenna2() = " << avBuf_p->antenna2() << endl;

    // Should not be able to reach here
    // (there is always at least one good baseline, even if flagged)
    if (nBln==0) 
        throw(AipsError("Could not finalize average, no good baselines"));
    // should not happen, they are normalized to FLT_MIN!
    else if (vbWtSum_p<=0.0) 
        throw(AipsError("Could not finalize average, weights are zero"));
    else
        throw(AipsError("Big problem in finalizeAverage!"));
  }
};

//----------------------------------------------------------------------------

void PlotMSVBAverager::initialize(vi::VisBuffer2& vb)
{
  // Initialize the averager

  if (prtlev()>2) cout << "  PMSVBA::initialize()" << endl;

  // Assign main meta info
  vi::VisBufferComponents2 components = vi::VisBufferComponents2::these
      ({VisBufferComponent2::FieldId, VisBufferComponent2::ObservationId, VisBufferComponent2::SpectralWindows, VisBufferComponent2::StateId});
  avBuf_p->copyComponents(vb, components, true, false);
  fieldid_ = vb.fieldId();
  spw_ = vb.spectralWindows();
  obsid_ = vb.observationId();
  stateid_ = vb.stateId();

  // Immutables:
  nChan_p = vb.nChannels();
  nCorr_p = vb.nCorrelations();

  if (blnAve_p)
    nBlnMax_p = 1;
  else if (antAve_p)
    nBlnMax_p = nAnt_p;
  else 
    nBlnMax_p = nAnt_p * (nAnt_p + 1)/2;

  blnOK_p.resize(nBlnMax_p);
  blnOK_p = false;
  blnWtSum_p.resize(nBlnMax_p);
  blnWtSum_p = 0.0;

  /*
  cout << "Shapes = " 
       << nCorr_p << " "
       << nChan_p << " "
       << nBlnMax_p << " "
       << ntrue(blnOK_p) << " "
       << endl;
  */

  // Set basic shapes
  avBuf_p->setShape(nCorr_p, nChan_p, nBlnMax_p, true);

  // Resize and fill in the antenna numbers for all rows
  avgAntenna1_.resize(nBlnMax_p);
  avgAntenna2_.resize(nBlnMax_p);
  if (blnAve_p) {
    avgAntenna1_.set(-1);
    avgAntenna2_.set(-1);
  }
  else if (antAve_p) {
    // Ensure zero or two cross-hands present
    //  (and set jcor_p)
    verifyCrosshands(vb);
    indgen(avgAntenna1_);
    avgAntenna2_.set(-1);
  }
  else {
    Int ibln = 0;
    for (Int iant1=0; iant1<nAnt_p; ++iant1) {
      for (Int iant2=iant1; iant2<nAnt_p; ++iant2) {
	avgAntenna1_(ibln) = iant1;
	avgAntenna2_(ibln) = iant2;
	++ibln;
      }
    }
  }

  // Resize and fill if larger
  Int nRows = vb.nRows();
  fieldid_.resize(nBlnMax_p, true);
  spw_.resize(nBlnMax_p, true);
  obsid_.resize(nBlnMax_p, true);
  stateid_.resize(nBlnMax_p, true);
  if (nBlnMax_p > nRows) fillIds(nRows);
  // Resize and initialize everything else
  avgTime_.resize(nBlnMax_p, false); 
  avgTime_.set(0.0);
  avgTimeInterval_.resize(nBlnMax_p, false); 
  avgTimeInterval_.set(0.0);
  avgScan_.resize(nBlnMax_p, false); 
  aveScan_p = vb.scan()(0);

  // All _rows_ assumed UNflagged, to start with
  //  (will refine later, if needed)
  avgFlagRow_.resize(nBlnMax_p, false); 
  avgFlagRow_.set(false);

  // all cells assumed flagged to start with
  avgFlagCube_.resize(nCorr_p,nChan_p, nBlnMax_p,false);
  avgFlagCube_.set(true);

  if (doUVW_p) {
    avgUvw_.resize(3,nBlnMax_p, false); 
    avgUvw_.set(0.0);
  }

  Complex czero(0.0);
  if (doVC_p) {
    avgVisCube_.resize(nCorr_p,nChan_p, nBlnMax_p,false);
    avgVisCube_.set(czero);
  }
  if (doMVC_p) {
    avgModelCube_.resize(nCorr_p,nChan_p, nBlnMax_p,false);
    avgModelCube_.set(czero);
  }
  if (doCVC_p) {
    avgCorrectedCube_.resize(nCorr_p,nChan_p, nBlnMax_p,false);
    avgCorrectedCube_.set(czero);
  }
  if (doFC_p) {
    avgFloatCube_.resize(nCorr_p,nChan_p, nBlnMax_p,false);
    avgFloatCube_.set(0.0);
  }

  if (doWC_p) {
    avgWeight_.resize(nCorr_p,nChan_p,nBlnMax_p, false); 
    avgWeight_.set(0.0f);
  }

  minTime_p = DBL_MAX;
  maxTime_p = -DBL_MAX;
  timeRef_p = vb.time()(0);

  // we are now initialized
  initialized_p=true;
};


//----------------------------------------------------------------------------
void PlotMSVBAverager::simpAccumulate (vi::VisBuffer2& vb)
{
// Accumulate a VisBuffer
// Input:
//    vb               const vi::VisBuffer2&   VisBuffer to accumulate
// Output to private data:
//    tStart_p         Double               Start time of current accumulation
//    nChan_p          Int                  No. of channels in the avg. buffer
//    avrow_p          Int                  Start row of current accumulation
//    avBuf_p          vi::CalVisBuffer2    Averaging buffer
//  

  if (prtlev()>2) cout << " PMSVBA::accumulate() " << endl;

  if (!initialized_p) initialize(vb);

  // Only accumulate VisBuffers with the same number of channels
  //  TBD: handle multiple shapes by separating ddis?
  if ((vb.nChannels() != nChan_p) || (vb.nCorrelations() != nCorr_p)) 
    throw(AipsError("PlotMSVBAverager: data shape does not conform"));

  Double vbWt(0.0);  // will accumulate this VBs total data weight
  const Float* wt;

  // Mutable vis cubes for accumulation
  Cube<Complex> accumVisCube;
  Cube<Complex> accumVisCubeModel;
  Cube<Complex> accumVisCubeCorrected;
  Cube<Float> accumVisCubeFloat;
  if (doVC_p) {
	accumVisCube.reference(vb.visCube());
	if (inCoh_p) convertToAP(accumVisCube);
  }
  if (doMVC_p) {
	accumVisCubeModel.reference(vb.visCubeModel());
	if (inCoh_p) convertToAP(accumVisCubeModel);
  }
  if (doCVC_p) {
	accumVisCubeCorrected.reference(vb.visCubeCorrected());
	if (inCoh_p) convertToAP(accumVisCubeCorrected);
  }
  if (doFC_p) {
	accumVisCubeFloat.reference(vb.visCubeFloat());
  }

  for (Int ibln=0; ibln<vb.nRows(); ++ibln) {
    // Calculate row from antenna numbers with the hash function.
    Int ant1 = vb.antenna1()(ibln);
    Int ant2 = vb.antenna2()(ibln);
    Int obln = baseline(ant1, ant2);

    // This baseline occurs in input, so preserve in output
    blnOK_p(obln) = true;

    wt = &vb.weight()(0,ibln);
    Float vbweight = *wt;
    // Ensure weights strictly positive
    // assumes chanIndepWt_p=true
    if (vbweight < FLT_MIN) vbweight = FLT_MIN;
    Double blnWt(0.0);

    for (Int chn=0; chn<vb.nChannels(); chn++) {
      for (Int cor=0; cor<nCorr_p; ++cor,++wt) {
	
	// Assume we won't accumulate anything in this cell
	//   (output is unflagged, input is flagged)
	Bool acc(false);

    IPosition flagPos(3, cor, chn, ibln);
	if (!vb.flagCube()(flagPos)) { // input UNflagged
	  // we will accumulate
	  acc=true;
	  if (avgFlagCube_(cor,chn,obln)) {  // output flagged
	    // This cell now NEWLY unflagged in output
	    avgFlagCube_(cor,chn,obln)=false;
	    // ...so zero the accumulators
	    if (doVC_p) 
	      avgVisCube_(cor,chn,obln) = 0.0;
	    if (doMVC_p) 
	      avgModelCube_(cor,chn,obln) = 0.0;
	    if (doCVC_p)
	      avgCorrectedCube_(cor,chn,obln) = 0.0;
	    if (doFC_p) 
	      avgFloatCube_(cor,chn,obln)=0.0;
	    if (doWC_p) 
	      avgWeight_(cor,chn,obln) = 0.0;
	  }
	}
	else  // input cell is flagged
	  // Only accumulate if output is also flagged
	  //  (yields average of flagged data if no unflagged data ever found)
	  if (avgFlagCube_(cor,chn,obln)) acc=true;

	// Accumulate this (cor,chn), if appropriate
	if (acc) {
	  if (doVC_p)
		avgVisCube_(cor,chn,obln) += 
			((vbweight) * accumVisCube(cor,chn,ibln) );
	  if (doMVC_p) 
		avgModelCube_(cor,chn,obln) += 
			((vbweight) * accumVisCubeModel(cor,chn,ibln) );
	  if (doCVC_p)
		avgCorrectedCube_(cor,chn,obln) +=
			((vbweight) * accumVisCubeCorrected(cor,chn,ibln) );
	  if (doFC_p)
		avgFloatCube_(cor,chn,obln)+=
			((vbweight) * accumVisCubeFloat(cor,chn,ibln) );
	  if (doWC_p) {
		avgWeight_(cor,chn,obln) += (vbweight);
	  }
	} // acc
	blnWt += vbweight;
      } // cor 
    } // chn

    // Don't let it be too large
    blnWt /= Double(nChan_p);
    vbWt += blnWt;
    
    blnWtSum_p(obln) += blnWt;
    
    // UVW
    if (doUVW_p && blnWt>0.0)
      for (uInt i=0;i<3;++i)
	avgUvw_(i,obln) += (vb.uvw()(i,ibln) * blnWt);
    
  } // ibln
  if (vbWt>0) {

    vbWtSum_p += vbWt;

    Double thisTime(vb.time()(0) - timeRef_p);
    Double thisInterval(vb.timeInterval()(0));

    minTime_p = std::min(minTime_p, (thisTime - thisInterval/2.0));
    maxTime_p = std::max(maxTime_p, (thisTime + thisInterval/2.0));

    aveTime_p += (thisTime * vbWt);
    aveInterval_p += vb.timeInterval()(0);
 
    Int thisScan = vb.scan()(0);
    if (aveScan_p != thisScan)
      aveScan_p = -1;

    // This doesn't work...
    //    Int thisField=vb.fieldId();
    //    if (avBuf_p->fieldId()!=thisField)
    //      avBuf_p->fieldId()=-1;

  }  

  //  cout << "Done with accumulate." << endl;


};

//----------------------------------------------------------------------------

void PlotMSVBAverager::antAccumulate (vi::VisBuffer2& vb)
{
// Accumulate a VisBuffer with per-antenna averaging
// Input:
//    vb               const vi::VisBuffer2&     VisBuffer to accumulate

  if (prtlev()>2) cout << " PMSVBA::antAccumulate() " << endl;

  if (!initialized_p) initialize(vb);

  // Only accumulate VisBuffers with the same number of channels
  //  TBD: handle multiple shapes by separating ddis?
  if ((vb.nChannels() != nChan_p) || (vb.nCorrelations() != nCorr_p)) 
    throw(AipsError("PlotMSVBAverager: data shape does not conform"));

  //  if (vb.spectralWindow()!=avBuf_p->spectralWindow())
  //    avBuf_p->spectralWindow()=-1;

  Double vbWt(0.0);  // will accumulate this VBs total data weight
  Cube<Float> wtsp;
  wtsp.resize();
  wtsp = vb.weightSpectrum();
  Float wt;

  // Mutable vis cubes for accumulation
  Cube<Complex> accumVisCube;
  Cube<Complex> accumVisCubeModel;
  Cube<Complex> accumVisCubeCorrected;
  Cube<Float> accumVisCubeFloat;
  if (doVC_p) {
	accumVisCube.reference(vb.visCube());
	if (inCoh_p) convertToAP(accumVisCube);
  }
  if (doMVC_p) {
	accumVisCubeModel.reference(vb.visCubeModel());
	if (inCoh_p) convertToAP(accumVisCubeModel);
  }
  if (doCVC_p) {
	accumVisCubeCorrected.reference(vb.visCubeCorrected());
	if (inCoh_p) convertToAP(accumVisCubeCorrected);
  }
  if (doFC_p) {
	accumVisCubeFloat.reference(vb.visCubeFloat());
  }

  for (Int ibln=0; ibln<vb.nRows(); ++ibln) {

    // The antennas in the baseline
    Vector<Int> oblnij(2);
    oblnij(0) = vb.antenna1()(ibln);
    oblnij(1) = vb.antenna2()(ibln);
    Int& obln_i(oblnij(0));
    Int& obln_j(oblnij(1));

    // These antennas occur (even if flagged) in input, so preserve in output
    blnOK_p(obln_i) = true;
    blnOK_p(obln_j) = true;
      
    Double blnWt(0.0);  // will accumulate this baseline's total data weight
    for (Int chn=0; chn<vb.nChannels(); chn++) {
      for (Int cor=0; cor<nCorr_p; ++cor) {
          wt = wtsp(cor, chn, ibln);
          if (wt < FLT_MIN) wt = FLT_MIN;
	// Assume we won't accumulate anything in this cell
	//   (output is unflagged, input is flagged)
	Bool acc_i(false),acc_j(false);

	// Consider accumulation according to state of flags
    IPosition flagPos(3, cor, chn, ibln);
	if (!vb.flagCube()(flagPos)) { // input UNflagged
	  // we will accumulate both ants
	  acc_i = acc_j = true;

	  // Zero accumulators if output cell currently flagged
	  for (Int ij=0; ij<2; ++ij) {
	    Int ia = oblnij(ij);
	    if (avgFlagCube_(cor,chn,ia)) {  // output flagged
	      // This cell now NEWLY unflagged in output
	      avgFlagCube_(cor,chn,ia) = false;
	      // ...so zero the accumulators
	      if (doVC_p) 
		avgVisCube_(cor,chn,ia) = 0.0;
	      if (doMVC_p) 
		avgModelCube_(cor,chn,ia) = 0.0;
	      if (doCVC_p)
		avgCorrectedCube_(cor,chn,ia) = 0.0;
	      if (doFC_p) 
		avgFloatCube_(cor,chn,ia)=0.0;
	      if (doWC_p)
		avgWeight_(cor,chn,ia) = 0.0;
	    }
	  }
	}
	else {        // input cell is flagged
	  // Only accumulate if output is also flagged:
	  //  (yields average of flagged data if no unflagged data ever found)
	  if (avgFlagCube_(cor,chn,obln_i)) acc_i = true;
	  if (avgFlagCube_(cor,chn,obln_j)) acc_j = true;
	}
	
	// Accumulate data, if appropriate
	if (acc_i) {
		if (doVC_p) avgVisCube_(cor,chn,obln_i) +=
				((wt) * accumVisCube(cor,chn,ibln));
		if (doMVC_p) avgModelCube_(cor,chn,obln_i) +=
				((wt) * accumVisCubeModel(cor,chn,ibln));
		if (doCVC_p) avgCorrectedCube_(cor,chn,obln_i) +=
				((wt) * accumVisCubeCorrected(cor,chn,ibln));
	  	if (doFC_p)  avgFloatCube_(cor,chn,obln_i) +=
               	          	((wt) * accumVisCubeFloat(cor,chn,ibln) );
		if (doWC_p)  avgWeight_(cor,chn,obln_i) += (wt);
	}
	if (acc_j) { 
		// NB: wt is implicitly indexed by cor, so index incoming data w/ cor
		Int jcor = jcor_p(cor);  // handle cross-hand swap
	  	if (doVC_p) avgVisCube_(jcor,chn,obln_j) +=
				((wt) * conj(accumVisCube(cor,chn,ibln)));
		if (doMVC_p) avgModelCube_(jcor,chn,obln_j) +=
				((wt) * conj(accumVisCubeModel(cor,chn,ibln)));
		if (doCVC_p) avgCorrectedCube_(jcor,chn,obln_j) +=
				((wt) * conj(accumVisCubeCorrected(cor,chn,ibln)));
	  	if (doFC_p)  avgFloatCube_(jcor,chn,obln_j) +=
				((wt) * accumVisCubeFloat(cor,chn,ibln));
		if (doWC_p)  avgWeight_(jcor,chn,obln_j) += (wt);
	}

	blnWt += (wt);
	
      } // cor 
    } // chn

    // Don't let it be too large
    blnWt /= Double(nChan_p);
    vbWt += blnWt;
    
    blnWtSum_p(obln_i) += blnWt;
    blnWtSum_p(obln_j) += blnWt;
    
    // UVW
    if (doUVW_p && blnWt>0.0)
      for (uInt i=0; i<3; ++i) {
	avgUvw_(i,obln_i) += (vb.uvw()(i,ibln) * blnWt);
	avgUvw_(i,obln_j) += (vb.uvw()(i,ibln) * blnWt);
      }      
  }

  if (vbWt>0) {

    vbWtSum_p += vbWt;

    Double thisTime(vb.time()(0) - timeRef_p);
    Double thisInterval(vb.timeInterval()(0));

    minTime_p = std::min(minTime_p, (thisTime - thisInterval/2.0));
    maxTime_p = std::max(maxTime_p, (thisTime + thisInterval/2.0));

    aveTime_p += (thisTime * vbWt);
    aveInterval_p += vb.timeInterval()(0);
 
    Int thisScan = vb.scan()(0);
    if (aveScan_p != thisScan)
      aveScan_p = -1;

    // This doesn't work...
    //    Int thisField=vb.fieldId();
    //    if (avBuf_p->fieldId()!=thisField)
    //      avBuf_p->fieldId()=-1;

  }  

};

//----------------------------------------------------------------------------

void PlotMSVBAverager::verifyCrosshands(vi::VisBuffer2& vb) {

  // Nominal values for jcor_p
  jcor_p.resize(nCorr_p);
  indgen(jcor_p);
  

  // Only perform cross-hand test if antenna-based averaging turned on
  if (antAve_p) {

    Vector<Int> corrs = vb.getCorrelationTypes();

    // Detect crosshands ids
    Vector<Int> chids(2,-1);
    Int ich(0);
    for (uInt icor=0; icor<jcor_p.nelements(); ++icor) {
      if (corrs(icor) == Stokes::RL ||
	  corrs(icor) == Stokes::LR ||
	  corrs(icor) == Stokes::XY ||
	  corrs(icor) == Stokes::YX) chids(ich++) = icor;
    }

    if (chids(0) > -1) {  // at least one ch detected
      if (chids(1) < 0) // second ch not detected, so abort
	throw(AipsError("Both cross-hands (or none) must be selected and present when Per Antenna averaging is selected."));
      else {  // Both cross-hands detected

	// Cross-hand difference must be 1  (indicates basis consistency)
	Int chdiff = abs(corrs(chids(1)) - corrs(chids(0)));
	if (chdiff != 1) throw(AipsError("Cross-hand basis inconsistency detected."));
      
	// swap cross-hands in jcor_p
	Int tmp = jcor_p(chids(0));
	jcor_p(chids(0)) = jcor_p(chids(1));
	jcor_p(chids(1)) = tmp;


      }
    }
  }

}

//----------------------------------------------------------------------------

Int PlotMSVBAverager::baseline(const Int& ant1, const Int& ant2)
{
// Compute row index in an accumulation interval for an
// interferometer index (ant1, ant2).
// Input:
//    ant1            const Int&      Antenna 1
//    ant2            const Int&      Antenna 2
// Output:
//    hashFunction    Int             Row offset in current accumulation
//

  if (blnAve_p)
    return 0;

  Int index;
  index = nAnt_p * ant1 - (ant1 * (ant1 - 1)) / 2 + ant2 - ant1;
  return index;
};

//----------------------------------------------------------------------------

void PlotMSVBAverager::convertToAP(Cube<Complex>& d) {
  Int n = d.nelements();
  Complex *c = d.data();
  Float a;
  for (Int i=0; i<n; ++i,++c) {
    a = abs(*c);
    *c = Complex(a,arg(*c));
  }
}

//----------------------------------------------------------------------------
void PlotMSVBAverager::fillIds(Int nrows) {
    // when you are resizing bigger than number of rows,
    // have to fill rest of vector or get seg fault from
    // invalid value (whatever is in memory; 
    // fill with last value in id vector
    for (Int i=nrows; i<nBlnMax_p; ++i) {
        fieldid_[i] = fieldid_[nrows-1];
        spw_[i] = spw_[nrows-1];
        obsid_[i] = obsid_[nrows-1];
        stateid_[i] = stateid_[nrows-1];
    }
}
} //# NAMESPACE CASA - END

