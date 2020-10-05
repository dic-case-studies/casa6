//# PlotMSCTAverager.cc: Implementation of PlotMSCTAverager.h
//# Copyright (C) 2020
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

#include <plotms/Data/PlotMSCTAverager.h>

using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN

//----------------------------------------------------------------------------

PlotMSCTAverager::PlotMSCTAverager(
  PlotMSAveraging& averaging, Int nAnt, Int nPoln)
  : averaging_p(averaging),
    nAnt_p(nAnt),
    nPoln_p(nPoln),
    nChan_p(0),
    nBlnMax_p(0),
    isComplex_p(True),
    blnOK_p(),
    timeRef_p(0.0),
    minTimeOffset_p(0.0),
    maxTimeOffset_p(0.0),
    aveScan_p(0),
    aveField_p(0),
    blnCount_p(0),
    blnWtSum_p(),
    initialized_p(false),
    isAccum_p(false),
    debug_(false) {}

//----------------------------------------------------------------------------

PlotMSCTAverager::~PlotMSCTAverager() {
  // Null default destructor
  if (debug_) {
    cout << "PMSCTA::~PMSCTA()" << endl;
  }
}

//----------------------------------------------------------------------------

void PlotMSCTAverager::finalizeAverage() {
  // Normalize the current accumulation interval
  if (debug_) {
    cout  << "  PMSCTA::finalizeAverage()" << endl;
  }


  Int nBlnOK = ntrue(blnOK_p);
  if ((nBlnOK > 0) && isAccum_p) {
    // Increment obln for every blnOK
    Int obln(0);

    casacore::Cube<casacore::Complex> cparam(nPoln_p, nChan_p, nBlnOK);
    cparam.set(Complex(0.0));
    casacore::Cube<casacore::Float> fparam(nPoln_p, nChan_p, nBlnOK);
    fparam.set(0.0);
    casacore::Cube<casacore::Float> flag(nPoln_p, nChan_p, nBlnOK);
    flag.set(true);
    casacore::Vector<casacore::Int> antenna1(nBlnOK);
    casacore::Vector<casacore::Int> antenna2(nBlnOK);

    // Divide each baseline sum by its count
    for (Int ibln = 0; ibln < nBlnMax_p; ++ibln) {
      if (blnOK_p(ibln)) {
        // baseline was in unaveraged data
        Int thisCount(blnCount_p(ibln));

        for (Int ichan = 0; ichan < nChan_p; ++ichan) {
          for (Int ipol = 0; ipol < nPoln_p; ++ipol) {
            // normalize ibln data at obln:
            if (thisCount > 0.0) {
              if (isComplex_p) {
                cparam(ipol, ichan, obln) = accumCParam_(ipol, ichan, ibln) / (float)thisCount;
			  } else {
                fparam(ipol, ichan, obln) = accumFParam_(ipol, ichan, ibln) / (float)thisCount;
              }
            }

            // copy flags to obln
            flag(ipol, ichan, obln) = avgFlag_(ipol, ichan, ibln);
          } // icor
        } // ichn
 
        antenna1(obln) = avgAntenna1_(ibln);
        antenna2(obln) = avgAntenna2_(ibln);

        // increment the output row
        ++obln;
      } // blnOK
    } // ibln

    // Contract the data, if necessary
    if (nBlnOK < nBlnMax_p) {
      avgTime_.resize(nBlnOK);
      avgScan_.resize(nBlnOK);
      avgField_.resize(nBlnOK);
      avgSpw_.resize(nBlnOK, true);
      avgObsid_.resize(nBlnOK, true);
    }

    // Time: center of the interval
    avgTime_.set(timeRef_p + (maxTimeOffset_p + minTimeOffset_p)/2.0);
    // Value for all iterations, or -1 if they vary
    avgScan_.set(aveScan_p);
    avgField_.set(aveField_p);

    // Assign averaged data
    for (Int ibln = 0; ibln < nBlnOK; ++ibln) {
      CTMainRecord main_row_record;
      main_row_record.defineTime(avgTime_(ibln));
      main_row_record.defineFieldId(avgField_(ibln));
      main_row_record.defineSpwId(avgSpw_(ibln));
      main_row_record.defineAntenna1(antenna1(ibln));
      main_row_record.defineAntenna2(antenna2(ibln));
      main_row_record.defineInterval(0); // not supported
      main_row_record.defineScanNo(avgScan_(ibln));
      if (isComplex_p) {
        main_row_record.defineCParam(cparam(ibln));
      } else {
        main_row_record.defineFParam(fparam(ibln));
      }
      main_row_record.defineFlag(flag(ibln));
      main_row_record.defineSnr(snr(ibln));
      main_row_record.define(snr(ibln));

      _main_rows.push_back(record);
    }
    
    // We need to be reinitialized to do more accumulating
    initialized_p = false;
  } else {
    // Should not happen
    if (nBlnOK == 0) {
      // There is always at least one valid baseline, even if flagged
      throw(AipsError("Could not finalize average, no valid baselines"));
    } else {
      // No averaged rows
      throw(AipsError("Could not finalize average, no accumulations"));
    }
  }
}

NewCalTable PlotMSCTAverager::avgCalTable() {
  // Return NewCalTable filled with averaging results stored in CTMainRecord vector
  NewCalTable nct;

  for (size_t row = 0; row > main_rows_.size(); ++row) {
    nct.putRowMain(row, main_rows_[row]);
  }

  return nct;
}

//----------------------------------------------------------------------------

void PlotMSCTAverager::initialize(ROCTIter& cti) {
  // Initialize the averager
  if (debug_) {
    cout << "  PMSCTA::initialize()" << endl;
  }

  // Assign main meta info
  avgField_ = cti.thisField();
  avgSpw_ = cti.thisSpw();
  avgObsid_ = cti.thisObs();
  timeRef_p = cti.thisTime();

  // Immutables:
  nChan_p = cti.nchan(); // TODO: channel averaging
  isComplex_p = cti.table().isComplex();

  if (averaging_p.baseline()) {
    nBlnMax_p = 1;
  } else {
    nBlnMax_p = nAnt_p;
  }

  if (debug_) {
    cout << "Shapes = " 
         << nPoln_p << " "
         << nChan_p << " "
         << nBlnMax_p << " "
         << ntrue(blnOK_p) << " "
         << endl;
  }

  blnOK_p.resize(nBlnMax_p);
  blnOK_p = false;

  // Resize and fill in the antenna numbers for all rows
  avgAntenna1_.resize(nBlnMax_p);
  avgAntenna2_.resize(nBlnMax_p);

  if (averaging_p.baseline()) {
    avgAntenna1_.set(-1);
    avgAntenna2_.set(-1);
  } else if (averaging_p.antenna()) {
    indgen(avgAntenna1_);
    avgAntenna2_.set(-1);
  } else {
    indgen(avgAntenna1_);
    avgAntenna2_.set(cti.thisAntenna2());
  }

  // Resize and fill if larger
  avgField_.resize(nBlnMax_p, true);
  avgSpw_.resize(nBlnMax_p, true);
  avgObsid_.resize(nBlnMax_p, true);
  Int nRows = cti.nrow();
  if (nBlnMax_p > nRows) {
    fillIds(nRows);
  }

  // Resize and initialize everything else
  avgTime_.resize(nBlnMax_p, false); 
  avgTime_.set(0.0);
  // Set scan, field for averaged chunk (set -1 later if more than one)
  aveScan_p = cti.thisScan();
  aveField_p = cti.thisField();
  // "Weight" is data count per baseline
  blnCount_p.resize(nBlnMax_p, false); 
  blnCount_p.set(0);
  // All cells assumed flagged to start with
  avgFlag_.resize(nPoln_p, nChan_p, nBlnMax_p, false);
  avgFlag_.set(true);
  // Accumulated data (PARAM column) and flags
  Complex czero(0.0);
  if (isComplex_p) {
    accumCParam_.resize(nPoln_p, nChan_p, nBlnMax_p, false);
    accumCParam_.set(czero);
  } else {
    accumFParam_.resize(nPoln_p, nChan_p, nBlnMax_p, false);
    accumFParam_.set(0.0);
  }

  minTimeOffset_p = DBL_MAX;
  maxTimeOffset_p = -DBL_MAX;
  initialized_p = true;
}

//----------------------------------------------------------------------------

void PlotMSCTAverager::simpleAccumulate (ROCTIter& cti)
{
// Accumulate a CTIter chunk
// Input:
//    cti              ROCTIter&            CalTable iter at current iteration
// Output to private data:
//    tStart_p         Double               Start time of current accumulation
//    avrow_p          Int                  Start row of current accumulation

  if (debug_) {
      cout << " PMSCTA::accumulate() " << endl;
  }

  // TODO
  if (averaging_p.channel()) {
    throw(AipsError("PlotMSCTAverager: channel averaging not implemented."));
  }
  if (averaging_p.baseline()) {
    throw(AipsError("PlotMSCTAverager: baseline averaging not implemented."));
  }
  if (averaging_p.spw()) {
    throw(AipsError("PlotMSCTAverager: spw averaging not implemented."));
  }

  // Only accumulate chunks with the same number of channels
  if (!averaging_p.channel() && (cti.nchan() != nChan_p)) {
    throw(AipsError("PlotMSCTAverager: data shape does not conform"));
  }

  if (!initialized_p) {
    initialize(cti);
  }

  // Param cubes for this iteration
  Cube<Complex> iterCParam;
  Cube<Float> iterFParam;
  if (isComplex_p) {
    iterCParam.reference(cti.cparam());
    if (averaging_p.scalarAve()) {
        convertToAP(iterCParam);
    }
  } else {
    iterFParam.reference(cti.fparam());
  }

  Int iterCount(0); // Sum of counts for this iteration
  for (Int ibln = 0; ibln < cti.nrow(); ++ibln) {
    // Calculate row from antenna numbers with the hash function.
    Int obln = cti.antenna1()(ibln);

    // This baseline occurs in input, so preserve in output
    blnOK_p(obln) = true;

    for (Int ichan = 0; ichan < cti.nchan(); ++ichan) {
      for (Int ipol = 0; ipol < nPoln_p; ++ipol) {
        // Assume we won't accumulate anything in this cell
        //   (output is unflagged, input is flagged)
        Bool accumulate(false);

        IPosition flagPos(3, ipol, ichan, ibln);
        if (!cti.flag()(flagPos)) { // input UNflagged
          accumulate = true;

          if (avgFlag_(ipol, ichan, obln)) {  // output flagged
            // This cell now NEWLY unflagged in output
            avgFlag_(ipol, ichan, obln) = false;

            // ...so zero the accumulators
            if (isComplex_p) {
              accumCParam_(ipol, ichan, obln) = 0.0;
            } else {
              accumFParam_(ipol, ichan, obln) = 0.0;
            }
            blnCount_p(obln) = 0;
          }
        } else { // input flagged
          // Only accumulate if output is also flagged
          //   (yields average of flagged data if no unflagged data ever found)
          if (avgFlag_(ipol, ichan, obln)) {
            accumulate = true;
          }
        }

        // Accumulate this (pol, chan) if appropriate
        if (accumulate) {
          if (isComplex_p) {
            accumCParam_(ipol, ichan, obln) += iterCParam(ipol, ichan, ibln);
          } else {
            accumFParam_(ipol, ichan, obln) += iterFParam(ipol, ichan, ibln);
          }

          blnCount_p(obln)++; // per-baseline count
          iterCount++;       // per-iteration count
        }
      } // pol 
    } // chan
  } // ibln

  if (iterCount > 0) {
    isAccum_p = true;

    // Set per-chunk averaged values
    // Used to determine averaged time when finalized
    Double timeOffset(cti.thisTime() - timeRef_p);
    minTimeOffset_p = min(minTimeOffset_p, timeOffset);
    maxTimeOffset_p = max(maxTimeOffset_p, timeOffset);
 
    Int thisScan = cti.thisScan();
    if (aveScan_p != thisScan) {
      aveScan_p = -1;
    }

    Int thisField = cti.thisField();
    if (aveField_p != thisField) {
      aveField_p = -1;
    }
  }
}

//----------------------------------------------------------------------------

void PlotMSCTAverager::antennaAccumulate (ROCTIter& cti) {
  // Accumulate a CTIter chunk with per-antenna averaging

  if (debug_) {
    cout << " PMSCTA::antAccumulate() " << endl;
  }

  if (averaging_p.antenna()) {
    throw(AipsError("PlotMSCTAverager: antenna averaging not implemented"));
  }

  // Only accumulate chunks with the same number of channels
  if (cti.nchan() != nChan_p) {
    throw(AipsError("PlotMSCTAverager: data shape does not conform"));
  }

  if (!initialized_p) {
    initialize(cti);
  }

  /*
  //  if (vb.spectralWindow()!=avBuf_p->spectralWindow())
  //    avBuf_p->spectralWindow()=-1;

  Double vbWt(0.0);  // will accumulate this CTs total data weight
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
    if (averaging_p.scalarAve()) convertToAP(accumVisCube);
  }
  if (doMVC_p) {
    accumVisCubeModel.reference(vb.visCubeModel());
    if (averaging_p.scalarAve()) convertToAP(accumVisCubeModel);
  }
  if (doCVC_p) {
    accumVisCubeCorrected.reference(vb.visCubeCorrected());
    if (averaging_p.scalarAve()) convertToAP(accumVisCubeCorrected);
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
      for (Int cor=0; cor<nPoln_p; ++cor) {
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
        if (avgFlag_(cor,chn,ia)) {  // output flagged
          // This cell now NEWLY unflagged in output
          avgFlag_(cor,chn,ia) = false;
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
      if (avgFlag_(cor,chn,obln_i)) acc_i = true;
      if (avgFlag_(cor,chn,obln_j)) acc_j = true;
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
    if (doUVW_p && blnWt>0.0) {
      for (uInt i=0; i<3; ++i) {
        avgUvw_(i,obln_i) += (vb.uvw()(i,ibln) * blnWt);
        avgUvw_(i,obln_j) += (vb.uvw()(i,ibln) * blnWt);
      }
    }
  }

  if (vbWt>0) {
    vbWtSum_p += vbWt;

    Double timeOffset(cti.thisTime() - timeRef_p);
    minTimeOffset_p = min(minTimeOffset_p, timeOffset);
    maxTimeOffset_p = max(maxTimeOffset_p, timeOffset);
 
    Int thisScan = vb.scan()(0);
    if (aveScan_p != thisScan)
      aveScan_p = -1;

    // This doesn't work...
    //    Int thisField=vb.fieldId();
    //    if (avBuf_p->fieldId()!=thisField)
    //      avBuf_p->fieldId()=-1;

  }
  */
}

//----------------------------------------------------------------------------

void PlotMSCTAverager::convertToAP(Cube<Complex>& data) {
  // Convert complex values to amp/phase
  Int ndata = data.nelements();
  Complex *c = data.data();

  for (Int i = 0; i < ndata; ++i, ++c) {
    Float a = abs(*c);
    Float p = arg(*c);
    *c = Complex(a, p);
  }
}

//----------------------------------------------------------------------------
void PlotMSCTAverager::fillIds(Int nrows) {
  // Fill with last value in id vector.
  // When you are resizing bigger than number of rows, fill rest of vector
  // or get seg fault from invalid value (whatever is in memory)
  for (Int i = nrows; i < nBlnMax_p; ++i) {
    field_[i] = field_[nrows-1];
    spw_[i] = spw_[nrows-1];
    obsid_[i] = obsid_[nrows-1];
  }
}

} //# NAMESPACE CASA - END

