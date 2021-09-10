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
    avgChan_p(false),
    nChanPerBin_p(1),
    nAvgChan_p(0),
    isComplex_p(True),
    isAntennaBased_p(True),
    blnOK_p(),
    timeRef_p(0.0),
    minTimeOffset_p(0.0),
    maxTimeOffset_p(0.0),
    avgScan_p(0),
    avgField_p(0),
    avgSpw_p(0),
    initialized_p(false),
    isAccum_p(false),
    debug_(false) {
}

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
    Int numchan(nchan());

    casacore::Cube<casacore::Complex> cparam(nPoln_p, numchan, nBlnOK);
    cparam.set(Complex(0.0));
    casacore::Cube<casacore::Float> fparam(nPoln_p, numchan, nBlnOK);
    fparam.set(0.0);
    casacore::Cube<casacore::Float> paramerr(nPoln_p, numchan, nBlnOK);
    paramerr.set(0.0);
    casacore::Cube<casacore::Bool> flag(nPoln_p, numchan, nBlnOK);
    flag.set(true);
    casacore::Cube<casacore::Float> snr(nPoln_p, numchan, nBlnOK);
    snr.set(0.0);
    casacore::Vector<casacore::Int> antenna1(nBlnOK);
    casacore::Vector<casacore::Int> antenna2(nBlnOK);
    if (avgChan_p) {
      avgFreq_.resize(nAvgChan_p);
    }

    // Divide each baseline sum by its weight/count
    for (Int ibln = 0; ibln < nBlnMax_p; ++ibln) {
      if (blnOK_p(ibln)) {
        // baseline was in unaveraged data
        for (Int ichan = 0; ichan < numchan; ++ichan) {
          for (Int ipol = 0; ipol < nPoln_p; ++ipol) {
            // normalize ibln data at obln:
            Float thisWt(accumWt_(ipol, ichan, ibln));
            if (thisWt > 0.0) {
              if (isComplex_p) {
                cparam(ipol, ichan, obln) = accumCParam_(ipol, ichan, ibln) / thisWt;
              } else {
                fparam(ipol, ichan, obln) = accumFParam_(ipol, ichan, ibln) / thisWt;
              }

              paramerr(ipol, ichan, obln) = accumParamErr_(ipol, ichan, ibln) / thisWt;
              snr(ipol, ichan, obln) = accumSnr_(ipol, ichan, ibln) / thisWt;

              if (avgChan_p) {
                avgFreq_(ichan) = accumFreq_(ichan) / nChanPerBin_p;
              }
            }

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
      avgSpw_.resize(nBlnOK);
      obsid_.resize(nBlnOK, true);
    }

    // Time: center of the interval
    avgTime_.set(timeRef_p + (maxTimeOffset_p + minTimeOffset_p)/2.0);
    // Value for all iterations, or -1 if they are combined
    avgScan_.set(avgScan_p);
    avgField_.set(avgField_p);
    avgSpw_.set(avgSpw_p);

    // Assign averaged data
    for (Int ibln = 0; ibln < nBlnOK; ++ibln) {
      CTMainRecord main_row_record;
      main_row_record.defineTime(avgTime_(ibln));
      main_row_record.defineFieldId(avgField_(ibln));
      main_row_record.defineSpwId(avgSpw_(ibln));
      main_row_record.defineAntenna1(antenna1(ibln));
      main_row_record.defineAntenna2(antenna2(ibln));
      main_row_record.defineInterval(0); // axis not supported
      main_row_record.defineScanNo(avgScan_(ibln));
      main_row_record.defineObsId(obsid_(ibln));

      if (isComplex_p) {
        main_row_record.defineCParam(cparam[ibln]);
      } else {
        main_row_record.defineFParam(fparam[ibln]);
      }
      main_row_record.defineParamerr(paramerr[ibln]);
      main_row_record.defineFlag(flag[ibln]);
      main_row_record.defineSnr(snr[ibln]);
      main_row_record.defineWeight(accumWt_[ibln]);

      main_rows_.push_back(main_row_record);
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

void PlotMSCTAverager::fillAvgCalTable(NewCalTable& caltab) {
  // Fill NewCalTable with averaging results stored in CTMainRecord vector
  // as rows in main table.
  // Cal table passed in should have subtables added.
  for (size_t row = 0; row < main_rows_.size(); ++row) {
    caltab.putRowMain(row, main_rows_[row]);
  }
}

//----------------------------------------------------------------------------

void PlotMSCTAverager::initialize(ROCTIter& cti) {
  // Initialize the averager
  if (debug_) {
    cout << "  PMSCTA::initialize()" << endl;
  }

  // Assign main meta info: some may be reset to -1 if multiple ids used for average
  avgField_p = cti.thisField();
  avgScan_p = cti.thisScan();
  avgSpw_p = cti.thisSpw();
  timeRef_p = cti.thisTime();
  Int obsid = cti.thisObs();

  if (averaging_p.antenna() && (cti.thisAntenna2() == -1)) {
    throw(AipsError("Antenna averaging invalid for pure antenna-based table"));
  }

  // Immutables:
  isComplex_p = cti.table().isComplex();

  // Set number of baselines
  if (averaging_p.baseline()) {
    nBlnMax_p = 1;
  } else if ((isAntennaBased_p) || averaging_p.antenna()) {
    nBlnMax_p = nAnt_p;
  } else {
    nBlnMax_p = nAnt_p * (nAnt_p + 1)/2;
  }

  // Set number of channels
  nChan_p = cti.nchan();
  // Set num chan and num chan per bin if averaging channels
  if (averaging_p.channel()) {
    avgChan_p = true;
    nChanPerBin_p = averaging_p.channelValue();

    if (nChanPerBin_p <= 1) {
        avgChan_p = false; // averaging 0 or 1 chan is no averaging
    }

    if (avgChan_p) {
      if (nChanPerBin_p > nChan_p) {
        nChanPerBin_p = nChan_p; // average all channels
      }

      nAvgChan_p = nChan_p / nChanPerBin_p;
    }
  }


  if (debug_) {
    cout << "Shapes = " 
         << nPoln_p << " "
         << nchan() << " "
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
    if (isAntennaBased_p) {
      indgen(avgAntenna1_);
      avgAntenna2_.set(cti.thisAntenna2());
    } else {
      // Form all combinations of baselines
      Int ibln = 0;
      for (Int iant1 = 0; iant1 < nAnt_p; ++iant1) {
        for (Int iant2 = iant1; iant2 < nAnt_p; ++iant2) {
	      avgAntenna1_(ibln) = iant1;
	      avgAntenna2_(ibln) = iant2;
	      ++ibln;
        }
      }
    } 
  }

  if (avgChan_p) {
    accumFreq_.resize(nAvgChan_p);
    accumFreq_.set(0.0);
  } else {
    avgFreq_ = cti.freq();
  }

  // Resize to number of baselines
  avgScan_.resize(nBlnMax_p);
  avgField_.resize(nBlnMax_p);
  avgSpw_.resize(nBlnMax_p);
  obsid_.resize(nBlnMax_p);
  obsid_.set(obsid);

  // Resize and initialize everything else
  avgTime_.resize(nBlnMax_p, false); 
  avgTime_.set(0.0);

  // Resize accumulated data (PARAM column) and flags
  Int numchan(nchan()); // depends on channel averaging
  avgFlag_.resize(nPoln_p, numchan, nBlnMax_p, false);
  avgFlag_.set(true); // All cells assumed flagged to start with
 
  Complex czero(0.0);
  if (isComplex_p) {
    accumCParam_.resize(nPoln_p, numchan, nBlnMax_p, false);
    accumCParam_.set(czero);
  } else {
    accumFParam_.resize(nPoln_p, numchan, nBlnMax_p, false);
    accumFParam_.set(0.0);
  }
  accumParamErr_.resize(nPoln_p, numchan, nBlnMax_p, false);
  accumParamErr_.set(0.0);
  accumSnr_.resize(nPoln_p, numchan, nBlnMax_p, false);
  accumSnr_.set(0.0);
  accumWt_.resize(nPoln_p, numchan, nBlnMax_p, false);
  accumWt_.set(0.0);

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

  if (!initialized_p) {
    initialize(cti);
  }

  // Only accumulate chunks with the same number of channels
  if (cti.nchan() != nChan_p) {
    throw(AipsError("PlotMSCTAverager: data shape does not conform"));
  }

  Cube<Complex> iterCParam;
  Cube<Float> iterFParam, iterParamErr, iterSnr, iterWt;
  Cube<Bool> iterFlag;
  Vector<Double> iterFreq;

  if (isComplex_p) {
    iterCParam.reference(cti.cparam());
    if (averaging_p.scalarAve()) {
        convertToAP(iterCParam);
    }
  } else {
    iterFParam.reference(cti.fparam());
  }
  iterParamErr.reference(cti.paramErr());
  iterSnr.reference(cti.snr());
  iterFlag.reference(cti.flag());
  iterFreq.reference(cti.freq());

  try {
    iterWt.reference(cti.wt());
  } catch (const AipsError& err) {
    casacore::IPosition datashape =
      (isComplex_p ? iterCParam.shape() : iterFParam.shape());
    iterWt.resize(datashape);
    iterWt.set(1.0);
  }

  Int accumCount(0);

  for (Int ibln = 0; ibln < cti.nrow(); ++ibln) {
    // Calculate output row from antenna1 numbers
    Int ant1 = cti.antenna1()(ibln);
    Int ant2 = cti.antenna2()(ibln);
    Int obln = baseline_index(ant1, ant2);

    // This baseline occurs in input, so preserve in output
    blnOK_p(obln) = true;

    for (Int ipol = 0; ipol < nPoln_p; ++ipol) {
      Int outchan(-1);
      for (Int ichan = 0; ichan < cti.nchan(); ++ichan) {
        // Assume we won't accumulate anything in this cell
        //   (output is unflagged, input is flagged)
        Bool accumulate(false);

        // Set outchan
        if (avgChan_p) {
          if ((ichan == 0) || ((ichan % nChanPerBin_p) == 0)) {
            // next outchan
            ++outchan;

            // Ignore remaining channels, incomplete bin
            if (outchan == nAvgChan_p) {
              break;
            }
          }
        } else {
          outchan = ichan; // use same channel when not averaging
        }

        IPosition inPos(3, ipol, ichan, ibln);
        IPosition outPos(3, ipol, outchan, obln);

        if (!iterFlag(inPos)) { // input UNflagged
          accumulate = true;

          if (avgFlag_(outPos)) {  // output flagged
            // This cell now NEWLY unflagged in output
            avgFlag_(outPos) = false;

            // ...so zero the accumulators
            if (isComplex_p) {
              accumCParam_(outPos) = 0.0;
            } else {
              accumFParam_(outPos) = 0.0;
            }
            accumParamErr_(outPos) = 0.0;
            accumSnr_(outPos) = 0.0;
            accumWt_(outPos) = 0.0;
            if (avgChan_p) {
              accumFreq_(outchan) = 0.0;
            }
          }
        } else { // input flagged
          // Only accumulate if output is also flagged
          //   (yields average of flagged data if no unflagged data ever found)
          if (avgFlag_(outPos)) {
            accumulate = true;
          }
        }

        // Accumulate this (pol, chan) if appropriate
        if (accumulate) {
          accumCount++;

          Float wt = iterWt(inPos);
          if (wt < FLT_MIN) {
            wt = FLT_MIN;
          }
          accumWt_(outPos) += wt;

          if (isComplex_p) {
            accumCParam_(outPos) += wt * iterCParam(inPos);
          } else {
            accumFParam_(outPos) += wt * iterFParam(inPos);
          }
          accumParamErr_(outPos) += wt * iterParamErr(inPos);
          accumSnr_(outPos) += wt * iterSnr(inPos);

          if (avgChan_p) {
            accumFreq_(outchan) += iterFreq(ichan);
          }
        }
      } // pol 
    } // chan
  } // ibln

  if (accumCount > 0) {
    isAccum_p = true;

    // Set per-chunk averaged values
    // Used to determine averaged time when finalized
    Double timeOffset(cti.thisTime() - timeRef_p);
    minTimeOffset_p = min(minTimeOffset_p, timeOffset);
    maxTimeOffset_p = max(maxTimeOffset_p, timeOffset);

    // Set multiple scans, fields, spws to -1 if combined
    if (avgScan_p != cti.thisScan()) {
      avgScan_p = -1;
    }
    if (avgField_p != cti.thisField()) {
      avgField_p = -1;
    }
    if (avgSpw_p != cti.thisSpw()) {
      avgSpw_p = -1;
    }
  }
}

//----------------------------------------------------------------------------

void PlotMSCTAverager::antennaAccumulate (ROCTIter& cti) {
  // Accumulate a CTIter chunk with per-antenna averaging

  if (debug_) {
    cout << " PMSCTA::antAccumulate() " << endl;
  }


  if (!initialized_p) {
    initialize(cti);
  }

  // Only accumulate chunks with the same number of channels
  if (cti.nchan() != nChan_p) {
    throw(AipsError("PlotMSCTAverager: data shape does not conform"));
  }

  // Param cubes for this iteration
  Cube<Complex> iterCParam;
  Cube<Float> iterFParam, iterParamErr, iterSnr, iterWt;
  Cube<Bool> iterFlag;
  Vector<Double> iterFreq;

  if (isComplex_p) {
    iterCParam.reference(cti.cparam());
    if (averaging_p.scalarAve()) {
        convertToAP(iterCParam);
    }
  } else {
    iterFParam.reference(cti.fparam());
  }
  iterParamErr.reference(cti.paramErr());
  iterSnr.reference(cti.snr());
  iterFlag.reference(cti.flag());
  iterFreq.reference(cti.freq());

  try {
    iterWt.reference(cti.wt());
  } catch (const AipsError& err) {
    casacore::IPosition datashape =
      (isComplex_p ? iterCParam.shape() : iterFParam.shape());
    iterWt.resize(datashape);
    iterWt.set(1.0);
  }

  Int accumCount(0);
  for (Int ibln = 0; ibln < cti.nrow(); ++ibln) {
    // The antennas in the baseline
    Vector<Int> oblnij(2);
    oblnij(0) = cti.antenna1()(ibln);
    oblnij(1) = cti.antenna2()(ibln);
    Int& obln_i(oblnij(0));
    Int& obln_j(oblnij(1));

    // These antennas occur (even if flagged) in input, so preserve in output
    blnOK_p(obln_i) = true;
    blnOK_p(obln_j) = true;
 
    for (Int ipol = 0; ipol < nPoln_p; ++ipol) {
      Int outchan(-1);
      for (Int ichan = 0; ichan < cti.nchan(); ichan++) {
        // Assume we won't accumulate anything in this cell
        //   (output is unflagged, input is flagged)
        Bool accum_i(false), accum_j(false);

        // Set outchan
        if (avgChan_p) {
          if ((ichan == 0) || ((ichan % nChanPerBin_p) == 0)) {
            // next outchan
            ++outchan;

            // Ignore remaining channels, incomplete bin
            if (outchan == nAvgChan_p) {
              break;
            }
          }
        } else {
          outchan = ichan; // use same channel when not averaging
        }

        // Consider accumulation according to state of flags
        IPosition inPos(3, ipol, ichan, ibln);
        if (!iterFlag(inPos)) { // input UNflagged
          // we will accumulate both ants
          accum_i = accum_j = true;

          // Zero accumulators if output cell currently flagged
          for (Int ij = 0; ij < 2; ++ij) {
            Int ia = oblnij(ij);
            IPosition iaPos(3, ipol, outchan, ia);
            if (avgFlag_(iaPos)) { // output flagged
              // This cell now NEWLY unflagged in output
              avgFlag_(iaPos) = false;
              // ...so zero the accumulators
              if (isComplex_p) {
                accumCParam_(iaPos) = 0.0;
              } else {
                accumFParam_(iaPos) = 0.0;
              }
              accumParamErr_(iaPos) = 0.0;
              accumSnr_(iaPos) = 0.0;
              accumWt_(iaPos) = 0.0;

              if (avgChan_p) {
                accumFreq_(iaPos) = 0.0;
              }
            }
          }
        } else { // input cell is flagged
          // Only accumulate if output is also flagged:
          //  (yields average of flagged data if no unflagged data ever found)
          if (avgFlag_(ipol, ichan, obln_i)) {
            accum_i = true;
          }
          if (avgFlag_(ipol, ichan, obln_j)) {
            accum_j = true;
          }
        }

        Float wt = iterWt(inPos);
        if (wt < FLT_MIN) {
          wt = FLT_MIN;
        }

        // Accumulate data, if appropriate
        if (accum_i) {
          IPosition ipos(3, ipol, outchan, obln_i);
          if (isComplex_p) {
            accumCParam_(ipos) += wt * iterCParam(inPos);
          } else {
            accumFParam_(ipos) += wt * iterFParam(inPos);
          }
          accumParamErr_(ipos) += wt * iterParamErr(inPos);
          accumSnr_(ipos) += wt * iterSnr(inPos);
          accumWt_(ipos) += wt;
        }

        if (accum_j) { 
          IPosition jpos(3, ipol, outchan, obln_j);
          if (isComplex_p) {
            accumCParam_(jpos) += wt * iterCParam(inPos);
          } else {
            accumFParam_(jpos) += wt * iterFParam(inPos);
          }
          accumParamErr_(jpos) += wt * iterParamErr(inPos);
          accumSnr_(jpos) += wt * iterSnr(inPos);
          accumWt_(jpos) += wt;
        }

        if (accum_i || accum_j) {
          accumCount++;
          if (avgChan_p) {
            accumFreq_(outchan) += iterFreq(ichan);
          }
        }
      } // chan 
    } // pol
  } // bln

  if (accumCount > 0) {
    isAccum_p = true;

    // average time of accumulated iterations
    Double timeOffset(cti.thisTime() - timeRef_p);
    minTimeOffset_p = min(minTimeOffset_p, timeOffset);
    maxTimeOffset_p = max(maxTimeOffset_p, timeOffset);
 
    if (avgScan_p != cti.thisScan()) {
      avgScan_p = -1;
    }
    if (avgField_p != cti.thisField()) {
      avgField_p = -1;
    }
    if (avgSpw_p != cti.thisSpw()) {
      avgSpw_p = -1;
    }
  }
  casacore::IPosition checkpos(3, 0,0,0);
}

//----------------------------------------------------------------------------

Int PlotMSCTAverager::baseline_index(const Int& ant1, const Int& ant2)
{
  // Compute row index in an accumulation interval for given ant1, ant2
  Int index;
  if (averaging_p.baseline()) {
    index = 0;
  } else if (isAntennaBased_p) {
    index = ant1;
  } else {
    index = nAnt_p * ant1 - (ant1 * (ant1 - 1)) / 2 + ant2 - ant1;
  }

  return index;
};

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

} //# NAMESPACE CASA - END

