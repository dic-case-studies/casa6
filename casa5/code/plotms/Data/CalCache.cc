//# CalCache.cc: Specialized PlotMSCache for filling CalTables
//# Copyright (C) 2009
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
//# $Id: $
#include <plotms/Data/CalCache.h>
#include <plotms/Data/PlotMSIndexer.h>
#include <plotms/Data/PlotMSAtm.h>
#include <plotms/PlotMS/PlotMS.h>
#include <plotms/Threads/ThreadCommunication.h>

#include <casa/OS/Timer.h>
#include <casa/OS/HostInfo.h>
#include <casa/OS/Memory.h>
#include <casa/Quanta/MVTime.h>
#include <casa/System/Aipsrc.h>
#include <casa/Utilities/Sort.h>
#include <casa/Arrays/ArrayMath.h>
#include <tables/Tables/Table.h>
#include <lattices/Lattices/ArrayLattice.h>
#include <synthesis/CalTables/CTColumns.h>
#include <synthesis/MeasurementComponents/VisCalGlobals.h>
#include <synthesis/MeasurementComponents/BPoly.h>
#include <synthesis/MeasurementComponents/GSpline.h>

using namespace casacore;

namespace casa {

// Define external CLIC solvers
#define cheb cheb_

extern "C" {
  void cheb(Int*, Double*, Double*, Int*);
}

CalCache::CalCache(PlotMSApp* parent):
  PlotMSCacheBase(parent),
  divZero_(False),
  ci_p(nullptr),
  wci_p(nullptr),
  basis_("unknown"),
  parsAreComplex_(False),
  msname_("")
{
}

CalCache::~CalCache() {}


String CalCache::polname(Int ipol) {
  if (polnRatio_) return "/";
  if (basis_=="Linear")
    return ( (ipol%2==0) ? String("X") : String("Y") );
  else if (basis_=="Circular")
    return ( (ipol%2==0) ? String("R") : String("L") );
  else { // "unknown", or antenna positions
    if (calType_=="KAntPos Jones") {
        switch(ipol) {
            case 0: return "X";
            case 1: return "Y";
            case 2: return "Z";
            default: return (String::toString(ipol));
        }
    } else {
        return ( String::toString(ipol) );
    }
  }
}

void CalCache::setFilename(String filename) { 
    filename_ = filename;
    Table tab(filename);
    calType_= tab.tableInfo().subType();

    if ((calType_=="T Jones") && (tab.keywordSet().isDefined("CAL_DESC"))) {
      throw AipsError(calType_ + " tables in the old cal table format are unsupported in plotms.");
    }
}

//*********************************
// protected method implementations

void CalCache::loadIt(vector<PMS::Axis>& loadAxes,
    vector<PMS::DataColumn>& /*loadData*/,
    ThreadCommunication* thread) {

  // this also sets calType_:
  setFilename(filename_);

  // Trap unsupported modes: averaging, transforms, poln ratio
  logLoad("Plotting a " + calType_ + " calibration table.");
  // Warn that averaging and transformations will be ignored
  if (averaging().anyAveraging())
    logWarn("CalCache::loadIt",
      "Averaging ignored: not supported for calibration tables");
  if (transformations().anyTransform())
    logWarn("CalCache::loadIt",
      "Transformations ignored: not supported for calibration tables");
  // poln ratio
  polnRatio_ = false;
  if (selection_.corr()=="/") {
    if ((calType_=="BPOLY") || (calType_[0] == 'T') ||
        ((calType_[0] == 'F') && !calType_.startsWith("Fringe"))) {
      throw(AipsError("Polarization ratio plots not supported for " + calType_ + " tables."));
    } else {
      polnRatio_ = true;
    }
  }

  antnames_.resize();
  stanames_.resize();
  antstanames_.resize();
  fldnames_.resize();
  positions_.resize();

  vector<PMS::DataColumn> loadData(loadAxes.size());
  for (uInt i=0; i<loadData.size(); ++i) { 
    loadData[i] = PMS::DEFAULT_DATACOLUMN;
  }

  if (calType_=="BPOLY") {
    loadBPoly(loadAxes, loadData, thread);
  } else if (calType_=="GSPLINE") {
    checkAxes(loadAxes);  // check for invalid axis before proceeding
    loadGSpline(loadAxes, loadData, thread);
  } else {
    loadNewCalTable(loadAxes, loadData, thread);
  }
}

// ======================== NewCalTable ==========================

void CalCache::loadNewCalTable(vector<PMS::Axis>& loadAxes,
    vector<PMS::DataColumn>& loadData, ThreadCommunication* thread) {
  // Get various names, properties from cal table
  TableLock lock(TableLock::AutoNoReadLocking);
  NewCalTable* ct = new NewCalTable(filename_, lock, Table::Old, Table::Plain);
  basis_ = ct->polBasis();
  parsAreComplex_ = ct->isComplex();
  ROCTColumns ctCol(*ct);
  antnames_ = ctCol.antenna().name().getColumn();
  stanames_ = ctCol.antenna().station().getColumn();
  antstanames_ = antnames_ + String("@") + stanames_;
  fldnames_ = ctCol.field().name().getColumn();
  positions_ = ctCol.antenna().position().getColumn();    
  nAnt_ = ctCol.antenna().nrow();

  // Apply selection to get selected cal table
  NewCalTable* selct = new NewCalTable();
  Vector<Vector<Slice> > chansel;
  Vector<Vector<Slice> > corrsel;
  selection_.apply(*ct, *selct, chansel, corrsel);

  Bool readonly(True); // no write access for loading cache
  setUpCalIter(*selct, readonly);
  countChunks(*ci_p, loadAxes, loadData, thread);
  loadCalChunks(*ci_p, loadAxes, thread);

  // delete NCT and iter to release table locks
  if (ct != nullptr) {
    delete ct;
    ct = nullptr;
  }
  if (selct != nullptr) {
    delete selct;
    selct = nullptr;
  }
  if (ci_p != nullptr) {
    delete ci_p;
    ci_p = nullptr;
  }
}

void CalCache::setUpCalIter(NewCalTable& selct, Bool readonly) {
  Int nsortcol(4);
  Block<String> columns(nsortcol);
  columns[0]="SCAN_NUMBER";
  columns[1]="FIELD_ID";
  columns[2]="SPECTRAL_WINDOW_ID";
  columns[3]="TIME";

  if (readonly) {
    // Readonly version, for caching
    ci_p = new ROCTIter(selct, columns);
    wci_p = nullptr;
  } else {
    // Writable, e.g. for flagging
    wci_p = new CTIter(selct, columns);
    ci_p = wci_p;  // const access
  }
}
      
void CalCache::countChunks(ROCTIter& ci,
    vector<PMS::Axis>& loadAxes,
    vector<PMS::DataColumn>& loadData,
    ThreadCommunication* thread) {
  // for NewCalTable
  if (thread!=nullptr) {
    thread->setStatus("Establishing cache size.  Please wait...");
    thread->setAllowedOperations(false,false,false);
  }

  // Count number of chunks.
  int chunk(0);
  ci.reset();
  while (!ci.pastEnd()) {
    ++chunk;
    ci.next0();
  }
  setCache(chunk, loadAxes, loadData);
}

void CalCache::loadCalChunks(ROCTIter& ci,
     const vector<PMS::Axis> loadAxes,
     ThreadCommunication* thread) {
  // for NewCalTable
  // permit cancel in progress meter:
  if(thread != nullptr)
    thread->setAllowedOperations(false,false,true);
  logLoad("Loading chunks......");

  Int chunk(0), lastscan(0), thisscan(0), lastspw(-1), thisspw(0);
  chshapes_.resize(4,nChunk_);
  goodChunk_.resize(nChunk_);
  goodChunk_.set(False);
  double progress;

  // Reset iterator
  ci.reset();
  while (!ci.pastEnd()) {
      // If a thread is given, check if the user canceled.
      if(thread != nullptr && thread->wasCanceled()) {
        dataLoaded_ = false;
        return;
      }

      // If a thread is given, update it.
      if(thread != nullptr && (nChunk_ <= (int)THREAD_SEGMENT ||
         chunk % THREAD_SEGMENT == 0)) {
          thread->setStatus("Loading chunk " + String::toString(chunk) +
              " / " + String::toString(nChunk_) + ".");
      }
      
      // Discern npar/nchan shape
      IPosition pshape(ci.flag().shape());

      // Determine data axis for determining param slice for npol
      String pol = selection_.corr();
      String paramAxis = toVisCalAxis(PMS::AMP);

      size_t nPol;
      if (polnRatio_) {  // length is for 1 poln, pick first one
          nPol = getParSlice(paramAxis, "R").length();
      } else {
          nPol = getParSlice(paramAxis, pol).length();
      }

      // Cache the data shapes
      chshapes_(0,chunk) = nPol;
      chshapes_(1,chunk) = pshape[1];
      chshapes_(2,chunk) = ci.nrow();
      chshapes_(3,chunk) = nAnt_;
      goodChunk_(chunk) = True;

      for(unsigned int i = 0; i < loadAxes.size(); i++) {
        loadCalAxis(ci, chunk, loadAxes[i], pol);
        // print atm stats once per scan
        if (loadAxes[i]==PMS::ATM || loadAxes[i]==PMS::TSKY) {
            thisscan = ci.thisScan();
            if (thisscan != lastscan) {
                printAtmStats(thisscan);
                lastscan = thisscan;
            }
            thisspw = ci.thisSpw();
            if (thisspw != lastspw) {
                uInt vectorsize = ( loadAxes[i]==PMS::ATM ?
                    (*atm_[chunk]).nelements() :
                    (*tsky_[chunk]).nelements());
                if (vectorsize==1) {
                    logWarn("load_cache", "Setting " + 
                        PMS::axis(loadAxes[i]) + " for spw " +
                        String::toString(thisspw) +
                        " to zero because it has only one channel.");
                }
                lastspw = thisspw;
            }
        }
      }

      chunk++;
      ci.next();

      // If a thread is given, update it.
      if (thread != nullptr && (nChunk_ <= (int)THREAD_SEGMENT ||
          chunk % THREAD_SEGMENT == 0)) {
          progress = ((double)chunk+1) / nChunk_;
          thread->setProgress((unsigned int)((progress * 100) + 0.5));
      }
  }

  if (divZero_) {
    logWarn("CalCache::loadIt", "Caught divide-by-zero exception in ratio plots; result(s) set to 1.0 and flagged");
  }
}

void CalCache::loadCalAxis(ROCTIter& cti, Int chunk, PMS::Axis axis, String pol) {
    // Load axis from NewCalTable
    Slice parSlice1 = Slice();
    Slice parSlice2 = Slice();
    if (PMS::axisNeedsCalSlice(axis)) {
        String calAxis = toVisCalAxis(axis);
        if (polnRatio_) {
            parSlice1 = getParSlice(calAxis, "R");
            parSlice2 = getParSlice(calAxis, "L");
        } else {
            parSlice1 = getParSlice(calAxis, pol);
        }
    }

    switch(axis) {
        case PMS::SCAN: // assumes scan unique
            scan_(chunk) = cti.thisScan();
            break;
        case PMS::FIELD:
            field_(chunk) = cti.thisField();
            break;
        case PMS::TIME: // assumes time unique 
            time_(chunk) = cti.thisTime();
            break;
        /*        
        case PMS::TIME_INTERVAL: // assumes timeInterval unique in VB
            timeIntr_(chunk) = cti.interval()(0); 
            break;
        */
        case PMS::SPW:
            spw_(chunk) = cti.thisSpw();
            break;
        case PMS::CHANNEL: 
            cti.chan(*chan_[chunk]);
            break;
        case PMS::FREQUENCY: {
            // TBD: Convert freq to desired frame
            cti.freq(*freq_[chunk]);
            (*freq_[chunk])/=1.0e9; // in GHz
            break;
        }
        /*
        case PMS::VELOCITY: {
            // Convert freq in the vb to velocity
            vbu_.toVelocity(*vel_[chunk], vb, transformations_.frame(),
            MVFrequency(transformations_.restFreqHz()),
            transformations_.veldef());
            (*vel_[chunk]) /= 1.0e3;  // in km/s
            break;
        }
        */
        case PMS::CORR: {
            corr_[chunk]->resize(chshapes_(0,chunk));
            if (pol=="" || pol=="RL" || pol=="XY") {
                indgen(*corr_[chunk]);
            } else if (pol== "R" || pol=="X") { 
                corr_[chunk]->resize(1);
                corr_[chunk]->set(0);
            } else if (pol== "L" || pol=="Y") { 
                corr_[chunk]->resize(1);
                corr_[chunk]->set(1);
            } else if (pol=="/") {
                corr_[chunk]->resize(1);
                corr_[chunk]->set(-1); // ???
            }
            break;
        }
        case PMS::ANTENNA1:
            *antenna1_[chunk] = cti.antenna1(); 
            break;
        case PMS::ANTENNA2:
            *antenna2_[chunk] = cti.antenna2(); 
            break;
        case PMS::BASELINE: {
            Vector<Int> a1(cti.antenna1());
            Vector<Int> a2(cti.antenna2());
            baseline_[chunk]->resize(cti.nrow());
            Vector<Int> bl(*baseline_[chunk]);
            for (Int irow=0;irow<cti.nrow();++irow) {
                if (a1(irow)<0) a1(irow)=chshapes_(3,0);
                if (a2(irow)<0) a2(irow)=chshapes_(3,0);
                bl(irow) = (chshapes_(3,0)+1)*a1(irow) -
                    (a1(irow)*(a1(irow) + 1))/2 + a2(irow);
            }
            break;
        }
        case PMS::ANTPOS: {
            if (!calType_.startsWith("KAntPos")) {
                throw(AipsError( "ANTPOS has no meaning for this table"));
            }

            Cube<Float> fArray = cti.fparam();
            *antpos_[chunk] = fArray(parSlice1, Slice(), Slice());
            break;
        }
        case PMS::GAMP:
        case PMS::AMP: {
            if (parsAreComplex()) {
                Cube<Complex> cArray = cti.cparam();
                if (polnRatio_) {
                    Array<Float> ampRatio = amplitude(cArray(parSlice1, Slice(),
                        Slice()) / cArray(parSlice2, Slice(), Slice()));
                    checkRatioArray(ampRatio, chunk);
                    *amp_[chunk] = ampRatio;
                } else {
                    *amp_[chunk] = amplitude(cArray(parSlice1, Slice(), Slice()));
                }
            } else {
                Cube<Float> fArray = cti.fparam();
                if (polnRatio_) {
                    if (calType_ == "Fringe Jones") { // subtract
                        *amp_[chunk] = fArray(parSlice1, Slice(), Slice()) -
                            fArray(parSlice2, Slice(), Slice());
                    } else {
                        Array<Float> ampRatio = fArray(parSlice1, Slice(), Slice()) /
                            fArray(parSlice2, Slice(), Slice());
                        checkRatioArray(ampRatio, chunk);
                        *amp_[chunk] = ampRatio;
                    }
                } else {        
                    *amp_[chunk] = fArray(parSlice1, Slice(), Slice());
                }
                if (calType_ == "F Jones") { // TEC table
                    (*amp_[chunk]) /= Float(1e+16);
                }
            }
            break;
        }
        case PMS::GPHASE:
        case PMS::PHASE: {
            if (parsAreComplex()) {
                Cube<Complex> cArray = cti.cparam();
                if (polnRatio_) {
                    Array<Float> phaseRatio = phase(cArray(parSlice1, Slice(),
                        Slice()) / cArray(parSlice2, Slice(), Slice()));
                    checkRatioArray(phaseRatio, chunk);
                    *pha_[chunk] = phaseRatio;
                } else {
                    *pha_[chunk] = phase(cArray(parSlice1, Slice(), Slice()));
                }
                (*pha_[chunk]) *= Float(180.0/C::pi);
            } else if (calType_ == "Fringe Jones") {
                Cube<Float> fArray = cti.fparam();
                if (polnRatio_) {
                    *pha_[chunk] = fArray(parSlice1, Slice(), Slice()) -
                        fArray(parSlice2, Slice(), Slice());
                } else {
                    *pha_[chunk] = fArray(parSlice1, Slice(), Slice());
                }
                (*pha_[chunk]) *= Float(180.0/C::pi);
            } else {
                throw(AipsError("phase has no meaning for this table"));
            }
            break;
        }
        case PMS::GREAL:   
        case PMS::REAL: {
            if (calType_ == "Fringe Jones") { // do not use float for this axis 
                throw(AipsError("real has no meaning for this table"));
            } else if (parsAreComplex()) {
                Cube<Complex> cArray = cti.cparam();
                if (polnRatio_) {
                    Array<Float> realRatio = real(cArray(parSlice1, Slice(),
                        Slice()) / cArray(parSlice2, Slice(), Slice()));
                    checkRatioArray(realRatio, chunk);
                    *real_[chunk] = realRatio;
                } else {        
                    *real_[chunk] = real(cArray(parSlice1, Slice(), Slice()));
                }
            } else {  // allow float for single dish cal tables
                Cube<Float> fArray = cti.fparam();
                if (polnRatio_) {
                    Array<Float> ampRatio = fArray(parSlice1, Slice(), Slice()) / fArray(parSlice2, Slice(), Slice());
                    checkRatioArray(ampRatio, chunk);
                    *real_[chunk] = ampRatio;
                } else {        
                    *real_[chunk] = fArray(parSlice1, Slice(), Slice());
                }
            }
            break;
        }
        case PMS::GIMAG:
        case PMS::IMAG: {
            if (parsAreComplex()) {
                Cube<Complex> cArray = cti.cparam();
                if (polnRatio_) {
                    Array<Float> imagRatio = imag(cArray(parSlice1, Slice(),
                        Slice()) / cArray(parSlice2, Slice(), Slice()));
                    checkRatioArray(imagRatio, chunk);
                    *imag_[chunk] = imagRatio;
                } else {        
                    *imag_[chunk] = imag(cArray(parSlice1, Slice(), Slice()));
                }
            } else {
                throw(AipsError("imag has no meaning for this table"));
            }
            break;
        }
        case PMS::DELAY:{
            if (!parsAreComplex()) {
                Cube<Float> fArray = cti.fparam();
                if (polnRatio_) {
                    Array<Float> delayRatio = fArray(parSlice1, Slice(), Slice())
                        - fArray(parSlice2, Slice(), Slice());
                    *par_[chunk] = delayRatio;
                } else {
                    *par_[chunk] = fArray(parSlice1, Slice(), Slice());
                }
            } else {
                throw(AipsError("delay has no meaning for this table"));
            }
            break;
        }
        case PMS::DELAY_RATE: {
            if (calType_.startsWith("Fringe") && !parsAreComplex()) {
                Cube<Float> fArray = cti.fparam();
                if (polnRatio_) {
                    Array<Float> delayRatio = fArray(parSlice1, Slice(), Slice())
                        - fArray(parSlice2, Slice(), Slice());
                    *par_[chunk] = delayRatio / 1.0e-12;
                } else {
                    *par_[chunk] = fArray(parSlice1, Slice(), Slice()) / 1.0e-12;
                }
            } else {
                throw(AipsError("delay rate has no meaning for this table"));
            }
            break;
        }
        case PMS::DISP_DELAY: {
            if (calType_.startsWith("Fringe") && !parsAreComplex()) {
                Cube<Float> fArray = cti.fparam();
                if (polnRatio_) {
                    Array<Float> dispDelayRatio = fArray(parSlice1, Slice(), Slice())
                        - fArray(parSlice2, Slice(), Slice());
                    // Divisor from PlotCal.cc
                    *par_[chunk] = dispDelayRatio / 1.334537;
                } else {
                    *par_[chunk] = fArray(parSlice1, Slice(), Slice()) / 1.334537;
                }
            } else {
                throw(AipsError("dispersive delay has no meaning for this table"));
            }
            break;
        }
        case PMS::OPAC: {
            if (!parsAreComplex() && calType_.contains("Opac")) {
                Cube<Float> fArray = cti.fparam();
                *par_[chunk] = fArray(parSlice1, Slice(), Slice());
            } else {
                throw(AipsError("opacity has no meaning for this table"));
            }
            break;
        }
        case PMS::SWP: {   // "SPGAIN" in plotcal
            if ( !parsAreComplex() && calType_.contains("EVLASWPOW")) {
                Cube<Float> fArray = cti.fparam();
                if (polnRatio_) {
                    Array<Float> swpRatio = fArray(parSlice1, Slice(), Slice()) /
                        fArray(parSlice2, Slice(), Slice());
                    checkRatioArray(swpRatio, chunk);
                    *par_[chunk] = swpRatio;
                } else {
                    *par_[chunk] = fArray(parSlice1, Slice(), Slice());
                }
            } else {
                throw(AipsError("SwPower has no meaning for this table"));
            }
            break;
        }
        case PMS::TSYS: {
            if ((!parsAreComplex()) &&
                (calType_.contains("EVLASWPOW") || calType_.contains("TSYS"))) {
                Cube<Float> fArray = cti.fparam();
                if (polnRatio_) {
                    Array<Float> tsysRatio = fArray(parSlice1, Slice(), Slice()) /
                        fArray(parSlice2, Slice(), Slice());
                    checkRatioArray(tsysRatio, chunk);
                    *par_[chunk] = tsysRatio;
                } else {
                    *par_[chunk] = fArray(parSlice1, Slice(), Slice());
                }
            } else {
                throw(AipsError("Tsys has no meaning for this table"));
            }
            break;
        }
        case PMS::SNR: {
            if (polnRatio_) {
                Array<Float> snrRatio = cti.snr()(parSlice1, Slice(), Slice()) /
                    cti.snr()(parSlice2, Slice(), Slice());
                checkRatioArray(snrRatio, chunk);
                *snr_[chunk] = snrRatio;
            } else {
                *snr_[chunk] = cti.snr()(parSlice1, Slice(), Slice());
            }
            break;
        }
        case PMS::TEC: {
            if (!parsAreComplex() && (calType_[0] == 'F') && (calType_ != "Fringe Jones")) {
                Cube<Float> fArray = cti.fparam();
                *tec_[chunk] = fArray(parSlice1, Slice(), Slice()) / 1e16;
            } else {
                throw(AipsError("TEC has no meaning for this table"));
            }
            break;
        }
        case PMS::FLAG: {
            if (polnRatio_) {
                *flag_[chunk] = cti.flag()(parSlice1, Slice(), Slice()) |
                    cti.flag()(parSlice2, Slice(), Slice());
            } else {
                *flag_[chunk] = cti.flag()(parSlice1, Slice(), Slice());
            }
            break;
        }
        /*
        case PMS::WT: {
            *wt_[chunk] = cti.weightMat();
            break;
        }
        case PMS::AZ0:
        case PMS::EL0: {
            Vector<Double> azel;
            cti.azel0Vec(cti.time()(0),azel);
            az0_(chunk) = azel(0);
            el0_(chunk) = azel(1);
            break;
        }
        case PMS::HA0: 
            ha0_(chunk) = cti.hourang(cti.time()(0))*12/C::pi;  // in hours
            break;
        case PMS::PA0: {
          pa0_(chunk) = cti.parang0(cti.time()(0))*180.0/C::pi; // in degrees
          if (pa0_(chunk)<0.0) pa0_(chunk)+=360.0;
          break;
        }
        */
        case PMS::ANTENNA: {
            antenna_[chunk]->resize(nAnt_);
            indgen(*antenna_[chunk]);
        break;
        }
        /*
        case PMS::AZIMUTH:
        case PMS::ELEVATION: {
            Matrix<Double> azel;
            cti.azelMat(cti.time()(0),azel);
            *az_[chunk] = azel.row(0);
            *el_[chunk] = azel.row(1);
            break;
        }
        case PMS::PARANG:
            *parang_[chunk] = cti.feed_pa(cti.time()(0))*(180.0/C::pi); //degrees
            break;
        case PMS::ROW: {
            *row_[chunk] = cti.rowIds();
            break;
        }
        */
        case PMS::OBSERVATION: {
          (*obsid_[chunk]).resize(1);
          *obsid_[chunk] = cti.thisObs();
          break;
        }
        case PMS::INTENT: {
          // metadata axis that always gets loaded so don't want to throw exception
          break;
        }
        case PMS::ATM:
        case PMS::TSKY: 
        case PMS::IMAGESB: {
          casacore::Int spw = cti.thisSpw();
          casacore::Int scan = cti.thisScan();
          casacore::Vector<casacore::Double> freqsGHz = cti.freq()/1e9;
          casacore::Vector<casacore::Double> curve(1, 0.0);
          if (axis == PMS::ATM) { 
              if (plotmsAtm_) {
                  plotmsAtm_->calcAtmTskyCurve(curve, spw, scan, freqsGHz);
              }
              *atm_[chunk] = curve;
          } else if (axis == PMS::TSKY) {
              if (plotmsAtm_) {
                  plotmsAtm_->calcAtmTskyCurve(curve, spw, scan, freqsGHz);
              }
              *tsky_[chunk] = curve;
          } else {
              if (plotmsAtm_) {
                  plotmsAtm_->calcImageCurve(curve, spw, scan, freqsGHz);
              }
              *imageSideband_[chunk] = curve;
          }
          break;
        }
        default:
          throw(AipsError("Axis choice not supported for Cal Tables"));
          break;
    }
}

void CalCache::flagToDisk(const PlotMSFlagging& flagging,
    Vector<Int>& flchunks, Vector<Int>& flrelids,
    Bool flag, PlotMSIndexer* indexer, int dataIndex ) {
  
  // Sort the flags by chunk:
  Sort sorter;
  sorter.sortKey(flchunks.data(),TpInt);
  sorter.sortKey(flrelids.data(),TpInt);
  Vector<uInt> order;
  uInt nflag;
  nflag = sorter.sort(order,flchunks.nelements());

  stringstream ss;

  // Make the VisIterator writable, with selection revised as appropriate
  NewCalTable* ct = new NewCalTable(filename_, Table::Update, Table::Plain);
  NewCalTable* selct = new NewCalTable();
  Vector<Vector<Slice> > chansel;
  Vector<Vector<Slice> > corrsel;
  selection_.apply(*ct, *selct, chansel, corrsel);

  Bool readonly(False); // write access for flagging
  setUpCalIter(*selct, readonly);
  ci_p->reset();

  Int iflag(0);
  for (Int ichk=0;ichk<nChunk_;++ichk) {
    if (ichk!=flchunks(order[iflag]) && !ci_p->pastEnd())
      // nothing to flag this chunk, just advance
      ci_p->next();
    else {
      // This chunk requires flag-setting
      Int ifl(iflag);
      
      // Get bits we need from the table
      Cube<Bool> ctflag;
      Vector<Int> channel,a1,a2;
      ci_p->flag(ctflag);
      ci_p->chan(channel);
      ci_p->antenna1(a1);
      ci_p->antenna2(a2);

      // Apply poln selection
      Int npar;
      String pol = selection_.corr();
      if (pol=="" || pol=="RL" || pol=="XY" || pol=="/") { // both axes
        npar = ctflag.shape()(0);
      } else { // poln selection using calParSlice
        String paramAxis = toVisCalAxis(PMS::FLAG);
        npar = getParSlice(paramAxis, pol).length();
      }    
      Int nchan = channel.nelements();
      Int nrow = ci_p->nrow();

      if (True) {
        Int currChunk=flchunks(order[iflag]);
        Double time=getTime(currChunk,0);
        Double cttime=ci_p->time()(0);
        Int spw=Int(getSpw(currChunk,0));
        Int ctspw=ci_p->thisSpw();
        Int field=Int(getField(currChunk,0));
        Int ctfld=ci_p->thisField();
        ss << "Time diff:  " << time-cttime << " " << time  << " " << cttime << "\n";
        ss << "Spw diff:   " << spw-ctspw   << " " << spw   << " " << ctspw  << "\n";
        ss << "Field diff: " << field-ctfld << " " << field << " " << ctfld  << "\n";
      }

      // Apply all flags in this chunk to this VB
      ifl=iflag;
      while (ifl<Int(nflag) && flchunks(order[ifl])==ichk) {
        Int currChunk=flchunks(order[ifl]);
        Int irel=flrelids(order[ifl]);
        Slice par1,chan,bsln;
        Slice par2 = Slice();

        // Set flag range on par axis:
        if (netAxesMask_[dataIndex](0) && !flagging.corrAll()) {
          // A specific single par
          if (pol=="" || pol=="RL" || pol=="XY") {  // flag both axes
            Int ipar=indexer->getIndex1000(currChunk,irel);
            par1 = Slice(ipar,1,1);
          } else if (polnRatio_) {
            par1 = getParSlice(toVisCalAxis(PMS::AMP), "R");
            par2 = getParSlice(toVisCalAxis(PMS::AMP), "L");
          } else {
            par1 = getParSlice(toVisCalAxis(PMS::AMP), pol);
          }
        } else {
          // all on par axis
          par1 = Slice(0,npar,1);
        }

        // Set Flag range on channel axis:
        if (netAxesMask_[dataIndex](1) && !flagging.channel()) {
          // A single specific channel
          Int ichan=indexer->getIndex0100(currChunk,irel);
          chan=Slice(ichan,1,1);
        } else {
          // Extend to all channels
          chan=Slice(0,nchan,1);
        }

        // Set Flags on the baseline axis:
        Int thisA1=Int(getAnt1(currChunk,indexer->getIndex0010(currChunk,irel)));
        Int thisA2=Int(getAnt2(currChunk,indexer->getIndex0010(currChunk,irel)));
        if (netAxesMask_[dataIndex](2) &&
            !flagging.antennaBaselinesBased() &&
            thisA1>-1 ) {
          // i.e., if baseline is an explicit data axis,
          //       full baseline extension is OFF
          //       and the first antenna in the selected point is > -1
          // Do some variety of detailed per-baseline flagging
          for (Int irow=0;irow<nrow;++irow) {
            if (thisA2>-1) {
              // match a baseline exactly
              if (a1(irow)==thisA1 &&
                  a2(irow)==thisA2) {
                ctflag(par1,chan,Slice(irow,1,1)) = flag;
                if (par2.length() > 0)
                  ctflag(par2,chan,Slice(irow,1,1)) = flag;
                break;  // found the one baseline, escape from for loop
              }
            } else {
              // either antenna matches the one specified antenna
              //  (don't break because there will be more than one)
              if (a1(irow)==thisA1 ||
                  a2(irow)==thisA1) {
                ctflag(par1,chan,Slice(irow,1,1)) = flag;
                if (par2.length() > 0)
                  ctflag(par2,chan,Slice(irow,1,1)) = flag;
              }
            }
          }
        } else {
          // Set flags for all baselines, because the plot
          //  is ordinarily implicit in baseline, we've turned on baseline
          //  extension, or we've avaraged over all baselines
          bsln=Slice(0,nrow,1);
          ctflag(par1,chan,bsln) = flag;
          if (par2.length() > 0)
            ctflag(par2,chan,bsln) = flag;
        }

      ++ifl;
      } // while
      
      // Put the flags back into the MS
      wci_p->setflag(ctflag);
      
      // Advance to the next vb
      if (!ci_p->pastEnd())
        ci_p->next();
      else
        // we are done, so escape chunk loop
        break;

      // step over the flags we've just done
      iflag=ifl;
      
      // Escape if we are already finished
      if (uInt(iflag)>=nflag) break;
    } // flaggable chunk
  } // ichk

  // Delete the NCTs and VisIter so lock is released
  if (ct != nullptr) {
    delete ct;
    ct = nullptr;
  }
  if (selct != nullptr) {
    delete selct;
    selct = nullptr;
  }
  if (wci_p != nullptr) {
    delete wci_p;
    wci_p = nullptr;
  }
  ci_p = nullptr;

  logFlag(ss.str());
}

// ======================== end NewCalTable ==========================

// ======================== CalTable ==========================

void CalCache::countChunks(Int nchunks, vector<PMS::Axis>& loadAxes,
    vector<PMS::DataColumn>& loadData, ThreadCommunication* thread) {
  // for CalTable
  if (thread!=nullptr) {
    thread->setStatus("Establishing cache size.  Please wait...");
    thread->setAllowedOperations(false,false,false);
  }
  setCache(nchunks, loadAxes, loadData);
}

void CalCache::setMSname(String msname) {
  // set msname_ (with path) if valid
  Path filepath(filename_);
  String path(filepath.dirName());
  if (!path.empty()) path += "/";
  msname_ = path + msname;
  if (msname.empty() || !Table::isReadable(msname_))
    throw(AipsError("Associated MS is not available, cannot plot solutions."));
}

void CalCache::getNamesFromMS() {
  // Set antenna and field names for Locate.
  MSMetaInfoForCal msmeta(msname_);
  msmeta.antennaNames(antnames_);
  antstanames_ = antnames_;
  msmeta.fieldNames(fldnames_);
  nAnt_ = msmeta.nAnt();
}

void CalCache::setUpLoad(ThreadCommunication* thread, Slice& parSlice) {
  // common setup for CalTables
  // permit cancel in progress meter:
  if(thread != nullptr)
    thread->setAllowedOperations(false,false,true);
  logLoad("Loading chunks......");

  // initial chunk info
  chshapes_.resize(4,nChunk_);
  goodChunk_.resize(nChunk_);
  goodChunk_.set(False);

  // get selected npol
  String polSelection(selection_.corr());
  String paramAxis(toVisCalAxis(PMS::AMP));
  // getParSlice() checks for valid axis and pol sel
  if (polnRatio_) {  // just pick one for length
    parSlice = getParSlice(paramAxis, "R");
  } else { 
    parSlice = getParSlice(paramAxis, polSelection);
  }
}

void CalCache::getCalDataAxis(PMS::Axis axis, Cube<Complex>& viscube,
    Int chunk) {
  // Get axes derived from calculated data cube; 
  // poln selection (parSlice) already applied to viscube
  switch(axis) {
    case PMS::GAMP:
    case PMS::AMP: {
      Cube<Float> ampcube = amplitude(viscube);
      if (polnRatio_) checkRatioArray(ampcube, chunk);
      *amp_[chunk] = ampcube;
      break;
    }
    case PMS::GPHASE:
    case PMS::PHASE: {
      Cube<Float> phasecube = phase(viscube);
      if (polnRatio_) checkRatioArray(phasecube, chunk);
      *pha_[chunk] = phasecube;
      (*pha_[chunk]) *= Float(180.0/C::pi);
      break;
    }
    case PMS::GREAL:
    case PMS::REAL: {
      Cube<Float> realcube = real(viscube);
      if (polnRatio_) checkRatioArray(realcube, chunk);
      *real_[chunk] = realcube;
      break;
    }
    case PMS::GIMAG:
    case PMS::IMAG: {
      Cube<Float> imagcube = imag(viscube);
      if (polnRatio_) checkRatioArray(imagcube, chunk);
      *imag_[chunk] = imagcube;
      break;
    }
    default:
      throw(AipsError("Axis choice not supported for Cal Tables"));
      break;
  }
}

// ======================== end CalTable ==========================

// ======================== BPOLY ==========================
void CalCache::loadBPoly(vector<PMS::Axis>& loadAxes,
    vector<PMS::DataColumn>& loadData, ThreadCommunication* thread) {
  // Load cache for a BPoly cal table
  BJonesPolyTable pt = BJonesPolyTable(filename_, Table::Update);

  // Set ms-related data
  ROCalDescColumns calDescCol(pt);
  String msname(calDescCol.msName()(0));
  setMSname(msname); // add path
  getNamesFromMS();  // field and antenna

  // Set number of antennas
  nAnt_ = pt.maxAntenna() + 1;

  // Create a B Jones table derived from BPOLY table
  NewCalTable* nct = virtualBPoly(pt);
  parsAreComplex_ = nct->isComplex();

  // Apply selection to get selected cal table
  NewCalTable* selct = new NewCalTable();
  Vector<Vector<Slice> > chansel;
  Vector<Vector<Slice> > corrsel;
  selection_.apply(*nct, *selct, chansel, corrsel);

  // Use NewCalTable implementation
  setUpCalIter(*selct, true);
  countChunks(*ci_p, loadAxes, loadData, thread);
  loadCalChunks(*ci_p, loadAxes, thread);
}

NewCalTable* CalCache::virtualBPoly(BJonesPolyTable& polyTable) {
  // Returns B Jones NewCalTable derived from BPOLY; based on BPoly::loadMemCalTable.

  // Generate a NCT in memory to hold the BPOLY as a B
  NewCalTable* nct = new NewCalTable("BPolyAsB.tmp", VisCalEnum::COMPLEX, "B Jones", msname_, false);

  // Frequency info from MS
  Vector<Vector<Double>> mschanfreq; // [nspw, nchan]
  getChanFreqsFromMS(mschanfreq);

  // Ensure sort on TIME, so CalSet is filled in order
  Block <String> sortCol(3);
  sortCol[0]="CAL_DESC_ID";
  sortCol[1]="TIME";
  sortCol[2]="ANTENNA1";
  polyTable.sort2(sortCol);

  Int nrows = polyTable.nRowMain();
  Int nDesc = polyTable.nRowDesc();

  // Spws to be calibrated by each caldesc
  Vector<Int> spwmap(nDesc,-1);
  for (Int idesc = 0; idesc < nDesc; ++idesc) {
    CalDescRecord calDescRec(polyTable.getRowDesc(idesc));
    Vector<Int> spwId;
    calDescRec.getSpwId(spwId);
    Int currSpw = spwId(0);
    spwmap(idesc) = currSpw;

    // Set SPW subtable freqs
    Vector<Double> freq = mschanfreq(currSpw);
    nct->setSpwFreqs(currSpw, freq);
  }
  
  // Solve arrays, so we can fill them
  Cube<Complex> cparam;
  Cube<Bool> flag;
  Cube<Float> err, snr;
  bool cubeInit(false);

  // Attach a calibration table columns accessor
  BJonesPolyMCol maincol(polyTable);

  for (Int row = 0; row < nrows; row++) {
    // Extract the polynomial coefficients in amplitude and phase
    Int nAmp = maincol.nPolyAmp().asInt(row);
    Int nPhase = maincol.nPolyPhase().asInt(row);
    Array<Double> ampCoeffArray, phaseCoeffArray;
    maincol.polyCoeffAmp().get(row, ampCoeffArray);
    maincol.polyCoeffPhase().get(row, phaseCoeffArray);

    Matrix<Double> ampCoeff(nAmp, 2);
    IPosition ampPos = ampCoeffArray.shape();
    ampPos = 0;
    for (Int k = 0; k < 2 * nAmp; k++) {
      ampPos.setLast(IPosition(1, k));
      ampCoeff(k % nAmp, k / nAmp) = ampCoeffArray(ampPos);
    };

    Matrix<Double> phaseCoeff(nPhase, 2);
    IPosition phasePos = phaseCoeffArray.shape();
    phasePos = 0;
    for (Int k = 0; k < 2 * nPhase; k++) {
      phasePos.setLast(IPosition(1, k));
      phaseCoeff(k % nPhase, k / nPhase) = phaseCoeffArray(phasePos);
    };

    // Get frequencies for this spw
    Int nPol(2);
    Int thisDesc = maincol.calDescId().asInt(row);
    Int thisSpw = spwmap(thisDesc);
    Vector<Double> freq = mschanfreq(thisSpw);
    Int nChan = freq.nelements();

    // Extract the valid domain for the polynomial
    Vector<Double> freqDomain(2);
    maincol.validDomain().get(row, freqDomain);
    Double x1 = freqDomain(0);
    Double x2 = freqDomain(1);

    Complex factor = maincol.scaleFactor().asComplex(row);
    Int thisAnt1 = maincol.antenna1().asInt(row);

    // Resize and initialize solve arrays
    if (!cubeInit) {
      cparam.resize(nPol, nChan, nAnt_);
      cparam.set(Complex(1.0));
      flag.resize(nPol, nChan, nAnt_);
      flag.set(true);
      err.resize(nPol, nChan, nAnt_);
      err.set(0.0);
      snr.resize(nPol, nChan, nAnt_);
      snr.set(1.0);
      cubeInit = true;
    }

    for (Int pol = 0; pol < 2; pol++) {
      Vector<Double> ac(ampCoeff.column(pol));
      Vector<Double> pc(phaseCoeff.column(pol));
      
      // Only do calculation if coeffs are non-zero
      if (anyNE(ac, Double(0.0)) || anyNE(pc, Double(0.0)) ) {
        for (Int chan = 0; chan < nChan; ++chan) {
          Double ampval(1.0), phaseval(0.0);
          // Calculate Cheby if freq in domain
          Double thisFreq(freq(chan));
          if ((thisFreq >= x1) && (thisFreq <= x2)) {
            ampval = getChebVal(ac, x1, x2, thisFreq);
            phaseval = getChebVal(pc, x1, x2, thisFreq);
            cparam(pol, chan, thisAnt1) = factor *
              Complex(exp(ampval)) * Complex(cos(phaseval), sin(phaseval));
            flag(pol, chan, thisAnt1) = false;
          } else {
            // Unflagged unit calibration for now
            cparam(pol, chan, thisAnt1) = Complex(1.0);
            flag(pol, chan, thisAnt1) = false;
          }
        } // chan
      }
    } // pol

    // Every nAnt rows, store the result
    if ((row + 1) % nAnt_ == 0) {
      Double thisTime = maincol.time().asdouble(row);
      Double thisInterval = maincol.interval().asdouble(row);
      Int thisField = maincol.fieldId().asInt(row);
      Int thisObs = maincol.obsId().asInt(row);
      Array<Int> refant = maincol.refAnt().get(row);
      IPosition first(refant.shape().size(), 0);
      Int thisRefant = refant(first);
      Vector<Int> ant1list; // NewCalTable generates antenna ids based on nrows

      nct->fillAntBasedMainRows(nAnt_, thisTime, thisInterval, thisField,
        thisSpw, thisObs, ant1list, thisRefant, cparam, flag, err, snr);

      // reset arrays next loop
      cubeInit = false;
    }
  }   // rows

  return nct;
}

void CalCache::getChanFreqsFromMS(Vector< Vector<Double> >& mschanfreqs) {
  // shape is (nspw, nchan)
  MeasurementSet ms(msname_);
  MSColumns mscol(ms);
  uInt nspw = mscol.spectralWindow().nrow();
  mschanfreqs.resize(nspw);
  for (uInt spw=0; spw<nspw; ++spw) {
    mschanfreqs(spw) = mscol.spectralWindow().chanFreq().get(spw);
  }
}

Double CalCache::getChebVal(const Vector<Double>& coeff, const Double& xinit,
  const Double& xfinal, const Double& x) {
// from synthesis/MeasurementComponents/BPoly.cc
// Compute a Chebyshev polynomial value using the CLIC library
// Input:
//    coeff       const Vector<Double>&       Chebyshev coefficients
//    xinit       const Double&               Domain start
//    xfinal      const Double&               Domain end
//    x           const Double&               x-ordinate
// Output:
//    getChebVal  Double                      Chebyshev polynomial value
//
  // Re-scale x-ordinate
  Double xcap = ((x - xinit) - (xfinal - x)) / (xfinal - xinit);

  // Compute polynomial
  Int deg = coeff.shape().asVector()(0);
  Vector<Double> val(deg);
  Bool check;
  Int checkval;
  cheb(&deg, &xcap, val.getStorage(check), &checkval);

  Double soly(0.0);
  for (Int mm = 0; mm < deg; mm++){
    soly += coeff[mm] * val[mm];
  }

  return soly;
}

// ======================== end BPOLY ==========================

// ======================== GSPLINE ==========================
void CalCache::loadGSpline(vector<PMS::Axis>& loadAxes,
    vector<PMS::DataColumn>& loadData, ThreadCommunication* thread) {
  GJonesSplineTable ct = GJonesSplineTable(filename_);
  GJonesSplineTable selct(ct);
  if (!selection_.timerange().empty()) {
    logWarn("PlotMS::load_cache",
        "Time selection not supported for GSPLINE calibration tables");
    selection_.setTimerange("");
  }
  // chansel not applicable, corrsel done with parSlice
  Vector<Vector<Slice> > chansel, corrsel;
  selection_.apply(ct, selct, chansel, corrsel);
  Vector<Int> selAnts = selection_.getSelectedAntennas1();

  ROGJonesSplineMCol mainCol(selct);
  ROCalDescColumns calDescCol(selct);
  String msname(calDescCol.msName()(0));
  setMSname(msname); // add path
  getNamesFromMS();  // field and antenna

  // count and load chunks
  Int nsample(1000); // make time samples to load cache
  countChunks(nsample, loadAxes, loadData, thread);
  loadCalChunks(mainCol, calDescCol, nsample, loadAxes, selAnts, thread);
}

void CalCache::loadCalChunks(ROGJonesSplineMCol& mcol, ROCalDescColumns& dcol,
    Int nsample, const vector<PMS::Axis> loadAxes, Vector<int>& selectedAnts,
    ThreadCommunication* thread) {
  // GSPLINE does not load per chunk or row as in other cal tables. 
  // Its "chunks" are per-time sample, generated from SPLINE_KNOTS_AMP/PHASE.
  // Therefore:
  //   Scalar columns are read once and all values plotted per time sample.
  //   Columns that normally have one value per chunk (field, spw, scan, time,
  //   and interval) will use the value in the column if all are the same; 
  //   else multiple values result in setting the value to -1 (for locate).
  //   The exception is field, which plots the field used for solutions.
  Slice parslice, parslice2;
  setUpLoad(thread, parslice);
  Int nPol = parslice.length();  // for chunk shapes
  if (polnRatio_) {
    parslice2 = getParSlice(toVisCalAxis(PMS::AMP), "L");
  }

  // load main table metadata once; use vector/value for every timestamp
  // field
  Vector<Int> fields(mcol.fieldId().getColumn());
  Int nItems = GenSort<Int>::sort(fields, Sort::Ascending, Sort::NoDuplicates);
  Int firstField(fields(0));
  Int fieldForSolve = (firstField < 0 ? 0 : firstField);
  Int fieldForCache = (nItems>1 ? -1 : firstField);
  // spw
  Array<Int> spws(dcol.spwId().getColumn());
  nItems = GenSort<Int>::sort(spws, Sort::Ascending, Sort::NoDuplicates);
  Int firstSpw(spws(IPosition(2,0,0)));
  Int spwForSolve = (firstSpw < 0 ? 0 : firstSpw);
  Int spwForCache = (nItems>1 ? -1 : firstSpw);
  // scan
  Vector<Int> scans(mcol.scanNo().getColumn());
  nItems = GenSort<Int>::sort(scans, Sort::Ascending, Sort::NoDuplicates);
  Int firstScan(scans(0));
  Int scanForSolve(firstScan < 0 ? 0 : firstScan);
  Int scanForCache = (nItems>1 ? -1 : firstScan);
  // antenna1
  Vector<Int> ant1(selectedAnts);  // use selected antenna1
  Int nSelAnts = ant1.size();
  Int nAnt(nSelAnts);
  if (nAnt == 0) {  // no ant1 selection, get from main table
    ant1 = mcol.antenna1().getColumn();
    nItems = GenSort<Int>::sort(ant1, Sort::Ascending, Sort::NoDuplicates);
    ant1.resize(nItems, true);
    nAnt = nItems;
  }

  // obsid
  Vector<Int> obsid(mcol.obsId().getColumn());
  nItems = GenSort<Int>::sort(obsid, Sort::Ascending, Sort::NoDuplicates);
  obsid.resize(nItems, true);
  Int obsForSolve(obsid(0) < 0 ? 0 : obsid(0));
  // feed1
  Vector<Int> feed1(mcol.feed1().getColumn());
  nItems = GenSort<Int>::sort(feed1, Sort::Ascending, Sort::NoDuplicates);
  feed1.resize(nItems, true);
  // interval
  Vector<Double> intervals(mcol.interval().getColumn());
  nItems = GenSort<Double>::sort(intervals, Sort::Ascending, Sort::NoDuplicates);
  intervals.resize(nItems, true);
  Double interval = (nItems>1 ? -1 : intervals(0));
  // poln
  Vector<Int> polns(parslice.length());
  if (polnRatio_) {
    polns = -1;
  } else {
    Int start(parslice.start()), inc(parslice.inc());
    indgen(polns, start, inc);
  }

  // set up gspline
  MSMetaInfoForCal msmeta(msname_);
  GJonesSpline* gspline = new GJonesSpline(msmeta);
  String mode(mcol.polyMode()(0));
  Record rec;  // for solving params
  rec.define("table", filename_);
  rec.define("apmode", mode);
  gspline->setApply(rec); // set up calbuffer for calcPar
  gspline->setSolve(rec); // set mode

  // get first row for chosen field
  uInt row;
  Vector<Int> fieldcol(mcol.fieldId().getColumn());
  for (row=0; row<fieldcol.size(); ++row)
    if (fieldcol(row) == fieldForSolve) break;

  MFrequency refFreq = mcol.refFreqMeas()(row)(IPosition(3,0,0,0));
  Double refFreqHz = refFreq.get("Hz").getValue();
  Vector<Double> freq(1, refFreqHz);    // for setMeta

  // Create 1000 time samples from spline knots: see PlotCal::virtualGSpline
  Vector<Double> splineKnots, times(nsample);
  // get splineKnots for that row
  if (mode.contains("AMP") || mode.contains("A&P")) {
    splineKnots = mcol.splineKnotsAmp()(row);
  } else if (mode.contains("PHAS") || mode.contains("A&P")) {
    splineKnots = mcol.splineKnotsPhase()(row);
  }
  // make samples based on spline knots
  Double dt((max(splineKnots)-min(splineKnots)) / Double(nsample));
  Double mintime(splineKnots(0) + (dt/2.0));
  for (Int sample=0; sample<nsample; sample++) {
    times(sample) = mintime + sample*dt;
  }

  // Now ready to load cache
  for (Int sample=0; sample<nsample; sample++) {
    // set time and field for calcPar
    Double time = times(sample);
    gspline->setMeta(obsForSolve, scanForSolve, time, spwForSolve, freq, fieldForSolve);
    gspline->doCalcPar();

    // Cache the data shapes
    chshapes_(0,sample) = nPol;
    chshapes_(1,sample) = 1; // nChan
    chshapes_(2,sample) = nAnt;
    chshapes_(3,sample) = 1; 
    goodChunk_(sample) = True;
    // load axes for each row
    for(unsigned int i = 0; i < loadAxes.size(); i++) {
      PMS::Axis axis = loadAxes[i];
      // slice viscube for poln selection
      Slicer slicer1(Slicer(parslice, Slice(), Slice()));
      if (PMS::axisIsData(axis)) { // amp, phase, real, imag 
        Cube<Complex> cpar, viscube;
        cpar = gspline->currCPar(); 
        viscube = cpar(slicer1);
        if (polnRatio_) {
          Slicer slicer2(Slicer(parslice2, Slice(), Slice()));
          viscube /= cpar(slicer2);
        }
        if (nAnt == nSelAnts)  // get selected rows 
          getSelectedCube(viscube, selectedAnts);
        getCalDataAxis(axis, viscube, sample);
      } else {
        switch(axis) {
          case PMS::SCAN:
            scan_(sample) = scanForCache;
            break;
          case PMS::FIELD:
            field_(sample) = fieldForCache;
            break;
          case PMS::TIME: 
            time_(sample) = time;
            break;
          case PMS::TIME_INTERVAL:
            timeIntr_[sample] = interval;
            break;
          case PMS::SPW:
            spw_(sample) = spwForCache;
            break;
          case PMS::CORR: {
            corr_[sample]->resize(nPol);
            String pol = selection_.corr();
            if (pol=="" || pol=="RL" || pol=="XY") {
              indgen(*corr_[sample]);
            } else if (pol=="/") {
              corr_[sample]->resize(1);
              corr_[sample]->set(-1);
            } else {  // R/X or L/Y
              Int poln = ((pol=="R" || pol=="X") ? 0 : 1);
              corr_[sample]->resize(1);
              corr_[sample]->set(poln);
            }
            break;
          }
          case PMS::ANTENNA1:
            *antenna1_[sample] = ant1;
            break;
          case PMS::ANTENNA:
            *antenna_[sample] = ant1;
            break;
          case PMS::FLAG: {
            Cube<Bool> parOK = gspline->currParOK();
            // OK=true means flag=false
            Cube<Bool> flagcube(!parOK(slicer1));
            if (polnRatio_) {
              Slicer slicer2(Slicer(parslice2, Slice(), Slice()));
              flagcube &= !parOK(slicer2);
            }
            if (nAnt == nSelAnts)  // get selected rows 
              getSelectedCube(flagcube, selectedAnts);
            *flag_[sample] = flagcube;
            break;
          }
          case PMS::ROW: {
            Vector<rownr_t> sampleRow(nAnt, sample);
            *row_[sample] = sampleRow;
            break;
          }
          case PMS::OBSERVATION:
            *obsid_[sample] = obsid;
            break;
          case PMS::FEED1:
            *feed1_[sample] = feed1;
            break;
          default:  // invalid axes weeded out in checkAxes
            break;
          }
      }
    }

    // If a thread is given, update it.
    if(thread != nullptr) {
      double progress = ((double)sample) / nsample;
      thread->setProgress((unsigned int)((progress * 100) + 0.5));
    }
  }
  if (divZero_)
    logWarn("CalCache::loadIt", "Caught divide-by-zero exception in ratio plots; result(s) set to 1.0 and flagged");
}

void CalCache::checkAxes(const vector<PMS::Axis>& loadAxes) {
  // trap user-requested axes that are invalid for GSPLINE
  for(unsigned int i = 0; i < loadAxes.size(); i++) {
    PMS::Axis axis(loadAxes[i]);
    switch (axis) {
      case PMS::SCAN:
      case PMS::FIELD:
      case PMS::TIME:
      case PMS::SPW:
      case PMS::TIME_INTERVAL:
      case PMS::CORR:
      case PMS::ANTENNA1:
      case PMS::ANTENNA:  // same as antenna1
      case PMS::AMP:
      case PMS::GAMP:
      case PMS::PHASE:
      case PMS::GPHASE:
      case PMS::REAL:
      case PMS::GREAL:
      case PMS::IMAG:
      case PMS::GIMAG:
      case PMS::FLAG:
      case PMS::OBSERVATION:
      case PMS::FEED1: {
        // allowed
        break;
      }
      case PMS::CHANNEL:
      case PMS::FREQUENCY:
      case PMS::SNR: {
        String msg("GSPLINE plotting does not support " + PMS::axis(axis) + " axis");
        if (axis==PMS::ANTENNA) msg += "; use ANTENNA1";
        throw(AipsError(msg));
        break;
      }
      case PMS::ANTENNA2:
      case PMS::BASELINE:
      case PMS::ROW:
      case PMS::DELAY:
      case PMS::DELAY_RATE:
      case PMS::DISP_DELAY:
      case PMS::OPAC:
      case PMS::SWP:
      case PMS::TSYS:
      case PMS::TEC:
      case PMS::INTENT: {
        String msg(PMS::axis(axis) + " has no meaning for this table");
        throw(AipsError(msg));
        break;
      }
      default:
        throw(AipsError("Axis choice not supported for Cal Tables"));
        break;
    }
  }
}

template<class T>
void CalCache::getSelectedCube(Cube<T>& inputCube, const Vector<Int>& selectedRows) {
  // replaces input cube with cube selected by rows in vector
  Cube<T> selectedCube;
  for (uInt irow=0; irow<selectedRows.size(); ++irow) {
    Slice rowSlice = Slice(selectedRows(irow));
    Slicer rowSlicer = Slicer(Slice(), Slice(), rowSlice);
    Cube<T> concatCube = concatenateArray(selectedCube, inputCube(rowSlicer));
    selectedCube.resize();
    selectedCube = concatCube;
  }
  inputCube.resize();
  inputCube = selectedCube;
}

// ======================== end GSPLINE ==========================


String CalCache::toVisCalAxis(PMS::Axis axis) {
    switch (axis) {
        // FLAG and SNR have same shape as AMP 
        // and should be sliced the same way
        case PMS::AMP:
        case PMS::GAMP:
        case PMS::FLAG:
        case PMS::SNR:
            if (calType_.contains("EVLASWP")) return "GAINAMP";
            if (calType_.contains("TSYS")) return "TSYS";
            if (calType_[0] == 'K' && !calType_.startsWith("KAntPos")) 
                return "DELAY";
            if (calType_ == "F Jones") return "TEC";
            if (calType_ == "Fringe Jones") return "DELAY";
            if (calType_ == "TOpac") return "OPAC";
            return "AMP";
            break;
        case PMS::PHASE:
        case PMS::GPHASE:
            return "PHASE";
            break;
        case PMS::REAL:
        case PMS::GREAL:
            return "REAL";
            break;
        case PMS::IMAG:
        case PMS::GIMAG:
            return "IMAG";
            break;
        case PMS::DELAY_RATE:
            return "RATE";
            break;
        case PMS::DISP_DELAY:
            return "DISP";
            break;
        default:
            return PMS::axis(axis);
            break;
    }
}

Slice CalCache::getParSlice(String axis, String polnSel) {
    Slice parSlice = Slice();
    try {
        parSlice = viscal::calParSlice(filename_, axis, polnSel);
    } catch(AipsError& err) {
        if (err.getMesg().contains("Unsupported value type")) {
            // Message a bit vague at top level, add some explanation
            String errMsg = err.getMesg() + ". Invalid axis or polarization selection for cal table type.";
            throw (AipsError(errMsg));
        } else if (calType_ == "M Mueller") {
            if (polnSel.empty()) {
                return Slice(0, 2, 1); // full selection
            }
            // include default "RL" when not set by user
            std::vector<casacore::String> valid_poln = {"R", "X", "RL"};

            for (auto& poln : valid_poln) {
                if (poln == polnSel) {
                    return Slice(0, 2, 1); // full selection
                }
            }
            throw (AipsError("Invalid polarization selection for cal table type."));
        } else { // unsupported cal type
            throw(AipsError(err));
        }
    }
    return parSlice;
}

void CalCache::checkRatioArray(Array<Float>& array, Int chunk) {
    // When array is result of division, check for division by zero (element is infinite)
    // Reset value to 1.0, flag, and set divZero_ for warning
    Cube<Float> ratioCube;
    ratioCube.reference(array);
    Cube<Bool> flags;
    flags.reference(*flag_[chunk]);

    IPosition cubeShape = ratioCube.shape();
    for (uInt i=0; i<cubeShape[0]; ++i) {
        for (uInt j=0; j<cubeShape[1]; ++j) {
            for (uInt k=0; k<cubeShape[2]; ++k) {
                if (isInf(ratioCube(i,j,k))) {
                    ratioCube(i,j,k) = 1.0;
                    flags(i,j,k) = True;
                    divZero_ = True;
                }
            }
        }
    }
}

}
