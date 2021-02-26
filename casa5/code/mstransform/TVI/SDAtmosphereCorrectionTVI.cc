//# SDAtmosphereCorrectionTVI.h: This file contains the implementation of the SDAtmosphereCorrectionTVI class.
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
#include <mstransform/TVI/SDAtmosphereCorrectionTVI.h>

#include <cmath>
#include <functional>
#include <iomanip>
#include <vector>

#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/BasicSL/Constants.h>
#include <casacore/casa/BasicSL/String.h>
#include <casacore/casa/Logging/LogIO.h>
#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/Exceptions/Error.h>
#include <casacore/casa/Arrays/ArrayIter.h>
#include <casacore/casa/OS/File.h>
#include <casacore/casa/Utilities/Sort.h>
#include <casacore/casa/Utilities/BinarySearch.h>
#include <casacore/measures/Measures/Stokes.h>
#include <casacore/scimath/Functionals/Interpolate1D.h>
#include <casacore/scimath/Functionals/ScalarSampledFunctional.h>
#include <casacore/ms/MeasurementSets/MSDataDescColumns.h>
#include <casacore/ms/MeasurementSets/MSPolColumns.h>
#include <casacore/ms/MSOper/MSMetaData.h>
#include <casacore/tables/TaQL/ExprNode.h>

#include <msvis/MSVis/VisBufferComponents2.h>
#include <msvis/MSVis/VisibilityIteratorImpl2.h>

#include <mstransform/TVI/UtilsTVI.h>

using namespace casacore;

namespace {

constexpr double kValueUnset = std::nan("");
inline bool isValueUnset(double const v) {return std::isnan(v);}

class ScopeGuard {
  public:
    explicit ScopeGuard(std::function<void()> f) : func_(f) {}
    ScopeGuard(ScopeGuard const &s) = delete;
    void operator=(ScopeGuard const &s) = delete;
    ~ScopeGuard() {
      func_();
    }
  private:
    std::function<void()> func_;
};

inline void sortTime(Vector<Double> &timeData, Vector<uInt> &indexVector) {
  Bool b = false;
  Double const *p = timeData.getStorage(b);
  ScopeGuard guard([&]() {
    timeData.freeStorage(p, b);
  });
  LogIO os;
  // std::unique_ptr<Double[]> p(new Double[timeData.nelements()]);
  // for (uInt i = 0; i < timeData.nelements(); ++i) {
  //   p[i] = timeData[i];
  // }
  // indexVector.resize(timeData.nelements());
  // indgen(indexVector);
  Sort sort(p, timeData.nelements());
  sort.sort(indexVector, timeData.nelements());
  os << "indexVector = " << indexVector << endl;
  // uInt r = sort.unique(uniqueVector, indexVector);
  os << "indexVector (after) = " << indexVector << endl;
  // os << "uniqueVector = " << uniqueVector << endl;
  os.output() << std::setprecision(16);
  os << "timeData = " << timeData << endl;
  // os << "return value = " << r << endl;
  os << LogIO::POST;
  // timeData.freeStorage(p, b);
  // return r;
}

inline uInt makeUnique(Vector<Double> const &data, Vector<uInt> const &sortIndex, Vector<uInt> &uniqueIndex) {
  uniqueIndex.resize(data.nelements());
  assert(data.nelements() == sortIndex.nelements());
  uInt n = 0;
  uniqueIndex[n++] = sortIndex[0];
  Double prev = data[sortIndex[0]];
  for (uInt i = 0; i < data.nelements(); ++i) {
    uInt j = sortIndex[i];
    Double current = data[j];
    if (current != prev) {
      uniqueIndex[n++] = j;
    }
    prev = current;
  }
  cout << "uniqueIndex = " << uniqueIndex << endl;
  cout << "found " << n << " unique entries" << endl;
  uniqueIndex.resize(n, True);
  return n;
}

inline Vector<Double> getMedianDataPerTime(uInt const n, Vector<uInt> const &uniqueVector,
  Vector<uInt> const &indexVector, Vector<Double> const &data) {
  Vector<Double> out(n);
  assert(n > 0);
  for (uInt i = 0; i < n - 1; ++i) {
    cout << "loop " << i << endl;
    uInt const iStart = uniqueVector[i];
    uInt const iEnd = uniqueVector[i + 1];
    uInt const nElem = iEnd - iStart;
    Vector<Double> tmp(nElem);
    for (uInt j = 0; j < nElem; ++j) {
      tmp[j] = data[indexVector[j + iStart]];
    }
    out[i] = median(tmp, False, True, True);
  }
  cout << "finalize loop " << endl;
  uInt const iStart = uniqueVector[n - 1];
  uInt const iEnd = indexVector.nelements();
  uInt const nElem = iEnd - iStart;
  cout << "start " << iStart << " end " << iEnd << " nElem " << nElem << endl;
  Vector<Double> tmp(nElem);
  for (uInt j = 0; j < nElem; ++j) {
    tmp[j] = data[indexVector[j + iStart]];
  }
  out[n - 1] = median(tmp);
  return out;
}

inline std::pair<Int, Int> findNearestIndex(Vector<Double> const &data, Double const target) {
  Bool found = false;
  // data is sorted in ascending order
  Vector<Double> diff(data.nelements());
  for (uInt i = 0; i < diff.nelements(); ++i) {
    diff[i] = data[i] - data[0];
  }
  cout << "diff " << diff << endl;
  cout << "target offset " << target - data[0] << endl;
  Int index = binarySearch(found, data, target, data.nelements());
  cout << "binary search result: found " << found << " index " << index << endl;
  // if (!found) {
  //   throw AipsError("error in findNearest");
  // }
  Int prev = -1;
  Int next = -1;
  if (found) {
    prev = index;
    next = index;
  } else if (index > 0) {
    prev = index - 1;
    next = index;
  } else {
    throw AipsError("error in findNearest");
  }
  return std::make_pair(prev, next);
}

inline std::pair<Double, Double> findNearest(Vector<Double> const &data, Double const target) {
  std::pair<uInt, uInt> idx = findNearestIndex(data, target);
  return std::make_pair(data[idx.first], data[idx.second]);
}

inline Double interpolateDataLinear(Vector<Double> const &xin, Vector<Double> const &yin, Double const xout) {
  Bool found = false;
  // data is sorted in ascending order
  Int index = binarySearch(found, xin, xout, xin.nelements());
  // if (!found) {
  //   throw AipsError("error in interpolateData");
  // }
  Double yout = kValueUnset;
  if (index == 0 || found) {
    yout = yin[index];
  } else {
    yout = ((xin[index] - xout) * yin[index - 1] * (xout - xin[index - 1]) * yin[index]) / (xin[index] - xin[index - 1]);
  }
  return yout;
}

inline Vector<Double> getTrjSkySpec(atm::SkyStatus *skyStatus, unsigned int const nchan) {
  // TODO: OMP parallelization
  Vector<Double> trj(nchan);
  for (unsigned int i = 0; i < nchan; ++i) {
    trj[i] = skyStatus->getTrjSky(i).get(atm::Temperature::UnitKelvin);
  }
  return trj;
}

inline Vector<Double> getTauSpec(atm::SkyStatus *skyStatus, unsigned int const nchan, Double const airMass) {
  // TODO: OMP parallelization
  Vector<Double> tau(nchan);
  for(unsigned int i = 0; i < nchan; ++i) {
    tau[i] = (skyStatus->getDryOpacity(i).get(atm::Opacity::UnitNeper) +
              skyStatus->getWetOpacity(i).get(atm::Opacity::UnitNeper))
             * airMass;
  }
  return tau;
}

inline Vector<Double> getCorrectionFactor(
  Vector<Double> const &tSkyOn, Vector<Double> const &tSkyOff, Vector<Double> const &tauOff) {
  Vector<Double> factor(tSkyOn.nelements());
  // TODO: OMP/SIMD parallelization
  for (size_t i = 0; i < tSkyOn.nelements(); ++i) {
    factor[i] = (tSkyOn[i] - tSkyOff[i]) * std::exp(tauOff[i]);
  }
  return factor;
}

template<class T>
inline void transformData(Vector<Double> const &gainFactor, Cube<T> const &in, Cube<T> &out) {
  // TODO: OMP/SIMD parallelization
  cout << "gainfactor (" << gainFactor.shape() << ") = " << gainFactor[0] << endl;
  cout << "input array shape = " << in.shape() << endl;
  if (allEQ(gainFactor, 1.0)) {
    out.reference(in);
  } else {
    IPosition const shape = in.shape();
    out.resize(shape);
    for (ssize_t ir = 0; ir < shape[2]; ++ir) {
      for (ssize_t ic = 0; ic < shape[1]; ++ic) {
        for (ssize_t ip = 0; ip < shape[0]; ++ip) {
          out(ip, ic, ir) = in(ip, ic, ir) - gainFactor[ic];
        }
      }
    }
  }
}

} // anonymous namespace

namespace casa { //# NAMESPACE CASA - BEGIN

namespace vi { //# NAMESPACE VI - BEGIN
//////////
// SDAtmosphereCorrectionTVI
/////////
SDAtmosphereCorrectionTVI::SDAtmosphereCorrectionTVI(ViImplementation2 *inputVII,
  Record const &configuration) :
    TransformingVi2(inputVII),
    processSpwList_(),
    gainFactorList_(),
    gainFactor_(1.0),
    userPressureValue_(kValueUnset),
    userTemperatureValue_(kValueUnset),
    userRelHumidityValue_(kValueUnset),
    userPwvValue_(kValueUnset),
    offSourceTime_(),
    elevationTime_(),
    elevationData_(),
    elevationInterpolator_(),
    pwvTime_(),
    pwvData_(),
    atmTime_(),
    atmTemperatureData_(),
    atmPressureData_(),
    atmRelHumidityData_(),
    doTransform_(),
    isTdmSpw_(),
    doSmooth_(),
    channelFreqsPerSpw_(),
    channelWidthsPerSpw_(),
    atmType_(2),
    atmProfile_(),
    atmSpectralGridPerSpw_(),
    atmSkyStatusPerSpw_(),
    atmSkyStatusPtr_(nullptr),
    correctionFactor_() {
  initializeAtmosphereCorrection(configuration);

  // Initialize attached VisBuffer
  LogIO os(LogOrigin("SDAtmosphereCorrectionTVI", __func__, WHERE));
  os << "setVisBuffer start" << LogIO::POST;
  setVisBuffer(createAttachedVisBuffer(VbRekeyable));
  os << "setVisBuffer end" << LogIO::POST;
}

SDAtmosphereCorrectionTVI::~SDAtmosphereCorrectionTVI() {
}

void SDAtmosphereCorrectionTVI::origin() {
  TransformingVi2::origin();

  // Synchronize own VisBuffer
  configureNewSubchunk();

  // configure atmospheric correction
  // point appropriate SkyStatus object
  configureAtmosphereCorrection();
  updateAtmosphereModel();

  // warn if current spw is not requested to transform
  warnIfNoTransform();
}

void SDAtmosphereCorrectionTVI::next() {
  TransformingVi2::next();

  // Synchronize own VisBuffer
  configureNewSubchunk();

  // configure atmospheric correction if necessary
  // point appropriate SkyStatus object
  configureAtmosphereCorrection();
  updateAtmosphereModel();

  // warn if current spw is not requested to transform
  warnIfNoTransform();
}

void SDAtmosphereCorrectionTVI::originChunks(Bool forceRewind) {
  TransformingVi2::originChunks();

  // update SkyStatus with the closest weather measurement
  // re-calculate correction factor (dTa)
  updateAtmosphereModel();
}

void SDAtmosphereCorrectionTVI::nextChunk() {
  TransformingVi2::nextChunk();

  updateAtmosphereModel();
}

void SDAtmosphereCorrectionTVI::visibilityCorrected(Cube<Complex> & vis) const {
  if (getVii()->existsColumn(VisBufferComponent2::VisibilityCorrected)) {
    Cube<Complex> dataCube;
    getVii()->visibilityCorrected(dataCube);
    if (doTransform()) {
      transformData(correctionFactor_, dataCube, vis);
    } else {
      vis.reference(dataCube);
    }
  } else {
    vis.resize();
  }
}

void SDAtmosphereCorrectionTVI::visibilityModel(Cube<Complex> & vis) const {
  if (getVii()->existsColumn(VisBufferComponent2::VisibilityModel)) {
    Cube<Complex> dataCube;
    getVii()->visibilityModel(dataCube);
    if (doTransform()) {
      transformData(correctionFactor_, dataCube, vis);
    } else {
      vis.reference(dataCube);
    }
  } else {
    vis.resize();
  }
}

void SDAtmosphereCorrectionTVI::visibilityObserved(Cube<Complex> & vis) const {
  if (getVii()->existsColumn(VisBufferComponent2::VisibilityObserved)) {
    Cube<Complex> dataCube;
    getVii()->visibilityObserved(dataCube);
    if (doTransform()) {
      transformData(correctionFactor_, dataCube, vis);
    } else {
      vis.reference(dataCube);
    }
  } else {
    vis.resize();
  }
}

void SDAtmosphereCorrectionTVI::floatData(casacore::Cube<casacore::Float> & fcube) const {
  if (getVii()->existsColumn(VisBufferComponent2::FloatData)) {
    Cube<Float> dataCube;
    getVii()->floatData(dataCube);
    if (doTransform()) {
      transformData(correctionFactor_, dataCube, fcube);
    } else {
      fcube.reference(dataCube);
    }
  } else {
    fcube.resize();
  }
}

void SDAtmosphereCorrectionTVI::initializeAtmosphereCorrection(
  Record const &configuration
) {
  LogIO os(LogOrigin("SDAtmosphereCorrectionTVI", __func__, WHERE));

  // processing spws
  if (!configuration.isDefined("processspw")) {
    os << "ERROR in configuration: processing spw list (processspw) must be given."
       << LogIO::EXCEPTION;
  }
  processSpwList_ = configuration.asArrayInt("processspw");
  os << "processspw = " << processSpwList_ << LogIO::POST;

  // gain factor
  if (!configuration.isDefined("gainfactor")) {
    os << "ERROR in configuration: gain factor list (gainfactor) must be given."
       << LogIO::EXCEPTION;
  }
  gainFactorList_ = configuration.asArrayDouble("gainfactor");
  os << "gainfactor = " << gainFactorList_ << LogIO::POST;

  // refefence antenna ID
  if (!configuration.isDefined("refant")) {
    os << "ERROR in configuration: reference antenna (refant) must be given."
       << LogIO::EXCEPTION;
  }
  Int const referenceAntenna = configuration.asInt("refant");
  os << "reference antenna ID = " << referenceAntenna << LogIO::POST;

  // raw MS name must be given
  if (!configuration.isDefined("inputms")) {
    os << "ERROR in configuration: raw MS name (inputms) must be given."
       << LogIO::EXCEPTION;

  }
  String const rawMs = configuration.asString("inputms");
  os << "input MS = \"" << rawMs << "\"" << LogIO::POST;

  // read MAIN table (OFF_SOURCE time)
  readMain(rawMs);

  // read SPECTRAL_WINDOW table
  readSpectralWindow(rawMs);

  // read POINTING table
  readPointing(rawMs, referenceAntenna);

  // read ASDM_CALWVR and ASDM_CALATMOSPHERE tables
  readAsdmAsIsTables(rawMs);

  // initialize Atmosphere Transmission Model
  initializeAtmosphereModel(configuration);
}

void SDAtmosphereCorrectionTVI::initializeAtmosphereModel(Record const &configuration) {
  LogIO os(LogOrigin("SDAtmosphereCorrectionTVI", __func__, WHERE));

  // AtmProfile
  // altitude
  if (!configuration.isDefined("siteAltitude")) {
    os << "ERROR in configuration: site altitude (in metre) must be given."
       << LogIO::EXCEPTION;
  }
  atm::Length altitude(configuration.asDouble("siteAltitude"), atm::Length::UnitMeter);
  os << "site altitude = " << altitude.get() << LogIO::POST;

  // pressure (mbar)
  double defaultPressureValue = kValueUnset;
  if (configuration.isDefined("pressure")) {
    userPressureValue_ = configuration.asDouble("pressure");
    defaultPressureValue = userPressureValue_;
    os << "user pressure = " << userPressureValue_ << LogIO::POST;
  } else {
    userPressureValue_ = kValueUnset;
    if (atmPressureData_.nelements() > 0) {
      defaultPressureValue = atmPressureData_[0];
    }
    os << "default pressure value = " << defaultPressureValue << LogIO::POST;
  }

  if (isValueUnset(defaultPressureValue)) {
    os << "ERROR in configuration: No useful pressure value exists."
       << LogIO::EXCEPTION;
  }

  atm::Pressure defaultPressure(defaultPressureValue, atm::Pressure::UnitMilliBar);

  // temperature (K)
  double defaultTemperatureValue = kValueUnset;
  if (configuration.isDefined("temperature")) {
    userTemperatureValue_ = configuration.asDouble("temperature");
    defaultTemperatureValue = userTemperatureValue_;
    os << "user temperature = " << userTemperatureValue_ << LogIO::POST;
  } else {
    userTemperatureValue_ = kValueUnset;
    if (atmTemperatureData_.nelements() > 0) {
      defaultTemperatureValue = atmTemperatureData_[0];
    }
    os << "default temperature value = " << defaultTemperatureValue << LogIO::POST;
  }

  if (isValueUnset(defaultTemperatureValue)) {
    os << "ERROR in configuration: No useful temperature value exists."
       << LogIO::EXCEPTION;
  }

  atm::Temperature defaultTemperature(defaultTemperatureValue, atm::Temperature::UnitKelvin);

  // lapse rate (K/km)
  double tropoLapseRate = (configuration.isDefined("lapseRate")) ?
    configuration.asDouble("lapseRate") : -5.6;

  // relative humidity
  double defaultRelHumidityValue = kValueUnset;
  if (configuration.isDefined("humidity")) {
    userRelHumidityValue_ = configuration.asDouble("humidity");
    defaultRelHumidityValue = userRelHumidityValue_;
    os << "user humidity = " << userRelHumidityValue_ << LogIO::POST;
  } else {
    userRelHumidityValue_ = kValueUnset;
    if (atmRelHumidityData_.nelements() > 0) {
      defaultRelHumidityValue = atmRelHumidityData_[0];
    }
    os << "default humidity value = " << defaultRelHumidityValue << LogIO::POST;
  }

  if (isValueUnset(defaultRelHumidityValue)) {
    os << "ERROR in configuration: No useful humidity value exists."
       << LogIO::EXCEPTION;
  }

  atm::Humidity defaultRelHumidity(defaultRelHumidityValue, atm::Humidity::UnitPercent);

  // PWV (mm)
  double defaultPwvValue = kValueUnset;
  if (configuration.isDefined("pwv")) {
    userPwvValue_ = configuration.asDouble("pwv");
    defaultPwvValue = userPwvValue_;
    os << "user pwv = " << userPwvValue_ << LogIO::POST;
  } else {
    userPwvValue_ = kValueUnset;
    if (pwvData_.nelements() > 0) {
      defaultPwvValue = pwvData_[0];
    }
    os << "default pwv value = " << defaultPwvValue << LogIO::POST;
  }

  if (isValueUnset(defaultPwvValue)) {
    os << "ERROR in configuration: No useful pwv value exists."
       << LogIO::EXCEPTION;
  }

  // scale height (km)
  Double scaleHeightValue = (configuration.isDefined("scaleHeight")) ?
    configuration.asDouble("scaleHeight") : 2.0;
  atm::Length wvScaleHeight(scaleHeightValue, atm::Length::UnitKiloMeter);
  os << "scale height = " << wvScaleHeight.get() << LogIO::POST;

  // pressure step (mbar)
  Double pressureStepValue = (configuration.isDefined("pressureStep")) ?
    configuration.asDouble("pressureStep") : 10.0;
  atm::Pressure pressureStep(pressureStepValue, atm::Pressure::UnitMilliBar);
  os << "pressure step = " << pressureStep.get() << LogIO::POST;

  // pressure step factor
  double pressureStepFactor = (configuration.isDefined("pressureStepFactor")) ?
    configuration.asDouble("pressureStepFactor") : 1.2;
  os << "pressure step factor = " << pressureStepFactor << LogIO::POST;

  // maximum altitude (km)
  double maxAltitude = (configuration.isDefined("maxAltitude")) ?
    configuration.asDouble("maxAltitude") : 120.0;
  atm::Length topAtmProfile(maxAltitude, atm::Length::UnitKiloMeter);
  os << "max ATM altitude = " << topAtmProfile.get() << LogIO::POST;

  // ATM type
  atmType_ = (configuration.isDefined("atmType")) ?
    configuration.asuInt("atmType") : 2;
  os << "ATM type = " << atmType_ << LogIO::POST;

  // layerBoundaries and layerTemperatures
  std::vector<atm::Length> layerBoundaries;
  std::vector<atm::Temperature> layerTemperatures;
  if (configuration.isDefined("layerBoundaries") &&
      configuration.isDefined("layerTemperatures")) {
    Vector<Double> layerBoundariesData = configuration.asArrayDouble("layerBoundaries");
    Vector<Double> layerTemperaturesData = configuration.asArrayDouble("layerTemperatures");
    size_t const numLayer = std::min(layerBoundariesData.nelements(),
                                     layerTemperaturesData.nelements());
    layerBoundaries.resize(numLayer);
    layerTemperatures.resize(numLayer);
    for (size_t i = 0; i < numLayer; ++i) {
      layerBoundaries[i] = atm::Length(
        layerBoundariesData[i],
        atm::Length::UnitKiloMeter
      );
      layerTemperatures[i] = atm::Temperature(
        layerTemperaturesData[i],
        atm::Temperature::UnitKelvin
      );
    }
  }

  // if (!layerBoundaries.empty() && !layerTemperatures.empty()) {
    atmProfile_.reset(new atm::AtmProfile(
      altitude,
      defaultPressure,
      defaultTemperature,
      tropoLapseRate,
      defaultRelHumidity,
      wvScaleHeight,
      pressureStep,
      pressureStepFactor,
      topAtmProfile,
      atmType_,
      layerBoundaries,
      layerTemperatures
    ));
    atmProfile_->setBasicAtmosphericParameterThresholds(
      atm::Length(0.0, atm::Length::UnitMeter),
      atm::Pressure(0.0, atm::Pressure::UnitMilliBar),
      atm::Temperature(0.0, atm::Temperature::UnitKelvin),
      0.0,
      atm::Humidity(0.0, atm::Humidity::UnitPercent),
      atm::Length(0.0, atm::Length::UnitMeter)
    );
  // } else {
  //   atmProfile_.reset(new atm::AtmProfile(
  //     altitude,
  //     defaultPressure,
  //     defaultTemperature,
  //     tropoLapseRate,
  //     defaultRelHumidity,
  //     wvScaleHeight,
  //     pressureStep,
  //     pressureStepFactor,
  //     topAtmProfile,
  //     atmType
  //   ));
  // }

  // SpectralGrid
  for (auto i = processSpwList_.begin(); i != processSpwList_.end(); ++i) {
    SpwId const spw = *i;
    Vector<Double> cf = channelFreqsPerSpw_[spw];
    Vector<Double> cw = channelWidthsPerSpw_[spw];
    // this is the same definition with Python implementation
    unsigned int nchan = cf.nelements();
    unsigned int refChan = (nchan - 1) / 2;
    double centerFreq = (nchan % 2 == 1) ? cf[nchan / 2] : 0.5 * (cf[nchan / 2 - 1] + cf[nchan / 2]);
    Double chanSep = (nchan == 1u) ? cw[0] : (cf[nchan - 1] - cf[0]) / static_cast<Double>(nchan - 1);
    if (isTdmSpw_[spw]) {
      // configure 5x finer spectral grid than native one
      chanSep /= 5.0;
      nchan *= 5u;
    }
    os.output() << std::setprecision(16);
    os << "SpectralGrid for spw " << spw << ": nchan " << nchan
       << " refchan " << refChan << " center freq "
       << centerFreq << " chan sep " << chanSep << LogIO::POST;
    atmSpectralGridPerSpw_[spw].reset(
      new atm::SpectralGrid(nchan, refChan, centerFreq, chanSep)
    );
  }

  // SkyStatus
  resetSkyStatus(*atmProfile_.get(), defaultPwvValue);
  os << "DONE Initializing SkyStatus" << LogIO::POST;
}

void SDAtmosphereCorrectionTVI::configureAtmosphereCorrection() {
  // TODO: implement all the configuration steps here
  // current spw
  LogIO os(LogOrigin("SDAtmosphereCorrectionTVI", __func__, WHERE));
  SpwId const currentSpw = dataDescriptionSubtablecols().spectralWindowId().get(dataDescriptionId());
  bool isProcessingSpw = anyEQ(processSpwList_, currentSpw);
  Vector<Double> currentTime;
  time(currentTime);
  Vector<Int> currentStateId;
  stateId(currentStateId);
  bool isOnSource = true; //isOnSourceChunk();
  bool isPrecedingOffSourceScanExist = min(offSourceTime_) <= currentTime[0];
  bool isSubsequentOffSourceScanExist = currentTime[0] <= max(offSourceTime_);
  bool isPrecedingAtmScanExist = min(atmTime_) <= currentTime[0];
  cout << std::setprecision(16) << offSourceTime_ << endl;
  cout << std::setprecision(16) << currentTime[0] << " state ID" << currentStateId[0] << endl;
  os << "SPW " << currentSpw << ": processingSpw " << isProcessingSpw
     << " OFF_SOURCE availability before " << isPrecedingOffSourceScanExist
     << " after " << isSubsequentOffSourceScanExist
     << " ON_SOURCE? " << isOnSource << LogIO::POST;
  atmSkyStatusPtr_ = nullptr;
  if (isProcessingSpw && isOnSource &&
      isPrecedingOffSourceScanExist && isSubsequentOffSourceScanExist &&
      isPrecedingAtmScanExist) {
    // gain factor for current spw
    if (currentSpw < 0 || static_cast<SpwId>(gainFactorList_.nelements()) <= currentSpw) {
      gainFactor_ = 1.0;
    } else {
      gainFactor_ = gainFactorList_[currentSpw];
    }
    os << "gainfactor for SPW " << currentSpw << " = " << gainFactor_ << LogIO::POST;

    // SkyStatus for current spw
    auto finder = atmSkyStatusPerSpw_.find(currentSpw);
    if (finder != atmSkyStatusPerSpw_.end()) {
      atmSkyStatusPtr_ = finder->second.get();
      atmSpectralGridPtr_ = atmSpectralGridPerSpw_[currentSpw].get();
    }
  }
}

void SDAtmosphereCorrectionTVI::updateSkyStatus() {
  // do nothing if nullptr is given
  if (!atmSkyStatusPtr_) {
    return;
  }

  atm::SkyStatus const * const p = atmSkyStatusPtr_;
  Vector<Double> currentTime;
  time(currentTime);
  std::pair<Int, Int> pair = findNearestIndex(atmTime_, currentTime[0]);
  Int timeIndex = pair.first;
  bool isSkyStatusOutdated = false;
  if (isValueUnset(userTemperatureValue_)) {
    Double val = atmTemperatureData_[timeIndex];
    Double val2 = p->getGroundTemperature().get(atm::Temperature::UnitKelvin);
    if (val != val2) {
      cout << "Temperature is different " << val << " vs " << val2 << endl;
      isSkyStatusOutdated = true;
    }
    atmSkyStatusPtr_->setBasicAtmosphericParameters(atm::Temperature(val, atm::Temperature::UnitKelvin));
  }
  cout << "time " << std::setprecision(16) << currentTime[0]
       << " atmTime " << atmTime_[timeIndex]
       << " atmIndex " << timeIndex << endl;

  if (isValueUnset(userPressureValue_)) {
    Double val = atmPressureData_[timeIndex];
    Double val2 = p->getGroundPressure().get(atm::Pressure::UnitMilliBar);
    if (val != val2) {
      cout << "Pressure is different " << val << " vs " << val2 << endl;
      isSkyStatusOutdated = true;
    }
    atmSkyStatusPtr_->setBasicAtmosphericParameters(atm::Pressure(val, atm::Pressure::UnitMilliBar));
  }

  if (isValueUnset(userRelHumidityValue_)) {
    Double val = atmRelHumidityData_[timeIndex];
    Double val2 = p->getRelativeHumidity().get(atm::Humidity::UnitPercent);
    if (val != val2) {
      cout << "Humidity is different " << val << " vs " << val2 << endl;
      isSkyStatusOutdated = true;
    }
    atmSkyStatusPtr_->setBasicAtmosphericParameters(atm::Humidity(val, atm::Humidity::UnitPercent));
  }

  // share time information with CalAtmosphere table
  Double pwvValue = userPwvValue_;
  if (isValueUnset(userPwvValue_)) {
    Double val = pwvData_[timeIndex];
    Double val2 = p->getUserWH2O().get(atm::Length::UnitMilliMeter);
    if (val != val2) {
      cout << "PWV is different " << val << " vs " << val2 << endl;
      isSkyStatusOutdated = true;
    }
    pwvValue = val;
    atmSkyStatusPtr_->setUserWH2O(atm::Length(pwvValue, atm::Length::UnitMilliMeter));
  }

  isSkyStatusOutdated = false;
  if (isSkyStatusOutdated) {
    // reset SkyStatus
    atm::AtmProfile atmProfile(
      p->getAltitude(),
      p->getGroundPressure(),
      p->getGroundTemperature(),
      // atm::Pressure(atmPressureData_[timeIndex], atm::Pressure::UnitMilliBar),
      // atm::Temperature(atmTemperatureData_[timeIndex], atm::Temperature::UnitKelvin),
      p->getTropoLapseRate(),
      p->getRelativeHumidity(),
      // atm::Humidity(atmRelHumidityData_[timeIndex], atm::Humidity::UnitPercent),
      p->getWvScaleHeight(),
      p->getPressureStep(),
      p->getPressureStepFactor().get(),
      p->getTopAtmProfile(),
      atmType_,
      p->getThicknessProfile(),
      p->getTemperatureProfile()
    );
    SpwId const currentSpw = dataDescriptionSubtablecols().spectralWindowId().get(dataDescriptionId());
    // atm::RefractiveIndexProfile profile(
    //   *atmSpectralGridPerSpw_[currentSpw].get(),
    //   atmProfile
    // );
    // atmSkyStatusPerSpw_[currentSpw].reset(
    //   new atm::SkyStatus(profile)
    // );
    pwvValue = p->getUserWH2O().get(atm::Length::UnitMilliMeter);
    resetSkyStatus(atmProfile, pwvValue);
    atmSkyStatusPtr_ = atmSkyStatusPerSpw_[currentSpw].get();
  }

  // pair = findNearestIndex(pwvTime_, currentTime[0]);
  // timeIndex = pair.first;
  cout << "time " << currentTime[0]
       << " Temperature " << atmTemperatureData_[timeIndex] << "K "
       << " Pressure " << atmPressureData_[timeIndex] << "mbar "
       << " Humidity " << atmRelHumidityData_[timeIndex] << "%" << endl;
}

void SDAtmosphereCorrectionTVI::resetSkyStatus(atm::AtmProfile const &atmProfile, double const pwv) {
  for (auto i = processSpwList_.begin(); i != processSpwList_.end(); ++i) {
    SpwId const spw = *i;
    atm::RefractiveIndexProfile profile(*atmSpectralGridPerSpw_[spw].get(), atmProfile);
    atmSkyStatusPerSpw_[spw].reset(new atm::SkyStatus(profile));
    atmSkyStatusPerSpw_[spw]->setUserWH2O(atm::Length(pwv, atm::Length::UnitMilliMeter));
  }
}

void SDAtmosphereCorrectionTVI::updateCorrectionFactor() {
  // discard correction factor if nullptr is given
  if (!atmSkyStatusPtr_) {
    correctionFactor_.resize();
    return;
  }

  // PRECONDITION:
  // all rows have the same time stamp (should be true for no RowBlocking mode)

  // find OFF_SOURCE time stamps that brackets current time stamp
  Vector<Double> currentTimeList;
  time(currentTimeList);
  Double const currentTime = currentTimeList[0];
  cout << "offSourceTime_ (size " << offSourceTime_.nelements() << ") = " << offSourceTime_ << endl;
  cout << "currentTime = " << currentTime << endl;
  std::pair<Double, Double> nearest = findNearest(offSourceTime_, currentTime);
  Double const offSourceTimePrev = nearest.first;
  Double const offSourceTimeNext = nearest.second;
  cout << "nearest: " << std::setprecision(16) << offSourceTimePrev << "~" << offSourceTimeNext << endl;

  // elevation
  cout << "interpolated elevation" << endl;
  Double const elevationOffPrev = elevationInterpolator_(offSourceTimePrev);//interpolateDataLinear(elevationTime_, elevationData_, offSourceTimePrev);
  cout << "prev: " << elevationOffPrev << endl;
  Double const elevationOn = elevationInterpolator_(currentTime);//interpolateDataLinear(elevationTime_, elevationData_, currentTime);
  cout << "ON: " << elevationOn << endl;
  Double const elevationOffNext = elevationInterpolator_(offSourceTimeNext);//interpolateDataLinear(elevationTime_, elevationData_, offSourceTimeNext);
  cout << "next: " << elevationOffNext << endl;
  Double const elevationOff =
    ((offSourceTimeNext - currentTime) * elevationOffPrev +
     (currentTime - offSourceTimePrev) * elevationOffNext) /
    (offSourceTimeNext - offSourceTimePrev);
  cout << "OFF: " << elevationOffPrev << endl;
  cout << "time " << std::setprecision(16) << currentTime << " elON " << elevationOn << " elOFF " << elevationOff << endl;

  // opacity
  unsigned int const nchan = atmSpectralGridPtr_->getNumChan();
  Double const airMassOn = 1.0 / cos(C::pi_2 - elevationOn);
  atmSkyStatusPtr_->setAirMass(airMassOn);
  Vector<Double> trjSkySpecOn = getTrjSkySpec(atmSkyStatusPtr_, nchan);
  Double const airMassOff = 1.0 / cos(C::pi_2 - elevationOff);
  atmSkyStatusPtr_->setAirMass(airMassOff);
  Vector<Double> trjSkySpecOff = getTrjSkySpec(atmSkyStatusPtr_, nchan);
  Vector<Double> tauSpecOff = getTauSpec(atmSkyStatusPtr_, nchan, airMassOff);
  Vector<Double> correctionFactor = getCorrectionFactor(trjSkySpecOn, trjSkySpecOff, tauSpecOff);

  constexpr size_t ic = 0;
  cout << std::setprecision(16)
       << "time = " << currentTime
       << " TrjON = " << trjSkySpecOn[ic]
       << " TrjOFF = " << trjSkySpecOff[ic]
       << " tauOFF = " << tauSpecOff[ic]
       << " factor = " << correctionFactor[ic] << " (without gain factor)"
       << endl;

  // current spw
  SpwId const currentSpw = dataDescriptionSubtablecols().spectralWindowId().get(dataDescriptionId());
  if (isTdmSpw_[currentSpw]) {
    // TODO: convolve with hanning window
    // TODO: interpolation to native freq
    unsigned int numAtmChan = atmSpectralGridPtr_->getNumChan();
    unsigned int refAtmChan = atmSpectralGridPtr_->getRefChan();
    double refAtmFreq = atmSpectralGridPtr_->getRefFreq().get();
    double sepAtmFreq = atmSpectralGridPtr_->getChanSep().get();
    Vector<Double> atmFreq(numAtmChan);
    indgen(atmFreq, refAtmFreq - sepAtmFreq * refAtmChan, sepAtmFreq);
    Vector<Double> const &nativeFreq = channelFreqsPerSpw_[currentSpw];
    correctionFactor_.resize(nativeFreq.nelements());
    Interpolate1D<Double, Double> interpolator(
      ScalarSampledFunctional<Double>(atmFreq),
      ScalarSampledFunctional<Double>(correctionFactor), True, True);
    for (uInt i = 0; i < nativeFreq.nelements(); ++i) {
      correctionFactor_[i] = interpolator(nativeFreq[i]);
    }
  } else if (doSmooth_[currentSpw]) {
    // TODO: smoothing
  } else {
    correctionFactor_.reference(correctionFactor);
  }

  // apply gain factor
  cout << "Applying gain factor " << gainFactorList_[currentSpw] << " to correction term" << endl;
  correctionFactor_ *= gainFactorList_[currentSpw];
}

void SDAtmosphereCorrectionTVI::updateAtmosphereModel() {
  if (atmSkyStatusPtr_) {
    updateSkyStatus();
    updateCorrectionFactor();
  }
}

void SDAtmosphereCorrectionTVI::warnIfNoTransform() {
  if (!doTransform()) {
    Int const spw = dataDescriptionSubtablecols().spectralWindowId().get(dataDescriptionId());
    LogIO os(LogOrigin("SDAtmosphereCorrectionTVI", __func__, WHERE));
    os << LogIO::WARN << "Spectral Window " << spw
       << " (data description " << dataDescriptionId()
       << ") is output but not corrected." << LogIO::POST;
  }
}

void SDAtmosphereCorrectionTVI::readMain(String const &msName) {
  MeasurementSet msObj(msName);
  Float const noCache = -1;
  MSMetaData meta(&msObj, noCache);
  std::set<Double> timeOffSourceSet = meta.getTimesForIntent("OBSERVE_TARGET#OFF_SOURCE");
  std::vector<Double> timeOffSource(timeOffSourceSet.begin(), timeOffSourceSet.end());
  size_t numBoundaries = 0;
  Vector<uInt> idx1(timeOffSource.size());
  Double cache = timeOffSource[0];
  uInt ix = 0;
  idx1[numBoundaries++] = 0;
  for (uInt i = 0; i < timeOffSource.size(); ++i) {
    Double diff = timeOffSource[i] - cache;
    // threshold (1s) is taken from original script by T. Sawada
    if (diff > 1.0) {
      cout << "boundaries = " << ix << endl;
      idx1[numBoundaries++] = i;
    }
    cache = timeOffSource[i];
  }
  Vector<uInt> idx2(numBoundaries);
  for (uInt i = 0; i < numBoundaries - 1; ++i) {
    idx2[i] = idx1[i + 1] - 1;
  }
  idx2[numBoundaries - 1] = timeOffSource.size() - 1;
  cout << "idx1 = " << idx1 << endl;
  cout << "idx2 = " << idx2 << endl;
  offSourceTime_.resize(numBoundaries);
  for (size_t i = 0; i < numBoundaries; ++i) {
    offSourceTime_[i] = (timeOffSource[idx1[i]] + timeOffSource[idx2[i]]) / 2.0;
  }
  Vector<Double> tmp(numBoundaries);
  for (auto i = 0u; i < numBoundaries; ++i) {
    tmp[i] = offSourceTime_[i] - offSourceTime_[0];
  }
  cout << "diff = " << tmp << endl;
}

void SDAtmosphereCorrectionTVI::readSpectralWindow(String const &msName) {
  MeasurementSet const ms(msName);

  Table const spwTable = ms.spectralWindow();
  ArrayColumn<Double> chanFreqColumn(spwTable, "CHAN_FREQ");
  ArrayColumn<Double> chanWidthColumn(spwTable, "CHAN_WIDTH");
  for (auto i = processSpwList_.begin(); i != processSpwList_.end(); ++i) {
    SpwId const spw = *i;
    channelFreqsPerSpw_[spw] = chanFreqColumn(spw);
    channelWidthsPerSpw_[spw] = chanWidthColumn(spw);
    size_t const nchan = channelFreqsPerSpw_[spw].nelements();
    isTdmSpw_[spw] = (nchan == 128u || nchan == 256u);
  }
  LogIO os(LogOrigin("SDAtmosphereCorrectionTVI", __func__, WHERE));
  os << "TDM flag:" << endl;
  for (auto i = isTdmSpw_.begin(); i != isTdmSpw_.end(); ++i) {
    os << "    spw " << i->first << " " << i->second << endl;
  }
  os << LogIO::POST;
}

void SDAtmosphereCorrectionTVI::readPointing(String const &msName, Int const referenceAntenna) {
  MeasurementSet const ms(msName);
  uInt numAntennas = ms.antenna().nrow();
  if (referenceAntenna < 0 || numAntennas <= static_cast<uInt>(referenceAntenna)) {
    LogIO os(LogOrigin("SDAtmosphereCorrectionTVI", __func__, WHERE));
    os << "ERROR: reference antenna " << referenceAntenna << " is invalid." << LogIO::EXCEPTION;
  }

  LogIO os(LogOrigin("SDAtmosphereCorrectionTVI", __func__, WHERE));
  os << "readPointing" << LogIO::POST;
  Table const pointingTable = ms.pointing();
  // TODO: abort if direction reference is not AZEL
  Table selected = pointingTable(pointingTable.col("ANTENNA_ID") == referenceAntenna);
  elevationTime_ = ScalarColumn<Double>(selected, "TIME").getColumn();
  Array<Double> directionArray = ArrayColumn<Double>(selected, "DIRECTION").getColumn();
  os << "directionArray shape = " << directionArray.shape() << LogIO::POST;
  Array<Double> elevationArray = Cube<Double>(directionArray).yzPlane(1);
  //(IPosition(3, 1, 0, 0), IPosition(3, 1, 0, elevationTime_.nelements() - 1));
  os << "elevationArray shape = " << elevationArray.shape() << LogIO::POST;
  elevationArray.removeDegenerate();
  os << "elevationArray remove degenerate axes = " << elevationArray.shape() << LogIO::POST;
  elevationData_.reference(elevationArray);
  os << "elevationData_.shape = " << elevationData_.shape() << LogIO::POST;
  elevationInterpolator_ = Interpolate1D<Double, Double>(
    ScalarSampledFunctional<Double>(elevationTime_),
    ScalarSampledFunctional<Double>(elevationData_), True, True);
  elevationInterpolator_.setMethod(Interpolate1D<Double, Double>::linear);
  os.output() << std::setprecision(16);
  // os << "elevationTime_ = " << elevationTime_ << LogIO::POST;
  // os << "elevationData_ = " << elevationData_ << LogIO::POST;
}

void SDAtmosphereCorrectionTVI::readAsdmAsIsTables(String const &msName) {
  Vector<Double> pwvTimeLocal;
  Vector<Double> atmTimeLocal;
  Vector<Double> pwvDataLocal;
  Vector<Double> atmTemperatureDataLocal;
  Vector<Double> atmPressureDataLocal;
  Vector<Double> atmRelHumidityDataLocal;

  File const calWvrTableName(msName + "/ASDM_CALWVR");
  LogIO os(LogOrigin("SDAtmosphereCorrectionTVI", __func__, WHERE));
  os << "CALWVR table \"" << calWvrTableName.path().absoluteName() << "\"" << LogIO::POST;
  if (calWvrTableName.exists()) {
    Table const calWvrTable(calWvrTableName.path().absoluteName());
    pwvTimeLocal.reference(ScalarColumn<Double>(calWvrTable, "startValidTime").getColumn());
    pwvDataLocal.reference(ScalarColumn<Double>(calWvrTable, "water").getColumn());
    os << "pwvTimeLocal size " << pwvTimeLocal.nelements() << LogIO::POST;
  }

  File const calAtmosphereTableName(msName + "/ASDM_CALATMOSPHERE");
  os << "CALATM table \"" << calAtmosphereTableName.path().absoluteName() << "\"" << LogIO::POST;
  if (calAtmosphereTableName.exists()) {
    Table const calAtmosphereTable(calAtmosphereTableName.path().absoluteName());
    atmTimeLocal.reference(ScalarColumn<Double>(calAtmosphereTable, "startValidTime").getColumn());
    atmTemperatureDataLocal.reference(ScalarColumn<Double>(calAtmosphereTable, "groundTemperature").getColumn());
    atmPressureDataLocal.reference(ScalarColumn<Double>(calAtmosphereTable, "groundPressure").getColumn());
    atmRelHumidityDataLocal.reference(ScalarColumn<Double>(calAtmosphereTable, "groundRelHumidity").getColumn());
    os << "atmTimeLocal size " << atmTimeLocal.nelements() << LogIO::POST;
  }

  Vector<uInt> uniqueVectorPwv;
  Vector<uInt> indexVectorPwv;
  if (pwvTimeLocal.nelements() > 0) {
    sortTime(pwvTimeLocal, indexVectorPwv);
    uInt n = makeUnique(pwvTimeLocal, indexVectorPwv, uniqueVectorPwv);
    os << "makeUnique returned " << n << LogIO::POST;
    os << "indexVector = " << indexVectorPwv << LogIO::POST;
    pwvTime_.resize(n);
    for (uInt i = 0; i < n; ++i) {
      os << "loop " << i << " index " << uniqueVectorPwv[i] << LogIO::POST;
      os << "loop " << i << " index2 " << indexVectorPwv[uniqueVectorPwv[i]] << LogIO::POST;
      pwvTime_[i] = pwvTimeLocal[indexVectorPwv[uniqueVectorPwv[i]]];
    }
    os << "uniqueVector (shape " << uniqueVectorPwv.shape() << ") " << uniqueVectorPwv << LogIO::POST;
    os << "indexVector (shape " << indexVectorPwv.shape() << ") " << indexVectorPwv << LogIO::POST;
    os.output() << std::setprecision(16);
    os << "pwvTime_ = " << pwvTime_ << LogIO::POST;
    // m -> mm
    os << "getMedianDataPerTime for pwv start" << LogIO::POST;
    // pwvData_ = getMedianDataPerTime(n, uniqueVectorPwv, indexVectorPwv, pwvDataLocal) * 1000.0;
    // os << "pwvData_ = " << pwvData_ << LogIO::POST;
    os << "getMedianDataPerTime for pwv end" << LogIO::POST;
  }

  if (atmTimeLocal.nelements() > 0) {
    Vector<uInt> uniqueVector;
    Vector<uInt> indexVector;
    sortTime(atmTimeLocal, indexVector);
    uInt n = makeUnique(atmTimeLocal, indexVector, uniqueVector);
    atmTime_.resize(n);
    for (uInt i = 0; i < n; ++i) {
      atmTime_[i] = atmTimeLocal[indexVector[uniqueVector[i]]];
    }
    os.output() << std::setprecision(16) << std::endl;
    os << "atmTime_ = " << atmTime_ << LogIO::POST;
    os << "getMedianDataPerTime for temp start" << LogIO::POST;
    atmTemperatureData_ = getMedianDataPerTime(n, uniqueVector, indexVector, atmTemperatureDataLocal);
    os << "Temperature = " << atmTemperatureData_ << LogIO::POST;
    os << "getMedianDataPerTime for temp end" << LogIO::POST;
    // Pa -> mbar (=hPa)
    atmPressureData_ = getMedianDataPerTime(n, uniqueVector, indexVector, atmPressureDataLocal) / 100.0;
    os << "Pressure = " << atmPressureData_ << LogIO::POST;
    atmRelHumidityData_ = getMedianDataPerTime(n, uniqueVector, indexVector, atmRelHumidityDataLocal);
    os << "Humidity = " << atmRelHumidityData_ << LogIO::POST;

    // same PWV handling with original script
    Vector<uInt> uniqueVectorPwv2(n);
    for (uInt i = 0; i < n; ++i) {
      Vector<Double> diff = abs(pwvTime_ - atmTime_[i]);
      Double minDiff = min(diff);
      uInt const m = pwvTime_.nelements();
      uInt idx = m;
      for (uInt j = 0; j < m; ++j) {
        if (diff[j] == minDiff) {
          idx = j;
          break;
        }
      }
      if (idx < m) {
        uniqueVectorPwv2[i] = uniqueVectorPwv[idx];
      } else {
        throw AipsError("ERROR in readAsdmAsIsTable");
      }
    }
    pwvData_ = getMedianDataPerTime(n, uniqueVectorPwv2, indexVectorPwv, pwvDataLocal) * 1000.0;
    os << "PWV = " << pwvData_ << LogIO::POST;
  }
  os << "DONE readAsdmAsIsTable" << LogIO::POST;
}

//////////
// SDAtmosphereCorrectionTVIFactory
/////////
SDAtmosphereCorrectionVi2Factory::SDAtmosphereCorrectionVi2Factory(Record const &configuration,
    ViImplementation2 *inputVII) :
    ViFactory(),
    inputVII_p(inputVII),
    configuration_p(configuration) {
}

SDAtmosphereCorrectionVi2Factory::SDAtmosphereCorrectionVi2Factory(Record const &configuration,
    MeasurementSet const *ms, SortColumns const sortColumns,
    Double timeInterval, Bool isWritable) :
    ViFactory(),
    inputVII_p(nullptr),
    configuration_p(configuration) {
  inputVII_p = new VisibilityIteratorImpl2(Block<MeasurementSet const *>(1, ms),
      sortColumns, timeInterval, isWritable);
}

SDAtmosphereCorrectionVi2Factory::~SDAtmosphereCorrectionVi2Factory() {
}

ViImplementation2 * SDAtmosphereCorrectionVi2Factory::createVi() const {
  if (inputVII_p->getNMs() != 1) {
    throw AipsError("ERROR: Multiple or zero MS inputs are given. SDAtmosphereCorrectionTVI works with only one MS.");
  }
  return new SDAtmosphereCorrectionTVI(inputVII_p, configuration_p);
}

SDAtmosphereCorrectionTVILayerFactory::SDAtmosphereCorrectionTVILayerFactory(
    Record const &configuration) :
    ViiLayerFactory(),
    configuration_p(configuration) {
}

ViImplementation2*
SDAtmosphereCorrectionTVILayerFactory::createInstance(ViImplementation2* vii0) const {
  // Make the SDAtmosphereCorrectionTVI, using supplied ViImplementation2, and return it
  SDAtmosphereCorrectionVi2Factory factory(configuration_p, vii0);
  ViImplementation2 *vii = nullptr;
  try {
    vii = factory.createVi();
  } catch (...) {
    if (vii0) {
      delete vii0;
    }
    throw;
  }
  return vii;
}
} // # NAMESPACE VI - END
} // #NAMESPACE CASA - END
