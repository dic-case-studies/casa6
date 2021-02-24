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

inline uInt sortTime(Vector<Double> &timeData, Vector<uInt> &uniqueVector, Vector<uInt> indexVector) {
  Bool b = false;
  Double *p = timeData.getStorage(b);
  ScopeGuard guard([&]() {
    timeData.putStorage(p, b);
  });
  Sort sort(p, timeData.nelements());
  return sort.unique(uniqueVector, indexVector);
}

inline Vector<Double> getMedianDataPerTime(uInt const n, Vector<uInt> const &uniqueVector,
  Vector<uInt> const &indexVector, Vector<Double> &data) {
  data.resize(n);
  assert(n > 0);
  for (uInt i = 0; i < n - 1; ++i) {
    uInt const iStart = uniqueVector[i];
    uInt const iEnd = uniqueVector[i + 1];
    uInt const nElem = iEnd - iStart;
    Vector<Double> tmp(nElem);
    for (uInt j = 0; j < nElem; ++j) {
      tmp[j] = data[indexVector[j + iStart]];
    }
    data[i] = median(tmp);
  }
  uInt const iStart = uniqueVector[n - 1];
  uInt const iEnd = indexVector.nelements();
  uInt const nElem = iEnd - iStart;
  Vector<Double> tmp(nElem);
  for (uInt j = 0; j < nElem; ++j) {
    tmp[j] = data[indexVector[j + iStart]];
  }
  data[n - 1] = median(tmp);
  return data;
}

inline std::pair<Double, Double> findNearest(Vector<Double> const &data, Double const target) {
  Bool found = false;
  // data is sorted in ascending order
  Int index = binarySearch(found, data, target, data.nelements());
  if (!found) {
    throw AipsError("error in findNearest");
  }
  Double prev = kValueUnset;
  Double next = kValueUnset;
  if (data[index] == target) {
    prev = target;
    next = target;
  } else {
    if (index == 0 || index == static_cast<Int>(data.nelements() - 1)) {
      throw AipsError("error in findNearest");
    }
    prev = data[index - 1];
    next = data[index];
  }
  return std::make_pair(prev, next);
}

inline Double interpolateDataLinear(Vector<Double> const &xin, Vector<Double> const &yin, Double const xout) {
  Bool found = false;
  // data is sorted in ascending order
  Int index = binarySearch(found, xin, xout, xin.nelements());
  if (!found) {
    throw AipsError("error in interpolateData");
  }
  Double yout = kValueUnset;
  if (index == 0 || xin[index] == xout) {
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

inline Vector<Double> getTauSpec(atm::SkyStatus *skyStatus, unsigned int const nchan) {
  // TODO: OMP parallelization
  Vector<Double> tau(nchan);
  for(unsigned int i = 0; i < nchan; ++i) {
    tau[i] = skyStatus->getDryOpacity(i).get(atm::Opacity::UnitNeper) +
             skyStatus->getWetOpacity(i).get(atm::Opacity::UnitNeper);
  }
  return tau;
}

inline Vector<Double> getCorrectionFactor(
  Vector<Double> const &tSkyOn, Vector<Double> const &tSkyOff, Vector<Double> const &tauOff) {
  Vector<Double> factor;
  // TODO: OMP/SIMD parallelization
  for (size_t i = 0; i < tSkyOn.nelements(); ++i) {
    factor[i] = (tSkyOn[i] - tSkyOff[i]) * std::exp(tauOff[i]);
  }
  return factor;
}

template<class T>
inline void transformData(Vector<Double> const gainFactor, Cube<T> const &in, Cube<T> &out) {
  // TODO: OMP/SIMD parallelization
  if (allEQ(gainFactor, 1.0)) {
    out.reference(in);
  } else {
    IPosition const shape = in.shape();
    out.resize(shape);
    for (ssize_t ir = 0; ir < shape[2]; ++ir) {
      for (ssize_t ic = 0; ic < shape[1]; ++ic) {
        for (ssize_t ip = 0; ip < shape[0]; ++ip) {
          out(ip, ic, ir) = gainFactor[ic] * in(ip, ic, ir);
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
    //elevationInterpolator_(),
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
    atmProfile_(),
    atmSpectralGridPerSpw_(),
    atmSkyStatusPerSpw_(),
    atmSkyStatusPtr_(nullptr),
    correctionFactor_() {
  initializeAtmosphereCorrection(configuration);

  // Initialize attached VisBuffer
  setVisBuffer(createAttachedVisBuffer(VbRekeyable));
}

SDAtmosphereCorrectionTVI::~SDAtmosphereCorrectionTVI() {
}

void SDAtmosphereCorrectionTVI::origin() {
  TransformingVi2::origin();

  // Synchronize own VisBuffer
  configureNewSubchunk();

  // configure atmospheric correction
  configureAtmosphereCorrection();

  // warn if current spw is not requested to transform
  warnIfNoTransform();
}

void SDAtmosphereCorrectionTVI::next() {
  TransformingVi2::next();

  // Synchronize own VisBuffer
  configureNewSubchunk();

  // configure atmospheric correction if necessary
  configureAtmosphereCorrection();

  // warn if current spw is not requested to transform
  warnIfNoTransform();
}

void SDAtmosphereCorrectionTVI::originChunks(Bool forceRewind) {
  TransformingVi2::originChunks();

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

  // gain factor
  if (!configuration.isDefined("gainfactor")) {
    os << "ERROR in configuration: gain factor list (gainfactor) must be given."
       << LogIO::EXCEPTION;
  }
  gainFactorList_ = configuration.asArrayDouble("gainfactor");

  // refefence antenna ID
  if (!configuration.isDefined("refant")) {
    os << "ERROR in configuration: reference antenna (refant) must be given."
       << LogIO::EXCEPTION;
  }
  Int const referenceAntenna = configuration.asInt("refant");

  // raw MS name must be given
  if (!configuration.isDefined("inputms")) {
    os << "ERROR in configuration: raw MS name (inputms) must be given."
       << LogIO::EXCEPTION;

  }
  String const rawMs = configuration.asString("inputms");

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

  // pressure (mbar)
  double defaultPressureValue = kValueUnset;
  if (configuration.isDefined("pressure")) {
    userPressureValue_ = configuration.asDouble("pressure");
    defaultPressureValue = userPressureValue_;
  } else {
    userPressureValue_ = kValueUnset;
    if (atmPressureData_.nelements() > 0) {
      defaultPressureValue = atmPressureData_[0];
    }
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
  } else {
    userTemperatureValue_ = kValueUnset;
    if (atmTemperatureData_.nelements() > 0) {
      defaultTemperatureValue = atmTemperatureData_[0];
    }
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
  } else {
    userRelHumidityValue_ = kValueUnset;
    if (atmRelHumidityData_.nelements() > 0) {
      defaultRelHumidityValue = atmRelHumidityData_[0];
    }
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
  } else {
    userPwvValue_ = kValueUnset;
    if (pwvData_.nelements() > 0) {
      defaultPwvValue = pwvData_[0];
    }
  }

  if (isValueUnset(defaultPwvValue)) {
    os << "ERROR in configuration: No useful pwv value exists."
       << LogIO::EXCEPTION;
  }

  // scale height (km)
  Double scaleHeightValue = (configuration.isDefined("scaleHeight")) ?
    configuration.asDouble("scaleHeight") : 2.0;
  atm::Length wvScaleHeight(scaleHeightValue, atm::Length::UnitKiloMeter);

  // pressure step (mbar)
  Double pressureStepValue = (configuration.isDefined("pressureStep")) ?
    configuration.asDouble("pressureStep") : 10.0;
  atm::Pressure pressureStep(pressureStepValue, atm::Pressure::UnitMilliBar);

  // pressure step factor
  double pressureStepFactor = (configuration.isDefined("pressureStepFactor")) ?
    configuration.asDouble("pressureStepFactor") : 1.2;

  // maximum altitude (km)
  double maxAltitude = (configuration.isDefined("maxAltitude")) ?
    configuration.asDouble("maxAltitude") : 120.0;
  atm::Length topAtmProfile(maxAltitude, atm::Length::UnitKiloMeter);

  // ATM type
  unsigned int atmType = (configuration.isDefined("atmType")) ?
    configuration.asuInt("atmType") : 2;

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

  if (!layerBoundaries.empty() && !layerTemperatures.empty()) {
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
      atmType,
      layerBoundaries,
      layerTemperatures
    ));
  } else {
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
      atmType
    ));
  }

  // SpectralGrid
  for (auto i = processSpwList_.begin(); i != processSpwList_.end(); ++i) {
    SpwId const spw = *i;
    Vector<Double> cf = channelFreqsPerSpw_[spw];
    Vector<Double> cw = channelWidthsPerSpw_[spw];
    size_t nchan = cf.nelements();
    Double chansep = (nchan == 1u) ? cw[0] : (cf[nchan - 1] - cf[0]) / static_cast<Double>(nchan - 1);
    if (isTdmSpw_[spw]) {
      // configure 5x finer spectral grid than native one
      chansep /= 5.0;
      nchan *= 5u;
    }
    atmSpectralGridPerSpw_[spw].reset(
      new atm::SpectralGrid(nchan, 0u, cf[0], chansep)
    );
  }

  // SkyStatus
  for (auto i = processSpwList_.begin(); i != processSpwList_.end(); ++i) {
    SpwId const spw = *i;
    atm::RefractiveIndexProfile profile(*atmSpectralGridPerSpw_[spw].get(), *atmProfile_.get());
    atm::SkyStatus *skyStatus = new atm::SkyStatus(profile);
    skyStatus->setUserWH2O(atm::Length(defaultPwvValue, atm::Length::UnitMilliMeter));
    atmSkyStatusPerSpw_[spw].reset(skyStatus);
  }
}

void SDAtmosphereCorrectionTVI::configureAtmosphereCorrection() {
  // TODO: implement all the configuration steps here
  // current spw
  SpwId const currentSpw = dataDescriptionSubtablecols().spectralWindowId().get(dataDescriptionId());
  bool isProcessingSpw = anyEQ(processSpwList_, currentSpw);
  Vector<Double> currentTime;
  time(currentTime);
  bool isPrecedingOffSourceScanExist = min(offSourceTime_) <= currentTime[0];
  bool isSubsequentOffSourceScanExist = currentTime[0] <= max(offSourceTime_);
  atmSkyStatusPtr_ = nullptr;
  if (isProcessingSpw && isPrecedingOffSourceScanExist && isSubsequentOffSourceScanExist) {
    // gain factor for current spw
    if (currentSpw < 0 || static_cast<SpwId>(gainFactorList_.nelements()) <= currentSpw) {
      gainFactor_ = 1.0;
    } else {
      gainFactor_ = gainFactorList_[currentSpw];
    }

    // SkyStatus for current spw
    auto finder = atmSkyStatusPerSpw_.find(currentSpw);
    if (finder != atmSkyStatusPerSpw_.end()) {
      atmSkyStatusPtr_ = finder->second.get();
      atmSpectralGridPtr_ = atmSpectralGridPerSpw_[currentSpw].get();
    }
  }

  updateAtmosphereModel();
}

void SDAtmosphereCorrectionTVI::updateSkyStatus(atm::SkyStatus *p) {
  // do nothing if nullptr is given
  if (!p) {
    return;
  }

  // TODO: configure appropriate index from time stamp
  // TODO: if no valid time stamp exists, make
  size_t timeIndex = 0;
  if (isValueUnset(userTemperatureValue_)) {
    Double temperatureValue = atmTemperatureData_[timeIndex];
    p->setBasicAtmosphericParameters(atm::Temperature(temperatureValue, atm::Temperature::UnitKelvin));
  }

  if (isValueUnset(userPressureValue_)) {
    Double pressureValue = atmPressureData_[timeIndex];
    p->setBasicAtmosphericParameters(atm::Pressure(pressureValue, atm::Pressure::UnitMilliBar));
  }

  if (isValueUnset(userRelHumidityValue_)) {
    Double relHumidityValue = atmRelHumidityData_[timeIndex];
    p->setBasicAtmosphericParameters(atm::Humidity(relHumidityValue, atm::Humidity::UnitPercent));
  }

  // TODO: configure appropriate index from time stamp
  timeIndex = 0;
  if (isValueUnset(userPwvValue_)) {
    Double pwvValue = pwvData_[timeIndex];
    p->setUserWH2O(atm::Length(pwvValue, atm::Length::UnitMilliMeter));
  }
}

void SDAtmosphereCorrectionTVI::updateCorrectionFactor(atm::SkyStatus *p) {
  // discard correction factor if nullptr is given
  if (!p) {
    correctionFactor_.resize();
    return;
  }

  // PRECONDITION:
  // all rows have the same time stamp (should be true for no RowBlocking mode)

  // find OFF_SOURCE time stamps that brackets current time stamp
  Vector<Double> currentTimeList;
  time(currentTimeList);
  Double const currentTime = currentTimeList[0];
  std::pair<Double, Double> nearest = findNearest(offSourceTime_, currentTime);
  Double const offSourceTimePrev = nearest.first;
  Double const offSourceTimeNext = nearest.second;

  // elevation
  Double const elevationOffPrev = interpolateDataLinear(elevationTime_, elevationData_, offSourceTimePrev);
  Double const elevationOn = interpolateDataLinear(elevationTime_, elevationData_, currentTime);
  Double const elevationOffNext = interpolateDataLinear(elevationTime_, elevationData_, offSourceTimeNext);
  Double const elevationOff =
    ((offSourceTimeNext - currentTime) * elevationOffPrev +
     (currentTime - offSourceTimePrev) * elevationOffNext) /
    (offSourceTimeNext - offSourceTimePrev);

  // opacity
  unsigned int const nchan = atmSpectralGridPtr_->getNumChan();
  Double const airMassOn = 1.0 / cos(C::pi_2 - elevationOn);
  atmSkyStatusPtr_->setAirMass(airMassOn);
  Vector<Double> trjSkySpecOn = getTrjSkySpec(atmSkyStatusPtr_, nchan);
  Double const airMassOff = 1.0 / cos(C::pi_2 - elevationOff);
  atmSkyStatusPtr_->setAirMass(airMassOff);
  Vector<Double> trjSkySpecOff = getTrjSkySpec(atmSkyStatusPtr_, nchan);
  Vector<Double> tauSpecOff = getTauSpec(atmSkyStatusPtr_, nchan) * airMassOff;
  Vector<Double> correctionFactor = getCorrectionFactor(trjSkySpecOn, trjSkySpecOff, tauSpecOff);

  // current spw
  SpwId const currentSpw = dataDescriptionSubtablecols().spectralWindowId().get(dataDescriptionId());
  if (isTdmSpw_[currentSpw]) {
    // TODO: interpolation to native freq
  } else if (doSmooth_[currentSpw]) {
    // TODO: smoothing
  } else {
    correctionFactor_.reference(correctionFactor);
  }
}

void SDAtmosphereCorrectionTVI::updateAtmosphereModel() {
  if (atmSkyStatusPtr_) {
    updateSkyStatus(atmSkyStatusPtr_);
    updateCorrectionFactor(atmSkyStatusPtr_);
  }
}

void SDAtmosphereCorrectionTVI::warnIfNoTransform() {
  if (!doTransform_[dataDescriptionId()]) {
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
  std::set<Double> timeOffSource = meta.getTimesForIntent("OBSERVE_TARGET#OFF_SOURCE");
  // midpoint of OFF_SOURCE scan
  size_t numBoundaries = 0;
  Vector<Double> boundaries(timeOffSource.size());
  auto iter = timeOffSource.begin();
  boundaries[numBoundaries++] = *iter;
  Double cache = *iter;
  for (auto i = iter; i != timeOffSource.end(); ++i) {
    Double diff = *i - cache;
    // threshold (1s) is taken from original script by T. Sawada
    if (diff > 1.0) {
      boundaries[numBoundaries++] = *i;
    }
    cache = *i;
  }
  boundaries[numBoundaries++] = cache;

  offSourceTime_.resize(numBoundaries - 1);
  for (size_t i = 0; i < numBoundaries - 1; ++i) {
    offSourceTime_[i] = (boundaries[i] + boundaries[i + 1]) / 2.0;
  }
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
}

void SDAtmosphereCorrectionTVI::readPointing(String const &msName, Int const referenceAntenna) {
  MeasurementSet const ms(msName);
  uInt numAntennas = ms.antenna().nrow();
  if (referenceAntenna < 0 || numAntennas <= static_cast<uInt>(referenceAntenna)) {
    LogIO os(LogOrigin("SDAtmosphereCorrectionTVI", __func__, WHERE));
    os << "ERROR: reference antenna " << referenceAntenna << " is invalid." << LogIO::EXCEPTION;
  }

  Table const pointingTable = ms.pointing();
  Table selected = pointingTable(pointingTable.col("ANTENNA_ID") == referenceAntenna);
  elevationTime_ = ScalarColumn<Double>(selected, "TIME").getColumn();
  Array<Double> elevationArray = ArrayColumn<Double>(selected, "DIRECTION")
    .getColumn()(IPosition(3, 1, 0, 0), IPosition(3, 1, 0, elevationTime_.nelements()));
  elevationArray.removeDegenerate();
  elevationData_.reference(elevationArray);
  // elevationInterpolator_ = Interpolate1D<Double, Double>(elevationTime_, elevationData_, True, True);
  // elevationInterpolator_.setMethod(Interpolate1D<Double, Double>::linear);
}

void SDAtmosphereCorrectionTVI::readAsdmAsIsTables(String const &msName) {
  Vector<Double> pwvTimeLocal;
  Vector<Double> atmTimeLocal;
  Vector<Double> pwvDataLocal;
  Vector<Double> atmTemperatureDataLocal;
  Vector<Double> atmPressureDataLocal;
  Vector<Double> atmRelHumidityDataLocal;

  File const calWvrTableName(msName + "/ASDM_CALWVR");

  if (calWvrTableName.exists()) {
    Table const calWvrTable(calWvrTableName.path().absoluteName());
    pwvTimeLocal.reference(ScalarColumn<Double>(calWvrTable, "startValidTime").getColumn());
    pwvDataLocal.reference(ScalarColumn<Double>(calWvrTable, "water").getColumn());
  }

  File const calAtmosphereTableName(msName + "/ASDM_CALATMOSPHERE");
  if (calAtmosphereTableName.exists()) {
    Table const calAtmosphereTable(calAtmosphereTableName.path().absoluteName());
    atmTimeLocal.reference(ScalarColumn<Double>(calAtmosphereTable, "startValidTime").getColumn());
    atmTemperatureDataLocal.reference(ScalarColumn<Double>(calAtmosphereTable, "groundTemperature").getColumn());
    atmPressureDataLocal.reference(ScalarColumn<Double>(calAtmosphereTable, "groundPressure").getColumn());
    atmRelHumidityDataLocal.reference(ScalarColumn<Double>(calAtmosphereTable, "groundRelHumidity").getColumn());
  }

  if (pwvTimeLocal.nelements() > 0) {
    Vector<uInt> uniqueVector;
    Vector<uInt> indexVector;
    uInt n = sortTime(pwvTimeLocal, uniqueVector, indexVector);
    pwvTime_.resize(n);
    for (uInt i = 0; i < n; ++i) {
      pwvTime_[i] = pwvTimeLocal[indexVector[uniqueVector[i]]];
    }
    // m -> mm
    pwvData_ = getMedianDataPerTime(n, uniqueVector, indexVector, pwvDataLocal) * 1000.0;
  }

  if (atmTimeLocal.nelements() > 0) {
    Vector<uInt> uniqueVector;
    Vector<uInt> indexVector;
    uInt n = sortTime(atmTimeLocal, uniqueVector, indexVector);
    atmTime_.resize(n);
    for (uInt i = 0; i < n; ++i) {
      atmTime_[i] = atmTimeLocal[indexVector[uniqueVector[i]]];
    }
    atmTemperatureData_ = getMedianDataPerTime(n, uniqueVector, indexVector, atmTemperatureDataLocal);
    // Pa -> mbar (=hPa)
    atmPressureData_ = getMedianDataPerTime(n, uniqueVector, indexVector, atmPressureDataLocal) / 100.0;
    atmRelHumidityData_ = getMedianDataPerTime(n, uniqueVector, indexVector, atmRelHumidityDataLocal);
  }
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
