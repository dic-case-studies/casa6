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

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <limits>
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
#include <casacore/ms/MeasurementSets/MSSpWindowColumns.h>
#include <casacore/ms/MeasurementSets/MSPolColumns.h>
#include <casacore/ms/MSOper/MSMetaData.h>
#include <casacore/tables/TaQL/ExprNode.h>

#include <msvis/MSVis/VisBufferComponents2.h>
#include <msvis/MSVis/VisibilityIteratorImpl2.h>

#include <mstransform/TVI/UtilsTVI.h>

using namespace casacore;

namespace {

constexpr float kNoCache = -1.0f;
constexpr double kValueUnset = std::numeric_limits<double>::quiet_NaN();
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
  Sort sort(p, timeData.nelements());
  sort.sort(indexVector, timeData.nelements());
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
  uniqueIndex.resize(n, True);
  return n;
}

inline Vector<Double> getMedianDataPerTime(uInt const n, Vector<uInt> const &uniqueVector,
  Vector<uInt> const &indexVector, Vector<Double> const &data) {
  Vector<Double> out(n);
  assert(n > 0);
  for (uInt i = 0; i < n - 1; ++i) {
    uInt const iStart = uniqueVector[i];
    uInt const iEnd = uniqueVector[i + 1];
    uInt const nElem = iEnd - iStart;
    Vector<Double> tmp(nElem);
    for (uInt j = 0; j < nElem; ++j) {
      tmp[j] = data[indexVector[j + iStart]];
    }
    out[i] = median(tmp, False, True, True);
  }
  uInt const iStart = uniqueVector[n - 1];
  uInt const iEnd = indexVector.nelements();
  uInt const nElem = iEnd - iStart;
  Vector<Double> tmp(nElem);
  for (uInt j = 0; j < nElem; ++j) {
    tmp[j] = data[indexVector[j + iStart]];
  }
  out[n - 1] = median(tmp, False, True, True);
  return out;
}

inline std::pair<Int, Int> findNearestIndex(Vector<Double> const &data, Double const target) {
  Bool found = false;
  // data is sorted in ascending order
  Vector<Double> diff(data.nelements());
  for (uInt i = 0; i < diff.nelements(); ++i) {
    diff[i] = data[i] - data[0];
  }
  Int index = binarySearch(found, data, target, data.nelements());
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

// inline Double interpolateDataLinear(Vector<Double> const &xin, Vector<Double> const &yin, Double const xout) {
//   Bool found = false;
//   // data is sorted in ascending order
//   Int index = binarySearch(found, xin, xout, xin.nelements());
//   Double yout = kValueUnset;
//   if (index == 0 || found) {
//     yout = yin[index];
//   } else {
//     yout = ((xin[index] - xout) * yin[index - 1] * (xout - xin[index - 1]) * yin[index]) / (xin[index] - xin[index - 1]);
//   }
//   return yout;
// }

// implementation of np.convolve(mode='same')
inline Vector<Double> convolve1DTriangle(Vector<Double> const &in) {
  // constexpr unsigned int kNumKernel = 3u;
  constexpr Double kKernelTriangle[] = {0.25, 0.5, 0.25};
  unsigned int const n = in.nelements();
  assert(n >= kNumKernel);
  Vector<Double> out(n, 0.0);
  // symmetric kernel
  out[0] = kKernelTriangle[0] * in[1] + kKernelTriangle[1] + in[0];
  for (unsigned int i = 1; i < n - 1; ++i) {
    out[i] = kKernelTriangle[0] * (in[i - 1] + in[i + 1]) + kKernelTriangle[1] + in[i];
  }
  out[n - 1] = kKernelTriangle[0] * in[n - 2] + kKernelTriangle[1] * in[n - 1];
  return out;
}

inline Vector<Double> convolve1DHanning(Vector<Double> const &in) {
  // normalized spectral response for Hanning window, FWHM=10
  constexpr unsigned int kNumKernel = 29u;
  constexpr unsigned int kIndexKernelCenter = kNumKernel / 2u;
  constexpr Double kKernelHanning[] = {
    -0.00098041, -0.00202866, -0.00265951, -0.00222265,
    0.00000000, 0.00465696, 0.01217214, 0.02260546, 0.03556241,
    0.05017949, 0.06519771, 0.07911911, 0.09042184, 0.09779662,
    0.10035898, 0.09779662, 0.09042184, 0.07911911, 0.06519771,
    0.05017949, 0.03556241, 0.02260546, 0.01217214, 0.00465696,
    0.00000000, -0.00222265, -0.00265951, -0.00202866, -0.00098041
  };
  unsigned int const n = in.nelements();
  assert(n >= kNumKernel);
  Vector<Double> out(n, 0.0);
  // symmetric kernel
  for (unsigned int i = 0u; i < kIndexKernelCenter; ++i) {
    out[i] = 0.0;
    for (unsigned int j = 0u; j < kIndexKernelCenter + i + 1; ++j) {
      unsigned int k = kIndexKernelCenter - i + j;
      out[i] += in[j] * kKernelHanning[k];
    }
  }
  for (unsigned int i = kIndexKernelCenter; i < n - kIndexKernelCenter; ++i) {
    out[i] = 0.0;
    for (unsigned int k = 0u; k < kNumKernel; ++k) {
      unsigned int j = i - kIndexKernelCenter + k;
      out[i] += in[j] * kKernelHanning[k];
    }
  }
  for (unsigned int i = n - kIndexKernelCenter; i < n; ++i) {
    out[i] = 0.0;
    for (unsigned int j = i - kIndexKernelCenter; j < n; ++j) {
      unsigned int k = j - (i - kIndexKernelCenter);
      out[i] += in[j] * kKernelHanning[k];
    }
  }
  return out;
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
  // cout << "gainfactor (" << gainFactor.shape() << ") = " << gainFactor[0] << endl;
  // cout << "input array shape = " << in.shape() << endl;
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
    configuration_(configuration),
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
    isTdmSpw_(),
    doSmooth_(),
    nchanBBPerSpw_(),
    channelFreqsPerSpw_(),
    channelWidthsPerSpw_(),
    currentTime_(kValueUnset),
    currentAtmTimeIndex_(-1),
    currentSpwId_(-1),
    atmType_(2),
    atmSkyStatusPerSpw_(),
    atmSkyStatusPtr_(nullptr),
    correctionFactor_() {

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

  updateCache();

  // update SkyStatus with the closest weather measurement
  // re-calculate correction factor (dTa)
  configureAtmosphereCorrection();
  updateAtmosphereModel();
}

void SDAtmosphereCorrectionTVI::next() {
  TransformingVi2::next();

  // Synchronize own VisBuffer
  configureNewSubchunk();

  updateCache();

  // update SkyStatus with the closest weather measurement
  // re-calculate correction factor (dTa)
  configureAtmosphereCorrection();
  updateAtmosphereModel();
}

void SDAtmosphereCorrectionTVI::originChunks(Bool forceRewind) {
  TransformingVi2::originChunks(forceRewind);

  // initialization
  initializeAtmosphereCorrection(configuration_);

  // setup some cache values
  updateCache();
  currentSpwId_ = dataDescriptionSubtablecols().spectralWindowId().get(dataDescriptionId());

  // initialize Atmosphere Transmission Model
  initializeAtmosphereModel(configuration_);


  // configure atmospheric correction
  // point appropriate SkyStatus object
  configureAtmosphereCorrection();
  updateAtmosphereModel();
}

void SDAtmosphereCorrectionTVI::nextChunk() {
  TransformingVi2::nextChunk();

  // setup some cache values
  updateCache();
  currentSpwId_ = dataDescriptionSubtablecols().spectralWindowId().get(dataDescriptionId());

  // configure atmospheric correction if necessary
  // point appropriate SkyStatus object
  configureAtmosphereCorrection();
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

  MSMetaData msmd(&ms(), kNoCache);
  std::set<uInt> allSpwIds = msmd.getSpwIDs();
  std::vector<uInt> nonProcessingSpws;
  nonProcessingSpws.reserve(allSpwIds.size());
  for (auto i = allSpwIds.begin(); i != allSpwIds.end(); ++i) {
    if (allNE(processSpwList_, static_cast<Int>(*i))) {
      nonProcessingSpws.push_back(*i);
    }
  }
  if (nonProcessingSpws.size() > 0) {
    os << LogIO::WARN << "SPWs ";
    for (auto i = nonProcessingSpws.begin(); i != nonProcessingSpws.end(); ++i) {
      os << *i << " ";
    }
    os << "are output but not corrected" << LogIO::POST;
  }

  // read MAIN table (OFF_SOURCE time)
  readMain(rawMs);

  // read SPECTRAL_WINDOW table
  readSpectralWindow(rawMs);

  // read POINTING table
  readPointing(rawMs, referenceAntenna);

  // read ASDM_CALWVR and ASDM_CALATMOSPHERE tables
  readAsdmAsIsTables(rawMs);

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
      defaultPressureValue = atmPressureData_[currentAtmTimeIndex_];
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
      defaultTemperatureValue = atmTemperatureData_[currentAtmTimeIndex_];
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
      defaultRelHumidityValue = atmRelHumidityData_[currentAtmTimeIndex_];
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
      defaultPwvValue = pwvData_[currentAtmTimeIndex_];
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
    // size_t const numLayer = std::min(layerBoundariesData.nelements(),
    //                                  layerTemperaturesData.nelements());
    if (layerBoundariesData.size() != layerTemperaturesData.size()) {
      os << "ERROR: list length of layerboundaries and layertemperature should be the same." << LogIO::EXCEPTION;
    }
    size_t const numLayer = layerBoundariesData.size();
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
  os << "user-defined layer:" << endl;
  if (layerBoundaries.size() > 0) {
    for (size_t i = 0; i < layerBoundaries.size(); ++i) {
      os << "  Height " << layerBoundaries[i].get(atm::Length::UnitKiloMeter)
         << " Temperature " << layerTemperatures[i].get(atm::Temperature::UnitKelvin)
         << endl;
    }
    os << LogIO::POST;
  } else {
    os << "NONE" << LogIO::POST;
  }

  atm::AtmProfile atmProfile(
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
  );

  // disable threshold so that any subtle changes in
  // temperature/pressure/humidity/pwv cause to udpate
  // the model
  atmProfile.setBasicAtmosphericParameterThresholds(
    atm::Length(0.0, atm::Length::UnitMeter),
    atm::Pressure(0.0, atm::Pressure::UnitMilliBar),
    atm::Temperature(0.0, atm::Temperature::UnitKelvin),
    0.0,
    atm::Humidity(0.0, atm::Humidity::UnitPercent),
    atm::Length(0.0, atm::Length::UnitMeter)
  );

  for (auto i = processSpwList_.begin(); i != processSpwList_.end(); ++i) {
    SpwId const spw = *i;
    Vector<Double> cf = channelFreqsPerSpw_[spw];
    Vector<Double> cw = channelWidthsPerSpw_[spw];
    unsigned int nchan = cf.nelements();

    // smoothing control
    Int polId = dataDescriptionSubtablecols().polarizationId().get(dataDescriptionId());
    Int numPol = polarizationSubtablecols().numCorr().get(polId);
    uInt numChanPerBB = nchanBBPerSpw_[spw];
    Int numChanPol = numChanPerBB * numPol;
    doSmooth_[spw] = (numChanPol == 256 || numChanPol == 8192);

    unsigned int refChan = (nchan - 1) / 2;
    double centerFreq = cf[refChan];
    Double chanSep = (nchan == 1u) ? cw[0] : (cf[nchan - 1] - cf[0]) / static_cast<Double>(nchan - 1);
    if (isTdmSpw_[spw]) {
      // configure 5x finer spectral grid than native one
      chanSep /= 5.0;
      nchan *= 5u;
      refChan = refChan * 5 + 2u;
    }

    // SpectralGrid
    os.output() << std::setprecision(16);
    os << "SpectralGrid for spw " << spw << ": nchan " << nchan
       << " refchan " << refChan << " center freq "
       << centerFreq << " chan sep " << chanSep << LogIO::POST;
    atm::SpectralGrid spectralGrid(nchan, refChan, centerFreq, chanSep);

    // SkyStatus
    os << "initializing SkyStatus for " << spw << LogIO::POST;
    atm::RefractiveIndexProfile profile(spectralGrid, atmProfile);
    atmSkyStatusPerSpw_[spw].reset(new atm::SkyStatus(profile));
    atmSkyStatusPerSpw_[spw]->setUserWH2O(atm::Length(defaultPwvValue, atm::Length::UnitMilliMeter));
  }
  os << "DONE Initializing SkyStatus" << LogIO::POST;

}

void SDAtmosphereCorrectionTVI::configureAtmosphereCorrection() {
  LogIO os(LogOrigin("SDAtmosphereCorrectionTVI", __func__, WHERE));
  bool isProcessingSpw = anyEQ(processSpwList_, currentSpwId_);
  bool isOnSource = isOnSourceChunk();
  bool isPrecedingOffSourceScanExist = min(offSourceTime_) <= currentTime_;
  bool isSubsequentOffSourceScanExist = currentTime_ <= max(offSourceTime_);
  bool isPrecedingAtmScanExist = min(atmTime_) <= currentTime_;
  // cout << std::setprecision(16) << offSourceTime_ << endl;
  // cout << std::setprecision(16) << currentTime_ << " state ID" << currentStateId[0] << endl;
  // os << "SPW " << currentSpwId_ << ": processingSpw " << isProcessingSpw
  //    << " OFF_SOURCE availability before " << isPrecedingOffSourceScanExist
  //    << " after " << isSubsequentOffSourceScanExist
  //    << " ON_SOURCE? " << isOnSource << LogIO::POST;
  atmSkyStatusPtr_ = nullptr;
  if (isProcessingSpw && isOnSource &&
      isPrecedingOffSourceScanExist && isSubsequentOffSourceScanExist &&
      isPrecedingAtmScanExist) {
    // gain factor for current spw
    if (currentSpwId_ < 0 || static_cast<SpwId>(gainFactorList_.nelements()) <= currentSpwId_) {
      gainFactor_ = 1.0;
    } else {
      gainFactor_ = gainFactorList_[currentSpwId_];
    }
    // os << "gainfactor for SPW " << currentSpwId_<< " = " << gainFactor_ << LogIO::POST;

    // SkyStatus for current spw
    auto finder = atmSkyStatusPerSpw_.find(currentSpwId_);
    if (finder != atmSkyStatusPerSpw_.end()) {
      atmSkyStatusPtr_ = finder->second.get();
    }
  }
}

void SDAtmosphereCorrectionTVI::updateSkyStatus() {
  LogIO os(LogOrigin("SDAtmosphereCorrectionTVI", __func__, WHERE));
  // os << "updateSkyStatus for SPW " << currentSpwId_ << LogIO::POST;

  // do nothing if nullptr is given
  if (!atmSkyStatusPtr_) {
    return;
  }

  Int timeIndex = currentAtmTimeIndex_;

  bool isSkyStatusOutdated = false;
  Double currentTemperatureValue = userTemperatureValue_;
  if (isValueUnset(currentTemperatureValue)) {
    Double val = atmTemperatureData_[timeIndex];
    Double val2 = atmSkyStatusPtr_->getGroundTemperature().get(atm::Temperature::UnitKelvin);
    if (val != val2) {
      // cout << "Temperature is different " << val << " vs " << val2 << endl;
      isSkyStatusOutdated = true;
    }
    currentTemperatureValue = val;
  }

  Double currentPressureValue = userPressureValue_;
  if (isValueUnset(currentPressureValue)) {
    Double val = atmPressureData_[timeIndex];
    Double val2 = atmSkyStatusPtr_->getGroundPressure().get(atm::Pressure::UnitMilliBar);
    if (val != val2) {
      // cout << "Pressure is different " << val << " vs " << val2 << endl;
      isSkyStatusOutdated = true;
    }
    currentPressureValue = val;
  }

  Double currentRelHumidityValue = userRelHumidityValue_;
  if (isValueUnset(currentRelHumidityValue)) {
    Double val = atmRelHumidityData_[timeIndex];
    Double val2 = atmSkyStatusPtr_->getRelativeHumidity().get(atm::Humidity::UnitPercent);
    if (val != val2) {
      // cout << "Humidity is different " << val << " vs " << val2 << endl;
      isSkyStatusOutdated = true;
    }
    currentRelHumidityValue = val;
  }

  if (isSkyStatusOutdated) {
    os << "updating SkyStatus for SPW " << currentSpwId_ << LogIO::POST;
    atmSkyStatusPtr_->setBasicAtmosphericParameters(
      atmSkyStatusPtr_->getAltitude(),
      atm::Pressure(currentPressureValue, atm::Pressure::UnitMilliBar),
      atm::Temperature(currentTemperatureValue, atm::Temperature::UnitKelvin),
      atmSkyStatusPtr_->getTropoLapseRate(),
      atm::Humidity(currentRelHumidityValue, atm::Humidity::UnitPercent),
      atmSkyStatusPtr_->getWvScaleHeight()
    );
  }

  // PYTHON IMPLEMENTATION
  // share time information with CalAtmosphere table
  if (isValueUnset(userPwvValue_)) {
    Double val = pwvData_[timeIndex];
    Double val2 = atmSkyStatusPtr_->getUserWH2O().get(atm::Length::UnitMilliMeter);
    if (val != val2) {
      // cout << "PWV is different " << val << " vs " << val2 << endl;
      os << "updaging PWV value for SPW " << currentSpwId_ << LogIO::POST;
      atmSkyStatusPtr_->setUserWH2O(atm::Length(val, atm::Length::UnitMilliMeter));
    }
  }

  // cout << "time " << currentTime_
  //      << " Temperature " << atmTemperatureData_[timeIndex] << "K "
  //      << " Pressure " << atmPressureData_[timeIndex] << "mbar "
  //      << " Humidity " << atmRelHumidityData_[timeIndex] << "%" << endl;

  // os << "DONE updateSkyStatus for SPW " << currentSpwId_ << LogIO::POST;
}

void SDAtmosphereCorrectionTVI::updateCorrectionFactor() {
  LogIO os(LogOrigin("SDAtmosphereCorrectionTVI", __func__, WHERE));
  // os << "updateCorrectionFactor for SPW " << currentSpwId_ << LogIO::POST;

  // discard correction factor if nullptr is given
  if (!atmSkyStatusPtr_) {
    correctionFactor_.resize();
    return;
  }

  // PRECONDITION:
  // all rows have the same time stamp (should be true for no RowBlocking mode)

  // find OFF_SOURCE time stamps that brackets current time stamp
  std::pair<Double, Double> nearest = findNearest(offSourceTime_, currentTime_);
  Double const offSourceTimePrev = nearest.first;
  Double const offSourceTimeNext = nearest.second;
  // cout << "nearest: " << std::setprecision(16) << offSourceTimePrev << "~" << offSourceTimeNext << endl;

  // elevation
  // cout << "interpolated elevation" << endl;
  Double const elevationOffPrev = elevationInterpolator_(offSourceTimePrev);
  // cout << "prev: " << elevationOffPrev << endl;
  Double const elevationOn = elevationInterpolator_(currentTime_);
  // cout << "ON: " << elevationOn << endl;
  Double const elevationOffNext = elevationInterpolator_(offSourceTimeNext);
  // cout << "next: " << elevationOffNext << endl;
  Double const elevationOff =
    ((offSourceTimeNext - currentTime_) * elevationOffPrev +
     (currentTime_ - offSourceTimePrev) * elevationOffNext) /
    (offSourceTimeNext - offSourceTimePrev);
  // cout << "OFF: " << elevationOffPrev << endl;
  // cout << "time " << std::setprecision(16) << currentTime_ << " elON " << elevationOn << " elOFF " << elevationOff << endl;

  // opacity
  unsigned int const nchan = atmSkyStatusPtr_->getNumChan();
  Double const airMassOn = 1.0 / cos(C::pi_2 - elevationOn);
  atmSkyStatusPtr_->setAirMass(airMassOn);
  Vector<Double> trjSkySpecOn = getTrjSkySpec(atmSkyStatusPtr_, nchan);
  Double const airMassOff = 1.0 / cos(C::pi_2 - elevationOff);
  atmSkyStatusPtr_->setAirMass(airMassOff);
  Vector<Double> trjSkySpecOff = getTrjSkySpec(atmSkyStatusPtr_, nchan);
  Vector<Double> tauSpecOff = getTauSpec(atmSkyStatusPtr_, nchan, airMassOff);
  Vector<Double> correctionFactor = getCorrectionFactor(trjSkySpecOn, trjSkySpecOff, tauSpecOff);

  // constexpr size_t ic = 0;
  // cout << std::setprecision(16)
  //      << "time = " << currentTime_
  //      << " TrjON = " << trjSkySpecOn[ic]
  //      << " TrjOFF = " << trjSkySpecOff[ic]
  //      << " tauOFF = " << tauSpecOff[ic]
  //      << " factor = " << correctionFactor[ic] << " (without gain factor)"
  //      << endl;

  // current spw
  // SpwId const currentSpw = dataDescriptionSubtablecols().spectralWindowId().get(dataDescriptionId());
  if (isTdmSpw_[currentSpwId_]) {
    // cout << "SPW " << currentSpwId_ << " is TDM " << endl;
    unsigned int numAtmChan = atmSkyStatusPtr_->getNumChan();
    unsigned int refAtmChan = atmSkyStatusPtr_->getRefChan();
    double refAtmFreq = atmSkyStatusPtr_->getRefFreq().get();
    double sepAtmFreq = atmSkyStatusPtr_->getChanSep().get();
    Vector<Double> atmFreq(numAtmChan);
    indgen(atmFreq, refAtmFreq - sepAtmFreq * refAtmChan, sepAtmFreq);
    Vector<Double> const &nativeFreq = channelFreqsPerSpw_[currentSpwId_];
    Vector<Double> smoothedCorrectionFactor = convolve1DHanning(correctionFactor);
    // atm frequency grid is 5x finer than native frequency grid
    correctionFactor_.resize(nativeFreq.nelements());
    for (uInt i = 0; i < nativeFreq.nelements(); ++i) {
      correctionFactor_[i] = smoothedCorrectionFactor[2 + i * 5];
    }
  } else if (doSmooth_[currentSpwId_]) {
    // cout << "SPW " << currentSpwId_ << " requires smoothing " << endl;
    correctionFactor_ = convolve1DTriangle(correctionFactor);
  } else {
    correctionFactor_.reference(correctionFactor);
  }

  // apply gain factor
  // cout << "Applying gain factor " << gainFactorList_[currentSpwId_] << " to correction term" << endl;
  correctionFactor_ *= gainFactorList_[currentSpwId_];

  // os << "DONE updateCorrectionFactor for SPW " << currentSpwId_ << LogIO::POST;
}

void SDAtmosphereCorrectionTVI::updateAtmosphereModel() {
  if (atmSkyStatusPtr_) {
    updateSkyStatus();
    updateCorrectionFactor();
  }
}

void SDAtmosphereCorrectionTVI::updateCache() {
  // non rowBlocking mode - time should be single value
  Vector<Double> timeData;
  time(timeData);
  currentTime_ = timeData[0];
  if (atmTime_.nelements() > 0) {
    std::pair<Int, Int> pair = findNearestIndex(atmTime_, currentTime_);
    currentAtmTimeIndex_ = pair.first;
    // cout << "updateCache: time " << std::setprecision(16) << currentTime_
    //      << " index " << currentAtmTimeIndex_ << endl;
  }
}

void SDAtmosphereCorrectionTVI::readMain(String const &msName) {
  MeasurementSet msObj(msName);
  MSMetaData meta(&msObj, kNoCache);
  std::set<Double> timeOffSourceSet = meta.getTimesForIntent("OBSERVE_TARGET#OFF_SOURCE");
  std::vector<Double> timeOffSource(timeOffSourceSet.begin(), timeOffSourceSet.end());
  size_t numBoundaries = 0;
  Vector<uInt> idx1(timeOffSource.size());
  Double cache = timeOffSource[0];
  idx1[numBoundaries++] = 0;
  for (uInt i = 0; i < timeOffSource.size(); ++i) {
    Double diff = timeOffSource[i] - cache;
    // PYTHON IMPLEMENTATION
    // threshold (1s) is taken from original script by T. Sawada
    if (diff > 1.0) {
      idx1[numBoundaries++] = i;
    }
    cache = timeOffSource[i];
  }
  Vector<uInt> idx2(numBoundaries);
  for (uInt i = 0; i < numBoundaries - 1; ++i) {
    idx2[i] = idx1[i + 1] - 1;
  }
  idx2[numBoundaries - 1] = timeOffSource.size() - 1;
  offSourceTime_.resize(numBoundaries);
  for (size_t i = 0; i < numBoundaries; ++i) {
    offSourceTime_[i] = (timeOffSource[idx1[i]] + timeOffSource[idx2[i]]) / 2.0;
  }
  Vector<Double> tmp(numBoundaries);
  for (auto i = 0u; i < numBoundaries; ++i) {
    tmp[i] = offSourceTime_[i] - offSourceTime_[0];
  }
}

void SDAtmosphereCorrectionTVI::readSpectralWindow(String const &msName) {
  MeasurementSet const ms(msName);

  Table const spwTable = ms.spectralWindow();
  auto const chanFreqColumn = spectralWindowSubtablecols().chanFreq();
  auto const chanWidthColumn = spectralWindowSubtablecols().chanWidth();
  for (auto i = processSpwList_.begin(); i != processSpwList_.end(); ++i) {
    SpwId const spw = *i;
    channelFreqsPerSpw_[spw] = chanFreqColumn(spw);
    channelWidthsPerSpw_[spw] = chanWidthColumn(spw);
  }
  LogIO os(LogOrigin("SDAtmosphereCorrectionTVI", __func__, WHERE));

  // number of channels per baseband
  MSMetaData meta(&ms, kNoCache);
  std::set<uInt> fdmSpws = meta.getFDMSpw();
  std::set<uInt> tdmSpws = meta.getTDMSpw();
  std::set<uInt> onSourceSpws = meta.getSpwsForIntent("OBSERVE_TARGET#ON_SOURCE");
  std::set<uInt> fdmTdmSpws;
  std::set_union(
    fdmSpws.begin(), fdmSpws.end(),
    tdmSpws.begin(), tdmSpws.end(),
    std::inserter(fdmTdmSpws, fdmTdmSpws.end())
  );
  std::set<uInt> scienceSpws;
  std::set_intersection(
    onSourceSpws.begin(), onSourceSpws.end(),
    fdmTdmSpws.begin(), fdmTdmSpws.end(),
    std::inserter(scienceSpws, scienceSpws.end())
  );
  std::map<Int, uInt> nchanPerBB;
  std::map<SpwId, Int> spwBBMap;
  auto const nchanColumn = spectralWindowSubtablecols().numChan();
  auto const nameColumn = spectralWindowSubtablecols().name();
  for (auto i = scienceSpws.begin(); i != scienceSpws.end(); ++i) {
    SpwId const spw = *i;
    String const name = nameColumn(spw);
    Int nchan = nchanColumn(spw);
    auto const pos = name.find("#BB_");
    auto const pos2 = name.find("#", pos + 1);
    if (pos != String::npos && pos + 4 < name.size()) {
      Int bb = String::toInt(name.substr(pos + 4, pos2 - pos - 4));
      // cout << "SPW " << spw << ": name \"" << name
      //      << "\" substr " << name.substr(pos + 4, pos2 - pos - 4)
      //      << " pos " << pos + 4 << " pos2 " << pos2
      //      << " bb " << bb << endl;
      spwBBMap[spw] = bb;
      auto j = nchanPerBB.find(bb);
      if (j == nchanPerBB.end()) {
        nchanPerBB[bb] = nchan;
      } else {
        nchanPerBB[bb] += nchan;
      }
    }
  }
  for (auto i = processSpwList_.begin(); i != processSpwList_.end(); ++i) {
    SpwId const spw = *i;
    auto j = spwBBMap.find(spw);
    if (j == spwBBMap.end()) {
      nchanBBPerSpw_[spw] = 0u;
    } else {
      nchanBBPerSpw_[spw] = nchanPerBB[j->second];
    }
    auto const nchan = nchanBBPerSpw_[spw];
    os << "nchan per BB for SPW " << spw << " is " << nchan << LogIO::POST;
    isTdmSpw_[spw] = (nchan == 128u || nchan == 256u);
    os << "SPW " << spw << " is " << ((isTdmSpw_[spw]) ? "TDM" : "FDM") << " like" << LogIO::POST;
  }

  os << "TDM flag:" << endl;
  os.output() << std::boolalpha;
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
  Table const pointingTable = ms.pointing();
  // TODO: abort if direction reference is not AZEL
  Table selected = pointingTable(pointingTable.col("ANTENNA_ID") == referenceAntenna);
  elevationTime_ = ScalarColumn<Double>(selected, "TIME").getColumn();
  Array<Double> directionArray = ArrayColumn<Double>(selected, "DIRECTION").getColumn();
  Array<Double> elevationArray = Cube<Double>(directionArray).yzPlane(1);
  elevationArray.removeDegenerate();
  elevationData_.reference(elevationArray);
  elevationInterpolator_ = Interpolate1D<Double, Double>(
    ScalarSampledFunctional<Double>(elevationTime_),
    ScalarSampledFunctional<Double>(elevationData_), True, True);
  elevationInterpolator_.setMethod(Interpolate1D<Double, Double>::linear);
  // os.output() << std::setprecision(16);
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
  }

  File const calAtmosphereTableName(msName + "/ASDM_CALATMOSPHERE");
  os << "CALATMOSPHERE table \"" << calAtmosphereTableName.path().absoluteName() << "\"" << LogIO::POST;
  if (calAtmosphereTableName.exists()) {
    Table const calAtmosphereTable(calAtmosphereTableName.path().absoluteName());
    atmTimeLocal.reference(ScalarColumn<Double>(calAtmosphereTable, "startValidTime").getColumn());
    atmTemperatureDataLocal.reference(ScalarColumn<Double>(calAtmosphereTable, "groundTemperature").getColumn());
    atmPressureDataLocal.reference(ScalarColumn<Double>(calAtmosphereTable, "groundPressure").getColumn());
    atmRelHumidityDataLocal.reference(ScalarColumn<Double>(calAtmosphereTable, "groundRelHumidity").getColumn());
  }

  Vector<uInt> uniqueVectorPwv;
  Vector<uInt> indexVectorPwv;
  if (pwvTimeLocal.nelements() > 0) {
    sortTime(pwvTimeLocal, indexVectorPwv);
    uInt n = makeUnique(pwvTimeLocal, indexVectorPwv, uniqueVectorPwv);
    pwvTime_.resize(n);
    for (uInt i = 0; i < n; ++i) {
      pwvTime_[i] = pwvTimeLocal[indexVectorPwv[uniqueVectorPwv[i]]];
    }
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
    atmTemperatureData_ = getMedianDataPerTime(n, uniqueVector, indexVector, atmTemperatureDataLocal);
    // Pa -> mbar (=hPa)
    atmPressureData_ = getMedianDataPerTime(n, uniqueVector, indexVector, atmPressureDataLocal) / 100.0;
    atmRelHumidityData_ = getMedianDataPerTime(n, uniqueVector, indexVector, atmRelHumidityDataLocal);

    // PYTHON IMPLEMENTATION
    // use time stamp in CALATMOSPHERE table
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
    // m -> mm
    pwvData_ = getMedianDataPerTime(n, uniqueVectorPwv2, indexVectorPwv, pwvDataLocal) * 1000.0;
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
