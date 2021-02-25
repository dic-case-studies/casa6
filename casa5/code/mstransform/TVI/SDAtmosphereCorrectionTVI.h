//# SDAtmosphereCorrectionTVI.h: Transforming VI for polarization averaging
//# Copyright (C) 1996,1997,1998,1999,2000,2001,2002,2003
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU Library General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the Implied warranty of MERCHANTABILITY or
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

#ifndef _SD_ATMOSPHERE_CORRECTION_TVI_H_
#define _SD_ATMOSPHERE_CORRECTION_TVI_H_

#include <casacore/casa/aips.h>

#include <atmosphere/ATM/ATMProfile.h>
#include <atmosphere/ATM/ATMSpectralGrid.h>
#include <atmosphere/ATM/ATMSkyStatus.h>

#include <msvis/MSVis/ViImplementation2.h>
#include <msvis/MSVis/TransformingVi2.h>
#include <msvis/MSVis/VisibilityIterator2.h>

#include <map>
#include <vector>

#include <casacore/measures/Measures/Stokes.h>
#include <casacore/casa/Containers/Record.h>
#include <casacore/ms/MeasurementSets/MSStateColumns.h>
#include <casacore/scimath/Functionals/Interpolate1D.h>

namespace casa { //# NAMESPACE CASA - BEGIN

namespace vi { // # NAMESPACE VI - BEGIN

//# forward decl

class VisBuffer2;

class ChannelSelector;
class ChannelSelectorCache;
typedef casacore::Vector<casacore::Vector<casacore::Slice> > ChannelSlicer;
class SpectralWindowChannelsCache;
class SpectralWindowChannels;
class SubtableColumns;

class SDAtmosphereCorrectionVi2Factory;

// <summary>
// SDAtmosphereCorrectionTVI
// </summary>

// <use visibility=export>

// <reviewed reviewer="" date="yyyy/mm/dd" tests="" demos="">
// </reviewed>

// <prerequisite>
//   <li> <linkto class="MSIter">MSIter</linkto>
//   <li> <linkto class="casacore::MeasurementSet">casacore::MeasurementSet</linkto>
//   <li> <linkto class="VisSet">VisSet</linkto>
// </prerequisite>
//
// <etymology>
// SDAtmosphereCorrectionTVI
// </etymology>
//
// <synopsis>
//
// </synopsis>
//
// <example>
// <code>
// //
// </code>
// </example>
//
// <motivation>
//
// </motivation>
//
// <thrown>
//    <li>
//    <li>
// </thrown>
//
// <todo asof="1997/05/30">
//   <li> cleanup the currently dual interface for visibilities and flags
//   <li> sort out what to do with weights when interpolating
// </todo>

class SDAtmosphereCorrectionTVI final : public TransformingVi2 {

public:
  using SpwId = casacore::Int;

  // Destructor

  virtual ~SDAtmosphereCorrectionTVI();

  // Report the the ViImplementation type
  //  (should be specialized in child classes)
  virtual casacore::String ViiType() const {
    return casacore::String("SDAtmosphereCorrection( ") + getVii()->ViiType() + " )";
  }

  // Methods to control and monitor subchunk iteration

  virtual void origin() override;
  virtual void next() override;
  virtual void originChunks(casacore::Bool forceRewind = false) override;
  virtual void nextChunk() override;

  // Return the visibilities as found in the casacore::MS, casacore::Cube (npol,nchan,nrow).

  virtual void visibilityCorrected(
  casacore::Cube<casacore::Complex> & vis) const;
  virtual void visibilityModel(casacore::Cube<casacore::Complex> & vis) const;
  virtual void visibilityObserved(
  casacore::Cube<casacore::Complex> & vis) const;

  // Return FLOAT_DATA as a casacore::Cube (npol, nchan, nrow) if found in the MS.

  virtual void floatData(casacore::Cube<casacore::Float> & fcube) const;

protected:

  // Constructor

  SDAtmosphereCorrectionTVI(ViImplementation2 * inputVi,
    casacore::Record const &configuration);

private:
  // initial configuratin of the correction
  void initializeAtmosphereCorrection(casacore::Record const &configuration);

  // initialize AtmosphereModel
  void initializeAtmosphereModel(casacore::Record const &configuration);

  // sync with the current chunk/subchunk
  void configureAtmosphereCorrection();
  void updateSkyStatus(atm::SkyStatus *p);
  void updateCorrectionFactor(atm::SkyStatus *p);
  void updateAtmosphereModel();

  // emit warning if no transform is applied to the current chunk
  void warnIfNoTransform();

  // read necessary data from MS
  void readMain(casacore::String const &msName);
  void readSpectralWindow(casacore::String const &msName);
  void readPointing(casacore::String const &msName, casacore::Int const referenceAntenna);
  void readAsdmAsIsTables(casacore::String const &msName);

  // check if transformation is necessary
  bool isOnSourceChunk() const {
    casacore::Vector<casacore::Int> stateIdList;
    stateId(stateIdList);
    static std::string const startstr("OBSERVE_TARGET#ON_SOURCE");
    casacore::String s = stateSubtablecols().obsMode().get(stateIdList[0]);
    std::cout << "state ID = " << stateIdList << " OBSMODE \"" << s << "\" match " << s.startsWith(startstr) << std::endl;
    return s.startsWith(startstr);
  }
  bool doTransform() const {return (atmSkyStatusPtr_ != nullptr);}

  // user inputs
  casacore::Vector<SpwId> processSpwList_;
  casacore::Vector<casacore::Double> gainFactorList_;
  casacore::Double gainFactor_;
  casacore::Double userPressureValue_;
  casacore::Double userTemperatureValue_;
  casacore::Double userRelHumidityValue_;
  casacore::Double userPwvValue_;

  // measurements recorded in MS
  casacore::Vector<casacore::Double> offSourceTime_;

  casacore::Vector<casacore::Double> elevationTime_;
  casacore::Vector<casacore::Double> elevationData_;
  // casacore::Interpolate1D<casacore::Double, casacore::Double> elevationInterpolator_;

  casacore::Vector<casacore::Double> pwvTime_;
  casacore::Vector<casacore::Double> pwvData_;

  casacore::Vector<casacore::Double> atmTime_;
  casacore::Vector<casacore::Double> atmTemperatureData_;
  casacore::Vector<casacore::Double> atmPressureData_;
  casacore::Vector<casacore::Double> atmRelHumidityData_;

  // Per data description flags whether transformation (i.e. atmospheric correction)
  // must be executed or not
  casacore::Vector<casacore::Bool> doTransform_;

  // spw related
  std::map<SpwId, bool> isTdmSpw_;
  std::map<SpwId, bool> doSmooth_;
  std::map<SpwId, casacore::Vector<casacore::Double> > channelFreqsPerSpw_;
  std::map<SpwId, casacore::Vector<casacore::Double> > channelWidthsPerSpw_;

  // ATM
  std::unique_ptr<atm::AtmProfile> atmProfile_;
  std::map<SpwId, std::unique_ptr<atm::SpectralGrid> > atmSpectralGridPerSpw_;
  atm::SpectralGrid *atmSpectralGridPtr_;
  std::map<SpwId, std::unique_ptr<atm::SkyStatus> > atmSkyStatusPerSpw_;
  atm::SkyStatus *atmSkyStatusPtr_;
  casacore::Vector<casacore::Double> correctionFactor_;

  friend SDAtmosphereCorrectionVi2Factory;
};


// <summary>
// A factory for generating ViImplementation2 for polarization averaging.
// </summary>
//
// <use visibility=export>
//
// <prerequisite>
//   <li> <linkto class="VisibilityIterator2:description">VisibilityIterator2</linkto>
// </prerequisite>
//
// <etymology>
// Factory for layered ViImplementation2 construction
// </etymology>
//
// <synopsis>
// PolAverageVi2Factory
// </synopsis>
//
// <motivation>
//
// </motivation>
//
// <example>
//
// </example>

class SDAtmosphereCorrectionVi2Factory final : public ViFactory {

public:
  // Constructor
  SDAtmosphereCorrectionVi2Factory(casacore::Record const &configuration,
    ViImplementation2 *inputVII);
  SDAtmosphereCorrectionVi2Factory(casacore::Record const &configuration,
    casacore::MeasurementSet const *ms, SortColumns const sortColumns,
    casacore::Double timeInterval, casacore::Bool isWritable);

  // Destructor
  ~SDAtmosphereCorrectionVi2Factory();

  ViImplementation2 * createVi() const;

private:
  ViImplementation2 *inputVII_p;
  casacore::Record configuration_p;
};

class SDAtmosphereCorrectionTVILayerFactory final : public ViiLayerFactory {

public:
  SDAtmosphereCorrectionTVILayerFactory(casacore::Record const &configuration);
  virtual ~SDAtmosphereCorrectionTVILayerFactory() {
  }
  ;

protected:

  ViImplementation2 * createInstance(ViImplementation2* vii0) const;

  casacore::Record configuration_p;

};

} //# NAMESPACE VI - END

} //# NAMESPACE CASA - END

#endif // _SD_ATMOSPHERE_CORRECTION_TVI_H_

