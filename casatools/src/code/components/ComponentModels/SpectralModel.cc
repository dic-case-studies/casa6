//# SpectralModel.cc:
//# Copyright (C) 1998,1999,2000,2002,2003
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
//# $Id: SpectralModel.cc 20652 2009-07-06 05:04:32Z Malte.Marquarding $

#include <components/ComponentModels/SpectralModel.h>
#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/Containers/RecordFieldId.h>
#include <casacore/casa/Containers/RecordInterface.h>
#include <casacore/casa/Exceptions/Error.h>
#include <casacore/casa/Arrays/IPosition.h>
#include <casacore/measures/Measures/MeasureHolder.h>
#include <casacore/measures/Measures/MFrequency.h>
#include <casacore/measures/Measures/MCFrequency.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/casa/Quanta/QuantumHolder.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/casa/Utilities/Assert.h>
#include <casacore/casa/Utilities/DataType.h>
#include <casacore/casa/BasicSL/String.h>
#include <casacore/casa/Logging/LogIO.h>
#include <casacore/casa/Logging/LogOrigin.h>
#include <iostream>

using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN

SpectralModel::SpectralModel()
  :itsRefFreq(Quantum<Double>(1, "GHz"), MFrequency::DEFAULT),
   itsFreqUnit("GHz"),
   itsFreqErr(0, itsFreqUnit)
{
  DebugAssert(SpectralModel::ok(), AipsError);
}

SpectralModel::SpectralModel(const MFrequency& refFreq, const Unit& freqUnit)
  :itsRefFreq(refFreq),
   itsFreqUnit(freqUnit),
   itsFreqErr(0, itsFreqUnit)
{
  DebugAssert(SpectralModel::ok(), AipsError);
}

SpectralModel::SpectralModel(const SpectralModel& other) 
  :RecordTransformable(),
   itsRefFreq(other.itsRefFreq),
   itsFreqUnit(other.itsFreqUnit),
   itsFreqErr(other.itsFreqErr)
{
  DebugAssert(SpectralModel::ok(), AipsError);
}

SpectralModel& SpectralModel::operator=(const SpectralModel& other) {
  if (this != &other) {
    itsRefFreq = other.itsRefFreq;
    itsFreqUnit = other.itsFreqUnit;
    itsFreqErr = other.itsFreqErr;
  }
  DebugAssert(SpectralModel::ok(), AipsError);
  return *this;
}

SpectralModel::~SpectralModel() {
}

const String& SpectralModel::ident() const {
  DebugAssert(SpectralModel::ok(), AipsError);
  static String typeString;
  typeString = ComponentType::name(type());
  return typeString;
}

const MFrequency& SpectralModel::refFrequency() const {
  DebugAssert(SpectralModel::ok(), AipsError);
  return itsRefFreq;
}

Double SpectralModel::refFreqInFrame(const MFrequency::Ref& elframe) const {
    const MFrequency& refFreq(refFrequency());
    Double nu0;
    Bool stupidTransform = (elframe.getType() == MFrequency::REST) ||  (elframe.getType() == MFrequency::N_Types) || (refFreq.getRef().getType() == MFrequency::REST) ||  (refFreq.getRef().getType() == MFrequency::N_Types);
  if (elframe != refFreq.getRef() && !stupidTransform) {
    try{
      nu0 = (MFrequency::Convert(refFreq, elframe))().getValue().getValue();
    }
    catch(...){
      throw(AipsError(String("Cannot transform frequency from ") +  MFrequency::showType(elframe.getType()) + String(" to ") + MFrequency::showType(refFreq.getRef().getType())));
    }
  } else {
    nu0 = refFreq.getValue().getValue();
  }

  return nu0;
}

void SpectralModel::setRefFrequency(const MFrequency& newRefFreq) {
  itsRefFreq = newRefFreq;
  DebugAssert(SpectralModel::ok(), AipsError);
}

const Unit& SpectralModel::frequencyUnit() const {  
  DebugAssert(SpectralModel::ok(), AipsError);
  return itsFreqUnit;
}

void SpectralModel::convertFrequencyUnit(const Unit& freqUnit) {
  itsFreqUnit = freqUnit;
  DebugAssert(SpectralModel::ok(), AipsError);
}

void SpectralModel::setRefFrequencyError(const Quantum<Double>& newRefFreqErr){
  if (badError(newRefFreqErr)) {
    LogIO logErr(LogOrigin("SpectralModel", "setRefFrequencyError"));
    logErr << "The errors must be positive values with units like Hz"
	   << LogIO::EXCEPTION;
  }
  itsFreqErr = newRefFreqErr;
  DebugAssert(SpectralModel::ok(), AipsError);
}

const Quantum<Double>& SpectralModel::refFrequencyError() const {
  return itsFreqErr;
}

void  SpectralModel::sample(Vector<Double>& scale, 
			    const Vector<MFrequency::MVType>& frequencies, 
			    const MFrequency::Ref& refFrame) const {
  DebugAssert(SpectralModel::ok(), AipsError);
  const uInt nSamples = frequencies.nelements();
  DebugAssert(scale.nelements() == nSamples, AipsError);
  
  for (uInt i = 0; i < nSamples; i++) {
    scale(i) = sample(MFrequency(frequencies(i), refFrame));
  }
}

Bool SpectralModel::fromRecord(String& errorMessage, 
				  const RecordInterface& record) {
  const String freqString("frequency");
  if (!record.isDefined(freqString)) {
    errorMessage += "The 'frequency' field does not exist\n";
    return false;
  }
  const RecordFieldId frequency(freqString);
  if (record.dataType(frequency) != TpRecord) {
    if ((record.dataType(frequency) == TpString) && 
	(record.shape(frequency) == IPosition(1,1)) &&
	(record.asString(frequency) == String("current"))) {
      return true;
    } else {
      errorMessage += "The 'frequency' field must be a record\n";
      errorMessage += "or the string 'current' (to use the current value)\n";
      return false;
    }
  }
  const Record& freqRecord = record.asRecord(frequency);
  MeasureHolder mh;
  if (!mh.fromRecord(errorMessage, freqRecord)) {
    errorMessage += "Could not parse the reference frequency\n";
    return false;
  }
  if (!mh.isMFrequency()) {
    errorMessage += "The reference frequency is not a frequency measure\n";
    return false;
  }
  SpectralModel::setRefFrequency(mh.asMFrequency());
  return true;
}

Bool SpectralModel::toRecord(String& errorMessage,
			     RecordInterface& record) const {
  DebugAssert(SpectralModel::ok(), AipsError);
  record.define(RecordFieldId("type"), ComponentType::name(type()));
  Record freqRecord;
  const MeasureHolder mh(refFrequency());
  if (!mh.toRecord(errorMessage, freqRecord)) {
    errorMessage += "Could not convert the reference frequency to a record\n";
    return false;
  }
  const String m0String("m0");
  if (freqRecord.isDefined(m0String)) {
    const RecordFieldId m0(m0String);
    if (freqRecord.dataType(m0) == TpRecord) {
      Record m0Rec = freqRecord.asRecord(m0);
      QuantumHolder qh;
      String error;
      if (qh.fromRecord(error, m0Rec) && qh.isQuantumDouble()) {
	Quantum<Double> q = qh.asQuantumDouble();
	q.convert(frequencyUnit());
	qh = QuantumHolder(q);
	if (qh.toRecord(error, m0Rec)) {
	  freqRecord.defineRecord(m0, m0Rec);
	}
      }
    }
  }
  record.defineRecord(RecordFieldId("frequency"), freqRecord);
  return true;
}

ComponentType::SpectralShape SpectralModel::
getType(String& errorMessage, const RecordInterface& record) {
  const String typeString("type");
  if (!record.isDefined(typeString)) {
    errorMessage += 
      String("\nThe record does not have a 'type' field.");
    return ComponentType::UNKNOWN_SPECTRAL_SHAPE;
  }
  const RecordFieldId type(typeString);
  if (record.dataType(type) != TpString) {
    errorMessage += String("\nThe 'type' field in the spectrum record,")
      + String("must be a String");
    return ComponentType::UNKNOWN_SPECTRAL_SHAPE;
  }      
  if (record.shape(type) != IPosition(1,1)) {
    errorMessage += String("The 'type' field, in the spectrum record,") + 
      String(" must have only 1 element\n");
    return ComponentType::UNKNOWN_SPECTRAL_SHAPE;
  }      
  const String& typeVal = record.asString(type);
  return ComponentType::spectralShape(typeVal);
}

Bool SpectralModel::ok() const {
  if (itsFreqUnit != Unit("GHz")) {
    LogIO logErr(LogOrigin("SpectralModel", "ok()"));
    logErr << LogIO::SEVERE << "The reference frequency has units with " 
	   << endl << " different dimensions than the Hz."
           << LogIO::POST;
    return false;
  }
  return true;
}

Bool SpectralModel::badError(const Quantum<Double>& quantum) {
  const UnitVal hz(1, "Hz");
  return !(quantum.check(hz)) || (quantum.getValue() < 0.0);
}

// Local Variables: 
// compile-command: "gmake SpectralModel"
// End: 

} //# NAMESPACE CASA - END

