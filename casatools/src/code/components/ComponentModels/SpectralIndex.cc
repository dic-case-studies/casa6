//# SpectralIndex.cc:
//# Copyright (C) 1998,1999,2000,2003
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
//# $Id: SpectralIndex.cc 21292 2012-11-28 14:58:19Z gervandiepen $

#include <components/ComponentModels/SpectralIndex.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Containers/RecordInterface.h>
#include <casacore/casa/Exceptions/Error.h>
#include <casacore/casa/Arrays/IPosition.h>
#include <casacore/casa/Logging/LogIO.h>
#include <casacore/casa/Logging/LogOrigin.h>
#include <casacore/casa/BasicMath/Math.h>
#include <casacore/measures/Measures/MFrequency.h>
#include <casacore/measures/Measures/MCFrequency.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/casa/Quanta/MVFrequency.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/casa/Utilities/Assert.h>
#include <casacore/casa/Utilities/DataType.h>
#include <casacore/casa/BasicSL/String.h>

using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN

SpectralIndex::SpectralIndex()
  :SpectralModel(),
   itsIndex(0.0),
   itsStokesIndex(Vector<Double>(1, 0.0)),
   itsError(0.0)
{
  DebugAssert(ok(), AipsError);
}

SpectralIndex::SpectralIndex(const MFrequency& refFreq, Double exponent)
  :SpectralModel(refFreq),
   itsIndex(exponent),
   itsStokesIndex(Vector<Double>(1, exponent)),
   itsError(0.0)
{
  DebugAssert(ok(), AipsError);
}

SpectralIndex::SpectralIndex(const SpectralIndex& other) 
  :SpectralModel(other),
   itsIndex(other.itsIndex),
   itsStokesIndex(other.itsStokesIndex),
   itsError(other.itsError)
{
  DebugAssert(ok(), AipsError);
}

SpectralIndex::~SpectralIndex() {
  DebugAssert(ok(), AipsError);
}

SpectralIndex& SpectralIndex::operator=(const SpectralIndex& other) {
  if (this != &other) {
    SpectralModel::operator=(other);
    itsIndex = other.itsIndex;
    itsStokesIndex.resize(other.itsStokesIndex.nelements());
    itsStokesIndex=other.itsStokesIndex;
    itsError = other.itsError;
  }
  DebugAssert(ok(), AipsError);
  return *this;
}

ComponentType::SpectralShape SpectralIndex::type() const {
  DebugAssert(ok(), AipsError);
  return ComponentType::SPECTRAL_INDEX;
}

const Double& SpectralIndex::index() const {
   DebugAssert(ok(), AipsError);
   return itsIndex;
}
const Vector<Double>& SpectralIndex::stokesIndex() const {
  DebugAssert(ok(), AipsError);
  return itsStokesIndex;

}
  void SpectralIndex::setIndex(const Double& newIndex) { 
  itsIndex = newIndex;
  itsStokesIndex.resize(1);
  itsStokesIndex[0]=newIndex;
  DebugAssert(ok(), AipsError);
}


  void SpectralIndex::setStokesIndex(const Vector<Double>& newIndex) { 
  itsIndex = newIndex[0];
  itsStokesIndex.resize();
  itsStokesIndex=newIndex;
  DebugAssert(ok(), AipsError);
}

Double SpectralIndex::sample(const MFrequency& centerFreq) const {
  DebugAssert(ok(), AipsError);
  const MFrequency& refFreq(refFrequency());
  const MFrequency::Ref& centerFreqFrame(centerFreq.getRef());
  Double nu0;
  if (centerFreqFrame != refFreq.getRef()) {
    nu0 = MFrequency::Convert(refFreq,centerFreqFrame)().getValue().getValue();
  } else {
    nu0 = refFreq.getValue().getValue();
  }
  if (nu0 <= 0.0) {
    throw(AipsError("SpectralIndex::sample(...) - "
		    "the reference frequency is zero or negative"));
  }
  const Double nu = centerFreq.getValue().getValue();
  return pow(nu/nu0, itsIndex);
}

void SpectralIndex::sampleStokes(
    const MFrequency& centerFreq, Vector<Double>& iquv
) const {
    Vector<Double> scale(4);
    if(itsStokesIndex.nelements()==1){
      scale.resize(1);
      scale(0)=sample(centerFreq);
      iquv(0) *= scale(0);
      return;
    }
    if(iquv.nelements() != 4)
      throw(AipsError("SpectralIndex::sampleStokes(...) 4 stokes parameters expected"));
    Double nu0=refFreqInFrame(centerFreq.getRef());
    
    if (nu0 <= 0.0) {
      throw(AipsError("SpectralIndex::sampleStokes(...) - "
		    "the reference frequency is zero or negative"));
    }
    const Double nu = centerFreq.getValue().getValue();
    iquv[0] *= pow(nu/nu0, itsStokesIndex[0]);
    //u and v get scaled so that fractional linear pol 
    // is scaled by pow(nu/nu0, alpha[1]);
    // as I is scaled by alpha[0]
    //easily shown then that u and q are scaled by pow(nu/nu0, alpha[0]+alpha[1])
    //similarly for scaling of fractional  circular pol
    

    for (uInt k=1; k < 3; ++k){
      iquv[k]*=pow(nu/nu0, itsStokesIndex[0]+itsStokesIndex[1]);
    }
    ///RM type rotation of linear pol. 
    if(itsStokesIndex[2] != 0.0){
      //angle= RM x lambda^2 : itsStokesIndex[2]=RM
      //if
      // Q'= cos(2*angle)*Q - sin(2*angle)*U
      // U'= sin(2*angle) *Q + cos(2*angle)*U
      // Q' and U' preserves fractional linear pol..and rotates it by angle
      Double twoalpha=2*itsStokesIndex[2]*C::c*C::c*(nu0*nu0-nu*nu)/(nu*nu*nu0*nu0);
      Double cos2alph=cos(twoalpha); 
      Double sin2alph=sin(twoalpha);
      //Q'
      Double tempQ=cos2alph*iquv[1]-sin2alph*iquv[2];
      //U'
      iquv[2]=sin2alph*iquv[1]+cos2alph*iquv[2];
      iquv[1]=tempQ;
    }
    iquv[3]*=pow(nu/nu0, itsStokesIndex[0]+itsStokesIndex[3]);
    
}
void SpectralIndex::sample(Vector<Double>& scale, 
			   const Vector<MFrequency::MVType>& frequencies, 
			   const MFrequency::Ref& refFrame) const {
  DebugAssert(ok(), AipsError);
  const uInt nSamples = frequencies.nelements();
  DebugAssert(scale.nelements() == nSamples, AipsError);

  const MFrequency& refFreq(refFrequency());
  Double nu0;
  if (refFrame != refFreq.getRef()) {
    nu0 = MFrequency::Convert(refFreq, refFrame)().getValue().getValue();
  } else {
    nu0 = refFreq.getValue().getValue();
  }
  if (nu0 <= 0.0) {
    throw(AipsError("SpectralIndex::sample(...) - "
		    "the reference frequency is zero or negative"));
  }

  Double nu;
  for (uInt i = 0; i < nSamples; i++) {
    nu = frequencies(i).getValue();
    scale(i) = pow(nu/nu0, itsIndex);
  }
}

void SpectralIndex::sampleStokes(
    Matrix<Double>& iquv, const Vector<MFrequency::MVType>& frequencies, 
    const MFrequency::Ref& refFrame
) const {
    ThrowIf(
        iquv.shape() != IPosition(2, frequencies.size(), 4),
        "Incorrect Matrix shape"
    );
    const auto nSamples = frequencies.nelements();
    const auto nu0 = refFreqInFrame(refFrame); 
    ThrowIf(nu0 <= 0.0, "the reference frequency is zero or negative");

    Double nu;
    //fractional linear and circular pol variation
    //u and v get scaled so that fractional linear pol 
    // is scaled by pow(nu/nu0, alpha[1]);
    // as I is scaled by alpha[0]
    //easily shown then that u and q are scaled by pow(nu/nu0, alpha[0]+alpha[1])
    //similarly for scaling of fractional  circular pol
    for (uInt i=0; i<nSamples; ++i) {
        nu = frequencies(i).getValue();
        iquv(i, 0) *= pow(nu/nu0, itsStokesIndex(0));
        iquv(i, 3) *= pow(nu/nu0, itsStokesIndex(0) + itsStokesIndex(3));
        iquv(i, 1) *= pow(nu/nu0, itsStokesIndex[0] + itsStokesIndex[1]);
        iquv(i, 2) *= pow(nu/nu0, itsStokesIndex[0] + itsStokesIndex[1]);
    }
    ///RM type rotation of linear pol. 
    if(itsStokesIndex[2] != 0.0) {
        //angle= RM x lambda^2 : itsStokesIndex[2]=RM
        //if
        // Q'= cos(2*angle)*Q - sin(2*angle)*U
        // U'= sin(2*angle) *Q + cos(2*angle)*U
        // Q' and U' preserves fractional linear pol..and rotates it by angle
        for (uInt i = 0; i < nSamples; i++) {
            nu = frequencies(i).getValue();
            Double twoalpha = 2*itsStokesIndex[2]*C::c*C::c*(nu0*nu0-nu*nu)/(nu*nu*nu0*nu0);
            Double cos2alph = cos(twoalpha); 
            Double sin2alph = sin(twoalpha);
            //Q'
            Double tempQ = cos2alph*iquv(i, 1) - sin2alph*iquv(i, 2);
            //U'
            iquv(i, 2) = sin2alph*iquv(i, 1) + cos2alph*iquv(i, 2);
            iquv(i, 1) = tempQ;
        }
    }
}
   
SpectralModel* SpectralIndex::clone() const {
  DebugAssert(ok(), AipsError);
  SpectralModel* tmpPtr = new SpectralIndex(*this);
  AlwaysAssert(tmpPtr != 0, AipsError);
  return tmpPtr;
}

uInt SpectralIndex::nParameters() const {
  DebugAssert(ok(), AipsError);
  return 1;
}

void SpectralIndex::setParameters(const Vector<Double>& newSpectralParms) {
  DebugAssert(newSpectralParms.nelements() == nParameters(), AipsError);
  itsIndex = newSpectralParms(0);
  DebugAssert(ok(), AipsError);
}

Vector<Double> SpectralIndex::parameters() const {
  DebugAssert(ok(), AipsError);
  return Vector<Double>(1, itsIndex);
}

void SpectralIndex::setErrors(const Vector<Double>& newSpectralErrs) {
  DebugAssert(newSpectralErrs.nelements() == nParameters(), AipsError);
  if (newSpectralErrs(0) < 0.0) {
    LogIO logErr(LogOrigin("SpectralIndex", "setErrors(...)"));
    logErr << "The errors must be non-negative."
	   << LogIO::EXCEPTION;
  }
  itsError = newSpectralErrs(0);
  DebugAssert(ok(), AipsError);
}

Vector<Double> SpectralIndex::errors() const {
  DebugAssert(ok(), AipsError);
  return Vector<Double>(1, itsError);
}

Bool SpectralIndex::fromRecord(String& errorMessage, 
 			       const RecordInterface& record) {
  if (!SpectralModel::fromRecord(errorMessage, record)) return false;
  if (!record.isDefined(String("index"))) {
    errorMessage += "The 'spectrum' record must have an 'index' field\n";
    return false;
  }
//
  {
     const RecordFieldId index("index");
     const IPosition shape(1,1);
     if (record.shape(index) != shape) {
       errorMessage += "The 'index' field must be a scalar\n";
       return false;
     }
     Double indexVal;
     switch (record.dataType(index)) {
     case TpDouble:
     case TpFloat:
     case TpInt:
       indexVal = record.asDouble(index);
       break;
     default:
       errorMessage += "The 'index' field must be a real number\n";
       return false;
     }
     setIndex(indexVal);
  }
  if(record.isDefined("stokesindex")){
    Vector<Double> tempstokes=record.asArrayDouble("stokesindex");
    if((tempstokes.nelements() != 4) && (tempstokes.nelements() != 1) ){
      errorMessage += "Stokes indices is not of length 1 or 4\n";
      return false;
    }
    itsStokesIndex.resize();
    itsStokesIndex=tempstokes;
  }
//
  {
      Vector<Double> errorVals(1, 0.0);
      if (record.isDefined("error")) {
        const RecordFieldId error("error");
        const IPosition shape(1,1);
        if (record.shape(error) != shape) {
          errorMessage += "The 'error' field must be a scalar\n";
          return false;
        }
        switch (record.dataType(error)) {
        case TpDouble:
        case TpFloat:
        case TpInt:
            errorVals[0] = record.asDouble(error);
          break;
        default:
          errorMessage += "The 'error' field must be a real number\n";
          return false;
        }
     }
//
     setErrors(errorVals);
  }
//
  DebugAssert(ok(), AipsError);
  return true;
}

Bool SpectralIndex::toRecord(String& errorMessage,
 			     RecordInterface& record) const {
  DebugAssert(ok(), AipsError);
  if (!SpectralModel::toRecord(errorMessage, record)) return false;
  record.define("index", index());
  if(itsStokesIndex.nelements() != 0)
    record.define("stokesindex", itsStokesIndex);
  record.define("error", errors()(0));
  return true;
}

Bool SpectralIndex::convertUnit(String& errorMessage,
				const RecordInterface& record) {
  const String fieldString("index");
  if (!record.isDefined(fieldString)) {
    return true;
  }
  const RecordFieldId field(fieldString);
  if (!((record.dataType(field) == TpString) && 
 	(record.shape(field) == IPosition(1,1)) &&
 	(record.asString(field) == ""))) {
    errorMessage += "The 'index' field must be an empty string\n";
    errorMessage += "(and not a vector of strings)\n";
    return false;
  }
  return true;
}

Bool SpectralIndex::ok() const {
  if (!SpectralModel::ok()) return false;
  if (refFrequency().getValue().getValue() <= 0.0) {
    LogIO logErr(LogOrigin("SpectralIndex", "ok()"));
    logErr << LogIO::SEVERE << "The reference frequency is zero or negative!" 
           << LogIO::POST;
    return false;
  }
  if (abs(itsIndex) > 100) {
    LogIO logErr(LogOrigin("SpectralIndex", "ok()"));
    logErr << LogIO::SEVERE << "The spectral index is greater than 100!" 
           << LogIO::POST;
    return false;
  }
  return true;
}

// Local Variables: 
// compile-command: "gmake SpectralIndex"
// End: 

} //# NAMESPACE CASA - END

