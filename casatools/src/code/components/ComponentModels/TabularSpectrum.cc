//# TabularSpectrum.cc:
//# Copyright (C) 2010
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
//# $Id: TabularSpectrum.cc 21292 2012-11-28 14:58:19Z gervandiepen $

#include <components/ComponentModels/TabularSpectrum.h>
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
#include <casacore/measures/Measures/MeasureHolder.h>
#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/Quanta/MVFrequency.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/casa/Utilities/Assert.h>
#include <casacore/casa/Utilities/DataType.h>
#include <casacore/casa/BasicSL/String.h>
#include <casacore/scimath/Mathematics/InterpolateArray1D.h>

using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN

TabularSpectrum::TabularSpectrum()
  :SpectralModel(),
   tabFreqVal_p(0),
   flux_p(0), ival_p(0),qval_p(0), uval_p(0), vval_p(0), referenceFreq_p(0.0), maxFreq_p(0.0), minFreq_p(0.0)
{
  freqRef_p=MFrequency::Ref(MFrequency::LSRK);

  DebugAssert(ok(), AipsError);
}

TabularSpectrum::TabularSpectrum(const MFrequency& refFreq,
                                 const Vector<MFrequency::MVType>& freq,
                                 const Vector<Flux<Double> >& flux,
                                 const MFrequency::Ref& refFrame)
  :SpectralModel(refFreq)
{

  Bool stupidTransform = (refFrame.getType() == MFrequency::REST) ||  (refFrame.getType() == MFrequency::N_Types) || (refFreq.getRef().getType() == MFrequency::REST) ||  (refFreq.getRef().getType() == MFrequency::N_Types);

  if (refFrame.getType() != refFreq.getRef().getType() && !stupidTransform) {
    referenceFreq_p = MFrequency::Convert(refFreq, refFrame)().getValue().get("Hz").getValue();
  } else {
    referenceFreq_p = refFreq.getValue().get("Hz").getValue();
  }
  setValues(freq, flux, refFrame);
  DebugAssert(ok(), AipsError);
}

TabularSpectrum::TabularSpectrum(const TabularSpectrum& other) 
  :SpectralModel(other)
{
  operator=(other);
  DebugAssert(ok(), AipsError);
}

TabularSpectrum::~TabularSpectrum() {
  DebugAssert(ok(), AipsError);
}

TabularSpectrum& TabularSpectrum::operator=(const TabularSpectrum& other) {
  if (this != &other) {
    SpectralModel::operator=(other);
    freqRef_p=other.freqRef_p;
    tabFreqVal_p.resize();
    tabFreqVal_p=other.tabFreqVal_p;
    flux_p.resize();
    flux_p=other.flux_p;
    referenceFreq_p=other.referenceFreq_p;
    maxFreq_p=other.maxFreq_p;
    minFreq_p=other.minFreq_p;
    ival_p=other.ival_p;
    qval_p=other.qval_p;
    uval_p=other.uval_p;
    vval_p=other.vval_p;
  }
  DebugAssert(ok(), AipsError);
  return *this;
}

ComponentType::SpectralShape TabularSpectrum::type() const {
  return ComponentType::TABULAR_SPECTRUM;
}

void TabularSpectrum::values(Vector<MFrequency::MVType>& freq, Vector<Flux<Double> >& flux) const {
   freq.resize(tabFreqVal_p.nelements());
   flux.resize(flux_p.nelements());
   flux=flux_p;
   for (uInt k=0; k < tabFreqVal_p.nelements(); ++k){
     freq(k)=MVFrequency(Quantity(tabFreqVal_p(k), "Hz"));
   }
}

void TabularSpectrum::setValues(const Vector<MFrequency::MVType>& frequencies, const Vector<Flux<Double> >& flux, const MFrequency::Ref& refFrame) { 
  if(flux.nelements() != frequencies.nelements()){
    throw(AipsError("frequencies length is not equal to flux length in TabularSpectrum::setValues"));
  }

  referenceFreq_p=refFreqInFrame(refFrame);

  freqRef_p=refFrame;
  tabFreqVal_p.resize(frequencies.nelements());
  flux_p.resize();
  flux_p=flux;
  ival_p.resize(frequencies.nelements());
  qval_p.resize(frequencies.nelements());
  uval_p.resize(frequencies.nelements());
  vval_p.resize(frequencies.nelements());

  for (uInt k=0; k < frequencies.nelements(); ++k){
    tabFreqVal_p(k)=frequencies(k).get("Hz").getValue();
    //IQUV
    flux_p(k).convertPol(ComponentType::STOKES);
    ival_p(k)=flux_p(k).value(Stokes::I).getValue();
    qval_p(k)=flux_p(k).value(Stokes::Q).getValue();
    uval_p(k)=flux_p(k).value(Stokes::U).getValue();
    vval_p(k)=flux_p(k).value(Stokes::V).getValue();
  }
  maxFreq_p=max(tabFreqVal_p);
  minFreq_p=min(tabFreqVal_p);
  //Just make sure the refVal_p is calculated
  this->setRefFrequency(refFrequency());
 
}
void TabularSpectrum::setRefFrequency(const MFrequency& newRefFreq) {
  SpectralModel::setRefFrequency(newRefFreq);
  referenceFreq_p=refFreqInFrame(freqRef_p);
  Vector<Double> xout(1, referenceFreq_p);
  Vector<Double> scale(1,0.0);
  refVal_p.resize(4);
  refVal_p=0.0;
  Vector<Vector<Double> > iquv(4);
  iquv[0].reference(ival_p);
  iquv[1].reference(qval_p);
  iquv[2].reference(uval_p);
  iquv[3].reference(vval_p);
  if(ival_p.nelements() < 1 || tabFreqVal_p.nelements() != ival_p.nelements())
    throw(AipsError("Values have to be set before referenceFrequency in TabularSpectrum"));
  for (uInt k=0; k < 4; ++k){
    InterpolateArray1D<Double, Double>::interpolate(scale, xout, tabFreqVal_p, iquv[k], InterpolateArray1D<Double, Double>::linear);
    refVal_p[k]=scale[0] != 0.0 ? scale[0] : max(iquv[k]);
    
  }
}

Double TabularSpectrum::sample(const MFrequency& centerFreq) const {
  const MFrequency& refFreq(refFrequency());
  const MFrequency::Ref& centerFreqFrame(centerFreq.getRef());
  Double nu;

  Bool stupidTransform = (centerFreqFrame.getType() == MFrequency::REST) ||  (centerFreqFrame.getType() == MFrequency::N_Types) || (freqRef_p.getType() == MFrequency::REST) ||  (freqRef_p.getType() == MFrequency::N_Types);
  if (centerFreqFrame.getType() != freqRef_p.getType() && !stupidTransform) {
    nu = MFrequency::Convert(centerFreq, freqRef_p)().getValue().get("Hz").getValue();
  } else {
    nu = refFreq.getValue().get("Hz").getValue();
  }
  if (nu < minFreq_p || nu > maxFreq_p) {
    throw(AipsError("TabularSpectrun::sample(...) - "
		    "the frequency requested out of range"));
  }

  Vector<Double> xout(1, referenceFreq_p);
  Vector<Double> scale(1,0.0);
  Double refy=refVal_p[0];
  xout[0]=nu;
  InterpolateArray1D<Double, Double>::interpolate(scale, xout, tabFreqVal_p, ival_p, InterpolateArray1D<Double, Double>::linear);
  

  if(refy !=0.0){
    return scale[0]/refy;
  }
  
  return 0.0 ;
}

void TabularSpectrum::sampleStokes(const MFrequency& centerFreq, Vector<Double>& retval) const {
  const MFrequency& refFreq(refFrequency());
  const MFrequency::Ref& centerFreqFrame(centerFreq.getRef());
  Double nu;
  retval.resize(4);
  retval.set(0.0);
  Bool stupidTransform = (centerFreqFrame.getType() == MFrequency::REST) ||  (centerFreqFrame.getType() == MFrequency::N_Types) || (freqRef_p.getType() == MFrequency::REST) ||  (freqRef_p.getType() == MFrequency::N_Types);
  if (centerFreqFrame.getType() != freqRef_p.getType() && !stupidTransform) {
    nu = MFrequency::Convert(centerFreq, freqRef_p)().getValue().get("Hz").getValue();
  } else {
    nu = refFreq.getValue().get("Hz").getValue();
  }
  if (nu < minFreq_p || nu > maxFreq_p) {
    throw(AipsError("TabularSpectrun::sample(...) - "
		    "the frequency requested out of range"));
  }

  Vector<Double> xout(1, referenceFreq_p);
  Vector<Double> scale(1,0.0);
  xout[0]=nu;
  Vector<Vector<Double> > iquv(4);
  iquv[0].reference(ival_p);
  iquv[1].reference(qval_p);
  iquv[2].reference(uval_p);
  iquv[3].reference(vval_p);
  for (uInt k=0; k < 4; ++k){
    InterpolateArray1D<Double, Double>::interpolate(scale, xout, tabFreqVal_p, iquv[k], InterpolateArray1D<Double, Double>::linear);
    retval(k)=scale(0);
  }
  
}

void TabularSpectrum::sample(Vector<Double>& scale, 
			   const Vector<MFrequency::MVType>& frequencies, 
			   const MFrequency::Ref& refFrame) const {
  const uInt nSamples = frequencies.nelements();

  MFrequency::Convert toThisFrame(refFrame, freqRef_p);
  Vector<Double> nu(frequencies.nelements());
  //try frame conversion only if it is not something stupid...
  //if it is then assume the frequencies are fine as is.
  Bool stupidTransform = (refFrame.getType() == MFrequency::REST) ||  (refFrame.getType() == MFrequency::N_Types) || (freqRef_p.getType() == MFrequency::REST) ||  (freqRef_p.getType() == MFrequency::N_Types);
  if ((refFrame.getType() != freqRef_p.getType()) && !stupidTransform) {
    for(uInt k=0; k < nSamples; ++k){
      nu(k) = toThisFrame(frequencies(k).getValue()).getValue().getValue();
    }
  } else {
    for(uInt k=0; k< nSamples; ++k){
      nu(k) = frequencies(k).getValue();
    }
  }
  /*  Vector<Double> xout(1, referenceFreq_p);
  Vector<Double> refVal(1,0.0);
  InterpolateArray1D<Double, Double>::interpolate(refVal, xout, tabFreqVal_p, ival_p, InterpolateArray1D<Double, Double>::linear);
  scale.resize(nSamples);
  */
  InterpolateArray1D<Double, Double>::interpolate(scale, nu, tabFreqVal_p, ival_p, InterpolateArray1D<Double, Double>::linear);
  if(refVal_p(0) !=0.0){
    for (uInt i = 0; i < nSamples; i++) {
      scale(i) = scale(i)/refVal_p(0);
    }
  }
  else{
    if(max(scale) != 0.0)
      scale /= max(scale);
  }

}

void TabularSpectrum::sampleStokes(
    Matrix<Double>& retvals, const Vector<MFrequency::MVType>& frequencies, 
    const MFrequency::Ref& refFrame
) const {
    ThrowIf(
        retvals.shape() != IPosition(2, frequencies.size(), 4),
        "Incorrect Matrix shape"
    );
    const auto nSamples = frequencies.nelements();
    retvals.set(0.0);
    MFrequency::Convert toThisFrame(refFrame, freqRef_p);
    Vector<Double> nu(frequencies.nelements());
    //try frame conversion only if it is not something stupid...
    //if it is then assume the frequencies are fine as is.
    Bool stupidTransform = (
        refFrame.getType() == MFrequency::REST)
        ||  (refFrame.getType() == MFrequency::N_Types)
        || (freqRef_p.getType() == MFrequency::REST)
        ||  (freqRef_p.getType() == MFrequency::N_Types
    );
    if ((refFrame.getType() != freqRef_p.getType()) && !stupidTransform) {
        for(uInt k=0; k < nSamples; ++k){
            nu(k) = toThisFrame(frequencies(k).getValue()).getValue().getValue();
        }
    }
    else {
        for(uInt k=0; k<nSamples; ++k) {
            nu(k) = frequencies(k).getValue();
        }
    }
    /*  Vector<Double> xout(1, referenceFreq_p);
    Vector<Double> refVal(1,0.0);
    InterpolateArray1D<Double, Double>::interpolate(refVal, xout, tabFreqVal_p, ival_p, InterpolateArray1D<Double, Double>::linear);
    scale.resize(nSamples);
    */
    Vector<Double> scaleone(nSamples);
    Vector<Vector<Double> > iquv(4);
    iquv[0].reference(ival_p);
    iquv[1].reference(qval_p);
    iquv[2].reference(uval_p);
    iquv[3].reference(vval_p);
    for (uInt k=0; k<4; ++k){
        InterpolateArray1D<Double, Double>::interpolate(
            scaleone, nu, tabFreqVal_p, iquv[k],
            InterpolateArray1D<Double, Double>::linear
        );
        for (uInt i=0; i<nSamples; ++i) {
	        retvals(i, k) = scaleone[i];  
        }
    }
}

SpectralModel* TabularSpectrum::clone() const {
  DebugAssert(ok(), AipsError);
  SpectralModel* tmpPtr = new TabularSpectrum(*this);
  AlwaysAssert(tmpPtr != 0, AipsError);
  return tmpPtr;
}

uInt TabularSpectrum::nParameters() const {
  return 0;
}

void TabularSpectrum::setParameters(const Vector<Double>& newSpectralParms) {
  AlwaysAssert(newSpectralParms.nelements() == nParameters(), AipsError);
}

Vector<Double> TabularSpectrum::parameters() const {
  return Vector<Double>(0);
}

void TabularSpectrum::setErrors(const Vector<Double>& newSpectralErrs) {
  AlwaysAssert(newSpectralErrs.nelements() == nParameters(), AipsError);
}

Vector<Double> TabularSpectrum::errors() const {
  return Vector<Double>(0);
}

Bool TabularSpectrum::fromRecord(String& errorMessage, 
 			       const RecordInterface& record) {
  if (!SpectralModel::fromRecord(errorMessage, record)) return false;
  //freqRef
  if (!record.isDefined(String("freqRef"))) {
    errorMessage += "The 'TabularSpectrum' record must have an 'freqRef' field\n";
    return false;
  }
  else{
    Record theTmpMF(record.asRecord("freqRef"));
    MeasureHolder mh;
    if(!mh.fromRecord(errorMessage, theTmpMF))
      return false;
    if(mh.isMFrequency())
      freqRef_p=(mh.asMFrequency().getRef());
    else
      return false;
  }


//tabFreqVal
if (!record.isDefined(String("tabFreqVal"))) {
    errorMessage += "The 'TabularSpectrum' record must have an 'tabFreqVal' field\n";
    return false;
  }
  else{
    tabFreqVal_p.resize();
    tabFreqVal_p=Vector<Double> (record.asArrayDouble("tabFreqVal"));
  }
////ival
 if (!record.isDefined(String("ival"))) {
   errorMessage += "The 'TabularSpectrum' record must have an 'ival' field\n";
    return false;
 }
 else{
    ival_p.resize();
    ival_p=Vector<Double> (record.asArrayDouble("ival"));
    
    qval_p=record.isDefined(String("qval")) ? Vector<Double> (record.asArrayDouble("qval")) : Vector<Double>(ival_p.nelements(), 0.0);
    uval_p=record.isDefined(String("uval")) ? Vector<Double> (record.asArrayDouble("uval")) : Vector<Double>(ival_p.nelements(), 0.0);
    vval_p=record.isDefined(String("vval")) ? Vector<Double> (record.asArrayDouble("vval")) : Vector<Double>(ival_p.nelements(), 0.0);
  }

//referenceFreq
 if (!record.isDefined(String("referenceFreq"))) {
   errorMessage += "The 'TabularSpectrum' record must have an 'referenceFreq' field\n";
    return false;
 }
 else{
   referenceFreq_p=record.asDouble("referenceFreq");
 }
 setRefFrequency(MFrequency(Quantity(referenceFreq_p, "Hz"), freqRef_p));
//maxFreq and minFreq
 if (!record.isDefined(String("maxFreq")) || !record.isDefined(String("minFreq")) ) {
   errorMessage += "The 'TabularSpectrum' record must have a 'maxFreq' and a 'minFreq' field\n";
   return false;
 }
 else{
   maxFreq_p=record.asDouble("maxFreq");
   minFreq_p=record.asDouble("minFreq");
 }

  return true;
}

Bool TabularSpectrum::toRecord(String& errorMessage,
 			     RecordInterface& record) const {
  DebugAssert(ok(), AipsError);
  if (!SpectralModel::toRecord(errorMessage, record)) return false;
  //save frame in a temporary MFrequency object
  MFrequency tmpMF(Quantity(0, "Hz"), freqRef_p);
  MeasureHolder mh(tmpMF);
  Record outRec;
  if(!mh.toRecord(errorMessage, outRec))
    return false;
  record.defineRecord("freqRef",outRec);
  record.define("tabFreqVal", tabFreqVal_p);
  record.define("ival", ival_p);
  record.define("qval", qval_p);
  record.define("uval", uval_p);
  record.define("vval", vval_p);
  record.define("referenceFreq", referenceFreq_p);
  record.define("maxFreq", maxFreq_p);
  record.define("minFreq", minFreq_p);
  return true;
}

Bool TabularSpectrum::convertUnit(String& errorMessage,
				const RecordInterface& record) {
  if (!record.isDefined("freqRef")){
    errorMessage+="Not a tabularSpectrum object";
    return false;
  }
  return true;
}

Bool TabularSpectrum::ok() const {
  if (!SpectralModel::ok()) return false;
  if (refFrequency().getValue().getValue() <= 0.0) {
    LogIO logErr(LogOrigin("TabularSpectrum", "ok()"));
    logErr << LogIO::SEVERE << "The reference frequency is zero or negative!" 
           << LogIO::POST;
    return false;
  }
  return true;
}

// Local Variables: 
// compile-command: "gmake SpectralIndex"
// End: 

} //# NAMESPACE CASA - END

