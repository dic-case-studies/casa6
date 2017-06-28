/*
 * ReaderInterface.h
 *
 *  Created on: Jan 5, 2016
 *      Author: nakazato
 */

#ifndef SINGLEDISH_FILLER_READERINTERFACE_H_
#define SINGLEDISH_FILLER_READERINTERFACE_H_

// std includes
#include <string>

#include <casacore/casa/Containers/Record.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/measures/Measures/Stokes.h>
#include <singledishfiller/Filler/AntennaRecord.h>
#include <singledishfiller/Filler/DataRecord.h>
#include <singledishfiller/Filler/FieldRecord.h>
#include <singledishfiller/Filler/FillerUtil.h>
#include <singledishfiller/Filler/ObservationRecord.h>
#include <singledishfiller/Filler/ProcessorRecord.h>
#include <singledishfiller/Filler/SourceRecord.h>
#include <singledishfiller/Filler/SpectralWindowRecord.h>

namespace casa { //# NAMESPACE CASA - BEGIN

// NonCopyable Mixin (CRTP)
template<class T>
class NonCopyable {
protected:
  NonCopyable() {
  }
  ~NonCopyable() {
  }
private:
  NonCopyable(NonCopyable const &) {
  }
  T &operator=(T const &) {
  }
};

class ReaderInterface: private NonCopyable<ReaderInterface> {
public:
  ReaderInterface(std::string const &name, bool const is_nro = false) :
    name_(name), is_nro_(is_nro) {
  }

  virtual ~ReaderInterface() {
  }

  std::string const &getName() const {
    return name_;
  }

  bool isNROData() {
    return is_nro_;
  }
  
  virtual int getNROArraySize() {
    return 0;
  }
  virtual int getNRONumBeam() {
    return 0;
  }
  virtual int getNRONumPol() {
    return 0;
  }
  virtual int getNRONumSpw() {
    return 0;
  }

  virtual int getNROArrayBeamId(int /* array_id */) {
    return -1;
  }
  virtual casacore::Stokes::StokesTypes getNROArrayPol(int /* array_id */) {
    return casacore::Stokes::Undefined;
  }
  virtual int getNROArraySpwId(int /* array_id */) {
    return -1;
  }

  virtual casacore::Bool isFloatData() const {
    return false;
  }

  virtual casacore::MDirection::Types getDirectionFrame() const {
    return casacore::MDirection::J2000;
  }

  void initialize() {
    initializeCommon();
    initializeSpecific();
  }

  void finalize() {
    finalizeCommon();
    finalizeSpecific();
  }

  // get number of rows for MAIN table
  virtual size_t getNumberOfRows() = 0;

  // to get OBSERVATION table
  // The method should return true if row entry is available.
  // If it return false, row will be invalid so it should not be used.
  virtual casacore::Bool getObservationRow(sdfiller::ObservationRecord &record) = 0;

  // to get ANTENNA table
  // The method should return true if row entry is available.
  // If it return false, row will be invalid so it should not be used.
  virtual casacore::Bool getAntennaRow(sdfiller::AntennaRecord &record) = 0;

  // to get PROCESSOR table
  // The method should return true if row entry is available.
  // If it return false, row will be invalid so it should not be used.
  virtual casacore::Bool getProcessorRow(sdfiller::ProcessorRecord &record) = 0;

  // to get SOURCE table
  // The method should return true if row entry is available.
  // If it return false, row will be invalid so it should not be used.
  virtual casacore::Bool getSourceRow(sdfiller::SourceRecord &row) = 0;

  // to get FIELD table
  // The method should return true if row entry is available.
  // If it return false, row will be invalid so it should not be used.
  virtual casacore::Bool getFieldRow(sdfiller::FieldRecord &row) = 0;

  // to get SPECTRAL WINDOW table
  // The method should return true if row entry is available.
  // If it return false, row will be invalid so it should not be used.
  virtual casacore::Bool getSpectralWindowRow(sdfiller::SpectralWindowRecord &row) = 0;

  // for DataAccumulator
  virtual casacore::Bool getData(size_t irow, sdfiller::DataRecord &record) = 0;

protected:
  virtual void initializeSpecific() = 0;
  virtual void finalizeSpecific() = 0;

  std::string const name_;
  bool const is_nro_;

private:
  // common initialization/finalization actions
  void initializeCommon() {
  }
  void finalizeCommon() {
  }
};

} //# NAMESPACE CASA - END

#endif /* SINGLEDISH_FILLER_READERINTERFACE_H_ */
