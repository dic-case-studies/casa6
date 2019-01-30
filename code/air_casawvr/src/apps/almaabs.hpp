/**
   Bojan Nikolic <b.nikolic@mrao.cam.ac.uk>, <bojan@bnikolic.co.uk>
   Initial version November 2009
   Maintained by ESO since 2013.

   \file almaabs.hpp

   Retrievals for ALMA based on absolute WVR measurements only
   
*/
#ifndef _LIBAIR_ALMAABS_HPP__
#define _LIBAIR_ALMAABS_HPP__

#include <vector>
#include <list>
#include <set>
#include <iostream>

#include <boost/scoped_ptr.hpp>
#include <boost/array.hpp>
#include <boost/ptr_container/ptr_list.hpp>

#include "alma_datastruct.h"
#include "antennautils.hpp"

namespace LibAIR2 {

  // Forward declarations
  struct ALMAWVRCharacter;
  class iALMAAbsRet;  
  struct ALMAResBase;
  struct ALMAContRes;
  class dTdLCoeffsBase;
  class InterpArrayData;
  class dTdLCoeffsSingleInterpolated;
  struct ALMARetOpts;

  /**
   */
  class ALMAAbsRet {
    boost::scoped_ptr<iALMAAbsRet> i;

    bool valid;

  public:
    /**
       \param TObs The observed sky temperatures

       \param el Elevation of observation (radian)

       \param WVRChar Characterisation of the WVR used to make this
       measurement
     */
    ALMAAbsRet(const std::vector<double> &TObs,
	       double el,
	       const ALMAWVRCharacter &WVRChar);

    virtual ~ALMAAbsRet();

    /** \brief Get retrieved results, return false if they are invalid
     */
    bool g_Res(ALMAResBase &res); 

  };

  std::ostream & operator<<(std::ostream &os, 
			    const  ALMAAbsInput &i);

  /** \brief Simple list of data inputs for retrieval
   */
  struct ALMAAbsInpL:
    public std::list<ALMAAbsInput> 
  {

  };

  std::ostream & operator<<(std::ostream &os, 
			    const  ALMAAbsInpL &i);

  /** Prepare input data for single retrieval from observation
      mid-point and for WVR on atenna 0 only
   */
  ALMAAbsInpL SimpleSingleI(const InterpArrayData &d);  

  /** 
      Retrieve for single antenna but multiple retrievals in time

      \param n Make n retrievals uniformly distributed in time across
      the data set
      
      \param states Consider only data with one of these state IDs
   */
  ALMAAbsInpL MultipleUniformI(const InterpArrayData &d,
			       size_t n,
			       const std::set<size_t> &states,
			       int refant);  

  /** Prepare input for mid-point of each new field (sequantially)

      \param time The stamps of each field ID 

      \param fb The segments corresponding to field boundaries
      
      \param states Consider only data with one of these state IDs
   */
  ALMAAbsInpL FieldMidPointI(const InterpArrayData &d,
			     const std::vector<std::pair<double, double> >  &fb,
			     const std::set<size_t> &states,
			     int refant);

  

  /**  Carry out the retrieval of coefficients form a list of inputs;
       remove the inputs which have zero Bayesian evidence from the list
   */
  boost::ptr_list<ALMAResBase> doALMAAbsRet(ALMAAbsInpL &il, 
					    std::vector<std::pair<double, double> > &fb,
					    LibAIR2::AntSet &problemAnts);
  

  /** \brief Calculate coefficients for phase correction from inputs
   */ 
  dTdLCoeffsBase *
  ALMAAbsProcessor(const ALMAAbsInpL &inp,
		   boost::ptr_list<ALMAResBase> &r);

  /** Single retrieval for mid point of the observation
   */
  dTdLCoeffsBase * 
  SimpleSingle(const InterpArrayData &d);  


  /** Single retrieval for mid point of the observation
   */
  dTdLCoeffsBase * 
  SimpleSingle(const InterpArrayData &d);

  /** Single retrieval at mid-point, but fitting also for the
      continuum
   */
  dTdLCoeffsBase * 
  SimpleSingleCont(const InterpArrayData &d, int refant);


  /** 
      Separate retrieval for each new field
  */
  dTdLCoeffsSingleInterpolated *
  SimpleMultiple(const std::vector<std::pair<double, double> > &fb,
		 boost::ptr_list<ALMAResBase> &r);

  
  /** Slightly experimental approach for a single-shot simple
      retrieveval of water vapour
   */
  void ALMAAbsContRetrieve(const std::vector<double> &TObs,
			   double el,
			   const ALMAWVRCharacter &WVRChar,
			   ALMAContRes &res,
			   const ALMARetOpts &opts);

  



}

#endif
