
#ifndef CPointingMethod_H
#define CPointingMethod_H

/*
 * ALMA - Atacama Large Millimeter casacore::Array
 * (c) European Southern Observatory, 2002
 * (c) Associated Universities Inc., 2002
 * Copyright by ESO (in the framework of the ALMA collaboration),
 * Copyright by AUI (in the framework of the ALMA collaboration),
 * All rights reserved.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY, without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA 02111-1307  USA
 * 
 * /////////////////////////////////////////////////////////////////
 * // WARNING!  DO NOT MODIFY THIS FILE!                          //
 * //  ---------------------------------------------------------  //
 * // | This is generated code!  Do not modify this file.       | //
 * // | Any changes will be lost when the file is re-generated. | //
 * //  ---------------------------------------------------------  //
 * /////////////////////////////////////////////////////////////////
 *
 * casacore::File CPointingMethod.h
 */

#ifndef __cplusplus
#error This is a C++ include file and cannot be used from plain C
#endif

#include <iostream>
#include <string>
#include <vector>
/**
  * A namespace to encapsulate the PointingMethod enumeration.
  */
#ifndef WITHOUT_ACS
#include <almaEnumerations_IFC.h>
#else

// This part mimics the behaviour of 
namespace PointingMethodMod
{
  //! PointingMethod.
  //!  [CalDM.CalPointing] Method of pointing measurement
  
  const char *const revision = "-1";
  const int version = 1;
  
  enum PointingMethod
  { 
    THREE_POINT /*!< Three-point scan */
     ,
    FOUR_POINT /*!< Four-point scan */
     ,
    FIVE_POINT /*!< Five-point scan */
     ,
    CROSS /*!< Cross scan */
     ,
    CIRCLE /*!< Circular scan */
     ,
    HOLOGRAPHY /*!<  */
     
  };
  typedef PointingMethod &PointingMethod_out;
} 
#endif

namespace PointingMethodMod {
	std::ostream & operator << ( std::ostream & out, const PointingMethod& value);
	std::istream & operator >> ( std::istream & in , PointingMethod& value );
}

/** 
  * A helper class for the enumeration PointingMethod.
  * 
  */
class CPointingMethod {
  public:
 
	/**
	  * Enumerators as strings.
	  */  
	
	static const std::string& sTHREE_POINT; /*!< A const string equal to "THREE_POINT".*/
	
	static const std::string& sFOUR_POINT; /*!< A const string equal to "FOUR_POINT".*/
	
	static const std::string& sFIVE_POINT; /*!< A const string equal to "FIVE_POINT".*/
	
	static const std::string& sCROSS; /*!< A const string equal to "CROSS".*/
	
	static const std::string& sCIRCLE; /*!< A const string equal to "CIRCLE".*/
	
	static const std::string& sHOLOGRAPHY; /*!< A const string equal to "HOLOGRAPHY".*/
	

	/**
	  * Return the major version number as an int.
	  * @return an int.
	  */
	  static int version() ;
	  
	  
	  /**
	    * Return the revision as a string.
	    * @return a string
	    *
	    */
	  static std::string revision() ;
	  
	  
     /**
       * Return the number of enumerators declared in PointingMethodMod::PointingMethod.
       * @return an unsigned int.
       */
       static unsigned int size() ;
       
       
    /**
      * Returns an enumerator as a string.
      * @param e an enumerator of PointingMethodMod::PointingMethod.
      * @return a string.
      */
	static std::string name(const PointingMethodMod::PointingMethod& e);
	
	/**
	  * Equivalent to the name method.
	  */
    static std::string toString(const PointingMethodMod::PointingMethod& f) { return name(f); }

	/** 
	  * Returns vector of  all the enumerators as strings. 
	  * The strings are stored in the vector in the same order than the enumerators are declared in the enumeration. 
	  * @return a vector of string.
	  */
     static const std::vector<std::string> names();	 
    
   	
   	// Create a PointingMethod enumeration object by specifying its name.
   	static PointingMethodMod::PointingMethod newPointingMethod(const std::string& name);
   	
   	/*! Return a PointingMethod's enumerator  given a string.
   	  * @param name the string representation of the enumerator.
   	 *  @return a PointingMethodMod::PointingMethod's enumerator.
   	 *  @throws a string containing an error message if no enumerator could be found for this name.
   	 */
 	static PointingMethodMod::PointingMethod literal(const std::string& name);
 	
    /*! Return a PointingMethod's enumerator given an unsigned int.
      * @param i the index of the enumerator in PointingMethodMod::PointingMethod.
      * @return a PointingMethodMod::PointingMethod's enumerator.
      * @throws a string containing an error message if no enumerator could be found for this integer.
      */
 	static PointingMethodMod::PointingMethod from_int(unsigned int i);	
 	

  private:
    /* Not Implemented.  This is a pure static class. */
    CPointingMethod();
    CPointingMethod(const CPointingMethod&);
    CPointingMethod& operator=(const CPointingMethod&);
    
    static std::string badString(const std::string& name) ;
  	static std::string badInt(unsigned int i) ;
  	
};
 
#endif /*!CPointingMethod_H*/
