//# FluxCalc_SS_JPL_Butler.h: class for getting the fluxes and angular
//# diameters of various Solar System objects.
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

#ifndef COMPONENTS_FLUXCALC_SS_JPL_BUTLER_H
#define COMPONENTS_FLUXCALC_SS_JPL_BUTLER_H

//# include directives
#include <casacore/casa/aips.h>
#include <casacore/casa/BasicSL/String.h>
#include <casacore/measures/Measures/MEpoch.h>
//#include <measures/Measures/MFrequency.h>
#include <casacore/measures/Measures/MDirection.h>
#include <components/ComponentModels/Flux.h>
#include <casacore/tables/Tables/ScalarColumn.h>

namespace casacore{

class MFrequency;
//class ScalarColumn<casacore::Double>;  There doesn't seem to be a way to forward
}

namespace casa { //# NAMESPACE CASA - BEGIN

//declare a template.

// <summary> 
// FluxCalc_SS_JPL_Butler: Compute flux densities and get angular diameters 
// for standard Solar System reference sources.
// </summary>

// <use visibility=export>

// <reviewed reviewer="" date="" tests="" demos="">

// <prerequisite>
// <li><linkto class="FluxStandard">FluxStandard</linkto> module
// </prerequisite>
//
// <etymology>
// From "flux density",
//      "Solar System",
//      "JPL" (ephemeris provides position and angular diameter),
//      and (Bryan) "Butler" (provides model to convert the above to a flux density).
// </etymology>
//
// <synopsis>
// This class organizes a set of functions to compute expected flux densities
// and angular diameters for several Solar System sources commonly used for
// flux calibration in (sub)mm astronomy.
// </synopsis>
//
// <example>
// <srcblock>
// </srcblock>
// </example>
//
// <motivation>
// Make available and encapsulate the Butler Solar System flux density models.
// </motivation>
//
// <todo asof="10/07/07">
// <li> Support more complicated (image) models.
// <li> Support a choice of models for a given object.
// <li> Refine the error estimate by looking at the extremes of the time and
//      frequency intervals.
// </todo>

class FluxCalc_SS_JPL_Butler
{
 public:
  // Default constructor; provided mainly so an array of these can be made more
  // easily.  (Not that there appears to be an immediate need for that.)
  FluxCalc_SS_JPL_Butler();

  // The more useful constructor.
  FluxCalc_SS_JPL_Butler(const casacore::String& objname, const casacore::MEpoch& time);
  
  ~FluxCalc_SS_JPL_Butler();

  // Chooses an object, and returns whether was successful.  It ignores the
  // case of objname.
  casacore::Bool setObj(const casacore::String& objname);
  
  void setTime(const casacore::MEpoch& time);
  //void setFreq(const casacore::MFrequency& freq);

  // These return whether or not the item has been set, and if has, copy their
  // item to the arg.
  casacore::Bool getName(casacore::String& output) const;
  casacore::Bool getTime(casacore::MEpoch& output) const;
  //casacore::Bool getFreq(casacore::MFrequency& output) const;

  // Returns the default direction (N. Pole) on failure.
  casacore::MDirection getDirection();

  // Sets angdiam to the source's angular diameter in radians and returns the
  // object's component type (i.e. DISK), or UNKNOWN_SHAPE on failure.
  ComponentType::Shape getShape(casacore::Double& angdiam);

  // Returns the distance from the object to the Sun in AU, or -1 if unknown.
  casacore::Double getHeliocentricDist();

  // returns the number of objects supported by this class.
  casacore::uInt n_known() const;
  
  // Compute the flux densities, their uncertainties, and the angular diameter
  // of the calibration source for a set of frequencies.  It will try to read
  // the right table in data/ephemerides/JPL-Horizons/.
  //
  // Return value:
  //    the model's shape, or UNKNOWN_SHAPE on failure.
  // Inputs:
  //    must be already set, or it returns UNKNOWN_SHAPE.
  //    report: If false, suppress logger messages like e.g. references for the
  //            temperature.  Intended for use in a loop, i.e. print the
  //            messages for the 1st spw but do not repeat them for the others.
  // Output args:
  //    value: the calculated flux.
  //    error: the estimated uncertainty of value.
  //    angdiam: angular diameter in radians
  ComponentType::Shape compute(casacore::Vector<Flux<casacore::Double> >& values,
                               casacore::Vector<Flux<casacore::Double> >& errors, casacore::Double& angdiam,
                               const casacore::Vector<casacore::MFrequency>& mfreqs,
			       const casacore::Bool report=true);

  // Same as the above, but for a set of sets of frequencies, i.e.
  // mfreqs[spw] is a set of frequencies for channels in spectral window spw.
  ComponentType::Shape compute(casacore::Vector<casacore::Vector<Flux<casacore::Double> > >& values,
                               casacore::Vector<casacore::Vector<Flux<casacore::Double> > >& errors,
                               casacore::Double& angdiam,
                               const casacore::Vector<casacore::Vector<casacore::MFrequency> >& mfreqs);
 private:
  enum KnownObjects {
    Mercury = 0,
    Venus,
    // Earth, // Too highly resolved
    Mars,
    Jupiter,
    Io,
    Ganymede,
    Europa,
    Callisto,
    // Saturn, // Modeling the rings is too complicated.
    Titan,
    Uranus,
    Neptune,
    Triton,
    Pluto,
    Ceres,
    Pallas,
    Vesta,
    Juno,
    Victoria,
    Davida,
    N_KNOWN
  };

  // Tries to set objnum_p to the KnownObject matching name_p.
  // Returns whether or not it found a match.
  casacore::Bool setObjNum();

  // Reads a JPL-Horizons ephemeris table to get temperature_p,
  // r_p (heliocentric distance), delta_p (geocentric distance),
  // and phang_p (phase angle).
  // Returns whether or not it was able to get the info.
  casacore::Bool readEphem();
  
  // Compute the flux densities assuming a uniform disk blackbody,
  // and their uncertainties.
  // It does not check whether everything is setup, since it assumes
  // compute() already took care of it!
  void compute_BB(casacore::Vector<Flux<casacore::Double> >& values,
                  casacore::Vector<Flux<casacore::Double> >& errors, const casacore::Double angdiam,
                  const casacore::Vector<casacore::MFrequency>& mfreqs);

  // Like compute_BB(), except it uses a graybody model where each frequency
  // has a corresponding temperature.
  void compute_GB(casacore::Vector<Flux<casacore::Double> >& values,
                  casacore::Vector<Flux<casacore::Double> >& errors, const casacore::Double angdiam,
                  const casacore::Vector<casacore::MFrequency>& mfreqs,
                  const casacore::Vector<casacore::Double>& temps);

  void compute_venus(casacore::Vector<Flux<casacore::Double> >& values,
		     casacore::Vector<Flux<casacore::Double> >& errors, const casacore::Double angdiam,
		     const casacore::Vector<casacore::MFrequency>& mfreqs);
  void compute_jupiter(casacore::Vector<Flux<casacore::Double> >& values,
                       casacore::Vector<Flux<casacore::Double> >& errors, const casacore::Double angdiam,
                       const casacore::Vector<casacore::MFrequency>& mfreqs);
  void compute_uranus(casacore::Vector<Flux<casacore::Double> >& values,
                      casacore::Vector<Flux<casacore::Double> >& errors, const casacore::Double angdiam,
                      const casacore::Vector<casacore::MFrequency>& mfreqs);
  void compute_neptune(casacore::Vector<Flux<casacore::Double> >& values,
                       casacore::Vector<Flux<casacore::Double> >& errors, const casacore::Double angdiam,
                       const casacore::Vector<casacore::MFrequency>& mfreqs);
  void compute_pluto(casacore::Vector<Flux<casacore::Double> >& values,
                     casacore::Vector<Flux<casacore::Double> >& errors, const casacore::Double angdiam,
                     const casacore::Vector<casacore::MFrequency>& mfreqs);

  // Uses objnum_p to look up a mean temperature, and uses that.
  // report: If false, suppress logger messages like references for the
  // temperature.  Intended for use in a loop, i.e. print the
  // messages for the 1st spw but do not repeat them for the others.
  // Returns whether or not it was successful.
  casacore::Bool compute_constant_temperature(casacore::Vector<Flux<casacore::Double> >& values,
				    casacore::Vector<Flux<casacore::Double> >& errors,
				    const casacore::Double angdiam,
				    const casacore::Vector<casacore::MFrequency>& mfreqs,
				    const casacore::Bool report=true);  

  // Find the row in mjd closest to time_p, and the rows just before and after
  // it, taking boundaries into account.
  casacore::Bool get_row_numbers(casacore::uInt& rowbef, casacore::uInt& rowclosest, casacore::uInt& rowaft,
		       const casacore::ScalarColumn<casacore::Double>& mjd);

  // Put a quadratic, linear, or nearest neighbor interpolation of colname into
  // val.  Returns whether or not it did it.
  // verbose: Send a message to the logger if nearest neighbor is used.
  static casacore::Bool get_interpolated_value(casacore::Double& val, const casacore::String& colname,
                                     const casacore::Table& tab, const casacore::uInt rowbef,
                                     const casacore::uInt rowclosest, const casacore::uInt rowaft,
                                     const casacore::Double f, const casacore::Double dt,
                                     const casacore::Double tp1mt0, const casacore::Double t0mtm1,
                                     const casacore::Bool verbose=true);

  // casacore::Data members which are initialized in the c'tor's initialization list:
  casacore::String     name_p;
  casacore::Bool       hasName_p;
  casacore::MEpoch     time_p;
  casacore::Bool       hasTime_p;
  // casacore::MFrequency freq_p;
  // casacore::Bool       hasFreq_p;
  casacore::Bool       hasEphemInfo_p;
  casacore::Unit       hertz_p;           // make it static and/or const?
  casacore::Bool       has_ra_p;
  casacore::Bool       has_dec_p;
  casacore::Bool       has_illu_p;

  // These are also initialized by the c'tor, but not in the initialization
  // list:
  FluxCalc_SS_JPL_Butler::KnownObjects objnum_p;
  casacore::Bool                                 hasObjNum_p;

  // casacore::Data members that are not initialized by the c'tor:
  casacore::Double temperature_p; // in K
  casacore::Double mean_rad_p;    // in AU (really!)
  casacore::Double r_p;           // heliocentric distance in AU
  casacore::Double delta_p;       // geocentric distance in AU
  casacore::Double phang_p;       // Phase angle in radians.
  casacore::Bool   has_r_p;
  casacore::Double illu_p;        // Illumination, as a fraction.
  casacore::Double ra_p;          // RA, in deg.
  casacore::Double dec_p;         // Declination, in deg.
};


} //# NAMESPACE CASA - END

#endif
