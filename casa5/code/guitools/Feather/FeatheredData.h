//# Copyright (C) 2005
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


#ifndef FEATHEREDDATA_H_
#define FEATHEREDDATA_H_

#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/aipstype.h>

namespace casa {

/**
 * casacore::Data structure class to collect related feather data
 * in one location.
 */

class FeatheredData {
public:
	FeatheredData();
	bool isEmpty() const;
	void setU( const casacore::Vector<casacore::Float>& xVal, const casacore::Vector<casacore::Float>& yVal );
	void setV( const casacore::Vector<casacore::Float>& xVal, const casacore::Vector<casacore::Float>& yVal );
	casacore::Vector<casacore::Float> getUX() const;
	casacore::Vector<casacore::Float> getUY() const;
	casacore::Vector<casacore::Float> getVX() const;
	casacore::Vector<casacore::Float> getVY() const;
	virtual ~FeatheredData();
private:
	casacore::Vector<casacore::Float> ux;
	casacore::Vector<casacore::Float> uy;
	casacore::Vector<casacore::Float> vx;
	casacore::Vector<casacore::Float> vy;
};

} /* namespace casa */
#endif /* FEATHEREDDATA_H_ */
