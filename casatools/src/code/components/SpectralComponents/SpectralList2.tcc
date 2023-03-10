//# SpectralList2.cc: Member templates for SpectralList
//# Copyright (C) 2001,2002
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
//# $Id: SpectralList2.tcc 19935 2007-02-27 05:07:40Z Malte.Marquarding $

//# Includes
#include <components/SpectralComponents/SpectralList.h>

#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Exceptions/Error.h>
#include <components/SpectralComponents/SpectralElement.h>

namespace casa { //# NAMESPACE CASA - BEGIN

//# Member templates
template <class MT>
void SpectralList::evaluate(casacore::Vector<MT> &y) const {
  for (casacore::uInt j=0; j<y.nelements(); j++) {
    if (list_p.nelements() > 0) y(j) = (*list_p[0])(j);
    else y(j) = 0;
  };
  for (casacore::uInt i=1; i<list_p.nelements(); i++) {
    for (casacore::uInt j=0; j<y.nelements(); j++) y(j) += (*list_p[i])(j);
  };
}

template <class MT>
void SpectralList::evaluate(casacore::Vector<MT> &y, const casacore::Vector<MT> &x) const {
  y.resize(x.nelements());
  for (casacore::uInt j=0; j<x.nelements(); j++) {
    if (list_p.nelements() > 0) y(j) = (*list_p[0])(x(j));
    else y(j) = 0;
  };
  for (casacore::uInt i=1; i<list_p.nelements(); i++) {
    for (casacore::uInt j=0; j<x.nelements(); j++) y(j) += (*list_p[i])(x(j));
  };
}

template <class MT>
void SpectralList::residual(casacore::Vector<MT> &y) const {
  for (casacore::uInt i=0; i<list_p.nelements(); i++) {
    for (casacore::uInt j=0; j<y.nelements(); j++) y(j) -= (*list_p[i])(j);
  };
}

template <class MT>
void SpectralList::residual(casacore::Vector<MT> &y, const casacore::Vector<MT> &x) const {
  if (x.nelements() != y.nelements()) {
    throw(casacore::AipsError("Unequal lengths in arguments SpectralList::residual"));
  };
  for (casacore::uInt i=0; i<list_p.nelements(); i++) {
    for (casacore::uInt j=0; j<x.nelements(); j++) y(j) -= (*list_p[i])(x(j));
  };
}

} //# NAMESPACE CASA - END

