// -*- C++ -*-
//# ConvFuncDiskCache.cc: Definition of the ConvFuncDiskCache class
//# Copyright (C) 1997,1998,1999,2000,2001,2002,2003
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
//# $Id$
#ifndef SYNTHESIS_CFDEFS_H
#define SYNTHESIS_CFDEFS_H
#include <casacore/casa/Arrays/Array.h>
#include <casacore/scimath/Mathematics/RigidVector.h>
#include <casacore/casa/aips.h>

namespace casa
{
  namespace CFDefs { //# NAMESPACE CASA - BEGIN
    typedef casacore::Array<casacore::Complex> CFType ;
    typedef casacore::Array<casacore::Double> CFTypeReal ;
    enum CACHETYPE {NOTCACHED=0,DISKCACHE, MEMCACHE};
    enum CFARRAYSHAPE {NXPOS=0,NYPOS,NWPOS,NPOLPOS,NBASEPOS,CFNDIM};
    typedef casacore::RigidVector<casacore::Int,2> MuellerElementType ;
    typedef casacore::Matrix<casacore::Int> MuellerMatrixType ;
    typedef casacore::Vector<casacore::Vector<casacore::Int> > PolMapType;

    enum CFCDirType {USERDIR=0, WORKDIR};

  } //# NAMESPACE CASA - END
}
#endif
