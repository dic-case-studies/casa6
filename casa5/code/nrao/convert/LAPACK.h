//# LAPACK.h: AIPS++ interface to LAPACK routines
//# Copyright (C) 1994
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

#ifndef NRAO_LAPACK_H
#define NRAO_LAPACK_H

//# Includes
#include <casa/aips.h>

#define NEED_FORTRAN_UNDERSCORES

#if defined(NEED_FORTRAN_UNDERSCORES)
 #define sgetrf sgetrf_
 #define dgetrf	dgetrf_
 #define cgetrf	cgetrf_
 #define zgetrf	zgetrf_
 #define sgetri	sgetri_
 #define dgetri	dgetri_
 #define cgetri	cgetri_
 #define zgetri	zgetri_
 #define ssolvx	ssolvx_
 #define dsolvx	dsolvx_
 #define csolvx	csolvx_
 #define zsolvx	zsolvx_
 #define sblda  sblda_
 #define cndnm  cndnm_
 #define MAIN   MAIN_
#endif

extern "C" {
// LU decomposition
    void sgetrf(int &, int &, float *, int &, int *, int &);
    void dgetrf(int &, int &, double *, int &, int *, int &);
    void cgetrf(int &, int &, void *, int &, int *, int &);
    void zgetrf(int &, int &, void *, int &, int *, int &);
               
// Inverse computation
    void sgetri(int &, float *, int &, const int *, float *, int &, int &);
    void dgetri(int &, double *, int &, const int *, double *, int &, int &);
    void cgetri(int &, void *, int &, const int *, void *, int &, int &);
    void zgetri(int &, void *, int &, const int *, void *, int &, int &);
               
// Solve AX=B  with errors
    void ssolvx(casacore::Int &, casacore::Int &, float *, casacore::Int &, const float *, casacore::Int &, 
		const casacore::Int *, float *, casacore::Int &, float *, casacore::Int &, float *, 
		float *, casacore::Int &, float *, float *, float *, casacore::Int *);
    void dsolvx(casacore::Int &, casacore::Int &, double *, casacore::Int &, const double *, casacore::Int &, 
	        const casacore::Int *, double *, casacore::Int &, double *, casacore::Int &, double *, 
		double *, casacore::Int &, double *, double *, double *, casacore::Int *);
    void csolvx(casacore::Int &, casacore::Int &, void *, casacore::Int &, const void *, casacore::Int &,
		const casacore::Int *, void *, casacore::Int &, void *, casacore::Int &, float *, 
		float *, casacore::Int &, float *, float *, void *, float *);
    void zsolvx(casacore::Int &, casacore::Int &, void *, casacore::Int &, const void *, casacore::Int &, 
		const casacore::Int *, void *, casacore::Int &, void *, casacore::Int &, double *, 
		double *, casacore::Int &, double *, double *, void *, double *);

// test routines from LAPACK rewrapped for easy C++ call
    void sblda(casacore::Int &, casacore::Int &, float &, float &, casacore::Int &, casacore::Int &, float *, int &);
    void cndnm(casacore::Int &, casacore::Int &, float *, casacore::Int &, float &, float &, int &);

    // called by fortran initialization, may not be needed on all machines,
    // in fact, this might cause problems on some machines.
    void MAIN_(); 
};

#endif
