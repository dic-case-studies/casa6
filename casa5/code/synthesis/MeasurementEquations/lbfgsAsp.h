//# lbfgsAsp.h: Using LBFGS to define the objective function for Asp
//# Copyright (C) 1996,1997,1998,1999,2000,2002
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
//# Correspondence concerning AIPS++ should be adressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//#
//# $Id$

#ifndef SYNTHESIS_LBFGSASP_H
#define SYNTHESIS_LBFGSASP_H

#include <ms/MeasurementSets/MeasurementSet.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/IPosition.h>
#include <images/Images/ImageInterface.h>
#include <images/Images/PagedImage.h>
#include <images/Images/TempImage.h>

#include <scimath/Mathematics/FFTServer.h>
#include <scimath/Functionals/Gaussian2D.h>


#include <Eigen/Core>

namespace casa { //# NAMESPACE CASA - BEGIN

class AspObjFunc
{
private:
  int n;
  int nX;
  int nY;
  casacore::Matrix<casacore::Float> itsMatDirty;
  casacore::Matrix<casacore::Complex> itsPsfFT;
  casacore::IPosition center;

public:
  AspObjFunc(const casacore::Matrix<casacore::Float>& dirty,
    const casacore::Matrix<casacore::Complex>& psf,
    const casacore::IPosition& positionOptimum) :
    itsMatDirty(dirty),
    itsPsfFT(psf),
    center(positionOptimum)
  {
    nX = itsMatDirty.shape()(0);
    nY = itsMatDirty.shape()(1);
  }

  ~AspObjFunc() = default;

  double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad)
  {
    double fx = 0.0;

    // x[0]: Amplitude, x[1]: scale
    casacore::Matrix<casacore::Float> Asp(nX, nY);
    casacore::Gaussian2D<casacore::Float> gbeam(x[0], center[0], center[1], x[1], 1, 0);
    for (int j = 0; j < nY; ++j)
    {
      for(int i = 0; i < nX; ++i)
      {
        int px = i - nX/2;
        int py = j - nY/2;
        Asp(i,j) = gbeam(px, py);
      }
    }

    casacore::Matrix<casacore::Complex> AspFT;
    casacore::FFTServer<casacore::Float,casacore::Complex> fft(itsMatDirty.shape());
    fft.fft0(AspFT, Asp);

    casacore::Matrix<casacore::Complex> cWork;
    cWork = AspFT * itsPsfFT;
    casacore::Matrix<casacore::Float> AspConvPsf(itsMatDirty.shape());
    fft.fft0(AspConvPsf, cWork, false);
    //std::cout << "AspConvPsf shape  " << AspConvPsf.shape() << std::endl;
    //fft.flip(AspConvPsf, false, false);

    // gradient. 0: amplitude; 1: scale
    casacore::Matrix<casacore::Float> GradAmp(nX, nY);
    casacore::Gaussian2D<casacore::Float> gbeamGradAmp(1, center[0], center[1], x[1], 1, 0);
    for (int j = 0; j < nY; ++j)
    {
      for(int i = 0; i < nX; ++i)
      {
        int px = i - nX/2;
        int py = j - nY/2;
        GradAmp(i,j) = (-2) * gbeamGradAmp(px, py);
      }
    }
    casacore::Matrix<casacore::Float> Grad0 = product(transpose(itsMatDirty), GradAmp);

    casacore::Matrix<casacore::Float> GradScale(nX, nY);
    for (int j = 0; j < nY; ++j)
    {
      for(int i = 0; i < nX; ++i)
        GradScale(i, j) = (-2)*2*(pow(i-center[0],2) + pow(j-center[1],2))*Asp(i,j)/pow(x[1],3);
    }
    casacore::Matrix<casacore::Float> Grad1 = product(transpose(itsMatDirty), GradScale);

    for (int j = 0; j < nY; ++j)
    {
      for(int i = 0; i < nX; ++i)
      {
        fx = fx + abs(double(itsMatDirty(i, j) - AspConvPsf(i,j)));
        grad[0] = grad[0] + Grad0(i,j);
        grad[1] = grad[1] + Grad1(i,j);
      }
    }
    std::cout << "fx " << fx << " AspConvPsf " << AspConvPsf(100,100) << std::endl;

    return fx;
  }

};

} //# NAMESPACE CASA - END

#endif