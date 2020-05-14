//# lbfgsAspZhang.h: Using LBFGS to define the objective function for Asp
//# whose computation is speed  up by Zhang 2016 approach.
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

#ifndef SYNTHESIS_LBFGSASPZHANG_H
#define SYNTHESIS_LBFGSASPZHANG_H

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

class AspZhangObjFunc
{
private:
  int n;
  int nX;
  int nY;
  unsigned int AspLen;
  casacore::Matrix<casacore::Float> itsMatDirty;
  std::vector<casacore::IPosition> center;

public:
  AspZhangObjFunc(const casacore::Matrix<casacore::Float>& dirty,
    const std::vector<casacore::IPosition>& positionOptimum) :
    itsMatDirty(dirty),
    center(positionOptimum)
  {
    nX = itsMatDirty.shape()(0);
    nY = itsMatDirty.shape()(1);
    AspLen = center.size();
  }

  ~AspZhangObjFunc() = default;

  double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad)
  {
    double fx = 0.0;
    casacore::Matrix<casacore::Float> AspSum(itsMatDirty.shape(), (casacore::Float)0.0);

    for (unsigned int k = 0; k < AspLen; k ++)
    {
      // generate a gaussian for each Asp in the Aspen set
      // x[0]: Amplitude0,       x[1]: scale0
      // x[2]: Amplitude1,       x[3]: scale1
      // x[i]: Amplitude(i/2), x[i+1]: scale(i/2)
      casacore::Matrix<casacore::Float> Asp(nX, nY);

      if (isnan(x[2*k])) // LBFGS encounters convergense issue
        return fx;

      casacore::Gaussian2D<casacore::Float> gbeam(x[2*k], center[k][0], center[k][1], x[2*k+1], 1, 0);
      for (int j = 0; j < nY; ++j)
      {
        for(int i = 0; i < nX; ++i)
        {
          int px = i - nX/2;
          int py = j - nY/2;
          Asp(i,j) = gbeam(px, py);
        }
      }


      // gradient. 0: amplitude; 1: scale
      // generate derivative of amplitude
      casacore::Matrix<casacore::Float> GradAmp(nX, nY);
      casacore::Gaussian2D<casacore::Float> gbeamGradAmp(1, center[k][0], center[k][1], x[2*k+1], 1, 0);
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

      // generate derivative of scale
      casacore::Matrix<casacore::Float> GradScale(nX, nY);
      for (int j = 0; j < nY; ++j)
      {
        for(int i = 0; i < nX; ++i)
          GradScale(i, j) = (-2)*2*(pow(i-center[k][0],2) + pow(j-center[k][1],2))*Asp(i,j)/pow(x[2*k+1],3);
      }
      casacore::Matrix<casacore::Float> Grad1 = product(transpose(itsMatDirty), GradScale);

      // generate objective function
      // returns the objective function value and gradient evaluated on x
      for (int j = 0; j < nY; ++j)
      {
        for(int i = 0; i < nX; ++i)
        {

          AspSum(i,j) = AspSum(i,j) + Asp(i,j);
          grad[2*k] = grad[2*k] + Grad0(i,j);
          grad[2*k+1] = grad[2*k+1] + Grad1(i,j);
        }
      }

    } // end of Aspen

    for (int j = 0; j < nY; ++j)
    {
      for(int i = 0; i < nX; ++i)
        fx = fx + pow(double(itsMatDirty(i, j) - AspSum(i,j)), 2);
    }

    return fx;
  }

};

} //# NAMESPACE CASA - END

#endif