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
#include <chrono>

namespace casa { //# NAMESPACE CASA - BEGIN

class AspObjFunc
{
private:
  int n;
  int nX;
  int nY;
  unsigned int AspLen;
  casacore::Matrix<casacore::Float> itsMatDirty;
  casacore::Matrix<casacore::Complex> itsPsfFT;
  std::vector<casacore::IPosition> center;

public:
  AspObjFunc(const casacore::Matrix<casacore::Float>& dirty,
    const casacore::Matrix<casacore::Complex>& psf,
    const std::vector<casacore::IPosition>& positionOptimum) :
    itsMatDirty(dirty),
    itsPsfFT(psf),
    center(positionOptimum)
  {
    nX = itsMatDirty.shape()(0);
    nY = itsMatDirty.shape()(1);
    AspLen = center.size();
  }

  ~AspObjFunc() = default;

  double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad)
  {
    double fx = 0.0;
    casacore::Matrix<casacore::Float> AspConvPsfSum(itsMatDirty.shape(), (casacore::Float)0.0);

    const int refi = nX/2;
    const int refj = nY/2;

    int minX = nX - 1;
    int maxX = 0;
    int minY = nY - 1;
    int maxY = 0;

    for (unsigned int k = 0; k < AspLen; k ++)
    {
      ////auto start = std::chrono::high_resolution_clock::now();

      if (isnan(x[2*k]) || x[2*k+1] <= 0) // LBFGS encounters convergense issue or scale < 0
      {
        std::cout << "nan? " << x[2*k] << " neg scale? " << x[2*k+1] << std::endl;
        fx = -999.0;
        return fx;
      }

      // generate a gaussian for each Asp in the Aspen set
      // x[0]: Amplitude0,       x[1]: scale0
      // x[2]: Amplitude1,       x[3]: scale1
      // x[2k]: Amplitude(k), x[2k+1]: scale(k+1)
      casacore::Matrix<casacore::Float> Asp(nX, nY);
      Asp = 0.0;

      const double sigma5 = 5 * x[2*k+1] / 2;
      /*const int minI = std::max(0, (int)(refi + center[k][0] - sigma5));
      const int maxI = std::min(nX-1, (int)(refi + center[k][0] + sigma5));
      const int minJ = std::max(0, (int)(refj + center[k][1] - sigma5));
      const int maxJ = std::min(nY-1, (int)(refj + center[k][1] + sigma5));*/
      const int minI = std::max(0, (int)(center[k][0] - sigma5));
      const int maxI = std::min(nX-1, (int)(center[k][0] + sigma5));
      const int minJ = std::max(0, (int)(center[k][1] - sigma5));
      const int maxJ = std::min(nY-1, (int)(center[k][1] + sigma5));

      casacore::Gaussian2D<casacore::Float> gbeam(1.0 / (sqrt(2*M_PI)*x[2*k+1])/*x[2*k]*/, center[k][0], center[k][1], x[2*k+1], 1, 0);
      for (int j = minJ; j <= maxJ; j++)
      {
        for (int i = minI; i <= maxI; i++)
        {
          //const int px = i - refi;
          //const int py = j - refj;
          const int px = i;
          const int py = j;
          Asp(i,j) = gbeam(px, py);
        }
      }

      ////auto stop = std::chrono::high_resolution_clock::now();
      ////auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
      ////std::cout << "minI " << minI << " maxI " << maxI << " minJ " << minJ << " maxJ " << maxJ << std::endl;
      ////std::cout << "LBFGS 1st Asp runtime " << duration.count() << " ms" << std::endl;

      ////start = std::chrono::high_resolution_clock::now();
      casacore::Matrix<casacore::Complex> AspFT;
      casacore::FFTServer<casacore::Float,casacore::Complex> fft(itsMatDirty.shape());
      fft.fft0(AspFT, Asp);

      casacore::Matrix<casacore::Complex> cWork;
      cWork = AspFT * itsPsfFT;
      casacore::Matrix<casacore::Float> AspConvPsf(itsMatDirty.shape(), (casacore::Float)0.0);
      fft.fft0(AspConvPsf, cWork, false);
      ////stop = std::chrono::high_resolution_clock::now();
      ////duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
      ////std::cout << "LBFGS 2st AspFFT runtime " << duration.count() << " ms" << std::endl;
      //std::cout << "AspConvPsf shape  " << AspConvPsf.shape() << std::endl;
      fft.flip(AspConvPsf, false, false); //genie need this?

      // gradient. 0: amplitude; 1: scale
      // generate derivatives of amplitude
      ////start = std::chrono::high_resolution_clock::now();
      Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> GradAmp = Eigen::MatrixXf::Zero(nX, nY);
      Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> GradScale = Eigen::MatrixXf::Zero(nX, nY);
      //casacore::Gaussian2D<casacore::Float> gbeamGradAmp(1, center[k][0], center[k][1], x[2*k+1], 1, 0);

      for (int j = minJ; j <= maxJ; j++)
      {
        for (int i = minI; i <= maxI; i++)
        {
          //const int px = i - refi;
          //const int py = j - refj;
          const int px = i;
          const int py = j;
          // generate derivatives of amplitude
          //GradAmp(i,j) = (-2) * gbeamGradAmp(px, py); // genie
          GradAmp(i,j) = (-2) * Asp(i,j) * sqrt(2*M_PI) * x[2*k+1]; // sanjay: -2*asp/Amp
          // generate derivative of scale
          //GradScale(i,j) = (-2)*2*(pow(i-center[k][0],2) + pow(j-center[k][1],2))*Asp(i,j)/pow(x[2*k+1],3); // genie
          GradScale(i,j) = 2 * (pow(i-center[k][0],2) + pow(j-center[k][1],2)) * Asp(i,j) / x[2*k+1]; //sanjay: 2*Asp*((x-xc)^2 + (y-yc)^2)/scale
        }
      }
      ////stop = std::chrono::high_resolution_clock::now();
      ////duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
      ////std::cout << "LBFGS 3rd Grad runtime " << duration.count() << " ms" << std::endl;

      ////start = std::chrono::high_resolution_clock::now();
      casacore::Bool ddel;
      const casacore::Float *dptr = itsMatDirty.getStorage(ddel);
      float *ddptr = const_cast<float*>(dptr);
      Eigen::MatrixXf M = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(ddptr, nX, nY);

      Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Grad0 = M * GradAmp;

      Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Grad1 = M * GradScale;
      itsMatDirty.freeStorage(dptr, ddel);
      ////stop = std::chrono::high_resolution_clock::now();
      ////duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
      ////std::cout << "LBFGS 4th product runtime " << duration.count() << " ms" << std::endl;

      ////start = std::chrono::high_resolution_clock::now();
      // generate objective function
      // returns the objective function value and gradient evaluated on x
      //std::cout << "before Asp# " << k << ": fx " << fx << " itsMatDirty " << itsMatDirty(0,0) << " AspConvPsf " << AspConvPsf(0,0) << std::endl;
      for (int j = minJ; j <= maxJ; j++)
      {
        for (int i = minI; i <= maxI; i++)
        {
          //fx = fx + abs(double(itsMatDirty(i, j) - AspConvPsf(i,j))); genie: seems wrong
          AspConvPsfSum(i,j) = AspConvPsfSum(i,j) + 0.1 * x[2*k] * AspConvPsf(i,j); //gain*optimumstrength*PsfConvAspen
          grad[2*k] = grad[2*k] + double(Grad0(i,j));
          grad[2*k+1] = grad[2*k+1] + double(Grad1(i,j));
        }
      }

      if (minI < minX)
        minX = minI;
      if (maxI > maxX)
        maxX = maxI;
      if (minJ < minY)
        minY = minJ;
      if (maxJ > maxY)
        maxY = maxJ;


      ////stop = std::chrono::high_resolution_clock::now();
      ////duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
      ////std::cout << "LBFGS 5th sum up runtime " << duration.count() << " ms" << std::endl;
    } // end of Aspen

    ////auto start = std::chrono::high_resolution_clock::now();

    for (int j = minY; j < maxY; ++j)
    {
      for(int i = minX; i < maxX; ++i)
      {
        //std::cout << "after Asp fx " << fx << " double " << double(pow(itsMatDirty(i,j) - AspConvPsfSum(i,j),2)) << " MatDirty " << itsMatDirty(i,j) << " PsfSum " << AspConvPsfSum(i,j) << std::endl;
        fx = fx + double(pow(itsMatDirty(i, j) - AspConvPsfSum(i,j),2));
      }
    }
    ////auto stop = std::chrono::high_resolution_clock::now();
    ////auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    ////std::cout << "minX " << minX << " maxX " << maxX << " minY " << minY << " maxY " << maxY << std::endl;
    ////std::cout << "LBFGS 6th runtime " << duration.count() << " ms" << std::endl;
    //std::cout << "after Asp fx " << fx << " AspConvPsfSum " << AspConvPsfSum(0,0) << std::endl;

    // memory used
    //std::cout << "Memory allocated in lbfgs " << double(casacore::HostInfo::memoryUsed()/1024) << " MB." << std::endl;

    return fx;
  }

};

} //# NAMESPACE CASA - END

#endif