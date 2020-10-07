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

#ifndef SYNTHESIS_LBFGSGAUSSIAN_H
#define SYNTHESIS_LBFGSGAUSSIAN_H

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

class GaussianObjFunc
{
private:
  int n;
  int nX;
  int nY;
  unsigned int AspLen;
  casacore::Matrix<casacore::Float> itsMatDirty;
  casacore::Matrix<casacore::Complex> itsPsfFT;
  std::vector<casacore::IPosition> center;
  casacore::Float itsGain;

public:
  GaussianObjFunc(const casacore::Matrix<casacore::Float>& dirty,
    const casacore::Matrix<casacore::Complex>& psf,
    const std::vector<casacore::IPosition>& positionOptimum,
    const casacore::Float gain) :
    itsMatDirty(dirty),
    itsPsfFT(psf),
    center(positionOptimum),
    itsGain(gain)
  {
    nX = itsMatDirty.shape()(0);
    nY = itsMatDirty.shape()(1);
    AspLen = center.size();
  }

  ~GaussianObjFunc() = default;

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

      //if (isnan(x[2*k]) || x[2*k+1] <= 0) // LBFGS encounters convergense issue or scale < 0
      if (isnan(x[k]) || x[k] <= 0) // LBFGS encounters convergense issue or scale < 0
      {
        std::cout << "nan? " << std::endl;//x[2*k] << " neg scale? " << x[2*k+1] << std::endl;
        fx = -999.0;
        return fx;
      }

      // generate a gaussian for each Asp in the Aspen set
      // x[0]: Amplitude0,       x[1]: scale0
      // x[2]: Amplitude1,       x[3]: scale1
      // x[2k]: Amplitude(k), x[2k+1]: scale(k+1)
      casacore::Matrix<casacore::Float> Asp(nX, nY);
      Asp = 0.0;
      casacore::Matrix<casacore::Float> dAsp(nX, nY);
      dAsp = 0.0;

      //const double sigma5 = 5 * x[2*k+1] / 2;
      const double sigma5 = 5 * x[k] / 2;
      const int minI = std::max(0, (int)(center[k][0] - sigma5));
      const int maxI = std::min(nX-1, (int)(center[k][0] + sigma5));
      const int minJ = std::max(0, (int)(center[k][1] - sigma5));
      const int maxJ = std::min(nY-1, (int)(center[k][1] + sigma5));

      //casacore::Gaussian2D<casacore::Float> gbeam(1.0 / (sqrt(2*M_PI)*x[2*k+1]), center[k][0], center[k][1], x[2*k+1], 1, 0);
      casacore::Gaussian2D<casacore::Float> gbeam(1.0 / (sqrt(2*M_PI)*x[k]), center[k][0], center[k][1], x[k], 1, 0);
      /*for (int j = minJ; j <= maxJ; j++)
      {
        for (int i = minI; i <= maxI; i++)
        {*/
      for (int j = 0; j <= nY-1; j++)
      {
        for (int i = 0; i <= nX-1; i++)
        {
          const int px = i;
          const int py = j;
          Asp(i,j) = gbeam(px, py); //verified by python
          dAsp(i,j)= Asp(i,j) * (((pow(i-center[k][0],2) + pow(j-center[k][1],2)) / pow(x[k],2) - 1) / x[k]); // verified by python
        }
      }
      std::cout << "peak(Asp) " << max(fabs(Asp)) << " scale " << x[k] << std::endl;
      std::cout << "peak(dAsp) " << max(fabs(dAsp)) << " dAsp(128,128) " << dAsp(128,128) << std::endl;
      //debug
      for (int j = 125; j <= 130; j++)
      {
        for (int i = 125; i <= 130; i++)
        {
          std::cout << " Asp(" << i << "," << j << ") = " << Asp(i,j) << std::endl;
          std::cout << " dAsp(" << i << "," << j << ") = " << dAsp(i,j) << std::endl;
        }
      }
      //

      casacore::Matrix<casacore::Complex> AspFT;
      casacore::FFTServer<casacore::Float,casacore::Complex> fft(itsMatDirty.shape());
      fft.fft0(AspFT, Asp);

      casacore::Matrix<casacore::Complex> cWork;
      cWork = AspFT * itsPsfFT;
      casacore::Matrix<casacore::Float> AspConvPsf(itsMatDirty.shape(), (casacore::Float)0.0);
      fft.fft0(AspConvPsf, cWork, false);
      fft.flip(AspConvPsf, false, false); //need this

      casacore::Matrix<casacore::Complex> dAspFT;
      fft.fft0(dAspFT, dAsp);
      casacore::Matrix<casacore::Complex> dcWork;
      dcWork = dAspFT * itsPsfFT;
      casacore::Matrix<casacore::Float> dAspConvPsf(itsMatDirty.shape(), (casacore::Float)0.0);
      fft.fft0(dAspConvPsf, dcWork, false);
      fft.flip(dAspConvPsf, false, false); //need this

      // gradient. 0: amplitude; 1: scale
      // generate derivatives of amplitude
      Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> GradAmp = Eigen::MatrixXf::Zero(nX, nY);
      Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> GradScale = Eigen::MatrixXf::Zero(nX, nY);

      //for (int j = minJ; j <= maxJ; j++)
      //{
      //  for (int i = minI; i <= maxI; i++)
      //  {
      /*for (int j = 0; j <= nY-1; j++)
      {
        for (int i = 0; i <= nX-1; i++)
        {
          const int px = i;
          const int py = j;
          // generate derivatives of amplitude
          //GradAmp(i,j) = (-2) * gbeamGradAmp(px, py); // genie
          //GradAmp(i,j) = (-2) * Asp(i,j) * sqrt(2*M_PI) * x[2*k+1]; // sanjay: -2*asp/Amp
          ////GradAmp(i,j) = (-2) * AspConvPsf(i,j); // try
          // generate derivative of scale
          //GradScale(i,j) = (-2)*2*(pow(i-center[k][0],2) + pow(j-center[k][1],2))*Asp(i,j)/pow(x[2*k+1],3); // genie
          //GradScale(i,j) = 2 * (pow(i-center[k][0],2) + pow(j-center[k][1],2)) * Asp(i,j) / x[2*k+1]; //sanjay: 2*Asp*((x-xc)^2 + (y-yc)^2)/scale
          ////GradScale(i,j) = (-2) * x[2*k] * ((pow(i-center[k][0],2) + pow(j-center[k][1],2)) / pow(x[2*k+1],2) - 1) * (AspConvPsf(i,j) / x[2*k+1]); //try
          GradScale(i,j) = (-2) * ((pow(i-center[k][0],2) + pow(j-center[k][1],2)) / pow(x[k],2) - 1) * (AspConvPsf(i,j) / x[k]);
        }
      }*/

      casacore::Bool ddel;
      const casacore::Float *dptr = itsMatDirty.getStorage(ddel);
      float *ddptr = const_cast<float*>(dptr);
      Eigen::MatrixXf Mdirty = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(ddptr, nX, nY);

      casacore::Bool ddelc;
      const casacore::Float *dptrc = AspConvPsf.getStorage(ddelc);
      float *ddptrc = const_cast<float*>(dptrc);
      Eigen::MatrixXf MAspConvPsf = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(ddptrc, nX, nY);

      casacore::Bool ddelc2;
      const casacore::Float *dptrc2 = dAspConvPsf.getStorage(ddelc2);
      float *ddptrc2 = const_cast<float*>(dptrc2);
      Eigen::MatrixXf MdAspConvPsf = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(ddptrc2, nX, nY);

      Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Mres;
      Mres = Mdirty - MAspConvPsf;
      //debug
      for (int j = 125; j <= 130; j++)
      {
        for (int i = 125; i <= 130; i++)
        {
          std::cout << " Mdirty(" << i << "," << j << ") = " << Mdirty(i,j) << std::endl;
          std::cout << " AspConvPsf(" << i << "," << j << ") = " << AspConvPsf(i,j) << std::endl;
          std::cout << " MAspConvPsf(" << i << "," << j << ") = " << MAspConvPsf(i,j) << std::endl;
          std::cout << " Mres(" << i << "," << j << ") = " << Mres(i,j) << std::endl;
        }
      }
      //
      std::cout << "verify: Mdirty(128,128) " << Mdirty(128,128) << std::endl;
      std::cout << "verify: MAspConvPsf(128,128) " << MAspConvPsf(128,128) << std::endl;
      std::cout << "verify: Mres(128,128) " << Mres(128,128) << std::endl;
      std::cout << "verify: MdAspConvPsf(128,128) " << MdAspConvPsf(128,128) << std::endl;
      ////Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Grad0 = M * GradAmp;
      ////Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Grad1 = M * GradScale;
      //Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Grad1;
      //Grad1 = Mres * GradScale; // this is wrong. Needs to be element by element multiplication
      //GradScale = (-2) * /*amp*/Mres.cwiseProduct(MdAspConvPsf); // this is still wrong. Needs to be element by element multiplication
      itsMatDirty.freeStorage(dptr, ddel);
      AspConvPsf.freeStorage(dptrc, ddelc);

      // generate objective function
      // returns the objective function value and gradient evaluated on x
      std::cout << "before Asp# " << k << ": fx " << fx << " simDirty " << itsMatDirty(128,128) << " AspConvPsf " << AspConvPsf(128,128) << std::endl;
      /*for (int j = minJ; j <= maxJ; j++)
      {
        for (int i = minI; i <= maxI; i++)
        {*/
      grad[k] = 0.0;
      std::cout << "before grad " << k << ": " << grad[k] << std::endl;
      for (int j = 0; j <= nY-1; j++)
      {
        for (int i = 0; i <= nX-1; i++)
        {
          ////AspConvPsfSum(i,j) = AspConvPsfSum(i,j) + /*itsGain * */x[2*k] * AspConvPsf(i,j); //gain*optimumstrength*PsfConvAspen
          AspConvPsfSum(i,j) = AspConvPsfSum(i,j) + AspConvPsf(i,j); //gain*optimumstrength*PsfConvAspen
          ////grad[2*k] = grad[2*k] + double(Grad0(i,j));
          ////grad[2*k+1] = grad[2*k+1] + double(Grad1(i,j));
          //grad[k] = grad[k] + double(Grad1(i,j));
          GradScale(i,j) = (-2) * /*amp*/Mres(i,j) * MdAspConvPsf(i,j);
          grad[k] = grad[k] + double(GradScale(i,j));
        }
      }
      std::cout << "verify: GradScale(128,128) " << GradScale(128,128) << std::endl;
      std::cout << "verify: Mres(127,128) " << Mres(127,128) << std::endl;
      std::cout << "verify: MdAspConvPsf(127,128) " << MdAspConvPsf(127,128) << std::endl;
      std::cout << "verify: GradScale(127,128) " << GradScale(127,128) << std::endl;
      std::cout << "verify: GradScale(128,127) " << GradScale(128,127) << std::endl;
      std::cout << "after grad " << k << ": " << grad[k] << std::endl;

      if (minI < minX)
        minX = minI;
      if (maxI > maxX)
        maxX = maxI;
      if (minJ < minY)
        minY = minJ;
      if (maxJ > maxY)
        maxY = maxJ;

    } // end of Aspen

    /*for (int j = minY; j < maxY; ++j)
    {
      for(int i = minX; i < maxX; ++i)
      {*/
    for (int j = 0; j <= nY-1; j++)
    {
      for (int i = 0; i <= nX-1; i++)
      {
        fx = fx + double(pow(itsMatDirty(i, j) - AspConvPsfSum(i,j),2));
      }
    }

    std::cout << "after simDirty fx " << fx << " AspConvPsfSum " << AspConvPsfSum(128,128) << std::endl;

    return fx;
  }

};

} //# NAMESPACE CASA - END

#endif