#ifndef SYNTHESIS_OBJFUNCCPU_H
#define SYNTHESIS_OBJFUNCCPU_H

#include <ms/MeasurementSets/MeasurementSet.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/IPosition.h>
#include <images/Images/ImageInterface.h>
#include <images/Images/PagedImage.h>
#include <images/Images/TempImage.h>

#include <scimath/Mathematics/FFTServer.h>
#include <scimath/Functionals/Gaussian2D.h>

#include <cuda_runtime.h>
#include "culbfgsb/culbfgsb.h"

#ifndef isnan
#define isnan(x) std::isnan(x)
#endif

namespace casa { //# NAMESPACE CASA - BEGIN

void objfunc_cpu(const casacore::Matrix<casacore::Float>& itsMatDirty,
    const casacore::Matrix<casacore::Complex>& itsPsfFT,
    const std::vector<casacore::IPosition>& center, 
    double* x, double& f, 
    double* fgrad, double** assist_buffer) 
{

  int nX;
  int nY;
  unsigned int AspLen;
  casacore::Matrix<casacore::Float> newResidual;
  casacore::Matrix<casacore::Float> AspConvPsf;
  casacore::Matrix<casacore::Float> dAspConvPsf;
  casacore::FFTServer<casacore::Float,casacore::Complex> fft;
  casacore::Matrix<casacore::Float> Asp;
  casacore::Matrix<casacore::Float> dAsp;
  *assist_buffer = new double[nX * nY];

  nX = itsMatDirty.shape()(0);
  nY = itsMatDirty.shape()(1);
  AspLen = center.size();
  newResidual.resize(nX, nY);
  AspConvPsf.resize(nX, nY);
  dAspConvPsf.resize(nX, nY);
  fft = casacore::FFTServer<casacore::Float,casacore::Complex>(itsMatDirty.shape());
  Asp.resize(nX, nY);
  dAsp.resize(nX, nY);

  f = 0;
  memset(fgrad, 0, 2 * sizeof(double));


    double amp = 1;

    const int refi = nX/2;
    const int refj = nY/2;

    int minX = nX - 1;
    int maxX = 0;
    int minY = nY - 1;
    int maxY = 0;

    // First, get the amp * AspenConvPsf for each Aspen to update the residual
    for (unsigned int k = 0; k < AspLen; k ++)
    {
        amp = x[2*k];
        double scale = x[2*k+1];
        std::cout << "f: amp " << amp << " scale " << scale << std::endl;

      if (isnan(amp) || scale < 0.4) // GSL scale < 0
      {
        //std::cout << "nan? " << amp << " neg scale? " << scale << std::endl;
        // If scale is small (<0.4), make it 0 scale to utilize Hogbom and save time
        scale = (scale = fabs(scale)) < 0.4 ? 0 : scale;
        //std::cout << "reset neg scale to " << scale << std::endl;

        if (scale <= 0)
          return;
      }

      // generate a gaussian for each Asp in the Aspen set
      // x[0]: Amplitude0,       x[1]: scale0
      // x[2]: Amplitude1,       x[3]: scale1
      // x[2k]: Amplitude(k), x[2k+1]: scale(k+1)
      //casacore::Matrix<casacore::Float> Asp(nX, nY);
      Asp = 0.0;
      dAsp = 0.0;

      const double sigma5 = 5 * scale / 2;
      const int minI = std::max(0, (int)(center[k][0] - sigma5));
      const int maxI = std::min(nX-1, (int)(center[k][0] + sigma5));
      const int minJ = std::max(0, (int)(center[k][1] - sigma5));
      const int maxJ = std::min(nY-1, (int)(center[k][1] + sigma5));

      if (minI < minX)
        minX = minI;
      if (maxI > maxX)
        maxX = maxI;
      if (minJ < minY)
        minY = minJ;
      if (maxJ > maxY)
        maxY = maxJ;

      for (int j = minJ; j <= maxJ; j++)
      {
        for (int i = minI; i <= maxI; i++)
        {
          const int px = i;
          const int py = j;

          Asp(i,j) = (1.0/(sqrt(2*M_PI)*fabs(scale)))*exp(-(pow(i-center[k][0],2) + pow(j-center[k][1],2))*0.5/pow(scale,2));
          dAsp(i,j)= Asp(i,j) * (((pow(i-center[k][0],2) + pow(j-center[k][1],2)) / pow(scale,2) - 1) / fabs(scale)); // verified by python
        }
      }

      casacore::Matrix<casacore::Complex> AspFT;
      fft.fft0(AspFT, Asp);
      casacore::Matrix<casacore::Complex> cWork;
      cWork = AspFT * itsPsfFT;
      fft.fft0(AspConvPsf, cWork, false);
      fft.flip(AspConvPsf, false, false); //need this

      // gradient. 0: amplitude; 1: scale
      // returns the gradient evaluated on x
      casacore::Matrix<casacore::Complex> dAspFT;
      fft.fft0(dAspFT, dAsp);
      casacore::Matrix<casacore::Complex> dcWork;
      dcWork = dAspFT * itsPsfFT;
      fft.fft0(dAspConvPsf, dcWork, false);
      fft.flip(dAspConvPsf, false, false); //need this
    } // end get amp * AspenConvPsf

    // reset grad to 0. This is important to get the correct optimization.
    double dA = 0.0;
    double dS = 0.0;

    // Update the residual using the current residual image and the latest Aspen.
    // Sanjay used, Res = OrigDirty - active-set aspen * Psf, in 2004, instead.
    // Both works but the current approach is simpler and performs well too.
    for (int j = minY; j < maxY; ++j)
    {
      for(int i = minX; i < maxX; ++i)
      {
        newResidual(i, j) = itsMatDirty(i, j) - amp * AspConvPsf(i, j);
        f = f + double(pow(newResidual(i, j), 2));

        // derivatives of amplitude
        dA += double((-2) * newResidual(i,j) * AspConvPsf(i,j));
        // derivative of scale
        dS += double((-2) * amp * newResidual(i,j) * dAspConvPsf(i,j));
      }
    }
    std::cout << "after f " << f << std::endl;

    fgrad[0] = dA;
    fgrad[1] = dS; 
}



// test CPU mode
//template <typename real>
double test_objfunc_cpu(const casacore::Matrix<casacore::Float>& itsMatDirty,
    const casacore::Matrix<casacore::Complex>& itsPsfFT,
    const std::vector<casacore::IPosition>& center) {

  // initialize LBFGSB option
  LBFGSB_CUDA_OPTION<double> lbfgsb_options;

  lbfgsbcuda::lbfgsbdefaultoption<double>(lbfgsb_options);
  lbfgsb_options.mode = LCM_NO_ACCELERATION;
  lbfgsb_options.eps_f = static_cast<double>(1e-8);
  lbfgsb_options.eps_g = static_cast<double>(1e-8);
  lbfgsb_options.eps_x = static_cast<double>(1e-8);
  lbfgsb_options.max_iteration = 1000;

  // initialize LBFGSB state
  LBFGSB_CUDA_STATE<double> state;
  memset(&state, 0, sizeof(state));
  double* assist_buffer_cpu = nullptr;

  double minimal_f = std::numeric_limits<double>::max();
  std::cout << "I'm called1" << std::endl;
  // setup callback function that evaluate function value and its gradient
  state.m_funcgrad_callback = [&itsMatDirty,&itsPsfFT, &center, &assist_buffer_cpu, &minimal_f](
                                  double* x, double& f, double* g,
                                  const cudaStream_t& stream,
                                  const LBFGSB_CUDA_SUMMARY<double>& summary) {
    std::cout << "I'm called2" << std::endl;
    objfunc_cpu(itsMatDirty, itsPsfFT, center, x, f, g, &assist_buffer_cpu); 

    minimal_f = fmin(minimal_f, f);
    return 0;
  };

  // initialize CPU buffers
  int N_elements = 2;

  double* x = new double[N_elements];
  double* g = new double[N_elements];

  double* xl = new double[N_elements];
  double* xu = new double[N_elements];

  // in this example, we don't have boundaries
  memset(xl, 0, N_elements * sizeof(xl[0]));
  memset(xu, 0, N_elements * sizeof(xu[0]));

  // initialize starting point
  x[0] = 1000.0;
  x[1] = 15.0;
  //double f_init = std::numeric_limits<double>::max();
  //dsscfg_cpu(itsMatDirty, itsPsfFT, center, g_nx, g_ny, x, f_init, nullptr, &assist_buffer_cpu, 'XS');

  // initialize number of bounds (0 for this example)
  int* nbd = new int[N_elements];
  memset(nbd, 0, N_elements * sizeof(nbd[0]));

  LBFGSB_CUDA_SUMMARY<double> summary;
  memset(&summary, 0, sizeof(summary));

  // call optimization
  lbfgsbcuda::lbfgsbminimize<double>(N_elements, state, lbfgsb_options, x, nbd,
                                   xl, xu, summary);

  // release allocated memory
  delete[] x;
  delete[] g;
  delete[] xl;
  delete[] xu;
  delete[] nbd;
  delete[] assist_buffer_cpu;

  return minimal_f;
}

} // end namespace casa

#endif // SYNTHESIS_OBJFUNCCPU_H
