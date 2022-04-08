#ifndef SYNTHESIS_OBJFUNCCPU_BK_H
#define SYNTHESIS_OBJFUNCCPU_BK_H

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

namespace casa { //# NAMESPACE CASA - BEGIN


//template <typename real>
void dsscfg_cpu(const casacore::Matrix<casacore::Float>& itsMatDirty,
  const casacore::Matrix<casacore::Complex>& itsPsfFT,
  const std::vector<casacore::IPosition>& center, 
  int const& nx, int const& ny, double* x, double& f,
                double* fgrad, double** assist_buffer, int task) {
  double hx = 1.0 / static_cast<double>(nx + 1);
  double hy = 1.0 / static_cast<double>(ny + 1);
  double area = 0.5 * hx * hy;
  //
  //     Compute the standard starting point if task = 'XS'.
  //
  if (task == 'XS') {
    double temp1 = 0.5;
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        double temp = static_cast<double>(std::min(j, ny - j)) * hy;
        int k = nx * j + i;
        x[k] =
            temp1 *
            sqrt(std::min(static_cast<double>(std::min(i, nx - i)) * hx, temp));
      }
    }
    *assist_buffer = new double[nx * ny * 4];
    std::cout << "itsMatDirty(100,100]) " << itsMatDirty(100,100) << std::endl;
    return;
  }
  //
  bool feval = task == 'F' || task == 'FG';
  bool geval = task == 'G' || task == 'FG';
  //
  //     Compute the function if task = 'F', the gradient if task = 'G', or
  //     both if task = 'FG'.
  //
  double fquad = 0.0;
  double fexp = 0.0;

  double* fgrad_r = NULL;
  double* fgrad_t = NULL;
  double* fgrad_l = NULL;
  double* fgrad_b = NULL;
  if (geval) {
    fgrad_r = *assist_buffer;
    fgrad_t = fgrad_r + nx * ny;
    fgrad_l = fgrad_t + nx * ny;
    fgrad_b = fgrad_l + nx * ny;

    memset(fgrad, 0, nx * ny * sizeof(double));
    memset(fgrad_r, 0, nx * ny * sizeof(double));
    memset(fgrad_t, 0, nx * ny * sizeof(double));
    memset(fgrad_l, 0, nx * ny * sizeof(double));
    memset(fgrad_b, 0, nx * ny * sizeof(double));
  }
  //
  //     Computation of the function and the gradient over the lower
  //     triangular elements.  The trapezoidal rule is used to estimate
  //     the integral of the exponential term.
  //
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      int k = nx * j + i;
      double v = 0.0;
      double vr = 0.0;
      double vt = 0.0;
      v = x[k];
      if (i != nx - 1) {
        vr = x[k + 1];
      }
      if (j != ny - 1) {
        vt = x[k + nx];
      }
      double dvdx = (vr - v) / hx;
      double dvdy = (vt - v) / hy;
      double expv = exp(v);
      double expvr = exp(vr);
      double expvt = exp(vt);
      if (feval) {
        fquad += dvdx * dvdx + dvdy * dvdy;
        fexp = fexp - 1.0 * (expv + expvr + expvt) / 3.0;
      }
      if (geval) {
        fgrad[k] += -dvdx / hx - dvdy / hy - 1.0 * expv / 3.0;
        if (i != nx - 1) {
          fgrad_r[k + 1] = dvdx / hx - 1.0 * expvr / 3.0;
        }
        if (j != ny - 1) {
          fgrad_t[k + nx] = dvdy / hy - 1.0 * expvt / 3.0;
        }
      }
      //
      //     Computation of the function and the gradient over the upper
      //     triangular elements.  The trapezoidal rule is used to estimate
      //     the integral of the exponential term.
      //
      double vb = 0.0;
      double vl = 0.0;

      if (j != 0) {
        vb = x[k - nx];
      }
      if (i != 0) {
        vl = x[k - 1];
      }

      dvdx = (v - vl) / hx;
      dvdy = (v - vb) / hy;
      double expvb = exp(vb);
      double expvl = exp(vl);
      expv = exp(v);
      if (feval) {
        fquad += dvdx * dvdx + dvdy * dvdy;
        fexp = fexp - 1.0 * (expvb + expvl + expv) / 3.0;
      }
      if (geval) {
        if (j != 0) {
          fgrad_b[k - nx] = -dvdy / hy - 1.0 * expvb / 3.0;
        }
        if (i != 0) {
          fgrad_l[k - 1] = -dvdx / hx - 1.0 * expvl / 3.0;
        }
        fgrad[k] += dvdx / hx + dvdy / hy - 1.0 * expv / 3.0;
      }
    }
  }

  //
  //     Scale the result.
  //
  if (feval) {
    f = 0.0;
    f = area * (.5 * fquad + fexp);
    std::cout << "dsscfg_cpu after f " << f << std::endl;
    std::cout << "center " << center[0] << std::endl;
  }
  if (geval) {
    for (int k = 0; k < nx * ny; ++k) {
      fgrad[k] =
          (fgrad[k] + fgrad_b[k] + fgrad_l[k] + fgrad_r[k] + fgrad_t[k]) * area;
    }
  }
  //
}



// test CPU mode
//template <typename real>
double test_dsscfg_cpu(const casacore::Matrix<casacore::Float>& itsMatDirty,
    const casacore::Matrix<casacore::Complex>& itsPsfFT,
    const std::vector<casacore::IPosition>& center) {
  const int g_nx = 4;
  const int g_ny = 4;

  // initialize LBFGSB option
  LBFGSB_CUDA_OPTION<double> lbfgsb_options;

  lbfgsbcuda::lbfgsbdefaultoption<double>(lbfgsb_options);
  lbfgsb_options.mode = LCM_NO_ACCELERATION;
  lbfgsb_options.eps_f = static_cast<double>(1e-3);
  lbfgsb_options.eps_g = static_cast<double>(1e-3);
  lbfgsb_options.eps_x = static_cast<double>(1e-3);
  lbfgsb_options.max_iteration = 10;

  // initialize LBFGSB state
  LBFGSB_CUDA_STATE<double> state;
  memset(&state, 0, sizeof(state));
  double* assist_buffer_cpu = nullptr;

  double minimal_f = std::numeric_limits<double>::max();
  // setup callback function that evaluate function value and its gradient
  state.m_funcgrad_callback = [&itsMatDirty,&itsPsfFT, &center, &assist_buffer_cpu, &minimal_f](
                                  double* x, double& f, double* g,
                                  const cudaStream_t& stream,
                                  const LBFGSB_CUDA_SUMMARY<double>& summary) {
    std::cout << "I'm called2" << std::endl;
    dsscfg_cpu(itsMatDirty, itsPsfFT, center, g_nx, g_ny, x, f, g, &assist_buffer_cpu, 'FG');
    if (summary.num_iteration % 100 == 0) {
      std::cout << "CPU iteration " << summary.num_iteration << " F: " << f
                << " x[0] " << x[0] << std::endl;    
    }

    minimal_f = fmin(minimal_f, f);
    return 0;
  };

  // initialize CPU buffers
  int N_elements = g_nx * g_ny;

  double* x = new double[N_elements];
  double* g = new double[N_elements];

  double* xl = new double[N_elements];
  double* xu = new double[N_elements];

  // in this example, we don't have boundaries
  memset(xl, 0, N_elements * sizeof(xl[0]));
  memset(xu, 0, N_elements * sizeof(xu[0]));

  // initialize starting point
  //double f_init = std::numeric_limits<double>::max();
  //dsscfg_cpu(itsMatDirty, itsPsfFT, center, g_nx, g_ny, x, f_init, nullptr, &assist_buffer_cpu, 'XS');
  double hx = 1.0 / static_cast<double>(g_nx + 1);
  double hy = 1.0 / static_cast<double>(g_ny + 1);
  //
  //     Compute the standard starting point if task = 'XS'.
  //

    double temp1 = 0.5;
    for (int j = 0; j < g_ny; ++j) {
      for (int i = 0; i < g_nx; ++i) {
        double temp = static_cast<double>(std::min(j, g_ny - j)) * hy;
        int k = g_nx * j + i;
        x[k] =
            temp1 *
            sqrt(std::min(static_cast<double>(std::min(i, g_nx - i)) * hx, temp));
      }
    }
    assist_buffer_cpu = new double[g_nx * g_ny * 4];
    std::cout << "itsMatDirty(100,100]) " << itsMatDirty(100,100) << std::endl;
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

/*#define INST_HELPER(real)                                                    \
  template void dsscfg_cpu<real>(const casacore::Matrix<casacore::Float>& itsMatDirty, \
                                 int const& nx, int const& ny, real* x,      \
                                 real& f, real* fgrad, real** assist_buffer, \
                                 int task);

INST_HELPER(float);
INST_HELPER(double);*/

} // end namespace casa

#endif // SYNTHESIS_OBJFUNCCPU_H
