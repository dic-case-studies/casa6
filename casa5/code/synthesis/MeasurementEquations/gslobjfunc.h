#ifndef SYNTHESIS_CPPSOLVERS_H
#define SYNTHESIS_CPPSOLVERS_H

#include <ms/MeasurementSets/MeasurementSet.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/IPosition.h>
#include <images/Images/ImageInterface.h>
#include <images/Images/PagedImage.h>
#include <images/Images/TempImage.h>

#include <scimath/Mathematics/FFTServer.h>
#include <scimath/Functionals/Gaussian2D.h>

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>

namespace casa { //# NAMESPACE CASA - BEGIN

class Params
{
public:
  Params(void *P)  {ParamsPtr = P;}

  //inline casacore::Matrix<casacore::Float>* dirty() {return static_cast<casacore::Matrix<casacore::Float> *>(ParamsPtr[0]);}
  //inline casacore::Matrix<casacore::Complex>* psfFT() {return (casacore::Matrix<casacore::Complex> *)(ParamsPtr[1]);}
  inline casacore::Matrix<casacore::Float>* dirty() {return (casacore::Matrix<casacore::Float> *)((int *) ParamsPtr)[0];}
  inline casacore::Matrix<casacore::Complex>* psfFT() {return (casacore::Matrix<casacore::Complex> *)((int *) ParamsPtr)[1];}
  //inline casacore::CountedPtr<casacore::Matrix<casacore::Float> > dirty() {return (casacore::Matrix<casacore::Float> *)((int *) ParamsPtr)[0];}
  //inline casacore::CountedPtr<casacore::Matrix<casacore::Complex> > psfFT() {return (casacore::Matrix<casacore::Complex> *)((int *) ParamsPtr)[1];}

private:
  void *ParamsPtr;
};

} // end namespace casa


namespace {

// objective fucntion
double my_f (const gsl_vector *x, void *params)
{

  casacore::Matrix<casacore::Float> *itsMatDirty;
  casacore::Matrix<casacore::Complex> *itsPsfFT;
  //casacore::CountedPtr<casacore::Matrix<casacore::Float> > itsMatDirty;
  //casacore::CountedPtr<casacore::Matrix<casacore::Complex> > itsPsfFT;

  // Extract the params from the void pointer...
  /*casa::Params MyP(params);
  itsMatDirty = MyP.dirty();
  itsPsfFT = MyP.psfFT();*/


  /*
  //
  // Make the model and residual image at the current location on the
  // Chisq surface
  //
  T ChiSq=0,A,S;
  int NX,NY, k;
  NX = (*Pj).size();

  k=0;
  //
  // Load only the active Asps from the minimization machine in the CC
  // list.
  //

  int m;
  m=0;
  for(int j=0;j<NX;j++)
    if ((*Pj)[j] == 1)
    {
      k=j;

      A  = (T)gsl_vector_get(v,m++);
      S  = (T)gsl_vector_get(v,m++);

      (*MComp)[k].setAmp(A);
      (*MComp)[k].setScale(S);
    }

  (*MImg).size(NX,NY);
  (*MImg).assign(0);

  (*MImg) += (*MComp);
  Convo.convolve((*MImg));
  MakeImages(*ModelDI, *DImg, *MImg, *ResImg);



  for (int i=EDGE;i<(NX-EDGE);i++)
    for (int j=EDGE;j<(NY-EDGE);j++)
      ChiSq += (*ResImg)(i,j)*(*ResImg)(i,j);
  */
  double ChiSq = 1e-4;
  return ChiSq;
}

// gradient
void my_df (const gsl_vector *x, void *params, gsl_vector *grad)
{

  casacore::Matrix<casacore::Float> *itsMatDirty;
  casacore::Matrix<casacore::Complex> *itsPsfFT;
  //casacore::CountedPtr<casacore::Matrix<casacore::Float> > itsMatDirty;
  //casacore::CountedPtr<casacore::Matrix<casacore::Complex> > itsPsfFT;

  // Extract the params from the void pointer
  casa::Params MyP(params);
  itsMatDirty = MyP.dirty();
  itsPsfFT = MyP.psfFT();
  //casacore::Matrix<casacore::Float>* itsMatDirty = new casacore::Matrix<casacore::Float>(*MyP.dirty());
  //casacore::Matrix<casacore::Complex>* itsPsfFT = new casacore::Matrix<casacore::Complex>(*MyP.psfFT());

  if (itsMatDirty == NULL)
  	std::cout << "It is null" << std::endl;
  else
  	std::cout << "It is not null " << (*itsMatDirty)(256,256) << std::endl;
  /*int nX = itsMatDirty->shape()(0);
  int nY = itsMatDirty->shape()(1);
  std::cout << " debug size " << nX << " " << nY << std::endl;*/
  //debug
  /*for (int j = 125; j < 130; j++)
  {
    for (int i = 125; i < 130; i++)
    {
      std::cout << "itsMatDirty(" << i << "," << j << ") = " << itsMatDirty[i,j] << std::endl;
      std::cout << "itsPsfFT(" << i << "," << j << ") = " << itsPsfFT[i,j] << std::endl;
    }
  }*/
  /*
  //
  // Recover the pointers to the parameters... :-| (Phew!)
  //
  Params<FTYPE> MyP(params);
  ResImg = MyP.ResImg();
  MComp  = MyP.ModComp();
  Pj     = MyP.CompList();
  ModelDI= MyP.ModelDI();
  PSF    = MyP.PSF();

  ResImg->size(Nx,Ny);

  k=0;
  NP=(*Pj).size();

  int m;
  m=0;
  for(int jj=0;jj<NP;jj++)
    if ((*Pj)[jj] == 1)
      {
	int XRange0,XRange1,YRange0,YRange1;
	Fluxon2D<T> tF;
	k=jj;

	dA = dS = 0.0;
	for (int i=0;i<Nx;i++)
	  for (int j=0;j<Ny;j++)
	  {
	    tF = (*MComp)[k];
	    //
	    // Derivatives are evaluated only for the region bound by (XRange0,YRange0)
	    // and (XRange1, YRange1).
	    // This range as of now is not strictly accurate for derivative computation
	    // (specially for dChiSq/dScale) - should be improved to take into account
	    // the multiplication with the parabolid for the derivative w.r.t scale.
	    //
	    tF.getRange(Nx,XRange0,XRange1,YRange0,YRange1);
	    if (tF.inRange(i,j,XRange0,XRange1,YRange0,YRange1))
	      {
		tF.evalDerivatives(i,j);

		tmp = (*ResImg)(i,j);
		dA  += tmp*tF.dF(tF.DA); // element by element
		dS  += tmp*tF.dF(tF.DS);
	      }
	  }

	gsl_vector_set(grad,m++  ,dA);
	gsl_vector_set(grad,m++,dS);

	k++;
      }*/
  gsl_vector_set(grad, 0, 1e-5);
  gsl_vector_set(grad, 1, 1e-5);
}

void my_fdf (const gsl_vector *v, void *params, double *f, gsl_vector *df)
{
  *f = my_f(v, params); // has to be double in GSL
  my_df(v, params, df);
}

void setupSolver(gsl_multimin_fdfminimizer **Solver,
			   const gsl_multimin_fdfminimizer_type *MinimizerType,
			   gsl_multimin_function_fdf *my_func,
			   int NAspen,
			   void *par[],
			   casacore::Matrix<casacore::Float> *dirty,
			   casacore::Matrix<casacore::Complex> *psfFT
			   )
{
  const int NParams = 2 * NAspen;

  if (*Solver)
  	gsl_multimin_fdfminimizer_free(*Solver);

  *Solver = gsl_multimin_fdfminimizer_alloc(MinimizerType, NParams);

  /*bool ddel;
  casacore::Float *dptr = dirty->getStorage(ddel);
  par[0] = (void *) dptr;
  bool ddel2;
  casacore::Float *dptr2 = psfFT->getStorage(ddel2);
  par[1] = (void *) dptr2;*/ // Kumar's fix but cause compile error

  par[0]         = (void *)dirty;
  par[1]         = (void *)psfFT;
  std::cout << "debug dirty = " << (*dirty)(256,256) << " " << dirty << std::endl;

  my_func->f      = &my_f;
  my_func->df     = &my_df;
  my_func->fdf    = &my_fdf;
  my_func->n      = NParams;
  my_func->params = par;
}

void debug_print(const gsl_multimin_fdfminimizer *s, const int k)
{
    const gsl_vector   *x          = NULL;
    const gsl_vector   *g          = NULL;

    std::cout << "At iteration k = " << k << std::endl;

    g = gsl_multimin_fdfminimizer_gradient(s);
    std::cout << "g = " << gsl_vector_get(g, 0) << " " << gsl_vector_get(g, 1) << std::endl;

    x = gsl_multimin_fdfminimizer_x(s);
    std::cout << "x = " << gsl_vector_get(x, 0) << " " << gsl_vector_get(x, 1) << std::endl;

    std::cout << "f(x) = " << gsl_multimin_fdfminimizer_minimum(s) << std::endl;
}

int findComponent(int NIter, gsl_multimin_fdfminimizer *s)
{
  int iter = 0;
  int status = 0;

  do
  {
	// Make the move!
	status = gsl_multimin_fdfminimizer_iterate(s);
	if (status == GSL_ENOPROG)
		gsl_multimin_fdfminimizer_restart(s);

	status = gsl_multimin_test_gradient(s->gradient, 1E-3);
    debug_print(s, iter);

	iter++;
  } while(status == GSL_CONTINUE && iter < NIter);

  return status;
}

//-----------example-------------//
	void helper()
    {
      double x = 5.0;
      double y = gsl_sf_bessel_J0 (x);
      std::cout << "GSL test header " << x << " " << y << std::endl;
    }


/*template<typename T>
class AspObjFunc : public cppoptlib::BoundedProblem<T> {
private:
  int n;
  int nX;
  int nY;
  unsigned int AspLen;
  casacore::Matrix<casacore::Float> itsMatDirty;
  casacore::Matrix<casacore::Complex> itsPsfFT;
  std::vector<casacore::IPosition> center;
  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> newResidual;

public:
    using typename cppoptlib::Problem<T>::TVector;
	using Superclass = cppoptlib::BoundedProblem<T>;

  AspObjFunc(const casacore::Matrix<casacore::Float>& dirty,
    const casacore::Matrix<casacore::Complex>& psf,
    const std::vector<casacore::IPosition>& positionOptimum, int dim) :
    itsMatDirty(dirty),
    itsPsfFT(psf),
    center(positionOptimum),
    Superclass(dim)
  {
    nX = itsMatDirty.shape()(0);
    nY = itsMatDirty.shape()(1);
    AspLen = center.size();
    newResidual = Eigen::MatrixXf::Zero(nX, nY);
  }

  ~AspObjFunc() = default;

    // objective function
    T value(const TVector &x) {
        T fx = 0.0;

	    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> AmpAspConvPsfSum = Eigen::MatrixXf::Zero(nX, nY);

	    const int refi = nX/2;
	    const int refj = nY/2;

	    int minX = nX - 1;
	    int maxX = 0;
	    int minY = nY - 1;
	    int maxY = 0;

        std::cout << "before fx " << fx << std::endl;
	    // First, get the amp * AspenConvPsf for each Aspen to update the residual
	    for (unsigned int k = 0; k < AspLen; k ++)
	    {

	      if (isnan(x[2*k]) || x[2*k+1] <= 0) // LBFGS encounters convergense issue or scale < 0
	      {
	        std::cout << "nan? " << x[2*k] << " neg scale? " << x[2*k+1] << std::endl;
	        //fx = -999.0;
	        //return fx;
	      }

	      // generate a gaussian for each Asp in the Aspen set
	      // x[0]: Amplitude0,       x[1]: scale0
	      // x[2]: Amplitude1,       x[3]: scale1
	      // x[2k]: Amplitude(k), x[2k+1]: scale(k+1)
	      casacore::Matrix<casacore::Float> Asp(nX, nY);
	      Asp = 0.0;

	      const double sigma5 = 5 * x[2*k+1] / 2;
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

	      / *for (int j = minJ; j <= maxJ; j++)
	      {
	        for (int i = minI; i <= maxI; i++)
	        {* /
	      for (int j = 0; j <= nY-1; j++)
	      {
	        for (int i = 0; i <= nX-1; i++)
	        {
	          const int px = i;
	          const int py = j;

	          Asp(i,j) = (1.0/(sqrt(2*M_PI)*fabs(x[2*k+1])))*exp(-(pow(i-center[k][0],2) + pow(j-center[k][1],2))*0.5/pow(x[2*k+1],2));
	        }
	      }
	      for (int j = 125; j < 130; j++)
	      {
	        for (int i = 125; i < 130; i++)
	        {
	          std::cout << "Asp(" << i << "," << j << ") = " << Asp(i,j) << std::endl;
	        }
	      }
	      std::cout << "peak(Asp) " << max(fabs(Asp)) << " amp " << x[2*k] << " scale " << x[2*k+1] << std::endl;

	      casacore::Matrix<casacore::Complex> AspFT;
	      casacore::FFTServer<casacore::Float,casacore::Complex> fft(itsMatDirty.shape());
	      fft.fft0(AspFT, Asp);

	      casacore::Matrix<casacore::Complex> cWork;
	      cWork = AspFT * itsPsfFT;
	      casacore::Matrix<casacore::Float> AspConvPsf(itsMatDirty.shape(), (casacore::Float)0.0);
	      fft.fft0(AspConvPsf, cWork, false);
	      fft.flip(AspConvPsf, false, false); //need this

	      casacore::Bool ddelc;
	      const casacore::Float *dptrc = AspConvPsf.getStorage(ddelc);
	      float *ddptrc = const_cast<float*>(dptrc);
	      Eigen::MatrixXf MAspConvPsf = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(ddptrc, nX, nY);

	      AmpAspConvPsfSum = AmpAspConvPsfSum + x[2*k] * MAspConvPsf; //optimumstrength*PsfConvAspen

	      //debug
	      for (int j = 125; j < 130; j++)
	      {
	        for (int i = 125; i < 130; i++)
	        {
	          std::cout << "AspConvPsf(" << i << "," << j << ") = " << AspConvPsf(i,j) << std::endl;
	          std::cout << "AmpAspConvPsfSum(" << i << "," << j << ") = " << AmpAspConvPsfSum(i,j) << std::endl;
	        }
	      }
	      AspConvPsf.freeStorage(dptrc, ddelc);
	    } // end get amp * AspenConvPsf

	    // update the residual
	    casacore::Bool ddel;
	    const casacore::Float *dptr = itsMatDirty.getStorage(ddel);
	    float *ddptr = const_cast<float*>(dptr);
	    Eigen::MatrixXf Mdirty = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(ddptr, nX, nY);
	    newResidual = Mdirty - AmpAspConvPsfSum;
	    itsMatDirty.freeStorage(dptr, ddel);

	    // generate objective function
	    // returns the objective function value
	    / *for (int j = minY; j < maxY; ++j)
	    {
	      for(int i = minX; i < maxX; ++i)
	      {* /
	    for (int j = 0; j <= nY-1; j++)
	    {
	      for (int i = 0; i <= nX-1; i++)
	      {
	        fx = fx + double(pow(newResidual(i, j), 2));
	      }
	    }
	    std::cout << "after fx " << fx << std::endl;

	    return fx;
    }

    // derivatives
    void gradient(const TVector &x, TVector &grad) {
        / *grad[0]  = -2 * (3 - x[0]) + 200 * (x[1] - x[0] * x[0]) * (-2 * x[0]);
        grad[1]  =                   200 * (x[1] - x[0] * x[0]);* /

        // gradient. 0: amplitude; 1: scale
	    // returns the gradient evaluated on x
	    for (unsigned int k = 0; k < AspLen; k ++)
	    {
	      casacore::Matrix<casacore::Float> Asp(nX, nY);
	      Asp = 0.0;
	      casacore::Matrix<casacore::Float> dAsp(nX, nY);
	      dAsp = 0.0;

	      const double sigma5 = 5 * x[2*k+1] / 2;
	      const int minI = std::max(0, (int)(center[k][0] - sigma5));
	      const int maxI = std::min(nX-1, (int)(center[k][0] + sigma5));
	      const int minJ = std::max(0, (int)(center[k][1] - sigma5));
	      const int maxJ = std::min(nY-1, (int)(center[k][1] + sigma5));

	      / *for (int j = minJ; j <= maxJ; j++)
	      {
	        for (int i = minI; i <= maxI; i++)
	        {* /
	      for (int j = 0; j <= nY-1; j++)
	      {
	        for (int i = 0; i <= nX-1; i++)
	        {
	          const int px = i;
	          const int py = j;

	          Asp(i,j) = (1.0/(sqrt(2*M_PI)*fabs(x[2*k+1])))*exp(-(pow(i-center[k][0],2) + pow(j-center[k][1],2))*0.5/pow(x[2*k+1],2));
	          dAsp(i,j)= Asp(i,j) * (((pow(i-center[k][0],2) + pow(j-center[k][1],2)) / pow(x[2*k+1],2) - 1) / fabs(x[2*k+1])); // verified by python
	        }
	      }
	      for (int j = 125; j < 130; j++)
	      {
	        for (int i = 125; i < 130; i++)
	        {
	          std::cout << "dAsp(" << i << "," << j << ") = " << dAsp(i,j) << std::endl;
	        }
	      }

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

	      casacore::Matrix<casacore::Float> GradAmp(itsMatDirty.shape(), (casacore::Float)0.0);
	      casacore::Matrix<casacore::Float> GradScale(itsMatDirty.shape(), (casacore::Float)0.0);

	      // reset grad to 0. This is important to get the correct optimization.
	      grad[2*k] = 0.0;
	      grad[2*k+1] = 0.0;
	      std::cout << "before grad " << 2*k << ": " << grad[2*k] << std::endl;
	      std::cout << "before grad " << 2*k+1 << ": " << grad[2*k+1] << std::endl;
	      / *for (int j = minJ; j <= maxJ; j++)
	      {
	        for (int i = minI; i <= maxI; i++)
	        {* /
	      for (int j = 0; j <= nY-1; j++)
	      {
	        for (int i = 0; i <= nX-1; i++)
	        {
	          // generate derivatives of amplitude
	          GradAmp(i,j) = (-2) * newResidual(i,j) * AspConvPsf(i,j);
	          // generate derivative of scale
	          GradScale(i,j) = (-2) * x[2*k] * newResidual(i,j) * dAspConvPsf(i,j);
	          grad[2*k] = grad[2*k] + double(GradAmp(i,j));
	          grad[2*k+1] = grad[2*k+1] + double(GradScale(i,j));
	        }
	      }

	      for (int j = 125; j < 130; j++)
	      {
	        for (int i = 125; i < 130; i++)
	        {
	          std::cout << "newResidual(" << i << "," << j << ") = " << newResidual(i,j) << std::endl;
	          std::cout << "dAspConvPsf(" << i << "," << j << ") = " << dAspConvPsf(i,j) << std::endl;
	          std::cout << "GradAmp(" << i << "," << j << ") = " << GradAmp(i,j) << std::endl;
	          std::cout << "GradScale(" << i << "," << j << ") = " << GradScale(i,j) << std::endl;
	        }
	      }
	      std::cout << "after grad " << 2*k << ": " << grad[2*k] << " amp " << x[2*k] << std::endl;
	      std::cout << "after grad " << 2*k+1 << ": " << grad[2*k+1] << " scale " << x[2*k+1] << std::endl;
	    } // end of derivatives
    }

};*/



} // end namespace

#endif