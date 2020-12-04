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

#include <Eigen/Core>

namespace casa { //# NAMESPACE CASA - BEGIN

class ParamObj
{
private:
  int nX;
  int nY;
  unsigned int AspLen;
  casacore::Matrix<casacore::Float> itsMatDirty;
  casacore::Matrix<casacore::Complex> itsPsfFT;
  std::vector<casacore::IPosition> center;
  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> newResidual;

public:
  ParamObj(const casacore::Matrix<casacore::Float>& dirty,
    const casacore::Matrix<casacore::Complex>& psf,
    const std::vector<casacore::IPosition>& positionOptimum) :
    itsMatDirty(dirty),
    itsPsfFT(psf),
    center(positionOptimum)
  {
    nX = itsMatDirty.shape()(0);
    nY = itsMatDirty.shape()(1);
    AspLen = center.size();
    newResidual = Eigen::MatrixXf::Zero(nX, nY);
  }

  ~ParamObj() = default;

  casacore::Matrix<casacore::Float>  getterDirty() { return itsMatDirty; }
  casacore::Matrix<casacore::Complex> getterPsfFT() { return itsPsfFT; }
  std::vector<casacore::IPosition> getterCenter() {return center;}
  unsigned int getterAspLen() { return AspLen; }
  int getterNX() { return nX; }
  int getterNY() { return nY; }
  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> getterRes() { return newResidual; }
  void setterRes(const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& res) { newResidual = res; }

};

} // end namespace casa


namespace {

// objective fucntion
double my_f (const gsl_vector *x, void *params)
{
    // retrieve params for GSL bfgs optimization
    casa::ParamObj *MyP = (casa::ParamObj *) params; //re-cast back to ParamObj to retrieve images
    casacore::Matrix<casacore::Float> itsMatDirty(MyP->getterDirty());
    casacore::Matrix<casacore::Complex> itsPsfFT(MyP->getterPsfFT());
    std::vector<casacore::IPosition> center = MyP->getterCenter();
    const unsigned int AspLen = MyP->getterAspLen();
    const int nX = MyP->getterNX();
    const int nY = MyP->getterNY();
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> newResidual = MyP->getterRes();

	double fx = 0.0;

	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> AmpAspConvPsfSum = Eigen::MatrixXf::Zero(nX, nY);

	const int refi = nX/2;
	const int refj = nY/2;

	int minX = nX - 1;
	int maxX = 0;
	int minY = nY - 1;
	int maxY = 0;

	//std::cout << "before fx " << fx << std::endl;

	// First, get the amp * AspenConvPsf for each Aspen to update the residual
	for (unsigned int k = 0; k < AspLen; k ++)
	{
      double amp = gsl_vector_get(x, 2*k);
      double scale = gsl_vector_get(x, 2*k+1);

	  if (isnan(amp) || scale < 0.4) // GSL scale < 0
	  {
	    std::cout << "nan? " << amp << " neg scale? " << scale << std::endl;
	    scale = (scale = fabs(scale)) < 0.4 ? 0.4 : scale;
	    std::cout << "reset neg scale to " << scale << std::endl;
	    //fx = -999.0;
	    //return fx;
	  }

	  // generate a gaussian for each Asp in the Aspen set
	  // x[0]: Amplitude0,       x[1]: scale0
	  // x[2]: Amplitude1,       x[3]: scale1
	  // x[2k]: Amplitude(k), x[2k+1]: scale(k+1)
	  casacore::Matrix<casacore::Float> Asp(nX, nY);
	  Asp = 0.0;

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
	  /*for (int j = 0; j <= nY-1; j++)
	  {
	    for (int i = 0; i <= nX-1; i++)
	    {*/
	      const int px = i;
	      const int py = j;

	      Asp(i,j) = (1.0/(sqrt(2*M_PI)*fabs(scale)))*exp(-(pow(i-center[k][0],2) + pow(j-center[k][1],2))*0.5/pow(scale,2));
	      //Asp(i,j) = (1.0/(2*M_PI*pow(scale,2)))*exp(-(pow(i-center[k][0],2) + pow(j-center[k][1],2))*0.5/pow(scale,2));
	    }
	  }
	  std::cout << " amp " << amp << " scale " << scale << std::endl;

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
	  Eigen::MatrixXf MAspConvPsf = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>(ddptrc, nX, nY);

	  AmpAspConvPsfSum = AmpAspConvPsfSum + amp * MAspConvPsf; //optimumstrength*PsfConvAspen

	  AspConvPsf.freeStorage(dptrc, ddelc);
	} // end get amp * AspenConvPsf

	// update the residual
	casacore::Bool ddel;
	const casacore::Float *dptr = itsMatDirty.getStorage(ddel);
	float *ddptr = const_cast<float*>(dptr);
	Eigen::MatrixXf Mdirty = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>(ddptr, nX, nY);
	newResidual = Mdirty - AmpAspConvPsfSum;
	itsMatDirty.freeStorage(dptr, ddel);

    // update newResidual back to the ParamObj
    MyP->setterRes(newResidual);

	// generate ChiSq
	// returns the objective function value
	for (int j = minY; j < maxY; ++j)
	{
	  for(int i = minX; i < maxX; ++i)
	  {
	/*for (int j = 0; j <= nY-1; j++)
	{
	  for (int i = 0; i <= nX-1; i++)
	  {*/
	    fx = fx + double(pow(newResidual(i, j), 2));
	  }
	}
	std::cout << "after fx " << fx << std::endl;

	return fx;

  /*
  // Make the model and residual image at the current location on the
  // Chisq surface

  // Load only the active Asps from the minimization machine in the CC
  // list. Res = OrigDirty - active-set aspen * Psf

  for (int i=EDGE;i<(NX-EDGE);i++)
    for (int j=EDGE;j<(NY-EDGE);j++)
      ChiSq += (*ResImg)(i,j)*(*ResImg)(i,j);
 */
}

// gradient
void my_df (const gsl_vector *x, void *params, gsl_vector *grad)
{
    casa::ParamObj *MyP = (casa::ParamObj *) params; //re-cast back to ParamObj to retrieve images
    casacore::Matrix<casacore::Float> itsMatDirty(MyP->getterDirty());
    casacore::Matrix<casacore::Complex> itsPsfFT(MyP->getterPsfFT());
    std::vector<casacore::IPosition> center = MyP->getterCenter();
    const unsigned int AspLen = MyP->getterAspLen();
    const int nX = MyP->getterNX();
    const int nY = MyP->getterNY();
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> newResidual = MyP->getterRes();

    // gradient. 0: amplitude; 1: scale
	// returns the gradient evaluated on x
	for (unsigned int k = 0; k < AspLen; k ++)
	{
	  casacore::Matrix<casacore::Float> Asp(nX, nY);
	  Asp = 0.0;
	  casacore::Matrix<casacore::Float> dAsp(nX, nY);
	  dAsp = 0.0;
	  double amp = gsl_vector_get(x, 2*k);
      double scale = gsl_vector_get(x, 2*k+1);

      if (isnan(amp) || scale < 0.4) // GSL scale < 0
	  {
	    std::cout << "grad: nan? " << amp << " neg scale? " << scale << std::endl;
	    scale = (scale = fabs(scale)) < 0.4 ? 0.4 : scale;
	    std::cout << "reset neg scale to " << scale << std::endl;
	  }

	  const double sigma5 = 5 * scale / 2;
	  const int minI = std::max(0, (int)(center[k][0] - sigma5));
	  const int maxI = std::min(nX-1, (int)(center[k][0] + sigma5));
	  const int minJ = std::max(0, (int)(center[k][1] - sigma5));
	  const int maxJ = std::min(nY-1, (int)(center[k][1] + sigma5));

	  for (int j = minJ; j <= maxJ; j++)
	  {
	    for (int i = minI; i <= maxI; i++)
	    {
	  /*for (int j = 0; j <= nY-1; j++)
	  {
	    for (int i = 0; i <= nX-1; i++)
	    {*/
	      const int px = i;
	      const int py = j;

	      Asp(i,j) = (1.0/(sqrt(2*M_PI)*fabs(scale)))*exp(-(pow(i-center[k][0],2) + pow(j-center[k][1],2))*0.5/pow(scale,2));
	      //Asp(i,j) = (1.0/(2*M_PI*pow(scale,2)))*exp(-(pow(i-center[k][0],2) + pow(j-center[k][1],2))*0.5/pow(scale,2));
	      dAsp(i,j)= Asp(i,j) * (((pow(i-center[k][0],2) + pow(j-center[k][1],2)) / pow(scale,2) - 1) / fabs(scale)); // verified by python
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
	  double dA = 0.0;
      double dS = 0.0;
	  //std::cout << "before grad " << 2*k << ": " << gsl_vector_get(grad, 2*k) << std::endl;
	  //std::cout << "before grad " << 2*k+1 << ": " << gsl_vector_get(grad, 2*k+1) << std::endl;
	  for (int j = minJ; j <= maxJ; j++)
	  {
	    for (int i = minI; i <= maxI; i++)
	    {
	  /*for (int j = 0; j <= nY-1; j++)
	  {
	    for (int i = 0; i <= nX-1; i++)
	    {*/
	      // generate derivatives of amplitude
	      GradAmp(i,j) = (-2) * newResidual(i,j) * AspConvPsf(i,j);
	      // generate derivative of scale
	      GradScale(i,j) = (-2) * amp * newResidual(i,j) * dAspConvPsf(i,j);
	      //grad[2*k] = grad[2*k] + double(GradAmp(i,j));
	      //grad[2*k+1] = grad[2*k+1] + double(GradScale(i,j));
	      dA += double(GradAmp(i,j));
	      dS += double(GradScale(i,j));
	    }
	  }

	  gsl_vector_set(grad, 2*k, dA); //M31
	  gsl_vector_set(grad, 2*k+1, dS); //M31
	  // the following scale up doesn't seem necessary since after opt scale is back to initial guess
	  //gsl_vector_set(grad, 2*k, dA * 1e2); //G55 test scale up
	  //gsl_vector_set(grad, 2*k+1, dS * 1e6); //G55 test scale up
	  //std::cout << "after grad " << 2*k << ": " << gsl_vector_get(grad, 2*k) << " amp " << amp << std::endl;
	  //std::cout << "after grad " << 2*k+1 << ": " << gsl_vector_get(grad, 2*k+1) << " scale " << scale << std::endl;
	} // end of derivatives


  /*
  // Recover the pointers to the parameters... :-| (Phew!)
  Params<FTYPE> MyP(params);
  ResImg = MyP.ResImg();
  ResImg->size(Nx,Ny);


	tF.evalDerivatives(i,j);

	tmp = (*ResImg)(i,j);
	dA  += tmp*tF.dF(tF.DA); // element by element
	dS  += tmp*tF.dF(tF.DS);

	gsl_vector_set(grad,m++  ,dA);
	gsl_vector_set(grad,m++,dS);

      }*/
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

  bool ddel;
  casacore::Float *dptr = dirty->getStorage(ddel);
  par[0] = (void *) dptr;
  bool ddel2;
  casacore::Complex *dptr2 = psfFT->getStorage(ddel2);
  par[1] = (void *) dptr2;

  //par[0]         = (void *)dirty;
  //par[1]         = (void *)psfFT;
  std::cout << "debug dirty = " << (*dirty)(256,256) << " " << dirty << " " << dptr << std::endl;

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

    //const unsigned int len = x->size;
    //std::cout << "x len " << len << std::endl;
	//for (unsigned int i = 0; i < len; i+=2)
	//{
	    g = gsl_multimin_fdfminimizer_gradient(s);
	    std::cout << "g = " << gsl_vector_get(g, 0) << " " << gsl_vector_get(g, 1) << std::endl;

	    x = gsl_multimin_fdfminimizer_x(s);
	    std::cout << "x = " << gsl_vector_get(x, 0) << " " << gsl_vector_get(x, 1) << std::endl;
	//}

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
	//std::cout << "debug: gsl status " << status << std::endl;
	/*if (status == GSL_ENOPROG) // 27: not making progress towards solution
		gsl_multimin_fdfminimizer_restart(s);*/
	if (status)
        break;


	status = gsl_multimin_test_gradient(s->gradient, 1E-3);
	//std::cout << "debug: grad status " << status << std::endl;
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

} // end namespace

#endif