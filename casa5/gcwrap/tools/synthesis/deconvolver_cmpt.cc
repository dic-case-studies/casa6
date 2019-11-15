
/***
 * Framework independent implementation file for deconvolver...
 *
 * Implement the deconvolver component here.
 * 
 * // TODO: WRITE YOUR DESCRIPTION HERE! 
 *
 * @author
 * @version 
 ***/

#include <iostream>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Containers/Record.h>
#include <casa/Exceptions/Error.h>
#include<casa/Utilities/Assert.h>
#include <deconvolver_cmpt.h>
#include <synthesis/MeasurementEquations/Deconvolver.h>
#include <casa/Logging/LogIO.h>
#include <scimath/Mathematics/GaussianBeam.h>
//#include <casa/System/PGPlotterNull.h>
//#include <graphics/Graphics/PGPlotterLocal.h>


using namespace std;
using namespace casacore;
using namespace casa;


using namespace casacore;
namespace casac {

deconvolver::deconvolver()
{

  itsDeconv = new Deconvolver();
  itsLog = new LogIO();
    //PGPlotterInterface *worker(0);
  //try {
    //worker = new PGPlotterLocal("/NULL");
  //} catch (AipsError x) {
    //worker = new PGPlotterNull("/NULL");
  //}
  //itsPlotter = new PGPlotter(worker);
  //itsDeconv->setPGPlotter(*itsPlotter);


}

deconvolver::~deconvolver()
{

  if(itsDeconv) delete itsDeconv;
  itsDeconv=0;
  if(itsLog) delete itsLog;
  itsLog=0;

}

bool
deconvolver::open(const std::string& dirty, const std::string& psf, bool warn)
{
  if(itsDeconv==0)
    itsDeconv = new Deconvolver();
  if(itsLog==0)
    itsLog = new LogIO();
  try {
    itsDeconv->open(String(dirty), String(psf), warn);
    return true;
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return false;
}

bool
deconvolver::reopen()
{
  try {

    itsDeconv->reopen();

    return true;
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return false;
}

bool
deconvolver::close()
{

  try {
    if(itsDeconv)
      delete itsDeconv;
    itsDeconv=new Deconvolver();
    //itsDeconv->setPGPlotter(*itsPlotter);
    return true;
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return false;
}

bool
deconvolver::done()
{

  return close();
}

bool
deconvolver::summary()
{
  try {

    return itsDeconv->summary();
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }

  return false;
}

bool
deconvolver::boxmask(const std::string& mask, const std::vector<int>& blc, const std::vector<int>& trc, const ::casac::variant& fillvalue, const ::casac::variant& outsidevalue)
{
  try {
    casacore::Quantity fill(1.0, "Jy");
    casacore::Quantity out(0.0, "Jy");
    if(String(fillvalue.toString()) != String("")){
      fill=casaQuantity(fillvalue);
    }
    if(String(outsidevalue.toString()) != String("")){
      out=casaQuantity(outsidevalue);
    }
    return itsDeconv->boxmask(String(mask), Vector<Int>(blc), Vector<Int>(trc),
			      fill, out);

  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }

  return false;
}


bool
deconvolver::regionmask(const std::string& mask, const ::casac::record& region,
			const ::casac::variant& boxes, const double value)
{
  Bool rstat(false);
  try {
    
    casacore::Record *pRegion = toRecord(region);
    if(pRegion->nfields() < 2){
      delete pRegion;
      pRegion=0;
    }
    casacore::Matrix<casacore::Quantity> quantBoxes;
    Matrix<Float> theBoxes;
    if((boxes.type() == ::casac::variant::DOUBLEVEC) ||
       (boxes.type() == ::casac::variant::INTVEC) ){
      Vector<Double> boxVector;
      if(boxes.type() == ::casac::variant::DOUBLEVEC)
	boxVector=Vector<Double>(boxes.toDoubleVec());
      else{
	std::vector<int> intvector=boxes.toIntVec();
	boxVector.resize(intvector.size());
	convertArray(boxVector, Vector<Int>(intvector));
	
      }
      Vector<Int> shape = boxes.arrayshape();
      theBoxes.resize(IPosition(shape));
      Vector<Double> localbox(boxVector);
      casacore::convertArray(theBoxes,localbox.reform(IPosition(shape)));
      quantBoxes.resize(theBoxes.shape());
      for (Int k=0; k < shape[1]; ++k){
	for(Int j=0; j < shape[0]; ++j){
	  quantBoxes(j,k)=casacore::Quantity(theBoxes(j,k), "pix");
	} 	    
      }
    }
    else if(boxes.type() == ::casac::variant::STRINGVEC){
      Vector<Int> shape = boxes.arrayshape();
      casacore::Vector<casacore::Quantity> theCorners;
      toCasaVectorQuantity(boxes, theCorners);
      if(casacore::product(shape) != Int(theCorners.nelements())){
	throw(AipsError("Error in converting list of strings to Quantities"));
      }
      if(shape.nelements()==1){
	// list is a single file so assuming every 4 elements is blc trc
	Int nrows=shape[0]/4;
	if(shape[0] != nrows*4)
	  throw(AipsError("list of blc trc has to have 4 elements each"));
	quantBoxes.resize(nrows,shape(0));
      }
      else if(shape.nelements()==2){
	quantBoxes.resize(shape);
      }
      else{
	throw(AipsError("Only dealing with 2D blc and trc's for now"));
      }
      Vector<casacore::Quantity> refQuantBoxes(quantBoxes.reform(IPosition(1, theCorners.nelements())));
      refQuantBoxes=theCorners;
      
      
    }
    
    
    rstat=itsDeconv->regionmask(String(mask), pRegion, quantBoxes, Float(value));
    if(pRegion != 0)
      delete pRegion;
    
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }



  return rstat;
}



bool
deconvolver::clipimage(const std::string& clippedimage, const std::string& inputimage, const ::casac::variant& threshold)
{

  try {
    casacore::Quantity thresh(0.0, "Jy");
    if(String(threshold.toString()) != String("")){
      thresh=casaQuantity(threshold);
    }

    return itsDeconv->clipimage(String(clippedimage), String(inputimage), 
				thresh);
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return false;
}

 casac::record*  deconvolver::fullclarkclean(const int niter, const double gain, const ::casac::variant& threshold, const std::string& model, const std::string& mask, const double cyclefactor){
   casac::record* rstat(0);
    try {
      casacore::Quantity thresh(0.0, "Jy");
      if(String(threshold.toString()) != String("")){
	thresh=casaQuantity(threshold);
      }
      
      Float maxResidual;
      Int iterations;
      Bool converged= itsDeconv->clarkclean(niter, gain, thresh, String(model), String(mask), maxResidual, iterations, cyclefactor);
      Record theRec;
      theRec.define("converged", converged);
      theRec.define("maxresidual", maxResidual);
      theRec.define("iterations", iterations);
      rstat = fromRecord(theRec);

    }
    catch  (AipsError x) {
      *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
      RETHROW(x);
    }
    
    return rstat;
  }


bool
deconvolver::clarkclean(const int niter, const double gain, const ::casac::variant& threshold, const bool displayprogress, const std::string& model, const std::string& mask, const int histbins, const std::vector<int>& psfpatchsize, const double maxextpsf, const double speedup, const int maxnumpix, const int maxnummajcycles, const int maxnummineriter)
{
  try {
    casacore::Quantity thresh(0.0, "Jy");
    if(String(threshold.toString()) != String("")){
      thresh=casaQuantity(threshold);
    }
    

    return itsDeconv->clarkclean(niter, gain, thresh, displayprogress, 
				 String(model), String(mask), histbins, 
				 Vector<Int>(psfpatchsize), maxextpsf, 
				 speedup, maxnumpix, maxnummajcycles, 
				 maxnummineriter);
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }

  return false;
}

std::string
deconvolver::dirtyname()
{

  try {


    return itsDeconv->dirtyname();
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }


  return "";
}

std::string
deconvolver::psfname()
{

  try {


    return itsDeconv->psfname();
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return "";
}

bool
deconvolver::make(const std::string& image)
{
  try {


    return itsDeconv->make(String(image));
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }

  return true;
}

bool
deconvolver::convolve(const std::string& convolvemodel, const std::string& model)
{
  try {

    return itsDeconv->convolve(String(convolvemodel), String(model));
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
 
  return true;
}

bool
deconvolver::makegaussian(const std::string& gaussianimage, const ::casac::variant& bmaj, const ::casac::variant& bmin, const ::casac::variant& bpa, const bool normalize)
{
  try {

    casacore::Quantity maj(1, "arcsec");
    casacore::Quantity min(1, "arcsec");
    casacore::Quantity pa(0, "deg");
    if(String(bmaj.toString()) != String("")){
      maj=casaQuantity(bmaj);
    }
    if(String(bmin.toString()) != String("")){
      min=casaQuantity(bmin);
    }
    if(String(bpa.toString()) != String("")){
      pa=casaQuantity(bpa);
    }
    GaussianBeam beam(maj, min, pa);
    return itsDeconv->makegaussian(String(gaussianimage), beam,
				   normalize);
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }

  return true;
}



bool
deconvolver::state()
{
  cout << "Not implemented " << endl;
    // TODO : IMPLEMENT ME HERE !
  return true;
}

bool
deconvolver::updatestate(const std::string& /*f*/, const std::string& /*method*/)
{

  cout << "Not implemented" << endl;
    // TODO : IMPLEMENT ME HERE !
  return true;
}


casac::record* 
deconvolver::clean(const std::string& algorithm, const int niter, const double gain, const ::casac::variant& threshold, const bool displayprogress, const std::string& model, const std::string& mask)
{
  casac::record* rstat(0);
  try {
    
    casacore::Quantity thresh(0.0, "Jy");
    if(String(threshold.toString()) != String("")){
      thresh=casaQuantity(threshold);
    }
    Float maxResidual=0.0;
    Int iterations=0;
    Bool converged=itsDeconv->clean(String(algorithm), niter, gain, thresh, 
				    displayprogress, String(model), String(mask), maxResidual, iterations);
    Record theRec;
    theRec.define("converged", converged);
    theRec.define("maxresidual", maxResidual);
    theRec.define("iterations", iterations);
    rstat = fromRecord(theRec);  
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  
  return rstat;
}


casac::record* 
deconvolver::naclean(const int niter, const double gain, const ::casac::variant& threshold, const std::string& model, const std::string& mask, const int masksupp, const std::string& memory, const double numsigma)
{
  casac::record* rstat(0);
  try {
    
    casacore::Quantity thresh(0.0, "Jy");
    if(String(threshold.toString()) != String("")){
      thresh=casaQuantity(threshold);
    }
    Float maxResidual=0.0;
    Int iterations=0;
    Int memType=2;
    if(casacore::String(memory)==casacore::String("none"))
      memType=0;
    else if (casacore::String(memory)==casacore::String("weak"))
      memType=1;
    else if(casacore::String(memory)==casacore::String("strong"))
      memType=3;
    Bool converged=itsDeconv->naclean(niter, gain, thresh, 
				      String(model), String(mask), masksupp, memType, Float(numsigma), maxResidual, iterations);
    Record theRec;
    theRec.define("converged", converged);
    theRec.define("maxresidual", maxResidual);
    theRec.define("iterations", iterations);
    rstat = fromRecord(theRec);  
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  
  return rstat;
}


bool
deconvolver::setscales(const std::string& scalemethod, const int nscales, const std::vector<double>& uservector)
{
  try {
    Vector<Float> vec(uservector.size());
    if(vec.nelements() != 0)
      convertArray(vec, Vector<Double>(uservector));
    return itsDeconv->setscales(String(scalemethod), nscales, vec);
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return false;
}

bool
deconvolver::ft(const std::string& model, const std::string& transform)
{
  try {
    

    return itsDeconv->ft(String(model), String(transform));
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return false;
}

bool
deconvolver::restore(const std::string& model, const std::string& image, const ::casac::variant& bmaj, const ::casac::variant& bmin, const ::casac::variant& bpa)
{
  try {

    GaussianBeam beam;
    if(toVectorString(bmaj.toStringVec()).nelements() != 0){
      casacore::Quantity bMaj(1, "arcsec");
      casacore::Quantity bMin(1, "arcsec");
      casacore::Quantity pa(0, "deg");
      bMaj=casaQuantity(bmaj);
      if(toVectorString(bmin.toStringVec()).nelements() != 0){
	bMin=casaQuantity(bmin);
      }
      if(toVectorString(bpa.toStringVec()).nelements() != 0){
	pa=casaQuantity(bpa);
      }
      beam=GaussianBeam(bMaj, bMin, pa); 
    } 
    return itsDeconv->restore(String(model), String(image), beam);
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }

  return false;
}

bool
deconvolver::residual(const std::string& model, const std::string& image)
{

  try {
    return itsDeconv->residual(String(model), String(image));
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }

  return false;
}

bool
deconvolver::smooth(const std::string& model, const std::string& image, const ::casac::variant& bmaj, const ::casac::variant& bmin, const ::casac::variant& bpa, const bool normalize)
{

  try {
    casacore::Quantity maj(1, "arcsec");
    casacore::Quantity min(1, "arcsec");
    casacore::Quantity pa(0, "deg");
    if(String(bmaj.toString()) != String("")){
      maj=casaQuantity(bmaj);
    }
    if(String(bmin.toString()) != String("")){
      min=casaQuantity(bmin);
    }
    if(String(bpa.toString()) != String("")){
      pa=casaQuantity(bpa);
    }
    

    GaussianBeam beam(maj, min, pa);
    return itsDeconv->smooth(String(model), String(image), beam,
			     normalize);
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }

  return false;
}

bool
deconvolver::mem(const std::string& entropy, const int niter, const ::casac::variant& sigma, const ::casac::variant& targetflux, const bool constrainflux, const bool displayprogress, const std::string& model, const std::string& prior, const std::string& mask, const bool imageplane)
{
  try {
    casacore::Quantity targflux(1.0, "Jy");
    if(String(targetflux.toString()) != String("")){
      targflux=casaQuantity(targetflux);
    }
    casacore::Quantity sig(1.0, "mJy");
    if(String(sigma.toString()) != String("")){
      sig=casaQuantity(sigma);
    }
    

    return itsDeconv->mem(String(entropy), niter, sig, targflux, constrainflux,
			  displayprogress, String(model), String(prior), 
			  String(mask), imageplane);
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return false;
}

bool
deconvolver::makeprior(const std::string& prior, const std::string& templateimage, const ::casac::variant& lowclipfrom, const ::casac::variant& lowclipto, const ::casac::variant& highclipfrom, const ::casac::variant& highclipto, const std::vector<int>& blc, const std::vector<int>& trc)
{

  try {
    casacore::Quantity lct(0.0, "Jy");
    if(String(lowclipto.toString()) != String("")){
      lct=casaQuantity(lowclipto);
    }
    casacore::Quantity lcf(0.0, "Jy");
    if(String(lowclipfrom.toString()) != String("")){
      lcf=casaQuantity(lowclipfrom);
    }
    casacore::Quantity hct(9.0e20, "Jy");
    if(String(highclipto.toString()) != String("")){
      hct=casaQuantity(highclipto);
    }
    casacore::Quantity hcf(9.0e20, "Jy");
    if(String(highclipfrom.toString()) != String("")){
      hcf=casaQuantity(highclipfrom);
    }
    
    return itsDeconv->makeprior(String(prior), String(templateimage), lcf, lct,
				hcf, hct, Vector<Int>(blc), Vector<Int>(trc));
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }    
  return false;
}


bool deconvolver::mtopen(const int ntaylor,
			 const std::vector<double>& scalevector,
			 const std::vector<std::string>& psfs)
{

  if(itsDeconv==0)
    itsDeconv = new Deconvolver();
  if(itsLog==0)
    itsLog = new LogIO();

  try {

    Vector<Float> vec(scalevector.size());
    if(vec.nelements() != 0)
      convertArray(vec, Vector<Double>(scalevector));

    Vector<String> psfvec(toVectorString(psfs));


    return itsDeconv->mtopen(ntaylor, vec, psfvec);

  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return false;
}



casac::record* 
deconvolver::mtclean(const std::vector<std::string>& residuals,
			  const std::vector<std::string>& models,
			  const int niter,
			  const double gain, 
			  const ::casac::variant& threshold, 
			  const bool displayprogress,
			  const std::string& mask)
{

  casac::record* rstat(0);

  try {

    Vector<String> resvec(toVectorString(residuals));
    Vector<String> modvec(toVectorString(models));

    casacore::Quantity thresh(0.0, "Jy");
    if(String(threshold.toString()) != String(""))
      thresh=casaQuantity(threshold);

    Float maxResidual=0.0;
    Int iterations=0;
    Bool converged=itsDeconv->mtclean(resvec, modvec, (Int)niter, (Float)gain, thresh, 
				      (Bool)displayprogress, String(mask), maxResidual, iterations );

    Record theRec;
    theRec.define("converged", converged);
    theRec.define("maxresidual", maxResidual);
    theRec.define("iterations", iterations);
    rstat = fromRecord(theRec);  

  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }

  return rstat;
}

bool deconvolver::mtrestore(const std::vector<std::string>& models,
			  const std::vector<std::string>& residuals,
			  const std::vector<std::string>& images,
			  const ::casac::variant& bmaj, 
			  const ::casac::variant& bmin, 
			  const ::casac::variant& bpa)
{

 try {

    Vector<String> modvec(toVectorString(models));
    Vector<String> resvec(toVectorString(residuals));
    Vector<String> imvec(toVectorString(images));

    casacore::Quantity maj(1, "arcsec");
    casacore::Quantity min(1, "arcsec");
    casacore::Quantity pa(0, "deg");
    if(String(bmaj.toString()) != String("")){
      maj=casaQuantity(bmaj);
    }
    if(String(bmin.toString()) != String("")){
      min=casaQuantity(bmin);
    }
    if(String(bpa.toString()) != String("")){
      pa=casaQuantity(bpa);
    }
    GaussianBeam beam(maj, min, pa);

    return itsDeconv->mtrestore(modvec, resvec, imvec, beam);

  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }

  return false;
}



bool deconvolver::mtcalcpowerlaw(const std::vector<std::string>& images,
                                 const std::vector<std::string>& residuals,
                                 const std::string& alphaname,
                                 const std::string& betaname,
                                 const ::casac::variant& threshold, 
                                 const bool calcerror)
{

 try {


    Vector<String> imvec(toVectorString(images));
    Vector<String> resvec(toVectorString(residuals));

    casacore::Quantity thresh(0.0, "Jy");
    if(String(threshold.toString()) != String(""))
      thresh=casaQuantity(threshold);
   
   return itsDeconv->mtcalcpowerlaw(imvec, resvec, alphaname, betaname, thresh, calcerror);

  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }


  return false;
}


} // casac namespace

