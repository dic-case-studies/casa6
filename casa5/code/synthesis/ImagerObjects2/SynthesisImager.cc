//# SynthesisImager.cc: Implementation of Imager.h
//# Copyright (C) 1997-2008
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This program is free software; you can redistribute it and/or modify it
//# under the terms of the GNU General Public License as published by the Free
//# Software Foundation; either version 2 of the License, or (at your option)
//# any later version.
//#
//# This program is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
//# more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with this program; if not, write to the Free Software Foundation, Inc.,
//# 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id$

#include <casa/Exceptions/Error.h>
#include <casa/iostream.h>
#include <casa/sstream.h>

#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/ArrayLogical.h>

#include <casa/Logging.h>
#include <casa/Logging/LogIO.h>
#include <casa/Logging/LogMessage.h>
#include <casa/Logging/LogSink.h>
#include <casa/Logging/LogMessage.h>
#include <casa/System/ProgressMeter.h>

#include <casa/OS/DirectoryIterator.h>
#include <casa/OS/File.h>
#include <casa/OS/HostInfo.h>
#include <casa/OS/Path.h>

#include <lattices/LRegions/LCBox.h>

#include <measures/Measures/MeasTable.h>

#include <ms/MeasurementSets/MSHistoryHandler.h>
#include <ms/MeasurementSets/MeasurementSet.h>
#include <ms/MSSel/MSSelection.h>


#include <synthesis/ImagerObjects/SIIterBot.h>


#include <synthesis/ImagerObjects/SynthesisUtilMethods.h>
#include <synthesis/ImagerObjects/SIImageStore.h>
#include <synthesis/ImagerObjects/SIImageStoreMultiTerm.h>

#include <synthesis/ImagerObjects2/SynthesisImager.h>

#include <synthesis/MeasurementEquations/ImagerMultiMS.h>
#include <msvis/MSVis/VisSetUtil.h>
#include <msvis/MSVis/VisImagingWeight.h>

#include <synthesis/TransformMachines2/GridFT.h>
#include <synthesis/TransformMachines2/WPConvFunc.h>
#include <synthesis/TransformMachines2/WProjectFT.h>
#include <synthesis/TransformMachines2/VisModelData.h>
//#include <synthesis/TransformMachines2/AWProjectFT.h>
#include <synthesis/TransformMachines2/MultiTermFTNew.h>
#include <synthesis/TransformMachines2/AWProjectWBFTNew.h>
#include <synthesis/TransformMachines2/AWConvFunc.h>
#include <synthesis/TransformMachines2/AWConvFuncEPJones.h>
#include <synthesis/TransformMachines2/NoOpATerm.h>
#include <synthesis/TransformMachines2/EVLAAperture.h>

#include <casadbus/viewer/ViewerProxy.h>
#include <casadbus/plotserver/PlotServerProxy.h>
//#include <casadbus/utilities/BusAccess.h>
//#include <casadbus/session/DBusSession.h>

#include <sys/types.h>
#include <unistd.h>
using namespace std;

using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN
namespace refim {//# namespace refactor imaging
  
using namespace casacore;
using namespace casa;
using namespace casacore;
using namespace casa::refim;
using namespace casacore;
using namespace casa::vi;
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  SynthesisImager::SynthesisImager() : itsMappers(SIMapperCollection()), writeAccess_p(true)
  {

     imwgt_p=VisImagingWeight("natural");
     imageDefined_p=false;
     useScratch_p=false;
     readOnly_p=true;

     mss_p=0;
     vi_p=0;

     facetsStore_p=-1;
     unFacettedImStore_p=NULL;

  }
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  SynthesisImager::~SynthesisImager() 
  {
    LogIO os( LogOrigin("SynthesisImager","destructor",WHERE) );
    os << LogIO::DEBUG1 << "SynthesisImager destroyed" << LogIO::POST;

    for (uInt k=0; k < mss_p.nelements(); ++k){
      delete mss_p[k];
    }
    //    cerr << "IN DESTR"<< endl;
    //    VisModelData::listModel(mss4vi_p[0]);
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  Bool SynthesisImager::selectData(const String& msname, 
				   const String& spw, 
				   const String& freqBeg, 
				   const String& freqEnd, 
				   const MFrequency::Types freqframe, 
				   const String& field, 
				   const String& antenna, 
				   const String& timestr,
				   const String& scan, 
				   const String& obs, 
				   const String& state,
				   const String& uvdist, 
				   const String& taql,
				   const Bool usescratch, 
				   const Bool readonly, 
				   const Bool incrModel)
  {
    SynthesisParamsSelect pars;
    pars.msname=msname;
    pars.spw=spw;
    pars.freqbeg=freqBeg;
    pars.freqend=freqEnd;
    pars.freqframe=freqframe;
    pars.field=field;
    pars.antenna=antenna;
    pars.timestr=timestr;
    pars.scan=scan;
    pars.obs=obs;
    pars.state=state;
    pars.uvdist=uvdist;
    pars.taql=taql;
    pars.usescratch=usescratch;
    pars.readonly=readonly;
    pars.incrmodel=incrModel;

    String err = pars.verify();

    if( err.length()>0 ) throw(AipsError("Invalid Selection parameters : " + err));

    selectData( pars );

    return true;
  }
  
  Bool SynthesisImager::selectData(const SynthesisParamsSelect& selpars)
  {
    LogIO os( LogOrigin("SynthesisImager","selectData",WHERE) );

    try
      {

    //Respect the readonly flag...necessary for multi-process access
    MeasurementSet thisms(selpars.msname, TableLock(TableLock::AutoNoReadLocking),
				selpars.readonly ? Table::Old : Table::Update);
    thisms.setMemoryResidentSubtables (MrsEligibility::defaultEligible());
    useScratch_p=selpars.usescratch;
    readOnly_p = selpars.readonly;
    //    cout << "**************** usescr : " << useScratch_p << "     readonly : " << readOnly_p << endl;
    //if you want to use scratch col...make sure they are there
    if(selpars.usescratch && !selpars.readonly){
      VisSetUtil::addScrCols(thisms, true, false, true, false);
      VisModelData::clearModel(thisms);
    }
    if(!selpars.incrmodel && !selpars.usescratch && !selpars.readonly)
    	VisModelData::clearModel(thisms, selpars.field, selpars.spw);

    os << "MS : " << selpars.msname << " | ";

    //Some MSSelection 
    //If everything is empty (which is valid) it will throw an exception..below
    //So make sure the main defaults are not empy i.e field and spw
    MSSelection thisSelection;
    if(selpars.field != ""){
      thisSelection.setFieldExpr(selpars.field);
      os << "Selecting on fields : " << selpars.field << " | " ;//LogIO::POST;
    }else
      thisSelection.setFieldExpr("*");
    if(selpars.spw != ""){
	thisSelection.setSpwExpr(selpars.spw);
	os << "Selecting on spw :"<< selpars.spw  << " | " ;//LogIO::POST;
    }else
      thisSelection.setSpwExpr("*");
    
    if(selpars.antenna != ""){
      Vector<String> antNames(1, selpars.antenna);
      // thisSelection.setAntennaExpr(MSSelection::nameExprStr( antNames));
      thisSelection.setAntennaExpr(selpars.antenna);
      os << "Selecting on antenna names : " << selpars.antenna << " | " ;//LogIO::POST;
	
    }            
    if(selpars.timestr != ""){
	thisSelection.setTimeExpr(selpars.timestr);
	os << "Selecting on time range : " << selpars.timestr << " | " ;//LogIO::POST;	
      }
    if(selpars.uvdist != ""){
      thisSelection.setUvDistExpr(selpars.uvdist);
      os << "Selecting on uvdist : " << selpars.uvdist << " | " ;//LogIO::POST;	
    }
    if(selpars.scan != ""){
      thisSelection.setScanExpr(selpars.scan);
      os << "Selecting on scan : " << selpars.scan << " | " ;//LogIO::POST;	
    }
    if(selpars.obs != ""){
      thisSelection.setObservationExpr(selpars.obs);
      os << "Selecting on Observation Expr : " << selpars.obs << " | " ;//LogIO::POST;	
    }
    if(selpars.state != ""){
      thisSelection.setStateExpr(selpars.state);
      os << "Selecting on Scan Intent/State : " << selpars.state << " | " ;//LogIO::POST;	
    }
    if(selpars.taql != ""){
	thisSelection.setTaQLExpr(selpars.taql);
	os << "Selecting via TaQL : " << selpars.taql << " | " ;//LogIO::POST;	
    }
    os << "[Opened " << (readOnly_p?"in readonly mode":(useScratch_p?"with scratch model column":"with virtual model column"))  << "]" << LogIO::POST;
    TableExprNode exprNode=thisSelection.toTableExprNode(&thisms);
    if(!(exprNode.isNull()))
      {
	mss_p.resize(mss_p.nelements()+1, false, true);
	mss_p[mss_p.nelements()-1]=new const  MeasurementSet(thisms(exprNode));
	os << "  NRows selected : " << (mss_p[mss_p.nelements()-1])->nrow() << LogIO::POST;
      }
    else{
      throw(AipsError("Selection for given MS "+selpars.msname+" is invalid"));
    }
    //We should do the select channel here for  the VI construction later
    //Need a cross check between channel selection and ms
    // replace below if/when viFrquencySelectionUsingChannels takes in a MSSelection
    // rather than the following gymnastics
    {
      Matrix<Int> chanlist = thisSelection.getChanList(mss_p[mss_p.nelements()-1]);
      
      IPosition shape = chanlist.shape();
      uInt nSelections = shape[0];
      Int spw,chanStart,chanEnd,chanStep,nchan;


      if(selpars.freqbeg==""){
    	  vi::FrequencySelectionUsingChannels channelSelector;
    	  //////////This is not implemented
    	  //channelSelector.add(thisSelection);
    	  fselections_p.add(channelSelector);

      }
      else{
    	  Quantity freq;
    	  Quantity::read(freq, selpars.freqbeg);
    	  Double lowfreq=freq.getValue("Hz");
    	  Quantity::read(freq, selpars.freqend);
    	  Double topfreq=freq.getValue("Hz");
    	  vi::FrequencySelectionUsingFrame channelSelector(selpars.freqframe);
    	  for(uInt k=0; k < nSelections; ++k)
    		  channelSelector.add(chanlist(k,0), lowfreq, topfreq);
    	  fselections_p.add(channelSelector);
      }


    }
    writeAccess_p=writeAccess_p && !selpars.readonly;
    createVisSet(writeAccess_p);

    /////// Remove this when the new vi/vb is able to get the full freq range.

    //// Set the data column on which to operate
    //    cout << "Using col : " << selpars.datacolumn << endl;
    if( selpars.datacolumn.contains("data") || selpars.datacolumn.contains("obs") ) 
      {datacol_p = refim::FTMachine::OBSERVED; }
    else if( selpars.datacolumn.contains("corr") )
      {datacol_p = refim::FTMachine::CORRECTED; }
    else { os << LogIO::WARN << "Invalid data column : " << datacol_p << ". Using corrected (or observed if corrected doesn't exist)" << LogIO::POST;  datacol_p = refim::FTMachine::CORRECTED; }

      }
    catch(AipsError &x)
      {
	throw( AipsError("Error in selectData() : "+x.getMesg()) );
      }

    return true;

  }

 

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  Bool SynthesisImager::defineImage(const String& imagename, const Int nx, const Int ny,
				    const Quantity& cellx, const Quantity& celly,
				    const String& stokes,
				    const MDirection& phaseCenter, 
				    const Int nchan,
				    const Quantity&freqStart,
				    const Quantity& freqStep, 
				    const Vector<Quantity>& restFreq,
				    const Int facets,
				    const String ftmachine, 
				    const Int nTaylorTerms,
				    const Quantity& refFreq,
				    const Projection& projection,
				    const Quantity& distance,
				    const MFrequency::Types& freqFrame,
				    const Bool trackSource, 
				    const MDirection& trackDir, 
				    const Bool overwrite,
				    const Float padding, 
				    const Bool useAutocorr, 
				    const Bool useDoublePrec, 
				    const Int wprojplanes, 
				    const String convFunc, 
				    const String startmodel,
				    // The extra params for WB-AWP
				    const Bool aTermOn,//    = true,
				    const Bool psTermOn,//   = true,
				    const Bool mTermOn,//    = false,
				    const Bool wbAWP,//      = true,
				    const String cfCache,//  = "",
				    const Bool doPointing,// = false,
				    const Bool doPBCorr,//   = true,
				    const Bool conjBeams,//  = true,
				    const Float computePAStep,         //=360.0
				    const Float rotatePAStep          //=5.0
				    )
{
  String err("");

  SynthesisParamsImage impars;
  impars.imageName=imagename;
  Vector<Int> ims(2);ims[0]=nx; ims[1]=ny;
  impars.imsize=ims;
  Vector<Quantity> cells(2); cells[0]=cellx, cells[1]=celly;
  impars.cellsize=cells;
  impars.stokes=stokes;
  impars.phaseCenter=phaseCenter;
  impars.nchan=nchan;
  impars.freqStart=freqStart;
  impars.freqStep=freqStep;
  impars.restFreq=restFreq;
  impars.nTaylorTerms=nTaylorTerms;
  impars.refFreq=refFreq;
  impars.projection=projection;
  impars.freqFrame=freqFrame;
  impars.overwrite=overwrite;
  impars.startModel=startmodel;

  err += impars.verify();

  SynthesisParamsGrid gridpars;
  gridpars.ftmachine=ftmachine;
  gridpars.distance=distance;
  gridpars.trackSource=trackSource;
  gridpars.trackDir=trackDir;
  gridpars.padding=padding;
  gridpars.facets=facets;
  gridpars.useAutoCorr=useAutocorr;
  gridpars.useDoublePrec=useDoublePrec;
  gridpars.wprojplanes=wprojplanes;
  gridpars.convFunc=convFunc;
  gridpars.aTermOn=aTermOn;
  gridpars.psTermOn=psTermOn;
  gridpars.mTermOn=mTermOn;
  gridpars.wbAWP=wbAWP;
  gridpars.cfCache=cfCache;
  gridpars.doPointing=doPointing;
  gridpars.doPBCorr=doPBCorr;
  gridpars.conjBeams=conjBeams;
  gridpars.computePAStep=computePAStep;
  gridpars.rotatePAStep=rotatePAStep;

  err += gridpars.verify();

  if( err.length()>0 ) throw(AipsError("Invalid Image/Gridding parameters : " + err));

  defineImage( impars, gridpars );

  return true;
}

  Bool SynthesisImager::defineImage(SynthesisParamsImage& impars, 
			   const SynthesisParamsGrid& gridpars)
  {

    LogIO os( LogOrigin("SynthesisImager","defineImage",WHERE) );
    if(mss_p.nelements() ==0)
      os << "SelectData has to be run before defineImage" << LogIO::EXCEPTION;

    CoordinateSystem csys;
    CountedPtr<refim::FTMachine> ftm, iftm;

    try
      {

	os << "Define image [" << impars.imageName << "]" << LogIO::POST;

	csys = impars.buildCoordinateSystem( vi_p );
	IPosition imshape = impars.shp();

	if( (itsMappers.nMappers()==0) || 
	    (impars.imsize[0]*impars.imsize[1] > itsMaxShape[0]*itsMaxShape[1]))
	  {
	    itsMaxShape=imshape;
	    itsMaxCoordSys=csys;
	  }

	/*
	os << "Define image  [" << impars.imageName << "] : nchan : " << impars.nchan 
	   //<< ", freqstart:" << impars.freqStart.getValue() << impars.freqStart.getUnit() 
	   << ", start:" << impars.start
	   <<  ", imsize:" << impars.imsize 
	   << ", cellsize: [" << impars.cellsize[0].getValue() << impars.cellsize[0].getUnit() 
	   << " , " << impars.cellsize[1].getValue() << impars.cellsize[1].getUnit() 
	   << LogIO::POST;
	*/
      }
    catch(AipsError &x)
      {
	os << "Error in building Coordinate System and Image Shape : " << x.getMesg() << LogIO::EXCEPTION;
      }

	
    try
      {
	createFTMachine(ftm, iftm, gridpars.ftmachine, impars.nTaylorTerms, gridpars.mType, 
			gridpars.facets, gridpars.wprojplanes,
			gridpars.padding,gridpars.useAutoCorr,gridpars.useDoublePrec,
			gridpars.convFunc,
			gridpars.aTermOn,gridpars.psTermOn, gridpars.mTermOn,
			gridpars.wbAWP,gridpars.cfCache,gridpars.doPointing,
			gridpars.doPBCorr,gridpars.conjBeams,
			gridpars.computePAStep,gridpars.rotatePAStep,
			gridpars.interpolation, impars.freqFrameValid);

      }
    catch(AipsError &x)
      {
	os << "Error in setting up FTMachine() : " << x.getMesg() << LogIO::EXCEPTION;
      }

    try
      {
	appendToMapperList(impars.imageName,  csys,  impars.shp(),
			   ftm, iftm,
			   gridpars.distance, gridpars.facets, impars.overwrite,
			   gridpars.mType, impars.nTaylorTerms);
	imageDefined_p=true;
      }
    catch(AipsError &x)
      {
	os << "Error in adding Mapper : "+x.getMesg() << LogIO::EXCEPTION;
      }

    // Set the model image for prediction -- Call an SIImageStore function that does the REGRIDDING.
    /*
      if( startmodel.length()>0 && !itsCurrentImages.null() )
      {
      os << "Setting " << startmodel << " as starting model for prediction " << LogIO::POST;
      itsCurrentImages->setModelImage( startmodel );
      }
    */
    return true;
  }
 
  Bool SynthesisImager::defineImage(CountedPtr<SIImageStore> imstor, 
				    const String& ftmachine)
  {
    CountedPtr<refim::FTMachine> ftm, iftm;

    // The following call to createFTMachine() uses the
    // following defaults
    //
    // facets=1, wprojplane=1, padding=1.0, useAutocorr=false, 
    // useDoublePrec=true, gridFunction=String("SF")
    //
    createFTMachine(ftm, iftm, ftmachine);
    
    Int id=itsMappers.nMappers();
    CoordinateSystem csys =imstor->residual()->coordinates();
    IPosition imshape=imstor->residual()->shape();
    Int nx=imshape[0], ny=imshape[1];
    if( (id==0) || (nx*ny > itsMaxShape[0]*itsMaxShape[1]))
      {
	itsMaxShape=imshape;
	itsMaxCoordSys=csys;
      }

    itsMappers.addMapper(  createSIMapper( "default", imstor, ftm, iftm, id ) );
    
    return true;
  }

   ///////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////
  void SynthesisImager::setComponentList(const ComponentList& cl, Bool sdgrid){
	  String cft="SimpleComponentFTMachine";
	  if(sdgrid)
		  cft="SimpCompGridFTMachine";
	  CountedPtr<SIMapper> sm=new SIMapper(cl, cft);
	  itsMappers.addMapper(sm);
	  ////itsMappers.addMapper(  createSIMapper( mappertype, imstor, ftm, iftm, id) );

  }
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //////////////////////Reset the Mapper
  ////////////////////
  void SynthesisImager::resetMappers(){
    ////reset code
	itsMappers=SIMapperCollection();
	unFacettedImStore_p=NULL;
  }
//////////////////////////////////////////////////////////////////
/////////////////////////////////////////////
  CountedPtr<SIImageStore> SynthesisImager::imageStore(const Int id)
  {
    if(facetsStore_p >1)
      {
	if(id==0)
	  {
	    return unFacettedImStore_p;
	  }
	else
	  {
	    return itsMappers.imageStore(facetsStore_p*facetsStore_p+id-1);
	  }
      }
    return itsMappers.imageStore(id);
  }


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  void SynthesisImager::executeMajorCycle(const Record& controlRecord)
  {
    LogIO os( LogOrigin("SynthesisImager","executeMajorCycle",WHERE) );


    Bool lastcycle=false;
    if( controlRecord.isDefined("lastcycle") )
      {
	controlRecord.get( "lastcycle" , lastcycle );
	//cout << "lastcycle : " << lastcycle << endl;
      }
    //else {cout << "No lastcycle" << endl;}

    os << "----------------------------------------------------------- Run ";
    if (lastcycle) os << "(Last) " ;
    os << "Major Cycle -------------------------------------" << LogIO::POST;
    
    try
      {    
	runMajorCycle(false, lastcycle);

	itsMappers.releaseImageLocks();

      }
    catch(AipsError &x)
      {
	throw( AipsError("Error in running Major Cycle : "+x.getMesg()) );
      }    

  }// end of executeMajorCycle
  //////////////////////////////////////////////
  /////////////////////////////////////////////

  void SynthesisImager::makePSF()
    {
      LogIO os( LogOrigin("SynthesisImager","makePSF",WHERE) );

      os << "----------------------------------------------------------- Make PSF ---------------------------------------------" << LogIO::POST;
    
      try
      {
	runMajorCycle(true, false);

    	  if(facetsStore_p >1)
	    {
	      setPsfFromOneFacet();
	    }

    	  itsMappers.releaseImageLocks();

      }
      catch(AipsError &x)
      {
    	  throw( AipsError("Error in making PSF : "+x.getMesg()) );
      }

    }


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void SynthesisImager::predictModel(){
    LogIO os( LogOrigin("SynthesisImager","predictModel ",WHERE) );
    
    {
      vi_p->originChunks();
      vi_p->origin();
      vi::VisBuffer2* vb=vi_p->getVisBuffer();
      itsMappers.initializeDegrid(*vb);
      for (vi_p->originChunks(); vi_p->moreChunks();vi_p->nextChunk())
	{
	  
	  for (vi_p->origin(); vi_p->more();vi_p->next())
	    {
	      
	      vb->setVisCubeModel(Cube<Complex>(vb->visCubeModel().shape(), Complex(0.0, 0.0)));
	      itsMappers.degrid(*vb, !useScratch_p);
	      if(vi_p->isWritable() && useScratch_p)
		vi_p->writeVisModel(vb->visCubeModel());
	    }
	}
      itsMappers.finalizeDegrid(*vb);
    }
    
  }// end of predictModel

   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  Bool SynthesisImager::weight(const String& type, const String& rmode,
			       const Quantity& noise, const Double robust,
			       const Quantity& fieldofview,
			       const Int npixels, const Bool multiField,
			       const String& filtertype, const Quantity& filterbmaj,
			       const Quantity& filterbmin, const Quantity& filterbpa   )
  {
    LogIO os(LogOrigin("SynthesisImager", "weight()", WHERE));

       try {
    	//Int nx=itsMaxShape[0];
    	//Int ny=itsMaxShape[1];
	 Quantity cellx=Quantity(itsMaxCoordSys.increment()[0], itsMaxCoordSys.worldAxisUnits()[0]);
	 Quantity celly=Quantity(itsMaxCoordSys.increment()[1], itsMaxCoordSys.worldAxisUnits()[1]);
	 os << LogIO::NORMAL // Loglevel INFO
	    << "Set imaging weights : " ; //<< LogIO::POST;
	 
	 if (type=="natural") {
	   os << LogIO::NORMAL // Loglevel INFO
	      << "Natural weighting" << LogIO::POST;
	   imwgt_p=VisImagingWeight("natural");
	 }
      else if (type=="radial") {
	os << "Radial weighting" << LogIO::POST;
    	  imwgt_p=VisImagingWeight("radial");
      }
      else{
    	  if(!imageDefined_p)
    		  throw(AipsError("Need to define image"));
    	  Int nx=itsMaxShape[0];
    	  Int ny=itsMaxShape[1];
    	  Quantity cellx=Quantity(itsMaxCoordSys.increment()[0], itsMaxCoordSys.worldAxisUnits()[0]);
    	  Quantity celly=Quantity(itsMaxCoordSys.increment()[1], itsMaxCoordSys.worldAxisUnits()[1]);
    	  if(type=="superuniform"){
    		  if(!imageDefined_p) throw(AipsError("Please define image first"));
    		  Int actualNpix=npixels;
    		  if(actualNpix <=0)
    			  actualNpix=3;
    		  os << LogIO::NORMAL // Loglevel INFO
    				  << "SuperUniform weighting over a square cell spanning ["
    				  << -actualNpix
    				  << ", " << actualNpix << "] in the uv plane" << LogIO::POST;
    		  imwgt_p=VisImagingWeight(*vi_p, rmode, noise, robust, nx,
    				  ny, cellx, celly, actualNpix,
    				  actualNpix, multiField);
    	  }
    	  else if ((type=="robust")||(type=="uniform")||(type=="briggs")) {
    		  if(!imageDefined_p) throw(AipsError("Please define image first"));
    		  Quantity actualFieldOfView(fieldofview);
    		  Int actualNPixels(npixels);
    		  String wtype;
    		  if(type=="briggs") {
    			  wtype = "Briggs";
    		  }
    		  else {
    			  wtype = "Uniform";
    		  }
    		  if(actualFieldOfView.get().getValue()==0.0&&actualNPixels==0) {
    			  actualNPixels=nx;
    			  actualFieldOfView=Quantity(actualNPixels*cellx.get("rad").getValue(),
    					  "rad");
    			  os << LogIO::NORMAL // Loglevel INFO
    					  << wtype
    					  << " weighting: sidelobes will be suppressed over full image"
    					  << LogIO::POST;
    		  }
    		  else if(actualFieldOfView.get().getValue()>0.0&&actualNPixels==0) {
    			  actualNPixels=nx;
    			  os << LogIO::NORMAL // Loglevel INFO
    					  << wtype
    					  << " weighting: sidelobes will be suppressed over specified field of view: "
    					  << actualFieldOfView.get("arcsec").getValue() << " arcsec" << LogIO::POST;
    		  }
    		  else if(actualFieldOfView.get().getValue()==0.0&&actualNPixels>0) {
    			  actualFieldOfView=Quantity(actualNPixels*cellx.get("rad").getValue(),
    					  "rad");
    			  os << LogIO::NORMAL // Loglevel INFO
    					  << wtype
    					  << " weighting: sidelobes will be suppressed over full image field of view: "
    					  << actualFieldOfView.get("arcsec").getValue() << " arcsec" << LogIO::POST;
    		  }
    		  else {
    			  os << LogIO::NORMAL // Loglevel INFO
    					  << wtype
    					  << " weighting: sidelobes will be suppressed over specified field of view: "
    					  << actualFieldOfView.get("arcsec").getValue() << " arcsec" << LogIO::POST;
    		  }
    		  os << LogIO::DEBUG1
    				  << "Weighting used " << actualNPixels << " uv pixels."
    				  << LogIO::POST;
    		  Quantity actualCellSize(actualFieldOfView.get("rad").getValue()/actualNPixels, "rad");

		  //		  cerr << "rmode " << rmode << " noise " << noise << " robust " << robust << " npixels " << actualNPixels << " cellsize " << actualCellSize << " multifield " << multiField << endl;
		  //		  Timer timer;
		  //timer.mark();
		  //Construct imwgt_p with old vi for now if old vi is in use as constructing with vi2 is slower 

		    ////This is slower by orders of magnitude as of 2014/06/25
		    imwgt_p=VisImagingWeight(*vi_p, rmode, noise, robust,
                                 actualNPixels, actualNPixels, actualCellSize,
                                 actualCellSize, 0, 0, multiField);

		  /*
		  if(rvi_p !=NULL){
		    imwgt_p=VisImagingWeight(*rvi_p, rmode, noise, robust,
                                 actualNPixels, actualNPixels, actualCellSize,
                                 actualCellSize, 0, 0, multiField);
		  }
		  else{
		    ////This is slower by orders of magnitude as of 2014/06/25
		    imwgt_p=VisImagingWeight(*vi_p, rmode, noise, robust,
                                 actualNPixels, actualNPixels, actualCellSize,
                                 actualCellSize, 0, 0, multiField);
		  }
		  */
		    //timer.show("After making visweight ");

    	  }
    	  else {
    		  //this->unlock();
    		  os << LogIO::SEVERE << "Unknown weighting " << type
    				  << LogIO::EXCEPTION;
    		  return false;
    	  }
      }
	 
	 //// UV-Tapering
	 //cout << "Taper type : " << filtertype << " : " << (filtertype=="gaussian") <<  endl;
	 if( filtertype == "gaussian" ) {
	   //	   os << "Setting uv-taper" << LogIO::POST;
	   imwgt_p.setFilter( filtertype,  filterbmaj, filterbmin, filterbpa );
	 }
	 
	 vi_p->useImagingWeight(imwgt_p);
	 
	 return true;
	 
       }
       catch(AipsError &x)
	 {
	   throw( AipsError("Error in Weighting : "+x.getMesg()) );
	 }
       
       
       return true;
  }
  
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////    Internal Functions start here.  These are not visible to the tool layer.
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////
  ////////////This should be called  at each defineimage
  CountedPtr<SIImageStore> SynthesisImager::createIMStore(String imageName, 
							  CoordinateSystem& cSys,
							  IPosition imShape, 
							  const Bool overwrite,
							  String mappertype,
							  uInt ntaylorterms,
							  Quantity distance,
							  uInt facets,
							  Bool useweightimage)
  {
    LogIO os( LogOrigin("SynthesisImager","createIMStore",WHERE) );

    CountedPtr<SIImageStore> imstor;

    try
      {
	// Prepare miscellaneous image information
	auto objectName = msc.field().name()(msc.fieldId()(0));
	///// misc info fpr ImageStore. This will go to the 'miscinfo' table keyword
	Record miscInfo;
	auto telescop=msc.observation().telescopeName()(0);
	miscInfo.define("INSTRUME", telescop);
	miscInfo.define("distance", distance.get("m").getValue());

	if( mappertype=="default" || mappertype=="imagemosaic" )
	  {
            imstor = std::make_shared<SIImageStore>(imageName, cSys, imShape, objectName,
                                                    miscInfo, overwrite,
                                                    (useweightimage || (mappertype=="imagemosaic")
                                                     ));
	  }
	else if (mappertype == "multiterm" )  // Currently does not support imagemosaic.
	  {
            // upcast with shared_ptr and then assign to CountedPtr<SIImageStore>
            std::shared_ptr<SIImageStore> multiTermStore =
                std::make_shared<SIImageStoreMultiTerm>(imageName, cSys, imShape,
                                                        objectName, miscInfo, facets,
                                                        overwrite, ntaylorterms, useweightimage);
            imstor = multiTermStore;
	  }
	else
	  {
	    throw(AipsError("Internal Error : Invalid mapper type in SynthesisImager::createIMStore"));
	  }

	// Get polRep from 'msc' here, and send to imstore. 
	StokesImageUtil::PolRep polRep(StokesImageUtil::CIRCULAR);
	Vector<String> polType=msc.feed().polarizationType()(0);
	if (polType(0)!="X" && polType(0)!="Y" &&  polType(0)!="R" && polType(0)!="L") {
	  os << LogIO::WARN << "Unknown stokes types in feed table: ["
	     << polType(0) << ", " << polType(1) << "]" << endl
	     << "Results open to question!" << LogIO::POST;
	}
	
	if (polType(0)=="X" || polType(0)=="Y") {
	  polRep=StokesImageUtil::LINEAR;
	  os << LogIO::DEBUG1 << "Preferred polarization representation is linear" << LogIO::POST;
	}
	else {
	  polRep=StokesImageUtil::CIRCULAR;
	  os << LogIO::DEBUG1 << "Preferred polarization representation is circular" << LogIO::POST;
	}
	/// end of reading polRep info
	
	///////// Send this info into ImageStore.
	imstor->setDataPolFrame(polRep);

      }
    catch(AipsError &x)
      {
	throw(AipsError("Error in createImStore : " + x.getMesg() ) );
      }
    
    
    return imstor;
  }
  
  CountedPtr<SIMapper> SynthesisImager::createSIMapper(String mappertype,  
							   CountedPtr<SIImageStore> imagestore,
						       CountedPtr<refim::FTMachine> ftmachine,
						       CountedPtr<refim::FTMachine> iftmachine,
						       uInt /*ntaylorterms*/)
  {
    LogIO os( LogOrigin("SynthesisImager","createSIMapper",WHERE) );
    
    CountedPtr<SIMapper> localMapper=NULL;

    try
      {
	
	if( mappertype == "default" || mappertype == "multiterm" )
	  {
	    localMapper = new SIMapper( imagestore, ftmachine, iftmachine );
	  }
	else if( mappertype == "imagemosaic") // || mappertype == "mtimagemosaic" )
	  {
	    localMapper = new SIMapperImageMosaic( imagestore, ftmachine, iftmachine );
	  }
	else
	  {
	    throw(AipsError("Unknown mapper type : " + mappertype));
	  }

      }
    catch(AipsError &x) {
	throw(AipsError("Error in createSIMapper : " + x.getMesg() ) );
      }
    return localMapper;
  }
  

  Block<CountedPtr<SIImageStore> > SynthesisImager::createFacetImageStoreList(
									      CountedPtr<SIImageStore> imagestore,
									      Int facets)
  {
    Block<CountedPtr<SIImageStore> > facetList( facets*facets );

    if( facets==1 ) { facetList[0] = imagestore;  return facetList; }

    // Remember, only the FIRST field in each run can have facets. So, check for this.
    if( ! unFacettedImStore_p.null() ) {
	throw( AipsError("A facetted image has already been set. Facets are supported only for the main (first) field. Please submit a feature-request if you need multiple facets for outlier fields as well. ") );
      }
    
    unFacettedImStore_p = imagestore;
    facetsStore_p = facets;
    
    // Note : facets : Number of facets on a side.
    // Note : facet : index from range(0, facets*facets)
    for (Int facet=0; facet< facets*facets; ++facet){
	facetList[facet] = unFacettedImStore_p->getSubImageStore(facet, facets);
      }
    
    return facetList;
  }

  void SynthesisImager::setPsfFromOneFacet()
  {

    if( unFacettedImStore_p.null() ){
	throw(AipsError("Internal Error in SynthesisImager : Setting PSF on Null unfacettedimage"));
      }

    // Copy the PSF from one facet to the center of the full image, to use for the minor cycle
    //
    // This code segment will work for single and multi-term 
    // - for single term, the index will always be 0, and SIImageStore's access functions know this.
    // - for multi-term, the index will be the taylor index and SIImageStoreMultiTerm knows this.
      {
	IPosition shape=(unFacettedImStore_p->psf(0))->shape();
	IPosition blc(4, 0, 0, 0, 0);
	IPosition trc=shape-1;
	for(uInt tix=0; tix<2 * unFacettedImStore_p->getNTaylorTerms() - 1; tix++)
	  {
	    TempImage<Float> onepsf((itsMappers.imageStore(0)->psf(tix))->shape(), 
				    (itsMappers.imageStore(0)->psf(tix))->coordinates());
	    onepsf.copyData(*(itsMappers.imageStore(0)->psf(tix)));
	    //now set the original to 0 as we have a copy of one facet psf
	    (unFacettedImStore_p->psf(tix))->set(0.0);
	    blc[0]=(shape[0]-(onepsf.shape()[0]))/2;
	    trc[0]=onepsf.shape()[0]+blc[0]-1;
	    blc[1]=(shape[1]-(onepsf.shape()[1]))/2;
	    trc[1]=onepsf.shape()[1]+blc[1]-1;
	    Slicer sl(blc, trc, Slicer::endIsLast);
	    SubImage<Float> sub(*(unFacettedImStore_p->psf(tix)), sl, true);
	    sub.copyData(onepsf);
	  }
      }

      //cout << "In setPsfFromOneFacet : sumwt : " << unFacettedImStore_p->sumwt()->get() << endl;

  }
  
  
  void SynthesisImager::appendToMapperList(String imagename,  
					   CoordinateSystem& csys, 
					   IPosition imshape,
					   CountedPtr<FTMachine>& ftm,
					   CountedPtr<FTMachine>& iftm,
					   Quantity distance, 
					   Int facets, 
					   const Bool overwrite,
					   String mappertype,
					   uInt ntaylorterms  )
    {
      LogIO log_l(LogOrigin("SynthesisImager", "appendToMapperList(ftm)"));
      //---------------------------------------------
      // Some checks..
      if(facets > 1 && itsMappers.nMappers() > 0)
	log_l << "Facetted image has to be first of multifields" << LogIO::EXCEPTION;
      
      AlwaysAssert( ( ( ! (ftm->name()=="MosaicFT" && mappertype=="imagemosaic") )  && 
      		      ( ! (ftm->name()=="AWProjectWBFTNew" && mappertype=="imagemosaic") )) ,
		    AipsError );
      //---------------------------------------------

      // Create the ImageStore object
      CountedPtr<SIImageStore> imstor;
      imstor = createIMStore(imagename, csys, imshape, overwrite,mappertype, ntaylorterms, distance,facets, iftm->useWeightImage() );

      // Create the Mappers
      if( facets<2 ) // One facet. Just add the above imagestore to the mapper list.
	{
	  itsMappers.addMapper(  createSIMapper( mappertype, imstor, ftm, iftm, ntaylorterms) );
	}
      else // This field is facetted. Make a list of reference imstores, and add all to the mapper list.
	{
	  // First, make sure that full images have been allocated before trying to make references.....
	  //	  if( ! imstor->checkValidity(true/*psf*/, true/*res*/,true/*wgt*/,true/*model*/,false/*image*/,false/*mask*/,true/*sumwt*/ ) ) 
	  //	    { throw(AipsError("Internal Error : Invalid ImageStore for " + imstor->getName())); }

	  // Make and connect the list.
	  Block<CountedPtr<SIImageStore> > imstorList = createFacetImageStoreList( imstor, facets );
	  for( uInt facet=0; facet<imstorList.nelements(); facet++)
	    {
	      CountedPtr<FTMachine> new_ftm, new_iftm;
	      if(facet==0){ new_ftm = ftm;  new_iftm = iftm; }
	      else{ new_ftm=ftm->cloneFTM();  new_iftm=iftm->cloneFTM(); }
	      itsMappers.addMapper(createSIMapper( mappertype, imstorList[facet], new_ftm, new_iftm, ntaylorterms));
	    }
	}

    }

  /////////////////////////
  /*
  Bool SynthesisImager::toUseWeightImage(CountedPtr<FTMachine>& ftm, String mappertype)
  {
    if( (ftm->name() == "GridFT" || ftm->name() == "WProjectFT")&&(mappertype!="imagemosaic") )  
      { return false; }
    else
      { return true; }
  }
  */

  // Make the FT-Machine and related objects (cfcache, etc.)
  void SynthesisImager::createFTMachine(CountedPtr<FTMachine>& theFT, 
					CountedPtr<FTMachine>& theIFT, 
					const String& ftname,
					const uInt nTaylorTerms,
					const String mType,
					const Int facets,            //=1
					//------------------------------
					const Int wprojplane,        //=1,
					const Float padding,         //=1.0,
					const Bool useAutocorr,      //=false,
					const Bool useDoublePrec,    //=true,
					const String gridFunction,   //=String("SF"),
					//------------------------------
					const Bool aTermOn,          //= true,
					const Bool psTermOn,         //= true,
					const Bool mTermOn,          //= false,
					const Bool wbAWP,            //= true,
					const String cfCache,        //= "",
					const Bool doPointing,       //= false,
					const Bool doPBCorr,         //= true,
					const Bool conjBeams,        //= true,
					const Float computePAStep,         //=360.0
					const Float rotatePAStep,          //=5.0
					const String interpolation,  //="linear"
					const Bool freqFrameValid, //=true
					const Int cache,             //=1000000000,
					const Int tile               //=16
					)

  {
    LogIO os( LogOrigin("SynthesisImager","createFTMachine",WHERE));

    if(ftname=="gridft"){
      if(facets >1){
	theFT=new refim::GridFT(cache, tile, gridFunction, mLocation_p, phaseCenter_p, padding, useAutocorr, useDoublePrec);
	theIFT=new refim::GridFT(cache, tile, gridFunction, mLocation_p, phaseCenter_p, padding, useAutocorr, useDoublePrec);

      }
      else{
	theFT=new refim::GridFT(cache, tile, gridFunction, mLocation_p, padding, useAutocorr, useDoublePrec);
	theIFT=new refim::GridFT(cache, tile, gridFunction, mLocation_p, padding, useAutocorr, useDoublePrec);
      }
    }
    else if(ftname== "wprojectft"){
      theFT=new refim::WProjectFT(wprojplane,  mLocation_p,
			   cache/2, tile, useAutocorr, padding, useDoublePrec);
      theIFT=new refim::WProjectFT(wprojplane,  mLocation_p,
			   cache/2, tile, useAutocorr, padding, useDoublePrec);
      CountedPtr<WPConvFunc> sharedconvFunc= new WPConvFunc();
      static_cast<refim::WProjectFT &>(*theFT).setConvFunc(sharedconvFunc);
      static_cast<refim::WProjectFT &>(*theFT).setConvFunc(sharedconvFunc);
    }
    else if ((ftname == "awprojectft") || (ftname== "mawprojectft") || (ftname == "protoft")) {
      createAWPFTMachine(theFT, theIFT, ftname, facets, wprojplane, 
			 padding, useAutocorr, useDoublePrec, gridFunction,
			 aTermOn, psTermOn, mTermOn, wbAWP, cfCache, 
			 doPointing, doPBCorr, conjBeams, computePAStep,
			 rotatePAStep, cache,tile);
    }
    else
      {
	throw( AipsError( "Invalid FTMachine name : " + ftname ) );
      }
    /* else if(ftname== "MosaicFT"){

       }*/

    ///////// Now, clone and pack the chosen FT into a MultiTermFT if needed.
    if( mType=="multiterm" )
      {
	AlwaysAssert( nTaylorTerms>=1 , AipsError );

	CountedPtr<refim::FTMachine> theMTFT = new refim::MultiTermFTNew( theFT , nTaylorTerms, true/*forward*/ );
	CountedPtr<refim::FTMachine> theMTIFT = new refim::MultiTermFTNew( theIFT , nTaylorTerms, false/*forward*/ );

	theFT = theMTFT;
	theIFT = theMTIFT;
      }


    ////// Now, set the SkyJones if needed, and if not internally generated.
    if( mType=="imagemosaic" && 
	(ftname != "awprojectft" && ftname != "mawprojectft" && ftname != "proroft") )
      {
	CountedPtr<SkyJones> vp = NULL;
	MSColumns msc(*mss_p[0]);
	Quantity parang(0.0,"deg");
	Quantity skyposthreshold(0.0,"deg");
	vp = new VPSkyJones(msc, true,  parang, BeamSquint::NONE,skyposthreshold);

	Vector<CountedPtr<SkyJones> > skyJonesList(1);
	skyJonesList(0) = vp;
	theFT->setSkyJones(  skyJonesList );
	theIFT->setSkyJones(  skyJonesList );

      }

    //// For mode=cubedata, set the freq frame to invalid..
    // get this info from buildCoordSystem
    Vector<Int> tspws(0);
    //theFT->setSpw( tspws, false );
    //theIFT->setSpw( tspws, false );
    theFT->setSpw( tspws, freqFrameValid );
    theIFT->setSpw( tspws, freqFrameValid );

    //// Set interpolation mode
    theFT->setFreqInterpolation( interpolation );
    theIFT->setFreqInterpolation( interpolation );

  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void SynthesisImager::createAWPFTMachine(CountedPtr<FTMachine>& theFT, CountedPtr<FTMachine>& theIFT, 
					   const String& ftmName,
					   const Int,// facets,            //=1
					   //------------------------------
					   const Int wprojPlane,        //=1,
					   const Float,// padding,         //=1.0,
					   const Bool,// useAutocorr,      //=false,
					   const Bool useDoublePrec,    //=true,
					   const String,// gridFunction,   //=String("SF"),
					   //------------------------------
					   const Bool aTermOn,          //= true,
					   const Bool psTermOn,         //= true,
					   const Bool mTermOn,          //= false,
					   const Bool wbAWP,            //= true,
					   const String cfCache,        //= "",
					   const Bool doPointing,       //= false,
					   const Bool doPBCorr,         //= true,
					   const Bool conjBeams,        //= true,
					   const Float computePAStep,   //=360.0
					   const Float rotatePAStep,    //=5.0
					   const Int cache,             //=1000000000,
					   const Int tile               //=16
					)

  {
    LogIO os( LogOrigin("SynthesisImager2","createAWPFTMachine",WHERE));

    if (wprojPlane<=1)
      {
	os << LogIO::NORMAL
	   << "You are using wprojplanes=1. Doing co-planar imaging (no w-projection needed)" 
	   << LogIO::POST;
	os << LogIO::NORMAL << "Performing WBA-Projection" << LogIO::POST; // Loglevel PROGRESS
      }
    if((wprojPlane>1)&&(wprojPlane<64)) 
      {
	os << LogIO::WARN
	   << "No. of w-planes set too low for W projection - recommend at least 128"
	   << LogIO::POST;
	os << LogIO::NORMAL << "Performing WBAW-Projection" << LogIO::POST; // Loglevel PROGRESS
      }

    CountedPtr<refim::ATerm> apertureFunction = createTelescopeATerm(*mss_p[0], aTermOn);
    CountedPtr<refim::PSTerm> psTerm = new refim::PSTerm();
    CountedPtr<refim::WTerm> wTerm = new refim::WTerm();
    
    //
    // Selectively switch off CFTerms.
    //
    if (aTermOn == false) {apertureFunction->setOpCode(refim::CFTerms::NOOP);}
    if (psTermOn == false) psTerm->setOpCode(refim::CFTerms::NOOP);

    //
    // Construct the CF object with appropriate CFTerms.
    //
    CountedPtr<refim::ConvolutionFunction> awConvFunc;
    //    awConvFunc = new AWConvFunc(apertureFunction,psTerm,wTerm, !wbAWP);
    if ((ftmName=="mawprojectft") || (mTermOn))
      awConvFunc = new refim::AWConvFuncEPJones(apertureFunction,psTerm,wTerm,wbAWP);
    else
      awConvFunc = new refim::AWConvFunc(apertureFunction,psTerm,wTerm,wbAWP);

    //
    // Construct the appropriate re-sampler.
    //
    CountedPtr<refim::VisibilityResamplerBase> visResampler;
    //    if (ftmName=="protoft") visResampler = new ProtoVR();
    //elsef
    visResampler = new refim::AWVisResampler();
    //    CountedPtr<VisibilityResamplerBase> visResampler = new VisibilityResampler();

    //
    // Construct and initialize the CF cache object.
    //
    CountedPtr<refim::CFCache> cfCacheObj = new refim::CFCache();
    cfCacheObj->setCacheDir(cfCache.data());
    cfCacheObj->initCache2();

    //
    // Finally construct the FTMachine with the CFCache, ConvFunc and
    // Re-sampler objects.  
    //
    Float pbLimit_l=1e-3;
    theFT = new refim::AWProjectWBFTNew(wprojPlane, cache/2, 
			      cfCacheObj, awConvFunc, 
			      visResampler,
			      /*true */doPointing, doPBCorr, 
			      tile, computePAStep, pbLimit_l, true,conjBeams,
			      useDoublePrec);

    Quantity rotateOTF(rotatePAStep,"deg");
    static_cast<refim::AWProjectWBFTNew &>(*theFT).setObservatoryLocation(mLocation_p);
    static_cast<refim::AWProjectWBFTNew &>(*theFT).setPAIncrement(Quantity(computePAStep,"deg"),rotateOTF);

    // theIFT = new AWProjectWBFT(wprojPlane, cache/2, 
    // 			       cfCacheObj, awConvFunc, 
    // 			       visResampler,
    // 			       /*true */doPointing, doPBCorr, 
    // 			       tile, computePAStep, pbLimit_l, true,conjBeams,
    // 			       useDoublePrec);

    // static_cast<AWProjectWBFT &>(*theIFT).setObservatoryLocation(mLocation_p);
    // static_cast<AWProjectWBFT &>(*theIFT).setPAIncrement(Quantity(computePAStep,"deg"),rotateOTF);

    theIFT = new refim::AWProjectWBFTNew(static_cast<refim::AWProjectWBFTNew &>(*theFT));

    //// Send in Freq info.
    //    vi_p->getFreqInSpwRange(
    os << "DO SOMETHING HERE for vi2/vb2 : Sending frequency selection information " <<  mssFreqSel_p  <<  " to AWP FTM." << LogIO::POST;
        theFT->setSpwFreqSelection( mssFreqSel_p );
        theIFT->setSpwFreqSelection( mssFreqSel_p );

  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  refim::ATerm* SynthesisImager::createTelescopeATerm(const MeasurementSet& ms, const Bool& isATermOn)
  {
    LogIO os(LogOrigin("SynthesisImager", "createTelescopeATerm",WHERE));
    
    if (!isATermOn) return new refim::NoOpATerm();
    
    MSObservationColumns msoc(ms.observation());
    String ObsName=msoc.telescopeName()(0);
    if ((ObsName == "EVLA") || (ObsName == "VLA"))
      return new refim::EVLAAperture();
    else
      {
	os << "Telescope name ('"+
	  ObsName+"') in the MS not recognized to create the telescope specific ATerm" 
	   << LogIO::WARN;
      }
    
    return NULL;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  // Do MS-Selection and set up vi/vb. 
  // Only this functions needs to know anything about the MS 
  void SynthesisImager::createVisSet(const Bool writeAccess)
  {
    LogIO os( LogOrigin("SynthesisImager","createVisSet",WHERE) );

    if(mss_p.nelements() != uInt(fselections_p.size()) && (fselections_p.size() !=0)){
      throw(AipsError("Discrepancy between Number of MSs and Frequency selections"));
    }
    vi_p=new vi::VisibilityIterator2(mss_p, vi::SortColumns(), writeAccess);
    if(fselections_p.size() !=0)
      vi_p->setFrequencySelection (fselections_p);
    //return *vi_p;
  }// end of createVisSet

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  void SynthesisImager::runMajorCycle(const Bool dopsf, 
				      const Bool savemodel)
  {
    LogIO os( LogOrigin("SynthesisImager","runMajorCycle",WHERE) );

    //    cout << "Savemodel : " << savemodel << "   readonly : " << readOnly_p << "   usescratch : " << useScratch_p << endl;

    Bool savemodelcolumn = savemodel && !readOnly_p && useScratch_p;
    Bool savevirtualmodel = savemodel && !readOnly_p && !useScratch_p;

    if( savemodelcolumn ) os << "Saving model column" << LogIO::POST;
    if( savevirtualmodel ) os << "Saving virtual model" << LogIO::POST;

    itsMappers.checkOverlappingModels("blank");

    {
    	vi_p->originChunks();
    	vi_p->origin();
    	vi::VisBuffer2* vb=vi_p->getVisBuffer();
    	if(!dopsf) itsMappers.initializeDegrid(*vb);
    	itsMappers.initializeGrid(*vb,dopsf);
    	for (vi_p->originChunks(); vi_p->moreChunks();vi_p->nextChunk())
    	{
    		for (vi_p->origin(); vi_p->more();vi_p->next())
    		{
    			if(!dopsf){
    				vb->setVisCubeModel(Cube<Complex>(vb->visCubeModel().shape(), Complex(0.0, 0.0)));
    				itsMappers.degrid(*vb, savevirtualmodel );
    				if( savemodelcolumn && vi_p->isWritable())
    					vi_p->writeVisModel(vb->visCubeModel());
    			}
    			itsMappers.grid(*vb, dopsf, datacol_p);
    		}
    	}
    	if(!dopsf) itsMappers.finalizeDegrid(*vb);
    	itsMappers.finalizeGrid(*vb, dopsf);
    }
    itsMappers.checkOverlappingModels("restore");

  }// end runMajorCycle

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}//# NAMESPACE refim - END
} //# NAMESPACE CASA - END

