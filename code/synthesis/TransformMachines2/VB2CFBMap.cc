//-*- C++ -*-
//# VB2CFBMap.h: Definition of the VB2CFBMap class
//# Copyright (C) 1997,1998,1999,2000,2001,2002,2003
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
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id$
//
#include <synthesis/TransformMachines/SynthesisMath.h>
#include <casa/Logging/LogIO.h>
#include <casa/Logging/LogSink.h>
#include <casa/Logging/LogOrigin.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Array.h>
#include <casa/Utilities/CountedPtr.h>
#include <synthesis/TransformMachines2/CFStore2.h>
#include <synthesis/TransformMachines2/ConvolutionFunction.h>
#include <synthesis/TransformMachines2/CFBuffer.h>
#include <msvis/MSVis/VisBuffer2.h>
#include <synthesis/TransformMachines2/VB2CFBMap.h>
#include <synthesis/TransformMachines2/Utils.h>
namespace casa{
  using namespace vi;
  namespace refim{
    Int mapAntIDToAntType(const casacore::Int& /*ant*/) {return 0;};

    VB2CFBMap::VB2CFBMap(): vb2CFBMap_p(), cfPhaseGrad_p(), baselineType_p(), vectorPhaseGradCalculator_p(), doPointing_p(false), cachedFieldId_p(-1), vbRow2BLMap_p(), vbRows_p(0), sigmaDev(), timer_p()
    {
      baselineType_p = new BaselineType();
      needsNewPOPG_p = false;
      needsNewFieldPG_p = false;
      totalCost_p=totalVB_p = 0.0;
      blNeedsNewPOPG_p.resize(0);
      blType_p=0;
      vectorPhaseGradCalculator_p.resize(0);
      // sigmaDev = SynthesisUtils::getenv("PO_SIGMADEV",3.0);
    };

    VB2CFBMap& VB2CFBMap::operator=(const VB2CFBMap& other)
    {
      if(this!=&other) 
	{
	  baselineType_p = other.baselineType_p;
	  cfPhaseGrad_p.assign(other.cfPhaseGrad_p);
	  vb2CFBMap_p.assign(other.vb2CFBMap_p);
	  vbRow2BLMap_p = other.vbRow2BLMap_p;
	  vectorPhaseGradCalculator_p.assign(other.vectorPhaseGradCalculator_p);
	  doPointing_p = other.doPointing_p;
	}
      return *this;
    };

    // void VB2CFBMap::setPhaseGradPerRow(const CountedPtr<PointingOffsets>& pointingOffset,
    // 				       const casacore::CountedPtr<CFBuffer>& cfb,
    // 				       const vi::VisBuffer2& vb,
    // 				       const int& row)
    // {
    //   // if (doPointing_p)
    //   // 	{
    //   // 	  if (phaseGradCalculator_p->needsNewPhaseGrad(pointingOffset, vb, 0))
    //   // 	    {
    //   // 	      phaseGradCalculator_p->ComputeFieldPointingGrad(pointingOffset,cfb,vb, 0);
    //   // 	      newPhaseGradComputed_p=true;
    //   // 	    }
    //   // 	}
    //   // else
    // 	{
    // 	  baselineType_p->phaseGradCalculator_p->ComputeFieldPointingGrad(pointingOffset,cfb,vb, 0);
    // 	}
    // 	{
    // 	  //cfPhaseGrad_p(row).assign(phaseGradCalculator_p->getFieldPointingGrad());
    // 	  cfPhaseGrad_p(row).reference(baselineType_p->phaseGradCalculator_p->field_phaseGrad_p);
    // 	}
    // }

    void VB2CFBMap::setBLNeedsNewPOPG(vector<int>& vbRow2BLMap_p)
    {
      blNeedsNewPOPG_p.resize(vbRow2BLMap_p.size(),false);
      vector<int> uniqueVbRow2BLMap_p = vbRow2BLMap_p;
      sort(uniqueVbRow2BLMap_p.begin(),uniqueVbRow2BLMap_p.end());
      auto itrBLMap = unique(uniqueVbRow2BLMap_p.begin(),uniqueVbRow2BLMap_p.end());
      uniqueVbRow2BLMap_p.erase(itrBLMap,uniqueVbRow2BLMap_p.end());
      
      for(int ii=0; ii < uniqueVbRow2BLMap_p.size(); ii++)
	{
	  // cerr << "unqVb2BL: "<< uniqueVbRow2BLMap_p[ii] << " ii " << ii << endl;
	  int idx = baselineType_p->returnIdx(vbRow2BLMap_p, uniqueVbRow2BLMap_p[ii]);
	  blNeedsNewPOPG_p[idx] = true;
	}

      
    }



    //______________________________________________________
    
    Int VB2CFBMap::makeVBRow2CFBMap(CFStore2& cfs,
				    const VisBuffer2& vb, 
				    const Quantity& dPA,
				    const Vector<Int>& /*dataChan2ImChanMap*/,
				    const Vector<Int>& /*dataPol2ImPolMap*/,
				    const CountedPtr<PointingOffsets>& pointingOffsets_p)
    //const Bool& /*doPointing*/)
    {
      //    VBRow2CFMapType& vbRow2CFMap_p,
      const Int nRow=vb.nRows(); 
      //UNUSED: nChan=dataChan2ImChanMap.nelements(), 
      //UNUSED: nPol=dataPol2ImPolMap.nelements();
      //    vbRow2CFMap_p.resize(nPol, nChan, nRow);
      vb2CFBMap_p.resize(nRow);
      cfPhaseGrad_p.resize(nRow);

      if(cachedFieldId_p != vb.fieldId()[0])
	{
	  needsNewFieldPG_p = true;
	  needsNewPOPG_p = true;
	  baselineType_p->setCachedGroups(false);
	  cachedFieldId_p = vb.fieldId()[0];
	  vectorPhaseGradCalculator_p.resize(0);
	}
      else if (vbRows_p == 0)
	{
	  vbRows_p = vb.nRows();
	  baselineType_p->setCachedGroups(false);
	  needsNewPOPG_p = true;
	}
      else if (vbRows_p != vb.nRows())
	{
	  vbRows_p = vb.nRows();
	  baselineType_p->setCachedGroups(false);
	  needsNewPOPG_p = true;
	}
      else
	{
	  baselineType_p->setCachedGroups(true);
	   needsNewPOPG_p = false;
	}
      

      Quantity pa(getPA(vb),"rad");
      //PolOuterProduct outerProduct;
      Int statusCode=CFDefs::MEMCACHE;

      // baselineType_p->setCacheGroups(vbRows_p, vb);
      // baselineType_p->setCachedGroups(true);
      baselineType_p->setDoPointing(doPointing_p);


      if(doPointing_p)
	{
	  // cerr<<"VB2CFBMap::makeVBRow2CFBMap DoP = T"<<endl;
	  Float A2R = 4.848137E-06;
	  Vector<Double> poIncrement = pointingOffsets_p->getIncrement();
	  // Double offsetDeviation = sigmaDev * A2R / sqrt(poIncrement[0]*poIncrement[0] + poIncrement[1]*poIncrement[1]);
	  baselineType_p->findAntennaGroupsM(vb,pointingOffsets_p, sigmaDev);

	  baselineType_p->makeVBRow2BLGMap(vb);
	  vbRow2BLMap_p = baselineType_p->getVBRow2BLMap();
	  setBLNeedsNewPOPG(vbRow2BLMap_p);



	  // cerr << "vbRow2BLMap_p" ;
	  // for(int i=0; i <= vbRow2BLMap_p.size(); i++)
	  //   cerr << vbRow2BLMap_p[i]<<" ";
	  // cerr <<endl;
	  // cerr << "baselineType->cachedGroups_p : "<< baselineType_p->cachedGroups_p <<endl;
	  if(baselineType_p->getCachedGroups())
	    {
	      unsigned int cachedPOSize_l = (baselineType_p->getCachedAntennaPO()).size();
	      unsigned int poSize_l = (pointingOffsets_p->pullPointingOffsets()).size();
	      Vector<Double> sumResPO_l;
	      double avgResPO_l = 0.0;
	      sumResPO_l.resize(2);
	      sumResPO_l[0]=0;
	      sumResPO_l[1]=0;
	      // cerr << "sumResPO_l "<< sumResPO_l << " poSize_l " << poSize_l << " cachedPOSize_l"<< cachedPOSize_l <<endl;
	      if(poSize_l == cachedPOSize_l)
		{
		  Vector < Vector <Double> > residualPointingOffsets_l = baselineType_p->getCachedAntennaPO() - pointingOffsets_p->pullPointingOffsets(); 
		  for(unsigned int ii=0; ii < poSize_l; ii++) 
		    sumResPO_l = sumResPO_l + residualPointingOffsets_l[ii];
		  avgResPO_l = sqrt(sumResPO_l[0]*sumResPO_l[0] + sumResPO_l[1]*sumResPO_l[1])/poSize_l;
		  // cerr << "avgResPO_l"<< avgResPO_l <<endl;
		  if(avgResPO_l*sqrt(poIncrement[0]*poIncrement[0] + poIncrement[1]*poIncrement[1]) >= sigmaDev*A2R)
		    {
		      baselineType_p->setCachedGroups(false);
		      baselineType_p->findAntennaGroupsM(vb,pointingOffsets_p, sigmaDev);
		      baselineType_p->setCachedAntennaPO(pointingOffsets_p->pullPointingOffsets());
		      baselineType_p->makeVBRow2BLGMap(vb);
		      vbRow2BLMap_p = baselineType_p->getVBRow2BLMap();
		      setBLNeedsNewPOPG(vbRow2BLMap_p);
		      if (!needsNewFieldPG_p)
			{
			  LogIO log_l(LogOrigin("VB2CFBMap", "makeVBRow2CFBMap[R&D]"));

			  log_l << "The average antenna Offset : " << avgResPO_l*sqrt(poIncrement[0]*poIncrement[0] + poIncrement[1]*poIncrement[1])/A2R 
				<< " arcsec, exceeds sigmaDev : " << sigmaDev 
				<< " Field ID = " << vb.fieldId()(0) << LogIO::POST;
			  
			}
		      cerr << "Avg of the Residual Pointing Offsets " << avgResPO_l*sqrt(poIncrement[0]*poIncrement[0] + poIncrement[1]*poIncrement[1])/A2R <<endl; 
		    }
		}
	      else
		{
		  baselineType_p->setCachedGroups(false);
		  needsNewPOPG_p = true;
		}
	      
	    }
	  baselineType_p->setCachedAntennaPO(pointingOffsets_p->pullPointingOffsets());
	}

      for (Int irow=0;irow<nRow;irow++)
	{
	  //
	  // Translate antenna ID to antenna type
	  //
	  Int ant1Type = mapAntIDToAntType(vb.antenna1()(irow)),
	    ant2Type = mapAntIDToAntType(vb.antenna2()(irow));
	  //
	  // Get the CFBuffer for the given PA and baseline catagorized
	  // by the two antenna types.  For homgeneous arrays, all
	  // baselines will map to a single antenna-type pair.
	  //
	  
	  CountedPtr<CFBuffer> cfb_l;
	  try
	    {
	      cfb_l = cfs.getCFBuffer(pa, dPA, ant1Type, ant2Type);
	      //cfb_l->show("From VRB: ");
	    }
	  catch (CFNotCached& x)
	    {
	      LogIO log_l(LogOrigin("VB2CFBMap", "makeVBRow2CFBMap[R&D]"));

	      log_l << "CFs not cached for " << pa.getValue("deg") 
		    << " deg, dPA = " << dPA.getValue("deg") 
		    << " Field ID = " << vb.fieldId()(0);
	      log_l << " Ant1Type, Ant2Type = " << ant1Type << "," << ant2Type << LogIO::POST;
	      statusCode=CFDefs::NOTCACHED;
	    }
	  
	  if (statusCode==CFDefs::NOTCACHED)
	    {
	      break;
	    }
	  else
	    {
	      // Set the phase grad for the CF per VB row
	      // setPhaseGradPerRow(pointingOffsets_p, cfb_l, vb, irow);
	      timer_p.mark();
	      // if (vbRows_p == 0)
	      //   {
	      //     vbRows_p = vb.nRows();
	      //     baselineType_p->setcachedGroups(false);
	      //   }
	      // else if (vbRows_p != vb.nRows())
	      //   {
	      //     vbRows_p = vb.nRows();
	      //     baselineType_p->setcachedGroups(false);
	      //   }
	      // else
	      //   baselineType_p->setcachedGroups_p(true);
	      // baselineType_p->setCacheGroups(vbRows_p, vb);
	      // baselineType_p->setDoPointing(doPointing_p);
	      // if(computeAntennaGroups_p)
	      // baselineType_p->findAntennaGroups(vb,pointingOffsets_p,sigmaDev);
	      cfPhaseGrad_p(irow).reference(setBLPhaseGrad(pointingOffsets_p, cfb_l, vb, irow, sigmaDev));
	      totalCost_p += timer_p.real();
	      totalVB_p++;
	      // Set the CFB per VB row
	      cfb_l->setPointingOffset(pointingOffsets_p->pullPointingOffsets());
	      vb2CFBMap_p(irow) = cfb_l;
	      vbRows_p = vb.nRows();
	      // if( needsNewPOPG_p ||  needsNewFieldPG_p)
	      // 	cerr << "NeedsNewPOPG_p " << needsNewPOPG_p << " NeedsNewFieldPG_p "<< needsNewFieldPG_p << " row " << myrow << endl;
		  
	    }
	}
      // {
      // 	double n=0;
      // 	for (int i=0;i<cfPhaseGrad_p.nelements();i++)
      // 	  n+=cfPhaseGrad_p[i].shape().product()*sizeof(casacore::Complex);
      // 	log_l << "Size of VB2CFBMap::cfPhaseGrad_p = " << n << " bytes" << LogIO::POST;
      // }
      return statusCode;
    }

    Matrix<Complex> VB2CFBMap::setBLPhaseGrad(const CountedPtr<PointingOffsets>& pointingOffsets_p ,
						 const CountedPtr<CFBuffer>& cfb,
						 const vi::VisBuffer2& vb,
						 const int& row,
					      const double& /*sigmaDev*/)
    {
      int myrow=row;
      int idx;
      if(doPointing_p)
	{
	  vector<int> uniqueVbRow2BLMap_p = vbRow2BLMap_p;
	  sort(uniqueVbRow2BLMap_p.begin(),uniqueVbRow2BLMap_p.end());
	  auto itrBLMap = unique(uniqueVbRow2BLMap_p.begin(),uniqueVbRow2BLMap_p.end());
	  uniqueVbRow2BLMap_p.erase(itrBLMap,uniqueVbRow2BLMap_p.end());
	  int maxVB2BLMap = *max_element(uniqueVbRow2BLMap_p.begin(),uniqueVbRow2BLMap_p.end());
	  if (myrow==0 )
	    {
	      blType_p = vbRow2BLMap_p[myrow];
	      blNeedsNewPOPG_p[myrow] = true;
	      // cerr <<"New POPG_p required for vb row : "<< myrow << " vbRow2BLMap_p[row] " << vbRow2BLMap_p[myrow] << endl;
	    }
	  // else if (blType_p != vbRow2BLMap_p[myrow])
	  //   {
	  //     needsNewPOPG_p=true;
	  //     blType_p = vbRow2BLMap_p[myrow];
	  //     // cerr <<"New POPG_p required for vb row : "<< myrow << " vbRow2BLMap_p[row] " << vbRow2BLMap_p[myrow] << endl;
	  //   }
	  // else
	  //   blNeedsNewPOPG_p[myrow]=false;
	    
	  // if (vectorPhaseGradCalculator_p.nelements() <= (unsigned int) vbRow2BLMap_p[row])
	  // auto itrBLMap_p = find(uniqueVbRow2BLMap_p.begin(),uniqueVbRow2BLMap_p.end(), vbRow2BLMap_p[row]);
	  // idx = distance(uniqueVbRow2BLMap_p.begin(), itrBLMap_p);
	  
	  // cerr<<"#######################" <<endl << "unqVBRow: ";
	  // for(int ii=0; ii<uniqueVbRow2BLMap_p.size();idx++)
	  //   cerr << uniqueVbRow2BLMap_p[ii] << " " ;
	  // cerr<<"#######################"<<endl;
	  // if (vectorPhaseGradCalculator_p.nelements() <= maxVB2BLMap + 1)
	  // vectorPhaseGradCalculator_p.resize(uniqueVbRow2BLMap_p.size()+1,true);


	  vectorPhaseGradCalculator_p.resize(maxVB2BLMap+1,true); // not zero counted hence the plus 1 - PJ
	    {
	       // cerr<<"vbRow2BLMap_p [row] doP=T "<< vbRow2BLMap_p[row] << " " <<row <<endl;
	      // vectorPhaseGradCalculator_p.resize(vbRow2BLMap_p[myrow]+1,true); // Revisit this.
	      for (unsigned int i=0;i<vectorPhaseGradCalculator_p.nelements(); i++)
		if (vectorPhaseGradCalculator_p[i].null())
		  {
		    vectorPhaseGradCalculator_p[i]=new PhaseGrad();
		    blNeedsNewPOPG_p[myrow] = true;
		  }
	    }


	  // if( baselineType_p->cachedGroups_p)
	  //   vectorPhaseGradCalculator_p[vbRow2BLMap_p[myrow]]->ComputeFieldPointingGrad(pointingOffsets_p,cfb,vb,0);
	  // else
	    // vectorPhaseGradCalculator_p[vbRow2BLMap_p[myrow]]->ComputeFieldPointingGrad(pointingOffsets_p,cfb,vb,row);
	  // vectorPhaseGradCalculator_p[vbRow2BLMap_p[myrow]]->needsNewPOPG_p = blNeedsNewPOPG_p[myrow];
	    // needsNewPOPG_p = true;
	    // vectorPhaseGradCalculator_p[idx]->needsNewPOPG_p = needsNewPOPG_p;
	    // vectorPhaseGradCalculator_p[idx]->needsNewFieldPG_p = needsNewFieldPG_p;
	    // vectorPhaseGradCalculator_p[idx]->ComputeFieldPointingGrad(pointingOffsets_p,cfb,vb,myrow);
	    // cerr << "maxVB2BLMap "<< maxVB2BLMap << " vbRow2BLMap_p[myrow]" << vbRow2BLMap_p[myrow] << " "  << myrow << " " << needsNewPOPG_p << endl;

	    vectorPhaseGradCalculator_p[vbRow2BLMap_p[myrow]]->needsNewPOPG_p = blNeedsNewPOPG_p[myrow];
	    vectorPhaseGradCalculator_p[vbRow2BLMap_p[myrow]]->needsNewFieldPG_p = needsNewFieldPG_p;
	    vectorPhaseGradCalculator_p[vbRow2BLMap_p[myrow]]->ComputeFieldPointingGrad(pointingOffsets_p,cfb,vb,myrow);

	  // vectorPhaseGradCalculator_p[vbRow2BLMap_p[myrow]]->needsNewPOPG_p = true;
	 
	}
      else
	{
	  needsNewPOPG_p = false;
	  myrow=0;
	  vbRow2BLMap_p.resize(vb.nRows(),0);
	  blNeedsNewPOPG_p.resize(vb.nRows(),false);
	  // vbRow2BLMap_p[0]=0;
	  if (vectorPhaseGradCalculator_p.nelements() <= (unsigned int) vbRow2BLMap_p[myrow])
	    {
	      //	    cerr<<"vbRow2BLMap_p [row] doP=F "<< vbRow2BLMap_p[myrow] << " " <<myrow <<endl;
	      vectorPhaseGradCalculator_p.resize(1);
	      if (vectorPhaseGradCalculator_p[myrow].null())
		vectorPhaseGradCalculator_p[myrow]=new PhaseGrad();
	    }
	  // vectorPhaseGradCalculator_p[vbRow2BLMap_p[myrow]]->ComputeFieldPointingGrad(pointingOffsets_p,cfb,vb,0);    
	  vectorPhaseGradCalculator_p[vbRow2BLMap_p[myrow]]->needsNewPOPG_p = needsNewPOPG_p;
	  vectorPhaseGradCalculator_p[vbRow2BLMap_p[myrow]]->needsNewFieldPG_p = needsNewFieldPG_p;
	  vectorPhaseGradCalculator_p[vbRow2BLMap_p[myrow]]->ComputeFieldPointingGrad(pointingOffsets_p,cfb,vb,myrow);
	  idx =0;
	}
      // if (needsNewPOPG_p)
	{

	  // vectorPhaseGradCalculator_p[vbRow2BLMap_p[myrow]]->needsNewFieldPG_p = needsNewFieldPG_p;
	  // vectorPhaseGradCalculator_p[vbRow2BLMap_p[myrow]]->ComputeFieldPointingGrad(pointingOffsets_p,cfb,vb,myrow);
	  needsNewPOPG_p = false;
	  needsNewFieldPG_p = false;
	}

      return  vectorPhaseGradCalculator_p[idx]->field_phaseGrad_p;
    
    };


  }
}
