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
#include <synthesis/TransformMachines2/PhaseGrad.h>
#include <msvis/MSVis/VisBuffer2.h>
#include <synthesis/TransformMachines2/VB2CFBMap.h>
namespace casa{
using namespace vi;
  namespace refim{
  Int mapAntIDToAntType(const casacore::Int& /*ant*/) {return 0;};

    VB2CFBMap::VB2CFBMap(): vb2CFBMap_p(), cfPhaseGrad_p(), phaseGradCalculator_p(),doPointing_p(false)
    {
      phaseGradCalculator_p = new PhaseGrad();
      newPhaseGradComputed_p = false;
    };

    VB2CFBMap& VB2CFBMap::operator=(const VB2CFBMap& other)
    {
      if(this!=&other) 
	{
	  phaseGradCalculator_p = other.phaseGradCalculator_p;
	  cfPhaseGrad_p.assign(other.cfPhaseGrad_p);
	  vb2CFBMap_p.assign(other.vb2CFBMap_p);
	  doPointing_p = other.doPointing_p;
	}
      return *this;
    };

    void VB2CFBMap::setPhaseGradPerRow(const CountedPtr<PointingOffsets>& pointingOffset,
				       const casacore::CountedPtr<CFBuffer>& cfb,
				       const vi::VisBuffer2& vb,
				       const int& row)
    {
      // if (doPointing_p)
      // 	{
      // 	  if (phaseGradCalculator_p->needsNewPhaseGrad(pointingOffset, vb, 0))
      // 	    {
      // 	      phaseGradCalculator_p->ComputeFieldPointingGrad(pointingOffset,cfb,vb, 0);
      // 	      newPhaseGradComputed_p=true;
      // 	    }
      // 	}
      // else
	{
	  phaseGradCalculator_p->ComputeFieldPointingGrad(pointingOffset,cfb,vb, 0);
	}
	{
	  //cfPhaseGrad_p(row).assign(phaseGradCalculator_p->getFieldPointingGrad());
	  cfPhaseGrad_p(row).reference(phaseGradCalculator_p->field_phaseGrad_p);
	}
    }


     Matrix<vector<int> > VB2CFBMap::findAntennaGroups(const vi::VisBuffer2& vb, 
						      const CountedPtr<PointingOffsets>& pointingOffsets_p, 
						      const double& sigmaDev)
    {
     
      const Int nRows = vb.nRows();
      vector< vector <double> > antDirPix_l = pointingOffsets_p->fetchAntOffsetToPix(vb, doPointing_p);
      MVDirection vbdir=vb.direction1()(0).getValue();	
      Vector<double> phaseDirPix_l = pointingOffsets_p->toPix(vb,vbdir,vbdir);
      int binsx = 5, binsy = 5;
      antennaGroups_p.resize(binsx,binsy);
      for (int ii=0; ii<binsx; ii++)
	for (int jj=0; jj<binsy;jj++)
	  antennaGroups_p(ii,jj).resize(0);
      Vector<double> pixSum, pixMean, pixAbsDev, pixMinDev, pixMaxDev, pixMedAbsDev;
      PO_DEBUG_P = SynthesisUtils::getenv("PO_DEBUG",0);
      // pixSum.resize(4);
      // pixMean.resize(4);
      pixAbsDev.resize(4);
      pixMinDev.resize(4);
      pixMaxDev.resize(4);
      pixMedAbsDev.resize(4);
      /* In order to compute the pointing deviation of the antennas to identify groups,
	 the following is done. We compute the absolute pointing deviation and the median
	 absolute deviation in pointing for the antenna given a visbuffer. once the groups
	 are identified we then put them on a matrix of vector ints containing the antenna
	 index so be used to figure out the baseline groups- PJ */
      for (int irow=0; irow<4; irow++) 
	{
	  pixAbsDev[irow] = 0;
	  pixMinDev[irow] = 0;
	  pixMaxDev[irow] = 0;
	}

      for (int irow=0; irow<nRows; irow++) 
	{
	  for (unsigned int ii=0; ii<phaseDirPix_l.size();ii++)
	    for (unsigned int jj=0; jj<phaseDirPix_l.size();jj++)
	      {
		pixAbsDev[ii+2*jj] += abs(antDirPix_l[ii+2*jj][irow] - phaseDirPix_l[ii]); 

		if(pixMinDev[ii+2*jj] >= abs(antDirPix_l[ii+2*jj][irow] - phaseDirPix_l[ii]))
		  pixMinDev[ii+2*jj]  = abs(antDirPix_l[ii+2*jj][irow] - phaseDirPix_l[ii]);
		else if (pixMaxDev[ii+2*jj] <= abs(antDirPix_l[ii+2*jj][irow] - phaseDirPix_l[jj]))
		  pixMaxDev[ii+2*jj]  = abs(antDirPix_l[ii+2*jj][irow] - phaseDirPix_l[ii]);
	      }
	}
      for (unsigned int ii=0; ii<pixAbsDev.size();ii++)
	{
	  pixAbsDev[ii] = pixAbsDev[ii]/nRows;
	  pixMedAbsDev[ii] = (pixMaxDev[ii]+pixMinDev[ii])/2 ;
	}
      // if (PO_DEBUG_P==1)
      // 	{
      // 	  for (int irow=0; irow<nRows; irow++)
      // 	    {

      // 	      cerr << "pixSigma : " << pixAbsDev << " medABSDev " << pixMedAbsDev << " pixMinDev " << pixMinDev << " pixMaxDev " << pixMaxDev << endl; 
      // 	    }
      // 	}
      
      for (int irow=0; irow<nRows; irow++) 
	{
	  int ii,jj,kk,ll;
	  if(pixMedAbsDev[0] == 0 )
	    ii = (antDirPix_l[0][irow] - phaseDirPix_l[0])/(sigmaDev);
	  else 
	    ii = (antDirPix_l[0][irow] - phaseDirPix_l[0])/(pixMedAbsDev[0]*sigmaDev);
	  if(pixMedAbsDev[1] == 0 )
	    jj = (antDirPix_l[1][irow] - phaseDirPix_l[1])/(sigmaDev);
	  else
	    jj = (antDirPix_l[1][irow] - phaseDirPix_l[1])/(pixMedAbsDev[0]*sigmaDev);
	  // cerr << " ii " << ii << " jj " << jj << " " << antennaGroups_p.shape() << " " << antennaGroups_p(ii,jj).shape()+1 << endl;
	  // if( ii < -1*int(phaseDirPix_l.size()))
	  //   ii = -1*(phaseDirPix_l.size());
	  // else if (ii > int(phaseDirPix_l.size()))
	  //   ii = phaseDirPix_l.size();
	      
	  // if( jj < -1*int(phaseDirPix_l.size()))
	  //   jj = -1*(phaseDirPix_l.size());
	  // else if (jj > int(phaseDirPix_l.size()))
	  //   jj = phaseDirPix_l.size();
		       
	       
	  antennaGroups_p(ii+int(binsx/2),jj+int(binsy/2)).push_back(vb.antenna1()[irow]);

	  if(pixMedAbsDev[2] == 0 )
	    kk = (antDirPix_l[2][irow] - phaseDirPix_l[0])/(sigmaDev);
	  else
	    kk = (antDirPix_l[2][irow] - phaseDirPix_l[0])/(pixMedAbsDev[2]*sigmaDev);

	  if(pixMedAbsDev[3] == 0 )
	    ll = (antDirPix_l[3][irow] - phaseDirPix_l[1])/(sigmaDev);
	  else
	    ll = (antDirPix_l[3][irow] - phaseDirPix_l[1])/(pixMedAbsDev[3]*sigmaDev);

	  // if( kk < -1*int(phaseDirPix_l.size()))
	  //   kk = -1*(phaseDirPix_l.size());
	  // else if (kk > int(phaseDirPix_l.size()))
	  //   kk = phaseDirPix_l.size();
	  // if( ll < -1*int(phaseDirPix_l.size()))
	  //   ll = -1*(phaseDirPix_l.size());
	  // else if (ll > int(phaseDirPix_l.size()))
	  //   ll = phaseDirPix_l.size();
	
	  antennaGroups_p(kk+int(binsx/2),ll+int(binsy/2)).push_back(vb.antenna1()[irow]);
	}

      LogIO log_l(LogOrigin("VB2CFBMap", "VB2CFMap::findAntennaGroups[R&D]"));
      for (int ii=0; ii<2; ii++) // This binning is a bit arbitrary and it looks like its time for some negative indices.
      	for (int jj=0; jj<2;jj++)
	  {
	 
	    std::sort(antennaGroups_p(ii,jj).begin(), antennaGroups_p(ii,jj).end());
	    auto last = std::unique(antennaGroups_p(ii,jj).begin(), antennaGroups_p(ii,jj).end());
	    antennaGroups_p(ii,jj).erase(last, antennaGroups_p(ii,jj).end());
	    for (unsigned int kk=0; kk<antennaGroups_p(ii,jj).size();kk++)
	
	      log_l << "Antenna " << std::to_string(antennaGroups_p(ii,jj)[kk]) << " is present in bin " << std::to_string(ii*2+jj) << LogIO::NORMAL<<LogIO::POST;
	    
	    
	  }

      return antennaGroups_p;



    }

    void VB2CFBMap::makeVBRow2BLGMap(const vi::VisBuffer2& vb) 
				    // const Matrix<vector<int> >& antennaGroups_p)
    {

      /* The goal here is take the antenna groups given a vb and retirn a baseline group
	 map for every VB row. This will then allow for the computation of the phase grad
	 once per phase grad per vb. */

      const Int nRows=vb.nRows(); 
      const Int nx = antennaGroups_p.ncolumn(), ny = antennaGroups_p.nrow();
      Int numAntGroups=0;
      vector<int> mapAntType1, mapAntType2;
      const Vector<Int> antenna1 = vb.antenna1(), antenna2 = vb.antenna2();

      mapAntGrp_p.resize(nx,ny);
      mapAntType1.resize(nRows);
      mapAntType2.resize(nRows);
      vbRow2BLMap_p.resize(nRows);

      for (int ii=0; ii < nx; ii++)
	for (int jj=0; jj < ny; jj++)      
	  {
	    mapAntGrp_p(ii,jj)=0;
	    if(antennaGroups_p(ii,jj).size() > 0)
	      {
		numAntGroups++;
		mapAntGrp_p(ii,jj) = numAntGroups;
	      }
	  }

      cout<<"mapAntGrp "<<mapAntGrp_p<<endl;

      mapBLGroup_p.resize(numAntGroups,numAntGroups);
      for (int ii=0; ii < numAntGroups; ii++)
	for (int jj=0; jj < numAntGroups; jj++)
	  mapBLGroup_p(ii,jj) = ii*numAntGroups+jj;

      cout<<"mapBLGrp "<<mapBLGroup_p<<endl;

      
      LogIO log_l(LogOrigin("VB2CFBMap", "makeVBRow2BLGMap[R&D]"));
      log_l << "Number of Baseline Groups found " << (numAntGroups + (numAntGroups*(numAntGroups-1))/2) << LogIO::POST;

      for (int kk=0; kk<nRows;kk++)
	{
	  for (int ii=0; ii < nx; ii++)
	    for (int jj=0; jj < ny; jj++)      
	      {
		int n=antennaGroups_p(ii,jj).size();
		if( n> 0)
		  for(int tt=0; tt < n; tt++)
		    {
		      if(antenna1[kk] == (antennaGroups_p(ii,jj))[tt])
			mapAntType1[kk] = mapAntGrp_p(ii,jj);
		      // else
		      // 	mapAntType2[kk] = -1;
			
		      if(antenna2[kk] == (antennaGroups_p(ii,jj))[tt])
			mapAntType1[kk] = mapAntGrp_p(ii,jj);
		      // else
		      // 	mapAntType2[kk] = -1;
		    }
	      }   

	}

      // cout << "Baseline Groups for vb Row are as follows : ";
      // for(int ii = 0; ii < mapAntType1.size(); ii++)
      // 	cout<<mapAntType1[ii]<<"  "<<mapAntType2[ii]<<" " <<ii <<endl;
      // cout<<endl<<"#############################################"<<endl;

     
      
      for (int kk=0; kk<nRows;kk++)
	if ((mapAntType1[kk]>=0) && (mapAntType2[kk]>=0)) vbRow2BLMap_p[kk] = mapBLGroup_p(mapAntType1[kk],mapAntType2[kk]);
	else vbRow2BLMap_p[kk] = -1;

      // for (int kk=0; kk<nRows();kk++)
      // 	for (int ii=0; ii < numAntGroups; ii++)
      // 	  for (int jj=ii; jj < numAntGroups; jj++)      
      // 	    {
      // 	      mapBLType.push_back(mapBLGroup(mapAntType1(kk),jj))
	      
      
      // return vbRow2BLMap_p;
      

    }


    
    Int VB2CFBMap::makeVBRow2CFBMap(CFStore2& cfs,
				    const VisBuffer2& vb, 
				    const Quantity& dPA,
				    const Vector<Int>& /*dataChan2ImChanMap*/,
				    const Vector<Int>& /*dataPol2ImPolMap*/,
				    const CountedPtr<PointingOffsets>& po_p)
				    //const Bool& /*doPointing*/)
    {
      //    VBRow2CFMapType& vbRow2CFMap_p,
      const Int nRow=vb.nRows(); 
      //UNUSED: nChan=dataChan2ImChanMap.nelements(), 
      //UNUSED: nPol=dataPol2ImPolMap.nelements();
      //    vbRow2CFMap_p.resize(nPol, nChan, nRow);
      vb2CFBMap_p.resize(nRow);
      cfPhaseGrad_p.resize(nRow);

      Quantity pa(getPA(vb),"rad");
      //PolOuterProduct outerProduct;
      Int statusCode=CFDefs::MEMCACHE;

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
	      setPhaseGradPerRow(po_p, cfb_l, vb, irow);

	      // Set the CFB per VB row
	      cfb_l->setPointingOffset(po_p->pullPointingOffsets());
	      vb2CFBMap_p(irow) = cfb_l;
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
}
}
