// -*- C++ -*-
//# BaselineType.cc: Implementation of the PhaseGrad class
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


/* 

The intent of BaselineType Object is two fold

1. Based on the type of the telescope and antenna define an enumerated baseline Type
2. This is also the location where PhaseGrad for the Pointing Table for a given BlType
calculated.

The first point is still an action time to finish. This will allow for CF production
for a truly heterogenous array, such as ALMA.
The second however is what the rest of the code here is accomplishing. We compute
the number of antenna groups that have a clustered pointing value based on nsigma
binning criterion from the user. We then parse the pointing Table and determine the
variance of the antennnas after subtracting the mean pointing direction. We then bin
the antennas based on their deviation from the mean into one of the bins. This defines
our antenna groups. We will then utilize nG + nG C 2 number of phase gradients during
imaging with doPointing=True


 */

#include <synthesis/TransformMachines2/BaselineType.h>
#include <casa/Logging/LogIO.h>
#include <casa/Logging/LogSink.h>
#include <casa/Logging/LogOrigin.h>
#include <scimath/Mathematics/Combinatorics.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/ArrayMath.h>

using namespace casacore;
namespace casa{
  using namespace refim;
  BaselineType::~BaselineType() 
  {
    vectorPhaseGradCalculator_p.resize(0);
  };
  
  // -----------------------------------------------------------------
  BaselineType::BaselineType(): phaseGradCalculator_p(), doPointing_p(true)
  {
    phaseGradCalculator_p = new PhaseGrad();
    newPhaseGradComputed_p = false;
  };

  BaselineType& BaselineType::operator=(const BaselineType& other)
  {
    if(this!=&other)
      {
        phaseGradCalculator_p = other.phaseGradCalculator_p;
	vectorPhaseGradCalculator_p = other.vectorPhaseGradCalculator_p;

      }
    return *this;
  
  };
  
  Matrix<Complex> BaselineType::setBLPhaseGrad(const CountedPtr<PointingOffsets>& pointingOffsets_p,
				       const CountedPtr<CFBuffer>& cfb,
				       const vi::VisBuffer2& vb,
				       const int& row)
  {
    int myrow=row;
    if(doPointing_p)
      {
	makeVBRow2BLGMap(vb);
      
	if (vectorPhaseGradCalculator_p.nelements() <= (unsigned int) vbRow2BLMap_p[row])
	  {
	    //  cerr<<"vbRow2BLMap_p [row] doP=T "<< vbRow2BLMap_p[row] << " " <<row <<endl;
	    vectorPhaseGradCalculator_p.resize(vbRow2BLMap_p[row]+1,true);
	    for (unsigned int i=0;i<vectorPhaseGradCalculator_p.nelements(); i++)
	      if (vectorPhaseGradCalculator_p[i].null())
		vectorPhaseGradCalculator_p[i]=new PhaseGrad();
	  }


      }
    else
      {
	myrow=0;
	vbRow2BLMap_p.resize(1);
	vbRow2BLMap_p[0]=0;
	if (vectorPhaseGradCalculator_p.nelements() <= (unsigned int) vbRow2BLMap_p[myrow])
	  {
	    //	    cerr<<"vbRow2BLMap_p [row] doP=F "<< vbRow2BLMap_p[myrow] << " " <<myrow <<endl;
	    vectorPhaseGradCalculator_p.resize(1);
	    if (vectorPhaseGradCalculator_p[myrow].null())
	      vectorPhaseGradCalculator_p[myrow]=new PhaseGrad();
	  }
      }
	
  
    vectorPhaseGradCalculator_p[vbRow2BLMap_p[myrow]]->ComputeFieldPointingGrad(pointingOffsets_p,cfb,vb,0);    
    
     return  vectorPhaseGradCalculator_p[vbRow2BLMap_p[myrow]]->field_phaseGrad_p;
    
  };


     Matrix<vector<int> > BaselineType::findAntennaGroups(const vi::VisBuffer2& vb, 
						      const CountedPtr<PointingOffsets>& pointingOffsets_p, 
						      const double& sigmaDev)
    {
      if (doPointing_p==false) return antennaGroups_p;
     
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
	  // if (PO_DEBUG_P==1)
	  //   {
	  //     cerr << " ii " << ii << " jj " << jj << " " << antennaGroups_p.shape() << " " << antennaGroups_p(ii,jj).shape()+1 << endl;
	  //     if( ii < -1*int(phaseDirPix_l.size()))
	  // 	ii = -1*(phaseDirPix_l.size());
	  //     else if (ii > int(phaseDirPix_l.size()))
	  // 	ii = phaseDirPix_l.size();
	      
	  //     if( jj < -1*int(phaseDirPix_l.size()))
	  // 	jj = -1*(phaseDirPix_l.size());
	  //     else if (jj > int(phaseDirPix_l.size()))
	  // 	jj = phaseDirPix_l.size();
	  //   }
		       
	       
	  antennaGroups_p(ii+int(binsx/2),jj+int(binsy/2)).push_back(vb.antenna1()[irow]);

	  if(pixMedAbsDev[2] == 0 )
	    kk = (antDirPix_l[2][irow] - phaseDirPix_l[0])/(sigmaDev);
	  else
	    kk = (antDirPix_l[2][irow] - phaseDirPix_l[0])/(pixMedAbsDev[2]*sigmaDev);

	  if(pixMedAbsDev[3] == 0 )
	    ll = (antDirPix_l[3][irow] - phaseDirPix_l[1])/(sigmaDev);
	  else
	    ll = (antDirPix_l[3][irow] - phaseDirPix_l[1])/(pixMedAbsDev[3]*sigmaDev);

	  // if (PO_DEBUG_P==1)
	  //   {
	  //     if( kk < -1*int(phaseDirPix_l.size()))
	  // 	kk = -1*(phaseDirPix_l.size());
	  //     else if (kk > int(phaseDirPix_l.size()))
	  // 	kk = phaseDirPix_l.size();
	  //     if( ll < -1*int(phaseDirPix_l.size()))
	  // 	ll = -1*(phaseDirPix_l.size());
	  //     else if (ll > int(phaseDirPix_l.size()))
	  // 	ll = phaseDirPix_l.size();
	  //   }

	  antennaGroups_p(kk+int(binsx/2),ll+int(binsy/2)).push_back(vb.antenna2()[irow]);
	}

      //LogIO log_l(LogOrigin("VB2CFBMap", "VB2CFMap::findAntennaGroups[R&D]"));
      for (int ii=0; ii<2; ii++) // This binning is a bit arbitrary and it looks like its time for some negative indices.
      	for (int jj=0; jj<2;jj++)
	  {
	 
	    std::sort(antennaGroups_p(ii,jj).begin(), antennaGroups_p(ii,jj).end());
	    auto last = std::unique(antennaGroups_p(ii,jj).begin(), antennaGroups_p(ii,jj).end());
	    antennaGroups_p(ii,jj).erase(last, antennaGroups_p(ii,jj).end());
	    // for (unsigned int kk=0; kk<antennaGroups_p(ii,jj).size();kk++)
	
	      // log_l << "Antenna " << std::to_string(antennaGroups_p(ii,jj)[kk]) << " is present in bin " << std::to_string(ii*2+jj) << LogIO::NORMAL<<LogIO::POST;
	    
	    
	  }

      return antennaGroups_p;



    }


  void BaselineType::makeVBRow2BLGMap(const vi::VisBuffer2& vb) 
  // const Matrix<vector<int> >& antennaGroups_p)
  {

    /* The goal here is take the antenna groups given a vb and return a baseline group
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

    // cout<<"mapAntGrp "<<mapAntGrp_p<<endl;

    mapBLGroup_p.resize(numAntGroups,numAntGroups);
    for (int ii=0; ii < numAntGroups; ii++)
      for (int jj=0; jj < numAntGroups; jj++)
	mapBLGroup_p(ii,jj) = ii*numAntGroups+jj;

    // cout<<"mapBLGrp "<<mapBLGroup_p<<endl;

      
    // LogIO log_l(LogOrigin("VB2CFBMap", "makeVBRow2BLGMap[R&D]"));
    // log_l << "Number of Baseline Groups found " << (numAntGroups + (numAntGroups*(numAntGroups-1))/2) << LogIO::POST;

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
  


  int BaselineType::numPhaseGrads(const int& nG)
  {
    if(nG<2)
      return nG;
    else
      return nG + Combinatorics::choose(nG,2);


  }

};
