//# FreqAxisTVI.h: This file contains the implementation of the FreqAxisTVI class.
//#
//#  CASA - Common Astronomy Software Applications (http://casa.nrao.edu/)
//#  Copyright (C) Associated Universities, Inc. Washington DC, USA 2011, All rights reserved.
//#  Copyright (C) European Southern Observatory, 2011, All rights reserved.
//#
//#  This library is free software; you can redistribute it and/or
//#  modify it under the terms of the GNU Lesser General Public
//#  License as published by the Free software Foundation; either
//#  version 2.1 of the License, or (at your option) any later version.
//#
//#  This library is distributed in the hope that it will be useful,
//#  but WITHOUT ANY WARRANTY, without even the implied warranty of
//#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//#  Lesser General Public License for more details.
//#
//#  You should have received a copy of the GNU Lesser General Public
//#  License along with this library; if not, write to the Free Software
//#  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//#  MA 02111-1307  USA
//# $Id: $

#include <mstransform/TVI/FreqAxisTVI.h>

using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN

namespace vi { //# NAMESPACE VI - BEGIN

//////////////////////////////////////////////////////////////////////////
// FreqAxisTVI class
//////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
FreqAxisTVI::FreqAxisTVI(	ViImplementation2 * inputVii,
							const Record &configuration):
							TransformingVi2 (inputVii)
{
	// Parse and check configuration parameters
	// Note: if a constructor finishes by throwing an exception, the memory
	// associated with the object itself is cleaned up — there is no memory leak.
	if (not parseConfiguration(configuration))
	{
		throw AipsError("Error parsing FreqAxisTVI configuration");
	}

	initialize();

	// Initialize attached VisBuffer
	setVisBuffer(createAttachedVisBuffer (VbRekeyable));

	return;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
FreqAxisTVI::~FreqAxisTVI()
{
	// The parent class destructor (~TransformingVi2) deletes the inner
	// ViImplementation2 object. However if it might have been already
	// deleted at the top level context
	// 2/8/2016 (jagonzal): As per request from George M. (via CAS-8220)
	// I allow TransformingVi2 destructor to delete its inner input VI;
	// This relies on the application layer that produces the inner VI not
	// deleting it which can be guaranteed when using the Factory pattern.
	// inputVii_p = NULL;

	return;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
Bool FreqAxisTVI::parseConfiguration(const Record &configuration)
{
	int exists = -1;
	Bool ret = true;

	// Parse spw selection (optional)
	exists = -1;
	exists = configuration.fieldNumber ("spw");
	if (exists >= 0)
	{
		configuration.get (exists, spwSelection_p);
		logger_p << LogIO::DEBUG1 << LogOrigin("FreqAxisTVI", __FUNCTION__)
				<< "spw selection is " << spwSelection_p
				<< LogIO::POST;
	}
	else
	{
		spwSelection_p = "*";
	}

	return ret;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void FreqAxisTVI::initialize()
{

    if (inputVii_p == nullptr)
        return;

    if (inputVii_p->msName()=="<noms>")
        // Handle "no-MS" case  (SimpleSimVi2 as base layer)
        formSelectedChanMap();
    else {

        // Get list of selected SPWs and channels
        MSSelection mssel;
        mssel.setSpwExpr(spwSelection_p);
        Matrix<Int> spwchan = mssel.getChanList(&(inputVii_p->ms()));
        logger_p << LogIO::DEBUG1 << LogOrigin("FreqAxisTVI", __FUNCTION__)
			        << "Selected SPW:Channels are " << spwchan << LogIO::POST;

        // Convert list of selected SPWs/Channels into a map
        spwInpChanIdxMap_p.clear();
        uInt nSelections = spwchan.shape()[0];
        Int spw,channelStart,channelStop,channelStep;
        for(uInt selection_i=0;selection_i<nSelections;selection_i++)
        {
            spw = spwchan(selection_i,0);
            channelStart = spwchan(selection_i,1);
            channelStop = spwchan(selection_i,2);
            channelStep = spwchan(selection_i,3);

            if (spwInpChanIdxMap_p.find(spw) == spwInpChanIdxMap_p.end())
            {
                spwInpChanIdxMap_p[spw].clear(); // Accessing the vector creates it
            }

            for (Int inpChan=channelStart;inpChan<=channelStop;inpChan += channelStep)
            {
                spwInpChanIdxMap_p[spw].push_back(inpChan);
            }
        }
    }

    return;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void FreqAxisTVI::formSelectedChanMap()
{
	// This triggers realization of the channel selection
	inputVii_p->originChunks();

	// Refresh map
	spwInpChanIdxMap_p.clear();

	for (Int ispw = 0; ispw < inputVii_p->nSpectralWindows(); ++ispw)
	{

		// TBD trap unselected spws with a continue

		Vector<Int> chansV;
		chansV.reference(inputVii_p->getChannels(0.0, -1, ispw, 0));

		Int nChan = chansV.nelements();
		if (nChan > 0)
		{
			spwInpChanIdxMap_p[ispw].clear(); // creates ispw's map
			for (Int ich = 0; ich < nChan; ++ich)
			{
				spwInpChanIdxMap_p[ispw].push_back(chansV[ich]); // accum into map
			}
		}
	} // ispw

	return;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void FreqAxisTVI::origin()
{
	// Drive underlying ViImplementation2
	getVii()->origin();

	// Synchronize own VisBuffer
	configureNewSubchunk();

	return;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void FreqAxisTVI::next()
{
	// Drive underlying ViImplementation2
	getVii()->next();

	// Synchronize own VisBuffer
	configureNewSubchunk();

	return;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
Bool FreqAxisTVI::existsColumn (VisBufferComponent2 id) const
{

	Bool ret;
	switch (id)
	{
		case VisBufferComponent2::WeightSpectrum:
		{
			ret = true;
			break;
		}
		case VisBufferComponent2::SigmaSpectrum:
		{
			ret = true;
			break;
		}
		default:
		{
			ret = getVii()->existsColumn(id);
			break;
		}
	}

	return ret;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
Vector<Int> FreqAxisTVI::getChannels (Double,Int,Int spectralWindowId,Int) const
{
	Vector<Int> ret(spwOutChanNumMap_p[spectralWindowId]);

	for (uInt chanIdx = 0; chanIdx<spwOutChanNumMap_p[spectralWindowId];chanIdx++)
	{
		ret(chanIdx) = chanIdx;
	}

	return ret;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void FreqAxisTVI::writeFlagRow (const Vector<Bool> & flag)
{
	getVii()->writeFlagRow(flag);
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void FreqAxisTVI::flagRow (Vector<Bool> & flagRow) const
{
	// Get flagCube from own VisBuffer
	const Cube<Bool> &flagCube = getVisBuffer()->flagCube();

	// Calculate output flagRow
	accumulateFlagCube(flagCube,flagRow);
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void FreqAxisTVI::weight (Matrix<Float> & weight) const
{
	if (weightSpectrumExists()) // Defined by each derived class or inner TVI
	{
		// Get flags and weightSpectrum from own VisBuffer
		const Cube<Bool> &flags = getVisBuffer()->flagCube();
		const Cube<Float> &weightSpectrum = getVisBuffer()->weightSpectrum();

		// Calculate output weight
		accumulateWeightCube(weightSpectrum,flags,weight);
	}
	else
	{
		getVii()->weight (weight);
	}

	return;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void FreqAxisTVI::sigma (Matrix<Float> & sigma) const
{
	if (sigmaSpectrumExists())
	{
		// Get flags and sigmaSpectrum from own VisBuffer
		const Cube<Bool> &flags = getVisBuffer()->flagCube();
		const Cube<Float> &sigmaSpectrum = getVisBuffer()->sigmaSpectrum();

		// Calculate output sigma
		accumulateWeightCube(sigmaSpectrum,flags,sigma);
	}
	else
	{
		getVii()->sigma (sigma);
	}

	return;
}

} //# NAMESPACE VI - END

} //# NAMESPACE CASA - END


