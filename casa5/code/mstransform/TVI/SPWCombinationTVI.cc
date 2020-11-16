//# SPWCombinationTVI.h: This file contains the implementation of the SPWCombinationTVI class.
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

#include <iomanip>
#include <mstransform/TVI/SPWCombinationTVI.h>
#include <casacore/ms/MeasurementSets/MSSpWindowColumns.h>

namespace casa { //# NAMESPACE CASA - BEGIN

namespace vi { //# NAMESPACE VI - BEGIN

//////////////////////////////////////////////////////////////////////////
// SPWCombinationTVI class
//////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
SPWCombinationTVI::SPWCombinationTVI(ViImplementation2 * inputVii) :
  FreqAxisTVI (inputVii)
{
std::cout << std::setprecision(10) << std::endl;
    initialize();

    return;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void SPWCombinationTVI::initialize()
{
    freqWidthChan_p = checkEqualWidth();

    spwInpChanOutFreqMap_p.clear();

std::cout << " bef npol " << std::endl;
    auto nPolIds = inputVii_p->nPolarizationIds();
    spwInpChanOutFreqMap_p.resize(nPolIds);
    spwInpChanOutMap_p.resize(nPolIds);
    //There are as many output SPWs as polarizations.
std::cout << " npol " << nPolIds << std::endl;
    for(ssize_t polId = 0; polId < nPolIds; polId++)
    {
        int outSpwId = polId;
        std::vector<double> thisOutSpwFreqs;
        std::vector<int> thisOutSpwChann;
        for(auto inpSpw : spwInpChanIdxMap_p)
        {
std::cout << " inpSpw " <<inpSpw.first << std::endl;
            spwInpChanOutFreqMap_p[outSpwId].clear();
            for(auto channel : inpSpw.second)
            {
                double channFreq = getChannelNominalFreq(inpSpw.first, channel);
                spwInpChanOutFreqMap_p[outSpwId][inpSpw.first].push_back(channFreq);
                thisOutSpwFreqs.push_back(channFreq);
std::cout << " chann freq " <<channFreq << std::endl;
            }
        }
        spwOutFirstFreq[outSpwId] = *std::min_element(thisOutSpwFreqs.begin(), thisOutSpwFreqs.end());
std::cout << " first freq " <<spwOutFirstFreq[outSpwId] << std::endl;
        for(auto inpSpw : spwInpChanIdxMap_p)
            for(auto freq : spwInpChanOutFreqMap_p[outSpwId][inpSpw.first])
            {
                spwInpChanOutMap_p[outSpwId][inpSpw.first].push_back(std::floor((freq - spwOutFirstFreq[outSpwId]) / freqWidthChan_p + 0.1));
                thisOutSpwChann.push_back(std::floor((freq - spwOutFirstFreq[outSpwId]) / freqWidthChan_p + 0.1));
            }
        spwOutChanNumMap_p[outSpwId] = *std::max_element(thisOutSpwChann.begin(), thisOutSpwChann.end()) + 1;
std::cout << " spwOutChanNumMap_p[outSpwId]   " << spwOutChanNumMap_p[outSpwId]<< std::endl;
    }
            
    return;
}

double SPWCombinationTVI::getChannelNominalFreq(size_t inpSpw, size_t channel) const
{
    //TODO: refactor with next one
    auto& inputSPWSubtablecols = inputVii_p->spectralWindowSubtablecols();

    auto& channelNominalFreqs = inputSPWSubtablecols.chanFreq();
//std::cout << " channelNominalFreqs " << channelNominalFreqs.getColumn() << std::endl;

    return channelNominalFreqs(inpSpw)(casacore::IPosition(1, channel));
}

double SPWCombinationTVI::checkEqualWidth() const
{
    double channelWidth = -1;
    auto& inputSPWSubtablecols = inputVii_p->spectralWindowSubtablecols();

    auto& channelWidths = inputSPWSubtablecols.chanWidth();

    for(auto& inpSpw : spwInpChanIdxMap_p)
    {
        if (channelWidth == -1)
            channelWidth = channelWidths(inpSpw.first)(casacore::IPosition(1, 0));
        auto spwChannelWidths = channelWidths(inpSpw.first); //This is a copy
        if(!std::all_of(spwChannelWidths.begin(), spwChannelWidths.end(), [&](double width){return width == channelWidth;}))
            throw casacore::AipsError("SPWs to combine have different widths");
    }
    return channelWidth;
}

void SPWCombinationTVI::origin()
{
    // Drive underlying ViImplementation2
    getVii()->origin();

    // Synchronize own VisBuffer
    configureNewSubchunk();

    return;
}

void SPWCombinationTVI::next()
{
    // Drive underlying ViImplementation2
    getVii()->next();

    // Synchronize own VisBuffer
    configureNewSubchunk();

    return;
}

void SPWCombinationTVI::flag(casacore::Cube<casacore::Bool>& flagCube) const
{
    getVii()->flag(flagCube);
}

void SPWCombinationTVI::flag(casacore::Vector<casacore::Cube<casacore::Bool>>& flagCubes) const
{
    getVii()->flag(flagCubes);
}

void SPWCombinationTVI::floatData(casacore::Cube<casacore::Float> & vis) const
{
    getVii()->floatData(vis);
}

void SPWCombinationTVI::floatData(casacore::Vector<casacore::Cube<casacore::Float>> & vis) const
{
    getVii()->floatData(vis);
}

void SPWCombinationTVI::visibilityObserved(casacore::Cube<casacore::Complex> & vis) const
{
    getVii()->visibilityObserved(vis);
}

void SPWCombinationTVI::visibilityObserved(casacore::Vector<casacore::Cube<casacore::Complex>> & vis) const
{
    getVii()->visibilityObserved(vis);
}

void SPWCombinationTVI::visibilityCorrected(casacore::Cube<casacore::Complex> & vis) const
{
    getVii()->visibilityCorrected(vis);
}

void SPWCombinationTVI::visibilityCorrected(casacore::Vector<casacore::Cube<casacore::Complex>> & vis) const
{
    getVii()->visibilityCorrected(vis);
}

void SPWCombinationTVI::visibilityModel(casacore::Cube<casacore::Complex> & vis) const
{
    getVii()->visibilityModel(vis);
}

void SPWCombinationTVI::visibilityModel(casacore::Vector<casacore::Cube<casacore::Complex>> & vis) const
{
    getVii()->visibilityModel(vis);
}

void SPWCombinationTVI::weight(casacore::Matrix<casacore::Float> &weight) const
{
    getVii()->weight(weight);
}

void SPWCombinationTVI::weight(casacore::Vector<casacore::Matrix<casacore::Float>> &weight) const
{
    getVii()->weight(weight);
}

void SPWCombinationTVI::weightSpectrum(casacore::Cube<casacore::Float> &weightSp) const
{
    getVii()->weightSpectrum(weightSp);
}

void SPWCombinationTVI::weightSpectrum(casacore::Vector<casacore::Cube<casacore::Float>> &weightSp) const
{
    getVii()->weightSpectrum(weightSp);
}

void SPWCombinationTVI::sigma(casacore::Matrix<casacore::Float> &sigma) const
{
    getVii()->sigma(sigma);
}

void SPWCombinationTVI::sigma(casacore::Vector<casacore::Matrix<casacore::Float>> &sigma) const
{
    getVii()->sigma(sigma);
}

void SPWCombinationTVI::sigmaSpectrum(casacore::Cube<casacore::Float> &sigmaSp) const
{
    getVii()->sigmaSpectrum(sigmaSp);
}

void SPWCombinationTVI::sigmaSpectrum(casacore::Vector<casacore::Cube<casacore::Float>> &sigmaSp) const
{
    getVii()->sigmaSpectrum(sigmaSp);
}

SPWCombinationTVIFactory::SPWCombinationTVIFactory (ViImplementation2 *inputVii)
 : inputVii_p(inputVii)
{
}

SPWCombinationTVIFactory::~SPWCombinationTVIFactory()
{
}

ViImplementation2 * SPWCombinationTVIFactory::createVi() const
{
    ViImplementation2* vii = new SPWCombinationTVI(inputVii_p);
    return vii;
}

SPWCombinationTVILayerFactory::SPWCombinationTVILayerFactory()
{
}

SPWCombinationTVILayerFactory::~SPWCombinationTVILayerFactory()
{
}

ViImplementation2* 
SPWCombinationTVILayerFactory::createInstance(ViImplementation2* vii0) const 
{
    // Make the SPWCombinationTVI, using supplied ViImplementation2, and return it
    ViImplementation2 *vii = new SPWCombinationTVI(vii0);
    return vii;
}

} //# NAMESPACE VI - END

} //# NAMESPACE CASA - END
