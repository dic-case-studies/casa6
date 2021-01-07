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
    
    // Where do the output SPWs indexes start. This assumes that
    // this TVI is not doing any reindexing in the output it produces.
    outSpwStartIdx_p = inputVii_p->nSpectralWindows();

    // TODO: If the input has not been "reindexed" there can be 
    // spwIds, polIds, ddIds which are not referenced in the main table
    // The proper way to do it is to filter all the input SPWs to those
    // ones for which there is real data.
    auto nPolIds = inputVii_p->nPolarizationIds();
    //There are as many output SPWs as polarizations.
std::cout << " npol " << nPolIds << std::endl;
    for(ssize_t polId = 0; polId < nPolIds; polId++)
    {
        int outSpwId = polId + outSpwStartIdx_p;
        std::vector<double> thisOutSpwFreqs;
        std::vector<int> thisOutSpwChann;
        spwInpChanOutFreqMap_p[outSpwId].clear();
        // TODO: Filter SPWs that do not appear together with this polarization
        // (they are not together in the same DDId)
        for(auto inpSpw : spwInpChanIdxMap_p)
        {
            for(auto channel : inpSpw.second)
            {
                double channFreq = getChannelNominalFreq(inpSpw.first, channel);
                spwInpChanOutFreqMap_p[outSpwId][inpSpw.first].push_back(channFreq);
                thisOutSpwFreqs.push_back(channFreq);
            }
        }
        spwOutFirstFreq_p[outSpwId] = *std::min_element(thisOutSpwFreqs.begin(), thisOutSpwFreqs.end());
        for(auto inpSpw : spwInpChanIdxMap_p)
        {
            for(auto freq : spwInpChanOutFreqMap_p[outSpwId][inpSpw.first])
            {
                //TODO: Check all the frequencies fall more or less in the same bins. With a bin tolerance (0.01% of the bin?)
                //TODO: Check that no frequencies overlap.
                spwInpChanOutMap_p[outSpwId][inpSpw.first].push_back(std::floor((freq - spwOutFirstFreq_p[outSpwId]) / freqWidthChan_p + 0.1));
                thisOutSpwChann.push_back(std::floor((freq - spwOutFirstFreq_p[outSpwId]) / freqWidthChan_p + 0.1));
            }
        }
        spwOutChanNumMap_p[outSpwId] = *std::max_element(thisOutSpwChann.begin(), thisOutSpwChann.end()) + 1;
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

casacore::Vector<double> SPWCombinationTVI::getFrequencies(double time, int frameOfReference,
                                                           int spectralWindowId, int msId) const
{
    if(spectralWindowId < outSpwStartIdx_p ||
        spectralWindowId > (outSpwStartIdx_p + spwInpChanOutFreqMap_p.size()))
        throw casacore::AipsError("SPWId out of valid range");

    auto inputSpwForThisOutputSpw = spwInpChanOutMap_p.at(spectralWindowId);
    std::vector<double> freqs(spwOutChanNumMap_p.at(spectralWindowId));
    for(auto inputSpw : inputSpwForThisOutputSpw)
    {
        auto innerFreqs = getVii()->getFrequencies(time, frameOfReference, inputSpw.first, msId);
        auto outChanThisInputSpw = spwInpChanOutMap_p.at(spectralWindowId).at(inputSpw.first)[0];
        std::copy(innerFreqs.begin(), innerFreqs.end(), freqs.begin() + outChanThisInputSpw);
    }
    return freqs;
}

void SPWCombinationTVI::origin()
{
    // Drive underlying ViImplementation2
    getVii()->origin();

    // Set structure parameters for this subchunk iteration
    setUpCurrentSubchunkStructure();

    // Synchronize own VisBuffer
    configureNewSubchunk();

    return;
}

void SPWCombinationTVI::next()
{
    // Drive underlying ViImplementation2
    getVii()->next();

    // Set structure parameters for this subchunk iteration
    setUpCurrentSubchunkStructure();

    // Synchronize own VisBuffer
    configureNewSubchunk();

    return;
}

void SPWCombinationTVI::setUpCurrentSubchunkStructure()
{
    auto& innerNRowsPerShape = getVii()->nRowsPerShape();
    getVii()->polarizationIds(currentSubchunkInnerPolIds_p);
    std::set<casacore::Int> uniquePolIDs;
    std::copy(currentSubchunkInnerPolIds_p.begin(),currentSubchunkInnerPolIds_p.end(),
              std::inserter(uniquePolIDs, uniquePolIDs.end()));

    ssize_t thisSubchunkNPolIds = uniquePolIDs.size();
    // This VisBuffer contains one single timestamp with "all" DDIds.
    // After SPW combination then number of rows is equal to the number of distinct DDIds,
    // i.e., the number of polarizations, which is also the number of distinct shapes
    nRowsPerShape_p.resize(thisSubchunkNPolIds);
    nChannelsPerShape_p.resize(thisSubchunkNPolIds);


    currentSubchunkSpwIds_p.resize(thisSubchunkNPolIds);
    std::iota(currentSubchunkSpwIds_p.begin(), currentSubchunkSpwIds_p.end(), outSpwStartIdx_p);

    getVii()->spectralWindows(currentSubchunkInnerSpwIds_p);

    // Set up the channels for shape
    size_t iShape=0;
    for(auto outSpw : spwOutChanNumMap_p)
    {
        nChannelsPerShape_p[iShape] = outSpw.second;
        iShape++;
    }

    nRowsPerShape_p = thisSubchunkNPolIds;

    //TODO: Check all the NRows are the same.
    //TODO: For several polarizations (nShapes)
}

casacore::rownr_t SPWCombinationTVI::nShapes() const
{
    return nRowsPerShape_p.size();
}

const casacore::Vector<casacore::rownr_t>& SPWCombinationTVI::nRowsPerShape() const
{
    return nRowsPerShape_p;
}

const casacore::Vector<casacore::Int>& SPWCombinationTVI::nChannelsPerShape() const
{
    return nChannelsPerShape_p;
}

void SPWCombinationTVI::spectralWindows(casacore::Vector<casacore::Int>& spws) const
{
    spws = currentSubchunkSpwIds_p;
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
    auto& visCubes = getVisBuffer()->visCubes();
    vis = visCubes[0];
}

void SPWCombinationTVI::visibilityObserved(casacore::Vector<casacore::Cube<casacore::Complex>> & vis) const
{
    // Get input VisBuffer and visibility observed cubes
    VisBuffer2 *vb = getVii()->getVisBuffer();
    auto& innerVisCubes = vb->visCubes();

    // Resize vis vector
    vis.resize(nShapes());

    size_t iShape = 0;
    // It is assumed that the VisBuffer contains unique metadata and timestamp except for DDId.
    // See checkSortingInner().
    for(auto& visCube : vis)
    {
        casacore::IPosition cubeShape(3, nCorrelationsPerShape()[iShape],
                            nChannelsPerShape_p[iShape], nRowsPerShape_p[iShape]);
        visCube.resize(cubeShape);
        ++iShape;
    }
    casacore::rownr_t inputRowsProcessed = 0;
    for(auto& inputCube : innerVisCubes)
    {
        // It is assumed that each input cube corresponds to 
        // an unique DDiD.
        // The case in which several DDiDs which have the same 
        // number of channels and polarizations have been merged in a single
        // cube with equal shape is not supported yet.
        // By construction of the rest of VB2 (VisibilityIteratorImpl2,
        // SimpleSimVI2, rest of TVIs) this doesn't happen yet anyway.
        casacore::Int thisCubePolId = currentSubchunkInnerPolIds_p[inputRowsProcessed];
        auto thisSpw = currentSubchunkInnerSpwIds_p[inputRowsProcessed];
        auto thisOutputSpw = thisCubePolId + outSpwStartIdx_p;

        auto outputChannel = spwInpChanOutMap_p.at(thisOutputSpw).at(thisSpw)[0];
        casacore::IPosition blcOutput(3, 0, outputChannel, 0);
        casacore::IPosition trcOutput(3, inputCube.shape()(0)-1, outputChannel + inputCube.shape()(1)-1, inputCube.shape()(2)-1);
        vis[thisCubePolId](blcOutput, trcOutput) = inputCube;
        inputRowsProcessed+=inputCube.shape()(2);
    }
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
