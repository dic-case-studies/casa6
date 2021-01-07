//# SPWCombinationTVI.h: This file contains the interface definition of the MSTransformManager class.
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

#ifndef SPWCombinationTVI_H_
#define SPWCombinationTVI_H_

// Base class
#include <mstransform/TVI/FreqAxisTVI.h>
#include <msvis/MSVis/VisibilityIterator2.h>


namespace casa { //# NAMESPACE CASA - BEGIN

namespace vi { //# NAMESPACE VI - BEGIN

//////////////////////////////////////////////////////////////////////////
// SPWCombinationTVI class
//////////////////////////////////////////////////////////////////////////

class SPWCombinationTVI : public FreqAxisTVI
{

public:

    // Constructor
    SPWCombinationTVI(ViImplementation2 * inputVii);

    // Report the the ViImplementation type
    virtual casacore::String ViiType() const override { return casacore::String("SPWCombination( ")+getVii()->ViiType()+" )"; };

    void origin() override;

    void next() override;

    // Structural changes
    casacore::rownr_t nShapes() const override;
    const casacore::Vector<casacore::rownr_t>& nRowsPerShape() const override;
    const casacore::Vector<casacore::Int>& nChannelsPerShape() const override;
    void spectralWindows(casacore::Vector<casacore::Int>& spws) const override;

    // Transformation of the flags
    void flag(casacore::Cube<casacore::Bool>& flagCube) const override;
    void flag(casacore::Vector<casacore::Cube<casacore::Bool>>& flagCubes) const override;

    // Transformation of float data
    void floatData(casacore::Cube<casacore::Float> & vis) const override;
    void floatData(casacore::Vector<casacore::Cube<casacore::Float>> & vis) const override;

    // Transformation of observed visibilities
    void visibilityObserved(casacore::Cube<casacore::Complex> & vis) const override;
    void visibilityObserved(casacore::Vector<casacore::Cube<casacore::Complex>> & vis) const override;

    // Transformation of corrected visibilities
    void visibilityCorrected(casacore::Cube<casacore::Complex> & vis) const override;
    void visibilityCorrected(casacore::Vector<casacore::Cube<casacore::Complex>> & vis) const override;

    // Transformation of model
    void visibilityModel(casacore::Cube<casacore::Complex> & vis) const override;
    void visibilityModel(casacore::Vector<casacore::Cube<casacore::Complex>> & vis) const override;

    // Transformation of weight 
    void weight(casacore::Matrix<casacore::Float> &weight) const override;
    void weight(casacore::Vector<casacore::Matrix<casacore::Float>> &weight) const override;
    void weightSpectrum(casacore::Cube<casacore::Float> &weightSp) const override;
    void weightSpectrum(casacore::Vector<casacore::Cube<casacore::Float>> &weightSp) const override;

    // Transformation of sigma
    void sigma(casacore::Matrix<casacore::Float> &sigma) const override;
    void sigma(casacore::Vector<casacore::Matrix<casacore::Float>> &sigma) const override;
    void sigmaSpectrum(casacore::Cube<casacore::Float> &sigmaSp) const override;
    void sigmaSpectrum(casacore::Vector<casacore::Cube<casacore::Float>> &sigmaSp) const override;

protected:

    void initialize();

    void setUpCurrentSubchunkStructure();

    double checkEqualWidth() const;

    double getChannelNominalFreq(size_t spw, size_t channel) const;

    double freqWidthChan_p;

    std::map<int, double> spwOutFirstFreq_p;

    // Variables that depend on the subchunk we are in
    casacore::Vector<casacore::rownr_t> nRowsPerShape_p;

    casacore::Vector<casacore::Int> nChannelsPerShape_p;

    casacore::Vector<casacore::Int> currentSubchunkSpwIds_p;

    casacore::Vector<casacore::Int> currentSubchunkInnerPolIds_p;

    casacore::Vector<casacore::Int> currentSubchunkInnerSpwIds_p;

    // vector[1..nPolId] map[key:inpSpwId] vector-> inp chann idx, same as spwInpChanIdxMap, value: outputChannel
    std::map<int, std::map<int, std::vector<int>>> spwInpChanOutMap_p;
    // vector[1..nPolId] map[key:inpSpwId] vector-> inp chann idx, same as spwInpChanIdxMap, value: outputChannel freq
    std::map<int, std::map<int, std::vector<double>>> spwInpChanOutFreqMap_p;
    // Start of te output SPW indexes.
    // This assumes that there is no reindexing performed by this TVI, i.e., the output
    // SPW is added to the list of input SPWs. With reindexing this would have been 0.
    casacore::Int outSpwStartIdx_p;
};

/*
 * Factory that allows the creation of SPWCombinationTVI classes.
 * This factory doesn't have any parameter to configure
 */
class SPWCombinationTVIFactory : public ViFactory
{

public:

    SPWCombinationTVIFactory(ViImplementation2 *inputVII);

    ~SPWCombinationTVIFactory();

protected:

    virtual vi::ViImplementation2 * createVi () const;

private:

    ViImplementation2 *inputVii_p;;
};

/*
 * Factory that allows the creation of SPWCombinationTVI classes
 * which act upon an underlying VI2
 * This factory doesn't have any parameter to configure
 */
class SPWCombinationTVILayerFactory : public ViiLayerFactory
{

public:

    SPWCombinationTVILayerFactory();

    virtual ~SPWCombinationTVILayerFactory();

protected:
  
    virtual ViImplementation2 * createInstance(ViImplementation2* vii0) const;
};


} //# NAMESPACE VI - END

} //# NAMESPACE CASA - END

#endif /* SPWCombinationTVI_H_ */
