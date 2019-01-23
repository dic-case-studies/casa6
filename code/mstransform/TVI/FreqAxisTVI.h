//# FreqAxisTVI.h: This file contains the interface definition of the MSTransformManager class.
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

#ifndef FreqAxisTVI_H_
#define FreqAxisTVI_H_

// Base class
#include <msvis/MSVis/TransformingVi2.h>

// VI/VB framework
#include <msvis/MSVis/VisBuffer2.h>
#include <msvis/MSVis/VisibilityIterator2.h>

// TVI framework
#include <mstransform/TVI/UtilsTVI.h>

// Measurement Set
#include <casacore/ms/MSSel/MSSelection.h>

// NOTE: See implementation include below


namespace casa { //# NAMESPACE CASA - BEGIN

namespace vi { //# NAMESPACE VI - BEGIN

//////////////////////////////////////////////////////////////////////////
// FreqAxisTVI class
//////////////////////////////////////////////////////////////////////////

template<class T> class FreqAxisTransformEngine; // Forward declaration
template<class T> class FreqAxisTransformEngine2; // Forward declaration

// <summary>
// FreqAxisTVI can be used as a base class for TVIs that operate on the frequency axis
// </summary>
//
// <synopsis>
// This class can be used as a convenience base class for TVIs that
// operate on the frequency axis.
// It provides a generic method to modify generic data cubes with a
// transformation provided by the derived class.
//
// The derived class has to initialize the following variable
// spwOutChanIdxMap_p: Map with the list of channels for each output SPW
//
// </synopsis>
//
class FreqAxisTVI : public TransformingVi2
{

public:

    // Lifecycle
    FreqAxisTVI(ViImplementation2 * inputVii);
    ~FreqAxisTVI();

    // Navigation methods
    virtual void origin ();
    virtual void next ();

    // General TVI info (common for all sub-classes)
    casacore::Bool existsColumn (VisBufferComponent2 id) const;
    casacore::Bool flagCategoryExists () const {return false;}

    // casacore::List of methods that should be implemented by derived classes
    // virtual void flag(casacore::Cube<casacore::Bool>& flagCube) const = 0;
    // virtual void floatData (casacore::Cube<casacore::Float> & vis) const = 0;
    // virtual void visibilityObserved (casacore::Cube<casacore::Complex> & vis) const = 0;
    // virtual void visibilityCorrected (casacore::Cube<casacore::Complex> & vis) const = 0;
    // virtual void visibilityModel (casacore::Cube<casacore::Complex> & vis) const = 0;
    // virtual void weightSpectrum(casacore::Cube<casacore::Float> &weightSp) const = 0;
    // virtual void sigmaSpectrum (casacore::Cube<casacore::Float> &sigmaSp) const = 0;
    // virtual casacore::Vector<casacore::Double> getFrequencies (	casacore::Double time, casacore::Int frameOfReference,casacore::Int spectralWindowId, casacore::Int msId) const = 0;
    // virtual void writeFlag (const casacore::Cube<casacore::Bool> & flagCube) = 0;

    // Common transformation for all sub-classes
    void writeFlagRow (const casacore::Vector<casacore::Bool> & flag);
    casacore::Vector<casacore::Int> getChannels (	casacore::Double time, casacore::Int frameOfReference,
            casacore::Int spectralWindowId, casacore::Int msId) const;
    void flagRow (casacore::Vector<casacore::Bool> & flagRow) const;
    void weight (casacore::Matrix<casacore::Float> & weight) const;
    void sigma (casacore::Matrix<casacore::Float> & sigma) const;

protected:

    // Method implementing main loop  (with auxiliary data)
    template <class T> void transformFreqAxis(	casacore::Cube<T> const &inputDataCube,
            casacore::Cube<T> &outputDataCube,
            FreqAxisTransformEngine<T> &transformer) const
    {
        // Re-shape output data cube
        outputDataCube.resize(getVisBuffer()->getShape(),false);

        // Get data shape for iteration
        const casacore::IPosition &inputShape = inputDataCube.shape();
        casacore::uInt nRows = inputShape(2);
        casacore::uInt nCorrs = inputShape(0);

        // Initialize input-output planes
        casacore::Matrix<T> inputDataPlane;
        casacore::Matrix<T> outputDataPlane;

        // Initialize input-output vectors
        casacore::Vector<T> inputDataVector;
        casacore::Vector<T> outputDataVector;

        for (casacore::uInt row=0; row < nRows; row++)
        {
            // Assign input-output planes by reference
            transformer.setRowIndex(row);
            inputDataPlane.reference(inputDataCube.xyPlane(row));
            outputDataPlane.reference(outputDataCube.xyPlane(row));

            for (casacore::uInt corr=0; corr < nCorrs; corr++)
            {
                // Assign input-output vectors by reference
                transformer.setCorrIndex(corr);
                inputDataVector.reference(inputDataPlane.row(corr));
                outputDataVector.reference(outputDataPlane.row(corr));

                // Transform data
                transformer.transform(inputDataVector,outputDataVector);
            }
        }

        return;
    }

    // Method implementing main loop  (with auxiliary data)
    template <class T> void transformFreqAxis2(	const casacore::IPosition &inputShape,
            FreqAxisTransformEngine2<T> &transformer,
            casacore::Int parallelCorrAxis=-1) const
    {
        casacore::uInt nRows = inputShape(2);
        if (parallelCorrAxis >= 0)
        {
            for (casacore::uInt row=0; row < nRows; row++)
            {
                transformer.setRowIndex(row);
                transformer.setCorrIndex(parallelCorrAxis);
                transformer.transform();
            }
        }
        else
        {
            casacore::uInt nCorrs = inputShape(0);
            for (casacore::uInt row=0; row < nRows; row++)
            {
                transformer.setRowIndex(row);

                for (casacore::uInt corr=0; corr < nCorrs; corr++)
                {
                    transformer.setCorrIndex(corr);

                    // jagonzal: Debug code
                    /*
					VisBuffer2 *vb = getVii()->getVisBuffer();
					if (vb->rowIds()(row)==0 and corr==0)
					{
						transformer.setDebug(True);
					}
					else
					{
						transformer.setDebug(False);
					}
                     */
                    transformer.transform();
                }
            }
        }

        return;
    }

    void initialize();

    // Form spwInpChanIdxMap_p via calls to underlying Vii
    void formChanMap();

    mutable casacore::LogIO logger_p;

    // Map with the list of channel indices for each input SPW
    std::map<casacore::Int,std::vector<casacore::Int> > spwInpChanIdxMap_p;

    // Map with the list of channel indices for each output SPW.
    // Note that this member needs to be initialized by each derived class
    std::map<casacore::Int,std::vector<casacore::Int> > spwOutChanIdxMap_p;
};

//////////////////////////////////////////////////////////////////////////
// FreqAxisTransformEngine class
//////////////////////////////////////////////////////////////////////////

template<class T> class FreqAxisTransformEngine
{

public:

    virtual void transform(	casacore::Vector<T> &,casacore::Vector<T> &) {};
    virtual void setRowIndex(casacore::uInt row) {row_p = row;}
    virtual void setCorrIndex(casacore::uInt corr) {corr_p = corr;}

protected:

    casacore::uInt row_p;
    casacore::uInt corr_p;

};

//////////////////////////////////////////////////////////////////////////
// FreqAxisTransformEngine2 class
//////////////////////////////////////////////////////////////////////////

template<class T> class FreqAxisTransformEngine2
{

public:

    FreqAxisTransformEngine2(DataCubeMap *inputData,DataCubeMap *outputData) :
        debug_p(false), rowIndex_p(0), corrIndex_p(0),
        inputData_p(inputData), outputData_p(outputData)
    {
    }

    void setRowIndex(casacore::uInt row)
    {
        rowIndex_p = row;
        inputData_p->setMatrixIndex(row);
        outputData_p->setMatrixIndex(row);

        return;
    }

    void setCorrIndex(casacore::uInt corr)
    {
        corrIndex_p = corr;
        inputData_p->setVectorIndex(corr);
        outputData_p->setVectorIndex(corr);

        return;
    }

    void setDebug(bool debug) { debug_p = debug;}

    virtual void transform() {}

protected:

    casacore::Bool debug_p;
    casacore::uInt rowIndex_p;
    casacore::uInt corrIndex_p;
    DataCubeMap *inputData_p;
    DataCubeMap *outputData_p;

};

} //# NAMESPACE VI - END

} //# NAMESPACE CASA - END

#endif /* FreqAxisTVI_H_ */

