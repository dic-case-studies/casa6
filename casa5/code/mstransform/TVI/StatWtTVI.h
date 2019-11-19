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

#ifndef STATWTTVI_H_
#define STATWTTVI_H_

#include <msvis/MSVis/TransformingVi2.h>

#include <casacore/casa/aipstype.h>
#include <casacore/casa/Arrays/ArrayIter.h>
#include <casacore/casa/BasicSL/Complex.h>
#include <casacore/ms/MSSel/MSSelection.h>
#include <casacore/scimath/StatsFramework/StatisticsAlgorithmFactory.h>

#include <msvis/MSVis/VisBuffer2.h>
#include <msvis/MSVis/VisibilityIterator2.h>
#include <mstransform/TVI/StatWtClassicalDataAggregator.h>
#include <mstransform/TVI/StatWtTypes.h>
#include <mstransform/TVI/UtilsTVI.h>
#include <stdcasa/variant.h>
#include <stdcasa/StdCasa/CasacSupport.h>

#include <map>
#include <vector>
#include <utility>
#include <set>
#include <memory>

namespace casa {

class StatWtVarianceAndWeightCalculator;

namespace vi {

class StatWtTVI : public TransformingVi2 {

public:

    static const casacore::String CHANBIN;

    using Baseline = std::pair<casacore::uInt, casacore::uInt>;

    /*
    struct ChanBin {
        casacore::uInt start = 0;
        casacore::uInt end = 0;

        bool operator<(const ChanBin& other) const {
            if (start < other.start) {
                return true;
            }
            if (start == other.start && end < other.end) {
                return true;
            }
            return false;
        }
    };
    */


    // The following fields are supported in the input configuration record
    // combine           String, if contains "corr", data will be aggregated
    //                   across correlations.
    // value in CHANBIN: Int or Quantity String, describes channel bin widths
    //                   in which to aggregate data within spectral windows
    //                   (spw boundaries are not crossed). If not supplied,
    //                   data for all channels in each spectral window are
    //                   aggregated.
    // minsamp:          Int, minimum number of samples required in an
    //                   aggregated set, if less than that, stats are not
    //                   computed and the data in the sample are flagged. If
    //                   not supplied, 2 is used.
    // statalg           String representing what statistics algorithm to use.
    //                   "cl", "ch", "f", "h".
    // maxiter           Int max number of iterations for Chauvenet algorithm
    // zscore            Double zscore for Chauvenet algorithm
    // center            String center for FitToHalf algorithm, "mean",
    //                   "median", or "zero"
    // lside             Bool side to use for FitToHalf algorithm, True means
    //                   <= center side.
    // fence             Double fence value for HingesFences algorithm
    // wtrange           Zero or two element Array<Double>. Specifies the range
    //                   of "good"
    //                   weight values. Data with weights computed to be outside
    //                   this range will be flagged. Both elements must be
    //                   non-negative. If zero length, all weights are
    //                   acceptable.
    // fitspw            String. MSSelection string representing channels to
    //                   exclude from weight computation.
    // datacolumn        String. Data column to use for computing weights.
    //                   Supports 'data' or 'corrected'. Minimum match, case
    //                   insensitive. If not provided. 'corrected' is used.
    // slidetimebin      Bool. If true, use a sliding window for binning in
    //                   time.
    // timebin           Double. Width of sliding time window. Not used if
    //                   doslidetime is not supplied or if doslidetime = false;
    StatWtTVI(
        ViImplementation2* inputVii, const casacore::Record &configuration
    );

    virtual ~StatWtTVI();

    virtual casacore::String ViiType() const {
        return casacore::String("StatWt( ") + getVii()->ViiType() + " )";
    };

    void initWeightSpectrum (const casacore::Cube<casacore::Float>& wtspec);

    void initSigmaSpectrum (const casacore::Cube<casacore::Float>& sigspec);

    void next();

    void origin();

    virtual void weightSpectrum(casacore::Cube<casacore::Float>& wtsp) const;

    virtual void sigmaSpectrum(casacore::Cube<casacore::Float>& sigmaSp) const;

    virtual void weight(casacore::Matrix<casacore::Float> & wtmat) const;
    
    virtual void sigma(casacore::Matrix<casacore::Float> & sigmaMat) const;

    virtual void flag(casacore::Cube<casacore::Bool>& flagCube) const;

    virtual void flagRow (casacore::Vector<casacore::Bool> & flagRow) const;

    void summarizeFlagging() const;

    void summarizeStats(
        casacore::Double& mean, casacore::Double& variance
    ) const;

    // Override unimplemented TransformingVi2 version
    void writeBackChanges(VisBuffer2* vb);

    // these are public so that class StatWt can call them. In general, other
    // clients shouldn't call them.
    static casacore::Double getTimeBinWidthInSec(
        const casacore::Quantity& binWidth
    );

    static void checkTimeBinWidth(casacore::Double binWidth);

    /*
    static casacore::Double getTimeBinWidthUsingInterval(
        const casacore::MeasurementSet *const ms, casacore::Int n
    );
    */

protected:

    void originChunks(casacore::Bool forceRewind);

    void nextChunk();
    
private:


    struct BaselineChanBin {
        Baseline baseline = std::make_pair(0, 0);
        casacore::uInt spw = 0;
        StatWtTypes::ChanBin chanBin;
        bool operator<(const BaselineChanBin& other) const {
            if (baseline < other.baseline) {
                return true;
            }
            if (baseline == other.baseline && spw < other.spw) {
                return true;
            }
            return baseline == other.baseline && spw == other.spw
                && chanBin < other.chanBin;
        };
    };
/*
    enum Column {
        // column(s) to use
        // DATA
        DATA,
        // CORRECTED_DATA
        CORRECTED,
        // CORRECTED_DATA - MODEL_DATA
        RESIDUAL,
        // DATA - MODEL_DATA
        RESIDUAL_DATA
    };
*/

    mutable casacore::Bool _weightsComputed = false;
    mutable std::shared_ptr<casacore::Bool> _mustComputeWtSp {};
    mutable casacore::Cube<casacore::Float> _newWtSp {};
    mutable casacore::Matrix<casacore::Float> _newWt {};
    mutable casacore::Cube<casacore::Bool> _newFlag {};
    mutable casacore::Vector<casacore::Bool> _newFlagRow {};
    // the vector represents separate correlations, there will be
    // only one element in the vector if _combineCorr is true
    mutable std::map<BaselineChanBin, casacore::Vector<casacore::Double>>
        _variancesOneShotProcessing {};
    // The key refers to the spw, the value vector refers to the
    // channel numbers within that spw that are the first, last channel pair
    // in their respective bins
    std::map<casacore::Int, std::vector<StatWtTypes::ChanBin>> _chanBins {};
    casacore::Int _minSamp = 2;
    casacore::Bool _combineCorr {false};
    std::shared_ptr<
        casacore::StatisticsAlgorithm<
            casacore::Double, casacore::Array<casacore::Float>::const_iterator,
            casacore::Array<casacore::Bool>::const_iterator,
            casacore::Array<casacore::Double>::const_iterator
        >
    > _statAlg {} ;
    std::shared_ptr<std::pair<casacore::Double, casacore::Double>> _wtrange {};
    // The _chanSelFlags key is the spw. The value is a Cube for convenience
    // for subchunk computations that require the same shaped cube of flags to
    // be applied. The dimension that counts is the second (zero-based 1) as it
    // has length equal to the number of channels in the spw. A value of True
    // indicates that the channel is "flagged", ie should not be used.
    std::map<casacore::uInt, casacore::Cube<casacore::Bool>> _chanSelFlags {};

    mutable size_t _nTotalPts = 0;
    mutable size_t _nNewFlaggedPts = 0;
    mutable size_t _nOrigFlaggedPts = 0;
    mutable StatWtTypes::Column _column = StatWtTypes::CORRECTED;
    mutable std::shared_ptr<
            std::map<casacore::uInt, std::pair<casacore::uInt, casacore::uInt>>
        > _samples {};
    mutable std::set<casacore::uInt> _processedRowIDs {};
    mutable std::vector<std::vector<casacore::Double>> _timeWindowWts {};
    mutable casacore::Cube<casacore::Double> _multiLoopWeights {};
    // if False, the a sliding time window is being used
    casacore::Bool _timeBlockProcessing = true;
    // we can process using classical VI/VB2 algorithm. Only happens if
    // we are not using a sliding time window and if we are not using an
    // integer number of time bins
    casacore::Bool _doOneShot = true;
    // for running time window, for each subchunk, map the rowID (in the MS)
    // to the row index in the chunk
    mutable std::map<casacore::uInt, casacore::uInt>
        _rowIDInMSToRowIndexInChunk {};
    // std::unique_ptr<casacore::Double> _slidingTimeWindowWidth {};
    // if defined means we are using a window width in seconds
    std::shared_ptr<casacore::Double> _binWidthInSeconds {};
    // if defined means we are using an integer number of timestamps for the
    // bin width
    std::unique_ptr<casacore::Int> _nTimeStampsInBin {};

    casacore::Bool _mustComputeSigma = casacore::False;
    casacore::Bool _updateWeight = casacore::True;
    casacore::Bool _noModel = casacore::False;

    std::shared_ptr<StatWtDataAggregator> _dataAggregator {};

    std::shared_ptr<
        casacore::ClassicalStatistics<casacore::Double,
        casacore::Array<casacore::Float>::const_iterator,
        casacore::Array<casacore::Bool>::const_iterator>
    > _wtStats {};

    std::unique_ptr<StatWtVarianceAndWeightCalculator> _varianceComputer;

    // idToChunksNeededByIDMap maps subchunkIDs to the range of subchunk IDs
    // they need. chunkNeededToIDsThatNeedChunkIDMap maps subchunk IDs that are
    // needed to the subchunkIDs that need them. min/max IDs (.first/.second)
    // in both cases
    void _limits(
        std::vector<std::pair<casacore::uInt, casacore::uInt>>& idToChunksNeededByIDMap,
        std::vector<std::pair<casacore::uInt, casacore::uInt>>& chunkNeededToIDsThatNeedChunkIDMap
    ) const;

    // returns True if this chunk has already been processed. This can happen
    // for the last chunk.
    casacore::Bool _checkFirstSubChunk(
        casacore::Int& spw, casacore::Bool& firstTime,
        const VisBuffer2 * const vb
    ) const;

    void _computeWeightSpectrumAndFlags() const;

    const casacore::Cube<casacore::Complex> _dataCube(
        const VisBuffer2 *const vb
    ) const;

    void _gatherAndComputeWeights() const;

    // sliding time bin window, and timebin was specified as a quantity
    void _gatherAndComputeWeightsMultiLoopProcessing() const;

    void _gatherAndComputeWeightsOneShotProcessing() const;

    // combines the flag cube with the channel selection flags (if any)
    casacore::Cube<casacore::Bool> _getResultantFlags(
        casacore::Cube<casacore::Bool>& chanSelFlagTemplate,
        casacore::Cube<casacore::Bool>& chanSelFlags,
        casacore::Bool& initChanSelFlags, casacore::Int spw,
        const casacore::Cube<casacore::Bool>& flagCube
    ) const;

    // CAS-12358
    void _logUsedChannels() const;

    void _computeVariancesOneShotProcessing(
        const std::map<BaselineChanBin, casacore::Cube<casacore::Complex>>& data,
        const std::map<BaselineChanBin, casacore::Cube<casacore::Bool>>& flags,
        const std::map<BaselineChanBin, casacore::Vector<casacore::Double>>& exposures
    ) const;

    void _computeWeightsMultiLoopProcessing(
        const casacore::Cube<casacore::Complex>& data,
        const casacore::Cube<casacore::Bool>& flags,
        const casacore::Vector<casacore::Double>& exposures,
        const std::vector<std::set<casacore::uInt>>& rowMap, casacore::uInt spw
    ) const;

    casacore::Bool _parseConfiguration(const casacore::Record &configuration);
	
    // swaps ant1/ant2 if necessary
    static Baseline _baseline(casacore::uInt ant1, casacore::uInt ant2);

    std::pair<
        casacore::Cube<casacore::Float>, casacore::Cube<casacore::Bool>
    > _getLowerLayerWtSpFlags(size_t& nOrigFlagged) const;

    void _setChanBinMap(casacore::Int binWidth);

    void _setChanBinMap(const casacore::Quantity& binWidth);

    void _setDefaultChanBinMap();

    void _clearCache();

    void _updateWtSpFlags(
        casacore::Cube<casacore::Float>& wtsp,
        casacore::Cube<casacore::Bool>& flags, casacore::Bool& checkFlags,
        const casacore::Slicer& slice, casacore::Float wt
    ) const;

    void _configureStatAlg(const casacore::Record& config);

    void _weightSpectrumFlagsOneShotProcessing(
        casacore::Cube<casacore::Float>& wtsp,
        casacore::Cube<casacore::Bool>& flagCube, casacore::Bool& checkFlags
     ) const;

    void _weightSpectrumFlagsMultiLoopProcessing(
        casacore::Cube<casacore::Float>& wtsp,
        casacore::Cube<casacore::Bool>& flagCube, casacore::Bool& checkFlags
    ) const;

    void _weightSingleChanBinOneShotProcessing(
        casacore::Matrix<casacore::Float>& wtmat, casacore::Int nrows
    ) const;

    void _weightSingleChanBinMultiLoopProcessing(
        casacore::Matrix<casacore::Float>& wtmat, casacore::Int nrows
    ) const;

};

}

}

#endif 

