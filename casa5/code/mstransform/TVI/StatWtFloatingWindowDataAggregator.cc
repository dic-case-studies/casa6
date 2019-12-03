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

#include <mstransform/TVI/StatWtFloatingWindowDataAggregator.h>

#include <scimath/StatsFramework/ClassicalStatistics.h>

#include <msvis/MSVis/ViImplementation2.h>
#include <mstransform/TVI/StatWtTypes.h>


#ifdef _OPENMP
#include <omp.h>
#endif

using namespace casacore;
using namespace std;

namespace casa {

namespace vi {

StatWtFloatingWindowDataAggregator::StatWtFloatingWindowDataAggregator(
    ViImplementation2 *const vii,
    // std::shared_ptr<Bool>& mustComputeWtSp,
    const std::map<casacore::Int, std::vector<StatWtTypes::ChanBin>>& chanBins,
    std::shared_ptr<map<uInt, pair<uInt, uInt>>>& samples,
    StatWtTypes::Column column, Bool noModel,
    const map<uInt, Cube<Bool>>& chanSelFlags, Bool combineCorr,
    shared_ptr<
        ClassicalStatistics<
            Double, Array<Float>::const_iterator,
            Array<Bool>::const_iterator
        >
    >& wtStats,
    shared_ptr<const pair<Double, Double>> wtrange,
    shared_ptr<const Double> binWidthInSeconds,
    shared_ptr<const Int> nTimeStampsInBin, const Bool timeBlockProcessing,
    shared_ptr<
        StatisticsAlgorithm<
            Double, Array<Float>::const_iterator,
            Array<Bool>::const_iterator, Array<Double>::const_iterator
        >
    >& statAlg
) : StatWtDataAggregator(
       vii, chanBins, samples, column, noModel, chanSelFlags, /* mustComputeWtSp,*/
       wtStats, wtrange, combineCorr, statAlg
    ), _binWidthInSeconds(binWidthInSeconds),
    _nTimeStampsInBin(nTimeStampsInBin),
    _timeBlockProcessing(timeBlockProcessing) {
        ThrowIf(
            ! (_binWidthInSeconds || _nTimeStampsInBin ),
            "Logic error: neither binWidthInSeconds "
            "nor nTimeStampsInBin has been specified"
        );
}

StatWtFloatingWindowDataAggregator::~StatWtFloatingWindowDataAggregator() {}

void StatWtFloatingWindowDataAggregator::aggregate() {
    // cout << __func__ << endl;

    //auto* vii = getVii();
    auto* vb = _vii->getVisBuffer();
    // map of rowID (the index of the vector) to rowIDs that should be used to
    // compute the weight for the key rowID
    std::vector<std::set<uInt>> rowMap;
    auto firstTime = True;
    // each subchunk is guaranteed to represent exactly one time stamp
    std::vector<Double> subChunkToTimeStamp;
    // Baseline-subchunk number pair to row index in that chunk
    std::map<std::pair<StatWtTypes::Baseline, uInt>, uInt> baselineSubChunkToIndex;
    //auto halfWidth = *_slidingTimeWindowWidth/2;
    // we build up the chunk data in the chunkData and chunkFlags cubes
    Cube<Complex> chunkData;
    Cube<Bool> chunkFlags;
    // chunkExposures needs to be a Cube rather than a Vector so the individual
    // values have one-to-one relationships to the corresponding data values
    // and so simplify dealing with the exposures which is fundamentally a
    // Vector
    // Cube<Double> chunkExposures;
    std::vector<Double> exposures;
   // cout << __FILE__ << " " << __LINE__ << endl;

    uInt subchunkStartRowNum = 0;
    auto initChanSelTemplate = True;
    Cube<Bool> chanSelFlagTemplate, chanSelFlags;
    // we cannot know the spw until inside the subchunk loop
    Int spw = -1;
    _rowIDInMSToRowIndexInChunk.clear();
    Slicer sl(IPosition(3, 0), IPosition(3, 1));
    auto slStart = sl.start();
    auto slEnd = sl.end();
    std::vector<std::pair<uInt, uInt>> idToChunksNeededByIDMap,
        chunkNeededToIDsThatNeedChunkIDMap;
    _limits(idToChunksNeededByIDMap, chunkNeededToIDsThatNeedChunkIDMap);
    cout << "idToChunksNeededByIDMap " << idToChunksNeededByIDMap << endl;
    cout << "chunkNeededToIDsThatNeedChunkIDMap " << chunkNeededToIDsThatNeedChunkIDMap << endl;
    uInt subChunkID = 0;
    // cout << __FILE__ << " " << __LINE__ << endl;

    for (_vii->origin(); _vii->more(); _vii->next(), ++subChunkID) {
        // cout << "subChunkID " << subChunkID << endl;
        // cout << __FILE__ << " " << __LINE__ << endl;

        if (_checkFirstSubChunk(spw, firstTime, vb)) {
            // cout << __FILE__ << " " << __LINE__ << endl;

            return;
        }
        // cout << "mustcomputewtsp " << _mustComputeWtSp << endl;
        if (! _mustComputeWtSp) {
            _mustComputeWtSp.reset(
                new Bool(
                    vb->existsColumn(VisBufferComponent2::WeightSpectrum)
                )
            );
            cout << "has set _mustComputeWtSp to " << _mustComputeWtSp
                    << " value " << *_mustComputeWtSp << endl;
        }
        // cout << __FILE__ << " " << __LINE__ << endl;

        _rowIDInMSToRowIndexInChunk[*vb->rowIds().begin()] = subchunkStartRowNum;
        const auto& ant1 = vb->antenna1();
        const auto& ant2 = vb->antenna2();
        // [nCorrs, nFreqs, nRows)
        const auto nrows = vb->nRows();
        // there is no guarantee a previous subchunk will be included,
        // eg if the timewidth is small enough
        // This is the first subchunk ID that should be used for averaging
        // grouping data for weight computation of the current subchunk ID.
        const auto firstChunkNeededByCurrentID = idToChunksNeededByIDMap[subChunkID].first;
        const auto lastChunkNeededByCurrentID = idToChunksNeededByIDMap[subChunkID].second;
        const auto firstChunkThatNeedsCurrentID = chunkNeededToIDsThatNeedChunkIDMap[subChunkID].first;
        //const auto lastChunkThatNeedsCurrentID = chunkNeededToIDsThatNeedChunkIDMap[subChunkID].second;
        /*
        cout << "firstChunkNeededByCurrentID " << firstChunkNeededByCurrentID << endl;
        cout << "firstChunkThatNeedsCurrentID " << firstChunkThatNeedsCurrentID << endl;
        cout << "idToChunksNeededByIDMap " << idToChunksNeededByIDMap << endl;
        cout << "chunkNeededToIDsThatNeedChunkIDMap " << chunkNeededToIDsThatNeedChunkIDMap << endl;
        */
        auto subchunkTime = vb->time()[0];
        auto rowInChunk = subchunkStartRowNum;
        pair<StatWtTypes::Baseline, uInt> mypair;
        mypair.second = subChunkID;

        // cout << __FILE__ << " " << __LINE__ << endl;

        for (Int row=0; row<nrows; ++row, ++rowInChunk) {
            // cout << __FILE__ << " " << __LINE__ << endl;

            // loop over rows in sub chunk, grouping baseline specific data
            // together
            const auto baseline = _baseline(ant1[row], ant2[row]);
            mypair.first = baseline;
            baselineSubChunkToIndex[mypair] = rowInChunk;
            std::set<uInt> neededRowNums;
            neededRowNums.insert(rowInChunk);
            // cout << __FILE__ << " " << __LINE__ << endl;
            if (subChunkID > 0) {

                auto s = min(
                    firstChunkNeededByCurrentID, firstChunkThatNeedsCurrentID
                );
                // cout << "subchunk start " << subchunk << endl;
                auto tpair = mypair;
                for (; s < subChunkID; ++s) {
                    // cout << __FILE__ << " " << __LINE__ << endl;

                    const auto myend = baselineSubChunkToIndex.end();
                    // cout << __FILE__ << " " << __LINE__ << endl;

                    tpair.second = s;
                    // cout << __FILE__ << " " << __LINE__ << endl;

                    const auto iter = baselineSubChunkToIndex.find(tpair);
                    auto found = iter != myend;
                    /*
                    if (ant1[row] == 0 && ant2[row] == 1) {
                        cout << "s=" << s << " found " << found << endl;
                    }
                    */
                    // cout << __FILE__ << " " << __LINE__ << endl;

                    if (found) {
                        // cout << "map " << baselineSubChunkToIndex << endl;
                        // cout << "test found " << *testFound << endl;
                        // cout << __FILE__ << " " << __LINE__ << endl;

                        const auto existingRowNum = iter->second;
                        // cout << __FILE__ << " " << __LINE__ << endl;

                        if (
                            s >= firstChunkNeededByCurrentID
                            && s <= lastChunkNeededByCurrentID
                        ) {
                           // cout << __FILE__ << " " << __LINE__ << endl;

                            // The subchunk data is needed for computation
                            // of the current subchunkID's weights
                            neededRowNums.insert(existingRowNum);
                            // cout << __FILE__ << " " << __LINE__ << endl;

                        }
                        // cout << "subchunk " << subchunk << endl;
                        // cout << "firstChunkThatNeedsCurrentID " << firstChunkThatNeedsCurrentID << endl;
                        if (
                            idToChunksNeededByIDMap[s].first <= subChunkID
                            && idToChunksNeededByIDMap[s].second >= subChunkID
                        ) {
                            // cout << __FILE__ << " " << __LINE__ << endl;
                            // cout << "rowMap size " << rowMap.size() << endl;
                            // cout << "existingRowNum " << existingRowNum << endl;
                            // cout << "rowInChunk " << rowInChunk << endl;
                            // if (rowMap.size() <= existingRowNum) {
                                // rowMap.resize
                            // }
                            rowMap[existingRowNum].insert(rowInChunk);
                            // cout << __FILE__ << " " << __LINE__ << endl;

                        }
                    }
                    // cout << __FILE__ << " " << __LINE__ << endl;

                }
            }
            // cout << __FILE__ << " " << __LINE__ << endl;
            rowMap.push_back(neededRowNums);
            // debug
            /*
            if (baseline == Baseline(4, 6)) {
                auto myrow = rowMap.size() - 1;
                cout << "row num " << myrow << " included rows " << rowMap[myrow] << endl;
            }
            */
            // cout << __FILE__ << " " << __LINE__ << endl;
            /*
            {
                if (baseline == StatWtTypes::Baseline(0,1)) {
                    auto debugPair = mypair;
                    for (uInt pp=0; pp<=subChunkID; ++pp) {
                        debugPair.second = pp;
                        auto myrownum = baselineSubChunkToIndex[debugPair];
                        cout << "pp " << pp << " pair " << debugPair << " myrownum "
                            << myrownum << endl;
                        cout << "rowMap[myrownum] " << rowMap[myrownum] << endl;
                    }
                }
            }
            */
        }
        // cout << __FILE__ << " " << __LINE__ << endl;
        const auto dataCube = _dataCube(vb);
        const auto resultantFlags = _getResultantFlags(
            chanSelFlagTemplate, chanSelFlags, initChanSelTemplate,
            spw, vb->flagCube()
        );
        const auto myExposures = vb->exposure().tovector();
        exposures.insert(
            exposures.end(), myExposures.begin(), myExposures.end()
        );
        const auto cubeShape = dataCube.shape();
        // Cube<Double> resultExposures(cubeShape);
        IPosition sliceStart(3, 0);
        auto sliceEnd = cubeShape - 1;
        // Slicer exposureSlice(sliceStart, sliceEnd, Slicer::endIsLast);
        // cout << __FILE__ << " " << __LINE__ << endl;

        /*
        for (uInt jj=0; jj<cubeShape[2]; ++jj) {
            sliceStart[2] = jj;
            sliceEnd[2] = jj;
            exposureSlice.setStart(sliceStart);
            exposureSlice.setEnd(sliceEnd);

            // set all exposures in the slice to the same value

            resultExposures(exposureSlice) = exposures[jj];
        }
        */

        // build up chunkData and chunkFlags one subchunk at a time
        // cout << __FILE__ << " " << __LINE__ << endl;

        if (chunkData.empty()) {
            chunkData = dataCube;
            chunkFlags = resultantFlags;
            //chunkExposures = resultExposures;
        }
        else {
            auto newShape = chunkData.shape();
            newShape[2] += nrows;
            chunkData.resize(newShape, True);
            chunkFlags.resize(newShape, True);
            //chunkExposures.resize(newShape, True);
            slStart[2] = subchunkStartRowNum;
            sl.setStart(slStart);
            slEnd = newShape - 1;
            sl.setEnd(slEnd);
            chunkData(sl) = dataCube;
            chunkFlags(sl) = resultantFlags;
            //chunkExposures(sl) = resultExposures;
        }
        subChunkToTimeStamp.push_back(subchunkTime);
        subchunkStartRowNum += nrows;


    }
    /*
    cout << "rowMap[0] " << rowMap[0] << endl;
    if (subChunkID > 9) {
        cout << "rowMap size " << rowMap.size() << endl;
        cout << "rowMap[550] " << rowMap[550] << endl;
    }
    // cout << "chunkexposres shape " << chunkExposures.shape() << endl;
    cout << __FILE__ << " " << __LINE__ << endl;
    */
    _computeWeightsMultiLoopProcessing(
        chunkData, chunkFlags, Vector<Double>(exposures), rowMap, spw
    );
    // cout << __FILE__ << " " << __LINE__ << endl;

}

void StatWtFloatingWindowDataAggregator::_computeWeightsMultiLoopProcessing(
    const Cube<Complex>& data, const Cube<Bool>& flags,
    const Vector<Double>& exposures, const std::vector<std::set<uInt>>& rowMap,
    uInt spw
) const {
    // cout << "data shape " << data.shape() << endl;
    // cout << "flags shape " << flags.shape() << endl;
    // cout << "exposures shape " << exposures.shape() << endl;
    // cout << "rowMap size " << rowMap.size() << endl;
    auto chunkShape = data.shape();
    const auto nActCorr = chunkShape[0];
    const auto ncorr = _combineCorr ? 1 : nActCorr;
    const auto& chanBins = _chanBins.find(spw)->second;
    // cout << __FILE__ << " " << __LINE__ << endl;

    _multiLoopWeights.resize(
        IPosition(3, ncorr, chanBins.size(), chunkShape[2]),
        False
    );
    const auto nRows = rowMap.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (size_t iRow=0; iRow<nRows; ++iRow) {
        IPosition chunkSliceStart(3, 0);
        auto chunkSliceLength = chunkShape;
        chunkSliceLength[2] = 1;
        Slicer chunkSlice(
            chunkSliceStart, chunkSliceLength, Slicer::endIsLength
        );
        auto chunkSliceEnd = chunkSlice.end();
        auto appendingSlice = chunkSlice;
        auto appendingSliceStart = appendingSlice.start();
        auto appendingSliceEnd = appendingSlice.end();
        auto intraChunkSlice = appendingSlice;
        auto intraChunkSliceStart = intraChunkSlice.start();
        auto intraChunkSliceEnd = intraChunkSlice.end();
        intraChunkSliceEnd[0] = nActCorr - 1;
        const auto& rowsToInclude = rowMap[iRow];
        auto dataShape = chunkShape;
        dataShape[2] = rowsToInclude.size();
        Cube<Complex> dataArray(dataShape);
        Cube<Bool> flagArray(dataShape);
        auto siter = rowsToInclude.begin();
        auto send = rowsToInclude.end();
        // cout << __FILE__ << " " << __LINE__ << endl;

        Vector<Double> exposureVector(rowsToInclude.size(), 0);
        uInt n = 0;
        // create an array with only the rows that should
        // be used in the computation of weights for the
        // current row
        for (; siter!=send; ++siter, ++n) {
            exposureVector[n] = exposures[*siter];
            appendingSliceStart[2] = n;
            appendingSlice.setStart(appendingSliceStart);
            appendingSliceEnd[2] = n;
            appendingSlice.setEnd(appendingSliceEnd);
            chunkSliceStart[2] = *siter;
            chunkSlice.setStart(chunkSliceStart);
            chunkSliceEnd[2] = *siter;
            chunkSlice.setEnd(chunkSliceEnd);
            dataArray(appendingSlice) = data(chunkSlice);
            flagArray(appendingSlice) = flags(chunkSlice);
        }
        // cout << __FILE__ << " " << __LINE__ << endl;

        // slice up for correlations and channel binning
        intraChunkSliceEnd[2] = dataShape[2] - 1;
        for (uInt corr=0; corr<ncorr; ++corr) {
            // cout << __FILE__ << " " << __LINE__ << endl;

            if (! _combineCorr) {
                intraChunkSliceStart[0] = corr;
                intraChunkSliceEnd[0] = corr;
            }
            auto citer = chanBins.begin();
            auto cend = chanBins.end();
            auto iChanBin = 0;
            for (; citer!=cend; ++citer, ++iChanBin) {
                // cout << __FILE__ << " " << __LINE__ << endl;

                intraChunkSliceStart[1] = citer->start;
                intraChunkSliceEnd[1] = citer->end;
                intraChunkSlice.setStart(intraChunkSliceStart);
                intraChunkSlice.setEnd(intraChunkSliceEnd);
                // cout << __FILE__ << " " << __LINE__ << endl;
                if (! _varianceComputer) {
                    cout << "variance computer not set" << endl;
                }
                _multiLoopWeights(corr, iChanBin, iRow)
                    = _varianceComputer->computeWeight(
                        dataArray(intraChunkSlice), flagArray(intraChunkSlice),
                        exposureVector, spw, exposures[iRow]
                    );
                // cout << __FILE__ << " " << __LINE__ << endl;

            }
            // cout << __FILE__ << " " << __LINE__ << endl;

        }
        // cout << __FILE__ << " " << __LINE__ << endl;

    }
    //cout << __FILE__ << " " << __LINE__ << endl;

}

void StatWtFloatingWindowDataAggregator::weightSpectrumFlags(
    Cube<Float>& wtsp, Cube<Bool>& flagCube, Bool& checkFlags,
    const Vector<Int>& ant1, const Vector<Int>& /* ant2 */,
    const Vector<Int>& spws, const Vector<Double>& /* exposures */,
    const Vector<uInt>& rowIDs
) const {

/*
void StatWtFloatingWindowDataAggregator::weightSpectrumFlags(
    Cube<Float>& wtsp, Cube<Bool>& flagCube, Bool& checkFlags
) const {
*/
    // fish out the rows relevant to this subchunk
    // Vector<uInt> rowIDs;
    // getRowIds(rowIDs);
    const auto start = _rowIDInMSToRowIndexInChunk.find(*rowIDs.begin());
    ThrowIf(
        start == _rowIDInMSToRowIndexInChunk.end(),
        "Logic Error: Cannot find requested subchunk in stored chunk"
    );
    // this is the row index in the chunk
    auto chunkRowIndex = start->second;
    // auto chunkRowEnd = chunkRowIndex + nRows();
    auto chunkRowEnd = chunkRowIndex + ant1.size();
    Slicer slice(IPosition(3, 0), flagCube.shape(), Slicer::endIsLength);
    auto sliceStart = slice.start();
    auto sliceEnd = slice.end();
    auto nCorrBins = _combineCorr ? 1 : flagCube.shape()[0];
    // Vector<Int> spws;
    // spectralWindows(spws);
    auto spw = *spws.begin();
    const auto& chanBins = _chanBins.find(spw)->second;
    auto subChunkRowIndex = 0;
    for (; chunkRowIndex < chunkRowEnd; ++chunkRowIndex, ++subChunkRowIndex) {
        sliceStart[2] = subChunkRowIndex;
        sliceEnd[2] = subChunkRowIndex;
        auto iChanBin = 0;
        for (const auto& chanBin: chanBins) {
            sliceStart[1] = chanBin.start;
            sliceEnd[1] = chanBin.end;
            auto corr = 0;
            for (; corr < nCorrBins; ++corr) {
                if (! _combineCorr) {
                    sliceStart[0] = corr;
                    sliceEnd[0] = corr;
                }
                slice.setStart(sliceStart);
                slice.setEnd(sliceEnd);
                _updateWtSpFlags(
                    wtsp, flagCube, checkFlags, slice,
                    _multiLoopWeights(corr, iChanBin, chunkRowIndex)
                );
            }
            ++iChanBin;
        }
    }
}

void StatWtFloatingWindowDataAggregator::_limits(
    std::vector<std::pair<uInt, uInt>>& idToChunksNeededByIDMap,
    std::vector<std::pair<uInt, uInt>>& chunkNeededToIDsThatNeedChunkIDMap
) const {
    // auto* vii = getVii();
    auto* vb = _vii->getVisBuffer();
    pair<uInt, uInt> p, q;
    uInt nTimes = _vii->nTimes();
    if (_nTimeStampsInBin) {
        // fixed number of time stamps specified
        if (_timeBlockProcessing) {
            // cout << "_nTimeStampsInBin && _timeBlockProcessing is True" << endl;
            // integer division
            uInt nBlocks = nTimes/(*_nTimeStampsInBin);
            // cout << "nTimes " << nTimes << endl;
            // cout << "*_nTimeStampsInBin " <<  *_nTimeStampsInBin << endl;
            if (nTimes % *_nTimeStampsInBin > 0) {
                ++nBlocks;
            }
            // cout << "nBlocks " << nBlocks << endl;
            uInt subChunkCount = 0;
            for (uInt blockCount = 0; blockCount < nBlocks; ++blockCount) {
                if ((subChunkCount + *_nTimeStampsInBin <= nTimes)) {
                    p.first = subChunkCount;
                    p.second = subChunkCount + *_nTimeStampsInBin - 1;
                }
                else {
                    // chunk upper edge
                    p.first = nTimes < (uInt)*_nTimeStampsInBin
                        ? 0 : nTimes - *_nTimeStampsInBin;
                    p.second = nTimes - 1;
                }
                q = p;
                for (uInt i=subChunkCount; i<=p.second; ++i, ++subChunkCount) {
                    idToChunksNeededByIDMap.push_back(p);
                    chunkNeededToIDsThatNeedChunkIDMap.push_back(q);
                }
            }
        }
        else {
            // sliding time window, fixed number of time stamps (timebin
            // specified as int
            const auto isEven = *_nTimeStampsInBin % 2 == 0;
            // integer division
            const uInt halfTimeBin = *_nTimeStampsInBin/2;
            const auto nBefore = isEven
                ? (halfTimeBin) : (*_nTimeStampsInBin - 1)/2;
            const auto nAfter = isEven ? nBefore + 1 : nBefore;
            // integer division
            // p.first is the first sub chunk needed by the current index.
            // p.second is the first sub chunk that needs the current index
            for (uInt i=0; i<nTimes; ++i) {
                if (i <= nBefore) {
                    p.first = 0;
                }
                else if (i >= nTimes - nAfter) {
                    p.first = nTimes - *_nTimeStampsInBin;
                }
                else {
                    p.first = i - nBefore;
                }
                if ((uInt)i >= nTimes - nAfter) {
                    p.second = nTimes - 1;
                }
                else {
                    p.second = i + nAfter;
                }
                if (i <= nAfter) {
                    q.first = 0;
                }
                else {
                    q.first = i - nAfter;
                }
                if (i + nAfter < nTimes) {
                    q.second = i + nAfter;
                }
                else {
                    q.second = nTimes - 1;
                }
                idToChunksNeededByIDMap.push_back(p);
                chunkNeededToIDsThatNeedChunkIDMap.push_back(q);
            }
        }
    }
    else {
        if (_timeBlockProcessing) {
            // shouldn't get in here
            ThrowCc("Logic error: shouldn't have gotten into this code block");
        }
        else {
            ThrowIf(
                ! _binWidthInSeconds,
                "Logic error: _binWidthInSeconds not defined"
            );
            auto halfBinWidth = *_binWidthInSeconds/2;
            vector<Double> subChunkTimes;
            uInt subChunkCount = 0;
            for (_vii->origin(); _vii->more(); _vii->next(), ++subChunkCount) {
                // all times in a subchunk are the same
                auto mytime = vb->time()[0];
                subChunkTimes.push_back(mytime);
                if (
                    subChunkCount == 0
                    || mytime - halfBinWidth <= subChunkTimes[0]
                ) {
                    p.first = 0;
                }
                else {
                    auto it = std::lower_bound(
                        subChunkTimes.cbegin(), subChunkTimes.cend(),
                        mytime - halfBinWidth
                    );
                    ThrowIf(
                        it == subChunkTimes.end(),
                        "Logic Error for std::lower_bound()"
                    );
                    p.first = std::distance(subChunkTimes.cbegin(), it);
                }
            }
            for (uInt subChunk=0; subChunk<=subChunkCount; ++subChunk) {
                // get upper bound
                auto upit = std::upper_bound(
                    subChunkTimes.cbegin(), subChunkTimes.cend(),
                    subChunkTimes[subChunk] - halfBinWidth
                );
                p.second = std::distance(subChunkTimes.cbegin(), upit);
            }
            q = p;
            idToChunksNeededByIDMap.push_back(p);
            chunkNeededToIDsThatNeedChunkIDMap.push_back(q);
        }
    }
}


/*
void StatWtFloatingWindowDataAggregator::aggregate() {
    // Drive NEXT LOWER layer's ViImpl to gather data into allvis:
    //  Assumes all sub-chunks in the current chunk are to be used
    //   for the variance calculation
    //  Essentially, we are sorting the incoming data into
    //   allvis, to enable a convenient variance calculation
    // cout << __FILE__ << " " << __LINE__ << endl;
    _variancesOneShotProcessing.clear();
    auto* vb = _vii->getVisBuffer();
    std::map<BaselineChanBin, Cube<Complex>> data;
    std::map<BaselineChanBin, Cube<Bool>> flags;
    std::map<BaselineChanBin, Vector<Double>> exposures;
    IPosition blc(3, 0);
    auto trc = blc;
    auto initChanSelTemplate = True;
    Cube<Bool> chanSelFlagTemplate, chanSelFlags;
    auto firstTime = True;
    // we cannot know the spw until we are in the subchunks loop
    Int spw = -1;
    // Vector<Double> exposureVector;
    // cout << __FILE__ << " " << __LINE__ << endl;

    for (_vii->origin(); _vii->more(); _vii->next()) {
        if (_checkFirstSubChunk(spw, firstTime, vb)) {
            return;
        }
        if (! _mustComputeWtSp) {
            _mustComputeWtSp.reset(
                new Bool(
                    vb->existsColumn(VisBufferComponent2::WeightSpectrum)
                )
            );
        }
        // cout << __FILE__ << " " << __LINE__ << endl;

        const auto& ant1 = vb->antenna1();
        const auto& ant2 = vb->antenna2();
        // [nCorr, nFreq, nRows)
        const auto& dataCube = _dataCube(vb);
        const auto& flagCube = vb->flagCube();
        const auto dataShape = dataCube.shape();
        const auto& exposureVector = vb->exposure();
        const auto nrows = vb->nRows();
        const auto npol = dataCube.nrow();
        const auto resultantFlags = _getResultantFlags(
            chanSelFlagTemplate, chanSelFlags, initChanSelTemplate,
            spw, flagCube
        );
        // cout << __FILE__ << " " << __LINE__ << endl;

        auto bins = _chanBins.find(spw)->second;
        BaselineChanBin blcb;
        blcb.spw = spw;
        IPosition dataCubeBLC(3, 0);
        auto dataCubeTRC = dataCube.shape() - 1;
        // cout << __FILE__ << " " << __LINE__ << endl;

        for (Int row=0; row<nrows; ++row) {
            dataCubeBLC[2] = row;
            dataCubeTRC[2] = row;
            blcb.baseline = _baseline(ant1[row], ant2[row]);
            auto citer = bins.cbegin();
            auto cend = bins.cend();
            for (; citer!=cend; ++citer) {
                dataCubeBLC[1] = citer->start;
                dataCubeTRC[1] = citer->end;
                blcb.chanBin.start = citer->start;
                blcb.chanBin.end = citer->end;
                auto dataSlice = dataCube(dataCubeBLC, dataCubeTRC);
                auto flagSlice = resultantFlags(dataCubeBLC, dataCubeTRC);
                if (data.find(blcb) == data.end()) {
                    data[blcb] = dataSlice;
                    flags[blcb] = flagSlice;
                    exposures[blcb] = Vector<Double>(1, exposureVector[row]);
                }
                else {
                    auto myshape = data[blcb].shape();
                    auto nplane = myshape[2];
                    auto nchan = myshape[1];
                    data[blcb].resize(npol, nchan, nplane+1, True);
                    flags[blcb].resize(npol, nchan, nplane+1, True);
                    exposures[blcb].resize(nplane+1, True);
                    trc = myshape - 1;
                    // because we've extended the cube by one plane since
                    // myshape was determined.
                    ++trc[2];
                    blc[2] = trc[2];
                    data[blcb](blc, trc) = dataSlice;
                    flags[blcb](blc, trc) = flagSlice;
                    exposures[blcb][trc[2]] = exposureVector[row];
                }
            }
        }
        // cout << __FILE__ << " " << __LINE__ << endl;

    }
    // cout << __FILE__ << " " << __LINE__ << endl;

    _computeVariancesOneShotProcessing(data, flags, exposures);
    // cout << __FILE__ << " " << __LINE__ << endl;

}
*/

/*
void StatWtFloatingWindowDataAggregator::_computeVariancesOneShotProcessing(
    const map<BaselineChanBin, Cube<Complex>>& data,
    const map<BaselineChanBin, Cube<Bool>>& flags,
    const map<BaselineChanBin, Vector<Double>>& exposures
) const {
    //cout << "exposures size " << exposures.size() << endl;
    //cout << "flags size " << flags.size() << endl;
    //cout << "data size " << data.size() << endl;
    auto diter = data.cbegin();
    auto dend = data.cend();
    const auto nActCorr = diter->second.shape()[0];
    const auto ncorr = _combineCorr ? 1 : nActCorr;
    // spw will be the same for all members
    const auto& spw = data.begin()->first.spw;
    std::vector<BaselineChanBin> keys(data.size());
    auto idx = 0;
    // cout << __FILE__ << " " << __LINE__ << endl;

    for (; diter!=dend; ++diter, ++idx) {
        const auto& blcb = diter->first;
        keys[idx] = blcb;
        _variancesOneShotProcessing[blcb].resize(ncorr);
    }
    // cout << __FILE__ << " " << __LINE__ << endl;

    auto n = keys.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (size_t i=0; i<n; ++i) {
        // cout << __FILE__ << " " << __LINE__ << endl;

        auto blcb = keys[i];
        // cout << __FILE__ << " " << __LINE__ << endl;

        auto dataForBLCB = data.find(blcb)->second;
        // cout << __FILE__ << " " << __LINE__ << endl;

        auto flagsForBLCB = flags.find(blcb)->second;
        // cout << __FILE__ << " " << __LINE__ << endl;

        auto exposuresForBLCB = exposures.find(blcb)->second;
        // cout << __FILE__ << " " << __LINE__ << endl;

        for (uInt corr=0; corr<ncorr; ++corr) {
            IPosition start(3, 0);
            auto end = dataForBLCB.shape() - 1;
            if (! _combineCorr) {
                start[0] = corr;
                end[0] = corr;
            }
            Slicer slice(start, end, Slicer::endIsLast);
            _variancesOneShotProcessing[blcb][corr]
                = _varianceComputer->computeVariance(
                    dataForBLCB(slice), flagsForBLCB(slice),
                    exposuresForBLCB, spw
                );
        }
        // cout << __FILE__ << " " << __LINE__ << endl;

    }
    //cout << __FILE__ << " " << __LINE__ << endl;

}
*/
/*

void StatWtFloatingWindowDataAggregator::weightSpectrumFlags(
    Cube<Float>& wtsp, Cube<Bool>& flagCube, Bool& checkFlags,
    const Vector<Int>& ant1, const Vector<Int>& ant2, const Vector<Int>& spws,
    const Vector<Double>& exposures
) const {
    //Vector<Int> ant1, ant2, spws;
    //Vector<Double> exposures;
    //antenna1(ant1);
    //antenna2(ant2);
    //spectralWindows(spws);
    //exposure(exposures);
    Slicer slice(IPosition(3, 0), flagCube.shape(), Slicer::endIsLength);
    auto sliceStart = slice.start();
    auto sliceEnd = slice.end();
    //auto nrows = nRows();
    auto nrows = ant1.size();
    for (uInt i=0; i<nrows; ++i) {
        sliceStart[2] = i;
        sliceEnd[2] = i;
        BaselineChanBin blcb;
        blcb.baseline = _baseline(ant1[i], ant2[i]);
        auto spw = spws[i];
        blcb.spw = spw;
        auto bins = _chanBins.find(spw)->second;
        for (const auto& bin: bins) {
            sliceStart[1] = bin.start;
            sliceEnd[1] = bin.end;
            blcb.chanBin = bin;
            auto variances = _variancesOneShotProcessing.find(blcb)->second;
            auto ncorr = variances.size();
            Vector<Double> weights = exposures[i]/variances;
            for (uInt corr=0; corr<ncorr; ++corr) {
                if (! _combineCorr) {
                    sliceStart[0] = corr;
                    sliceEnd[0] = corr;
                }
                slice.setStart(sliceStart);
                slice.setEnd(sliceEnd);
                _updateWtSpFlags(
                    wtsp, flagCube, checkFlags, slice, weights[corr]
                );
            }
        }
    }
}
*/

}

}
