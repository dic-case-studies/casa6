//# StatWtTVI.cc: This file contains the implementation of the StatWtTVI class.
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

#include <mstransform/TVI/StatWtFloatingWindowAggregator.h>

#include <scimath/StatsFramework/ClassicalStatistics.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace casacore;
using namespace std;

namespace casa {

namespace vi {

StatWtFloatingWindowAggregator::StatWtFloatingWindowAggregator(
    ViImplementation2 *const vii,
    std::shared_ptr<const Bool> mustComputeWtSp,
    const std::map<casacore::Int, std::vector<StatWtTypes::ChanBin>>& chanBins,
    std::shared_ptr<map<uInt, pair<uInt, uInt>>> samples,
    StatWtTypes::Column column, Bool noModel,
    const map<uInt, Cube<Bool>>& chanSelFlags, Bool combineCorr,
    shared_ptr<
        ClassicalStatistics<
            Double, Array<Float>::const_iterator,
            Array<Bool>::const_iterator
        >
    > wtStats,
    std::shared_ptr<
        const std::pair<casacore::Double, casacore::Double>
    > wtrange
) : StatWtDataAggregator(
       chanBins, samples, column, noModel, chanSelFlags, mustComputeWtSp,
       wtStats, wtrange
    ),
    _vii(vii), _combineCorr(combineCorr) {}

StatWtFloatingWindowAggregator::~StatWtFloatingWindowAggregator() {}

void StatWtFloatingWindowAggregator::aggregate() {
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

void StatWtFloatingWindowAggregator::_computeVariancesOneShotProcessing(
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

void StatWtFloatingWindowAggregator::weightSpectrumFlags(
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

}

}
