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

#include <mstransform/TVI/StatWtClassicalDataAggregator.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace casacore;

namespace casa {

namespace vi {

StatWtClassicalDataAggregator::StatWtClassicalDataAggregator(
    ViImplementation2 *const vii,
    std::shared_ptr<casacore::Bool> mustComputeWtSp,
    const std::map<casacore::Int, std::vector<StatWtTypes::ChanBin>> chanBins
) : StatWtDataAggregator(chanBins), _vii(vii), _mustComputeWtSp(mustComputeWtSp)
    {}

void StatWtClassicalDataAggregator::aggregate() {
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

}

}
