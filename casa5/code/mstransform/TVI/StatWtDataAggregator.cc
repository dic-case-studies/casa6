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

#include <mstransform/TVI/StatWtDataAggregator.h>

#include <mstransform/TVI/StatWtVarianceAndWeightCalculator.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace casacore;

namespace casa {

namespace vi {

StatWtDataAggregator::StatWtDataAggregator(
    const std::map<casacore::Int, std::vector<StatWtTypes::ChanBin>>& chanBins
) : _chanBins(chanBins) {}

StatWtDataAggregator::Baseline StatWtTVI::_baseline(uInt ant1, uInt ant2) {
    return Baseline(min(ant1, ant2), max(ant1, ant2));
}

Bool StatWtDataAggregator::_checkFirstSubChunk(
    Int& spw, Bool& firstTime, const VisBuffer2 * const vb
) const {
    if (! firstTime) {
        // this chunk has already been checked, it has not
        // been processed previously
        return False;
    }
    const auto& rowIDs = vb->rowIds();
    if (_processedRowIDs.find(rowIDs[0]) == _processedRowIDs.end()) {
        // haven't processed this chunk
        _processedRowIDs.insert(rowIDs[0]);
        // the spw is the same for all subchunks, so it only needs to
        // be set once
        spw = *vb->spectralWindows().begin();
        if (_samples->find(spw) == _samples->end()) {
            (*_samples)[spw].first = 0;
            (*_samples)[spw].second = 0;
        }
        firstTime = False;
        return False;
    }
    else {
        // this chunk has been processed, this can happen at the end
        // when the last chunk is processed twice
        return True;
    }
}

const Cube<Complex> StatWtDataAggregator::_dataCube(
    const VisBuffer2 *const vb
) const {
    // cout << __func__ << endl;
    switch (_column) {
    case CORRECTED:
        return vb->visCubeCorrected();
    case DATA:
        return vb->visCube();
    case RESIDUAL:
        if (_noModel) {
            return vb->visCubeCorrected();
        }
        else {
            return vb->visCubeCorrected() - vb->visCubeModel();
        }
    case RESIDUAL_DATA:
        if(_noModel) {
            return vb->visCube();
        }
        else {
            return vb->visCube() - vb->visCubeModel();
        }
    default:
        ThrowCc("Logic error: column type not handled");
    }
}

Cube<Bool> StatWtDataAggregator::_getResultantFlags(
    Cube<Bool>& chanSelFlagTemplate, Cube<Bool>& chanSelFlags,
    Bool& initTemplate, Int spw, const Cube<Bool>& flagCube
) const {
    if (_chanSelFlags.find(spw) == _chanSelFlags.cend()) {
        // no selection of channels to ignore
        return flagCube;
    }
    if (initTemplate) {
        // this can be done just once per chunk because all the rows
        // in the chunk are guaranteed to have the same spw
        // because each subchunk is guaranteed to have a single
        // data description ID.
        chanSelFlagTemplate = _chanSelFlags.find(spw)->second;
        initTemplate = False;
    }
    auto dataShape = flagCube.shape();
    chanSelFlags.resize(dataShape, False);
    auto ncorr = dataShape[0];
    auto nrows = dataShape[2];
    IPosition start(3, 0);
    IPosition end = dataShape - 1;
    Slicer sl(start, end, Slicer::endIsLast);
    for (uInt corr=0; corr<ncorr; ++corr) {
        start[0] = corr;
        end[0] = corr;
        for (Int row=0; row<nrows; ++row) {
            start[2] = row;
            end[2] = row;
            sl.setStart(start);
            sl.setEnd(end);
            chanSelFlags(sl) = chanSelFlagTemplate;
        }
    }
    return flagCube || chanSelFlags;
}

}

}
