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

#ifndef STATWTCLASSICALDATAAGGREGATOR_H_
#define STATWTCLASSICALDATAAGGREGATOR_H_

#include <mstransform/TVI/StatWtDataAggregator.h>
#include <msvis/MSVis/TransformingVi2.h>

#include <mstransform/TVI/StatWtTVI.h>

#include <map>

namespace casa {

namespace vi {

// Used by StatWt. Developers should not directly use this API.
// Aggregates data allowing a sliding window inside the chunk, which allows
// aggregating data that span subchunks. This algorithm is used when timebin is
// an int and or when slidetimewindow is true.

class StatWtFloatingWindowAggregator: public StatWtDataAggregator {

public:
    
    StatWtFloatingWindowAggregator() = delete;

    // out of necessity, the passed in pointer like variables are shared with
    // the caller.
    StatWtFloatingWindowAggregator(
        ViImplementation2 *const vii,
        std::shared_ptr<const casacore::Bool> mustComputeWtSp,
        const std::map<
            casacore::Int, std::vector<StatWtTypes::ChanBin>
        >& chanBins,
        std::shared_ptr<
            std::map<casacore::uInt, std::pair<casacore::uInt, casacore::uInt>>
        > samples,
        StatWtTypes::Column column, casacore::Bool noModel,
        const std::map<casacore::uInt, casacore::Cube<casacore::Bool>>&
            chanSelFlags,
        casacore::Bool combineCorr,
        std::shared_ptr<
            casacore::ClassicalStatistics<casacore::Double,
            casacore::Array<casacore::Float>::const_iterator,
            casacore::Array<casacore::Bool>::const_iterator>
        > wtStats,
        std::shared_ptr<
            const std::pair<casacore::Double, casacore::Double>
        > wtrange
    );

    ~StatWtFloatingWindowAggregator();

    // aggregates the data and computes the weights
    void aggregate();

    void weightSpectrumFlags(
        casacore::Cube<casacore::Float>& wtsp,
        casacore::Cube<casacore::Bool>& flagCube, casacore::Bool& checkFlags,
        const casacore::Vector<casacore::Int>& ant1,
        const casacore::Vector<casacore::Int>& ant2,
        const casacore::Vector<casacore::Int>& spws,
        const casacore::Vector<casacore::Double>& exposures
     ) const;

private:
    /*
    struct BaselineChanBin {
        vi::StatWtTVI::Baseline baseline = std::make_pair(0, 0);
        casacore::uInt spw = 0;
        vi::StatWtTypes::ChanBin chanBin;
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

    mutable std::map<BaselineChanBin, casacore::Vector<casacore::Double>>
            _variancesOneShotProcessing {};

    ViImplementation2 *const _vii;
    const casacore::Bool _combineCorr;

    void _computeVariancesOneShotProcessing(
        const std::map<BaselineChanBin, casacore::Cube<casacore::Complex>>& data,
        const std::map<BaselineChanBin, casacore::Cube<casacore::Bool>>& flags,
        const std::map<BaselineChanBin, casacore::Vector<casacore::Double>>& exposures
    ) const;
    */
};

}

}

#endif 

