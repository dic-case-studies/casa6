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
// Aggregates data according to the fundamental VI/VB2 methodology. That is,
// only one subchunk is loaded into memory and aggregated at a time.
// This algorithm is used when timebin is a time quantity (not an int) and
// when slidetimewindow is false.

class StatWtClassicalDataAggregator: public StatWtDataAggregator {

public:
    
    StatWtClassicalDataAggregator() = delete;

    // out of necessity, the passed in pointer like variables are shared with
    // the caller.
    StatWtClassicalDataAggregator(
        ViImplementation2 *const vii,
        std::shared_ptr<casacore::Bool> mustComputeWtSp,
        const std::map<
            casacore::Int, std::vector<StatWtTypes::ChanBin>
        > chanBins
    );

    // aggregates the data and computes the weights
    void aggregate();

private:

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
    std::shared_ptr<casacore::Bool> _mustComputeWtSp;
};

}

}

#endif 

