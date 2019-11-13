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

#ifndef STATWTDATAAGGREGATOR_H_
#define STATWTDATAAGGREGATOR_H_

#include <msvis/MSVis/VisBuffer2.h>
#include <mstransform/TVI/StatWtTypes.h>

#include <map>

namespace casa {

namespace vi {

// Pure virtual base class of various statwt data aggregators. Only StatWt
// needs to use this; developers should not use this code directly.

class StatWtDataAggregator {

public:
    
    using Baseline = std::pair<casacore::uInt, casacore::uInt>;

    StatWtDataAggregator() = delete;

    StatWtDataAggregator(
        const std::map<
            casacore::Int, std::vector<StatWtTypes::ChanBin>
        >& chanBins
    );

    virtual ~StatWtDataAggregator();

    // aggregates the data and computes the weights
    virtual void aggregate() = 0;
    
protected:

    const std::map<casacore::Int, std::vector<StatWtTypes::ChanBin>> _chanBins;

    // swaps ant1/ant2 if necessary
    static Baseline _baseline(casacore::uInt ant1, casacore::uInt ant2);

    // returns True if this chunk has already been processed. This can happen
    // for the last chunk.
    casacore::Bool _checkFirstSubChunk(
        casacore::Int& spw, casacore::Bool& firstTime,
        const vi::VisBuffer2 * const vb
    ) const;

    const casacore::Cube<casacore::Complex> _dataCube(
        const VisBuffer2 *const vb
    ) const;

    // combines the flag cube with the channel selection flags (if any)
    casacore::Cube<casacore::Bool> _getResultantFlags(
        casacore::Cube<casacore::Bool>& chanSelFlagTemplate,
        casacore::Cube<casacore::Bool>& chanSelFlags,
        casacore::Bool& initChanSelFlags, casacore::Int spw,
        const casacore::Cube<casacore::Bool>& flagCube
    ) const;

};

}

}

#endif 

