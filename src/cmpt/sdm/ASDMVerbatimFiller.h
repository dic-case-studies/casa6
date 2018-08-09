#ifndef _ASDMVERBATIMFILLER_H_
#define _ASDMVERBATIMFILLER_H_
#include "ASDMTableBase.h"
#include "ASDMTables.h"

#include <set>
#include <alma/ASDM/ASDM.h>

using namespace std;
using namespace asdm;

namespace casac {
    class ASDMVerbatimFiller {
      public:
        virtual ~ASDMVerbatimFiller();
        //  ASDMVerbatimFiller(casacore::MS* ms_p, const set<const ASDM_TABLE_BASE*>& table); 
        ASDMVerbatimFiller(casacore::MS* ms_p, const set<ASDM_TABLE_BASE*>& table); 
        void fill(const ASDM& asdm);
      private:
        set<ASDM_TABLE_BASE*> table_;
        ASDMVerbatimFiller();
    };
}

#endif // _ASDMVERBATIMFILLER_H_

