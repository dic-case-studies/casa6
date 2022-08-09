#include <sstream>
#include <alma/ASDM/Frequency.h>
#include <alma/Enumtcl/Enum.hpp>
#include <alma/Enumtcl/CorrelationMode.h>
#include <alma/Enumtcl/SpectralResolutionType.h>
#include <alma/Enumtcl/TimeSampling.h>
#include <alma/Enumtcl/AtmPhaseCorrection.h>
#include <alma/ASDMBinaries/SDMDataViews.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>

namespace asdm {
    class ExecBlockRow;
    class SpectralWindowRow;
    class ASDM;
    class SysPowerRow;
    class MainRow;
    class ScanRow;
}

namespace sdmbin {
    class SDMBinData;
}

namespace casac {
    class ASDM2MSFiller;
    class UvwCoords;
    class MSFlagEval;
}
    
namespace FrequencyReferenceCodeMod {
    class Frequency;
}
