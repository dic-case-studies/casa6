#ifndef BDF2ASDMSTMANINDEX
#define BDF2ASDMSTMANINDEX
#include <alma/Enumerations/CPrimitiveDataType.h>
#include <casacore/casa/Containers/Block.h>
#include <asdmstman/AsdmIndex.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/ArrayUtil.h>
#include <casacore/casa/Arrays/ArrayLogical.h>
#include <casacore/casa/IO/AipsIO.h>
#include <casacore/casa/Containers/BlockIO.h>
#include <map>

#include <stdint.h>
#include <algorithm>

namespace casac {

    class BDF2AsdmStManIndex {
      public:
        BDF2AsdmStManIndex();
        BDF2AsdmStManIndex(const std::vector<std::string>& bdfNames, bool isBigEndian, const string& fname);
        virtual ~BDF2AsdmStManIndex();
        void init (const std::vector<std::string>& bdfNames, bool isBigEndian, const string& fname);
        void setNumberOfDataDescriptions(unsigned int numberOfDataDescriptions);
        void						done();
        void						clearIndexes();
        void						clearAutoIndexes();
        void						clearCrossIndexes();
        void appendAutoIndex(unsigned int             iDD,
                             const string&		bdfName,
                             unsigned int		nBl,
                             unsigned int		nSpw,
                             unsigned int		nChan,
                             unsigned int		nPol,
                             unsigned int		stepBl,
                             unsigned int		stepSpw,
                             const std::vector<double>&	scaleFactors,
                             uint64_t	         	fileOffset,
                             uint32_t                 spwOffset);

        void appendWVRIndex(unsigned int              iDD,
                            const string&		bdfName,
                            unsigned int		nBl,
                            unsigned int		nSpw,
                            unsigned int		nChan,
                            unsigned int		nPol,
                            unsigned int		stepBl,
                            unsigned int		stepSpw,
                            const std::vector<double>&	scaleFactors,
                            uint64_t	         	fileOffset,
                            uint32_t                  spwOffset);

        void appendCrossIndex(unsigned int            iDD,
                              const string&		bdfName,
                              unsigned int		nBl,
                              unsigned int		nSpw,
                              unsigned int		nChan,
                              unsigned int		nPol,
                              unsigned int		stepBl,
                              unsigned int		stepSpw,
                              const std::vector<double>&	scaleFactors,
                              uint64_t		fileOffset,
                              uint32_t                spwOffset,
                              PrimitiveDataTypeMod::PrimitiveDataType       dataType);

        void	dumpAutoCross();
        void	dumpCrossAuto();

        static int version();

      private:
        uint32_t                      numberOfDataDescriptions;
        casacore::Block<casacore::String>	bdfNames;
        casacore::String			fname;
        std::map<std::string, int>	s2i_m;
        std::vector<std::vector<casa::AsdmIndex> >	autoIndexes_vv;
        std::vector<std::vector<casa::AsdmIndex> >	crossIndexes_vv;
        std::vector<casa::AsdmIndex>  allIndexes_v;
        bool fileAttached;
        casacore::AipsIO			aio;
        uint64_t			MSMainRowNumber;
    };
}

#endif
