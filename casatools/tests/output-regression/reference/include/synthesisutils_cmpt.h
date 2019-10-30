#ifndef _SYNTHESISUTILS_XML_SYNTHESISUTILS_CMPT_
#define _SYNTHESISUTILS_XML_SYNTHESISUTILS_CMPT_
/******************** generated by xml-casa (v2) from synthesisutils.xml ************
********************* 446f8d72a8d4a8d7938e884151ff05b6 *****************************/

#include <vector>
#include <string>
#include <complex>
#include <stdcasa/record.h>
#include <casaswig_types.h>
#include <casa/Exceptions/Error.h>
#include <synthesisutils_forward.h>

#include <synthesisimstore_cmpt.h>


using namespace std;

namespace casac {

  class  synthesisutils  {
    public:

      synthesisutils();
      record* contdatapartition(const record& _selpars=initialize_record(""""), int _npart=int(1));
      record* cubedatapartition(const record& _selpars=initialize_record(""""), int _npart=int(1), const variant& _fstart=variant( ), const variant& _fend=variant( ), const string& _frame=string("LSRK"));
      record* cubeimagepartition(const record& _impars=initialize_record(""""), int _npart=int(1));
      record* cubedataimagepartition(const record& _selpars=initialize_record(""""), const record& _incsys=initialize_record(""""), int _npart=int(1), int _nchannel=int(1));
      record* checkselectionparams(const record& _selpars=initialize_record(""""));
      record* checkimageparams(const record& _impars=initialize_record(""""));
      record* checkgridparams(const record& _gridpars=initialize_record(""""));
      record* updateimpars(const record& _impars=initialize_record(""""));
      int getOptimumSize(int _size=int(100));
      bool done();

        ~synthesisutils( );

    private:

#include <synthesisutils_private.h>


      // --- declarations of static parameter defaults ---
    public:

  };

}

#endif