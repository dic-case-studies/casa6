#ifndef _IMAGEMETADATA_XML_IMAGEMETADATA_CMPT_
#define _IMAGEMETADATA_XML_IMAGEMETADATA_CMPT_
/******************** generated by xml-casa (v2) from imagemetadata.xml *************
********************* a6949a9b048f9c174974cdafb7d2b5c3 *****************************/

#include <vector>
#include <string>
#include <complex>
#include <stdcasa/record.h>
#include <casaswig_types.h>
#include <casa/Exceptions/Error.h>
#include <imagemetadata_forward.h>


using namespace std;

namespace casac {

  class  imagemetadata  {
    public:

      imagemetadata();
      bool add(const string& _key=string(""), const variant& _value=variant( ));
      bool close();
      bool done();
      variant* get(const string& _key=string(""));
      record* list(bool _verbose=bool(true));
      bool open(const string& _infile=string(""""));
      bool remove(const string& _key=string(""), const variant& _value=variant( ));
      bool set(const string& _key=string(""), const variant& _value=variant( ));

        ~imagemetadata( );

    private:

#include <imagemetadata_private.h>


      // --- declarations of static parameter defaults ---
    public:

  };

}

#endif
