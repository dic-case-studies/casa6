#ifndef _REGIONMANAGER_XML_REGIONMANAGER_CMPT_
#define _REGIONMANAGER_XML_REGIONMANAGER_CMPT_
/******************** generated by xml-casa (v2) from regionmanager.xml *************
********************* 4ff09691e40ce125872be37332fdf213 *****************************/

#include <vector>
#include <string>
#include <complex>
#include <stdcasa/record.h>
#include <casaswig_types.h>
#include <casa/Exceptions/Error.h>
#include <regionmanager_forward.h>


using namespace std;

namespace casac {

  class  regionmanager  {
    public:

      regionmanager();
      string absreltype(int _absrelvalue=int(0));
      record* box(const std::vector<double>& _blc=std::vector<double>({0}), const std::vector<double>& _trc=std::vector<double>({-1}), const std::vector<double>& _inc=std::vector<double>({1}), const string& _absrel=string("abs"), bool _frac=bool(false), const string& _comment=string(""));
      record* frombcs(const record& _csys=initialize_record(""""), const std::vector<int>& _shape=std::vector<int>({0}), const string& _box=string(""), const string& _chans=string(""), const string& _stokes=string(""), const string& _stokescontrol=string("a"), const variant& _region=variant( ));
      record* complement(const variant& _region=variant( ), const string& _comment=string(""));
      record* concatenation(const variant& _box=variant( ), const variant& _regions=variant( ), const string& _comment=string(""));
      bool deletefromtable(const string& _tablename=string(""), const string& _regionname=string(""));
      record* difference(const record& _region1=initialize_record(""""), const record& _region2=initialize_record(""""), const string& _comment=string(""));
      bool done();
      std::vector<int> selectedchannels(const string& _specification=string(""), const std::vector<int>& _shape=std::vector<int>({0}));
      record* fromtextfile(const string& _filename=string(""), const std::vector<int>& _shape=std::vector<int>({0}), const record& _csys=initialize_record(""""));
      record* fromtext(const string& _text=string(""), const std::vector<int>& _shape=std::vector<int>({1}), const record& _csys=initialize_record(""""));
      record* fromfiletorecord(const string& _filename=string(""), bool _verbose=bool(true), const string& _regionname=string(""));
      bool tofile(const string& _filename=string(""), const record& _region=initialize_record(""""));
      string fromrecordtotable(const string& _tablename=string(""), const variant& _regionname=variant( ), const record& _regionrec=initialize_record(""""), bool _asmask=bool(false), bool _verbose=bool(true));
      record* fromtabletorecord(const string& _tablename=string(""), const variant& _regionname=variant( ), bool _verbose=bool(true));
      record* intersection(const variant& _regions=variant( ), const string& _comment=string(""));
      bool ispixelregion(const record& _region=initialize_record(""""));
      bool isworldregion(const record& _region=initialize_record(""""));
      std::vector<std::string> namesintable(const string& _tablename=string(""));
      bool setcoordinates(const record& _csys=initialize_record(""""));
      record* makeunion(const variant& _regions=variant( ), const string& _comment=string(""));
      record* wbox(const variant& _blc=variant( ), const variant& _trc=variant( ), const std::vector<int>& _pixelaxes=std::vector<int>({-1}), const record& _csys=initialize_record(""""), const string& _absrel=string("abs"), const string& _comment=string(""));
      record* wpolygon(const variant& _x=variant( ), const variant& _y=variant( ), const std::vector<int>& _pixelaxes=std::vector<int>({-1}), const record& _csys=initialize_record(""""), const string& _absrel=string("abs"), const string& _comment=string(""));

        ~regionmanager( );

    private:

#include <regionmanager_private.h>


      // --- declarations of static parameter defaults ---
    public:

  };

}

#endif
