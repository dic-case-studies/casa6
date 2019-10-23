#ifndef _QUANTA_XML_QUANTA_CMPT_
#define _QUANTA_XML_QUANTA_CMPT_
/******************** generated by xml-casa (v2) from quanta.xml ********************
********************* 78b1178759f40bb6be6a580df9a7dfa8 *****************************/

#include <vector>
#include <string>
#include <complex>
#include <stdcasa/record.h>
#include <casaswig_types.h>
#include <casa/Exceptions/Error.h>
#include <quanta_forward.h>


using namespace std;

namespace casac {

  class  quanta  {
    public:

      quanta();
      record* convertfreq(const variant& _v=variant( ), const string& _outunit=string("Hz"));
      record* convertdop(const variant& _v=variant( ), const string& _outunit=string("km/s"));
      record* quantity(const variant& _v=variant( ), const string& _unitname=string(""));
      std::vector<double> getvalue(const variant& _v=variant( ));
      string getunit(const variant& _v=variant( ));
      record* canonical(const variant& _v=variant( ));
      record* canon(const variant& _v=variant( ));
      record* convert(const variant& _v=variant( ), const variant& _outunit=variant( ));
      bool define(const string& _name=string(""), const variant& _v=variant( ));
      string map(const string& _v=string("all"));
      record* maprec(const string& _v=string("all"));
      bool fits();
      std::vector<std::string> angle(const variant& _v=variant( ), int _prec=int(0), const std::vector<std::string>& _form=std::vector<std::string>({}), bool _showform=bool(false));
      std::vector<std::string> time(const variant& _v=variant( ), int _prec=int(0), const std::vector<std::string>& _form=std::vector<std::string>({}), bool _showform=bool(false));
      record* add(const variant& _v=variant( ), const variant& _a=variant( ));
      record* sub(const variant& _v=variant( ), const variant& _a=variant( ));
      record* mul(const variant& _v=variant( ), const variant& _a=variant( ));
      record* div(const variant& _v=variant( ), const variant& _a=variant( ));
      record* neg(const variant& _v=variant( ));
      record* norm(const variant& _v=variant( ), double _a=double(-0.5));
      bool le(const variant& _v=variant( ), const variant& _a=variant( ));
      bool lt(const variant& _v=variant( ), const variant& _a=variant( ));
      bool eq(const variant& _v=variant( ), const variant& _a=variant( ));
      bool ne(const variant& _v=variant( ), const variant& _a=variant( ));
      bool gt(const variant& _v=variant( ), const variant& _a=variant( ));
      bool ge(const variant& _v=variant( ), const variant& _a=variant( ));
      record* sin(const variant& _v=variant( ));
      record* cos(const variant& _v=variant( ));
      record* tan(const variant& _v=variant( ));
      record* asin(const variant& _v=variant( ));
      record* acos(const variant& _v=variant( ));
      record* atan(const variant& _v=variant( ));
      record* atan2(const variant& _v=variant( ), const variant& _a=variant( ));
      record* abs(const variant& _v=variant( ));
      record* ceil(const variant& _v=variant( ));
      record* floor(const variant& _v=variant( ));
      record* log(const variant& _v=variant( ));
      record* log10(const variant& _v=variant( ));
      record* exp(const variant& _v=variant( ));
      record* sqrt(const variant& _v=variant( ));
      bool compare(const variant& _v=variant( ), const variant& _a=variant( ));
      bool check(const string& _v=string(""));
      bool checkfreq(const variant& _cm=variant( ));
      record* pow(const variant& _v=variant( ), int _a=int(1));
      record* constants(const string& _v=string("pi"));
      bool isangle(const variant& _v=variant( ));
      record* totime(const variant& _v=variant( ));
      record* toangle(const variant& _v=variant( ));
      record* splitdate(const variant& _v=variant( ));
      string tos(const variant& _v=variant( ), int _prec=int(9));
      string type();
      bool done(bool _kill=bool(false));
      record* unit(const variant& _v=variant( ), const string& _unitname=string(""));
      bool isquantity(const variant& _v=variant( ));
      bool setformat(const string& _t=string(""), const string& _v=string("F"));
      string getformat(const string& _t=string(""));
      string formxxx(const variant& _v=variant( ), const string& _format=string("dms"), int _prec=int(2));

        ~quanta( );

    private:

#include <quanta_private.h>


      // --- declarations of static parameter defaults ---
    public:

  };

}

#endif
