#ifndef _UTILS_XML_UTILS_CMPT_
#define _UTILS_XML_UTILS_CMPT_
/******************** generated by xml-casa (v2) from utils.xml *********************
********************* 7d4ba21f55d1669ed174a0de42d6e3c1 *****************************/

#include <vector>
#include <string>
#include <complex>
#include <stdcasa/record.h>
#include <casaswig_types.h>
#include <casa/Exceptions/Error.h>
#include <utils_forward.h>


using namespace std;

namespace casac {

  class  utils  {
    public:

      utils();
      bool verify(const record& _input=initialize_record(""""), const variant& _xmldescriptor=variant( ), bool _throwexecpt=bool(false));
      bool setconstraints(const variant& _xmldescriptor=variant( ));
      bool verifyparam(const record& _param=initialize_record(""""));
      variant* expandparam(const string& _name=string(""), const variant& _value=variant( ));
      record* torecord(const string& _input=string(""));
      string toxml(const record& _input=initialize_record(""""), bool _asfile=bool(false), const string& _filename=string("recordas.xml"));
      string getrc(const string& _rcvar=string(""));
      bool removetable(const std::vector<std::string>& _tablenames=std::vector<std::string>({}));
      record* tableinfo(const string& _tablename=string(""));
      std::vector<std::string> lockedtables();
      record* hostinfo();
      string c_exception();
      void c_exception_clear();
      string _crash_reporter_initialize(const string& _crashDumpDirectory=string(""), const string& _crashDumpPosterApplication=string(""), const string& _crashPostingUrl=string(""), const string& _logFile=string(""));
      bool _trigger_segfault(int _faultType=int(0));
      bool initialize(const std::vector<std::string>& _default_path=std::vector<std::string>({}));
      std::vector<std::string> defaultpath();
      bool setpath(const std::vector<std::string>& _dirs=std::vector<std::string>({}));
      std::vector<std::string> getpath();
      void clearpath();
      string resolve(const string& _path=string(""));
      record* registry();
      record* services();
      void shutdown();
      std::vector<int> version();
      string version_desc();
      string version_info();
      string version_string();
      bool compare_version(const string& _comparitor=string(""), const std::vector<int>& _vec=std::vector<int>({}));

        ~utils( );

    private:

#include <utils_private.h>


      // --- declarations of static parameter defaults ---
    public:

  };

}

#endif
