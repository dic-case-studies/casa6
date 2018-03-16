#ifndef _MS_XML_MS_CMPT_
#define _MS_XML_MS_CMPT_
/******************** generated by xml-casa (v2) from ms.xml ************************
********************* ec64e6cb974b440d59e96f53eb826a52 *****************************/

#include <vector>
#include <string>
#include <complex>
#include <stdcasa/record.h>
#include <casaswig_types.h>
#include <casa/Exceptions/Error.h>
#include <ms_forward.h>

#include <table_cmpt.h>

#include <msmetadata_cmpt.h>


using namespace std;

namespace casac {

  class  ms  {
    public:

      ms();
      bool open(const string& _thems=string(""), bool _nomodify=bool(true), bool _lock=bool(false), bool _check=bool(false));
      bool reset();
      bool close();
      bool done();
      string name();
      bool iswritable();
      int nrow(bool _selected=bool(false));
      record* getdata(const std::vector<std::string>& _items=std::vector<std::string>({}), bool _ifraxis=bool(false), int _ifraxisgap=int(0), int _increment=int(1), bool _average=bool(false));
      bool putdata(const record& _items=initialize_record(""""));
      bool fromfits(const string& _msfile=string(""), const string& _fitsfile=string(""), bool _nomodify=bool(true), bool _lock=bool(false), int _obstype=int(0), const string& _host=string(""), bool _forcenewserver=bool(false), const string& _antnamescheme=string("old"));
      bool fromfitsidi(const string& _msfile=string(""), const string& _fitsfile=string(""), bool _nomodify=bool(true), bool _lock=bool(false), int _obstype=int(0));
      bool tofits(const string& _fitsfile=string(""), const string& _column=string("corrected"), const variant& _field=variant( ), const variant& _spw=variant( ), const variant& _baseline=variant( ), const string& _time=string(""), const variant& _scan=variant( ), const variant& _uvrange=variant( ), const string& _taql=string(""), bool _writesyscal=bool(false), bool _multisource=bool(false), bool _combinespw=bool(false), bool _writestation=bool(false), bool _padwithflags=bool(false), bool _overwrite=bool(false));
      bool listfits(const string& _fitsfile=string(""));
      string asdmref(const string& _abspath=string(""));
      bool concatenate(const string& _msfile=string(""), const variant& _freqtol=variant( ), const variant& _dirtol=variant( ), float _weightscale=float(1.), int _handling=int(0), const string& _destmsfile=string(""), bool _respectname=bool(false));
      bool testconcatenate(const string& _msfile=string(""), const variant& _freqtol=variant( ), const variant& _dirtol=variant( ), bool _respectname=bool(false));
      bool virtconcatenate(const string& _msfile=string(""), const string& _auxfilename=string(""), const variant& _freqtol=variant( ), const variant& _dirtol=variant( ), float _weightscale=float(1.), bool _respectname=bool(true));
      bool createmultims(const string& _outputTableName=string(""), const std::vector<std::string>& _tables=std::vector<std::string>({}), const std::vector<std::string>& _subtables=std::vector<std::string>({}), bool _nomodify=bool(true), bool _lock=bool(false), bool _copysubtables=bool(false), const std::vector<std::string>& _omitsubtables=std::vector<std::string>({}));
      bool ismultims();
      bool split(const string& _outputms=string(""), const variant& _field=variant( ), const variant& _spw=variant( ), const std::vector<int>& _step=std::vector<int>({1}), const variant& _baseline=variant( ), const variant& _timebin=variant( ), const string& _time=string(""), const variant& _scan=variant( ), const variant& _uvrange=variant( ), const string& _taql=string(""), const string& _whichcol=string("DATA"), const variant& _tileshape=variant( ), const variant& _subarray=variant( ), const string& _combine=string(""), const string& _correlation=string(""), const string& _intent=string(""), const string& _obs=string(""));
      bool partition(const string& _outputms=string(""), const variant& _field=variant( ), const variant& _spw=variant( ), const variant& _baseline=variant( ), const variant& _timebin=variant( ), const string& _time=string(""), const variant& _scan=variant( ), const variant& _uvrange=variant( ), const string& _taql=string(""), const string& _whichcol=string("DATA"), const variant& _tileshape=variant( ), const variant& _subarray=variant( ), const string& _combine=string(""), const string& _intent=string(""), const string& _obs=string(""));
      record* summary(bool _verbose=bool(false), const string& _listfile=string(""), bool _listunfl=bool(false), double _cachesize=double(50), bool _overwrite=bool(false), bool _wantreturn=bool(true));
      record* getscansummary();
      record* getspectralwindowinfo();
      std::vector<std::string> getreferencedtables();
      variant* getfielddirmeas(const string& _dircolname=string("PHASE_DIR"), int _fieldid=int(0), double _time=double(0), const string& _format=string("measure"));
      bool listhistory();
      bool writehistory(const string& _message=string(""), const string& _parms=string(""), const string& _origin=string("MSHistoryHandler::addMessage()"), const string& _msname=string(""), const string& _app=string("ms"));
      record* statistics(const string& _column=string(""), const string& _complex_value=string(""), bool _useflags=bool(true), bool _useweights=bool(false), const string& _spw=string(""), const string& _field=string(""), const string& _baseline=string(""), const string& _uvrange=string(""), const string& _time=string(""), const string& _correlation=string(""), const string& _scan=string(""), const string& _intent=string(""), const string& _array=string(""), const string& _obs=string(""), const string& _reportingaxes=string(""), bool _timeaverage=bool(false), const string& _timebin=string("0s"), const string& _timespan=string(""), double _maxuvwdistance=double(0.0));
      record* statisticsold(const string& _column=string(""), const string& _complex_value=string(""), bool _useflags=bool(true), const string& _spw=string(""), const string& _field=string(""), const string& _baseline=string(""), const string& _uvrange=string(""), const string& _time=string(""), const string& _correlation=string(""), const string& _scan=string(""), const string& _array=string(""), const string& _obs=string(""));
      record* range(const std::vector<std::string>& _items=std::vector<std::string>({}), bool _useflags=bool(true), int _blocksize=int(10));
      bool lister(const string& _options=string(""), const string& _datacolumn=string("data"), const string& _field=string(""), const string& _spw=string(""), const string& _antenna=string(""), const string& _timerange=string(""), const string& _correlation=string(""), const string& _scan=string(""), const string& _feed=string(""), const string& _array=string(""), const string& _observation=string(""), const string& _uvrange=string(""), const string& _average=string(""), bool _showflags=bool(false), const string& _msselect=string(""), int _pagerows=int(50), const string& _listfile=string(""));
      casac::msmetadata* metadata(float _cachesize=float(50));
      bool msselect(const record& _items=initialize_record(""""), bool _onlyparse=bool(false));
      record* msselectedindices();
      record* msseltoindex(const string& _vis=string(""), const variant& _spw=variant( ), const variant& _field=variant( ), const variant& _baseline=variant( ), const variant& _time=variant( ), const variant& _scan=variant( ), const variant& _uvrange=variant( ), const variant& _observation=variant( ), const variant& _polarization=variant( ), const string& _taql=string(""));
      bool selectinit(int _datadescid=int(0), bool _reset=bool(false));
      bool select(const record& _items=initialize_record(""""));
      bool selecttaql(const string& _msselect=string(""));
      bool selectchannel(int _nchan=int(1), int _start=int(0), int _width=int(1), int _inc=int(1));
      bool selectpolarization(const std::vector<std::string>& _wantedpol=std::vector<std::string>({}));
      record* statwt2(const string& _combine=string(""), const variant& _timebin=variant( ), bool _slidetimebin=bool(false), const variant& _chanbin=variant( ), int _minsamp=int(2), const string& _statalg=string("classic"), double _fence=double(-1), const string& _center=string("mean"), bool _lside=bool(true), double _zscore=double(-1), int _maxiter=int(-1), const string& _excludechans=string(""), const std::vector<double>& _wtrange=std::vector<double>({}), bool _preview=bool(false), const string& _datacolumn=string("corrected"));
      bool statwt(bool _dorms=bool(false), bool _byantenna=bool(true), bool _sepacs=bool(true), const variant& _fitspw=variant( ), const variant& _fitcorr=variant( ), const string& _combine=string(""), const variant& _timebin=variant( ), int _minsamp=int(3), const variant& _field=variant( ), const variant& _spw=variant( ), const variant& _antenna=variant( ), const string& _timerange=string(""), const variant& _scan=variant( ), const string& _intent=string(""), const variant& _array=variant( ), const string& _correlation=string(""), const string& _obs=string(""), const string& _datacolumn=string("corrected_data"));
      bool regridspw(const string& _outframe=string("LSRK"), const string& _mode=string("chan"), double _restfreq=double(-3E30), const string& _interpolation=string("LINEAR"), double _start=double(-3E30), double _center=double(-3E30), double _bandwidth=double(-1.), double _chanwidth=double(-1.), bool _hanning=bool(true));
      bool cvel(const string& _mode=string("channel"), int _nchan=int(-1), const variant& _start=variant( ), const variant& _width=variant( ), const string& _interp=string("linear"), const variant& _phasec=variant( ), const variant& _restfreq=variant( ), const string& _outframe=string(""), const string& _veltype=string("radio"), bool _hanning=bool(true));
      bool hanningsmooth(const string& _datacolumn=string("corrected"));
      std::vector<double> cvelfreqs(const std::vector<int>& _spwids=std::vector<int>({0}), const std::vector<int>& _fieldids=std::vector<int>({0}), const string& _obstime=string(""), const string& _mode=string("channel"), int _nchan=int(-1), const variant& _start=variant( ), const variant& _width=variant( ), const variant& _phasec=variant( ), const variant& _restfreq=variant( ), const string& _outframe=string(""), const string& _veltype=string("radio"), bool _verbose=bool(true));
      bool contsub(const string& _outputms=string(""), const variant& _fitspw=variant( ), int _fitorder=int(1), const string& _combine=string(""), const variant& _spw=variant( ), const variant& _unionspw=variant( ), const variant& _field=variant( ), const variant& _scan=variant( ), const string& _intent=string(""), const string& _correlation=string(""), const string& _obs=string(""), const string& _whichcol=string("CORRECTED_DATA"));
      bool continuumsub(const variant& _field=variant( ), const variant& _fitspw=variant( ), const variant& _spw=variant( ), const variant& _solint=variant( ), int _fitorder=int(0), const string& _mode=string("subtract"));
      bool uvsub(bool _reverse=bool(false));
      bool addephemeris(int _id=int(-1), const string& _ephemerisname=string(""), const string& _comment=string(""), const variant& _field=variant( ));
      bool timesort(const string& _newmsname=string(""));
      bool sort(const string& _newmsname=string(""), const std::vector<std::string>& _columns=std::vector<std::string>({}));
      bool iterinit(const std::vector<std::string>& _columns=std::vector<std::string>({}), double _interval=double(0.0), int _maxrows=int(0), bool _adddefaultsortcolumns=bool(true));
      bool iterorigin();
      bool iternext();
      bool iterend();
      record* ngetdata(const std::vector<std::string>& _items=std::vector<std::string>({}), bool _ifraxis=bool(false), int _ifraxisgap=int(0), int _increment=int(1), bool _average=bool(false));
      bool niterinit(const std::vector<std::string>& _columns=std::vector<std::string>({}), double _interval=double(0.0), int _maxrows=int(0), bool _adddefaultsortcolumns=bool(true));
      bool niterorigin();
      bool niternext();
      bool niterend();
      int nrowold(bool _selected=bool(false));
      record* rangeold(const std::vector<std::string>& _items=std::vector<std::string>({}), bool _useflags=bool(true), int _blocksize=int(10));
      bool selectinitold(int _datadescid=int(0), bool _reset=bool(false));
      bool selectold(const record& _items=initialize_record(""""));
      bool selecttaqlold(const string& _msselect=string(""));
      bool selectchannelold(int _nchan=int(1), int _start=int(0), int _width=int(1), int _inc=int(1));
      bool selectpolarizationold(const std::vector<std::string>& _wantedpol=std::vector<std::string>({}));
      record* getdataold(const std::vector<std::string>& _items=std::vector<std::string>({}), bool _ifraxis=bool(false), int _ifraxisgap=int(0), int _increment=int(1), bool _average=bool(false));
      bool putdataold(const record& _items=initialize_record(""""));
      bool iterinitold(const std::vector<std::string>& _columns=std::vector<std::string>({}), double _interval=double(0.0), int _maxrows=int(0), bool _adddefaultsortcolumns=bool(true));
      bool iteroriginold();
      bool iternextold();
      bool iterendold();
      bool continuumsubold(const variant& _field=variant( ), const variant& _fitspw=variant( ), const variant& _spw=variant( ), const variant& _solint=variant( ), int _fitorder=int(0), const string& _mode=string("subtract"));

        ~ms( );

    private:

#include <ms_private.h>


      // --- declarations of static parameter defaults ---
    public:

  };

}

#endif
