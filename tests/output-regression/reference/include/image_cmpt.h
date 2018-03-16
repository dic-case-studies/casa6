#ifndef _IMAGE_XML_IMAGE_CMPT_
#define _IMAGE_XML_IMAGE_CMPT_
/******************** generated by xml-casa (v2) from image.xml *********************
********************* a02e4b3e0b3c4c5f3f8f771db3b7b111 *****************************/

#include <vector>
#include <string>
#include <complex>
#include <stdcasa/record.h>
#include <casaswig_types.h>
#include <casa/Exceptions/Error.h>
#include <image_forward.h>

#include <coordsys_cmpt.h>

#include <componentlist_cmpt.h>


using namespace std;

namespace casac {

  class  image  {
    public:

      image();
      casac::image* newimage(const string& _infile=string(""));
      casac::image* newimagefromfile(const string& _infile=string(""));
      casac::image* imagecalc(const string& _outfile=string(""), const string& _pixels=string(""), bool _overwrite=bool(false), const string& _imagemd=string(""));
      casac::image* collapse(const string& _function=string(""), const variant& _axes=variant( ), const string& _outfile=string(""), const variant& _region=variant( ), const string& _box=string(""), const string& _chans=string(""), const string& _stokes=string(""), const string& _mask=string(""), bool _overwrite=bool(false), bool _stretch=bool(false));
      casac::image* decimate(const string& _outfile=string(""), int _axis=int(0), int _factor=int(1), const string& _method=string("copy"), const variant& _region=variant( ), const string& _mask=string(""), bool _overwrite=bool(false), bool _stretch=bool(false));
      bool dohistory(bool _enable=bool(True));
      casac::image* imageconcat(const string& _outfile=string(""), const variant& _infiles=variant( ), int _axis=int(-1), bool _relax=bool(false), bool _tempclose=bool(true), bool _overwrite=bool(false), bool _reorder=bool(false));
      bool fromarray(const string& _outfile=string(""), const variant& _pixels=variant( ), const record& _csys=initialize_record(""""), bool _linear=bool(false), bool _overwrite=bool(false), bool _log=bool(true));
      bool fromascii(const string& _outfile=string(""), const string& _infile=string(""), const std::vector<int>& _shape=std::vector<int>({-1}), const string& _sep=string(":"), const record& _csys=initialize_record(""""), bool _linear=bool(false), bool _overwrite=bool(false));
      bool fromcomplist(const string& _outfile=string(""), const std::vector<int>& _shape=std::vector<int>({}), const variant& _cl=variant( ), const record& _csys=initialize_record(""""), bool _overwrite=bool(false), bool _log=bool(true), bool _cache=bool(true));
      bool fromfits(const string& _outfile=string(""), const string& _infile=string(""), int _whichrep=int(0), int _whichhdu=int(0), bool _zeroblanks=bool(false), bool _overwrite=bool(false));
      bool fromimage(const string& _outfile=string(""), const string& _infile=string(""), const variant& _region=variant( ), const variant& _mask=variant( ), bool _dropdeg=bool(false), bool _overwrite=bool(false));
      bool fromshape(const string& _outfile=string(""), const std::vector<int>& _shape=std::vector<int>({0}), const record& _csys=initialize_record(""""), bool _linear=bool(false), bool _overwrite=bool(false), bool _log=bool(true), const string& _type=string("f"));
      bool maketestimage(const string& _outfile=string(""), bool _overwrite=bool(false));
      casac::image* deviation(const string& _outfile=string(""), const variant& _region=variant( ), const string& _mask=string(""), bool _overwrite=bool(false), bool _stretch=bool(false), const std::vector<int>& _grid=std::vector<int>({1,1}), const variant& _anchor=variant( ), const variant& _xlength=variant( ), const variant& _ylength=variant( ), const string& _interp=string("cubic"), const string& _stattype=string("sigma"), const string& _statalg=string("classic"), double _zscore=double(-1), int _maxiter=int(-1));
      casac::image* adddegaxes(const string& _outfile=string(""), bool _direction=bool(false), bool _spectral=bool(false), const string& _stokes=string(""), bool _linear=bool(false), bool _tabular=bool(false), bool _overwrite=bool(false), bool _silent=bool(false));
      bool addnoise(const string& _type=string("normal"), const std::vector<double>& _pars=std::vector<double>({0.0,1.0}), const variant& _region=variant( ), bool _zero=bool(false), const std::vector<int>& _seeds=std::vector<int>({}));
      casac::image* convolve(const string& _outfile=string(""), const variant& _kernel=variant( ), double _scale=double(-1.0), const variant& _region=variant( ), const variant& _mask=variant( ), bool _overwrite=bool(false), bool _stretch=bool(false));
      record* boundingbox(const variant& _region=variant( ));
      casac::image* boxcar(const string& _outfile=string(""), const variant& _region=variant( ), const variant& _mask=variant( ), int _axis=int(-1), int _width=int(2), bool _drop=bool(true), const string& _dmethod=string("copy"), bool _overwrite=bool(false), bool _stretch=bool(false));
      string brightnessunit();
      bool calc(const string& _pixels=string(""), bool _verbose=bool(True));
      bool calcmask(const string& _mask=string(""), const string& _name=string(""), bool _asdefault=bool(true));
      bool close();
      casac::image* continuumsub(const string& _outline=string(""), const string& _outcont=string("continuumsub.im"), const variant& _region=variant( ), const std::vector<int>& _channels=std::vector<int>({-1}), const string& _pol=string(""), int _fitorder=int(0), bool _overwrite=bool(false));
      record* convertflux(const variant& _value=variant( ), const variant& _major=variant( ), const variant& _minor=variant( ), const string& _type=string("Gaussian"), bool _topeak=bool(true), int _channel=int(-1), int _polarization=int(-1));
      casac::image* convolve2d(const string& _outfile=string(""), const std::vector<int>& _axes=std::vector<int>({0,1}), const string& _type=string("gaussian"), const variant& _major=variant( ), const variant& _minor=variant( ), const variant& _pa=variant( ), double _scale=double(-1), const variant& _region=variant( ), const variant& _mask=variant( ), bool _overwrite=bool(false), bool _stretch=bool(false), bool _targetres=bool(false), const record& _beam=initialize_record(""""));
      casac::coordsys* coordsys(const std::vector<int>& _axes=std::vector<int>({-1}));
      record* coordmeasures(const std::vector<double>& _pixel=std::vector<double>({-1}), const string& _dframe=string("cl"), const string& _sframe=string("cl"));
      record* decompose(const variant& _region=variant( ), const variant& _mask=variant( ), bool _simple=bool(false), double _threshold=double(-1), int _ncontour=int(11), int _minrange=int(1), int _naxis=int(2), bool _fit=bool(true), double _maxrms=double(-1), int _maxretry=int(-1), int _maxiter=int(256), double _convcriteria=double(0.0001), bool _stretch=bool(false));
      record* deconvolvecomponentlist(const record& _complist=initialize_record(""""), int _channel=int(-1), int _polarization=int(-1));
      record* deconvolvefrombeam(const variant& _source=variant( ), const variant& _beam=variant( ));
      record* beamforconvolvedsize(const variant& _source=variant( ), const variant& _convolved=variant( ));
      record* commonbeam();
      bool remove(bool _done=bool(false), bool _verbose=bool(true));
      bool removefile(const string& _file=string(""));
      bool done(bool _remove=bool(false), bool _verbose=bool(true));
      bool fft(const string& _real=string(""), const string& _imag=string(""), const string& _amp=string(""), const string& _phase=string(""), const std::vector<int>& _axes=std::vector<int>({-1}), const variant& _region=variant( ), const variant& _mask=variant( ), bool _stretch=bool(false), const string& _complex=string(""));
      record* findsources(int _nmax=int(20), double _cutoff=double(0.1), const variant& _region=variant( ), const variant& _mask=variant( ), bool _point=bool(true), int _width=int(5), bool _negfind=bool(false));
      record* fitprofile(const string& _box=string(""), const variant& _region=variant( ), const string& _chans=string(""), const string& _stokes=string(""), int _axis=int(-1), const variant& _mask=variant( ), int _ngauss=int(1), int _poly=int(-1), const string& _estimates=string(""), int _minpts=int(1), bool _multifit=bool(false), const string& _model=string(""), const string& _residual=string(""), const string& _amp=string(""), const string& _amperr=string(""), const string& _center=string(""), const string& _centererr=string(""), const string& _fwhm=string(""), const string& _fwhmerr=string(""), const string& _integral=string(""), const string& _integralerr=string(""), bool _stretch=bool(false), bool _logresults=bool(true), const variant& _pampest=variant( ), const variant& _pcenterest=variant( ), const variant& _pfwhmest=variant( ), const variant& _pfix=variant( ), const variant& _gmncomps=variant( ), const variant& _gmampcon=variant( ), const variant& _gmcentercon=variant( ), const variant& _gmfwhmcon=variant( ), const std::vector<double>& _gmampest=std::vector<double>({0.0}), const std::vector<double>& _gmcenterest=std::vector<double>({0.0}), const std::vector<double>& _gmfwhmest=std::vector<double>({0.0}), const variant& _gmfix=variant( ), const string& _spxtype=string(""), const std::vector<double>& _spxest=std::vector<double>({}), const std::vector<bool>& _spxfix=std::vector<bool>({}), const variant& _div=variant( ), const string& _spxsol=string(""), const string& _spxerr=string(""), const string& _logfile=string(""), bool _append=bool(true), const variant& _pfunc=variant( ), const std::vector<double>& _goodamprange=std::vector<double>({0.0}), const std::vector<double>& _goodcenterrange=std::vector<double>({0.0}), const std::vector<double>& _goodfwhmrange=std::vector<double>({0.0}), const variant& _sigma=variant( ), const string& _outsigma=string(""), const std::vector<int>& _planes=std::vector<int>({}));
      record* fitcomponents(const string& _box=string(""), const variant& _region=variant( ), const variant& _chans=variant( ), const string& _stokes=string(""), const variant& _mask=variant( ), const std::vector<double>& _includepix=std::vector<double>({-1}), const std::vector<double>& _excludepix=std::vector<double>({-1}), const string& _residual=string(""), const string& _model=string(""), const string& _estimates=string(""), const string& _logfile=string(""), bool _append=bool(true), const string& _newestimates=string(""), const string& _complist=string(""), bool _overwrite=bool(false), bool _dooff=bool(false), double _offset=double(0.0), bool _fixoffset=bool(false), bool _stretch=bool(false), const variant& _rms=variant( ), const variant& _noisefwhm=variant( ), const string& _summary=string(""));
      bool fromrecord(const record& _record=initialize_record(""""), const string& _outfile=string(""));
      variant* getchunk(const std::vector<int>& _blc=std::vector<int>({-1}), const std::vector<int>& _trc=std::vector<int>({-1}), const std::vector<int>& _inc=std::vector<int>({1}), const std::vector<int>& _axes=std::vector<int>({-1}), bool _list=bool(false), bool _dropdeg=bool(false), bool _getmask=bool(false));
      variant* getregion(const variant& _region=variant( ), const std::vector<int>& _axes=std::vector<int>({-1}), const variant& _mask=variant( ), bool _list=bool(false), bool _dropdeg=bool(false), bool _getmask=bool(false), bool _stretch=bool(false));
      record* getprofile(int _axis=int(-1), const string& _function=string("mean"), const variant& _region=variant( ), const string& _mask=string(""), const string& _unit=string(""), bool _stretch=bool(false), const string& _spectype=string("default"), const variant& _restfreq=variant( ), const string& _frame=string(""), const string& _logfile=string(""));
      record* getslice(const std::vector<double>& _x=std::vector<double>({}), const std::vector<double>& _y=std::vector<double>({}), const std::vector<int>& _axes=std::vector<int>({0,1}), const std::vector<int>& _coord=std::vector<int>({-1}), int _npts=int(0), const string& _method=string("linear"));
      casac::image* hanning(const string& _outfile=string(""), const variant& _region=variant( ), const variant& _mask=variant( ), int _axis=int(-10), bool _drop=bool(true), bool _overwrite=bool(false), bool _async=bool(false), bool _stretch=bool(false), const string& _dmethod=string("copy"));
      std::vector<bool> haslock();
      record* histograms(const std::vector<int>& _axes=std::vector<int>({-1}), const variant& _region=variant( ), const variant& _mask=variant( ), int _nbins=int(25), const std::vector<double>& _includepix=std::vector<double>({-1}), bool _cumu=bool(false), bool _log=bool(false), bool _stretch=bool(false));
      std::vector<std::string> history(bool _list=bool(true));
      bool insert(const string& _infile=string(""), const variant& _region=variant( ), const std::vector<double>& _locate=std::vector<double>({-1}), bool _verbose=bool(false));
      bool isopen();
      bool ispersistent();
      bool lock(bool _writelock=bool(false), int _nattempts=int(0));
      bool makecomplex(const string& _outfile=string(""), const string& _imag=string(""), const variant& _region=variant( ), bool _overwrite=bool(false));
      std::vector<std::string> maskhandler(const string& _op=string("default"), const std::vector<std::string>& _name=std::vector<std::string>({}));
      record* miscinfo();
      bool modify(const record& _model=initialize_record(""""), const variant& _region=variant( ), const variant& _mask=variant( ), bool _subtract=bool(true), bool _list=bool(true), bool _stretch=bool(false));
      record* maxfit(const variant& _region=variant( ), bool _point=bool(true), int _width=int(5), bool _negfind=bool(false), bool _list=bool(true));
      casac::image* moments(const std::vector<int>& _moments=std::vector<int>({0}), int _axis=int(-10), const variant& _region=variant( ), const variant& _mask=variant( ), const std::vector<std::string>& _method=std::vector<std::string>({}), const std::vector<int>& _smoothaxes=std::vector<int>({-1}), const variant& _smoothtypes=variant( ), const std::vector<double>& _smoothwidths=std::vector<double>({0.0}), const std::vector<double>& _includepix=std::vector<double>({-1}), const std::vector<double>& _excludepix=std::vector<double>({-1}), double _peaksnr=double(3.0), double _stddev=double(0.0), const string& _doppler=string("RADIO"), const string& _outfile=string(""), const string& _smoothout=string(""), bool _overwrite=bool(false), bool _drop=bool(true), bool _stretch=bool(false), bool _async=bool(false));
      string name(bool _strippath=bool(false));
      bool open(const string& _infile=string(""""), bool _cache=bool(true));
      casac::image* pad(const string& _outfile=string(""), int _npixels=int(1), double _value=double(0), bool _padmask=bool(false), bool _overwrite=bool(false), const variant& _region=variant( ), const string& _box=string(""), const string& _chans=string(""), const string& _stokes=string(""), const string& _mask=string(""), bool _stretch=bool(false), bool _wantreturn=bool(true));
      casac::image* crop(const string& _outfile=string(""), const std::vector<int>& _axes=std::vector<int>({}), bool _overwrite=bool(false), const variant& _region=variant( ), const string& _box=string(""), const string& _chans=string(""), const string& _stokes=string(""), const string& _mask=string(""), bool _stretch=bool(false), bool _wantreturn=bool(true));
      record* pixelvalue(const std::vector<int>& _pixel=std::vector<int>({-1}));
      bool putchunk(const variant& _pixels=variant( ), const std::vector<int>& _blc=std::vector<int>({-1}), const std::vector<int>& _inc=std::vector<int>({1}), bool _list=bool(false), bool _locking=bool(true), bool _replicate=bool(false));
      bool putregion(const variant& _pixels=variant( ), const variant& _pixelmask=variant( ), const variant& _region=variant( ), bool _list=bool(false), bool _usemask=bool(true), bool _locking=bool(true), bool _replicate=bool(false));
      casac::image* rebin(const string& _outfile=string(""), const std::vector<int>& _bin=std::vector<int>({}), const variant& _region=variant( ), const variant& _mask=variant( ), bool _dropdeg=bool(false), bool _overwrite=bool(false), bool _async=bool(false), bool _stretch=bool(false), bool _crop=bool(false));
      casac::image* regrid(const string& _outfile=string(""), const std::vector<int>& _shape=std::vector<int>({-1}), const record& _csys=initialize_record(""""), const std::vector<int>& _axes=std::vector<int>({-1}), const variant& _region=variant( ), const variant& _mask=variant( ), const string& _method=string("linear"), int _decimate=int(10), bool _replicate=bool(false), bool _doref=bool(true), bool _dropdeg=bool(false), bool _overwrite=bool(false), bool _force=bool(false), bool _asvelocity=bool(false), bool _async=bool(false), bool _stretch=bool(false));
      casac::image* transpose(const string& _outfile=string(""), const variant& _order=variant( ));
      casac::image* rotate(const string& _outfile=string(""), const std::vector<int>& _shape=std::vector<int>({-1}), const variant& _pa=variant( ), const variant& _region=variant( ), const variant& _mask=variant( ), const string& _method=string("cubic"), int _decimate=int(0), bool _replicate=bool(false), bool _dropdeg=bool(false), bool _overwrite=bool(false), bool _stretch=bool(false));
      bool rotatebeam(const variant& _angle=variant( ));
      bool rename(const string& _name=string(""), bool _overwrite=bool(false));
      bool replacemaskedpixels(const variant& _pixels=variant( ), const variant& _region=variant( ), const variant& _mask=variant( ), bool _update=bool(false), bool _list=bool(false), bool _stretch=bool(false));
      record* beamarea(int _channel=int(-1), int _polarization=int(-1));
      record* restoringbeam(int _channel=int(-1), int _polarization=int(-1));
      casac::image* sepconvolve(const string& _outfile=string(""), const std::vector<int>& _axes=std::vector<int>({-1}), const std::vector<std::string>& _types=std::vector<std::string>({""}), const variant& _widths=variant( ), double _scale=double(-1), const variant& _region=variant( ), const variant& _mask=variant( ), bool _overwrite=bool(false), bool _stretch=bool(false));
      bool set(const variant& _pixels=variant( ), int _pixelmask=int(-1), const variant& _region=variant( ), bool _list=bool(false));
      bool setbrightnessunit(const string& _unit=string(""));
      bool setcoordsys(const record& _csys=initialize_record(""""));
      bool sethistory(const string& _origin=string(""), const std::vector<std::string>& _history=std::vector<std::string>({}));
      bool setmiscinfo(const record& _info=initialize_record(""""));
      std::vector<int> shape();
      bool setrestoringbeam(const variant& _major=variant( ), const variant& _minor=variant( ), const variant& _pa=variant( ), const record& _beam=initialize_record(""""), bool _remove=bool(false), bool _log=bool(true), int _channel=int(-1), int _polarization=int(-1), const string& _imagename=string(""));
      record* statistics(const std::vector<int>& _axes=std::vector<int>({-1}), const variant& _region=variant( ), const variant& _mask=variant( ), const std::vector<double>& _includepix=std::vector<double>({-1}), const std::vector<double>& _excludepix=std::vector<double>({-1}), bool _list=bool(false), bool _force=bool(false), bool _disk=bool(false), bool _robust=bool(false), bool _verbose=bool(false), bool _stretch=bool(false), const string& _logfile=string(""), bool _append=bool(true), const string& _algorithm=string("classic"), double _fence=double(-1), const string& _center=string("mean"), bool _lside=bool(true), double _zscore=double(-1), int _maxiter=int(-1), const string& _clmethod=string("auto"));
      bool twopointcorrelation(const string& _outfile=string(""), const variant& _region=variant( ), const variant& _mask=variant( ), const std::vector<int>& _axes=std::vector<int>({-1}), const string& _method=string("structurefunction"), bool _overwrite=bool(false), bool _stretch=bool(false));
      casac::image* subimage(const string& _outfile=string(""), const variant& _region=variant( ), const variant& _mask=variant( ), bool _dropdeg=bool(false), bool _overwrite=bool(false), bool _list=bool(true), bool _stretch=bool(false), bool _wantreturn=bool(true), const std::vector<int>& _keepaxes=std::vector<int>({}));
      record* summary(const string& _doppler=string("RADIO"), bool _list=bool(true), bool _pixelorder=bool(true), bool _verbose=bool(false));
      bool tofits(const string& _outfile=string(""), bool _velocity=bool(false), bool _optical=bool(true), int _bitpix=int(-32), double _minpix=double(1), double _maxpix=double(-1), const variant& _region=variant( ), const variant& _mask=variant( ), bool _overwrite=bool(false), bool _dropdeg=bool(false), bool _deglast=bool(false), bool _dropstokes=bool(false), bool _stokeslast=bool(true), bool _wavelength=bool(false), bool _airwavelength=bool(false), bool _async=bool(false), bool _stretch=bool(false), bool _history=bool(true));
      bool toASCII(const string& _outfile=string(""), const variant& _region=variant( ), const variant& _mask=variant( ), const string& _sep=string(":"), const string& _format=string("%e"), double _maskvalue=double(-999), bool _overwrite=bool(false), bool _stretch=bool(false));
      record* torecord();
      string type();
      record* topixel(const variant& _value=variant( ));
      record* toworld(const variant& _value=variant( ), const string& _format=string("n"), bool _dovelocity=bool(True));
      bool unlock();
      casac::image* newimagefromarray(const string& _outfile=string(""), const variant& _pixels=variant( ), const record& _csys=initialize_record(""""), bool _linear=bool(false), bool _overwrite=bool(false), bool _log=bool(true));
      casac::image* newimagefromfits(const string& _outfile=string(""), const string& _infile=string(""), int _whichrep=int(0), int _whichhdu=int(0), bool _zeroblanks=bool(false), bool _overwrite=bool(false));
      casac::image* newimagefromimage(const string& _infile=string(""), const string& _outfile=string(""), const variant& _region=variant( ), const variant& _mask=variant( ), bool _dropdeg=bool(false), bool _overwrite=bool(false));
      casac::image* newimagefromshape(const string& _outfile=string(""), const std::vector<int>& _shape=std::vector<int>({0}), const record& _csys=initialize_record(""""), bool _linear=bool(false), bool _overwrite=bool(false), bool _log=bool(true), const string& _type=string("f"));
      casac::image* pbcor(const variant& _pbimage=variant( ), const string& _outfile=string(""), bool _overwrite=bool(false), const string& _box=string(""), const variant& _region=variant( ), const string& _chans=string(""), const string& _stokes=string(""), const string& _mask=string(""), const string& _mode=string("divide"), float _cutoff=float(-1.0), bool _stretch=bool(false));
      casac::image* pv(const string& _outfile=string(""), const variant& _start=variant( ), const variant& _end=variant( ), const variant& _center=variant( ), const variant& _length=variant( ), const variant& _pa=variant( ), const variant& _width=variant( ), const string& _unit=string("arcsec"), bool _overwrite=bool(false), const variant& _region=variant( ), const string& _chans=string(""), const string& _stokes=string(""), const string& _mask=string(""), bool _stretch=bool(false), bool _wantreturn=bool(true));
      variant* makearray(double _v=double(0.0), const std::vector<int>& _shape=std::vector<int>({0}));
      bool isconform(const string& _other=string(""));

        ~image( );

    private:

#include <image_private.h>


      // --- declarations of static parameter defaults ---
    public:

  };

}

#endif
