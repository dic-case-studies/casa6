##################### generated by xml-casa (v2) from calanalysis.xml ###############
##################### 6027fdcee2daf9c77e67911b3056537d ##############################
from __future__ import absolute_import 
from .__casac__ import calanalysis as _calanalysis
from .platform import str_encode as _str_encode
from .platform import str_decode as _str_decode
from .typecheck import validator as _pc
from .coercetype import coerce as _coerce


class calanalysis:
    ### self
    def __init__(self, *args, **kwargs):
        """Construct a calibration analysis tool.
        """
        self._swigobj = kwargs.get('swig_object',None)
        if self._swigobj is None:
            self._swigobj = _calanalysis()

    def open(self, caltable=''):
        """This member function opens a calibration table.
        """
        schema = {'caltable': {'type': 'cReqPath', 'coerce': _coerce.expand_path}}
        doc = {'caltable': caltable}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _open_result = self._swigobj.open(_pc.document['caltable'])
        return _open_result

    def close(self):
        """This member function closes a calibration table.
        """
        _close_result = self._swigobj.close()
        return _close_result

    def calname(self):
        """This member function returns calibration table name.
        """
        _calname_result = self._swigobj.calname()
        return _calname_result

    def msname(self):
        """This member function returns the name of the MS that created this calibration
        table.
        """
        _msname_result = self._swigobj.msname()
        return _msname_result

    def viscal(self):
        """This member function returns the type of calibration table ('B', 'G', 'T',
        etc.).
        """
        _viscal_result = self._swigobj.viscal()
        return _viscal_result

    def partype(self):
        """This member function returns the parameter column type in the calibration table
        ('Complex' or 'Float').
        """
        _partype_result = self._swigobj.partype()
        return _partype_result

    def polbasis(self):
        """This member function returns the polarization basis in the calibration table
        ('L' for linear or 'C' for circular).
        """
        _polbasis_result = self._swigobj.polbasis()
        return _polbasis_result

    def numfield(self):
        """This member function returns the number of fields in the calibration table.
        """
        _numfield_result = self._swigobj.numfield()
        return _numfield_result

    def field(self, name=True):
        """This member function returns the fields in the calibration table.
        """
        schema = {'name': {'type': 'cBool'}}
        doc = {'name': name}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _field_result = self._swigobj.field(_pc.document['name'])
        return _field_result

    def numantenna(self):
        """This member function returns the number of antennas in the calibration table.
        """
        _numantenna_result = self._swigobj.numantenna()
        return _numantenna_result

    def numantenna1(self):
        """This member function returns the number of antenna 1s in the calibration table.
        """
        _numantenna1_result = self._swigobj.numantenna1()
        return _numantenna1_result

    def numantenna2(self):
        """This member function returns the number of antenna 2s in the calibration table.
        """
        _numantenna2_result = self._swigobj.numantenna2()
        return _numantenna2_result

    def antenna(self, name=True):
        """This member function returns the antennas in the calibration table.
        """
        schema = {'name': {'type': 'cBool'}}
        doc = {'name': name}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _antenna_result = self._swigobj.antenna(_pc.document['name'])
        return _antenna_result

    def antenna1(self, name=True):
        """This member function returns the antenna 1s in the calibration table.
        """
        schema = {'name': {'type': 'cBool'}}
        doc = {'name': name}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _antenna1_result = self._swigobj.antenna1(_pc.document['name'])
        return _antenna1_result

    def antenna2(self, name=True):
        """This member function returns the antenna 2s in the calibration table.
        """
        schema = {'name': {'type': 'cBool'}}
        doc = {'name': name}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _antenna2_result = self._swigobj.antenna2(_pc.document['name'])
        return _antenna2_result

    def numfeed(self):
        """This member function returns the number of feeds in the calibration table.
        """
        _numfeed_result = self._swigobj.numfeed()
        return _numfeed_result

    def feed(self):
        """This member function returns the feeds in the calibration table.
        """
        _feed_result = self._swigobj.feed()
        return _feed_result

    def numtime(self):
        """This member function returns the number of times in the calibration table.
        """
        _numtime_result = self._swigobj.numtime()
        return _numtime_result

    def time(self):
        """This member function returns the times (in MJD seconds) in the calibration
        table.
        """
        _time_result = self._swigobj.time()
        return _time_result

    def numspw(self):
        """This member function returns the number of spectral windows in the calibration
        table.
        """
        _numspw_result = self._swigobj.numspw()
        return _numspw_result

    def spw(self, name=True):
        """This member function returns the spectral windows in the calibration table.
        """
        schema = {'name': {'type': 'cBool'}}
        doc = {'name': name}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _spw_result = self._swigobj.spw(_pc.document['name'])
        return _spw_result

    def numchannel(self):
        """This member function returns the number of channels per spectral window in the
        calibration table.
        """
        _numchannel_result = self._swigobj.numchannel()
        return _numchannel_result

    def freq(self):
        """This member function returns the frequencies per spectral window in the
        calibration table.
        """
        _freq_result = self._swigobj.freq()
        return _freq_result

    def get(self, field=[ ], antenna=[ ], timerange=[ ], spw=[ ], feed=[ ], axis='TIME', ap='AMPLITUDE', norm=False, unwrap=False, jumpmax=float(0.0)):
        """This member function returns the calibration data.
        """
        schema = {'field': {'type': 'cVariant'}, 'antenna': {'type': 'cVariant'}, 'timerange': {'type': 'cVariant'}, 'spw': {'type': 'cVariant'}, 'feed': {'type': 'cVariant'}, 'axis': {'type': 'cStr'}, 'ap': {'type': 'cStr'}, 'norm': {'type': 'cBool'}, 'unwrap': {'type': 'cBool'}, 'jumpmax': {'type': 'cFloat', 'coerce': _coerce.to_float}}
        doc = {'field': field, 'antenna': antenna, 'timerange': timerange, 'spw': spw, 'feed': feed, 'axis': axis, 'ap': ap, 'norm': norm, 'unwrap': unwrap, 'jumpmax': jumpmax}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _get_result = self._swigobj.get(_pc.document['field'], _pc.document['antenna'], _pc.document['timerange'], _pc.document['spw'], _pc.document['feed'], _str_encode(_pc.document['axis']), _str_encode(_pc.document['ap']), _pc.document['norm'], _pc.document['unwrap'], _pc.document['jumpmax'])
        return _get_result

    def fit(self, field=[ ], antenna=[ ], timerange=[ ], spw=[ ], feed=[ ], axis='TIME', ap='AMPLITUDE', norm=False, unwrap=False, jumpmax=float(0.0), order='AVERAGE', type='LSQ', weight=False):
        """This member function returns the calibration data and fits along the
        non-iteration axis.
        """
        schema = {'field': {'type': 'cVariant'}, 'antenna': {'type': 'cVariant'}, 'timerange': {'type': 'cVariant'}, 'spw': {'type': 'cVariant'}, 'feed': {'type': 'cVariant'}, 'axis': {'type': 'cStr'}, 'ap': {'type': 'cStr'}, 'norm': {'type': 'cBool'}, 'unwrap': {'type': 'cBool'}, 'jumpmax': {'type': 'cFloat', 'coerce': _coerce.to_float}, 'order': {'type': 'cStr'}, 'type': {'type': 'cStr'}, 'weight': {'type': 'cBool'}}
        doc = {'field': field, 'antenna': antenna, 'timerange': timerange, 'spw': spw, 'feed': feed, 'axis': axis, 'ap': ap, 'norm': norm, 'unwrap': unwrap, 'jumpmax': jumpmax, 'order': order, 'type': type, 'weight': weight}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _fit_result = self._swigobj.fit(_pc.document['field'], _pc.document['antenna'], _pc.document['timerange'], _pc.document['spw'], _pc.document['feed'], _str_encode(_pc.document['axis']), _str_encode(_pc.document['ap']), _pc.document['norm'], _pc.document['unwrap'], _pc.document['jumpmax'], _str_encode(_pc.document['order']), _str_encode(_pc.document['type']), _pc.document['weight'])
        return _fit_result

