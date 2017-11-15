##################### generated by xml-casa (v2) from calanalysis.xml ###############
##################### 6027fdcee2daf9c77e67911b3056537d ##############################
from __future__ import absolute_import 
from .__casac__ import calanalysis as _calanalysis
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
        return self._swigobj.open(_pc.document['caltable'])

    def close(self):
        """This member function closes a calibration table.
        """
        return self._swigobj.close()

    def calname(self):
        """This member function returns calibration table name.
        """
        return self._swigobj.calname()

    def msname(self):
        """This member function returns the name of the MS that created this calibration
        table.
        """
        return self._swigobj.msname()

    def viscal(self):
        """This member function returns the type of calibration table ('B', 'G', 'T',
        etc.).
        """
        return self._swigobj.viscal()

    def partype(self):
        """This member function returns the parameter column type in the calibration table
        ('Complex' or 'Float').
        """
        return self._swigobj.partype()

    def polbasis(self):
        """This member function returns the polarization basis in the calibration table
        ('L' for linear or 'C' for circular).
        """
        return self._swigobj.polbasis()

    def numfield(self):
        """This member function returns the number of fields in the calibration table.
        """
        return self._swigobj.numfield()

    def field(self, name=True):
        """This member function returns the fields in the calibration table.
        """
        schema = {'name': {'type': 'cBool'}}
        doc = {'name': name}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.field(_pc.document['name'])

    def numantenna(self):
        """This member function returns the number of antennas in the calibration table.
        """
        return self._swigobj.numantenna()

    def numantenna1(self):
        """This member function returns the number of antenna 1s in the calibration table.
        """
        return self._swigobj.numantenna1()

    def numantenna2(self):
        """This member function returns the number of antenna 2s in the calibration table.
        """
        return self._swigobj.numantenna2()

    def antenna(self, name=True):
        """This member function returns the antennas in the calibration table.
        """
        schema = {'name': {'type': 'cBool'}}
        doc = {'name': name}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.antenna(_pc.document['name'])

    def antenna1(self, name=True):
        """This member function returns the antenna 1s in the calibration table.
        """
        schema = {'name': {'type': 'cBool'}}
        doc = {'name': name}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.antenna1(_pc.document['name'])

    def antenna2(self, name=True):
        """This member function returns the antenna 2s in the calibration table.
        """
        schema = {'name': {'type': 'cBool'}}
        doc = {'name': name}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.antenna2(_pc.document['name'])

    def numfeed(self):
        """This member function returns the number of feeds in the calibration table.
        """
        return self._swigobj.numfeed()

    def feed(self):
        """This member function returns the feeds in the calibration table.
        """
        return self._swigobj.feed()

    def numtime(self):
        """This member function returns the number of times in the calibration table.
        """
        return self._swigobj.numtime()

    def time(self):
        """This member function returns the times (in MJD seconds) in the calibration
        table.
        """
        return self._swigobj.time()

    def numspw(self):
        """This member function returns the number of spectral windows in the calibration
        table.
        """
        return self._swigobj.numspw()

    def spw(self, name=True):
        """This member function returns the spectral windows in the calibration table.
        """
        schema = {'name': {'type': 'cBool'}}
        doc = {'name': name}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.spw(_pc.document['name'])

    def numchannel(self):
        """This member function returns the number of channels per spectral window in the
        calibration table.
        """
        return self._swigobj.numchannel()

    def freq(self):
        """This member function returns the frequencies per spectral window in the
        calibration table.
        """
        return self._swigobj.freq()

    def get(self, field=[ ], antenna=[ ], timerange=[ ], spw=[ ], feed=[ ], axis='TIME', ap='AMPLITUDE', norm=False, unwrap=False, jumpmax=float(0.0)):
        """This member function returns the calibration data.
        """
        schema = {'field': {'type': 'cVariant'}, 'antenna': {'type': 'cVariant'}, 'timerange': {'type': 'cVariant'}, 'spw': {'type': 'cVariant'}, 'feed': {'type': 'cVariant'}, 'axis': {'type': 'cStr'}, 'ap': {'type': 'cStr'}, 'norm': {'type': 'cBool'}, 'unwrap': {'type': 'cBool'}, 'jumpmax': {'type': 'cFloat', 'coerce': _coerce.to_float}}
        doc = {'field': field, 'antenna': antenna, 'timerange': timerange, 'spw': spw, 'feed': feed, 'axis': axis, 'ap': ap, 'norm': norm, 'unwrap': unwrap, 'jumpmax': jumpmax}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.get(_pc.document['field'], _pc.document['antenna'], _pc.document['timerange'], _pc.document['spw'], _pc.document['feed'], _pc.document['axis'], _pc.document['ap'], _pc.document['norm'], _pc.document['unwrap'], _pc.document['jumpmax'])

    def fit(self, field=[ ], antenna=[ ], timerange=[ ], spw=[ ], feed=[ ], axis='TIME', ap='AMPLITUDE', norm=False, unwrap=False, jumpmax=float(0.0), order='AVERAGE', type='LSQ', weight=False):
        """This member function returns the calibration data and fits along the
        non-iteration axis.
        """
        schema = {'field': {'type': 'cVariant'}, 'antenna': {'type': 'cVariant'}, 'timerange': {'type': 'cVariant'}, 'spw': {'type': 'cVariant'}, 'feed': {'type': 'cVariant'}, 'axis': {'type': 'cStr'}, 'ap': {'type': 'cStr'}, 'norm': {'type': 'cBool'}, 'unwrap': {'type': 'cBool'}, 'jumpmax': {'type': 'cFloat', 'coerce': _coerce.to_float}, 'order': {'type': 'cStr'}, 'type': {'type': 'cStr'}, 'weight': {'type': 'cBool'}}
        doc = {'field': field, 'antenna': antenna, 'timerange': timerange, 'spw': spw, 'feed': feed, 'axis': axis, 'ap': ap, 'norm': norm, 'unwrap': unwrap, 'jumpmax': jumpmax, 'order': order, 'type': type, 'weight': weight}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.fit(_pc.document['field'], _pc.document['antenna'], _pc.document['timerange'], _pc.document['spw'], _pc.document['feed'], _pc.document['axis'], _pc.document['ap'], _pc.document['norm'], _pc.document['unwrap'], _pc.document['jumpmax'], _pc.document['order'], _pc.document['type'], _pc.document['weight'])

