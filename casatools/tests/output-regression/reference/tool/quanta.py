##################### generated by xml-casa (v2) from quanta.xml ####################
##################### 78b1178759f40bb6be6a580df9a7dfa8 ##############################
from __future__ import absolute_import 
from .__casac__ import quanta as _quanta
from .platform import str_encode as _str_ec
from .platform import str_decode as _str_dc
from .platform import dict_encode as _dict_ec
from .platform import dict_decode as _dict_dc
from .platform import dict_encode as _quant_ec
from .platform import dict_decode as _quant_dc
from .platform import encode as _any_ec
from .platform import decode as _any_dc
from .typecheck import validator as _pc
from .coercetype import coerce as _coerce


class quanta:
    ### self
    def __init__(self, *args, **kwargs):
        """Create a quanta tool on the specified host (or by default the
        host you are running on).
        """
        self._swigobj = kwargs.get('swig_object',None)
        if self._swigobj is None:
            self._swigobj = _quanta()

    def convertfreq(self, v=[ ], outunit='Hz'):
        """convertfreq converts a frequency quantity to another unit.
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'outunit': {'type': 'cStr'}}
        doc = {'v': v, 'outunit': outunit}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _convertfreq_result = _dict_dc(self._swigobj.convertfreq(_any_ec(_pc.document['v']), _str_ec(_pc.document['outunit'])))
        return _convertfreq_result

    def convertdop(self, v=[ ], outunit='km/s'):
        """convertfreq converts a velocity quantity to another unit. Units are either
        velocity or dimensionless.
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'outunit': {'type': 'cStr'}}
        doc = {'v': v, 'outunit': outunit}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _convertdop_result = _dict_dc(self._swigobj.convertdop(_any_ec(_pc.document['v']), _str_ec(_pc.document['outunit'])))
        return _convertdop_result

    def quantity(self, v=[ ], unitname=''):
        """quantity makes a quantity from a string, or from a value and a
        string. Note that a function unit exists which is a synonym for
        quantity. If only a string is given, it can be a scalar string.
        The result will be a scalar quantity.
        
        
        If a numeric value and a unit string
        are given, the numeric value can be any numeric type, and can also be
        a vector of numeric values.  print qa.map() to get a list of recognized units.
        'd' is usually days, but can be degrees (see example).
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'unitname': {'type': 'cStr'}}
        doc = {'v': v, 'unitname': unitname}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _quantity_result = _dict_dc(self._swigobj.quantity(_any_ec(_pc.document['v']), _str_ec(_pc.document['unitname'])))
        return _quantity_result

    def getvalue(self, v=[ ]):
        """getvalue returns the internal value of a quantity. It also can handle an array
        of quantities.
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _getvalue_result = self._swigobj.getvalue(_any_ec(_pc.document['v']))
        return _getvalue_result

    def getunit(self, v=[ ]):
        """getunit returns the internal unit string of a quantity
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _getunit_result = _str_dc(self._swigobj.getunit(_any_ec(_pc.document['v'])))
        return _getunit_result

    def canonical(self, v=[ ]):
        """canonical (with alias canon) gets the canonical value of a quantity
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _canonical_result = _dict_dc(self._swigobj.canonical(_any_ec(_pc.document['v'])))
        return _canonical_result

    def canon(self, v=[ ]):
        """canon gets the canonical value of a quantity
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _canon_result = _dict_dc(self._swigobj.canon(_any_ec(_pc.document['v'])))
        return _canon_result

    def convert(self, v=[ ], outunit=[ ]):
        """convert converts a quantity to another unit. If no output unit given,
        conversion is to canonical units
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'outunit': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v, 'outunit': outunit}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _convert_result = _dict_dc(self._swigobj.convert(_any_ec(_pc.document['v']), _any_ec(_pc.document['outunit'])))
        return _convert_result

    def define(self, name, v=[ ]):
        """define defines the name and value of a user defined unit
        """
        schema = {'name': {'type': 'cStr'}, 'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'name': name, 'v': v}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _define_result = self._swigobj.define(_str_ec(_pc.document['name']), _any_ec(_pc.document['v']))
        return _define_result

    def map(self, v='all'):
        """map lists the known mapping of units and constants. It has a single argument,
        which can be a coded string (no-case, minimax match):
        begin{description}
        item[all] all of the following units (not constants): also the default
        item[Prefix] known decimal prefixes
        item[SI] known SI units
        item[Customary] a set of customary units known to programs
        item[User] units defined by the user
        item[Constants] known constants (note: only 'const', 'Const', 'constants'
        and 'Constants' recognised).
        end{description}
        """
        schema = {'v': {'type': 'cStr'}}
        doc = {'v': v}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _map_result = _str_dc(self._swigobj.map(_str_ec(_pc.document['v'])))
        return _map_result

    def maprec(self, v='all'):
        """maprec returns a record with the known mapping of units and constants. It has a single argument,
        which can be a coded string (no-case, minimax match):
        begin{description}
        item[all] all of the following units (not constants): also the default
        item[Prefix] known decimal prefixes
        item[SI] known SI units
        item[Customary] a set of customary units known to programs
        item[User] units defined by the user
        end{description}
        """
        schema = {'v': {'type': 'cStr'}}
        doc = {'v': v}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _maprec_result = _dict_dc(self._swigobj.maprec(_str_ec(_pc.document['v'])))
        return _maprec_result

    def fits(self):
        """fits defines some unit names used in reading and writing FITS files.
        """
        _fits_result = self._swigobj.fits()
        return _fits_result

    def angle(self, v=[ ], prec=int(0), form=[  ], showform=False):
        """angle converts an angle quantity to a formatted string. The formatting
        information is a precision (0 is default, 6 includes +-ddd.mm.ss) and a
        string array of codes (no-case, minimax match):
        Codes include:
        begin{description}
        item[clean] delete leading/trailing superfluous separators
        item[no_d] do not show degrees part
        item[no_dm] do not show degrees and minutes part
        item[dig2] show only 2 digits of degrees in angle format
        item[time] show as time (hh:mm:ss.ttt) rather than as angle
        end{description}
        If a multi-dimensional value is given for the value $v$, the returned value
        is a string vector of a length equal to last dimension. Each string has a
        number of fields equal to the number of elements in all earlier
        dimensions. If the {em showform} is $T$, each vector element is surrounded
        by a pair of square brackets if there is more than one entry, and fields are
        separated by a ','.
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'prec': {'type': 'cInt'}, 'form': {'type': 'cStrVec', 'coerce': [_coerce.to_list,_coerce.to_strvec]}, 'showform': {'type': 'cBool'}}
        doc = {'v': v, 'prec': prec, 'form': form, 'showform': showform}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _angle_result = [_str_dc(_x) for _x in self._swigobj.angle(_any_ec(_pc.document['v']), _pc.document['prec'], [_str_ec(_x) for _x in _pc.document['form']], _pc.document['showform'])]
        return _angle_result

    def time(self, v=[ ], prec=int(0), form=[  ], showform=False):
        """time converts a time quantity to a formatted string. The formatting
        information is a precision (0 is default, 6 includes hh.mm.ss) and a
        string array of codes (no-case, minimax match):
        Codes include:
        begin{description}
        item[clean] delete leading/trailing superfluous separators
        item[no_d] do not show hours part
        item[no_dm] do not show hours and minutes part
        item[ymd] include a date as yyyy/mm/dd (date is by default not shown)
        item[dmy] include a date as ddMMMyyyy (date is by default not shown)
        item[mjd] include a date as Modified Julian Day (date is by default not shown)
        item[fits] include a date and show time in FITS format: le from OS
        item[angle] show in angle (dd.mm.ss.ttt) rather than time format
        item[day] prefix day-of-week to output
        item[local] show local time rather than UTC (add timezone offset)
        item[no_time] suppress printing of time part
        end{description}
        If a multi-dimensional value is given for the value $v$, the returned value
        is a string vector of a length equal to last dimension. Each string has a
        number of fields equal to the number of elements in all earlier
        dimensions. If the {em showform} is $T$, each vector element is surrounded
        by a pair of square brackets if there is more than one entry, and fields are
        separated by a ','.
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'prec': {'type': 'cInt'}, 'form': {'type': 'cStrVec', 'coerce': [_coerce.to_list,_coerce.to_strvec]}, 'showform': {'type': 'cBool'}}
        doc = {'v': v, 'prec': prec, 'form': form, 'showform': showform}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _time_result = [_str_dc(_x) for _x in self._swigobj.time(_any_ec(_pc.document['v']), _pc.document['prec'], [_str_ec(_x) for _x in _pc.document['form']], _pc.document['showform'])]
        return _time_result

    def add(self, v=[ ], a=[ ]):
        """add adds two quantities
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'a': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v, 'a': a}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _add_result = _dict_dc(self._swigobj.add(_any_ec(_pc.document['v']), _any_ec(_pc.document['a'])))
        return _add_result

    def sub(self, v=[ ], a=[ ]):
        """sub subtracts two quantities
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'a': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v, 'a': a}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _sub_result = _dict_dc(self._swigobj.sub(_any_ec(_pc.document['v']), _any_ec(_pc.document['a'])))
        return _sub_result

    def mul(self, v=[ ], a=[ ]):
        """mul multiplies two quantities
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'a': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v, 'a': a}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _mul_result = _dict_dc(self._swigobj.mul(_any_ec(_pc.document['v']), _any_ec(_pc.document['a'])))
        return _mul_result

    def div(self, v=[ ], a=[ ]):
        """div divides two quantities
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'a': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v, 'a': a}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _div_result = _dict_dc(self._swigobj.div(_any_ec(_pc.document['v']), _any_ec(_pc.document['a'])))
        return _div_result

    def neg(self, v=[ ]):
        """neg negates a quantity
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _neg_result = _dict_dc(self._swigobj.neg(_any_ec(_pc.document['v'])))
        return _neg_result

    def norm(self, v=[ ], a=float(-0.5)):
        """norm normalise angles in interval of $2pi$ radians. The default interval is
        from -0.5 to +0.5 of a full interval (i.e. from -180 to +180 degrees). The
        lower end of the interval can be set as a fraction of $2pi$
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'a': {'type': 'cFloat', 'coerce': _coerce.to_float}}
        doc = {'v': v, 'a': a}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _norm_result = _dict_dc(self._swigobj.norm(_any_ec(_pc.document['v']), _pc.document['a']))
        return _norm_result

    def le(self, v=[ ], a=[ ]):
        """le compares two quantities for less than or equal.
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'a': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v, 'a': a}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _le_result = self._swigobj.le(_any_ec(_pc.document['v']), _any_ec(_pc.document['a']))
        return _le_result

    def lt(self, v=[ ], a=[ ]):
        """lt compares two quantities for less than.
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'a': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v, 'a': a}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _lt_result = self._swigobj.lt(_any_ec(_pc.document['v']), _any_ec(_pc.document['a']))
        return _lt_result

    def eq(self, v=[ ], a=[ ]):
        """eq compares two quantities for equality.
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'a': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v, 'a': a}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _eq_result = self._swigobj.eq(_any_ec(_pc.document['v']), _any_ec(_pc.document['a']))
        return _eq_result

    def ne(self, v=[ ], a=[ ]):
        """ne compares two quantities for non equality.
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'a': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v, 'a': a}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _ne_result = self._swigobj.ne(_any_ec(_pc.document['v']), _any_ec(_pc.document['a']))
        return _ne_result

    def gt(self, v=[ ], a=[ ]):
        """gt compares two quantities for greater than.
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'a': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v, 'a': a}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _gt_result = self._swigobj.gt(_any_ec(_pc.document['v']), _any_ec(_pc.document['a']))
        return _gt_result

    def ge(self, v=[ ], a=[ ]):
        """ge  compares two quantities for greater than or equal.
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'a': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v, 'a': a}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _ge_result = self._swigobj.ge(_any_ec(_pc.document['v']), _any_ec(_pc.document['a']))
        return _ge_result

    def sin(self, v=[ ]):
        """sin gives sine of angle quantity
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _sin_result = _dict_dc(self._swigobj.sin(_any_ec(_pc.document['v'])))
        return _sin_result

    def cos(self, v=[ ]):
        """cos gives cosine of angle quantity
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _cos_result = _dict_dc(self._swigobj.cos(_any_ec(_pc.document['v'])))
        return _cos_result

    def tan(self, v=[ ]):
        """tan gives tangent of angle quantity
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _tan_result = _dict_dc(self._swigobj.tan(_any_ec(_pc.document['v'])))
        return _tan_result

    def asin(self, v=[ ]):
        """asin gives arcsine of non-dimensioned quantity
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _asin_result = _dict_dc(self._swigobj.asin(_any_ec(_pc.document['v'])))
        return _asin_result

    def acos(self, v=[ ]):
        """acos gives arccosine of non-dimensioned quantity
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _acos_result = _dict_dc(self._swigobj.acos(_any_ec(_pc.document['v'])))
        return _acos_result

    def atan(self, v=[ ]):
        """atan gives arctangent of non-dimensioned quantity
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _atan_result = _dict_dc(self._swigobj.atan(_any_ec(_pc.document['v'])))
        return _atan_result

    def atan2(self, v=[ ], a=[ ]):
        """atan gives arctangent of two non-dimensioned quantity
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'a': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v, 'a': a}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _atan2_result = _dict_dc(self._swigobj.atan2(_any_ec(_pc.document['v']), _any_ec(_pc.document['a'])))
        return _atan2_result

    def abs(self, v=[ ]):
        """abs gives absolute value of quantity
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _abs_result = _dict_dc(self._swigobj.abs(_any_ec(_pc.document['v'])))
        return _abs_result

    def ceil(self, v=[ ]):
        """ceil gives ceiling value of quantity
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _ceil_result = _dict_dc(self._swigobj.ceil(_any_ec(_pc.document['v'])))
        return _ceil_result

    def floor(self, v=[ ]):
        """floor gives flooring value of quantity
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _floor_result = _dict_dc(self._swigobj.floor(_any_ec(_pc.document['v'])))
        return _floor_result

    def log(self, v=[ ]):
        """log gives natural logarithm of dimensionless quantity
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _log_result = _dict_dc(self._swigobj.log(_any_ec(_pc.document['v'])))
        return _log_result

    def log10(self, v=[ ]):
        """log10 gives logarithm of dimensionless quantity
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _log10_result = _dict_dc(self._swigobj.log10(_any_ec(_pc.document['v'])))
        return _log10_result

    def exp(self, v=[ ]):
        """exp gives exponential value of dimensionless quantity
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _exp_result = _dict_dc(self._swigobj.exp(_any_ec(_pc.document['v'])))
        return _exp_result

    def sqrt(self, v=[ ]):
        """sqrt gives square root of quantity with only even powered dimensions
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _sqrt_result = _dict_dc(self._swigobj.sqrt(_any_ec(_pc.document['v'])))
        return _sqrt_result

    def compare(self, v=[ ], a=[ ]):
        """compare compares the dimensionality of units of two qauntities
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'a': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v, 'a': a}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _compare_result = self._swigobj.compare(_any_ec(_pc.document['v']), _any_ec(_pc.document['a']))
        return _compare_result

    def check(self, v):
        """check checks if the argument has a properly defined unit string
        """
        schema = {'v': {'type': 'cStr'}}
        doc = {'v': v}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _check_result = self._swigobj.check(_str_ec(_pc.document['v']))
        return _check_result

    def checkfreq(self, cm=[ ]):
        """checkfreq checks if the argument has a properly defined frequency interpretable
        unit string
        """
        schema = {'cm': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'cm': cm}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _checkfreq_result = self._swigobj.checkfreq(_any_ec(_pc.document['cm']))
        return _checkfreq_result

    def pow(self, v=[ ], a=int(1)):
        """pow raises a quantity to an integer power
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'a': {'type': 'cInt'}}
        doc = {'v': v, 'a': a}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _pow_result = _dict_dc(self._swigobj.pow(_any_ec(_pc.document['v']), _pc.document['a']))
        return _pow_result

    def constants(self, v='pi'):
        """constants gets a named constant quantity. Names (no-case, minimax) are:
        
        pi    3.14..                    3.14159
        ee    2.71..                    2.71828
        c     light vel.                2.99792e+08 m/s
        G     grav. const               6.67259e-11 N.m2/kg2
        h     Planck const              6.62608e-34 J.s
        HI    HI line                   1420.41 MHz
        R     gas const                 8.31451 J/K/mol
        NA    Avogadro number           6.02214e+23 mol-1
        e     electron charge           1.60218e-19 C
        mp    proton mass               1.67262e-27 kg
        mp_me mp/me                     1836.15
        mu0   permeability vac.         1.25664e-06 H/m
        eps0  permittivity vac.         1.60218e-19 C
        k     Boltzmann const           1.38066e-23 J/K
        F     Faraday const             96485.3 C/mol
        me    electron mass             9.10939e-31 kg
        re    electron radius           2.8179e-15 m
        a0    Bohr's radius             5.2918e-11 m
        R0    solar radius              6.9599e+08 m
        k2    IAU grav. const^2         0.000295912 AU3/d2/S0
        """
        schema = {'v': {'type': 'cStr'}}
        doc = {'v': v}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _constants_result = _dict_dc(self._swigobj.constants(_str_ec(_pc.document['v'])))
        return _constants_result

    def isangle(self, v=[ ]):
        """isangle checks if the argument is a valid angle/time quantity.
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _isangle_result = self._swigobj.isangle(_any_ec(_pc.document['v']))
        return _isangle_result

    def totime(self, v=[ ]):
        """totime converts an angle quantity (or a time) to a time quantity
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _totime_result = _dict_dc(self._swigobj.totime(_any_ec(_pc.document['v'])))
        return _totime_result

    def toangle(self, v=[ ]):
        """toangle converts a time quantity (or an angle) to an angle quantity
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _toangle_result = _dict_dc(self._swigobj.toangle(_any_ec(_pc.document['v'])))
        return _toangle_result

    def splitdate(self, v=[ ]):
        """splitdate splits a date/time quantity into a record with constituent fields
        like year, yearday, month etc. All fields will be integer (to enable use as
        index and easy personal formatting), with the exception of the {em s} field
        which is a double float. See the example for the fields returned.
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _splitdate_result = _dict_dc(self._swigobj.splitdate(_any_ec(_pc.document['v'])))
        return _splitdate_result

    def tos(self, v=[ ], prec=int(9)):
        """tos converts a quantity to a string with the precision defined with
        the {em setformat('prec')} (which defaults to 9). If the optional
        {em prec} argument is set to an integer value greater than 1, that
        precision is used in the conversion
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'prec': {'type': 'cInt'}}
        doc = {'v': v, 'prec': prec}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _tos_result = _str_dc(self._swigobj.tos(_any_ec(_pc.document['v']), _pc.document['prec']))
        return _tos_result

    def type(self):
        """type will return the tool name.
        """
        _type_result = _str_dc(self._swigobj.type())
        return _type_result

    def done(self, kill=False):
        """Currently, this method is an NOP.
        """
        schema = {'kill': {'type': 'cBool'}}
        doc = {'kill': kill}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _done_result = self._swigobj.done(_pc.document['kill'])
        return _done_result

    def unit(self, v=[ ], unitname=''):
        """unit makes a quantity from a string, or from a value and a string.
        Note that unit is a synonym for quantity (see example there).
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'unitname': {'type': 'cStr'}}
        doc = {'v': v, 'unitname': unitname}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _unit_result = _dict_dc(self._swigobj.unit(_any_ec(_pc.document['v']), _str_ec(_pc.document['unitname'])))
        return _unit_result

    def isquantity(self, v=[ ]):
        """Checks if the operand is a correct quantity
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}}
        doc = {'v': v}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _isquantity_result = self._swigobj.isquantity(_any_ec(_pc.document['v']))
        return _isquantity_result

    def setformat(self, t='', v='F'):
        """
        """
        schema = {'t': {'type': 'cStr'}, 'v': {'type': 'cStr'}}
        doc = {'t': t, 'v': v}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _setformat_result = self._swigobj.setformat(_str_ec(_pc.document['t']), _str_ec(_pc.document['v']))
        return _setformat_result

    def getformat(self, t=''):
        """getformat returns the current format value set for the different
        format possibilities. See the
        setformat function for the
        different format type descriptions. The known types are: 
        prec, aprec, tprec, long, lat, len, dtime, elev, auto, vel, freq,
        dop, unit.
        """
        schema = {'t': {'type': 'cStr'}}
        doc = {'t': t}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _getformat_result = _str_dc(self._swigobj.getformat(_str_ec(_pc.document['t'])))
        return _getformat_result

    def formxxx(self, v=[ ], format='dms', prec=int(2)):
        """form.xxx (xxx can be lat, long, len, vel, freq, dtime, unit) will format the
        input into a string using the global format information set by setformat().
        """
        schema = {'v': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'format': {'type': 'cStr'}, 'prec': {'type': 'cInt'}}
        doc = {'v': v, 'format': format, 'prec': prec}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _formxxx_result = _str_dc(self._swigobj.formxxx(_any_ec(_pc.document['v']), _str_ec(_pc.document['format']), _pc.document['prec']))
        return _formxxx_result

