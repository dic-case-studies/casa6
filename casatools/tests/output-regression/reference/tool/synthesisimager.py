##################### generated by xml-casa (v2) from synthesisimager.xml ###########
##################### 31f256e80795a946d4bd76b14412671d ##############################
from __future__ import absolute_import 
from .__casac__ import synthesisimager as _synthesisimager
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
from .synthesisimstore import synthesisimstore as _wrap_synthesisimstore

class synthesisimager:
    ### self
    def __init__(self, *args, **kwargs):
        """This is used to construct {tt synthesisimager} tool.
        """
        self._swigobj = kwargs.get('swig_object',None)
        if self._swigobj is None:
            self._swigobj = _synthesisimager()

    def selectdata(self, selpars={ }):
        """
        """
        schema = {'selpars': {'type': 'cDict'}}
        doc = {'selpars': selpars}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _selectdata_result = self._swigobj.selectdata(_dict_ec(_pc.document['selpars']))
        return _selectdata_result

    def tuneselectdata(self):
        """
        """
        _tuneselectdata_result = _dict_dc(self._swigobj.tuneselectdata())
        return _tuneselectdata_result

    def defineimage(self, impars={ }, gridpars={ }):
        """
        """
        schema = {'impars': {'type': 'cDict'}, 'gridpars': {'type': 'cDict'}}
        doc = {'impars': impars, 'gridpars': gridpars}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _defineimage_result = self._swigobj.defineimage(_dict_ec(_pc.document['impars']), _dict_ec(_pc.document['gridpars']))
        return _defineimage_result

    def setdata(self, msname='', spw='', freqbeg='', freqend='', freqframe='LSRK', field='', antenna='', timestr='', scan='', obs='', state='', uvdist='', taql='', usescratch=False, readonly=False, incrmodel=False):
        """Select data from one MS. Call this function in succession if there are
        multiple MSs.
        """
        schema = {'msname': {'type': 'cStr'}, 'spw': {'type': 'cStr'}, 'freqbeg': {'type': 'cStr'}, 'freqend': {'type': 'cStr'}, 'freqframe': {'type': 'cStr'}, 'field': {'type': 'cStr'}, 'antenna': {'type': 'cStr'}, 'timestr': {'type': 'cStr'}, 'scan': {'type': 'cStr'}, 'obs': {'type': 'cStr'}, 'state': {'type': 'cStr'}, 'uvdist': {'type': 'cStr'}, 'taql': {'type': 'cStr'}, 'usescratch': {'type': 'cBool'}, 'readonly': {'type': 'cBool'}, 'incrmodel': {'type': 'cBool'}}
        doc = {'msname': msname, 'spw': spw, 'freqbeg': freqbeg, 'freqend': freqend, 'freqframe': freqframe, 'field': field, 'antenna': antenna, 'timestr': timestr, 'scan': scan, 'obs': obs, 'state': state, 'uvdist': uvdist, 'taql': taql, 'usescratch': usescratch, 'readonly': readonly, 'incrmodel': incrmodel}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _setdata_result = self._swigobj.setdata(_str_ec(_pc.document['msname']), _str_ec(_pc.document['spw']), _str_ec(_pc.document['freqbeg']), _str_ec(_pc.document['freqend']), _str_ec(_pc.document['freqframe']), _str_ec(_pc.document['field']), _str_ec(_pc.document['antenna']), _str_ec(_pc.document['timestr']), _str_ec(_pc.document['scan']), _str_ec(_pc.document['obs']), _str_ec(_pc.document['state']), _str_ec(_pc.document['uvdist']), _str_ec(_pc.document['taql']), _pc.document['usescratch'], _pc.document['readonly'], _pc.document['incrmodel'])
        return _setdata_result

    def setimage(self, imagename='', nx=int(128), ny=int(-1), cellx=[ ], celly=[ ], stokes='I', phasecenter=[ ], nchan=int(-1), freqstart=[ ], freqstep=[ ], restfreq=[ ], facets=int(1), ftmachine='gridft', ntaylorterms=int(1), reffreq=[ ], projection='SIN', distance=[ ], freqframe='LSRK', tracksource=False, trackdir=[ ], overwrite=True, padding=float(1.0), useautocorr=False, usedoubleprec=True, wprojplanes=int(1), convfunc='SF', startmodel='', aterm=True, psterm=True, mterm=False, wbawp=True, cfcache='', dopointing=False, dopbcorr=True, conjbeams=True, computepastep=float(360.0), rotatepastep=float(5.0)):
        """Define the image coordinate systems and shapes.
        """
        schema = {'imagename': {'type': 'cStr'}, 'nx': {'type': 'cInt'}, 'ny': {'type': 'cInt'}, 'cellx': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'celly': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'stokes': {'type': 'cStr'}, 'phasecenter': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'nchan': {'type': 'cInt'}, 'freqstart': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'freqstep': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'restfreq': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'facets': {'type': 'cInt'}, 'ftmachine': {'type': 'cStr'}, 'ntaylorterms': {'type': 'cInt'}, 'reffreq': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'projection': {'type': 'cStr'}, 'distance': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'freqframe': {'type': 'cStr'}, 'tracksource': {'type': 'cBool'}, 'trackdir': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'overwrite': {'type': 'cBool'}, 'padding': {'type': 'cFloat', 'coerce': _coerce.to_float}, 'useautocorr': {'type': 'cBool'}, 'usedoubleprec': {'type': 'cBool'}, 'wprojplanes': {'type': 'cInt'}, 'convfunc': {'type': 'cStr'}, 'startmodel': {'type': 'cStr'}, 'aterm': {'type': 'cBool'}, 'psterm': {'type': 'cBool'}, 'mterm': {'type': 'cBool'}, 'wbawp': {'type': 'cBool'}, 'cfcache': {'type': 'cStr'}, 'dopointing': {'type': 'cBool'}, 'dopbcorr': {'type': 'cBool'}, 'conjbeams': {'type': 'cBool'}, 'computepastep': {'type': 'cFloat', 'coerce': _coerce.to_float}, 'rotatepastep': {'type': 'cFloat', 'coerce': _coerce.to_float}}
        doc = {'imagename': imagename, 'nx': nx, 'ny': ny, 'cellx': cellx, 'celly': celly, 'stokes': stokes, 'phasecenter': phasecenter, 'nchan': nchan, 'freqstart': freqstart, 'freqstep': freqstep, 'restfreq': restfreq, 'facets': facets, 'ftmachine': ftmachine, 'ntaylorterms': ntaylorterms, 'reffreq': reffreq, 'projection': projection, 'distance': distance, 'freqframe': freqframe, 'tracksource': tracksource, 'trackdir': trackdir, 'overwrite': overwrite, 'padding': padding, 'useautocorr': useautocorr, 'usedoubleprec': usedoubleprec, 'wprojplanes': wprojplanes, 'convfunc': convfunc, 'startmodel': startmodel, 'aterm': aterm, 'psterm': psterm, 'mterm': mterm, 'wbawp': wbawp, 'cfcache': cfcache, 'dopointing': dopointing, 'dopbcorr': dopbcorr, 'conjbeams': conjbeams, 'computepastep': computepastep, 'rotatepastep': rotatepastep}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _setimage_result = self._swigobj.setimage(_str_ec(_pc.document['imagename']), _pc.document['nx'], _pc.document['ny'], _any_ec(_pc.document['cellx']), _any_ec(_pc.document['celly']), _str_ec(_pc.document['stokes']), _any_ec(_pc.document['phasecenter']), _pc.document['nchan'], _any_ec(_pc.document['freqstart']), _any_ec(_pc.document['freqstep']), _any_ec(_pc.document['restfreq']), _pc.document['facets'], _str_ec(_pc.document['ftmachine']), _pc.document['ntaylorterms'], _any_ec(_pc.document['reffreq']), _str_ec(_pc.document['projection']), _any_ec(_pc.document['distance']), _str_ec(_pc.document['freqframe']), _pc.document['tracksource'], _any_ec(_pc.document['trackdir']), _pc.document['overwrite'], _pc.document['padding'], _pc.document['useautocorr'], _pc.document['usedoubleprec'], _pc.document['wprojplanes'], _str_ec(_pc.document['convfunc']), _str_ec(_pc.document['startmodel']), _pc.document['aterm'], _pc.document['psterm'], _pc.document['mterm'], _pc.document['wbawp'], _str_ec(_pc.document['cfcache']), _pc.document['dopointing'], _pc.document['dopbcorr'], _pc.document['conjbeams'], _pc.document['computepastep'], _pc.document['rotatepastep'])
        return _setimage_result

    def setweighting(self, type='natural', rmode='norm', noise=[ ], robust=float(0.0), fieldofview=[ ], npixels=int(0), multifield=False, uvtaper=[  ]):
        """
        """
        schema = {'type': {'type': 'cStr'}, 'rmode': {'type': 'cStr'}, 'noise': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'robust': {'type': 'cFloat', 'coerce': _coerce.to_float}, 'fieldofview': {'type': 'cVariant', 'coerce': [_coerce.to_variant]}, 'npixels': {'type': 'cInt'}, 'multifield': {'type': 'cBool'}, 'uvtaper': {'type': 'cStrVec', 'coerce': [_coerce.to_list,_coerce.to_strvec]}}
        doc = {'type': type, 'rmode': rmode, 'noise': noise, 'robust': robust, 'fieldofview': fieldofview, 'npixels': npixels, 'multifield': multifield, 'uvtaper': uvtaper}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _setweighting_result = self._swigobj.setweighting(_str_ec(_pc.document['type']), _str_ec(_pc.document['rmode']), _any_ec(_pc.document['noise']), _pc.document['robust'], _any_ec(_pc.document['fieldofview']), _pc.document['npixels'], _pc.document['multifield'], [_str_ec(_x) for _x in _pc.document['uvtaper']])
        return _setweighting_result

    def makepsf(self):
        """
        """
        _makepsf_result = self._swigobj.makepsf()
        return _makepsf_result

    def predictmodel(self):
        """
        """
        _predictmodel_result = self._swigobj.predictmodel()
        return _predictmodel_result

    def drygridding(self, cflist=[ '' ]):
        """
        """
        schema = {'cflist': {'type': 'cStrVec', 'coerce': [_coerce.to_list,_coerce.to_strvec]}}
        doc = {'cflist': cflist}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _drygridding_result = self._swigobj.drygridding([_str_ec(_x) for _x in _pc.document['cflist']])
        return _drygridding_result

    def fillcfcache(self, cflist=[ '' ], ftmname='', cfcpath='', pstermon=False, atermon=True, conjbeams=True):
        """
        """
        schema = {'cflist': {'type': 'cStrVec', 'coerce': [_coerce.to_list,_coerce.to_strvec]}, 'ftmname': {'type': 'cStr'}, 'cfcpath': {'type': 'cStr'}, 'pstermon': {'type': 'cBool'}, 'atermon': {'type': 'cBool'}, 'conjbeams': {'type': 'cBool'}}
        doc = {'cflist': cflist, 'ftmname': ftmname, 'cfcpath': cfcpath, 'pstermon': pstermon, 'atermon': atermon, 'conjbeams': conjbeams}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _fillcfcache_result = self._swigobj.fillcfcache([_str_ec(_x) for _x in _pc.document['cflist']], _str_ec(_pc.document['ftmname']), _str_ec(_pc.document['cfcpath']), _pc.document['pstermon'], _pc.document['atermon'], _pc.document['conjbeams'])
        return _fillcfcache_result

    def reloadcfcache(self):
        """
        """
        _reloadcfcache_result = self._swigobj.reloadcfcache()
        return _reloadcfcache_result

    def executemajorcycle(self, controls={ }):
        """
        """
        schema = {'controls': {'type': 'cDict'}}
        doc = {'controls': controls}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _executemajorcycle_result = self._swigobj.executemajorcycle(_dict_ec(_pc.document['controls']))
        return _executemajorcycle_result

    def makepb(self):
        """
        """
        _makepb_result = self._swigobj.makepb()
        return _makepb_result

    def makesdimage(self):
        """
        """
        _makesdimage_result = self._swigobj.makesdimage()
        return _makesdimage_result

    def makesdpsf(self):
        """
        """
        _makesdpsf_result = self._swigobj.makesdpsf()
        return _makesdpsf_result

    def getimstore(self, id=int(0)):
        """
        """
        schema = {'id': {'type': 'cInt'}}
        doc = {'id': id}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _getimstore_result = _wrap_synthesisimstore(swig_object=self._swigobj.getimstore(_pc.document['id']))
        return _getimstore_result

    def getcsys(self):
        """
        """
        _getcsys_result = _dict_dc(self._swigobj.getcsys())
        return _getcsys_result

    def updatenchan(self):
        """
        """
        _updatenchan_result = self._swigobj.updatenchan()
        return _updatenchan_result

    def getweightdensity(self):
        """
        """
        _getweightdensity_result = self._swigobj.getweightdensity()
        return _getweightdensity_result

    def setweightdensity(self):
        """
        """
        _setweightdensity_result = self._swigobj.setweightdensity()
        return _setweightdensity_result

    def done(self):
        """
        """
        _done_result = self._swigobj.done()
        return _done_result

