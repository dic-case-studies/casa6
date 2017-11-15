##################### generated by xml-casa (v2) from synthesisutils.xml ############
##################### 446f8d72a8d4a8d7938e884151ff05b6 ##############################
from __future__ import absolute_import 
from CASAtools.__casac__ import synthesisutils as _synthesisutils
from CASAtools.typecheck import validator as _pc
from CASAtools.coercetype import coerce as _coerce
from CASAtools.synthesisimstore import synthesisimstore as _wrap_synthesisimstore

class synthesisutils:
    ### self
    def __init__(self, *args, **kwargs):
        """This is used to construct {tt synthesisutils} tool.
        """
        self._swigobj = kwargs.get('swig_object',None)
        if self._swigobj is None:
            self._swigobj = _synthesisutils()

    def contdatapartition(self, selpars={ }, npart=int(1)):
        """
        """
        schema = {'selpars': {'type': 'cDict'}, 'npart': {'type': 'cInt'}}
        doc = {'selpars': selpars, 'npart': npart}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.contdatapartition(_pc.document['selpars'], _pc.document['npart'])

    def cubedatapartition(self, selpars={ }, npart=int(1), fstart=[ ], fend=[ ], frame='LSRK'):
        """returns a dictionary with data spectral parttiion that maps  data  to  nparts
        of the input range frequency... usually to be used for doing data selection
        when imaging a cube from fstart to fend in npart subcubes
        """
        schema = {'selpars': {'type': 'cDict'}, 'npart': {'type': 'cInt'}, 'fstart': {'type': 'cVariant'}, 'fend': {'type': 'cVariant'}, 'frame': {'type': 'cStr'}}
        doc = {'selpars': selpars, 'npart': npart, 'fstart': fstart, 'fend': fend, 'frame': frame}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.cubedatapartition(_pc.document['selpars'], _pc.document['npart'], _pc.document['fstart'], _pc.document['fend'], _pc.document['frame'])

    def cubeimagepartition(self, impars={ }, npart=int(1)):
        """
        """
        schema = {'impars': {'type': 'cDict'}, 'npart': {'type': 'cInt'}}
        doc = {'impars': impars, 'npart': npart}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.cubeimagepartition(_pc.document['impars'], _pc.document['npart'])

    def cubedataimagepartition(self, selpars={ }, incsys={ }, npart=int(1), nchannel=int(1)):
        """
        """
        schema = {'selpars': {'type': 'cDict'}, 'incsys': {'type': 'cDict'}, 'npart': {'type': 'cInt'}, 'nchannel': {'type': 'cInt'}}
        doc = {'selpars': selpars, 'incsys': incsys, 'npart': npart, 'nchannel': nchannel}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.cubedataimagepartition(_pc.document['selpars'], _pc.document['incsys'], _pc.document['npart'], _pc.document['nchannel'])

    def checkselectionparams(self, selpars={ }):
        """
        """
        schema = {'selpars': {'type': 'cDict'}}
        doc = {'selpars': selpars}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.checkselectionparams(_pc.document['selpars'])

    def checkimageparams(self, impars={ }):
        """
        """
        schema = {'impars': {'type': 'cDict'}}
        doc = {'impars': impars}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.checkimageparams(_pc.document['impars'])

    def checkgridparams(self, gridpars={ }):
        """
        """
        schema = {'gridpars': {'type': 'cDict'}}
        doc = {'gridpars': gridpars}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.checkgridparams(_pc.document['gridpars'])

    def updateimpars(self, impars={ }):
        """
        """
        schema = {'impars': {'type': 'cDict'}}
        doc = {'impars': impars}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.updateimpars(_pc.document['impars'])

    def getOptimumSize(self, size=int(100)):
        """
        """
        schema = {'size': {'type': 'cInt'}}
        doc = {'size': size}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.getOptimumSize(_pc.document['size'])

    def done(self):
        """
        """
        return self._swigobj.done()

