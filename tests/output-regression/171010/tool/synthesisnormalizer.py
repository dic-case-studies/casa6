##################### generated by xml-casa (v2) from synthesisnormalizer.xml #######
##################### e7758c1bfca4799c53d8fddadd53d3e4 ##############################
from __future__ import absolute_import 
from CASAtools.__casac__ import synthesisnormalizer as _synthesisnormalizer
from CASAtools.typecheck import validator as _pc
from CASAtools.coercetype import coerce as _coerce
from CASAtools.synthesisimstore import synthesisimstore as _wrap_synthesisimstore

class synthesisnormalizer:
    ### self
    def __init__(self, *args, **kwargs):
        """This is used to construct {tt synthesisnormalizer} tool.
        """
        self._swigobj = kwargs.get('swig_object',None)
        if self._swigobj is None:
            self._swigobj = _synthesisnormalizer()

    def setupnormalizer(self, normpars={ }):
        """
        """
        schema = {'normpars': {'type': 'cDict'}}
        doc = {'normpars': normpars}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.setupnormalizer(_pc.document['normpars'])

    def gatherweightdensity(self):
        """
        """
        return self._swigobj.gatherweightdensity()

    def scatterweightdensity(self):
        """
        """
        return self._swigobj.scatterweightdensity()

    def gatherpsfweight(self):
        """
        """
        return self._swigobj.gatherpsfweight()

    def gatherresidual(self):
        """
        """
        return self._swigobj.gatherresidual()

    def dividepsfbyweight(self):
        """
        """
        return self._swigobj.dividepsfbyweight()

    def normalizeprimarybeam(self):
        """
        """
        return self._swigobj.normalizeprimarybeam()

    def divideresidualbyweight(self):
        """
        """
        return self._swigobj.divideresidualbyweight()

    def dividemodelbyweight(self):
        """
        """
        return self._swigobj.dividemodelbyweight()

    def multiplymodelbyweight(self):
        """
        """
        return self._swigobj.multiplymodelbyweight()

    def scattermodel(self):
        """
        """
        return self._swigobj.scattermodel()

    def getimstore(self):
        """
        """
        return _wrap_synthesisimstore(swig_object=self._swigobj.getimstore())

    def setimstore(self, imstore=None):
        """
        """
        schema = {'imstore': {'type': 'csynthesisimstoreTool'}}
        doc = {'imstore': imstore}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.setimstore(_pc.document['imstore if imstore is None else imstore._swigobj'])

    def done(self):
        """
        """
        return self._swigobj.done()

