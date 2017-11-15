##################### generated by xml-casa (v2) from utils.xml #####################
##################### 1cf6407697ca1942b37e37ea830645c3 ##############################
from __future__ import absolute_import 
from .__casac__ import utils as _utils
from .typecheck import validator as _pc
from .coercetype import coerce as _coerce


class utils:
    ### self
    def __init__(self, *args, **kwargs):
        """
        """
        self._swigobj = kwargs.get('swig_object',None)
        if self._swigobj is None:
            self._swigobj = _utils()

    def verify(self, input={ }, xmldescriptor=[ ], throwexecpt=False):
        """
        """
        schema = {'input': {'type': 'cDict'}, 'xmldescriptor': {'type': 'cVariant'}, 'throwexecpt': {'type': 'cBool'}}
        doc = {'input': input, 'xmldescriptor': xmldescriptor, 'throwexecpt': throwexecpt}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.verify(_pc.document['input'], _pc.document['xmldescriptor'], _pc.document['throwexecpt'])

    def setconstraints(self, xmldescriptor=[ ]):
        """
        """
        schema = {'xmldescriptor': {'type': 'cVariant'}}
        doc = {'xmldescriptor': xmldescriptor}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.setconstraints(_pc.document['xmldescriptor'])

    def verifyparam(self, param={ }):
        """
        """
        schema = {'param': {'type': 'cDict'}}
        doc = {'param': param}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.verifyparam(_pc.document['param'])

    def expandparam(self, name='', value=[ ]):
        """
        """
        schema = {'name': {'type': 'cStr'}, 'value': {'type': 'cVariant'}}
        doc = {'name': name, 'value': value}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.expandparam(_pc.document['name'], _pc.document['value'])

    def torecord(self, input):
        """
        """
        schema = {'input': {'type': 'cStr'}}
        doc = {'input': input}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.torecord(_pc.document['input'])

    def toxml(self, input={ }, asfile=False, filename='recordas.xml'):
        """
        """
        schema = {'input': {'type': 'cDict'}, 'asfile': {'type': 'cBool'}, 'filename': {'type': 'cStr'}}
        doc = {'input': input, 'asfile': asfile, 'filename': filename}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.toxml(_pc.document['input'], _pc.document['asfile'], _pc.document['filename'])

    def getrc(self, rcvar=''):
        """
        """
        schema = {'rcvar': {'type': 'cStr'}}
        doc = {'rcvar': rcvar}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.getrc(_pc.document['rcvar'])

    def removetable(self, tablenames=[  ]):
        """
        """
        schema = {'tablenames': {'type': 'cStrVec', 'coerce': _coerce.to_strvec}}
        doc = {'tablenames': tablenames}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.removetable(_pc.document['tablenames'])

    def tableinfo(self, tablename=''):
        """Currently this only returns the pid of the process locking the table (lockpid), if the lock
        is permanent (lockperm), and the status (lockstatus) -- 'not in use', 'open', 'read', 'write',
        or 'unknown'. However, the hope is that this will eventually return a complete description of
        the table.
        
        """
        schema = {'tablename': {'type': 'cStr'}}
        doc = {'tablename': tablename}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.tableinfo(_pc.document['tablename'])

    def lockedtables(self):
        """
        """
        return self._swigobj.lockedtables()

    def hostinfo(self):
        """
        """
        return self._swigobj.hostinfo()

    def c_exception(self):
        """Returns detailed information from the last CASA C++ exception (i.e., AipsError).  The
        exception message and the stack trace (mangled; use the shell's c++filt to demangle)
        from the last CASA C++ exception.  The information is from the last one generated
        and may not represent an exception from the last action; c_exception_clear can be
        used to remove stale information.  The information's exception might also
        have been caught in the C++ code and not have been translated into a Python-level
        exception.
        
        """
        return self._swigobj.c_exception()

    def c_exception_clear(self):
        """Clears the CASA C++ exception information.  This allows the user to be sure that
        information retrieved using c_exception is not from an exception in the
        distant past.
        
        """
        return self._swigobj.c_exception_clear()

    def _crash_reporter_initialize(self, crashDumpDirectory, crashDumpPosterApplication, crashPostingUrl, logFile):
        """Initializes the crash reporter which will generate a crash report if casapy
        crashes.  For reporter purposes a crash is the reception of an signal by
        casapy which would normally result in the program being terminated.  This includes
        segfaults, aborts, etc., plus any unhandled C++ exceptions (C++ generates an
        abort signal for unhandled exceptions).  This method is intended for use by the
        casapy infrastructure and should not be called by other code or by users; however,
        the call will only install the crash reporter the first time it is called so any
        subsequent calls should be no-ops.  Returns true if initialization occurred and
        false if the crash reporter was stubbed out (i.e., symbol UseCrashReporter was
        not defined).
        
        """
        schema = {'crashDumpDirectory': {'type': 'cStr'}, 'crashDumpPosterApplication': {'type': 'cStr'}, 'crashPostingUrl': {'type': 'cStr'}, 'logFile': {'type': 'cStr'}}
        doc = {'crashDumpDirectory': crashDumpDirectory, 'crashDumpPosterApplication': crashDumpPosterApplication, 'crashPostingUrl': crashPostingUrl, 'logFile': logFile}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj._crash_reporter_initialize(_pc.document['crashDumpDirectory'], _pc.document['crashDumpPosterApplication'], _pc.document['crashPostingUrl'], _pc.document['logFile'])

    def _trigger_segfault(self, faultType=int(0)):
        """This triggers a segfault for testing the crash reporter.  Obviously you
        shouldn't call this unless that's what you want.  It's in here for
        development/debugging purposes and ought to be removed before you see this.
        
        """
        schema = {'faultType': {'type': 'cInt'}}
        doc = {'faultType': faultType}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj._trigger_segfault(_pc.document['faultType'])

    def initialize(self, default_path):
        """returns true if initalization was performed; returns false if initialization was already done
        """
        schema = {'default_path': {'type': 'cStrVec', 'coerce': _coerce.to_strvec}}
        doc = {'default_path': default_path}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.initialize(_pc.document['default_path'])

    def defaultpath(self):
        """Returns the default data path. This path is used unless the user has set the current path to something else using the setpath function.
        """
        return self._swigobj.defaultpath()

    def setpath(self, dirs=[  ]):
        """Sets the data path to the specified list of directories. Returns true if all directories were added
        returns false otherwise.
        """
        schema = {'dirs': {'type': 'cStrVec', 'coerce': _coerce.to_strvec}}
        doc = {'dirs': dirs}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.setpath(_pc.document['dirs'])

    def getpath(self):
        """Returns the list of directories that are currently in the data path.
        """
        return self._swigobj.getpath()

    def clearpath(self):
        """Removes all directories from the data path.
        """
        return self._swigobj.clearpath()

    def resolve(self, path=''):
        """If the provided path already represents a file or a directory, it is returned. If it does not,
        this function tries to find a complete path by matching up this partial directory with the
        elements of the data path.
        """
        schema = {'path': {'type': 'cStr'}}
        doc = {'path': path}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.resolve(_pc.document['path'])

    def version(self):
        """Returns a four element vector representing the version (major, minor, patch and feature).
        """
        return self._swigobj.version()

    def version_desc(self):
        """The descriptive string describes a particular packaged version. During a development
        cycle there are different sorts of packaged distributions. For example, a development
        version ("DEV") or a release version ("REL").
        """
        return self._swigobj.version_desc()

    def version_info(self):
        """Returns a description string that includes the version information and the descriptive string..
        """
        return self._swigobj.version_info()

    def version_string(self):
        """Returns a description string that includes the version information and the descriptive string..
        """
        return self._swigobj.version_string()

    def compare_version(self, comparitor, vec):
        """Returns a description string that includes the version information and the descriptive string..
        """
        schema = {'comparitor': {'type': 'cStr'}, 'vec': {'type': 'cIntVec', 'coerce': _coerce.to_intvec}}
        doc = {'comparitor': comparitor, 'vec': vec}
        assert _pc.validate(doc,schema), str(_pc.errors)
        return self._swigobj.compare_version(_pc.document['comparitor'], _pc.document['vec'])

