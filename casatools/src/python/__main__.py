from __future__ import absolute_import
import sys
from .config import build
import os as __os

for flag in sys.argv:
    ## qmake has difficuties with quoted paths, e.g. '-I/some/path/for/grpc/include/files'
    ## if such quoting is required (e.g. for paths with spaces in them), then probably
    ## another -quoted flag should be introduced to change the default behavior
    if flag == '--grpc-link':
        print(' '.join(build['build.flags.link.grpc']))
    if flag == '--grpc-compile':
        print(' '.join(build['build.flags.compile.grpc']))
    if flag == '--grpc-protocpp':
        print(build['build.compiler.protocpp'])
    if flag == '--grpc-protopy':
        print(build['build.compiler.protopy'])
    if flag == '--grpc-libpath':
        print(' '.join(map(lambda y: y[2:],filter(lambda x: x.startswith('-L'), build['build.flags.link.grpc']))))
    if flag == '--proto-registrar':
        print(__os.path.join(__os.path.dirname(__os.path.abspath(__file__)),'__casac__','proto','registrar.proto'))
    if flag == '--proto-shutdown':
        print(__os.path.join(__os.path.dirname(__os.path.abspath(__file__)),'__casac__','proto','shutdown.proto'))
    if flag == '--proto-ping':
        print(__os.path.join(__os.path.dirname(__os.path.abspath(__file__)),'__casac__','proto','ping.proto'))
    if flag == '--compiler-cc':
        print(build['build.compiler.cc'])
    if flag == '--compiler-cxx':
        print(build['build.compiler.cxx'])
    if flag == '--compiler-xml':
        print(build['build.compiler.xml-casa'])
    if flag == '--update-user-data' or flag.startswith("--update-user-data="):
        import casatools as _ct
        from subprocess import Popen, PIPE

        def _execute(cmd):
            with open(__os.devnull,"w") as fnull:
                popen = Popen(cmd, stdout=PIPE, stderr=fnull, universal_newlines=True)
            for stdout_line in iter(popen.stdout.readline, ""):
                yield stdout_line
            popen.stdout.close()
            return_code = popen.wait()
            if return_code:
                raise subprocess.CalledProcessError(return_code, cmd)

        from casatools import ctuser as cfg
        _user_data = __os.path.expanduser("~/.casa/data")
        if flag.startswith("--update-user-data="):
            _user_data = flag[19:]
        elif hasattr(cfg,'rundata'):
            _user_data = __os.path.expanduser(cfg.rundata)

        if not __os.path.exists(_user_data):
            try:
                __os.makedirs(_user_data)
            except:
                sys.exit('error, could not create %s' % _user_data)
        elif __os.path.isfile(_user_data):
            sys.exit('error, ~/.casa/data exists and is a file instead of a directory')

        try:
            print("attempting to install/update runtime data in %s" % _user_data)
            for line in _execute([ 'rsync', '-avz', 'rsync://casa-rsync.nrao.edu/casa-data', _user_data ]):
                sys.stdout.write(".")
                sys.stdout.flush( )
            print("")
        except: sys.exit('data fetch failed...')

    if flag == '--help':
        print("--compiler-cc\t\t\tpath to C compiler used to build casatools")
        print("--compiler-cxx\t\t\tpath to C++ compiler used to build casatools")
        print("--compiler-xml\t\t\tpath to compiler used to generate bindings from XML")
        print("--grpc-compile\t\t\tflags to compile C++ source files")
        print("--grpc-link\t\t\tflags to use to link gRPC C++ applications")
        print("--grpc-libpath\t\t\tpath to gRPC C++ libraries")
        print("--grpc-protocpp\t\t\tpath to C++ protoc compiler")
        print("--grpc-protopy\t\t\tpath to Python protoc compiler")
        print("--proto-registrar\t\tpath to registrar protobuf spec")
        print("--proto-shutdown\t\tpath to shutdown protobuf spec")
        print("--proto-ping\t\t\tpath to ping protobuf spec")
        print("--update-user-data[=dir]\tinstall or update measures data in ~/.casa/data")
