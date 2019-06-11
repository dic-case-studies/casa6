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
    if flag == '--compiler-cc':
        print(build['build.compiler.cc'])
    if flag == '--compiler-cxx':
        print(build['build.compiler.cxx'])
    if flag == '--help':
        print("--compiler-cc\t\tpath to C compiler used to build casatools")
        print("--compiler-cxx\t\tpath to C++ compiler used to build casatools")
        print("--grpc-compile\t\tflags to compile C++ source files")
        print("--grpc-link\t\tflags to use to link gRPC C++ applications")
        print("--grpc-libpath\t\tpath to gRPC C++ libraries")
        print("--grpc-protocpp\t\tpath to C++ protoc compiler")
        print("--grpc-protopy\t\tpath to Python protoc compiler")
        print("--proto-registrar\tpath to registrar protobuf spec")
        print("--proto-shutdown\tpath to shutdown protobuf spec")
