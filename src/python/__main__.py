from __future__ import absolute_import
import sys
from .config import build

for flag in sys.argv:
    if flag == '--grpc-link':
        print(' '.join(map(lambda x: repr(x),build['build.flags.link.grpc'])))
    if flag == '--grpc-compile':
        print(' '.join(map(lambda x: repr(x),build['build.flags.compile.grpc'])))
    if flag == '--grpc-protoc':
        print(build['build.compiler.protoc'])
