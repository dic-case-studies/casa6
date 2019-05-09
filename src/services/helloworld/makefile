
###
### this is just a rudamentary test set up for one specific python-version
### right now to ensure that the gRPC/protobuf install and the casacore changes (for
### service list) works properly...
### 
ifeq ($(shell uname),Linux)
CC=g++
PROTOC=../../../build/lib.linux-x86_64-3.6/casatools/__casac__/bin/protoc
GRPCPLUGIN=../../../build/lib.linux-x86_64-3.6/casatools/__casac__/bin/grpc_cpp_plugin
LIBDIR=../../../build/lib.linux-x86_64-3.6/casatools/__casac__/lib
CFLAGS=-I../../../build/lib.linux-x86_64-3.6/casatools/__casac__/include
PROTODIR=../protos
loader_origin=\$$ORIGIN
else
CC=clang++
PROTOC=../../../build/lib.macosx-10.12-x86_64-3.6/casatools/__casac__/bin/protoc
GRPCPLUGIN=../../../build/lib.macosx-10.12-x86_64-3.6/casatools/__casac__/bin/grpc_cpp_plugin
LIBDIR=../../../build/lib.macosx-10.12-x86_64-3.6/casatools/__casac__/lib
CFLAGS=-I../../../build/lib.macosx-10.12-x86_64-3.6/casatools/__casac__/include
PROTODIR=../protos
loader_origin=@loader_path
endif

%.o: %.cc
	$(CC) -std=c++11 $(CFLAGS) -g -c $< -o $@

vpath %.proto $(PROTODIR)

all: client_future client server

client: helloworld.pb.o helloworld.grpc.pb.o client.o
	$(CC) -o $@ $^ -L$(LIBDIR) -Wl,-rpath,$(loader_origin)/$(LIBDIR) -L$(LIBDIR) -lgrpc++ -lgrpc -Wl,-rpath,$(loader_origin)/$(LIBDIR) -lprotobuf

client_future: helloworld.pb.o helloworld.grpc.pb.o client_future.o
	$(CC) -g -o $@ $^ -L$(LIBDIR) -Wl,-rpath,$(loader_origin)/$(LIBDIR) -L$(LIBDIR) -lgrpc++ -lgrpc -Wl,-rpath,$(loader_origin)/$(LIBDIR) -lprotobuf -lpthread

server: helloworld.pb.o helloworld.grpc.pb.o server.o
	$(CC) -o $@ $^ -L$(LIBDIR) -Wl,-rpath,$(loader_origin)/$(LIBDIR) -L$(LIBDIR) -lgrpc++ -lgrpc -Wl,-rpath,$(loader_origin)/$(LIBDIR) -lprotobuf

.PRECIOUS: %.grpc.pb.cc
%.grpc.pb.cc: %.proto
	$(PROTOC) -I$(PROTODIR) --grpc_out=. --plugin=protoc-gen-grpc=$(GRPCPLUGIN) $<

.PRECIOUS: %.pb.cc
%.pb.cc: %.proto
	$(PROTOC) -I$(PROTODIR) --cpp_out=. $<

.PRECIOUS: %.grpc.pb.h
%.grpc.pb.h: %.proto
	$(PROTOC) -I$(PROTODIR) --grpc_out=. --plugin=protoc-gen-grpc=$(GRPCPLUGIN) $<

.PRECIOUS: %.pb.h
%.pb.h: %.proto
	$(PROTOC) -I$(PROTODIR) --cpp_out=. $<

clean:
	rm -f *.o *.pb.cc *.pb.h client client_future server
