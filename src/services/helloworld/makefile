###
### this is just a rudamentary test set up for one specific architecture/python-version
### right now to ensure that the gRPC/protobuf install and the casacore changes (for
### service list) works properly...
### 
PROTOC=../../../build/lib.macosx-10.12-x86_64-3.6/CASAtools/__casac__/bin/protoc
GRPCPLUGIN=../../../build/lib.macosx-10.12-x86_64-3.6/CASAtools/__casac__/bin/grpc_cpp_plugin
LIBDIR=../../../build/lib.macosx-10.12-x86_64-3.6/CASAtools/__casac__/lib
CFLAGS=-I../../../build/lib.macosx-10.12-x86_64-3.6/CASAtools/__casac__/include
PROTODIR=../protos

%.o: %.cc
	clang++ -std=c++11 $(CFLAGS) -c $< -o $@

vpath %.proto $(PROTODIR)

all: client server

client: helloworld.pb.o helloworld.grpc.pb.o client.o
	clang++ -o $@ $^ -L$(LIBDIR) -Wl,-rpath,@loader_path/$(LIBDIR) -L$(LIBDIR) -lgrpc++ -lgrpc -Wl,-rpath,@loader_path/$(LIBDIR) -lprotobuf

server: helloworld.pb.o helloworld.grpc.pb.o server.o
	clang++ -o $@ $^ -L$(LIBDIR) -Wl,-rpath,@loader_path/$(LIBDIR) -L$(LIBDIR) -lgrpc++ -lgrpc -Wl,-rpath,@loader_path/$(LIBDIR) -lprotobuf

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
	rm -f *.o *.pb.cc *.pb.h client server
