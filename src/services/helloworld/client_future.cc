/*
 *
 * Copyright 2015 gRPC authors.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <iostream>
#include <memory>
#include <string>
#include <future>
#include <unistd.h>
#include <functional>

#include <grpc++/grpc++.h>

#ifdef BAZEL_BUILD
#include "examples/protos/helloworld.grpc.pb.h"
#else
#include "helloworld.grpc.pb.h"
#endif

#include "casatools/Data/Opt.h"
using casatools::opt;

using grpc::Channel;
using grpc::ClientContext;
using grpc::Status;
using helloworld::HelloRequest;
using helloworld::HelloReply;
using helloworld::Greeter;

class GreeterClient {
 public:
  GreeterClient(std::shared_ptr<Channel> channel)
      : stub_(Greeter::NewStub(channel)) {}

    // Assembles the client's payload, sends it and presents the response back
    // from the server.
    casatools::OptionValue<std::string> SayHello(const std::string& user) {
    std::cout << "entering SayHello..." << std::endl;
    // Data we are sending to the server.
    HelloRequest request;
    request.set_name(user);

    // Container for the data we expect from the server.
    HelloReply reply;

    // Context for the client. It could be used to convey extra information to
    // the server and/or tweak certain RPC behaviors.
    ClientContext context;

    // The actual RPC.

    std::cout << "SayHello sleeping..." << std::endl;
    sleep(10);
    std::cout << "SayHello sending request..." << std::endl;
    Status status = stub_->SayHello(&context, request, &reply);
    std::cout << "SayHello received result..." << std::endl;

    // Act upon its status.
    if (status.ok()) {
      std::cout << "SayHello result = " << reply.message( ) << std::endl;
      return opt<std::string>::some(reply.message());
    } else {
      std::cout << "SayHello no result to retrieve..." << std::endl;
      return opt<std::string>::none( );
    }
  }

 private:
  std::unique_ptr<Greeter::Stub> stub_;
};

int main(int argc, char** argv) {
  // Instantiate the client. It requires a channel, out of which the actual RPCs
  // are created. This channel models a connection to an endpoint (in this case,
  // localhost at port 50051). We indicate that the channel isn't authenticated
  // (use of InsecureChannelCredentials()).
  GreeterClient greeter(grpc::CreateChannel(
      "localhost:50051", grpc::InsecureChannelCredentials()));
  std::string user("world");

  auto reply = std::async( std::launch::async, &GreeterClient::SayHello, &greeter, user );
  std::cout << "\tentering loop" << std::endl;
  while ( reply.wait_for(std::chrono::milliseconds(10)) != std::future_status::ready ) {
      std::cout << "\t...waiting... " << std::endl;
      sleep(1);
  }

  // calling <future>.get(...) more than once will result in exceptions...
  auto result = reply.get( );
  if ( result.has_value( ) )
      std::cout << "\tGreeter received: " << result.get( ) << std::endl;
  else
      std::cout << "\tGreeter request failed..." << std::endl;

  return 0;
}
