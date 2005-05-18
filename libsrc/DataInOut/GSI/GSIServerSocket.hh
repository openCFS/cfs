// Definition of the GSI::ServerSocket class

#ifndef GSI_SERVERSOCKET
#define GSI_SERVERSOCKET

#include "GSISocket.hh"
#include "GSISocketException.hh"
#include "GSITypeDefs.hh"

namespace GridlibSocketInterface
{

class ServerSocket : public Socket
{
 public:

  ServerSocket ( int32 port, int32 timeout = -1 );
  ServerSocket ( int32 timeout = -1 )
    : Socket(timeout) {};

  virtual ~ServerSocket();

  int32 accept ( ServerSocket&, int32 single = 0 ) throw (SocketException);
  //  void accept ( ServerSocket&) throw (SocketException);
};
 
}

#endif //GSI_SERVERSOCKET
