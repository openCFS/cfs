// Definition of the GSI::ClientSocket class

#ifndef GSI_CLIENTSOCKET
#define GSI_CLIENTSOCKET

#include "GSISocket.hh"
#include "GSITypeDefs.hh"

namespace GridlibSocketInterface
{

class ClientSocket : public Socket
{
 public:

  ClientSocket ( const std::string& host, const int32 port );
  virtual ~ClientSocket(){};
};
 
}

#endif //GSI_CLIENTSOCKET
