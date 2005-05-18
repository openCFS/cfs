// GSI::SocketException class


#ifndef GSI_SOCKETEXCEPTION
#define GSI_SOCKETEXCEPTION

#include <string>

#include "GSIException.hh"

namespace GridlibSocketInterface
{

class SocketException : public Exception
{
 public:
  SocketException ( std::string s ) : Exception(s) {};
  ~SocketException (){};
};
 
}

#endif //GSI_SOCKETEXCEPTION
