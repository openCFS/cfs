// Definition of the GSIClientSocket class

#ifndef GSIClientSocket_class
#define GSIClientSocket_class

#include "GSISocket.hh"


class GSIClientSocket : public GSISocket
{
 public:

  GSIClientSocket ( const std::string& host, const int port );
  virtual ~GSIClientSocket(){};
};


#endif
