// Definition of the GSIServerSocket class

#ifndef GSIServerSocket_class
#define GSIServerSocket_class

#include "GSISocket.hh"
#include "GSISocketException.hh"

class GSIServerSocket : public GSISocket
{
 public:

  GSIServerSocket ( int port, int timeout = -1 );
  GSIServerSocket ( int timeout = -1 )
    : GSISocket(timeout) {};

  virtual ~GSIServerSocket();

  int accept ( GSIServerSocket&, int single = 0 ) throw (GSISocketException);
  //  void accept ( GSIServerSocket&) throw (GSISocketException);
};


#endif
