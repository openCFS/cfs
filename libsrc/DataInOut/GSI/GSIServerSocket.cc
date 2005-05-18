// Implementation of the GSI::ServerSocket class

#include "GSIServerSocket.hh"

namespace GridlibSocketInterface
{

ServerSocket::ServerSocket ( int32 port, int32 timeout )
  : Socket(timeout)
{
  if ( ! Socket::create() )
    {
         Socket::cleanup();
         //      throw SocketException ( "Could not create server socket." );
    }

  if ( ! Socket::bind ( port ) )
    {
         Socket::cleanup();
         //      throw SocketException ( "Could not bind to port." );
    }

  if ( ! Socket::listen() )
    {
        Socket::cleanup();
        //      throw SocketException ( "Could not listen to socket." );
    }

}

ServerSocket::~ServerSocket()
{
}


int32 ServerSocket::accept ( ServerSocket& sock, int32 single ) throw (SocketException)
{
  int32 ret;
  
  if ( (ret = Socket::accept ( sock )) < 0 )
    {
      throw SocketException ( "Could not accept socket." );
    }

  if (single) 
  {
      Socket::shutdownServer();
  }
  
  return ret;
}
/*
void ServerSocket::accept ( ServerSocket& sock ) throw (SocketException) 
{
  if ( ! Socket::accept ( sock ) )
  {
      throw SocketException ( "Could not accept socket." );
  }
}
*/
 
}
