// Implementation of the GSIServerSocket class

#include "GSIServerSocket.hh"
#include "GSISocketException.hh"


GSIServerSocket::GSIServerSocket ( int port, int timeout )
  : GSISocket(timeout)
{
  if ( ! GSISocket::create() )
    {
         GSISocket::cleanup();
         //      throw GSISocketException ( "Could not create server socket." );
    }

  if ( ! GSISocket::bind ( port ) )
    {
         GSISocket::cleanup();
         //      throw GSISocketException ( "Could not bind to port." );
    }

  if ( ! GSISocket::listen() )
    {
        GSISocket::cleanup();
        //      throw GSISocketException ( "Could not listen to socket." );
    }

}

GSIServerSocket::~GSIServerSocket()
{
}


int GSIServerSocket::accept ( GSIServerSocket& sock, int single ) throw (GSISocketException)
{
  int ret;
  
  if ( (ret = GSISocket::accept ( sock )) < 0 )
    {
      throw GSISocketException ( "Could not accept socket." );
    }

  if (single) 
  {
      GSISocket::shutdownServer();
  }
  
  return ret;
}
/*
void GSIServerSocket::accept ( GSIServerSocket& sock ) throw (GSISocketException) 
{
  if ( ! GSISocket::accept ( sock ) )
  {
      throw GSISocketException ( "Could not accept socket." );
  }
}
*/
