// Implementation of the GSIClientSocket class

#include <iostream>

#include "GSIClientSocket.hh"
#include "GSISocketException.hh"


GSIClientSocket::GSIClientSocket ( const std::string& host, const int port )
{
  
  if ( ! GSISocket::create() )
    {
        GSISocket::cleanup();
        
        //    throw GSISocketException ( "Could not create client socket." );
    }

  if ( ! GSISocket::connect ( host, port ) )
    {
        //      throw GSISocketException ( "Could not bind to port." );

        std::cout << "Wixenkacke! Connect geht nicht!";
        
        GSISocket::cleanup();
    }

}

