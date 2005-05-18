// Implementation of the GSI::ClientSocket class

#include <iostream>

#include "GSIClientSocket.hh"

namespace GridlibSocketInterface
{

ClientSocket::ClientSocket ( const std::string& host, const int32 port )
{
  
  if ( ! Socket::create() )
    {
        Socket::cleanup();

        //    throw SocketException ( "Could not create client socket." );
    }

  if ( ! Socket::connect ( host, port ) )
    {
        //      throw SocketException ( "Could not bind to port." );

        std::cout << "Wixenkacke! Connect geht nicht!";
        
        Socket::cleanup();
    }

}
 
}
