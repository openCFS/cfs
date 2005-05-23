#ifdef USE_RCSID
static const char RCSid_GSIServerSocket[] = "$Id$";
#endif

/*----------------------------------------------------------------------
  |
  |
  | $Log$
  | Revision 1.4  2005/05/23 22:15:31  ahauck
  | Corrected indentation of all cc/hh-files according to xemacs' 'gnu' style
  | and according to our coding rules, so all TABs were replaced by spaces.
  |
  | Revision 1.3  2005/05/18 19:26:03  strieben
  | Upgraded GSI library to newest available version.
  |
  | Revision 1.1.1.1  2004/08/31 15:53:00  simon
  | Initial GSI import
  |
  |
  +---------------------------------------------------------------------*/


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
