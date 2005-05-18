#ifdef USE_RCSID
static const char RCSid_GSIClientSocket[] = "$Id$";
#endif

/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.3  2005/05/18 19:26:02  strieben
| Upgraded GSI library to newest available version.
|
| Revision 1.1.1.1  2004/08/31 15:53:00  simon
| Initial GSI import
|
|
+---------------------------------------------------------------------*/

#include <iostream>

#include "GSIClientSocket.hh"

namespace GridlibSocketInterface
{

ClientSocket::ClientSocket (const std::string& host,
                            const int32 port)
{
  
  if ( ! Socket::create() )
    {
        Socket::cleanup();
    }

  if ( ! Socket::connect ( host, port ) )
    {
        std::cout << "Wixenkacke! Connect geht nicht!";
        
        Socket::cleanup();
    }

}
 
}
