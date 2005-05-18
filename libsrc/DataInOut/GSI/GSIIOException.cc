#ifdef USE_RCSID
static const char RCSid_GSIIOException[] = "$Id$";
#endif

/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.2  2005/05/18 19:26:03  strieben
| Upgraded GSI library to newest available version.
|
| Revision 1.1.1.1  2004/08/31 15:53:00  simon
| Initial GSI import
|
|
+---------------------------------------------------------------------*/

#include <cstring>

#include "GSIIOException.hh"

namespace GridlibSocketInterface
{

  std::string IOException::GetErrorString()
  {
    char buf[256];
	strncpy(buf, strerror(errno_), 256);
    return buf;
  }
  
}
