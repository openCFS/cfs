#ifdef USE_RCSID
static const char RCSid_GSIIOException[] = "$Id$";
#endif

/*----------------------------------------------------------------------
  |
  |
  | $Log$
  | Revision 1.3  2005/05/23 22:15:31  ahauck
  | Corrected indentation of all cc/hh-files according to xemacs' 'gnu' style
  | and according to our coding rules, so all TABs were replaced by spaces.
  |
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
