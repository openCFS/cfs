#include <cstring>

#include "GSIIOException.hh"

namespace GridlibSocketInterface
{

  std::string IOException::GetErrorString()
  {
    char buf[256];
    strerror_r(errno_, buf, 256);
    return buf;
  }
  
}
