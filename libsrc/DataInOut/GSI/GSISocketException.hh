// SocketException class


#ifndef GSISocketException_class
#define GSISocketException_class

#include <string>

#include "GSIException.hh"

class GSISocketException : public GSIException
{
 public:
  GSISocketException ( std::string s ) : GSIException(s) {};
  ~GSISocketException (){};
};

#endif
