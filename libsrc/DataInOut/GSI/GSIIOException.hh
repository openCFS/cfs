// IOException class


#ifndef GSIIOException_class
#define GSIIOException_class

#include <string>

#include "GSIException.hh"

class GSIIOException : public GSIException
{
 public:
  GSIIOException ( std::string s ) : GSIException(s) {};
  ~GSIIOException (){};
};

#endif
