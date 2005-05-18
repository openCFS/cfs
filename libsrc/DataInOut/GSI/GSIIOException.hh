// GSI::IOException class


#ifndef GSI_IOEXCEPTION
#define GSI_IOEXCEPTION

#include <string>

#include "GSITypeDefs.hh"
#include "GSIException.hh"

namespace GridlibSocketInterface
{

class IOException : public Exception
{
 public:
  IOException ( std::string s ) : Exception(s) {};
  ~IOException (){};

  int32 GetErrno() { return errno_; }
  void SetErrno(int32 err) { errno_ = err; }
  std::string GetErrorString();

 private:
   int32 errno_;
};
 
}

#endif //GSI_IOEXCEPTION
