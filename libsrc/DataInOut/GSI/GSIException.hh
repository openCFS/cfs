// GSI::Exception class

#ifndef GSI_EXCEPTION
#define GSI_EXCEPTION

#include <string>

#include "GSITypeDefs.hh"

namespace GridlibSocketInterface
{

class Exception
{
 public:
  Exception ( std::string s) :  m_s(s) { };
  ~Exception (){};

  std::string GetDescription() { return m_s; }
  void SetDescription(const std::string& descr) { m_s = descr; }

 protected:

  std::string m_s;
};
 
}

#endif //GSI_EXCEPTION
