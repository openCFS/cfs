// Exception class

#ifndef GSIException_class
#define GSIException_class

#include <string>

class GSIException
{
 public:
  GSIException ( std::string s) :  m_s(s) { };
  ~GSIException (){};

  std::string description() { return m_s; }

 protected:

  std::string m_s;
};

#endif
