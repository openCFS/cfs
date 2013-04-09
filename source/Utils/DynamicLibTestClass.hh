#ifndef DYNAMIC_LIB_TEST_H
#define DYNAMIC_LIB_TEST_H

#include "Utils/DynamicLibTestClassIface.hh"

// Test example for loading classes from a DLL.

namespace CoupledField
{

  class DynamicLibTestClass : public DynamicLibTestClassIface
  {
  public:
    DynamicLibTestClass(void (*delObj)(void*));
    virtual ~DynamicLibTestClass(void);

    virtual void TestIt(std::string& outStr,
                        const std::string& inStr);
  };
  
}

#endif
