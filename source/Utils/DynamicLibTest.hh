#ifndef DYNAMIC_LIB_TEST_H
#define DYNAMIC_LIB_TEST_H

#include "Utils/DynamicObject.hh"

// Test example for loading classes from a DLL.

namespace CoupledField
{

  class DynamicLibTestClass : public DynamicObject
  {
  public:
    DynamicLibTestClass(void (*delObj)(void*));
    virtual ~DynamicLibTestClass(void);

    virtual void TestIt(std::string& outStr,
                        const std::string& inStr);
  };
  
}

#endif
