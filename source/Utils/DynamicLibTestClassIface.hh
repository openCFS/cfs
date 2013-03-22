#ifndef DYNAMIC_LIB_TEST_IFACE_H
#define DYNAMIC_LIB_TEST_IFACE_H

#include "Utils/DynamicObject.hh"

// Test example for loading classes from a DLL.

namespace CoupledField
{

  class DynamicLibTestClassIface : public DynamicObject
  {
  public:
    DynamicLibTestClassIface(void (*delObj)(void*));
    virtual ~DynamicLibTestClassIface(void);

    virtual void TestIt(std::string& outStr,
                        const std::string& inStr) = 0;
  };
  
}

#endif
