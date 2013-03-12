#include <cstring>

#include <iostream>

#include "DynamicLibTestClass.hh"

namespace CoupledField
{

  DynamicLibTestClass::DynamicLibTestClass(void (*delObj)(void*))
    : DynamicLibTestClassIface(delObj) 
  {
    std::cout << "Creating DynamicLibTestClass..." << std::endl;
  }
  
  DynamicLibTestClass::~DynamicLibTestClass(void)
  {
    std::cout << "Destructing DynamicLibTestClass..." << std::endl;
  }
  

  void DynamicLibTestClass::TestIt(std::string& outStr,
                                   const std::string& inStr)
  {
    outStr = "GEORG";
    
    std::cout << "Got inStr: " << inStr << std::endl;
    std::cout << "Returning outStr: " << outStr << std::endl;
  }

  // The Dynamic library should also contain the following external C 
  // function definitions
  // Factory methods for generating and deleting objects
  // from this DLL
  extern "C"
  {
    void deleteObject(void* obj) {
      delete reinterpret_cast<DynamicObject*>(obj);
    }

    void* loadObject(const char* name, int argc, void** argv) {
      if(std::strncmp(name,
                      "DynamicLibTestClass",
                      strlen(name) < 19 
                 ? strlen(name) 
                 : 19) == 0) {
        return new DynamicLibTestClass(deleteObject);
      }

      return NULL;
    }
  }
  
}
