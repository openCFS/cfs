#include <iostream>
#include <cstdlib>

#include "DynamicLibrary.hh"

namespace CoupledField
{
  DynamicLibrary::DynamicLibrary(void* objFile)
    : _objFile(objFile)
  {
  }
  
  DynamicLibrary::~DynamicLibrary(void)
  {
    dlclose(_objFile);
  }

  DynamicObject*
  DynamicLibrary::newObject(const char* name, int argc, void** argv)
  {
    // If there is no valid library, return null
    if(_objFile == 0) {
      return 0;
    }

    typedef void* (*FuncPtrType)(const char*, int, void**); // a function pointer type
    
    // Get the loadObject() function.  If it doesn't exist, return NULL.
    // http://www.mr-edd.co.uk/blog/supressing_gcc_warnings
#ifdef __GNUC__
__extension__
#endif
    FuncPtrType loadSym = reinterpret_cast<FuncPtrType>( dlsym(_objFile, "loadObject") );
    if(loadSym == 0) {
      return 0;
    }

    // Load a new instance of the requested class, and return it
    void* obj = 0;
    obj = loadSym(name, argc, argv);
    return reinterpret_cast<DynamicObject*>(obj);
  }
  
  DynamicLibrary*
  DynamicLoader::loadObjectFile(const char* file, int flags)
  {
    void* objFile = dlopen(file, flags);
    if(objFile == 0) {
      return 0;
    }
    return new DynamicLibrary(objFile);
  }
}
