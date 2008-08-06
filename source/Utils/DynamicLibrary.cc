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
    // Get the loadObject() function.  If it doesn't exist, return NULL.
    void* loadSym = dlsym(_objFile, "loadObject");
    if(loadSym == 0) {
      return 0;
    }
    // Load a new instance of the requested class, and return it
    void* obj = 
      ((void* (*)(const char*, int, void**))(loadSym))(name, argc, argv);
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
