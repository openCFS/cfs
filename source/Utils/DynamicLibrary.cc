#include <iostream>
#include <cstdlib>

#include "DynamicLibrary.hh"

namespace CoupledField
{
#ifndef _WIN32
  DynamicLibrary::DynamicLibrary(void* objFile)
    : objFile_(objFile)
  {
  }
#else
  DynamicLibrary::DynamicLibrary(HMODULE hmod)
    : hmod_(hmod)
  {
  }
#endif
  
  DynamicLibrary::~DynamicLibrary(void)
  {
#ifndef _WIN32
    dlclose(objFile_);
#else
    FreeLibrary(hmod_);
#endif
  }

  DynamicObject*
  DynamicLibrary::newObject(const char* name, int argc, void** argv)
  {
#ifndef _WIN32
    // If there is no valid library, return null
    if(objFile_ == 0) {
      return 0;
    }

    typedef void* (*FuncPtrType)(const char*, int, void**); // a function pointer type

    // Get the loadObject() function.  If it doesn't exist, return NULL.
    // http://www.mr-edd.co.uk/blog/supressing_gcc_warnings
#ifdef __GNUC__
__extension__
#endif
    FuncPtrType loadSym = reinterpret_cast<FuncPtrType>( dlsym(objFile_, "loadObject") );
#else
    // If there is no valid library, return null
    if(hmod_ == 0) {
      return 0;
    }

    typedef void* (__cdecl *FuncPtrType)(const char*, int, void**); // a function pointer type
    FuncPtrType loadSym = (FuncPtrType) GetProcAddress(hmod_, "loadObject");
#endif

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
#ifndef _WIN32
    void* objFile = dlopen(file, flags);
    if(objFile == 0) {
      return 0;
    }
    return new DynamicLibrary(objFile);
#else
    HMODULE hmod = LoadLibrary(file);
    if (hmod == NULL) 
    {
      return 0;
    }
    return new DynamicLibrary(hmod);
#endif
  }
}
