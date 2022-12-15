#ifndef DYNAMICLIBRARY_H
#define DYNAMICLIBRARY_H

#ifndef _WIN32
#include <dlfcn.h>
#else
#define NOMINMAX
#include <windows.h>
#endif

namespace CoupledField
{

  class DynamicObject;
  
  // The DynamicLibrary class.
  // This class will be created by the dynamic loader, 
  // and can be used to create objects out of a loaded library
  class DynamicLibrary
  {
  protected:
#ifndef _WIN32
    // The handle to the shared library that was opened
    void* objFile_;
#else
    HMODULE hmod_;
#endif

    // Since an instance of DynamicLibrary manages lifetime of an open 
    // library, it is important to make sure that the object isn't 
    // copied.
    DynamicLibrary(const DynamicLibrary&) {}
    DynamicLibrary& operator=(const DynamicLibrary&) { return *this; }
    
#ifndef _WIN32
    // Creates a new library, with the object file handle that is passed 
    // in. Protected so that only the DynamicLoader can create an 
    // instance (since it is declared friend.
    DynamicLibrary(void* objFile);
#else
    // Creates a new library, with the object file handle that is passed 
    // in. Protected so that only the DynamicLoader can create an 
    // instance (since it is declared friend.
    DynamicLibrary(HMODULE hmod);
#endif
  public:
    // Destructor, closes the open shared library
    ~DynamicLibrary(void);
    
    // Creates a new instance of the named class, or returns NULL is the 
    // class isn't found. The array of void*'s takes the arguments for   
    // the class's constructor, with argc giving the number of arguments.
    DynamicObject* newObject(const char* name, int argc, void** argv);
    
    // Declare the DynamicLoader as a friend, so it can
    // instantiate DynamicLibrary.
    friend class DynamicLoader;
    
  };

  // The dynamic loader class, used for loading DynamicLibraries.
  class DynamicLoader
  {
  public:
    // Loads a DynamicLibrary, given the shared library file
    // "file", with the dlopen flags supplied.
    static DynamicLibrary* loadObjectFile(const char* file, int flags);
  };

}

#endif
