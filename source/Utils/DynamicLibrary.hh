#ifndef DYNAMICLIBRARY_H
#define DYNAMICLIBRARY_H

#include <dlfcn.h>

namespace CoupledField
{

  class DynamicObject;
  
  // The DynamicLibrary class.
  // This class will be created by the dynamic loader, 
  // and can be used to create objects out of a loaded library
  class DynamicLibrary
  {
  protected:
    // The handle to the shared library that was opened
    void* _objFile;
    // Since an instance of DynamicLibrary manages lifetime of an open 
    // library, it is important to make sure that the object isn't 
    // copied.
    DynamicLibrary(const DynamicLibrary&) {}
    DynamicLibrary& operator=(const DynamicLibrary&) {}
    
    // Creates a new library, with the object file handle that is passed 
    // in. Protected so that only the DynamicLoader can create an 
    // instance (since it is declared friend.
    DynamicLibrary(void* objFile);
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
