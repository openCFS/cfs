#include <iostream>
#include <string>

#include "DynamicLibrary.hh"
#include "DynamicLibTestClass.hh"

using namespace CoupledField;

int main(int argc, char** argv) 
{
#ifndef _WIN32
  int dynLoadFlags = RTLD_NOW;
  std::string dllName = "libDynamicLibTestClass.so";
#else
  int dynLoadFlags = 0;
  std::string dllName = "libDynamicLibTestClass.dll";
#endif

  DynamicLibrary* dynLibrary = 
    DynamicLoader::loadObjectFile(dllName.c_str(), dynLoadFlags);

  if(!dynLibrary) 
  {
    std::cerr << "Couldn't load the dynamic library '" << dllName << "'.\n";
    return 1;
  }
  
  // Load a DynamicLibTestClass object
  DynamicLibTestClassIface* dynLibTestClass = 
    dynamic_cast<DynamicLibTestClassIface*>(dynLibrary->newObject("DynamicLibTestClass", 0, NULL));

  if(dynLibTestClass == NULL) {
    std::cerr << "Couldn't create DynamicLibTestClass object." << std::endl;
    return 1;
  }
  
  std::string output;
  
  dynLibTestClass->TestIt(output, "FRANZ");
  
  std::cout << "TestIt returned: " << output << std::endl;
  
  // Get rid of the binlibIface_ object
  dynLibTestClass->deleteSelf();
  dynLibTestClass = NULL;
  
  // Close the dynamic library
  delete dynLibrary;
  
  return 0;
  
}


