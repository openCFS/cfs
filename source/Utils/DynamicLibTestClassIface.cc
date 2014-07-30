#include <cstring>

#include <iostream>

#include "DynamicLibTestClassIface.hh"

namespace CoupledField
{
  DynamicLibTestClassIface::DynamicLibTestClassIface(void (*delObj)(void*))
    : DynamicObject(delObj) 
  {
    std::cout << "Creating DynamicLibTestClassIface..." << std::endl;
  }
  
  DynamicLibTestClassIface::~DynamicLibTestClassIface(void)
  {
    std::cout << "Destructing DynamicLibTestClassIface..." << std::endl;
  }
}
