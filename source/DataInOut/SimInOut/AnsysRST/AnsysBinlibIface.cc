//AnsysBinlibIface.cpp
#include "AnsysBinlibIface.hh"

namespace CoupledField
{

  AnsysBinlibIface::AnsysBinlibIface(void (*delObj)(void*))
    : DynamicObject(delObj)
  {
  }

  AnsysBinlibIface::~AnsysBinlibIface(void)
  {
  }

}
