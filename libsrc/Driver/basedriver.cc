#include <fstream>
#include <iostream>
#include <string>

#include "basedriver.hh"

namespace CoupledField
{

template<class Dim>
BaseDriver<Dim> :: BaseDriver(Domain<Dim> * adomain)
{
#ifdef TRACE
  (*trace) << "entering BaseDriver::BaseDriver" << std::endl;
#endif
  ptdomain_ = adomain;

    // read info Should we save first,second derivatives or not
//  ptdomain->GetInFile()->ReadOutputOptions(SaveDer1,SaveDer2);

}

template<class Dim>
BaseDriver<Dim> :: ~BaseDriver()
{
#ifdef TRACE
  (*trace) << "entering BaseDriver::~BaseDriver" << std::endl;
#endif

}

}
