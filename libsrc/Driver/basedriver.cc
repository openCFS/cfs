#include <fstream>
#include <iostream>
#include <string>

#include "basedriver.hh"

namespace CoupledField
{

BaseDriver :: BaseDriver(Domain * adomain)
{
#ifdef TRACE
  (*trace) << "entering BaseDriver::BaseDriver" << std::endl;
#endif
  ptdomain_ = adomain;

    // read info Should we save first,second derivatives or not
//  ptdomain->GetInFile()->ReadOutputOptions(SaveDer1,SaveDer2);

}

BaseDriver :: ~BaseDriver()
{
#ifdef TRACE
  (*trace) << "entering BaseDriver::~BaseDriver" << std::endl;
#endif

}

}
