#include <fstream>
#include <iostream>
#include <string>

#include "basedriver.hh"
#include <DataInOut/writeresults.hh>
namespace CoupledField
{

BaseDriver :: BaseDriver(Domain * adomain)
{
#ifdef TRACE
  (*trace) << "entering BaseDriver::BaseDriver" << std::endl;
#endif
  ptdomain_ = adomain;

  nummeshes_=0;
    // read info Should we save first,second derivatives or not
//  ptdomain->GetInFile()->ReadOutputOptions(SaveDer1,SaveDer2);

}

BaseDriver :: ~BaseDriver()
{
#ifdef TRACE
  (*trace) << "entering BaseDriver::~BaseDriver" << std::endl;
#endif

}

void BaseDriver :: PrintSeqMeshes()
{
#ifdef TRACE
  (*trace) << "entering TransientDriver :: PrintSeqMeshes" << std::endl;
#endif

   if (nummeshes_) {
     ptMeshes_->OpenFile(nummeshes_);
   }

   Integer level=0;
   ptMeshes_->WriteGrid(level);

   nummeshes_++;
}

}
