#include <iostream>
#include <fstream>

#include "getetrahedral.hh"

namespace CoupledField
{

void GeTetrahedral::SetIntPoints()
{
#ifdef TRACE
  (*trace) << "entering GeTetrahedral::SetIntPoints" << std::endl;
#endif

  switch(IntegType)
   {
   case GaussOrder3:
      break;

   default:
     std::cerr << "Integration type " << IntegType
           << " is not implemented \n" << std::endl; exit(-1);
   }
}

void GeTetrahedral::SetTransformFncAtIntPoints()
{
#ifdef TRACE
  (*trace) << "entering GeTetrahedral::SetTransformFncAtIntPoints" << std::endl;
#endif
  

}

void GeTetrahedral::SetDerTransformFncAtIntPoints()
{
#ifdef TRACE
  (*trace) << "entering GeTetrahedral::SetDShapeFnc" << std::endl;
#endif
   
}

} // end of namespace

