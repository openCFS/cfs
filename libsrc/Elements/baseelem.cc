#include <iostream>
#include <fstream>

#include "baseelem.hh"
 
namespace CoupledField
{

 BaseElem :: BaseElem()
{
#ifdef TRACE
  (*trace) << "entering BaseElem::BaseElem" << std::endl;
#endif
  ;
}
 
BaseElem :: ~BaseElem()
{
#ifdef TRACE
  (*trace) << "entering BaseElem::~BaseElem" << std::endl;
#endif
 
  ;
}

}
