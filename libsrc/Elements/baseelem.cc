#include <iostream>
#include <fstream>
#include <string>

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

enum IntegrationType BaseElem::String2EnumIntegrationType(const Char * inttype)
{
#ifdef TRACE
  (*trace) << "entering BaseElem::String2EnumIntegrationType" << std::endl;
#endif

 enum IntegrationType result=null;
 if (!std::strcmp(inttype,"GaussOrder2")) result=GaussOrder2;
 if (!std::strcmp(inttype,"GaussOrder3")) result=GaussOrder3;
 if (!std::strcmp(inttype,"GaussOrder4")) result=GaussOrder4;
 if (!std::strcmp(inttype,"GaussOrder5")) result=GaussOrder5;
 if (!std::strcmp(inttype,"GaussOrder7")) result=GaussOrder7;

 if (result==null) Error("Check your config file. Your integration type is wrong",__FILE__,__LINE__);

 return result;
}

}
