#include <fstream>
#include <iostream>
#include <string>

#include "basepde.hh"
 
namespace CoupledField
{

BasePDE::BasePDE(Material * aMatFile, FileType * aInFile, WriteResults<Point2D> * aOutFile, TimeFunc *aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering BasePDE::BasePDE" << std::endl;
#endif

  MatFile    = aMatFile;
  InFile     = aInFile;
  OutFile    = aOutFile;
  ptTimeFunc = aptTimeFunc;
}

void BasePDE::SetStepData()
{
#ifdef TRACE
  (*trace) << "entering BasePDE::SetStepData" << std::endl;
#endif

  SetMatrixFactors();
}

} // end of namespace
