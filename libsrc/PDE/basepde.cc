#include <fstream>
#include <iostream>
#include <string>

#include "basepde.hh"
 
namespace CoupledField
{

BasePDE::BasePDE(AbstractAlgebraicSys * aptalgsys, Material * aMatFile, FileType * aInFile, WriteResults<Point2D> * aOutFile, TimeFunc *aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering BasePDE::BasePDE" << std::endl;
#endif

  MatFile_    = aMatFile;
  InFile_     = aInFile;
  OutFile_    = aOutFile;
  ptTimeFunc_ = aptTimeFunc;
  ptalgsys_   = aptalgsys;

}

} // end of namespace
