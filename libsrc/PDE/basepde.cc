#include <fstream>
#include <iostream>
#include <string>

#include "basepde.hh"
 
namespace CoupledField
{

BasePDE::BasePDE(AbstractAlgebraicSys * aptalgsys, Material * aMatFile, FileType * aInFile, WriteResults * aOutFile, TimeFunc *aptTimeFunc)
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

void BasePDE::SetAlgSys_id(const Integer as_sysid)
{
 AS_sysid_ = as_sysid;
}

} // end of namespace
