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
  OutFile3d_=NULL;
  ptTimeFunc_ = aptTimeFunc;
  ptalgsys_   = aptalgsys;

}

BasePDE::BasePDE(AbstractAlgebraicSys * aptalgsys, Material * aMatFile, FileType * aInFile, WriteResults<Point3D> * aOutFile, TimeFunc *aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering BasePDE::BasePDE" << std::endl;
#endif

  MatFile_    = aMatFile;
  InFile_     = aInFile;
  OutFile3d_    = aOutFile;
  OutFile_=NULL;
  ptTimeFunc_ = aptTimeFunc;
  ptalgsys_   = aptalgsys;

}

void BasePDE::SetAlgSys_id(const Integer as_sysid)
{
 AS_sysid_ = as_sysid;
}

} // end of namespace
