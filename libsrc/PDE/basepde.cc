#include <fstream>
#include <iostream>
#include <string>

#include "basepde.hh"
#include "conffile.hh"
 
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
  as_sysid_ = as_sysid;
}

void BasePDE::ReadBCs(const std::string eq)
{
#ifdef TRACE
  (*trace) << " entering BasePDE::ReadBCs " << std::endl;
#endif

  conf->getliststr("homogeneous_dirichlet",bcs_hd_,eq); 
  conf->getliststr("inhomogeneous_dirichlet",bcs_id_,eq);

  Integer i;

  val_id_.resize(bcs_id_.size());

  for(i=0; i<bcs_id_.size(); i++)
    conf->get(bcs_id_[i],val_id_[i],eq,"bc_conditions","inhomogeneous_dirichlet"); 
}

Integer BasePDE::GetNumRestraints(BCs* ptBCs, const Integer level)
{
#ifdef TRACE
  (*trace) << "entering BasePDE::GetNumRestraints" << std::endl;
#endif
    
  Integer res=0;
  Integer i;
  for (i=0; i<bcs_hd_.size(); i++)
    {
      res+=ptBCs->GetNumNodesLevel(bcs_hd_[i],level);
    }

  for (i=0; i<bcs_id_.size(); i++)
    {
      res+=ptBCs->GetNumNodesLevel(bcs_id_[i],level);
    }

  return res;
}

BasePDE::~BasePDE()
{
;
}

} // end of namespace
