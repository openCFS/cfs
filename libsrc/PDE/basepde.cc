#include <fstream>
#include <iostream>
#include <string>

#include "basepde.hh"
#include <DataInOut/conffile.hh>
 
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

  //standard paramewter for solver
  eps_         = 1.0e-8;
  dampiter_    = 0.7;
  maxnumiter_  = 100;
  solvertype_  = RealDirect;
  precondtype_ = ID;
  numeqcoarse_ = 200;
  coarsealpha_ = 0.1;
 
}

void BasePDE::SpecifySolver(Integer &solvertype, Integer &precondtype, Double &eps, Double &dampiter, 
			    Integer &maxnumiter, Integer &numeqcoarse, Double &coarsealpha )
{
#ifdef TRACE
  (*trace) << "entering  BasePDE::SpecifySolver" << std::endl;
#endif

  //set values to default;
  eps = eps_;
  dampiter = dampiter_;
  maxnumiter  = maxnumiter_;
  solvertype  = solvertype_;
  precondtype = precondtype_;
  numeqcoarse = numeqcoarse_;
  coarsealpha = coarsealpha_;

  //if values are defined in conf-file, take these
  conf->ifget("eps",eps,pdename_); // relative accuracy in the precond. energy
  conf->ifget("dampiter",dampiter,pdename_); // damping parameter for Jacobi, SSOR
  conf->ifget("maxnumit",maxnumiter,pdename_); // maximal number of iterations
  conf->ifget("solvertype",solvertype,pdename_); // solver
  conf->ifget("precondtype", precondtype,pdename_); //preconditioner
  conf->ifget("numeqcoarse",numeqcoarse,pdename_); // number of equation for coarsing
  conf->ifget("coarsealpha",coarsealpha,pdename_); // coarsing parameter for AMG

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
  //val_id_=new Integer[bcs_id_.size()];

  for(i=0; i<bcs_id_.size(); i++)
    {
    conf->get(bcs_id_[i],val_id_[i],eq,"bc_conditions","inhomogeneous_dirichlet");
    }
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
#ifdef TRACE
  (*trace) << " entering BasePDE::~BasePDE() " << std::endl;
#endif

}

} // end of namespace
