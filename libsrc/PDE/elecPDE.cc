#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "elecPDE.hh"
#include <DataInOut/Unverg/outUnverg.hh>
#include <Forms/forms_header.hh>
#include <Estimator/spaceerror.hh>
 
namespace CoupledField
{

ElecPDE::ElecPDE(Grid * aptgrid, BCs *aptbcs, Material *ptMaterial, TimeFunc *aptTimeFunc, 
		     FileType *aptFileType, WriteResults *aptOut)
:BasePDE(aptgrid,aptbcs,ptMaterial,aptFileType,aptOut,aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering ElecPDE::ElecPDE " << std::endl;
#endif

  dofspernode_=1;

  size_=ptgrid_->GetMaxnumnodes(0);
  sol_.Resize(size_);
  sol_.Init(0);

  SetMatrixFactors();
}


void ElecPDE::SetMatrixFactors()
{
#ifdef TRACE
  (*trace) << "entering ElecPDE::SetMatrixFactors" << std::endl;
#endif
 
  matrix_factor_[0] = 0.0;  // factor for stiffness matrix
  matrix_factor_[1] = 0.0;  // factor for damping matrix
  matrix_factor_[2] = 0.0;  // factor for convection matrix
  matrix_factor_[3] = 0.0;  // factor for mass matrix
}

void ElecPDE::SpecifyMatrices(Integer &matrixtype, Integer * matrixsystype, Integer &graphtype, 
				Integer &numdofpernode, Integer &numdirichlets, Integer &numconstraints)
{
#ifdef TRACE
  (*trace) << "entering ElecPDE::SpecifyMatrices" << std::endl;
#endif

  matrixtype = RSCALAR; 
  graphtype  = NODEGRAPH; 

  matrixsystype[0] = SYSTEM;      // memory for the system matrix
  matrixsystype[1] = 0;           // memory for the stiffness matrix
  matrixsystype[2] = 0;           // memory for the damping matrix
  matrixsystype[3] = 0;           // memory for the convection matrix
  matrixsystype[4] = 0;           // memory for the mass matrix


  numdofpernode  = dofspernode_;
  numdirichlets  = GetNumRestraints(actlevel_);
  numconstraints = 0;
}

void ElecPDE::SetAlgSys(const Integer as_sysid)
{
  as_sysid_ = as_sysid;

  //allocate according algebraic system
  algsys_ = new StandardSystem();

  //set solver parameters  
  SetSolverParameters();

  //set the graph type used for the system matrices
  Integer numnode = ptgrid_->GetMaxnumnodes(actlevel_);
  SetupMatrixGraph(numnode,NODEGRAPH);

  //allocate the necessary matrices as well as solver and preconditioner
  CreateMatrices_Solver();
}


void ElecPDE::SolveStepStatic(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering ElecPDE::SolveStepStatic" << std::endl;
#endif

  Integer update = 0;
  Integer job = 1;

  Double * ptsol;

  //compute and assemble element matrices
  SetupMatrices(level);

  //account for bcs
  SetBCs(level,update,0);

  algsys_->CalcPrecond(job);
  algsys_->Solve();

  ptsol = algsys_->GetSolutionVal();

  // save solution
  Vector<Double> transsol(ptgrid_->GetMaxnumnodes(level), ptsol);
  sol_=transsol;
}

void ElecPDE::CalcCoeff(Vector<Double> & coeff)
{
  if (!MatFile_) Error("You didn't specialize material file. Use option -m");

  coeff.Resize(subdoms_.size());

  Integer i, matnum;
  for (i=0; i<subdoms_.size(); i++)
    {
      conf->get(subdoms_[i],matnum,"list_subdomains");

      Double dielectr;
      MatFile_->ReadDielectricTerms(dielectr,matnum); 

      coeff[i]=dielectr;
    }
}

void ElecPDE::WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering ElecPDE::WriteResultsInFile" << std::endl;
#endif

  Integer step=0;
  Double time=0;

  // write results
  if (OutFile_->IsGMV())
    OutFile_->WriteSolution(sol_,step,time,"electric_potential");
  else
    OutFile_->WriteSolution(sol_,step,time,"electric potential");

}

ElecPDE::~ElecPDE()
{
 ;
}

} // end of namespace


