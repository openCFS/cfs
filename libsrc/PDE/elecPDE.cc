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

ElecPDE::ElecPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc, FileType *aptFileType, 
		 WriteResults *aptOut)
:BasePDE(aptgrid, aptbcs, aptFileType, aptOut, aptTimeFunc)
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


void ElecPDE::DiscreteParamsPDE()
{
#ifdef TRACE
  (*trace) << "entering ElecPDE::DiscreteParamsPDE" << std::endl;
#endif

  MatrixType_   = RSCALAR;
  GraphType_    = NODEGRAPH; 
  SystemMatrix_ = TRUE;

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

#ifdef DEBUG
  std::string matFileName = "solMat";
  OutFile_->WriteSolMatrix(ptgrid_, level, sol_, matFileName);
#endif
}

void ElecPDE::WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering ElecPDE::WriteResultsInFile" << std::endl;
#endif

  Integer step=0;
  Double time=0;

  // write results
  OutFile_->WriteSolution(sol_, step, time, "electric_potential");

}

ElecPDE::~ElecPDE()
{
 ;
}

} // end of namespace


