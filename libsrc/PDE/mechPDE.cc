#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include <Estimator/actimeerror.hh>
#include <DataInOut/Unverg/outUnverg.hh>
#include <DataInOut/GMV/outGMV.hh>
#include <Forms/forms_header.hh>
#include <Estimator/spaceerror.hh>

#include "mechPDE.hh" 

namespace CoupledField
{

MechPDE::MechPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc, FileType *aptFileType, 
		 WriteResults *aptOut)
:BasePDE(aptgrid, aptbcs, aptFileType, aptOut, aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::MechPDE " << std::endl;
#endif

  SetMatrixFactors();
}

void MechPDE::DiscreteParamsPDE()
{
#ifdef TRACE
  (*trace) << "entering MechPDE ::DiscreteParamsPDE" << std::endl;
#endif

  MatrixType_   = RBLOCK;
  GraphType_    = NODEGRAPH; 
  SystemMatrix_ = TRUE;

}


void MechPDE::SetMatrixFactors()
{
#ifdef TRACE
  (*trace) << "entering MechPDE::SetMatrixFactors" << std::endl;
#endif

  matrix_factor_[0] = 0.0;  // factor for stiffness matrix
  matrix_factor_[1] = 0.0;  // factor for damping matrix
  matrix_factor_[2] = 0.0;  // factor for convection matrix
  matrix_factor_[3] = 0.0;  // factor for mass matrix
}


void MechPDE::SolveStepStatic(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::SolveStepStatic" << std::endl;
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
  Vector<Double> transsol(ptgrid_->GetMaxnumnodes(level)*dofspernode_, ptsol);
  disp_=transsol;
}


void MechPDE::WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering MechPDE::WriteResultsInFile" << std::endl;
#endif

  Integer laststepcalc=0;
  Double  lasttimecalc=0;

  if (OutFile_->IsGMV())
    {
      OutFile_->WriteSolution(disp_,laststepcalc,lasttimecalc,"displacement", dofspernode_);
    }
  else
    {
      OutFile_->WriteSolution(disp_,laststepcalc,lasttimecalc,"displacement",  dofspernode_);
    }
}
}
