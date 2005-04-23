#include <fstream>
#include <iostream>
#include <string>

#include "solveStepElec.hh"
#include "assemble.hh"

#include "Utils/nodestoresol.hh"
#include "PDE/StdPDE.hh"

namespace CoupledField {

  SolveStepElec::SolveStepElec(StdPDE& apde) : StdSolveStep(apde)
  {

    ENTER_FCN( "SolveStepElec::SolveStepElec" );
  }


  // ======================================================
  // Solve Step Static SECTION  
  // ======================================================

void SolveStepElec:: PreStepStatic(const Integer kstep, const Double asteptime,
			     const Integer level, const Boolean reset)
{
  ENTER_FCN( "SolveStepElec::PreStepStatic" );

  if (isIterCoupled_)     
    algsys_->InitSol();
  
  if (geoUpdate_)
    {
      algsys_->InitRHS();
      algsys_->InitSol();
      algsys_->InitMatrix();
      assemble_->SetReassemble();   
    }
}

void SolveStepElec::PostStepStatic(const Integer kstep, const Double asteptime,
			     const Integer level)
{
  ENTER_FCN( "SolveStepElec::PostStepStatic" );

  if (isIterCoupled_)
    (*iterCoupledCounter_)++;


#ifdef ADAPTGRID
  if (flags->CalcErrorMap_)
    {
      Double         totalErr;
      ElemStoreSol<Double>  Sol_Mesh;
      Vector<Double> solVec;
      
      ptError_=new SpaceErrorEstimator();

      ptError_->Init(this);
      
      sol_->TransformNodeSolution(Sol_Mesh,sol_,eqnData_,ptgrid_);

      sol_->GetCompleteVector(solVec);

      //  solVec.Resize(sol_.size());
      //       int i;
      //       for (i=0; i<sol_.size(); i++)
      // 	solVec[i]=sol_[0][i];

      ptError_->CalcErrorMap(solVec,subdoms_,ptgrid_,errorMap_,totalErr,level);
      
      std::cout << " total error of calculation:: " << totalErr << std::endl;
      *data << errorMap_ << std::endl;
    }
#endif
}



  //   Default Destructor
  // **********************
  SolveStepElec::~SolveStepElec() {

    ENTER_FCN( "SolveStepElec::~SolveStepElec" );
 
  }

} // end of namespace

