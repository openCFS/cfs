#include <fstream>
#include <iostream>
#include <string>

#include "staticdriver.hh"
#include "DataInOut/GMV/outGMV.hh"
#include "CoupledPDE/basecoupledpde.hh"
#include "General/environment.hh"
#include "PDE/basePDE.hh"

namespace CoupledField {

  // ***************
  //   Constructor
  // ***************
  StaticDriver::StaticDriver(Domain * adomain) : SingleDriver(adomain) {
    ENTER_FCN( "StaticDriver::StaticDriver" );
  }


  // **********************
  //   Default destructor
  // **********************
  StaticDriver::~StaticDriver() {
    ENTER_FCN( "StaticDriver::~StaticDriver" );
  }


  // *****************
  //   Solve problem
  // *****************
  void StaticDriver::SolveProblem() {
    ENTER_FCN( "StaticDriver::SolveProblem" );

    Integer level=0;
    Integer pdenumber  = 0;

    if (PrintGridOnly) {
      ptdomain_->PrintGrid(level);
      exit(0);
    }

    const Integer nstep = 1;
    Double  steptime = 0;
    Boolean reset = FALSE;
  
    if (ptdomain_->GetNumPDE() <= 1) {
      ptdomain_->GetPDE(pdenumber)->WriteGeneralPDEdefines();
      ptdomain_->GetPDE(pdenumber)->SolveStepStatic(nstep, steptime, level,
						    reset);   
      ptdomain_->GetPDE(pdenumber)->PostProcess(level);
      ptdomain_->PrintGrid(level);
      ptdomain_->GetPDE(pdenumber)->WriteResultsInFile();
    }
    else {
      ptdomain_->GetCoupledPDE()->WriteGeneralPDEdefines();
      ptdomain_->GetCoupledPDE()->SolveStepStatic(nstep, steptime, level,
						  reset);
      ptdomain_->PrintGrid(level);
      ptdomain_->GetCoupledPDE()->WriteResultsInFile();
    }
  }

} // end of namespace
