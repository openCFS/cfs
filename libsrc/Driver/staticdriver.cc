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
  StaticDriver::StaticDriver(Domain * adomain, 
			     Integer stepOffset,
			     Double timeOffset,
			     std::string driverTag,
			     Boolean isPartOfSequence) 
    : SingleDriver(adomain, stepOffset, timeOffset, 
		   driverTag, isPartOfSequence) {
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
    
    // if PDEs are not set explicitly, 
    // get all from the domain
    if (pdes_.GetSize() == 0) {
      pdes_.Resize(ptdomain_->GetNumPDE());
      for (Integer iPDE=0; iPDE<ptdomain_->GetNumPDE(); iPDE++)
	pdes_[iPDE] = ptdomain_->GetPDE(iPDE);
    }
    
    // initialize pdes only, if this driver
    // is not part of multiSequence driver
    StdVector<std::string> tags;
    if (! isPartOfSequence_) {
      tags.Resize(pdes_.GetSize());
      tags.Init("anyTag");
      ptdomain_->InitPDEs(1,tags);
      Info->StartProgress ("Starting to solve problem", FALSE);
    }
    

    // Initialize 'TimeStepping'
    const Integer nstep = 1;
    Double  steptime = 0.0;
    Boolean reset = FALSE;
    
    // Solve problem
    if (pdes_.GetSize() <= 1) {

      // branch for single PDE
      pdes_[pdenumber]->WriteGeneralPDEdefines();
      pdes_[pdenumber]->SolveStepStatic(nstep, steptime, level,
						    reset);   
      pdes_[pdenumber]->PostProcess(level);

      // if multiSequence is performed, the ms-driver
      // writes out the grid one time
      if (! isPartOfSequence_)
	ptdomain_->PrintGrid(level);
      
      pdes_[pdenumber]->WriteResultsInFile(stepOffset_, timeOffset_);
    }
    else {
      // branch for coupled PDEs
      ptdomain_->GetCoupledPDE()->WriteGeneralPDEdefines();

      // define which PDEs participate in solving process
      ptdomain_->GetCoupledPDE()->DefineSolvingPDEs(pdes_);
      ptdomain_->GetCoupledPDE()->SolveStepStatic(nstep, steptime, level,
						  reset);
      // if multiSe
      // quence is performed, the ms-driver
      // writes out the grid one time
      if (! isPartOfSequence_)
	ptdomain_->PrintGrid(level);

      ptdomain_->GetCoupledPDE()->WriteResultsInFile(stepOffset_, timeOffset_);
    }
  }

} // end of namespace
