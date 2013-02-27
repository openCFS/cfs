// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "IterCoupledPDE.hh"


#include "PDE/SinglePDE.hh"
#include "Driver/SolveSteps/IterSolveStep.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/CoefFunction/CoefFunctionAccumulator.hh"

namespace CoupledField
{

  IterCoupledPDE::IterCoupledPDE(StdVector<StdPDE*> & PDEs,
                                 PtrParamNode paramNode) 
  : BasePDE( paramNode ) {

    PDEs_       = PDEs;
    myParam_ = paramNode;
    infoNode_ = info->Get("PDE")->Get("iterCoupledPDE", ParamNode::APPEND);

    // Loop over single PDEs
    for( UInt i = 0; i < PDEs_.GetSize(); ++i ) {
      singlePDEs_.Push_back(dynamic_cast<SinglePDE*>(PDEs_[i]));
    }
    
    numPDEs_ = singlePDEs_.GetSize();
    solveStep_ = NULL;

    // Concatenate PDE name strings for output in info-file
    pdename_ = "CoupledPDE: ";
    for (UInt actPDE=0; actPDE < PDEs.GetSize()-1; actPDE++)
      pdename_ += PDEs[actPDE] -> GetName() + "+";
    pdename_ += PDEs[PDEs.GetSize()-1] -> GetName();

    
    // Create IterSolveStep instance
    IterSolveStep * solveStep =new IterSolveStep( *this, myParam_, infoNode_ ); 
    solveStep->Init();
    solveStep_ = solveStep;
  }


  IterCoupledPDE::~IterCoupledPDE() {
    // delete solveStep-object
    delete solveStep_;

  }

  PtrCoefFct IterCoupledPDE::GetCouplingCoefFct( SolutionType type,
                                                 shared_ptr<EntityList>  list,
                                                 const std::string& pdeName,
                                                 bool& updateGeo ) {
    
    // directly pass the query to the IterSolveStep instance
    return solveStep_->GetCouplingCoefFct( type, list, pdeName, updateGeo );
    
  }

  void IterCoupledPDE::WriteRestart() {

    for (UInt actPDE=0; actPDE < PDEs_.GetSize(); actPDE++)
      PDEs_[actPDE]->WriteRestart( );
  }


  void IterCoupledPDE::ReadRestart(UInt &startStep)  {

    for (UInt actPDE=0; actPDE < PDEs_.GetSize(); actPDE++)
      PDEs_[actPDE]->ReadRestart(startStep);
  }


  void IterCoupledPDE::WriteResultsInFile(const UInt kstep,
                                          const Double asteptime ) {

    for (UInt i=0; i<PDEs_.GetSize(); i++)
      PDEs_[i]->WriteResultsInFile(kstep, asteptime );
  }


  void IterCoupledPDE::ToInfo(PtrParamNode base) {
    std::cerr << "only test\n";
  }


  void IterCoupledPDE::WriteGeneralPDEdefines()
  {

    for (UInt i=0; i<PDEs_.GetSize(); i++)
      PDEs_[i]->WriteGeneralPDEdefines();
  }


  BaseSolveStep * IterCoupledPDE::GetSolveStep()
  {
    return solveStep_;
  }
  
  void IterCoupledPDE::UpdateToSolStrategy() {

    for (UInt i=0; i<PDEs_.GetSize(); i++)
      PDEs_[i]->UpdateToSolStrategy();
  }


} // end of namespace
