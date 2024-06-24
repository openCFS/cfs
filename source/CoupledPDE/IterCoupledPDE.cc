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

  IterCoupledPDE::IterCoupledPDE( StdVector<SinglePDE*>& singlePdes,
                                  StdVector<DirectCoupledPDE*>& cplPdes,
                                  PtrParamNode paramNode,
                                  PtrParamNode infoNode,
                                  shared_ptr<SimState> simState,
                                  Domain* domain) 
  : BasePDE( paramNode, infoNode, simState, domain ) {

    
    myParam_ = paramNode;
    infoNode_ = infoNode->Get("PDE")->Get("iterCoupledPDE", ParamNode::APPEND);

    // Initially we only store the single PDEs and the DirectCoupledPDEs.
    singlePDEs_ = singlePdes;
    coupledPDEs_ = cplPdes;
    
    solveStep_ = NULL;

    // Concatenate PDE name strings for output in info-file
    pdename_ = "CoupledPDE: ";
    for (UInt actPDE=0; actPDE < singlePdes.GetSize()-1; actPDE++)
      pdename_ += singlePdes[actPDE] -> GetName() + "+";
    pdename_ += singlePdes[singlePdes.GetSize()-1] -> GetName();

    
    // Create IterSolveStep instance
    IterSolveStep * solveStep =
        new IterSolveStep( *this, myParam_, infoNode_ ); 
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

  void IterCoupledPDE::GetUpdateGeoForPDE( SolutionType type,
                                           shared_ptr<EntityList>  list,
                                           const std::string& pdeName,
                                           bool& updateGeo ) {
    
    // directly pass the query to the IterSolveStep instance
    return solveStep_->GetUpdateGeoForPDE( type, list, pdeName, updateGeo );
    
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
