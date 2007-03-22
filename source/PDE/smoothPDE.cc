// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "Forms/forms_header.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Driver/assemble.hh"
#include "Driver/stdSolveStep.hh"
#include "smoothPDE.hh" 

namespace CoupledField {

  SmoothPDE::SmoothPDE(Grid * aptgrid, ParamNode* paramNode )
    :SinglePDE(aptgrid, paramNode ) {
    ENTER_FCN( "SmoothPDE::SmoothPDE" );
    
    pdename_ = "smooth";
    pdematerialclass_ = MECHANIC;
    firstTurn_ = true;
    
    // No time step algorithm for this PDE
    isAlwaysStatic_ = true;

  
  

    
    method_ = "mechanic";

    //is a nonlinear PDE, since in each iteration, we have to setup the matrices new!
    //nonLin_ = true;
  }


  void SmoothPDE::DefineIntegrators() {
    ENTER_FCN( "SmoothPDE::DefineIntegerators" );
    
    // Get problem geometry and PDE subtype
    myParam_->Get("subType", subType_ );
    
    std::string probGeo;
    param->Get("domain")->Get("geometryType", probGeo );

    // Weigthening factors for smoothening 
    factor_.Resize(numElems_);
    for (UInt i=0; i<factor_.GetSize(); i++)
      factor_[i] = 1.0;

    //transform the type
    SubTensorType type;
    String2Enum(subType_,type);

    for (UInt actSD = 0; actSD < subdoms_.GetSize(); actSD++) {
      
      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( subdoms_[actSD] );
      
      // ==============  add stiffness ======================================
      
      BaseMaterial* actSDMat = materials_[subdoms_[actSD]];
      
      BaseForm * bilinearStiff;
      bilinearStiff = new SmoothInt(actSDMat, type, true );
      
      BiLinFormContext * stiffContext = 
        new BiLinFormContext( bilinearStiff, STIFFNESS );
      stiffContext->SetPtPdes(this, this);
      stiffContext->SetResults( results_[0], results_[0],
                                actSDList, actSDList );  
      assemble_->AddBiLinearForm( stiffContext );
      
      // Give result to equation numbering class
      eqnMap_->AddResult( *results_[0], actSDList );
    }
  }


  void SmoothPDE::DefineSolveStep() {
    ENTER_FCN( "SmoothPDE::DefineSolveStep" );
    
    solveStep_ = new StdSolveStep(*this); 
  }


  void SmoothPDE::InitCoupling(PDECoupling * coupling) {
    ENTER_FCN( "SmoothPDE::Initcoupling" );

    isIterCoupled_ = true;
    ptCoupling_   = coupling; 
    StdVector<UInt> * nodes = NULL;
    UInt dofspernode  = results_[0]->dofNames.GetSize();

    // input couplings
    for (UInt i=0; i<ptCoupling_->GetNumInputCouplings(); i++) {
        
      // check for input mechanic displacement
      if (ptCoupling_->GetInputQuantity(i) == MECH_DISPLACEMENT) {
        numCouplingBcs_ += (dofspernode * ptCoupling_->GetInputNumNodes(i));
      }
    }

    // output couplings
    for (UInt i=0; i<ptCoupling_->GetNumOutputCouplings(); i++) {
      // check for output displacement
      if (ptCoupling_->GetOutputQuantity(i) == SMOOTH_DISPLACEMENT) {
        ptCoupling_->CreateCouplingVector(i,isComplex_); 
      }
    }
    
  }



  void SmoothPDE::CalcOutputCoupling() {
    ENTER_FCN( "SmoothPDE::CalcOutputCoupling" );
    
    SolutionType quantity;
    StdVector<UInt> * couplingnodes;
    CFSVector * values;
    
    // at first, check if this PDE is iterative coupled
    if (isIterCoupled_ == false)
      return;
    
    // loop over all output coupling quantities
    for (UInt i=0; i<ptCoupling_->GetNumOutputCouplings(); i++) {
      quantity = ptCoupling_->GetOutputQuantity(i);

      switch(ptCoupling_->GetOutputType(i)) {
        
      case NODE:
        
        ptCoupling_->GetOutputNodes(i, couplingnodes);
        ptCoupling_->GetOutputValues(i, values);
        
        if (quantity == SMOOTH_DISPLACEMENT) {
          sol_->NodeSolutionToCoupling(*values,*couplingnodes);
        }
        
        break;
        
      case ELEM:
        EXCEPTION("No Element coupling output");
      }
      
    }
  }
  
 


  bool SmoothPDE::HasOutput(SolutionType output) {
    ENTER_FCN( "SmoothPDE::HasOutput" );
    
    if (output == SMOOTH_DISPLACEMENT)
      return true;
    
    return false;
  }
  
  void SmoothPDE::CalcResults( shared_ptr<BaseResult> res ) {
    ENTER_FCN( "SmoothPDE::CalcResults" );
    
    switch (res->GetResultInfo()->resultType ) {
    case SMOOTH_DISPLACEMENT:
      if( isComplex_ ) {
        ExtractResult<Complex>( res, sol_ );
      } else {
        ExtractResult<Double>( res, sol_ );
      }
      break;
    default:
      Warning( "Resulttype not computable by smoothing PDE",
               __FILE__, __LINE__ );
    }
  }
    
  void SmoothPDE::DefineAvailResults() {
    ENTER_FCN( "SmoothPDE::DefineAvailResults" );

    // =====================================================================
    // set solution information
    // =====================================================================
    shared_ptr<ResultInfo> res1(new ResultInfo);
    shared_ptr<AnsatzFct> fct(new LagrangeFct);
    res1->resultType = SMOOTH_DISPLACEMENT;
    if( dim_ == 3 ) {
      res1->dofNames = "x", "y", "z";
    } else {
      if( isaxi_ ) {
        res1->dofNames = "r", "z";
      } else {
        res1->dofNames = "x", "y";
      }
    }
    res1->unit = "m";
    res1->definedOn = ResultInfo::NODE;
    res1->fctType = fct;
    res1->entryType = ResultInfo::VECTOR;
    results_.Push_back( res1 );
    availResults_.insert( res1 );
      
  }

 

} // end namespace CoupledField
