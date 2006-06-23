#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "DataInOut/Unverg/outUnverg.hh"
#include "DataInOut/GMV/outGMV.hh"
#include "Forms/forms_header.hh"
#include "Estimator/spaceerror.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "Driver/assemble.hh"
#include "Driver/stdSolveStep.hh"
#include "smoothPDE.hh" 

namespace CoupledField
{

  SmoothPDE::SmoothPDE(Grid * aptgrid, TimeFunc *aptTimeFunc, WriteResults *aptOut )
    :SinglePDE(aptgrid, aptOut, aptTimeFunc)
  {
    ENTER_FCN( "SmoothPDE::SmoothPDE" );
  
    pdename_ = "smooth";
    pdematerialclass_ = MECHANIC;
    firstTurn_ = true;

    // No time step algorithm for this PDE
    isAlwaysStatic_ = true;

  
    // Set flag that the geometry has changed
  
    // Get problem geometry and PDE subtype
    params->Get( "subType", subType_, pdename_ );
    std::string probGeo;
    params->Get( "type", probGeo, "geometry" );

    // =====================================================================
    // set solution information
    // =====================================================================
    shared_ptr<ResultDof> res1(new ResultDof);
    shared_ptr<AnsatzFct> fct(new LagrangeFct);
    res1->resultType = SMOOTH_DISPLACEMENT;
    if( dim_ == 3 ) {
      res1->dofNames = "ux", "uy", "uz";
    } else {
      res1->dofNames = "ux", "uy";
    }
    res1->definedOn = ResultDof::NODE;
    res1->fctType = fct;
    results_.Push_back( res1 );
      
    method_ = "mechanic";

    //is a nonlinear PDE, since in each iteration, we have to setup the matrices new!
    //nonLin_ = true;
  }


  void SmoothPDE::DefineIntegrators()
  {
    ENTER_FCN( "SmoothPDE::DefineIntegerators" );

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


  void SmoothPDE::DefineSolveStep()
  {
    ENTER_FCN( "SmoothPDE::DefineSolveStep" );
    
    solveStep_ = new StdSolveStep(*this); 
  }


  void SmoothPDE::InitCoupling(PDECoupling * coupling)
  {
    ENTER_FCN( "SmoothPDE::Initcoupling" );

    isIterCoupled_ = true;
    ptCoupling_   = coupling; 
    StdVector<UInt> * nodes = NULL;
    UInt dofspernode  = results_[0]->dofNames.GetSize();

    // input couplings
    for (UInt i=0; i<ptCoupling_->GetNumInputCouplings(); i++)
      {
        
        // check for input mechanic displacement
        if (ptCoupling_->GetInputQuantity(i) == MECH_DISPLACEMENT)
          {
            numCouplingBcs_ += (dofspernode * ptCoupling_->GetInputNumNodes(i));
          }
      }

    // output couplings
    for (UInt i=0; i<ptCoupling_->GetNumOutputCouplings(); i++)
      {
        // check for output displacement
        if (ptCoupling_->GetOutputQuantity(i) == SMOOTH_DISPLACEMENT)
          {
            ptCoupling_->CreateCouplingVector(i,isComplex_); 
          }
      }

  }



  void SmoothPDE::CalcOutputCoupling()
  {
    ENTER_FCN( "SmoothPDE::CalcOutputCoupling" );

    SolutionType quantity;
    StdVector<UInt> * couplingnodes;
    CFSVector * values;

    // at first, check if this PDE is iterative coupled
    if (isIterCoupled_ == false)
      return;

    // loop over all output coupling quantities
    for (UInt i=0; i<ptCoupling_->GetNumOutputCouplings(); i++)
      {
        quantity = ptCoupling_->GetOutputQuantity(i);

        switch(ptCoupling_->GetOutputType(i))
          {

          case NODE:
          
            ptCoupling_->GetOutputNodes(i, couplingnodes);
            ptCoupling_->GetOutputValues(i, values);

            if (quantity == SMOOTH_DISPLACEMENT)
              {
                sol_->NodeSolutionToCoupling(*values,*couplingnodes);
              }
          
            break;

          case ELEM:
            Error("No Element coupling output", __FILE__,__LINE__);
          }

      }
  }

  void SmoothPDE::WriteResultsInFile(const UInt kstep,
                                     const Double asteptime,
                                     UInt stepOffset,
                                     Double timeOffset)
  {
    ENTER_FCN( "SmoothPDE::WriteResultsInFile" );
  
    NodeStoreSol<Double> * solTransient;
  
    UInt actStep = kstep + stepOffset;
    Double actTime = asteptime + timeOffset;

    if (analysistype_ == STATIC ||
        analysistype_ == HARMONIC|| analysistype_==MULTIHARMONIC) {
      if (saveSol_) {
        solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);
        outFile_->WriteNodeSolutionTransient(*solTransient, actStep, actTime);
      }
    }
    else
      Error("SmoothPDE: Only static and transient results can be written out",
            __FILE__, __LINE__);
  }


  void SmoothPDE::WriteHistoryInFile(const UInt kstep,
                                     const Double asteptime,
                                     UInt stepOffset,
                                     Double timeOffset)
  {
    ENTER_FCN( "SmoothPDE::WriteHistoryInFile" );
  
  }

  bool SmoothPDE::HasOutput(SolutionType output)
  {
    ENTER_FCN( "SmoothPDE::HasOutput" );
  
    if (output == SMOOTH_DISPLACEMENT)
      return true;

    return false;
  }

  void SmoothPDE::ReadStoreResults()
  {
    ENTER_FCN( "SmoothPDE::ReadStoreRestuls");

    // Construct vectors for restricted parameter search
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    std::string quantity;
    // *****************************
    // Determine nodal results
    // ***************************** 
    StdVector<std::string> nodeValues;
    Enum2String(SMOOTH_DISPLACEMENT, quantity);
    keyVec  = pdename_, "storeResults", "nodeResults", "region";
    attrVec = "", "", "type";
    valVec = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, nodeValues);
    if (nodeValues.GetSize() > 0) {
      saveSol_ = true;
      hasOutput_ = true;
    }

  }


} // end namespace CoupledField
