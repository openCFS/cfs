// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "Forms/smoothInt.hh"
#include "Forms/smoothNLInt.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Driver/assemble.hh"
#include "Driver/solveStepSmooth.hh"
#include "smoothPDE.hh" 
#include "pseudoTS.hh" 

namespace CoupledField {

  SmoothPDE::SmoothPDE(Grid * aptgrid, PtrParamNode paramNode )
    :SinglePDE(aptgrid, paramNode ) {
    
    pdename_ = "smooth";
    pdematerialclass_ = MECHANIC;
    firstTurn_ = true;
    
    // No time step algorithm for this PDE
    //isAlwaysStatic_ = true;

    // Get problem geometry and PDE subtype
    myParam_->GetValue("subType", subType_ );

    std::string probGeo;
    param->Get("domain")->GetValue("geometryType", probGeo );

    // Set number of degrees of freedom and
    // ensure that subtype fits to problem geometry
    if ( subType_ == "3d" && probGeo == "3d" ) {
      stressDim_ = 6;
      Info->PrintF("", "=== 3D PROBLEM\n");
    }
    else if ( subType_ == "axi" && probGeo == "axi" ) {
      isaxi_ = true;
      stressDim_ = 4;
      Info->PrintF("", "=== AXISYSMMETRIC PROBLEM\n");
    }
    else if ( subType_ == "planeStrain" && probGeo == "plane" ) {
      stressDim_ = 3;
      Info->PrintF("", "=== PLANE STRAIN PROBLEM\n");
    }
      
    else if ( subType_ == "planeStress" && probGeo == "plane" ) {
      stressDim_ = 3;
      Info->PrintF("", "=== PLANE STRESS PROBLEM\n");
    }

    else {
      EXCEPTION( "Subtype '" <<  subType_ << "' of PDE '"
                 <<  pdename_ <<  "' does not fit to problem  geometry '"
                 << probGeo << "'"; );
    }

    method_ = "mechanic";

    //is a nonlinear PDE, since in each iteration, we have to setup the matrices new!
    //nonLin_ = true;
    
    //elastWeight_="byArea";
    elastWeight_          = myParam_->Get("smoothing")->As<std::string>();
    characteristicLength_ = myParam_->Get("size")->As<Double>();
    exponent_             = myParam_->Get("exponent")->As<Double>();

    //    elastWeight_="byStrain";
//    characteristicLength_ = 1.0;
//    exponent_ = -1.3;
    
  }


  void SmoothPDE::DefineIntegrators() {
	  bool coordUpdate = false;
	  
    // Get problem geometry and PDE subtype
    myParam_->GetValue("subType", subType_ );
    
    std::string probGeo;
    param->Get("domain")->GetValue("geometryType", probGeo );

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

      Double A_elem, AMax_elem=0.0, AMin_elem=0.0;
      StdVector<Elem*> elemssd;
      Matrix<Double> ptCoord;

      ptgrid_->GetElems(elemssd,subdoms_[actSD]);

      ptgrid_->GetElemNodesCoord( ptCoord,elemssd[0]->connect);
      A_elem = elemssd[0]->ptElem->CalcVolume(ptCoord, isaxi_);
      AMax_elem=A_elem; AMin_elem=A_elem;

      for(UInt i=1; i< elemssd.GetSize(); i++){
        ptgrid_->GetElemNodesCoord( ptCoord,elemssd[i]->connect);
        A_elem = elemssd[i]->ptElem->CalcVolume(ptCoord, isaxi_);

        if(A_elem>AMax_elem)
          AMax_elem=A_elem;
        if(A_elem<AMin_elem)
          AMin_elem=A_elem;
      }

      // ==============  add stiffness ======================================
      
      BaseMaterial* actSDMat = materials_[subdoms_[actSD]];
      
      BaseForm * bilinearStiff = NULL;
      if (elastWeight_=="byArea"){
        bilinearStiff = new SmoothInt(actSDMat, type, coordUpdate );
      }
      else if (elastWeight_=="byStrain"){
        bilinearStiff = new SmoothNLInt(actSDMat, elastFactors_, couplingNodes_,
                                        AMax_elem, AMin_elem,
                                        characteristicLength_, 
                                        exponent_,
                                        type, coordUpdate);
        bilinearStiff->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
      }

      
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
	  //solveStep_ = new StdSolveStep(*this); 
	  solveStep_ = new SolveStepSmooth(*this); 
  }

  void SmoothPDE::InitNonLin() {
	  if (elastWeight_=="byStrain")
		  nonLin_=true;
	  else
		  nonLin_ = false;
  }

  void SmoothPDE::InitStabParams(){
    Integer paramSize;
    //the size of the matrix with the elastic weighting factors are currently
    //Num2DElements/Num3DElements long.
    //Bec currently i have no idea how to access local elem numbers in the integrator
    paramSize = ptgrid_->GetNumElems();

    for(Integer i=0; i<paramSize; i++) {
      //elastFactors_[i][4]=(-1.0);
      elastFactors_[i][3]=(-1.0);
    }
  }

  void SmoothPDE::InitCoupling(PDECoupling * coupling)
  {

    isIterCoupled_ = true;
    ptCoupling_   = coupling; 
    // StdVector<UInt> * nodes = NULL;
    Vector<Double> globCoord;
    StdVector<UInt> * nodes = NULL;

    UInt dofspernode  = results_[0]->dofNames.GetSize();

    // input couplings
    for (UInt i=0; i<ptCoupling_->GetNumInputCouplings(); i++)
      {
        
        // check for input mechanic displacement
        if (ptCoupling_->GetInputQuantity(i) == MECH_DISPLACEMENT)
          {
            numCouplingBcs_ += (dofspernode * ptCoupling_->GetInputNumNodes(i));

            ptCoupling_->GetInputNodes(i, nodes);
            couplingNodes_.Resize(nodes->GetSize(),2);
            couplingNodes_.Init();
            for ( UInt j = 0; j < nodes->GetSize(); j++ ) {
              ptgrid_->GetNodeCoordinate( globCoord, (*nodes)[j] );

              couplingNodes_[j][0]=globCoord[0];
              couplingNodes_[j][1]=globCoord[1];
              //couplingNodes_[j][2]=globCoord[2];
            }
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
        else if (ptCoupling_->GetOutputQuantity(i) == SMOOTH_VELOCITY)
          {
            ptCoupling_->CreateCouplingVector(i,isComplex_); 
          }
        else if (ptCoupling_->GetOutputQuantity(i) == GRID_VELOCITY)
          {
            ptCoupling_->CreateCouplingVector(i,isComplex_);
          }
      }

  } 

  void SmoothPDE::InitTimeStepping()
  {
    // timestepping formulation
    TS_alg_ = new PseudoTS( algsys_);
  }

  void SmoothPDE::CalcOutputCoupling() {
    
    SolutionType quantity;
    StdVector<UInt> * couplingnodes;
    SingleVector * values;
    
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
        else if (quantity == SMOOTH_VELOCITY) {
        	solDeriv1_.SetAlgSysVector(getDeriv(FIRST_DERIV));     
        	solDeriv1_.NodeSolutionToCoupling((*values),*couplingnodes);
        }
        else if (quantity == GRID_VELOCITY) {
        	solDeriv1_.SetAlgSysVector(getDeriv(FIRST_DERIV));
        	solDeriv1_.NodeSolutionToCoupling((*values),*couplingnodes);
        }

        break;
        
      case ELEM:
        EXCEPTION("No Element coupling output");
      }
      
    }
  }
  
 


  bool SmoothPDE::HasOutput(SolutionType output) {
    
	    if (output == SMOOTH_DISPLACEMENT)
	      return true;
	    else if (output == SMOOTH_VELOCITY)
	      return true;
	    else if (output == GRID_VELOCITY)
	      return true;

    return false;
  }
  
  void SmoothPDE::CalcResults( shared_ptr<BaseResult> res ) {
    
    switch (res->GetResultInfo()->resultType ) {
    case SMOOTH_DISPLACEMENT:
      if( isComplex_ ) {
        ExtractResult<Complex>( res, sol_ );
      } else {
        ExtractResult<Double>( res, sol_ );
      }
      break;
    case SMOOTH_VELOCITY:
      ExtractDerivResult( res, FIRST_DERIV );
      break;
    case GRID_VELOCITY:
        ExtractDerivResult( res, FIRST_DERIV );
      break;
    default:
      WARN( "Resulttype not computable by smoothing PDE" );
    }
  }
  
  
  // ***********************************************************************
  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void SmoothPDE::ReadSpecialResults() {

      solDeriv1_.SetNumSolutions(1);
      solDeriv1_.SetNumNodes(numPDENodes_);
      solDeriv1_.SetSolutionTypeName(SMOOTH_DISPLACEMENT, SMOOTH_VELOCITY);
      solDeriv1_.SetNumDofs(dim_);
      solDeriv1_.SetResults( results_ );
      solDeriv1_.SetPtrEQNData( eqnMap_.get(),ptgrid_); 
      solDeriv1_.SetRegions( subdoms_ );
      solDeriv1_.Init();
  }

  
  void SmoothPDE::DefineAvailResults() {
	  //the size of the matrix with the stabilization parameters is currently
      //Num2DElements/Num3DElements long. 
      //Bec currently i have no idea how to access local elem numbers in the integrator
	  elastFactors_.Resize(ptgrid_->GetNumElems()+1,4);
      elastFactors_.InitValue(-1.0);    	    	

	  
    // =====================================================================
    // set solution information
    // =====================================================================
	  StdVector<std::string> dispDofNames;
	  // === SMOOTH DISPLACEMENT ===
    shared_ptr<ResultInfo> res1(new ResultInfo);
    shared_ptr<AnsatzFct> fct(new LagrangeFct);
    res1->resultType = SMOOTH_DISPLACEMENT;
    if( dim_ == 3 ) {
    	dispDofNames = "x", "y", "z";
    } else {
      if( isaxi_ ) {
    	  dispDofNames = "r", "z";
      } else {
    	  dispDofNames = "x", "y";
      }
    }
    res1->dofNames = dispDofNames;
    res1->unit = "m";
    res1->definedOn = ResultInfo::NODE;
    res1->fctType = fct;
    res1->entryType = ResultInfo::VECTOR;
    results_.Push_back( res1 );
    availResults_.insert( res1 );

    // === SMOOTH VELOCITY ===
    shared_ptr<ResultInfo> vel(new ResultInfo);
    vel->resultType = SMOOTH_VELOCITY;
    vel->dofNames = dispDofNames;
    vel->unit = "m/s";
    vel->entryType = res1->entryType;
    vel->definedOn = res1->definedOn;
    vel->fctType = res1->fctType;
    availResults_.insert( vel );
  }

 

} // end namespace CoupledField
