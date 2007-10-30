// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <math.h>

#include "heatCondPDE.hh"

#include "Forms/laplaceInt.hh"
#include "Forms/linHeatCondInt.hh"
#include "Forms/massInt.hh"
#include "Forms/linNeumannInt.hh"

#include "PDE/trapezoidal.hh"

#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/ansatzFct.hh"

#include "Driver/stdSolveStep.hh"
#include "Driver/assemble.hh"


namespace CoupledField {


  // ======================================================
  // SET SOLUTION INFORMATION
  // ======================================================
  HeatCondPDE::HeatCondPDE(Grid * aptgrid, ParamNode* paramNode )
    :SinglePDE( aptgrid, paramNode ) {

    pdename_          = "heatConduction";
    pdematerialclass_ = THERMIC;
    maxTimeDerivOrder_ = 1;

    nonLin_    = false;
    isElectroCoupled_ = false;
    isMechCoupled_ = false;
    
    InitialCondition_=0.0;
	
  }

  
  
  void HeatCondPDE::SetInitialCondition() {

	try {
		// fetch paramnodes for initial condition
		myParam_->Get("InitialCondition",InitialCondition_);

		//std::cerr << "\n Initial Temperature : " << InitialCondition_ << std::endl;


		if (!isDirectCoupled_ ) {
			if (isComplex_ ) {
				solVec_->Init(Complex(InitialCondition_, 0.0) );
			} else {
				solVec_->Init( InitialCondition_);
			}
			sol_->SetAlgSysDataPointer(solVec_->GetSize(), 
			                    dynamic_cast<Vector<Double>&>(*solVec_).GetPointer() );
			
			// to test -----------------------------------------------------------
			//Vector< Double > h;
			//h.Init();
			//sol_->GetGlobalSolVector(HEAT_TEMPERATURE, h);
			//std::cout << "\n Initial solution vector = " << h.Serialize() << std::endl;
			// to test -----------------------------------------------------------
		}

		
		isSetInitialCondition_ = true;


	}
	catch(Exception & ex ) {
		InitialCondition_=0.0;
		//std::cerr << "\n Initial Temperature : " << InitialCondition_ << std::endl;
	}

  }
  
  
  void HeatCondPDE::ReadSpecialBCs() {

    // First, delete all of the "normal" boundary conditions
    inBcs_.Clear();

    // fetch paramnodes for inbc
    StdVector<ParamNode*> inbcNodes = 
      myParam_->Get("bcsAndLoads")->GetList("neumannInhom");
    
    std::string inName, inDof, inValue, inPhase,  inType;
    std::string inHtc, inTSolid, inTFluid;
    
    // iterate over all parameter nodes
    for( UInt i = 0; i < inbcNodes.GetSize(); i++ ) {
      try {
        inDof = "";
        inbcNodes[i]->Get( "name", inName );
        inbcNodes[i]->Get( "entityType", inType ); 
        inbcNodes[i]->Get( "value", inValue );
        inbcNodes[i]->Get( "dof", inDof );
        inbcNodes[i]->Get( "phase", inPhase );
        inbcNodes[i]->Get( "heatTransferCoefficient", inHtc );
        inbcNodes[i]->Get( "tempSolid", inTSolid );
        inbcNodes[i]->Get( "tempFluid", inTFluid );
        
        // Create inhomogeneous heat Neumann boundary condition
        shared_ptr<InhomHeatNeumannBc> actBc ( new InhomHeatNeumannBc );
        EntityList::ListType listType;
        EntityList::String2Enum( inType, listType );
        shared_ptr<EntityList> actList =
          ptgrid_->GetEntityList( listType,
                                  inName, EntityList::NAMED_NODES );
        actBc->entities = actList;
        actBc->result = results_[0];
        actBc->eqnMap = eqnMap_;
        if( inDof  == "" ) {
          actBc->dof = 1;
        } else {
          actBc->dof = results_[0]->GetDofIndex( inDof );
        }
      actBc->value = inValue;
      actBc->phase = inPhase;
      actBc->htc = inHtc;
      actBc->tSolid = inTSolid;
      actBc->tFluid = inTFluid;
      
      // add definition
      heatInBcs_.Push_back( actBc );

      } catch (Exception & ex ) {
        RETHROW_EXCEPTION( ex, "Can not create inhomogeneous Neumann conditions on '"
                           << inName << "'" );
      }
    }
  }

  void HeatCondPDE::DefineIntegrators() {

    Double density, heatCapacity, thermalConductivity;
    Double coeffmass, coeffstiff;
    Matrix<Double> thermalConductivityTensor;

    //type of geometry-------------------------------------------------------------
	std::string geometryType;
	param->Get("domain")->Get("geometryType", geometryType );

	// convert to tensor type
	SubTensorType tensorType;
	if (geometryType == "3d") {
		tensorType = FULL;
	} else if (geometryType == "plane") {
		tensorType = PLANE_STRAIN;
	} else if (geometryType == "axi") {
		tensorType = AXI;
		isaxi_ = true;
	}
	//-----------------------------------------------------------------------------
    
    for (UInt actSD = 0; actSD < subdoms_.GetSize(); actSD++) {
      
      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( subdoms_[actSD] );

      BaseMaterial * actMat = materials_[subdoms_[actSD]];
      actMat->GetScalar(density,DENSITY,REAL);
      actMat->GetScalar(heatCapacity,HEAT_CAPACITY,REAL);
      //actMat->GetScalar(thermalConductivity,HEAT_CONDUCTIVITY,REAL);

      
      //check how is comming the heat conductivity-----------------------------------
      BaseForm *bilinearStiff = NULL;
      
      if(actMat->IsSet(HEAT_CONDUCTIVITY)){
    	  // stiffness integrator for isotropic material
          actMat->GetScalar(thermalConductivity,HEAT_CONDUCTIVITY,REAL);
          
          coeffstiff = thermalConductivity;
          bilinearStiff = new LaplaceInt(coeffstiff,isaxi_, true );
      }
      else{
          bilinearStiff = new linHeatCondInt( actMat, tensorType, true );
      }
      //-----------------------------------------------------------------------------
      
      // stiffness integrator
//      coeffstiff = thermalConductivity;
//      BaseForm * bilinearStiff = new LaplaceInt(coeffstiff,isaxi_, true );
      BiLinFormContext * stiffContext = 
        new BiLinFormContext(bilinearStiff, STIFFNESS );

      stiffContext->SetPtPdes(this, this);
      
      stiffContext->SetResults( results_[0], results_[0],
                                 actSDList, actSDList );

      // Finally add the standard integrators - stiffness
      assemble_->AddBiLinearForm( stiffContext );
      
	  coeffmass = density * heatCapacity;
		
      if (isElectroCoupled_ == false && isMechCoupled_ == false ) {
			// mass integrator

			BaseForm * bilinearMass = new MassInt(coeffmass, 1, isaxi_, true );
			BiLinFormContext * massContext = new BiLinFormContext(bilinearMass, MASS );

			massContext->SetPtPdes(this, this);
			massContext->SetResults(results_[0], results_[0], actSDList,
					actSDList );

			// Finally add the standard integrators - mass
			assemble_->AddBiLinearForm(massContext );
		} else {
			// damping integrator
			BaseForm * bilinearDamp = new MassInt(coeffmass, 1, isaxi_, true );
			BiLinFormContext * dampContext = new BiLinFormContext(bilinearDamp, DAMPING );

			dampContext->SetPtPdes(this, this);
			dampContext->SetResults(results_[0], results_[0], actSDList,
					actSDList );

			// Finally add the standard integrators - mass
			assemble_->AddBiLinearForm(dampContext );

		}

      
      // Finally add the standard integrators
      //assemble_->AddBiLinearForm( stiffContext );
      //assemble_->AddBiLinearForm( massContext );
      
      // Give result to equation numbering class
      eqnMap_->AddResult( *results_[0], actSDList );
    }

    // Neumann boundary condition
    try{
	    for( UInt iBc = 0; iBc < heatInBcs_.GetSize(); iBc++ ) {
	      
	      // get current Bc
	      InhomHeatNeumannBc const & actBc = *heatInBcs_[iBc];
	
	      // Because the value of the Neumann BC can be determined by either
	      //  'value' or
	      //  'heatTransferCoefficient', 'tSolid' and 'tFluid'
	      // we convert everything to doubles and look at the values
	      std::string value  = actBc.value;
	
	      // If 'htc' is different from zero, we assume the user wants this definition
	      if (actBc.htc  != "0.0") {
	        value = actBc.htc +  " * ( " + actBc.tSolid + "- " + actBc.tFluid + ")";
	        
	        Info->PrintF( pdename_, "For inhomogeneous Neumann BC use \
	        			\n  heat transfer coefficient: %s \n  TempSolid:  %s \n  TempFluid: %s\n\n",
	                    actBc.htc.c_str(), actBc.tSolid.c_str() , actBc.tFluid.c_str() );
	      }
	
	      LinearSurfForm *neumannBC = new LinNeumannInt( value,
	                                                     actBc.phase, 
	                                                     HEAT_CONDUCTIVITY,
	                                                     isaxi_ );
	      neumannBC->SetVoluInfo( materials_ );
	      
	      LinearFormContext * neumannContext =
	        new LinearFormContext( neumannBC );
	      neumannContext->SetPtPde( this );
	      neumannContext->SetResult( actBc.result, actBc.entities );
	      
	      assemble_->AddLinearForm( neumannContext );
	      
	      // Give result to equation numbering class
	      eqnMap_->AddResult( *actBc.result, actBc.entities );
	    }
    }
    catch(Exception & ex){
		RETHROW_EXCEPTION(ex, "Could not calculate Neumann boundary condition"
				<<" in HeatCondPDE" );
    }
    
  }

  void HeatCondPDE::DefineSolveStep() {

    solveStep_ = new StdSolveStep(*this);

  }

  // ======================================================
  // TIME STEPPING SECTION
  // ======================================================

  void HeatCondPDE::InitTimeStepping() {
    
    // Until now no effective mass formulation in the trapezoidal 
    //  integration scheme is implemented!
    TS_alg_ = new Trapezoidal( algsys_ );

  }

  // ======================================================
  // COUPLING SECTION
  // ======================================================

  void HeatCondPDE::InitCoupling(PDECoupling * Coupling) {
    
    isIterCoupled_ = true;
    ptCoupling_   = Coupling;
    
    // Intialize the memory of the coupling values
    for (UInt i=0; i<ptCoupling_->GetNumOutputCouplings(); i++) {
      if (ptCoupling_->GetOutputQuantity(i) == HEAT_TEMPERATURE) {
        ptCoupling_->CreateCouplingVector(i,isComplex_);
        
        // now since we need a incremental formulation, 
        //  initialize some necessary vectors
        isIncrFormulation_ = true;
      }
    }
  }
  

  void HeatCondPDE::CalcOutputCoupling() {

   //  UInt dof;
//     SolutionType quantity;
//     StdVector<Elem*> * couplingElems = NULL;
//     StdVector<UInt> * couplingNodes = NULL;
//     CFSVector * temp_values = NULL;
  
//     // at first, check if this PDE is iterative coupled
//     if (isIterCoupled_ == false)
//       return;

    // loop over all output coupling quantities
//     for (UInt i=0; i<ptCoupling_->GetNumOutputCouplings(); i++) {
      
//       quantity = ptCoupling_->GetOutputQuantity(i);
//       ptCoupling_->GetOutputValues(i, temp_values);

//       // hard coded cast, since coupling is only possible with
//       // real valued entries
//       Vector<Double> * values = dynamic_cast<Vector<Double>*>(temp_values);

//       switch(ptCoupling_->GetOutputType(i)) {
        
//       case NODE:

//         ptCoupling_->GetOutputNodes(i, couplingNodes);
//         ptCoupling_->GetOutputElements(i, couplingElems);
//         dof = ptCoupling_->GetOutputDof(i);    
  
//         break;
//       case ELEM:
//         EXCEPTION("No Element coupling output", __FILE__,__LINE__);
//       }
//     }

  }



  void HeatCondPDE::CalcResults( shared_ptr<BaseResult> result ) {
    
    switch (result->GetResultInfo()->resultType ) {
    case HEAT_TEMPERATURE:
      if( isComplex_ ) {
        ExtractResult<Complex>( result, sol_ );
      } else {
        ExtractResult<Double>( result, sol_ );
      }
      break;
      
    case HEAT_RHS_LOAD:
      if( isComplex_ ) {
        ExtractRhsResult<Complex>( result, results_[0] );
      } else {
        ExtractRhsResult<Double>( result, results_[0] );
      }
      break;
      
    default:   
      Warning( "Resulttype not computable by thermic PDE",
               __FILE__, __LINE__ );
    }
  }



  void HeatCondPDE::DefineAvailResults() {
    
    // === TEMPERATURE ===
    shared_ptr<ResultInfo> res1(new ResultInfo);
    shared_ptr<AnsatzFct> fct(new LagrangeFct);
    res1->resultType = HEAT_TEMPERATURE;
    res1->dofNames = "";
    res1->unit = "K";
    res1->definedOn = ResultInfo::NODE;
    res1->entryType = ResultInfo::SCALAR;
    res1->fctType = fct;
    results_.Push_back( res1 );
    availResults_.insert( res1 );
      
    // === TEMPERATURE RHS ===
    shared_ptr<ResultInfo> rhs(new ResultInfo);
    rhs->resultType = HEAT_RHS_LOAD;
    rhs->dofNames = "";
    rhs->unit = "?";
    rhs->definedOn = ResultInfo::NODE;
    rhs->entryType = ResultInfo::SCALAR;
    rhs->fctType = fct;
    
    availResults_.insert( rhs );

  }
  
  void HeatCondPDE::SetElectroCoupling()
  {
    isElectroCoupled_ = true;

  }
  
  void HeatCondPDE::SetMechCoupling()
  {
    isMechCoupled_ = true;

  }

} // end of namespace CoupledField
