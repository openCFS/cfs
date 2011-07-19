// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
#include <boost/lexical_cast.hpp>

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
#include "Forms/linearForm.hh"
#include "Forms/abcInt.hh"

#include "PDE/trapezoidal.hh"

#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/ansatzFct.hh"
#include "Domain/domain.hh"

#include "Driver/stdSolveStep.hh"
#include "Driver/assemble.hh"
#include "CoupledPDE/pdecoupling.hh"
#include "Optimization/Design/DesignSpace.hh"


namespace CoupledField {


// ======================================================
// SET SOLUTION INFORMATION
// ======================================================
HeatCondPDE::HeatCondPDE(Grid * aptgrid, PtrParamNode paramNode )
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
    myParam_->GetValue("InitialCondition", InitialCondition_);

    //std::cerr << "\n Initial Temperature : " << InitialCondition_ << std::endl;


    if (!isDirectCoupled_ ) {
      if (isComplex_ ) {
        Vector<Complex> & solComplex = dynamic_cast<Vector<Complex>& >(*solVec_); 
        solComplex.Init(Complex(InitialCondition_, 0.0) );
      } else {
        Vector<Double> & solReal = dynamic_cast<Vector<Double>& >(*solVec_);
        solReal.Init( InitialCondition_);
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

  } catch(Exception & ex ) {
    InitialCondition_=0.0;
    //std::cerr << "\n Initial Temperature : " << InitialCondition_ << std::endl;
  }
}


void HeatCondPDE::ReadSpecialBCs() {

  // First, delete all of the Neumann boundary conditions object
  //inBcs_.Clear();

  // fetch paramnodes for Robin boundary condition
  ParamNodeList rbcNodes =
    myParam_->Get("bcsAndLoads")->GetList("robin");

  std::string myDof, myName, myType, myHTC, myBulkTemp;

  // iterate over all parameter nodes
  for( UInt i = 0; i < rbcNodes.GetSize(); i++ ) {
    try {
      myDof.clear();
      rbcNodes[i]->GetValue( "name", myName );
      rbcNodes[i]->GetValue( "entityType", myType );
      rbcNodes[i]->GetValue( "heatTransferCoefficient", myHTC );
      rbcNodes[i]->GetValue( "bulkTemperature", myBulkTemp );

      // Create robin boundary condition
      shared_ptr<RobinBc> actBc ( new RobinBc );

      EntityList::ListType listType = EntityList::listType.Parse(myType);
      shared_ptr<EntityList> actList =
        ptgrid_->GetEntityList( listType, myName, EntityList::REGION );

      actBc->entities = actList;
      actBc->result = results_[0];
      actBc->eqnMap = eqnMap_;
      if( myDof.empty() ) {
        actBc->dof = 1;
      } else {
        actBc->dof = results_[0]->GetDofIndex( myDof );
      }
      actBc->HTC = myHTC;
      actBc->bulkTemp = myBulkTemp;

      // add definition boundary condition list
      robinBcs_.Push_back( actBc );

    } catch (Exception & ex ) {
      RETHROW_EXCEPTION( ex, "Can not create Cauchy boundary condition on '"
          << myName << "'!" );
    }
  }
}

void HeatCondPDE::DefineIntegrators()
{
  Double density, heatCapacity;
  Matrix<Double> thermalConductivityTensor;

  //type of geometry
  std::string geometryType;
  param->Get("domain")->GetValue("geometryType", geometryType );

  // convert to tensor type
  SubTensorType tensorType = FULL;
  if (geometryType == "plane") {
    tensorType = PLANE_STRAIN;
  } else if (geometryType == "axi") {
    tensorType = AXI;
    isaxi_ = true;
  }
  else
    EXCEPTION("subtensortype not implemented");

  // loop over all subdomains
  for (UInt actSD = 0; actSD < subdoms_.GetSize(); actSD++) {

    // create new entity list
    shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
    actSDList->SetRegion( subdoms_[actSD] );

    BaseMaterial * actMat = materials_[subdoms_[actSD]];
    actMat->GetScalar(density,DENSITY,Global::REAL);
    actMat->GetScalar(heatCapacity,HEAT_CAPACITY,Global::REAL);
    //actMat->GetScalar(thermalConductivity,HEAT_CONDUCTIVITY,Global::REAL);

    // ====================================================================
    // stiffness integrator
    // ====================================================================

    BaseForm *bilinearStiff = NULL;
    if( actMat->IsSet(HEAT_CONDUCTIVITY) )
    {
      Double coeffstiff(0.0);
      // stiffness integrator for isotropic material
      actMat->GetScalar(coeffstiff,HEAT_CONDUCTIVITY,Global::REAL);

      bilinearStiff = new LaplaceInt(coeffstiff, isaxi_, true );

      Info->PrintF( pdename_, "Assemble Laplace integrator with multiplicative factor %6.1f.\n", coeffstiff );
    }
    
    if( actMat->IsSet(HEAT_CONDUCTIVITY_TENSOR) )
    {
      bilinearStiff = new linHeatCondInt( actMat, tensorType, true );
    }
    assert(bilinearStiff != NULL);

    BiLinFormContext * stiffContext = new BiLinFormContext(bilinearStiff, STIFFNESS );

    stiffContext->SetPtPdes(this, this);
    stiffContext->SetResults( results_[0], results_[0], actSDList, actSDList );

    // Finally add the standard integrators - stiffness
    assemble_->AddBiLinearForm( stiffContext );


    // ====================================================================
    // mass integrator
    // ====================================================================
    Double coeffmass(density * heatCapacity);

    if (isElectroCoupled_ == false && isMechCoupled_ == false ) {

      BaseForm * bilinearMass = new MassInt(coeffmass, 1, isaxi_, true );
      BiLinFormContext * massContext = new BiLinFormContext(bilinearMass, MASS );

      massContext->SetPtPdes(this, this);
      massContext->SetResults(results_[0], results_[0], actSDList,
          actSDList );

      // Finally add the standard integrators - mass
      assemble_->AddBiLinearForm(massContext );

      Info->PrintF( pdename_,
          "Assemble Mass integrator with multiplicative factor %6.1f.\n", coeffmass );
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

    // Give result to equation numbering class
    eqnMap_->AddResult( *results_[0], actSDList );
  }

  // ======================================================================
  // Neumann boundary condition
  // ======================================================================
  try{
    for( UInt iBc = 0; iBc < inBcs_.GetSize(); iBc++ ) {

      // get current Bc
      InhomNeumannBc const & actBc = *inBcs_[iBc];

      LinearSurfForm *neumannBC =
          new LinNeumannInt( actBc.value, actBc.phase,
                             NO_SOLUTION_TYPE, HEAT_CONDUCTIVITY, isaxi_ );

      neumannBC->SetVoluInfo( materials_ );

      LinearFormContext * neumannContext =
        new LinearFormContext( neumannBC );
      neumannContext->SetPtPde( this );
      neumannContext->SetResult( actBc.result, actBc.entities );

      assemble_->AddLinearForm( neumannContext );

      // Give result to equation numbering class
      eqnMap_->AddResult( *actBc.result, actBc.entities );
    }
  } catch(Exception & ex){
    RETHROW_EXCEPTION(ex, "Could not assemble Neumann boundary condition"
        <<" in HeatCondPDE!" );
  }

  // ======================================================================
  // Robin boundary condition
  // ======================================================================
  try{
    for( UInt iBc = 0; iBc < robinBcs_.GetSize(); iBc++ ) {

      // get current Bc
      RobinBc const & actBc = *robinBcs_[iBc];

      // integrator due to bulk temperature in the surrounding fluid
      //   i.e. convection surface heat flow vector
      Double bTemp = boost::lexical_cast<Double>( actBc.bulkTemp );
      if ( bTemp != 0.0) {

        LinearForm *bulkTempInt = new VolumeSrcInt( actBc.HTC, isaxi_ );

        LinearFormContext * bulkTempContext = new LinearFormContext( bulkTempInt );
        bulkTempContext->SetPtPde( this );
        bulkTempContext->SetResult( actBc.result, actBc.entities );
        assemble_->AddLinearForm( bulkTempContext );

        Info->PrintF( pdename_,
            "Assemble RHS integrator with bulk temperature %6.1f.\n", bTemp );
      } else {
        Info->PrintF( pdename_,
            "Bulk temperature is zero, no linear form on RHS!\n" );
      }

      // integrator due to temperature of the solid, where there is only heat conduction
      //   i.e. convection surface conductivity matrix
      Double hTC = boost::lexical_cast<Double>( actBc.HTC );
      HeatFluxInt *fluxInt = new HeatFluxInt( actBc.HTC, isaxi_ );

      BiLinFormContext * fluxContext = new BiLinFormContext( fluxInt, STIFFNESS );
      fluxContext->SetPtPdes( this, this );
      fluxContext->SetResults( actBc.result, actBc.result,
          actBc.entities, actBc.entities );
      assemble_->AddBiLinearForm( fluxContext );

      // Give result to equation numbering class
      eqnMap_->AddResult( *actBc.result, actBc.entities );

      Info->PrintF( pdename_,
          "Assemble additional integrator for stiffness matrix with HTC %6.1f.\n", hTC );
    }
  } catch(Exception & ex){
    RETHROW_EXCEPTION(ex, "Could not assemble Cauchy boundary condition"
        <<" in HeatCondPDE!" );
  }

  // ======================================================================
  // RHS source values
  // ======================================================================

  // fetch paramnodes for RHS source
  ParamNodeList rhsValuesNodes =
    myParam_->Get("bcsAndLoads")->GetList("rhsValues");

  bool isharmonic;
  std::string rhsFileId = "default";

  try {
    // iterate over all parameter nodes
    for( UInt i = 0; i < rhsValuesNodes.GetSize(); i++ )
    {
      std::string rhsRegion(rhsValuesNodes[i]->Get("region")->As<std::string>());
      rhsValuesNodes[i]->GetValue("isharmonic", isharmonic, ParamNode::EX);
      //rhsFileId = rhsValuesNodes[i]->Get("inputId")->As<std::string>();

      if ( isharmonic ) {
        Info->PrintF( pdename_,
            "Use same RHS source vector for region '%s' for all timesteps.\n\n",
            rhsRegion.c_str() );
      } else {
        Info->PrintF( pdename_,
            "Use RHS source vector for region '%s' from acoustic results.\n\n",
            rhsRegion.c_str() );

      }

      // get material density
      BaseMaterial * actMat = materials_[ ptgrid_->GetRegion().Parse(rhsRegion) ];
      actMat->GetScalar(density,DENSITY,Global::REAL);

      linAcouPowerSourceInt* sourceRHSInt = new linAcouPowerSourceInt( isaxi_,
          isharmonic, rhsFileId, rhsRegion, density );

      LinearFormContext* sourceRHSContext = new LinearFormContext( sourceRHSInt );
      sourceRHSContext->SetPtPde( this );

      shared_ptr<ElemList> rhsElemList( new ElemList(ptgrid_ ) );
      rhsElemList->SetRegion( ptgrid_->GetRegion().Parse(rhsRegion) );
      sourceRHSContext->SetResult( results_[0], rhsElemList );
      assemble_->AddLinearForm( sourceRHSContext );
    }
  }
  catch(Exception & ex)
  {
    RETHROW_EXCEPTION(ex, "Could not assemble RHS source integrator"
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

  PtrParamNode systemNode = FindLinearSystem(pdename_);
  // Until now no effective mass formulation in the trapezoidal
  //  integration scheme is implemented!
  TS_alg_ = new Trapezoidal( algsys_, systemNode );

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


void HeatCondPDE::CalcOutputCoupling()
{
  EXCEPTION("CalcOutputCoupling not implemented!");
  
  //     UInt dof;
  //     SolutionType quantity;
  //     StdVector<Elem*> * couplingElems = NULL;
  //     StdVector<UInt> * couplingNodes = NULL;
  //     SingleVector * temp_values = NULL;
  //     // at first, check if this PDE is iterative coupled
  //     if (isIterCoupled_ == false)
  //       return;
  //     // loop over all output coupling quantities
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
  //         EXCEPTION("No Element coupling output");
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
    
  case OPT_RESULT_1:
  case OPT_RESULT_2:
  case OPT_RESULT_3:
  case OPT_RESULT_4:
  case OPT_RESULT_5:
  case OPT_RESULT_6:
  case OPT_RESULT_7:
  case OPT_RESULT_8:
  case OPT_RESULT_9:
    // design should work, this is checked in AvailabeResults()
    domain->GetErsatzMaterial()->ExtractResults(result, isComplex_);
    break;

  default:
    WARN( "Resulttype not computable by thermic PDE" );
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

void HeatCondPDE::SetElectroCoupling() {

  isElectroCoupled_ = true;
}

void HeatCondPDE::SetMechCoupling() {

  isMechCoupled_ = true;
}

} // end of namespace CoupledField
