// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
#include <iostream>
#include <string>
#include <utility>

#include "CoupledPDE/pdecoupling.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/WriteInfo.hh"
#include "Domain/ansatzFct.hh"
#include "Domain/domain.hh"
#include "Domain/entityList.hh"
#include "Domain/grid.hh"
#include "Domain/resultInfo.hh"
#include "Driver/assemble.hh"
#include "Driver/formsContext.hh"
#include "Driver/stdSolveStep.hh"
#include "Forms/abcInt.hh"
#include "Forms/baseForm.hh"
#include "Forms/laplaceInt.hh"
#include "Forms/linHeatCondInt.hh"
#include "Forms/linNeumannInt.hh"
#include "Forms/linSurfForm.hh"
#include "Forms/linearForm.hh"
#include "Forms/massInt.hh"
#include "Forms/nLinHeatInt.hh"
#include "Forms/nonConformingInt.hh"
#include "General/Enum.hh"
#include "General/defs.hh"
#include "General/exception.hh"
#include "MatVec/SingleVector.hh"
#include "MatVec/matrix.hh"
#include "MatVec/vector.hh"
#include "Materials/baseMaterial.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/StdPDE.hh"
#include "PDE/eqnMap.hh"
#include "PDE/trapezoidal.hh"
#include "Utils/basenodestoresol.hh"
#include "Utils/nodestoresol.hh"
#include "Utils/result.hh"
#include "boost/lexical_cast.hpp"
#include "boost/tokenizer.hpp"
#include "heatCondPDE.hh"
#include "math.h"


namespace CoupledField {

  DECLARE_LOG(heatcondpde)
  DEFINE_LOG(heatcondpde, "pde.heatcond")

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
    // fetch paramnodes for initial condition (if not present, leave)
    if ( !myParam_->Has("InitialCondition") ) {
      return;
    }
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

    // read volume force definition
    ReadRegionLoads();

}


  // ****************************
  //  Initialize Nonlinearities
  // ****************************
  void HeatCondPDE::InitNonLin() {

    SinglePDE::InitNonLin();

    //now do PDE specifics
    if ( nonLin_ ) {
      //we perform a total formulation (no incremental one)
      totalFormulation_ = true;
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
//  else
//    EXCEPTION("subtensortype not implemented");

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

    //get possible nonlinearities defined in this region
    //regionID = subdoms_[actSD]
    StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[subdoms_[actSD]]; 
    
    if ( nonLinTypes.Find(NLHEAT_CONDUCTIVITY) != -1 ) {
      // informs material that approx./interpol. for heat conductivity is needed
      std::cout << "Do NL Cond" << std::endl;

      actMat->NeedApproxMatCurve( HEAT_CONDUCTIVITY );

      BaseForm *nlBilinearStiff = new nlinHeatStiffInt( actMat, tensorType, false );

      nlBilinearStiff->SetNonLinMethod( nonLinMethod_ );      
      nlBilinearStiff->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));

      BiLinFormContext *stiffContext = new BiLinFormContext(nlBilinearStiff, STIFFNESS );
      
      stiffContext->SetPtPdes(this, this);
      stiffContext->SetResults( results_[0], results_[0], actSDList, actSDList );
      
      // Finally add the standard integrators - stiffness
      assemble_->AddBiLinearForm( stiffContext );

//       // nonlinear RHS linearform!!
//       LinearForm * rhsNL;
//       rhsNL = new nLinHeat_linFormInt( actMat, isaxi_, false );
//       rhsNL->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));

//       LinearFormContext * rhsContext = new LinearFormContext( rhsNL );
//       rhsContext->SetPtPde( this );
//       rhsContext->SetResult( results_[0], actSDList );
//       assemble_->AddLinearForm( rhsContext );
    }
    else {
      if( actMat->IsSet(HEAT_CONDUCTIVITY) ) {

        BaseForm *bilinearStiff; 
        Double coeffstiff(0.0);
        // stiffness integrator for isotropic material
        actMat->GetScalar(coeffstiff,HEAT_CONDUCTIVITY,Global::REAL);
        bilinearStiff = new LaplaceInt(coeffstiff, isaxi_, true );
        
        Info->PrintF( pdename_, "Assemble Laplace integrator with multiplicative factor %6.1f.\n", coeffstiff );
        BiLinFormContext * stiffContext = new BiLinFormContext(bilinearStiff, STIFFNESS );
        
        stiffContext->SetPtPdes(this, this);
        stiffContext->SetResults( results_[0], results_[0], actSDList, actSDList );
      
        // Finally add the standard integrators - stiffness
        assemble_->AddBiLinearForm( stiffContext );
      }
      else if( actMat->IsSet(HEAT_CONDUCTIVITY_TENSOR) ) {
        BaseForm *bilinearStiff; 
        bilinearStiff = new linHeatCondInt( actMat, tensorType, true );
        BiLinFormContext * stiffContext = new BiLinFormContext(bilinearStiff, STIFFNESS );
        
        stiffContext->SetPtPdes(this, this);
        stiffContext->SetResults( results_[0], results_[0], actSDList, actSDList );
      
        // Finally add the standard integrators - stiffness
        assemble_->AddBiLinearForm( stiffContext );
      }

//       if ( nonLin_ == true ) {
//           // we need nonlinear RHS linearform for both linear and nonlinear
//           // subdomains; just in case of material nonlinearity!
//         LinearForm * rhsNL;
//         rhsNL = new nLinHeat_linFormInt( actMat, isaxi_, false );
//         rhsNL->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
        
//         LinearFormContext * rhsContext = new LinearFormContext( rhsNL );
//         rhsContext->SetPtPde( this );
//         rhsContext->SetResult( results_[0], actSDList );
//         assemble_->AddLinearForm( rhsContext );
//      }
      
    }


    // ====================================================================
    // mass integrator
    // ====================================================================
    Double coeffmass(density * heatCapacity);

    if (isElectroCoupled_ == false && isMechCoupled_ == false ) {

      if ( nonLinTypes.Find( NLHEAT_CAPACITY) != -1 ) {
        // informs material that approx./interpol. for heat conductivity is needed
        actMat->NeedApproxMatCurve( HEAT_CAPACITY );

        BaseForm *nlBilinearMass = new nlinHeatMassInt( actMat, tensorType, false );

        nlBilinearMass->SetNonLinMethod( nonLinMethod_ );      
        nlBilinearMass->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
        
        BiLinFormContext *massContext = new BiLinFormContext(nlBilinearMass, MASS );
        
        massContext->SetPtPdes(this, this);
        massContext->SetResults( results_[0], results_[0], actSDList, actSDList );
        
        // Finally add the standard integrators - stiffness
        assemble_->AddBiLinearForm( massContext );
      }
      else {
        BaseForm * bilinearMass = new MassInt(coeffmass, 1, isaxi_, true );
        BiLinFormContext * massContext = new BiLinFormContext(bilinearMass, MASS );
        
        massContext->SetPtPdes(this, this);
        massContext->SetResults(results_[0], results_[0], actSDList,
                                actSDList );
        
        // Finally add the standard integrators - mass
        assemble_->AddBiLinearForm(massContext );
        
        Info->PrintF( pdename_,
                      "Assemble Mass integrator with multiplicative factor %6.1f.\n", coeffmass );
      }
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

    // =======================================================================
    // Integrators for NonConforming Interfaces
    // =======================================================================
    
    // Get index of LAGRANGE_MULT result, just in case... who knows...
    UInt lmResultIdx = 0;
    for(UInt i=0, n=results_.GetSize(); i<n; i++) {
      if(results_[i]->resultType == LAGRANGE_MULT) {
        lmResultIdx = i;
        break;
      }
    }
    LOG_DBG2(heatcondpde) << "NonMatching: Index of LAGRANGE_MULT result: "
                     << lmResultIdx;
    
    for( UInt i = 0; i < ncIFaces_.GetSize(); i++ ) {
      
      // get regionId of Lagrangian surface
      StdVector<std::string> keyVec, attrVec, valVec;
      std::string slaveSide;
      std::string ncIfaceName = ptgrid_->GetRegion().ToString(ncIFaces_[i]);

      PtrParamNode ncIfaceListNode;
      ncIfaceListNode = param->Get("domain")->Get("ncInterfaceList");

      slaveSide = ncIfaceListNode->
          GetByVal("ncInterface", "name",
                   ncIfaceName)->Get("slaveSide")->As<std::string>();

      // Part 1: Define integrator M(Psi, Lambda) on
      //         non-conforming interface
      LOG_DBG2(heatcondpde) << "NonMatching: Defining nonconforming integrator"
                        << " for M on interface '"
                        << ptgrid_->GetRegion().ToString(ncIFaces_[i]) << "'.";
      shared_ptr<ElemList> actNcList( new ElemList(ptgrid_ ) );
      actNcList->SetRegion( ncIFaces_[i] );
      
      NonConformingInt * ncInt = 
        new NonConformingInt( 1, isaxi_ );

      NcBiLinFormContext * stiffIntDescr = 
     	  new NcBiLinFormContext( ncInt , STIFFNESS );

      // Force assembling of M(Psi, Lambda)^T
      stiffIntDescr->SetCounterPart( true );

      stiffIntDescr->SetPtPdes(this, this);
      stiffIntDescr->SetResults( results_[0], results_[lmResultIdx],
                                 actNcList, actNcList );
      
      assemble_->AddBiLinearForm( stiffIntDescr );


      // Part 2: Define integrator D(Psi, Lambda) on
      //         Lagrangian surface
      LOG_DBG2(heatcondpde) << "NonMatching: Defining mass integrator"
                        << " for D on interface '"
                        << ptgrid_->GetRegion().ToString(ncIFaces_[i]) << "'.";
      shared_ptr<SurfElemList> actSDList( new SurfElemList(ptgrid_ ) );      
      actSDList->SetRegion( ptgrid_->GetRegion().Parse(slaveSide));

      // D(Psi, Lambda) has the form of a standard mass
      // integrator with factor 1.0
      MassInt * dMatInt = new MassInt( 1.0, 1, isaxi_ );
      BiLinFormContext * dMatContext = 
        new BiLinFormContext( dMatInt, STIFFNESS );

      // Force assembling of D(Psi, Lambda)^T
      dMatContext->SetCounterPart( true );
      dMatContext->SetPtPdes( this, this );
      dMatContext->SetResults( results_[0], results_[lmResultIdx],
                               actSDList, actSDList );
      
      assemble_->AddBiLinearForm( dMatContext );

      // Give result LAGRANGE_MULT to equation numbering class
      eqnMap_->AddResult( *results_[lmResultIdx], actSDList );
    }


  // ======================================================================
  // RHS source values
  // ======================================================================

  // Add integrators for region loads
  DefineRegionLoadIntegrators(regionLoads_);

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


void HeatCondPDE::DefineRegionLoadIntegrators(std::map<RegionIdType, RegionLoad>& regionLoads, StdVector<LinearFormContext*>* linForms){
    VolumeSrcInt * volSrcInt;

    std::map<RegionIdType, RegionLoad>::iterator loadIt = regionLoads.begin();
    for( loadIt = regionLoads.begin(); loadIt != regionLoads.end(); loadIt++ ) {
      volSrcInt = (*loadIt).second.GetSrcScalarIntegrator();

      // Create new element list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( loadIt->first );
      LinearFormContext * volSrcContext =
        new LinearFormContext( volSrcInt );
      volSrcContext->SetPtPde(this);
      volSrcContext->SetResult( results_[0], actSDList );
      if(linForms != NULL){
        linForms->Push_back(volSrcContext);
      }else{
        assemble_->AddLinearForm( volSrcContext );
      }

      (*loadIt).second.ToInfo(infoNode_->Get("regionLoad"));

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


    // ===================================
    // Check for non-conforming interfaces
    // ===================================
    StdVector<std::string> ncIfaceNames, ncIfaceNamesForPDE;
    StdVector<RegionIdType> ncIfaceIds;
    
    LOG_DBG2(heatcondpde) << "NonMatching: Checking if nonconforming "
                      << "interfaces of PDE exist in domain.";

    PtrParamNode heatcondpdeNCIfaceListNode;
    heatcondpdeNCIfaceListNode = param->GetByVal("sequenceStep", std::string("index"), sequenceStep_)
    ->Get("pdeList/heatConduction/ncInterfaceList", ParamNode::PASS);
    
    if(!heatcondpdeNCIfaceListNode)
      return;

    PtrParamNode domainNCIfaceListNode;
    domainNCIfaceListNode = param->Get("domain")->Get("ncInterfaceList", ParamNode::PASS);

    if(!domainNCIfaceListNode)
    {
      EXCEPTION("No nonmatching interfaces have been specified in domain!");
    }

    ParamNodeList pdeNCIfaceNodes;
    pdeNCIfaceNodes = heatcondpdeNCIfaceListNode->GetList("ncInterface");

    for (UInt i = 0; i < pdeNCIfaceNodes.GetSize(); i++) {
      std::string pdeIfaceName = pdeNCIfaceNodes[i]->Get("name")->As<std::string>();

      PtrParamNode domainIfaceNode = domainNCIfaceListNode
          ->GetByVal("ncInterface", "name", pdeIfaceName, ParamNode::PASS);
      if(!domainIfaceNode)
      {
        LOG_DBG2(heatcondpde) << "NonMatching: Nonconforming "
        << "interface '" << ncIfaceNames[i]
                                         << "' does not exist in domain.";

        EXCEPTION( "ncInterface referenced from PDE not defined in domain!");
      }

      ncIfaceNamesForPDE.Push_back(pdeIfaceName);
    }
    ptgrid_->GetRegion().Parse(ncIfaceNamesForPDE, ncIfaceIds);

    for (UInt i = 0; i < ncIfaceIds.GetSize(); i++) {
      ncIFaces_.Push_back(ncIfaceIds[i]);
    }

    // In the case of the presence of non-conforming interfaces,
    // a second resultdof object has to be created, which describes the 
    // Lagrange multiplier
    if( ncIFaces_.GetSize() > 0 ) {
      LOG_DBG2(heatcondpde) << "NonMatching: Defining new ResultDof Lagrange.";
      shared_ptr<ResultInfo> lagr ( new ResultInfo );
      lagr->resultType = LAGRANGE_MULT;
      lagr->dofNames = "l";
      lagr->fctType = results_[0]->fctType;
      lagr->definedOn = results_[0]->definedOn;
      results_.Push_back( lagr );
    } 

}

void HeatCondPDE::SetElectroCoupling() {

  isElectroCoupled_ = true;
}

void HeatCondPDE::SetMechCoupling() {

  isMechCoupled_ = true;
}

} // end of namespace CoupledField
