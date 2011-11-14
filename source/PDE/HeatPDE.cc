// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
#include <boost/lexical_cast.hpp>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <math.h>

#include <boost/tokenizer.hpp>

#include "HeatPDE.hh"

#include "General/defs.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Domain/CoefFunction/CoefFunction.hh"
#include "Utils/StdVector.hh"
#include "Driver/SolveSteps/SolveStepElec.hh"
#include "CoupledPDE/PDECoupling.hh"
#include "Driver/Assemble.hh"
#include "Utils/ApproxData.hh"
#include "Utils/SmoothSpline.hh"
#include "Materials/Models/Hysteresis.hh"
#include "Materials/Models/Preisach.hh"
#include "FeBasis/H1/FeSpaceH1.hh"
#include "FeBasis/H1/FeSpaceH1Nodal.hh"
#include "FeBasis/FeFunctions.hh"
#include "DataInOut/WriteInfo.hh"

#include "Driver/TimeSchemes/Trapezoidal.hh"

#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/Assemble.hh"
#include "CoupledPDE/PDECoupling.hh"

//new integrator concept
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"

//new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"

namespace CoupledField {

  DECLARE_LOG(heatcondpde)
  DEFINE_LOG(heatcondpde, "pde.heatcond")

// ======================================================
// SET SOLUTION INFORMATION
// ======================================================
HeatPDE::HeatPDE(Grid * aptgrid, PtrParamNode paramNode )
:SinglePDE( aptgrid, paramNode ) {

  pdename_           = "heatConduction";
  pdematerialclass_  = THERMIC;
  maxTimeDerivOrder_ = 1;
  nonLin_            = false;
  InitialCondition_  = 0.0;

}



void HeatPDE::SetInitialCondition() {

  REFACTOR;

//   try {
//     // fetch paramnodes for initial condition (if not present, leave)
//     if ( !myParam_->Has("InitialCondition") ) {
//       return;
//     }
//     myParam_->GetValue("InitialCondition", InitialCondition_);

//     //std::cerr << "\n Initial Temperature : " << InitialCondition_ << std::endl;


//     if (!isDirectCoupled_ ) {
//       if (isComplex_ ) {
//         Vector<Complex> & solComplex = dynamic_cast<Vector<Complex>& >(*solVec_); 
//         solComplex.Init(Complex(InitialCondition_, 0.0) );
//       } else {
//         Vector<Double> & solReal = dynamic_cast<Vector<Double>& >(*solVec_);
//         solReal.Init( InitialCondition_);
//       }
//       sol_->SetAlgSysDataPointer(solVec_->GetSize(),
//           dynamic_cast<Vector<Double>&>(*solVec_).GetPointer() );

//       // to test -----------------------------------------------------------
//       //Vector< Double > h;
//       //h.Init();
//       //sol_->GetGlobalSolVector(HEAT_TEMPERATURE, h);
//       //std::cout << "\n Initial solution vector = " << h.Serialize() << std::endl;
//       // to test -----------------------------------------------------------
//     }

//     isSetInitialCondition_ = true;

//   } catch(Exception & ex ) {
//     InitialCondition_=0.0;
//     //std::cerr << "\n Initial Temperature : " << InitialCondition_ << std::endl;
//   }
}


void HeatPDE::ReadSpecialBCs() {

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
      //      actBc->eqnMap = eqnMap_;
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

  std::map<SolutionType, shared_ptr<FeSpace> >
  HeatPDE::CreateFeSpaces(const std::string& formulation, PtrParamNode infoNode) {
    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    if(formulation == "default" || formulation == "H1"){
      PtrParamNode potSpaceNode = infoNode->Get("heatTemperature");
      crSpaces[HEAT_TEMPERATURE] =
        FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1);
      crSpaces[HEAT_TEMPERATURE]->Init();
    }else{
      EXCEPTION("The formulation " << formulation << "of heat PDE is not known!");
    }
    return crSpaces;
  }


  // ****************************
  //  Initialize Nonlinearities
  // ****************************
  void HeatPDE::InitNonLin() {

    SinglePDE::InitNonLin();

    //now do PDE specifics
    if ( nonLin_ ) {
      //we perform a total formulation (no incremental one)
      totalFormulation_ = true;
    }

  }


void HeatPDE::DefineIntegrators() {

  RegionIdType actRegion;
  BaseMaterial * actSDMat = NULL;  

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

  // Define integrators for "standard" materials
  std::map<RegionIdType, BaseMaterial*>::iterator it;
  shared_ptr<FeSpace> mySpace = feFunctions_[HEAT_TEMPERATURE]->GetFeSpace();
  
  for ( it = materials_.begin(); it != materials_.end(); it++ ) {
    
    // Set current region and material
    actRegion = it->first;
    actSDMat = it->second;
    
    // Get current region name
    std::string regionName = ptgrid_->GetRegion().ToString(actRegion);
    
    // create new entity list
    shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
    actSDList->SetRegion( actRegion );
    
    // --- Set the FE ansatz for the current region ---
    PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
    std::string polyId = curRegNode->Get("polyId")->As<std::string>();
    std::string integId = curRegNode->Get("integId")->As<std::string>();
    mySpace->SetRegionApproximation(actRegion, polyId,integId);
    
    // ====================================================================
    // stiffness integrator
    // ====================================================================

    // --- standard real-valued stiffness integrator ---
    shared_ptr<CoefFunction > curCoef = 
      actSDMat->GetCoefFunction(HEAT_CONDUCTIVITY, tensorType, Global::REAL);
    
    
    BDBInt< GradientOperator ,FeH1,Double,Double >* stiffInt;
    stiffInt = new BDBInt<GradientOperator,FeH1,Double,Double >(curCoef,1.0 );

    BiLinFormContext * stiffIntDescr =
      new BiLinFormContext(stiffInt, STIFFNESS );
    
    feFunctions_[HEAT_TEMPERATURE]->AddEntityList( actSDList );

    //stiffIntDescr->SetPtPdes(this, this);
    stiffIntDescr->SetEntities( actSDList, actSDList );
    stiffIntDescr->SetFeFunctions(feFunctions_[HEAT_TEMPERATURE],feFunctions_[HEAT_TEMPERATURE]);
    stiffInt->SetFeSpace( feFunctions_[HEAT_TEMPERATURE]->GetFeSpace());
    
    assemble_->AddBiLinearForm( stiffIntDescr );
    bdbInts_[actRegion] = stiffInt;
      
    // ====================================================================
    // mass integrator
    // ====================================================================
    Double density, heatCapacity, massFactor;
    actSDMat->GetScalar(density,DENSITY,Global::REAL);
    actSDMat->GetScalar(heatCapacity,HEAT_CAPACITY,Global::REAL);
    massFactor = density * heatCapacity;

    BBInt<IdentityOperator,FeH1,Double> *massInt;
    massInt = new BBInt<IdentityOperator, FeH1, Double>(massFactor);
    massInt->SetFeSpace( feFunctions_[HEAT_TEMPERATURE]->GetFeSpace() );

    BiLinFormContext *massContext =  new BiLinFormContext(massInt, MASS );
     
    massContext->SetEntities( actSDList, actSDList );
    massContext->SetFeFunctions( feFunctions_[HEAT_TEMPERATURE],feFunctions_[HEAT_TEMPERATURE]);
    assemble_->AddBiLinearForm( massContext );
  }

  // ======================================================================
  // Neumann boundary condition
  // ======================================================================
//   try{
//     for( UInt iBc = 0; iBc < inBcs_.GetSize(); iBc++ ) {

//       // get current Bc
//       InhomNeumannBc const & actBc = *inBcs_[iBc];

//       LinearSurfForm *neumannBC =
//           new LinNeumannInt( actBc.value, actBc.phase,
//                              NO_SOLUTION_TYPE, HEAT_CONDUCTIVITY, isaxi_ );

//       neumannBC->SetVoluInfo( materials_ );

//       LinearFormContext * neumannContext =
//         new LinearFormContext( neumannBC );
//       neumannContext->SetPtPde( this );
//       neumannContext->SetResult( actBc.result, actBc.entities );

//       assemble_->AddLinearForm( neumannContext );

//       // Give result to equation numbering class
//       eqnMap_->AddResult( *actBc.result, actBc.entities );
//     }
//   } catch(Exception & ex){
//     RETHROW_EXCEPTION(ex, "Could not assemble Neumann boundary condition"
//         <<" in HeatPDE!" );
//   }

//   // ======================================================================
//   // Robin boundary condition
//   // ======================================================================
//   try{
//     for( UInt iBc = 0; iBc < robinBcs_.GetSize(); iBc++ ) {

//       // get current Bc
//       RobinBc const & actBc = *robinBcs_[iBc];

//       // integrator due to bulk temperature in the surrounding fluid
//       //   i.e. convection surface heat flow vector
//       Double bTemp = boost::lexical_cast<Double>( actBc.bulkTemp );
//       if ( bTemp != 0.0) {

//         LinearForm *bulkTempInt = new VolumeSrcInt( actBc.HTC, isaxi_ );

//         LinearFormContext * bulkTempContext = new LinearFormContext( bulkTempInt );
//         bulkTempContext->SetPtPde( this );
//         bulkTempContext->SetResult( actBc.result, actBc.entities );
//         assemble_->AddLinearForm( bulkTempContext );

//         Info->PrintF( pdename_,
//             "Assemble RHS integrator with bulk temperature %6.1f.\n", bTemp );
//       } else {
//         Info->PrintF( pdename_,
//             "Bulk temperature is zero, no linear form on RHS!\n" );
//       }

//       // integrator due to temperature of the solid, where there is only heat conduction
//       //   i.e. convection surface conductivity matrix
//       Double hTC = boost::lexical_cast<Double>( actBc.HTC );
//       HeatFluxInt *fluxInt = new HeatFluxInt( actBc.HTC, isaxi_ );

//       BiLinFormContext * fluxContext = new BiLinFormContext( fluxInt, STIFFNESS );
//       fluxContext->SetPtPdes( this, this );
//       fluxContext->SetResults( actBc.result, actBc.result,
//           actBc.entities, actBc.entities );
//       assemble_->AddBiLinearForm( fluxContext );

//       // Give result to equation numbering class
//       eqnMap_->AddResult( *actBc.result, actBc.entities );

//       Info->PrintF( pdename_,
//           "Assemble additional integrator for stiffness matrix with HTC %6.1f.\n", hTC );
//     }
//   } catch(Exception & ex){
//     RETHROW_EXCEPTION(ex, "Could not assemble Cauchy boundary condition"
//         <<" in HeatPDE!" );
//   }

//     // =======================================================================
//     // Integrators for NonConforming Interfaces
//     // =======================================================================
    
//     // Get index of LAGRANGE_MULT result, just in case... who knows...
//     UInt lmResultIdx = 0;
//     for(UInt i=0, n=results_.GetSize(); i<n; i++) {
//       if(results_[i]->resultType == LAGRANGE_MULT) {
//         lmResultIdx = i;
//         break;
//       }
//     }
//     LOG_DBG2(heatcondpde) << "NonMatching: Index of LAGRANGE_MULT result: "
//                      << lmResultIdx;
    
//     for( UInt i = 0; i < ncIFaces_.GetSize(); i++ ) {
      
//       // get regionId of Lagrangian surface
//       StdVector<std::string> keyVec, attrVec, valVec;
//       std::string slaveSide;
//       std::string ncIfaceName = ptgrid_->GetRegion().ToString(ncIFaces_[i]);

//       PtrParamNode ncIfaceListNode;
//       ncIfaceListNode = param->Get("domain")->Get("ncInterfaceList");

//       slaveSide = ncIfaceListNode->
//           GetByVal("ncInterface", "name",
//                    ncIfaceName)->Get("slaveSide")->As<std::string>();

//       // Part 1: Define integrator M(Psi, Lambda) on
//       //         non-conforming interface
//       LOG_DBG2(heatcondpde) << "NonMatching: Defining nonconforming integrator"
//                         << " for M on interface '"
//                         << ptgrid_->GetRegion().ToString(ncIFaces_[i]) << "'.";
//       shared_ptr<ElemList> actNcList( new ElemList(ptgrid_ ) );
//       actNcList->SetRegion( ncIFaces_[i] );
      
//       NonConformingInt * ncInt = 
//         new NonConformingInt( 1, isaxi_ );

//       NcBiLinFormContext * stiffIntDescr = 
//      	  new NcBiLinFormContext( ncInt , STIFFNESS );

//       // Force assembling of M(Psi, Lambda)^T
//       stiffIntDescr->SetCounterPart( true );

//       stiffIntDescr->SetPtPdes(this, this);
//       stiffIntDescr->SetResults( results_[0], results_[lmResultIdx],
//                                  actNcList, actNcList );
      
//       assemble_->AddBiLinearForm( stiffIntDescr );


//       // Part 2: Define integrator D(Psi, Lambda) on
//       //         Lagrangian surface
//       LOG_DBG2(heatcondpde) << "NonMatching: Defining mass integrator"
//                         << " for D on interface '"
//                         << ptgrid_->GetRegion().ToString(ncIFaces_[i]) << "'.";
//       shared_ptr<SurfElemList> actSDList( new SurfElemList(ptgrid_ ) );      
//       actSDList->SetRegion( ptgrid_->GetRegion().Parse(slaveSide));

//       // D(Psi, Lambda) has the form of a standard mass
//       // integrator with factor 1.0
//       MassInt * dMatInt = new MassInt( 1.0, 1, isaxi_ );
//       BiLinFormContext * dMatContext = 
//         new BiLinFormContext( dMatInt, STIFFNESS );

//       // Force assembling of D(Psi, Lambda)^T
//       dMatContext->SetCounterPart( true );
//       dMatContext->SetPtPdes( this, this );
//       dMatContext->SetResults( results_[0], results_[lmResultIdx],
//                                actSDList, actSDList );
      
//       assemble_->AddBiLinearForm( dMatContext );

//       // Give result LAGRANGE_MULT to equation numbering class
//       eqnMap_->AddResult( *results_[lmResultIdx], actSDList );
//     }


  // TO BE DONE

//   try {
//     // iterate over all parameter nodes
//     for( UInt i = 0; i < rhsValuesNodes.GetSize(); i++ )
//     {
//       std::string rhsRegion(rhsValuesNodes[i]->Get("region")->As<std::string>());
//       rhsValuesNodes[i]->GetValue("isharmonic", isharmonic, ParamNode::EX);
//       //rhsFileId = rhsValuesNodes[i]->Get("inputId")->As<std::string>();

//       if ( isharmonic ) {
//         Info->PrintF( pdename_,
//             "Use same RHS source vector for region '%s' for all timesteps.\n\n",
//             rhsRegion.c_str() );
//       } else {
//         Info->PrintF( pdename_,
//             "Use RHS source vector for region '%s' from acoustic results.\n\n",
//             rhsRegion.c_str() );

//       }

//       // get material density
//       BaseMaterial * actMat = materials_[ ptgrid_->GetRegion().Parse(rhsRegion) ];
//       actMat->GetScalar(density,DENSITY,Global::REAL);

//       linAcouPowerSourceInt* sourceRHSInt = new linAcouPowerSourceInt( isaxi_,
//           isharmonic, rhsFileId, rhsRegion, density );

//       LinearFormContext* sourceRHSContext = new LinearFormContext( sourceRHSInt );
//       sourceRHSContext->SetPtPde( this );

//       shared_ptr<ElemList> rhsElemList( new ElemList(ptgrid_ ) );
//       rhsElemList->SetRegion( ptgrid_->GetRegion().Parse(rhsRegion) );
//       sourceRHSContext->SetResult( results_[0], rhsElemList );
//       assemble_->AddLinearForm( sourceRHSContext );
//     }
//   }
//   catch(Exception & ex)
//   {
//     RETHROW_EXCEPTION(ex, "Could not assemble RHS source integrator"
//         <<" in HeatPDE" );
//   }

}


  LinearFormContext* HeatPDE::CreateRhsLinearForm(SolutionType rhsType,
                                                      shared_ptr<CoefFunction > rhsCoef){
    LinearFormContext * mContext = NULL;
    switch(rhsType){
    case HEAT_SOURCE_DENSITY:
      BUIntegrator<IdentityOperator,FeH1,Double>* curInt;
      curInt = new BUIntegrator<IdentityOperator,FeH1,Double>(1.0,rhsCoef);

      mContext = new LinearFormContext(curInt);
      mContext->SetFeFunction( feFunctions_[HEAT_TEMPERATURE]);
      curInt->SetFeSpace( feFunctions_[HEAT_TEMPERATURE]->GetFeSpace());
      break;
    default:
      Exception("Right hand side quantity not known for heatPDE");
    }
    return mContext;
  }



void HeatPDE::DefineSolveStep() {

  solveStep_ = new StdSolveStep(*this);

}

// ======================================================
// TIME STEPPING SECTION
// ======================================================

void HeatPDE::InitTimeStepping() {

  PtrParamNode systemNode = FindLinearSystem(pdename_);
  // Until now no effective mass formulation in the trapezoidal
  //  integration scheme is implemented!
  TS_alg_ = new Trapezoidal( algsys_, systemNode );

}

void HeatPDE::CalcResults( shared_ptr<BaseResult> res ) {

  switch (res->GetResultInfo()->resultType ) {
  case HEAT_TEMPERATURE:
    feFunctions_[HEAT_TEMPERATURE]->ExtractResult( res );
    break;

//   case HEAT_RHS_LOAD:
//     if( isComplex_ ) {
//       ExtractRhsResult<Complex>( result, results_[0] );
//     } else {
//       ExtractRhsResult<Double>( result, results_[0] );
//     }
//     break;
    
  default:
    WARN( "Resulttype not computable by thermic PDE" );
  }
}



void HeatPDE::DefinePrimaryResults() {

  // === TEMPERATURE ===
  shared_ptr<ResultInfo> res1( new ResultInfo);
  res1->resultType = HEAT_TEMPERATURE;
    
  res1->dofNames = "";
  res1->unit = "K";
  res1->definedOn = ResultInfo::NODE;
  res1->entryType = ResultInfo::SCALAR;
  feFunctions_[HEAT_TEMPERATURE]->SetResultInfo(res1);
  results_.Push_back( res1 );
  availResults_.insert( res1 );


  // === TEMPERATURE RHS ===
  shared_ptr<ResultInfo> rhs ( new ResultInfo );
  rhs->resultType = HEAT_RHS_LOAD;
  rhs->dofNames = "";
  rhs->unit = "?";
  rhs->definedOn = results_[0]->definedOn;
  rhs->entryType = ResultInfo::SCALAR;
  availResults_.insert( rhs );
  postProcResults_[HEAT_RHS_LOAD] = HEAT_TEMPERATURE;
  

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
      lagr->definedOn = results_[0]->definedOn;
      results_.Push_back( lagr );
    } 

}

} // end of namespace CoupledField
