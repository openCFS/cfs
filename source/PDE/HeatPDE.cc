#include <boost/lexical_cast.hpp>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>

#include <typeinfo>

#include "HeatPDE.hh"

#include "General/defs.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Domain/CoefFunction/CoefFunction.hh"
#include "Domain/CoefFunction/CoefFunctionApprox.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefFunctionOpt.hh"
#include "Utils/StdVector.hh"

#include "Driver/Assemble.hh"
#include "Utils/ApproxData.hh"
#include "Utils/SmoothSpline.hh"
#include "Materials/Models/Hysteresis.hh"
#include "Materials/Models/Preisach.hh"
#include "FeBasis/H1/FeSpaceH1Nodal.hh"
#include "FeBasis/FeFunctions.hh"

#include "Driver/TimeSchemes/TimeSchemeGLM.hh"

#include "Driver/SolveSteps/StdSolveStep.hh"


#include "Driver/Assemble.hh"

//new integrator concept
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/BiLinWrappedLinForm.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/KXInt.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/LinForms/SingleEntryInt.hh"

//new postprocessing concept
#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/Results/ResultFunctor.hh"

namespace CoupledField {

  DECLARE_LOG(heatcondpde)
  DEFINE_LOG(heatcondpde, "pde.heatcond")

// ======================================================
// SET SOLUTION INFORMATION
// ======================================================
HeatPDE::HeatPDE(Grid * aptgrid, PtrParamNode paramNode,
                 PtrParamNode infoNode,
                 shared_ptr<SimState> simState, Domain* domain)
:SinglePDE( aptgrid, paramNode, infoNode, simState, domain ) {

  pdename_           = "heatConduction";
  pdematerialclass_  = THERMIC;
  nonLin_            = false;
  
  //! Always use updated Lagrangian formulation 
  updatedGeo_        = true;
}



void HeatPDE::ReadSpecialBCs() {

  // First, delete all of the Neumann boundary conditions object
  //inBcs_.Clear();

  // fetch paramnodes for Robin boundary condition
  ParamNodeList rbcNodes = myParam_->Get("bcsAndLoads")->GetList("robin");

  std::string myDof, myName, myType, myHTC, myBulkTemp;

  // iterate over all parameter nodes
  for( UInt i = 0; i < rbcNodes.GetSize(); i++ ) {
    try {
      rbcNodes[i]->GetValue( "name", myName );
      rbcNodes[i]->GetValue( "entityType", myType );
      rbcNodes[i]->GetValue( "heatTransferCoefficient", myHTC );
      rbcNodes[i]->GetValue( "bulkTemperature", myBulkTemp );

      // Create robin boundary condition
      shared_ptr<RobinBc> actBc ( new RobinBc );

      
      shared_ptr<EntityList> actList =
          ptGrid_->GetEntityList( EntityList::SURF_ELEM_LIST, myName); 

      actBc->entities = actList;
      actBc->result = results_[0];
      //      actBc->eqnMap = eqnMap_;
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
        FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1, ptGrid_);
      crSpaces[HEAT_TEMPERATURE]->Init(solStrat_);
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
 //  nonLinTotalFormulation_ = true;

  }


void HeatPDE::DefineIntegrators() {

  RegionIdType actRegion;
  BaseMaterial * actSDMat = NULL;  

  // convert to tensor type
  SubTensorType tensorType = FULL;
  if ( ptGrid_->GetDim() == 2 ) {
    if ( ptGrid_->IsAxi() ) {
      tensorType = AXI;
      isaxi_ = true;
    } else {
      tensorType = PLANE;
    }
  }

  // Define integrators for "standard" materials
  std::map<RegionIdType, BaseMaterial*>::iterator it;
  shared_ptr<FeSpace> mySpace = feFunctions_[HEAT_TEMPERATURE]->GetFeSpace();
  
  for ( it = materials_.begin(); it != materials_.end(); it++ ) {
    
    // Set current region and material
    actRegion = it->first;
    actSDMat = it->second;
    
    // Get current region name
    std::string regionName = ptGrid_->GetRegion().ToString(actRegion);
    
    // create new entity list
    shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
    actSDList->SetRegion( actRegion );
    
    // --- Set the FE ansatz for the current region ---
    PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
    std::string polyId = curRegNode->Get("polyId")->As<std::string>();
    std::string integId = curRegNode->Get("integId")->As<std::string>();
    mySpace->SetRegionApproximation(actRegion, polyId,integId);
    
    // pass entitylist of fespace / fefunction
    shared_ptr<BaseFeFunction> feFunc = feFunctions_[HEAT_TEMPERATURE];
    feFunc->AddEntityList( actSDList );

    // ====================================================================
    // stiffness integrator
    // ====================================================================
  
    //get possible nonlinearities defined in this region
    StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion]; 
    if ( nonLinTypes.Find(NLHEAT_CONDUCTIVITY) != -1 ) {
      // informs material that approx./interpol. for heat conductivity is needed
      //BaseBOperator * bOp = new IdentityOperator<FeH1>();
      PtrCoefFct heatCoef = this->GetCoefFct(HEAT_TEMPERATURE);
      PtrCoefFct condNL = actSDMat->GetScalCoefFncNonLin( HEAT_CONDUCTIVITY, Global::REAL, heatCoef);
                                      
      // create stiffness integrator
      BaseBDBInt* stiffInt = NULL;
      if( dim_ == 2 ) {
        stiffInt = new BBInt<>(new GradientOperator<FeH1,2>(), condNL,
                               1.0, updatedGeo_ );
      } else {
        stiffInt = new BBInt<>(new GradientOperator<FeH1,3>(), condNL,
                               1.0, updatedGeo_ );
      }
      stiffInt->SetName("StiffnessIntegrator-NL");

      BiLinFormContext* stiffContext = new BiLinFormContext(stiffInt, STIFFNESS );
      stiffContext->SetEntities( actSDList, actSDList );
      stiffContext->SetFeFunctions( feFunc, feFunc );

      assemble_->AddBiLinearForm( stiffContext );
      bdbInts_[actRegion] = stiffInt;

      // =================================
      //  Nonlinear RHS-integrator
      // =================================
//      LinearForm * rhsNlinForm = new KXIntegrator<Double>(stiffInt, -1.0,
//                                                          feFunc );
//      rhsNlinForm->SetName("RHSNonLinFormHeatStiff");
//      LinearFormContext * rhsNlinContext =
//          new LinearFormContext( rhsNlinForm );
//      rhsNlinContext->SetEntities( actSDList );
//      rhsNlinContext->SetFeFunction( feFunc );
//      assemble_->AddLinearForm( rhsNlinContext );
    }
    else {
      // --- linear real-valued stiffness integrator ---
      shared_ptr<CoefFunction > curCoef = actSDMat->GetTensorCoefFnc( HEAT_CONDUCTIVITY_TENSOR, tensorType, Global::REAL );

      // when we do optimization we wrap the original CoefFunction. Don't check for region to handle dim-1 pressure on dim elements
      if(domain->GetDesign(false) != NULL)
      {
        CoefFunctionOpt* tmpFnc = new CoefFunctionOpt(domain->GetDesign(), curCoef, this); // takes double and complex
        curCoef.reset(tmpFnc);
      }

      BaseBDBInt* stiffInt = NULL;
      if( dim_ == 2 ) {
        stiffInt = new BDBInt<>(new GradientOperator<FeH1,2>(), curCoef,1.0, updatedGeo_ );
      } else {
        stiffInt = new BDBInt<>(new GradientOperator<FeH1,3>(), curCoef,1.0, updatedGeo_ );
      }
      stiffInt->SetName("HeatConductivity");
      
      // the integrator has a coef function but for the optimization case the opt coef needs to know also the integrator
      if(domain->GetDesign(false) != NULL)
        dynamic_pointer_cast<CoefFunctionOpt>(curCoef)->SetForm(stiffInt);

      BiLinFormContext* stiffIntDescr = new BiLinFormContext(stiffInt, STIFFNESS );
          
      //stiffIntDescr->SetPtPdes(this, this);
      stiffIntDescr->SetEntities( actSDList, actSDList );
      stiffIntDescr->SetFeFunctions(feFunc,feFunc);
      stiffInt->SetFeSpace( feFunctions_[HEAT_TEMPERATURE]->GetFeSpace());
      
      assemble_->AddBiLinearForm( stiffIntDescr );
      bdbInts_[actRegion] = stiffInt;

//      if ( nonLinTypes.Find(NLHEAT_CONDUCTIVITY) ||
//                nonLinTypes.Find(NLHEAT_CAPACITY) != -1 ) {
//        // === Additional RHS integrator in case of Non-linearity ===
//        LinearForm * rhsNlinForm = new KXIntegrator<Double>(stiffInt, -1.0,
//                                                            feFunc );
//        rhsNlinForm->SetName("RHSNonLinFormHeatStiff-Lin");
//        LinearFormContext * rhsNlinContext =
//            new LinearFormContext( rhsNlinForm );
//        rhsNlinContext->SetEntities( actSDList );
//        rhsNlinContext->SetFeFunction( feFunc );
//        assemble_->AddLinearForm( rhsNlinContext );
//      }
    }

    // ====================================================================
    // mass integrator
    // ====================================================================

    if ( nonLinTypes.Find(NLHEAT_CAPACITY) != -1 ) {
      // Factor for mass matrix: density * heatCapacity
      PtrCoefFct heatCoef = this->GetCoefFct(HEAT_TEMPERATURE);
      PtrCoefFct density = actSDMat->GetScalCoefFnc( DENSITY, Global::REAL );
      //BaseBOperator * bOp = new IdentityOperator<FeH1>();
      PtrCoefFct capNL = 
          actSDMat->GetScalCoefFncNonLin( HEAT_CAPACITY, Global::REAL, heatCoef );

      PtrCoefFct nlMassCoeff =
          CoefFunction::Generate(mp_, Global::REAL,
                                 CoefXprBinOp(mp_,  capNL, density, CoefXpr::OP_MULT ) );

      // create stiffness integrator
      BaseBDBInt* massIntNL = NULL;
      if( dim_ == 2 ) {
        massIntNL = new BBInt<>(new IdentityOperator<FeH1,2>(), nlMassCoeff,
                                1, updatedGeo_ );
      } else {
        massIntNL = new BBInt<>(new IdentityOperator<FeH1,3>(), nlMassCoeff,
                                1, updatedGeo_ );
      }
      massIntNL->SetName("MassIntegrator-NL");

      BiLinFormContext * massNLContext =
        new BiLinFormContext(massIntNL, DAMPING );
      massNLContext->SetEntities( actSDList, actSDList );
      massNLContext->SetFeFunctions( feFunc, feFunc );

      assemble_->AddBiLinearForm( massNLContext );
      bdbInts_[actRegion] = massIntNL;
    }
    else {
      // Linear case
      // Factor for mass matrix: density * heatCapacity
      PtrCoefFct density = actSDMat->GetScalCoefFnc( DENSITY, Global::REAL );
      PtrCoefFct heatCapacity =
          actSDMat->GetScalCoefFnc( HEAT_CAPACITY, Global::REAL );
      PtrCoefFct massFactor =
          CoefFunction::Generate(mp_, Global::REAL,
                                 CoefXprBinOp( mp_, density, heatCapacity,
                                               CoefXpr::OP_MULT ) );

      BiLinearForm *massInt = NULL;
      if(dim_==2)
        massInt = new BBInt<>(new IdentityOperator<FeH1,2,1,Double>(), massFactor,1.0, updatedGeo_ );
      else
        massInt = new BBInt<>(new IdentityOperator<FeH1,3,1,Double>(), massFactor,1.0, updatedGeo_ );

      massInt->SetName("MassIntegrator");
      massInt->SetFeSpace( feFunctions_[HEAT_TEMPERATURE]->GetFeSpace() );
      
      BiLinFormContext *massContext =  new BiLinFormContext(massInt, DAMPING );
      
      massContext->SetEntities( actSDList, actSDList );
      massContext->SetFeFunctions( feFunctions_[HEAT_TEMPERATURE],feFunctions_[HEAT_TEMPERATURE]);
      assemble_->AddBiLinearForm( massContext );
    }
  }

  // ===============
  //  electric power density
  // ===============
  LOG_DBG(heatcondpde) << "Reading electric power densities";

  shared_ptr<BaseFeFunction> myFct = feFunctions_[HEAT_TEMPERATURE];
  StdVector<std::string> dispDofNames = myFct->GetResultInfo()->dofNames;
  StdVector<shared_ptr<EntityList> > ent;
  StdVector<PtrCoefFct > coef;
  LinearForm * lin = NULL;

  ReadRhsExcitation( "elecPowerDensity", dispDofNames, ResultInfo::SCALAR, isComplex_,
                      ent, coef, updatedGeo_ );
  for( UInt i = 0; i < ent.GetSize(); ++i ) {
    // check type of entitylist
    if (ent[i]->GetType() == EntityList::NODE_LIST) {
      EXCEPTION("Electric power density must be defined on elements")
    }

    if( dim_ == 2) {
      if(isComplex_) {
        lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1,2>(),
                                          Complex(1.0), coef[i], updatedGeo_ );
      } else {
        lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1,2>(),
                                         1.0, coef[i], updatedGeo_ );
      }
    } else  {
      if(isComplex_) {
        lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1,3>(),
                                          Complex(1.0), coef[i], updatedGeo_ );
      } else {
        lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1,3>(),
                                         1.0, coef[i], updatedGeo_ );
      }
    }
    lin->SetName("ElectricPowerDensityInt");
    LinearFormContext *ctx = new LinearFormContext( lin );
    ctx->SetEntities( ent[i] );
    ctx->SetFeFunction(myFct);
    assemble_->AddLinearForm(ctx);
    myFct->AddEntityList(ent[i]);
  } // for

  // ======================================================================
  // Heat flux boundary condition
  // ======================================================================
  ReadRhsExcitation( "heatFlux", dispDofNames, ResultInfo::SCALAR, isComplex_,
                      ent, coef, updatedGeo_ );
  for( UInt i = 0; i < ent.GetSize(); ++i ) {
    // check type of entitylist
    if (ent[i]->GetType() == EntityList::NODE_LIST) {
      EXCEPTION("Heat flux must be defined on elements")
    }

    if( dim_ == 2) {
      if(isComplex_) {
        lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1,2>(),
                                          Complex(1.0), coef[i], updatedGeo_ );
      } else {
        lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1,2>(),
                                         1.0, coef[i], updatedGeo_ );
      }
    } else  {
      if(isComplex_) {
        lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1,3>(),
                                          Complex(1.0), coef[i], updatedGeo_ );
      } else {
        lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1,3>(),
                                         1.0, coef[i], updatedGeo_ );
      }
    }
    lin->SetName("HeatFluxInt");
    LinearFormContext *ctx = new LinearFormContext( lin );
    ctx->SetEntities( ent[i] );
    ctx->SetFeFunction(myFct);
    assemble_->AddLinearForm(ctx);
    myFct->AddEntityList(ent[i]);
  } // for


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
//       std::string ncIfaceName = ptGrid_->GetRegion().ToString(ncIFaces_[i]);

//       PtrParamNode ncIfaceListNode;
//       ncIfaceListNode = domain_->GetParamRoot()->Get("domain")->Get("ncInterfaceList");

//       slaveSide = ncIfaceListNode->
//           GetByVal("ncInterface", "name",
//                    ncIfaceName)->Get("slaveSide")->As<std::string>();

//       // Part 1: Define integrator M(Psi, Lambda) on
//       //         non-conforming interface
//       LOG_DBG2(heatcondpde) << "NonMatching: Defining nonconforming integrator"
//                         << " for M on interface '"
//                         << ptGrid_->GetRegion().ToString(ncIFaces_[i]) << "'.";
//       shared_ptr<ElemList> actNcList( new ElemList(ptGrid_ ) );
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
//                         << ptGrid_->GetRegion().ToString(ncIFaces_[i]) << "'.";
//       shared_ptr<SurfElemList> actSDList( new SurfElemList(ptGrid_ ) );      
//       actSDList->SetRegion( ptGrid_->GetRegion().Parse(slaveSide));

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
//       BaseMaterial * actMat = materials_[ ptGrid_->GetRegion().Parse(rhsRegion) ];
//       actMat->GetScalar(density,DENSITY,Global::REAL);

//       linAcouPowerSourceInt* sourceRHSInt = new linAcouPowerSourceInt( isaxi_,
//           isharmonic, rhsFileId, rhsRegion, density );

//       LinearFormContext* sourceRHSContext = new LinearFormContext( sourceRHSInt );
//       sourceRHSContext->SetPtPde( this );

//       shared_ptr<ElemList> rhsElemList( new ElemList(ptGrid_ ) );
//       rhsElemList->SetRegion( ptGrid_->GetRegion().Parse(rhsRegion) );
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

void HeatPDE::DefineNcIntegrators() {
  StdVector< NcInterfaceInfo >::iterator ncIt = ncInterfaces_.Begin(),
                                         endIt = ncInterfaces_.End();
  for ( ; ncIt != endIt; ++ncIt ) {
    switch (ncIt->type) {
    case NC_MORTAR:
      if (dim_ == 2)
        DefineMortarCoupling<2,1>(HEAT_TEMPERATURE, *ncIt);
      else
        DefineMortarCoupling<3,1>(HEAT_TEMPERATURE, *ncIt);
      break;
    case NC_NITSCHE:
      if (dim_ == 2)
        DefineNitscheCoupling<2,1>(HEAT_TEMPERATURE, *ncIt );
      else
        DefineNitscheCoupling<3,1>(HEAT_TEMPERATURE, *ncIt );
      break;
    default:
      EXCEPTION("Unknown type of ncInterface");
      break;
    }
  }
}

void HeatPDE::DefineRhsLoadIntegrators() {
  
  LOG_TRACE(heatcondpde) << "Defining rhs load integrators for thermal PDE";

    // Get FESpace and FeFunction of electric potential
    shared_ptr<BaseFeFunction> myFct = feFunctions_[HEAT_TEMPERATURE];
    shared_ptr<FeSpace> mySpace = myFct->GetFeSpace();

    StdVector<shared_ptr<EntityList> > ent;
    StdVector<PtrCoefFct > coef;
    LinearForm * lin = NULL;
    StdVector<std::string> dofNames;

    
    // @Manfred: Is there an equivalent to total charge (= nodal values)
    // for the thtermal PDE? This would go here.
//    // =========================
//    //  Charges (volume, nodal)
//    // =========================
//    LOG_DBG(elecpde) << "Reading charges";
//    ReadRhsExcitation( "charge", dofNames, ResultInfo::SCALAR, 
//                       isComplex_, ent, coef );
//
//    for( UInt i = 0; i < ent.GetSize(); ++i ) {
//      // check type of entitylist
//      if (ent[i]->GetType() == EntityList::NODE_LIST) {
//
//        // ---------------
//        //  Nodal Charges 
//        // ---------------
//        // Nodal charge must be constant
//        if( coef[i]->GetDependency() == CoefFunction::GENERAL ) {
//          EXCEPTION("Nodal charges must not be spatial dependent");
//        }
//
//        UInt numNodes = ent[i]->GetSize();
//        if( numNodes > 1 ) {
//          // Here we would divide the nodal force by the number of nodes
//          // in the list, in order to ensure that the whole force corresponds
//          // to the prescribed value. However, this requires modification of 
//          // the expressions of the coefficient functions, which depends on real/harm
//          // and the specific type (const, timefreq, variable).
//          WARN("The chareg value will not be divided by the number of nodes and thus "
//              << "depends on the number of nodes" );
//        }
//
//        lin = new SingleEntryInt(coef[i]);
//        lin->SetName("NodalChargeInt");
//        LinearFormContext *ctx = new LinearFormContext( lin );
//        ctx->SetEntities( ent[i] );
//        ctx->SetFeFunction(myFct);
//        assemble_->AddLinearForm(ctx);
//      } else {
//        // --------------------------
//        //  Surface / Volume Charges 
//        // --------------------------
//        EXCEPTION("Not yet implemented");
//
//        // Same issue here as above: We need to "divide" the total force by the
//        // area / volume to get the force density.
//      }
//    } // for

    bool coefUpdateGeo = true;
    // =====================
    //  HEAT SOURCE DENSITY
    // =====================
    LOG_DBG(heatcondpde) << "Reading heat source density";
    
    ReadRhsExcitation( "heatSourceDensity", dofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("Heat source density must be defined on elements")
      }
      EntityIterator it = ent[i]->GetIterator();
      it.Begin();
      
      if(isComplex_) {
        lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1>(), Complex(1.0), coef[i], coefUpdateGeo);
      } else  {
        lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1>(), 1.0, coef[i], coefUpdateGeo);
      }
      lin->SetName("HeatSourceDensityInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
    } // for

    // ========================
    //  HEAT SOURCE(nodal)
    // ========================
    LOG_DBG(heatcondpde) << "Reading heat source values";

    coefUpdateGeo = false;

    ReadRhsExcitation( "heatSource", dofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo);
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // assume that we have elem list due to specification of a region instead of named nodes in xml file
      if (ent[i]->GetType() != EntityList::NODE_LIST && ent[i]->GetType() != EntityList::ELEM_LIST) {
        EXCEPTION("Heat source must be defined on nodes!")
      }

      UInt numNodes = ent[i]->GetSize();
      if( numNodes > 1 ) {
        Global::ComplexPart part = Global::REAL;
        coef[i] = CoefFunction::Generate(mp_, part, CoefXprVecScalOp(mp_, coef[i], boost::lexical_cast<std::string>(numNodes), CoefXpr::OP_DIV) );
      }

      lin = new SingleEntryInt(coef[i]);
      lin->SetName("NodalHeatInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
    }

    // ========================
    //  DESIGN DEPENDENT HEAT SOURCE
    // ========================
    LOG_DBG(heatcondpde) << "Reading heat source values (design dependent)";

//    ReadRhsExcitation( "designDependentHeatSource", dofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo);
    ReadRhsExcitation( "designDependentHeatSource", dofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo);
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // assume that we have elem list due to specification of a region instead of named nodes in xml file
      if (ent[i]->GetType() != EntityList::NODE_LIST && ent[i]->GetType() != EntityList::ELEM_LIST) {
        EXCEPTION("Design dependent heat source must be defined on nodes!")
      }

      if(domain->GetDesign(false) != NULL)
      {
        CoefFunctionOpt* tmpFnc = new CoefFunctionOpt(domain->GetDesign(), coef[i], this); // takes double and complex
        coef[i].reset(tmpFnc);
//        coef[i]->SetConservative(true);
//        rhsFeFunctions_[HEAT_TEMPERATURE]->AddLoadCoefFunction(coef[i],ent[i]);
      }

      lin = new SingleEntryInt(coef[i]);
      lin->SetName("HeatConductivity");

      BiLinWrappedLinForm* linWrapped = new BiLinWrappedLinForm(lin,false); // have to wrapp lin form since CoefFunctionOpt expects a bilinform

      // the integrator has a coef function but for the optimization case the opt coef needs to know also the integrator
      if(domain->GetDesign(false) != NULL)
        dynamic_pointer_cast<CoefFunctionOpt>(coef[i])->SetForm(linWrapped);

      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities(ent[i]);
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
    }
}


void HeatPDE::DefineSolveStep() {

  solveStep_ = new StdSolveStep(*this);

}

// ======================================================
// TIME STEPPING SECTION
// ======================================================

void HeatPDE::InitTimeStepping() {

  // Until now no effective mass formulation in the trapezoidal
  //  integration scheme is implemented!
  Double gamma = 0.5; 
  GLMScheme * scheme = new Trapezoidal(gamma);

  if ( nonLinTotalFormulation_ ) {
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(scheme, 0)); //, nlType) );
    feFunctions_[HEAT_TEMPERATURE]->SetTimeScheme(myScheme);
  }
  else {
    TimeSchemeGLM::NonLinType nlType = (nonLin_)? TimeSchemeGLM::INCREMENTAL : TimeSchemeGLM::NONE;
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(scheme, 0, nlType) );
    feFunctions_[HEAT_TEMPERATURE]->SetTimeScheme(myScheme);
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
  res1->SetFeFunction(feFunctions_[HEAT_TEMPERATURE]);
  DefineFieldResult( feFunctions_[HEAT_TEMPERATURE], res1 );


  // -----------------------------------
  //  Define xml-names of Dirichlet BCs
  // -----------------------------------
  idbcSolNameMap_[HEAT_TEMPERATURE] = "temperature";
  
  // === TEMPERATURE RHS ===
//  shared_ptr<ResultInfo> rhs ( new ResultInfo );
//  rhs->resultType = HEAT_RHS_LOAD;
//  rhs->dofNames = "";
//  rhs->unit = "?";
////  rhs->definedOn = results_[0]->definedOn;
//  rhs->definedOn = ResultInfo::NODE;
//  rhs->entryType = ResultInfo::SCALAR;
//  availResults_.insert( rhs );
//  rhsFeFunctions_[HEAT_TEMPERATURE]->SetResultInfo(rhs);
//  DefineFieldResult( rhsFeFunctions_[HEAT_TEMPERATURE], rhs );

}

void HeatPDE::DefinePostProcResults() {
  shared_ptr<BaseFeFunction> feFct = feFunctions_[HEAT_TEMPERATURE];

  if ( analysistype_ != STATIC ) {
    // === TEMPERATURE D1===
    shared_ptr<ResultInfo> heatD1( new ResultInfo);
    heatD1->resultType = HEAT_TEMPERATURE_D1;

    heatD1->dofNames = "";
    heatD1->unit = "K/s";
    heatD1->definedOn = ResultInfo::NODE;
    heatD1->entryType = ResultInfo::SCALAR;
    availResults_.insert( heatD1 );
    DefineTimeDerivResult( HEAT_TEMPERATURE_D1, 1, HEAT_TEMPERATURE );
  }
  
  // === TEMPERATURE RHS ===
  shared_ptr<ResultInfo> rhs ( new ResultInfo );
  rhs->resultType = HEAT_RHS_LOAD;
  rhs->dofNames = "";
  rhs->unit = "K";
  rhs->definedOn = ResultInfo::NODE;
  rhs->entryType = ResultInfo::SCALAR;
  rhsFeFunctions_[HEAT_TEMPERATURE]->SetResultInfo(rhs);
  DefineFieldResult( rhsFeFunctions_[HEAT_TEMPERATURE], rhs );

  // === HEAT FLUX DENSITY ===
  shared_ptr<ResultInfo> flux ( new ResultInfo );
  flux->resultType = HEAT_FLUX_DENSITY;
  flux->SetVectorDOFs(dim_, isaxi_);
  flux->unit = "W/m^2";
  flux->definedOn = ResultInfo::ELEMENT;
  flux->entryType = ResultInfo::VECTOR;
  shared_ptr<CoefFunctionFormBased> fluxFunc;
  if( isComplex_ ) {
    fluxFunc.reset(new CoefFunctionFlux<Complex>(feFct, flux, Complex(-1.0)));
  } else {
    fluxFunc.reset(new CoefFunctionFlux<Double>(feFct, flux, -1.0));
  }
  DefineFieldResult( fluxFunc, flux );
  stiffFormCoefs_.insert(fluxFunc);

  // optimization results are provided in DesignSpace::ExtractResults()
  // copied from MechPDE
  // === MECH_PSEUDO_DENISTY ===
  shared_ptr<ResultInfo> mpd(new ResultInfo);
  mpd->resultType = MECH_PSEUDO_DENSITY;
  mpd->entryType = ResultInfo::SCALAR;
  mpd->definedOn = ResultInfo::ELEMENT;
  mpd->dofNames = "";
  mpd->fromOptimization = true;
  DefineFieldResult(shared_ptr<FeFunction<double> >(new FeFunction<double>(NULL)), mpd); // the fe-function is only a dummy

  // === PHYSICAL_PSEUDO_DENISTY ===
  shared_ptr<ResultInfo> ppd(new ResultInfo);
  ppd->resultType = PHYSICAL_PSEUDO_DENSITY;
  ppd->entryType = ResultInfo::SCALAR;
  ppd->definedOn = ResultInfo::ELEMENT;
  ppd->dofNames = "";
  ppd->fromOptimization = true;
  DefineFieldResult(shared_ptr<FeFunction<double> >(new FeFunction<double>(NULL)), ppd);
}

} // end of namespace CoupledField
