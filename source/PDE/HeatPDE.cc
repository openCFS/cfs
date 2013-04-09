#include <boost/lexical_cast.hpp>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>

#include "HeatPDE.hh"

#include "General/defs.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Domain/CoefFunction/CoefFunction.hh"
#include "Domain/CoefFunction/CoefFunctionApprox.hh"
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
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/KXInt.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"

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
      tensorType = PLANE_STRAIN;
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
      std::cout << "Do NL Cond" << std::endl;

      BaseBOperator * bOp = new IdentityOperator<FeH1>();
      PtrCoefFct condNL = 
          actSDMat->GetScalCoefFncNonLin( HEAT_CONDUCTIVITY, Global::REAL, feFunc, bOp);
                                      
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

      BiLinFormContext * stiffContext =
        new BiLinFormContext(stiffInt, STIFFNESS );
      stiffContext->SetEntities( actSDList, actSDList );
      stiffContext->SetFeFunctions( feFunc, feFunc );

      assemble_->AddBiLinearForm( stiffContext );
      bdbInts_[actRegion] = stiffInt;

      // =================================
      //  Nonlinear RHS-integrator
      // =================================
      LinearForm * rhsNlinForm = new KXIntegrator<Double>(stiffInt, -1.0, 
                                                          feFunc );
        rhsNlinForm->SetName("RHSNonLinFormHeatStiff");
        LinearFormContext * rhsNlinContext =
            new LinearFormContext( rhsNlinForm );
        rhsNlinContext->SetEntities( actSDList );
        rhsNlinContext->SetFeFunction( feFunc );
        assemble_->AddLinearForm( rhsNlinContext );
    }
    else {
      // --- linear real-valued stiffness integrator ---
      shared_ptr<CoefFunction > curCoef = 
        actSDMat->GetTensorCoefFnc( HEAT_CONDUCTIVITY, tensorType, 
                                    Global::REAL );

      BaseBDBInt* stiffInt = NULL;
      if( dim_ == 2 ) {
        stiffInt = new BDBInt<>(new GradientOperator<FeH1,2>(), curCoef,1.0, updatedGeo_ );
      } else {
        stiffInt = new BDBInt<>(new GradientOperator<FeH1,3>(), curCoef,1.0, updatedGeo_ );
      }
      stiffInt->SetName("StiffnessIntegrator");
      
      BiLinFormContext * stiffIntDescr =
        new BiLinFormContext(stiffInt, STIFFNESS );
          
      //stiffIntDescr->SetPtPdes(this, this);
      stiffIntDescr->SetEntities( actSDList, actSDList );
      stiffIntDescr->SetFeFunctions(feFunc,feFunc);
      stiffInt->SetFeSpace( feFunctions_[HEAT_TEMPERATURE]->GetFeSpace());
      
      assemble_->AddBiLinearForm( stiffIntDescr );
      bdbInts_[actRegion] = stiffInt;
    }

    // ====================================================================
    // mass integrator
    // ====================================================================

    // Factor for mass matrix: density * heatCapacity
    PtrCoefFct density = actSDMat->GetScalCoefFnc( DENSITY, Global::REAL );
    PtrCoefFct heatCapacity = 
        actSDMat->GetScalCoefFnc( HEAT_CAPACITY, Global::REAL );
    PtrCoefFct massFactor =
        CoefFunction::Generate(mp_, Global::REAL,
                               CoefXprBinOp( mp_, density, heatCapacity, 
                                             CoefXpr::OP_MULT ) );
    

    if ( nonLinTypes.Find(NLHEAT_CAPACITY) != -1 ) {
      // informs material that approx./interpol. for heat conductivity is needed
      std::cout << "Do NL Capacity" << std::endl;
//      
      BaseBOperator * bOp = new IdentityOperator<FeH1>();
      PtrCoefFct capNL = 
          actSDMat->GetScalCoefFncNonLin( HEAT_CAPACITY, Global::REAL,
                                          feFunc, bOp);

      PtrCoefFct nlMassCoeff = 
          CoefFunction::Generate(mp_, Global::REAL,
                                 CoefXprBinOp(mp_,  capNL, density, CoefXpr::OP_MULT ) );

      // create stiffness integrator
      BaseBDBInt* massIntNL = NULL;
      if( dim_ == 2 ) {
        massIntNL = new BBInt<>(new IdentityOperator<FeH1,2>(), nlMassCoeff, 
                                1.0, updatedGeo_ );
      } else {
        massIntNL = new BBInt<>(new IdentityOperator<FeH1,3>(), nlMassCoeff, 
                                1.0, updatedGeo_ );
      }
      massIntNL->SetName("MassIntegrator-NL");

      BiLinFormContext * massNLContext =
        new BiLinFormContext(massIntNL, DAMPING );
      massNLContext->SetEntities( actSDList, actSDList );
      massNLContext->SetFeFunctions( feFunc, feFunc );

      assemble_->AddBiLinearForm( massNLContext );
      bdbInts_[actRegion] = massIntNL;

      // =================================
      //  Nonlinear RHS-integrator
      // =================================
      LinearForm * rhsNlinForm = new KXIntegrator<Double>(massIntNL, -1.0, feFunc );
      rhsNlinForm->SetName("RHSNonLinFormHeatMass");
      LinearFormContext * rhsNlinContext =
        new LinearFormContext( rhsNlinForm );
      rhsNlinContext->SetEntities( actSDList );
      rhsNlinContext->SetFeFunction( feFunc );
      assemble_->AddLinearForm( rhsNlinContext );
    }
    else {
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
    
    ReadRhsExcitation( "heatSourceDensity", dofNames, 
                       ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("Heat source density must be defined on elements")
      }
      EntityIterator it = ent[i]->GetIterator();
      it.Begin();
      
      if(isComplex_) {
        lin = new BUIntegrator<IdentityOperator<FeH1>, Complex>(Complex(1.0), coef[i], coefUpdateGeo);
      } else  {
        lin = new BUIntegrator<IdentityOperator<FeH1>, Double>(1.0, coef[i], coefUpdateGeo);
      }
      lin->SetName("HeatSourceDensityInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
    } // for
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
  //TS_alg_ = new Trapezoidal( algsys_, olasNode_ );
  shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(GLMScheme::TRAPEZOIDAL, 0) );

  feFunctions_[HEAT_TEMPERATURE]->SetTimeScheme(myScheme);


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
  shared_ptr<ResultInfo> rhs ( new ResultInfo );
  rhs->resultType = HEAT_RHS_LOAD;
  rhs->dofNames = "";
  rhs->unit = "?";
  rhs->definedOn = results_[0]->definedOn;
  rhs->entryType = ResultInfo::SCALAR;
  availResults_.insert( rhs );
  rhsFeFunctions_[HEAT_TEMPERATURE]->SetResultInfo(rhs);
  DefineFieldResult( rhsFeFunctions_[HEAT_TEMPERATURE], rhs );

  // ===================================
  // Check for non-conforming interfaces
  // ===================================
  StdVector<std::string> ncIfaceNames, ncIfaceNamesForPDE;
    StdVector<RegionIdType> ncIfaceIds;
    
    LOG_DBG2(heatcondpde) << "NonMatching: Checking if nonconforming "
                      << "interfaces of PDE exist in domain.";

    PtrParamNode heatcondpdeNCIfaceListNode;
    heatcondpdeNCIfaceListNode = domain_->GetParamRoot()->GetByVal("sequenceStep", std::string("index"), sequenceStep_)
    ->Get("pdeList/heatConduction/ncInterfaceList", ParamNode::PASS);
    
    if(!heatcondpdeNCIfaceListNode)
      return;

    PtrParamNode domainNCIfaceListNode;
    domainNCIfaceListNode = domain_->GetParamRoot()->Get("domain")->Get("ncInterfaceList", ParamNode::PASS);

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
    ptGrid_->GetRegion().Parse(ncIfaceNamesForPDE, ncIfaceIds);

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
