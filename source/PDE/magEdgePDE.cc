// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>

#include "magEdgePDE.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Utils/Coil.hh"
#include "Utils/SmoothSpline.hh"
#include "Utils/LinInterpolate.hh"
#include "Forms/curlfieldop.hh"
#include "Forms/curlCurlEdgeInt.hh"
#include "Forms/nLincurlCurlEdgeInt.hh"
#include "Forms/massEdgeInt.hh"
#include "Forms/linearForm.hh"
#include "trapezoidal.hh"
#include "CoupledPDE/pdecoupling.hh"
#include "Domain/ansatzFct.hh"
#include "Driver/assemble.hh"
#include "Utils/coordSystem.hh"
#include "Elements/fespaceHCurlHi.hh"
#include "Elements/HCurlElems.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "OLAS/algsys/standardsys.hh"

#ifdef USE_SCRIPTING
#include "DataInOut/Scripting/cfsmessenger.hh"
#endif

namespace CoupledField {

// declare class specific logging stream
DECLARE_LOG(magEdgePde)
DEFINE_LOG(magEdgePde, "magEdgePde")


  // **************
  //  Constructor
  // **************
  MagEdgePDE::MagEdgePDE( Grid * aptgrid, PtrParamNode paramNode )
    :SinglePDE( aptgrid, paramNode ) {

    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_          = "magneticEdge";
    pdematerialclass_ = ELECTROMAGNETIC;
    maxTimeDerivOrder_ = 1;

    // check if we have a 3d setup
    bool is3d = param->Get("domain")->Get("geometryType")->As<std::string>() == "3d";
    if ( !is3d )
      EXCEPTION("MagEdgePDE is just implemented for 3D setups!");
  
    // read in solution strategy for this PDE
    PtrParamNode  solStratNode = myParam_->Get("solutionStrategy",ParamNode::PASS);
    if( solStratNode )
      solStrategy_  = SolStrategyEnum.Parse(solStratNode->As<std::string>());

  }


  // *************
  //  Destructor
  // *************
  MagEdgePDE::~MagEdgePDE() {
  }


  // *********************************
  //  Re-Initialize PDE
  // *********************************
  void MagEdgePDE::SetSolutionStep(UInt stepNum ) {
    solStep_ = stepNum;
    if( solStrategy_ == STRAT_TWO_LEVEL && solStep_ == 2 ) {

      // First, delete all stuff
      // => Take code of destructors of SinglePDE and StdPDE
      // => make sure to delete the algebraic system as well



      // delete algebraic system
      // delete Assemble class
      // delete Time-Stepping class
      // delete timestepping *** not allowed here
      // delete solVec_, rhsVec_
      // delete materials_



      // Then, re-initialize the PDE (maybe we can adapt the init-method)
      // (called from Domain::InitPDEs())
      // - define analysis type
      // - read in regions 
      // - call DefineFEFunctions() *** Make sure in second step, to
      // - call DefineAvailResults()
      // - create new Standardsystem (hard coded -> change top SBM!)
      // - create Assemble class
      // - create new NodeStoreSol obejct and solution vector
      // - call ReadMaterialData()
      // - call ReadDampingInformation()
      // - call InitNonLin()
      // - determine "usePenalty_"
      // - call ReadBCs()
      // - call ReaSpecialBCs()
      // - call InitTimeStepping()
      // - call DefineIntegrators()
      // - Loop over all FSpaces and initialize (map equations)
      // - Initialize StoreSolution classes (numUnknowns, feSpace, regions)
      // - Initialize solPrev_ objects
      // - Initialize solVec_ and rhsVec_ objects
      // - Initialize TimeStepping algorithm
      // - call ReadStoreResults()
      // - call ReadSpecialResults()
      // - call DefineSolveStep() **** Prevent in 2nd step ****
      // - call PreparePDE4Computation()


      // Things performed in DefineAlgSys() (called from Domain::InitPDEs())
      // - call ReadOLASParams()  *** Make sure to use iterative system
      // - call olas->GraphSetupInit()
      // - call olas->RegisterPDE() *** may have to be adjusted for multiple spaces
      // - call solveStep->SetPDEId() *** may have to be adjusted
      // - call algssys->AssembleInit()
      // - call assemble_->SetupMatrixGraph()
      // - call algsys->AssembleDone()
      // - call algsys->GetReordering()
      // - pass reordering to FeFunctions / FeSpaces
      // - collect number of dirichletBCs from Spaces
      // - call algsys->SetNumDirichletBCs
      // - call SinglePDE::CreateMatricesSolver()
      // - call SetInitialCondition()
      // - Pass the algebraic system to all FeFunctions

      // Thins performed in StdPDE::CreateMatricesSolver:
      // - algsys->CreateLinSys()
      // - algsys->InitMatrix()
      // - algsys->CreateSolver
      // - algsys->CreatePrecond
      // - algsys->InitRHS()
      // -algsys->InitSol()


      // Alternative Approach is to keep as many objects as possible:
      // - re-initialize FeSpace and FeFunctions
      // - delete algebraic system and re-create it
      // - re-initialize StoreSol-objects (check, if really needed)
      // - re-size solPvec_, solVec_, rhsVec_ etc.
      // - Maybe re-initialize timsteppeing
      // - Create new algebraic system (and pass it to solvestep)
      // 

      // ----------------------------------------
      
      // delete old algebraic system and create new one
      REFACTOR;
      
//      std::cerr << "algsys_ was : " << algsys_;
//      delete algsys_;
//      algsys_ = NULL;
//      std::cerr << "algsys_ after delete is  " << algsys_ << std::endl;
//      algsys_ = new StandardSystem(FindLinearSystem("magneticEdge2"));
//      std::cerr << "algsys_ is : " << algsys_;
//      // Get parameter and report object of OLAS
//      olasParams_ = algsys_->GetOLASParams();
//      olasReport_ = algsys_->GetOLASReport();
//
//      // Obtain unique pde identifier
//      pdeId_ = algsys_->ObtainPDEId( pdename_ );
//      assemble_->SetAlgSys(algsys_);
//      assemble_->ResetMatrixReassembly();
//      
//      // Initialize FeSpace objects
//      numPdeEquations_ = 0;
//      numPdeUnknowns_ = 0;
//      numPdeInHomDirBc_ = 0;
//      //so this implementation is just for the magnetic potential
//      shared_ptr<FeSpace> actSpace = feFunctions_[MAG_POTENTIAL]->GetFeSpace();
//
//      // IMPORTANT: Set 2nd step in two-Level scheme
//      actSpace->SetStrategy(STRAT_TWO_LEVEL, 2);
//      //actSpace->Finalize();
//      actSpace->PreCalcShapeFncs();
//      numPdeEquations_ += actSpace->GetNumEquations();
//      numPdeUnknowns_ += actSpace->GetNumFreeEquations();
//      numPdeInHomDirBc_ += actSpace->GetNumInhomDirichletBc();
//      
//      // Re-set store-solution and other vector objects
//      sol_->Init();
//      solVec_->Resize(numPdeUnknowns_);
//      solVec_->Init();
//      rhsVec_->Resize(numPdeUnknowns_);
//      rhsVec_->Init();
//      sol_->SetAlgSysDataPointer(solVec_->GetSize(),
//                                         dynamic_cast<Vector<Double>&>(*solVec_).GetPointer() );
//      
//      DefineAlgSys();
//      //CreateMatrices_Solver();
    }
  }
  
  // **********************************
  //  Read special boundary conditions
  // **********************************

  void MagEdgePDE::ReadSpecialBCs() {


    // --------------------------------------------------------------------
    //   Get information about coils and open files for measurement coils
    // --------------------------------------------------------------------
    ReadCoils();

    // -----------------------------
    // Check for permanent magnets
    // -----------------------------
    ReadMagnets();

  }


  // ****************************
  //  Initialize Nonlinearities
  // ****************************
  void MagEdgePDE::InitNonLin() {

    nonLin_ = false;
    isHysteresis_ = false;

    // Check, if "nonLinList" is present
    PtrParamNode nonLinListNode = myParam_->Get("nonLinList", ParamNode::PASS );
    if( !nonLinListNode)
      return;

    // Get nonlinear types
    ParamNodeList nonLinNodes = nonLinListNode->GetChildren();
    for( UInt i = 0; i < nonLinNodes.GetSize(); i++ ) {

      std::string actTypeString = nonLinNodes[i]->GetName();
      std::string actId = nonLinNodes[i]->Get("id")->As<std::string>();

      NonLinType actType;
      String2Enum( actTypeString, actType );

      // check type
      if( actType == PERMEABILITY ) {
        nonLin_ = true;
      }
      if( actType == HYSTERESIS ) {
        isHysteresis_ = true;
      }
      nonLinIdType_[actId] = actType;
    }

    // Run over all region and set entry in "regionNonLinId"
    ParamNodeList regionNodes =
      myParam_->Get("regionList")->GetChildren();

    RegionIdType actRegionId;
    std::string actRegionName, actNonLinId;

    if( regionNodes.GetSize() > 0 ) {
      Info->PrintF( pdename_, "Non-linearity in following region(s)\n" );
    }
    for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {

      // get data
      regionNodes[i]->GetValue( "name", actRegionName );
      regionNodes[i]->GetValue( "nonLinId", actNonLinId );

      if( actNonLinId == "" )
        continue;

      actRegionId = ptgrid_->GetRegion().Parse( actRegionName );

      // Check nonLinId was already registerd
      if( nonLinIdType_.find( actNonLinId) == nonLinIdType_.end() ) {
        EXCEPTION( "NonLinearity with id '" << actNonLinId
                   << "' was not defined in 'nonLinList'" );
      }

      regionNonLinId_[actRegionId] = actNonLinId;
      regionNonLinType_[actRegionId] = nonLinIdType_[actNonLinId];

      // Log to info file
      std::string nonLinString;
      Enum2String( nonLinIdType_[actNonLinId], nonLinString );
      Info->PrintF( pdename_, " %s: %s\n", actRegionName.c_str(),
                    nonLinString.c_str() );

    }

    // set nonlinearity flag only, if any region references
    // a nonlinearity at all
    if( regionNonLinId_.size() > 0 ) {
      nonLin_ = true;
    }

    // Here we need in addition the nonLinMethod_ for the definition
    // of the integrators
    nonLinMethod_ = FIXEDPOINT;
    PtrParamNode nonLinNode = myParam_->Get("nonLinear", ParamNode::PASS );
    if( nonLinNode ) {
      std::string methodString;
      nonLinNode->GetValue(  "method", methodString, ParamNode::PASS );
      nonLinMethod_ = NonLinMethodTypeEnum.Parse(methodString);
    }
  }


  // *****************************
  //   Definition of Integrators
  // *****************************
  void MagEdgePDE::DefineIntegrators() {

    RegionIdType actRegion;
    BaseMaterial * actMat = NULL;

    // flag for updatedLagrange formulation
    bool upLagrangeForm = true;

    
    // initially, check for regularization factor
    Double regularizationFactor = 1e-6;
    myParam_->GetValue("penaltyFactor", regularizationFactor, ParamNode::PASS);
    
    std::cerr << "regularization factor is " << regularizationFactor << std::endl;
    
    
    //==============================================================
    //begin new implementation
    //==============================================================
    for(UInt iRegion = 0; iRegion < subdoms_.GetSize() ; iRegion ++){
      actRegion = subdoms_[iRegion];
      actMat    = materials_[actRegion];
      std::string regionName = ptgrid_->GetRegion().ToString(actRegion);

      shared_ptr<BaseFeFunction> feFunc = feFunctions_[MAG_POTENTIAL];
      shared_ptr<FeSpace> feSpace = feFunc->GetFeSpace();

      //obtain integration method and order of current subdomain
      IntScheme::IntegMethod curMethod;
      Matrix<Integer> order;
      feSpace->GetIntegration(actRegion,curMethod,order);

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( actRegion );
      // Switch, if region is linear / nonlinear
      if ( regionNonLinType_[actRegion] != NO_NONLINEARITY ) {
       // =================================
       //  Nonlinear Stiffness Integrator
       // =================================
       nLinCurlCurlEdgeInt* curlcurlNL =
           new nLinCurlCurlEdgeInt( actMat, upLagrangeForm );


       curlcurlNL->SetIntegration(ptgrid_->GetIntegrationScheme(),
                                  curMethod,order[0][0]);
       curlcurlNL->SetNonLinMethod( nonLinMethod_ );
       curlcurlNL->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
       curlcurlNL->SetFeSpace( feSpace );
       BiLinFormContext * stiffContext = new BiLinFormContext( curlcurlNL, STIFFNESS );
       stiffContext->SetEntities( actSDList, actSDList );
       stiffContext->SetFeFunctions( feFunc, feFunc );
       assemble_->AddBiLinearForm( stiffContext );

       // we have to save the bilinear form for later usage
       nlinBilinForms_[actRegion] = curlcurlNL;

       // =================================
       //  Nonlinear RHS-integrator
       // =================================
       nLinMagEdge_linFormInt* rhsSource
       = new nLinMagEdge_linFormInt( actMat, upLagrangeForm);
       rhsSource->SetIntegration(ptgrid_->GetIntegrationScheme(),
                                 curMethod,order[0][0]);
       rhsSource->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
       rhsSource->SetFeSpace( feSpace );
       LinearFormContext * rhsContext =
           new LinearFormContext( rhsSource );
       rhsContext->SetEntities( actSDList );
       rhsContext->SetFeFunction( feFunc );
       assemble_->AddLinearForm( rhsContext );

      } else {

       // ***************************************
       // LINEAR PART
       // ***************************************

       // ===============================
       //  Standard Stiffness Integrator
       // ===============================
       CurlCurlEdgeInt* curlcurl =
           new CurlCurlEdgeInt( actMat, upLagrangeForm);
       curlcurl->SetIntegration(ptgrid_->GetIntegrationScheme(),
                                curMethod,order[0][0]);
       curlcurl->SetFeSpace( feSpace );

       BiLinFormContext * stiffContext =
           new BiLinFormContext(curlcurl, STIFFNESS );
       stiffContext->SetEntities( actSDList, actSDList );
       stiffContext->SetFeFunctions( feFunc, feFunc );
       assemble_->AddBiLinearForm( stiffContext );
       linBilinForms_[actRegion] = curlcurl;

       // === Additional RHS integrator in case of Non-linearity ===
       if ( nonLin_ == true ) {
         nLinMagEdge_linFormInt* rhsSource =
             new nLinMagEdge_linFormInt( actMat, upLagrangeForm );
         rhsSource->SetIntegration(ptgrid_->GetIntegrationScheme(),
                                   curMethod,order[0][0]);
         rhsSource->SetFeSpace( feSpace );
         rhsSource->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));

         LinearFormContext * rhsContext =
             new LinearFormContext( rhsSource );
         rhsContext->SetEntities( actSDList );
         rhsContext->SetFeFunction( feFunc );
         assemble_->AddLinearForm( rhsContext );
       }
     } // END OF NOLIN / LIN PART

      // ============================
      // Standard Mass Matrix
      // ============================
      Double conductivity = 0.0; // , maxPerm = 0.0; // TODO: Check if this is still needed
      bool scaleByEdgeSize = false;
      materials_[actRegion]->GetScalar(conductivity,MAG_CONDUCTIVITY,Global::REAL);

      if ( conductivity < 1e-10 || analysistype_ == STATIC ) {
        Double perm;
        // get tensor of permeability and determine max. value
        materials_[actRegion]->GetScalar( perm, MAG_PERMEABILITY, Global::REAL );
        scaleByEdgeSize = true;
        //regularizationFactor = 1e-6;
        conductivity =  regularizationFactor / perm;
      }

      MassEdgeInt *massInt =
          new MassEdgeInt(conductivity, scaleByEdgeSize, upLagrangeForm );
      massInt->SetIntegration(ptgrid_->GetIntegrationScheme(),
                              curMethod,order[0][0]);
      massInt->SetFeSpace( feSpace );

      BiLinFormContext * massContext;
      if ( analysistype_ == STATIC) {
        // we have to guarantee, that we add some mass to
        // the curl-curl-matrix!!
        massContext =  new BiLinFormContext(massInt, STIFFNESS );
      } else {
        massContext =
            new BiLinFormContext(massInt, MASS );
      }
      massContext->SetEntities( actSDList, actSDList );
      massContext->SetFeFunctions( feFunc, feFunc );
      assemble_->AddBiLinearForm( massContext );

      feFunc->AddEntityList( actSDList );
      // ============================
      // COIL INTEGRATORS
      // ============================
      // If this subdomain is a coil we have to do special things
      for ( UInt coil = 0; coil < coilDef_.GetSize(); coil++ ) {
        if ( actRegion == coilRegionId_[coil] ) {
          std::string factor = coilDef_[coil]->value_ + "/" +
            lexical_cast<std::string>(coilDef_[coil]->windingCrossSection_);

          VolForceInt *coilSource3d =
              new LinearEdgeSrcInt ( 3, coilDef_[coil]->phase_,
                                     isaxi_ );

          StdVector<std::string> currDensity(3);
          currDensity[0] = factor + "*" + lexical_cast<std::string>(coilDef_[coil]->locFlowDir_[0]);
          currDensity[1] = factor + "*" + lexical_cast<std::string>(coilDef_[coil]->locFlowDir_[1]);
          currDensity[2] = factor + "*" + lexical_cast<std::string>(coilDef_[coil]->locFlowDir_[2]);
          coilSource3d->SetVolForceVector( currDensity,
                                           coilDef_[coil]->flowCoordSys_,
                                           true, 1.0 );
          LinearFormContext * coilContext =
              new LinearFormContext( coilSource3d );
          coilSource3d->SetFeSpace(feSpace );
          coilSource3d->SetIntegration( ptgrid_->GetIntegrationScheme(),
                                        curMethod,order[0][0] );
          coilContext->SetEntities( actSDList );
          coilContext->SetFeFunction( feFunc );
          assemble_->AddLinearForm( coilContext );
        }
      }

        
    }


        // ============================
        // PERMANENT MAGNETS
        // ============================
        //
//      // check, if this subdomain is a permanent magnet
//      for ( UInt perm = 0; perm < magnetsDomain_.GetSize(); perm++ ) {
//        if ( actRegion == magnetsDomain_[perm] ) {
//          EXCEPTION("Currently magnetic 3D with edge elements do not support permanent magnets");
//
//          Vector<Double> magnetization(dim_);
//          magnetization[0] = magnetsOriX_[perm];
//          magnetization[1] = magnetsOriY_[perm];
//          magnetization[2] = magnetsOriZ_[perm];
//
//          // Get reluctivity for this domain and perform consistency check
//          Double reluctivity;
//          actMat->GetScalar(reluctivity,MAG_RELUCTIVITY,Global::REAL);
//
//          std::string fncname = "none";
//          LinearForm *permSource =
//            new MagPerm3DInt(magnetization, reluctivity,
//                             isaxi_, upLagrangeForm );
//
//          LinearFormContext * permContext =
//            new LinearFormContext( permSource );
//          permContext->SetPtPde( this );
//          permContext->SetResult( results_[0], actSDList );
//          assemble_->AddLinearForm( permContext );
//        }
//      }
//    }
  }

  void MagEdgePDE::DefineSolveStep()
  {
    SolveStepMagEdge *magSolveStep = new SolveStepMagEdge(*this); 
    
   magSolveStep->SetSolStrategy( solStrategy_ );
    solveStep_ = magSolveStep; 
  }


  // ======================================================
  // TIME-STEPPING SECTION
  // ======================================================

  void MagEdgePDE::InitTimeStepping() {
    PtrParamNode systemNode = FindLinearSystem(pdename_);
    TS_alg_ = new Trapezoidal( algsys_, systemNode );
  }


  
  void MagEdgePDE::ReadOlasParams( std::string sysName ) {
    std::string name = pdename_;
    if (solStrategy_ == STRAT_TWO_LEVEL && solStep_ == 2 ) {
      name +="2";
    }
    StdPDE::ReadOlasParams(name);
  }


  // ***************
  //   CalcResults
  // ***************
  void MagEdgePDE::CalcResults ( shared_ptr<BaseResult> res ) {

    switch (res->GetResultInfo()->resultType ) {
      case MAG_RHS_LOAD:
        if( isComplex_ ) {
          ExtractRhsResult<Complex>( res, results_[0] );
        } else {
          ExtractRhsResult<Double>( res, results_[0] );
        }
        break;

      case MAG_FLUX_DENSITY:
        if( isComplex_ ) {
          CalcFluxDensity<Complex>( res );
        } else {
          CalcFluxDensity<Double>( res );
        }
        break;
        
      case MAG_POTENTIAL:
        if( isComplex_ ) {
          CalcVecPotential<Complex>( res );
        } else {
          CalcVecPotential<Double>( res );
        }
        break;
      
//      case MAG_POTENTIAL_DIV:
//        if( isComplex_ ) {
//          CalcMagPotentialDiv<Complex>( res );
//        } else {
//          CalcMagPotentialDiv<Double>( res );
//        }
//        break;
//      
//      case MAG_EDDY_CURRENT:
//         if( isComplex_ ) {
//           CalcEddyCurrent<Complex>( res );
//         } else {
//           CalcEddyCurrent<Double>( res );
//         }
//         break;
//         
      case MAG_ELEM_PERMEABILITY:
        if( isComplex_ ) {
          CalcPermeability<Complex>( res );
        } else {
          CalcPermeability<Double>( res );
        }
        break;
         
      case MAG_ENERGY:
        if( isComplex_ ) {
          CalcEnergy<Complex>( res );
        } else {
          CalcEnergy<Double>( res );
        }
        break;

      default:
        WARN( "Resulttype not computable by magnetic PDE" );
    }

  }
  
  void MagEdgePDE::CalcSpecialResults() {
    
#ifndef USE_INTERPOLATION
    EXCEPTION("Special Results can just be calculated with INTERPOLATION active");
#else
    // fetch fluxlist
    for( UInt iList = 0; iList < calcFlux_.GetSize(); ++iList ) {

      FluxAtPoints & actList = calcFlux_[iList];

      // open file
      std::ofstream  out(actList.fileName.c_str(),std::ios::out );

      out << "# Flux density\n";
      out << "# #elem \t x\t y \t z \t x(Vs/Am) \ty(Vs/Am( \t z(Vs/Am)\txi \teta \tzeta\n";

      // Loop over all points 
      Vector<Double> locCoord, flux;
      for( UInt iPoint = 0; iPoint < actList.points.GetSize(); iPoint++) { 
        const Elem * ptElem = NULL;
        const Vector<Double> & globCoord =actList.points[iPoint].coord;
        if( globCoord.GetSize() == 0) continue;
        ptElem = ptgrid_->GetElemAtGlobalCoord( globCoord, locCoord );
        if( !ptElem ) {
          std::string warnStr = "Cold not find element at position " 
              + actList.points[iPoint].coord.ToString(); 
          Warning( warnStr.c_str());
        } else {
          LocPoint lp(locCoord);
          CalcFluxDensityAtIP(ptElem, lp, flux);

          // write to file
          out << ptElem->elemNum << "\t";
          out << globCoord[0] << "\t" << globCoord[1] << "\t" << globCoord[2] << "\t";
          out << flux[0] << "\t" << flux[1] << "\t" << flux[2] << "\t"
              << locCoord[0] << "\t" << locCoord[1] << "\t" << locCoord[2] << std::endl;
        }

      }
      // close file
      out.close();
    }
#endif
  }



  template<class TYPE>
  void MagEdgePDE::CalcFluxDensity( shared_ptr<BaseResult> result ) {
     Vector<TYPE> elemFlux(dim_);

    // fetch result object and convert data
    Result<TYPE> &  actSol =
      dynamic_cast<Result<TYPE>&>(*result);
    Vector<TYPE> & actVal = actSol.GetVector();
    actVal.Resize( actSol.GetEntityList()->GetSize() * dim_ );

    // loop over elements
    EntityIterator it = actSol.GetEntityList()->GetIterator();
    for ( it.Begin(); !it.IsEnd(); it++ ) {
      LocPoint lp;
      const Elem * el = it.GetElem();
      lp.coord = Elem::shapes[el->type].midPointCoord;
      CalcFluxDensityAtIP( el, lp, elemFlux );
      for( UInt iDof = 0; iDof < dim_; iDof++ ) {
        actVal[it.GetPos()*dim_ + iDof] = elemFlux[iDof];
      }
    }
  }

  template<class TYPE>
  void MagEdgePDE::CalcVecPotential( shared_ptr<BaseResult> result ) {
    Vector<TYPE> elemVecPot(dim_);

    // fetch result object and convert data
    Result<TYPE> &  actSol =
        dynamic_cast<Result<TYPE>&>(*result);
    Vector<TYPE> & actVal = actSol.GetVector();
    actVal.Resize( actSol.GetEntityList()->GetSize() * dim_ );

    // loop over elements
    EntityIterator it = actSol.GetEntityList()->GetIterator();
    for ( it.Begin(); !it.IsEnd(); it++ ) {
      LocPoint lp;
      const Elem * el = it.GetElem();
      lp.coord = Elem::shapes[el->type].midPointCoord;
      CalcVecPotentialAtIP( it, lp, elemVecPot);
      for( UInt iDof = 0; iDof < dim_; iDof++ ) {
        actVal[it.GetPos()*dim_ + iDof] = elemVecPot[iDof];
      }
    }
    }


  template<class TYPE>
  void MagEdgePDE::CalcEddyCurrent( shared_ptr<BaseResult> result ) {
    EXCEPTION("Not yet adapted to new implementation")

//    // fetch result object and convert data
//    Result<TYPE> &  actSol =
//      dynamic_cast<Result<TYPE>&>(*result);
//    EntityIterator it = actSol.GetEntityList()->GetIterator();
//    Vector<TYPE> & actVal = actSol.GetVector();
//    UInt jEddyDofs = 3;
//    actVal.Resize( actSol.GetEntityList()->GetSize() * jEddyDofs );
//
//    Vector<TYPE> jEddyElem;
//    for ( it.Begin(); !it.IsEnd(); it++ ) {
//      CalcEddyCurrentAtIP( it, 0, jEddyElem );
//      for( UInt iDof = 0; iDof < jEddyDofs; iDof++ ) {
//        actVal[it.GetPos()*jEddyDofs + iDof] = jEddyElem[iDof];
//      }
//    }
  }
  
  template<class TYPE>
  void MagEdgePDE::CalcEnergy( shared_ptr<BaseResult> result ) {

    Result<TYPE> &  actSol = 
        dynamic_cast<Result<TYPE>&>(*result);      
    EntityIterator regionIt = actSol.GetEntityList()->GetIterator();

    // resize vector
    Vector<TYPE> & actVal = actSol.GetVector();
    actVal.Resize( actSol.GetEntityList()->GetSize() );

    Vector<TYPE> help;
    TYPE energy;
    Matrix<Double> elemmat;
    
    BaseForm *bilinear_stiff = NULL;
    // Loop over regions
    for( regionIt.Begin(); !regionIt.IsEnd(); regionIt++ ) {
      if ( regionNonLinType_[regionIt.GetRegion()] != NO_NONLINEARITY ) {
        nlinBilinForms_[regionIt.GetRegion()]->SetNonLinMethod(FIXEDPOINT);
        bilinear_stiff = nlinBilinForms_[regionIt.GetRegion()];
      } else {
        bilinear_stiff = linBilinForms_[regionIt.GetRegion()];
      }
      
      ElemList actSDList(ptgrid_ );
      actSDList.SetRegion( regionIt.GetRegion() );
      EntityIterator elemIt = actSDList.GetIterator();
      // Loop over elements
      energy = 0;
      Vector<TYPE> magPot;
      for( elemIt.Begin(); !elemIt.IsEnd(); elemIt++ ) {

        bilinear_stiff->CalcElementMatrix( elemmat, elemIt, elemIt );

        sol_->GetElemSolution(magPot, elemIt );
        help =  elemmat * magPot;
        energy += 0.5 * (help * magPot);

      }  
      actVal[regionIt.GetPos()] = energy;
    }
    
    
    
    
    //    EXCEPTION("Not yet adapted to new implementation")
//    Matrix<Double> elemmat, ptCoord, massMat;
//
//    StdVector<UInt> connecth;
//    StdVector<Integer> Eqns;  
//    Vector<TYPE> help;
//    BaseForm *bilinear_stiff = NULL; 
//
//    // get result
//    Result<TYPE> &  actSol = 
//      dynamic_cast<Result<TYPE>&>(*result);      
//      EntityIterator regionIt = actSol.GetEntityList()->GetIterator();
//
//      // resize vector
//      Vector<TYPE> & actVal = actSol.GetVector();
//      actVal.Resize( actSol.GetEntityList()->GetSize() );
//
//      // Loop over regions
//      for( regionIt.Begin(); !regionIt.IsEnd(); regionIt++ ) {
//
//
//        // Create stiffness integrator
//        if ( regionNonLinType_[regionIt.GetRegion()] != NO_NONLINEARITY ) {
//
////          bilinear_stiff = new nLinCurlCurlEdgeInt( materials_[regionIt.GetRegion()],
////                                                    true );
////
////          bilinear_stiff->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
////
////          // VERY IMPORTANT: Set nonlinear-method "fixpoint", as otherwise also
////          // the Frechet part of the stiffness is calculated!
////          bilinear_stiff->SetNonLinMethod( "fixPoint" );
//        } else {
//          bilinear_stiff = new CurlCurlEdgeInt( materials_[regionIt.GetRegion()], true);
//
//        }
//        ElemList actSDList(ptgrid_ );
//        actSDList.SetRegion( regionIt.GetRegion() );
//        EntityIterator elemIt = actSDList.GetIterator();
//
//        TYPE energy = 0.0;
//        Vector<TYPE> magvecpot;   
//        for ( elemIt.Begin(); !elemIt.IsEnd(); elemIt++ ) {
//          sol_->GetElemSolution(magvecpot, elemIt);
//          bilinear_stiff->CalcElementMatrix(elemmat, elemIt, elemIt);
//          help = elemmat * magvecpot;
//          TYPE temp = 0.5 * (help * magvecpot);
//          energy += temp;
//        }
//        
//        actVal[regionIt.GetPos()] = energy;
//        delete bilinear_stiff;   
//      }
  }


  // ******************************************************
  //   Query parameter object for information about coils
  // ******************************************************
  void MagEdgePDE::ReadCoils() {

    // Check if the element "coils" is present at all.
    // Otherwise leave
    PtrParamNode coilNode = myParam_->Get( "coils", ParamNode::PASS );
    if ( !coilNode )
      return;

    // Get single coil nodes
    ParamNodeList coilNodes = coilNode->GetChildren();

    // Trigger reading in of definitions
    if( coilNodes.GetSize() > 0 ) {
      Info->PrintF( pdename_, "Using the following coils:\n" );
      for( UInt i = 0; i < coilNodes.GetSize(); i++ ) {

        // get region name of actual coil
        std::string regionName = coilNodes[i]->Get("name")->As<std::string>();
        RegionIdType regionId = ptgrid_->GetRegion().Parse( regionName );

        coilRegionId_.Push_back( regionId );
        coilDef_.Push_back( shared_ptr<Coil>( new Coil( regionId,
                                                        coilNodes[i], ptgrid_) ) );
        Info->PrintCoil( *coilDef_.Last(), analysistype_ );
      }
    }
  }


  // ********************************************************
  //   Query parameter object for information about magnets
  // ********************************************************
  void MagEdgePDE::ReadMagnets() {

    // Check if the element "magnets" is present at all.
    // Otherwise leave
    PtrParamNode magnetNode = myParam_->Get( "magnets", ParamNode::PASS );
    if ( !magnetNode )
      return;

    // Get single magnet nodes
    ParamNodeList magnetNodes = magnetNode->GetChildren();

    // trigger definition of magnets
    if( magnetNodes.GetSize() > 0 ) {
      Info->PrintF( pdename_,
              "Found permanent magnets in the following regions:\n" );

      Double tmpDir;
      for( UInt i = 0; i < magnetNodes.GetSize(); i++ ) {

        // get region name of actual magnet
        std::string regionName = magnetNodes[i]->Get("name")->As<std::string>();
        RegionIdType regionId = ptgrid_->GetRegion().Parse( regionName );

        magnetsDomain_.Push_back( regionId );

        // read orientation
        magnetNodes[i]->GetValue( "orientX", tmpDir );
        magnetsOriX_.Push_back( tmpDir );

        magnetNodes[i]->GetValue( "orientY", tmpDir );
        magnetsOriY_.Push_back( tmpDir );

        magnetNodes[i]->GetValue( "orientZ", tmpDir );
        magnetsOriZ_.Push_back( tmpDir );

        // report name to logfile
        Info->PrintF( pdename_, " %s\n", regionName.c_str());
      }
      Info->PrintF( "", "\n" );
    }
  }


  void MagEdgePDE::DefineAvailResults() {

    StdVector<std::string> vecComponents;
    vecComponents = "x", "y", "z";

    // === MAGNETIC VECTOR POTENTIAL ===
    shared_ptr<ResultInfo> res1(new ResultInfo);
//    shared_ptr<AnsatzFct> fct(new NedelecFct);
    res1->resultType = MAG_POTENTIAL;
    res1->dofNames = "";
    res1->unit = "Vs/m";
    res1->definedOn = ResultInfo::EDGE;
    res1->entryType = ResultInfo::SCALAR;
//    res1->fctType = fct;
    
    results_.Push_back( res1 );
    availResults_.insert( res1 );
    
    feFunctions_[MAG_POTENTIAL]->SetResultInfo(res1);
    feFunctions_[MAG_POTENTIAL]->SetPDE(shared_ptr<MagEdgePDE>(this));
    
    // === MAGNETIC FLUX DENSITY ===
    shared_ptr<ResultInfo> flux(new ResultInfo);
    flux->resultType = MAG_FLUX_DENSITY;
    flux->dofNames = vecComponents;
    flux->unit = "Vs/m^2";
    flux->definedOn = ResultInfo::ELEMENT;
    flux->entryType = ResultInfo::VECTOR;
//    flux->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( flux );
    postProcResults_[MAG_FLUX_DENSITY] = MAG_POTENTIAL;

    // === MAGNETIC VECTOR POTENTIAL ===
    shared_ptr<ResultInfo> magPot(new ResultInfo);
    magPot->resultType = MAG_POTENTIAL;
    magPot->dofNames = vecComponents;
    magPot->unit = "A/m^2";
    magPot->definedOn = ResultInfo::ELEMENT;
    magPot->entryType = ResultInfo::VECTOR;
    availResults_.insert( magPot );
    
//    // === EDDY CURRENT DENSITY ===
//    shared_ptr<ResultInfo> eddy(new ResultInfo);
//    eddy->resultType = MAG_EDDY_CURRENT;
//    eddy->dofNames = vecComponents;
//    eddy->unit = "A/m^2";
//    eddy->definedOn = ResultInfo::ELEMENT;
//    eddy->entryType = ResultInfo::VECTOR;
//    eddy->fctType = shared_ptr<ConstFct>(new ConstFct() );
//    availResults_.insert( eddy );
//    
//    // === DIVERGENCE OF MAGNETIC POTENTIAL ===
//    shared_ptr<ResultInfo> div(new ResultInfo);
//    div->resultType = MAG_POTENTIAL_DIV;
//    div->dofNames = "";
//    div->unit = "Vs/m^3";
//    div->definedOn = ResultInfo::ELEMENT;
//    div->entryType = ResultInfo::SCALAR;
//    div->fctType = shared_ptr<ConstFct>(new ConstFct() );
//    availResults_.insert( div );
    
    // === PERMEABILITY  ===
    shared_ptr<ResultInfo> perm(new ResultInfo);
    perm->resultType = MAG_ELEM_PERMEABILITY;
    perm->dofNames = "";
    perm->unit = "Vs/Am";
    perm->definedOn = ResultInfo::ELEMENT;
    perm->entryType = ResultInfo::SCALAR;
    availResults_.insert( perm );

    // === MAGNETIC ENERGY ===
    shared_ptr<ResultInfo> energy(new ResultInfo);
    energy->resultType = MAG_ENERGY;
    energy->dofNames = "";
    energy->unit = "Ws";
    energy->definedOn = ResultInfo::REGION;
    energy->entryType = ResultInfo::SCALAR;
    availResults_.insert( energy );
    postProcResults_[MAG_ENERGY] = MAG_POTENTIAL;
  }

  void MagEdgePDE::ReadSpecialResults() {

    // check, if flux density should be written at special points

    PtrParamNode node; 
    node = myParam_->Get("storeResults")->Get("interpolate",ParamNode::PASS);
    if (!node) return;
    ParamNodeList partList = node->GetList("part");

    // loop over all parts
    for( UInt iPart = 0; iPart < partList.GetSize(); ++iPart ) {
      PtrParamNode  actPartNode = partList[iPart];
      ParamNodeList listNodes = actPartNode->GetList("list");
      calcFlux_.Push_back(FluxAtPoints());
      FluxAtPoints & actFlux = calcFlux_.Last();
      actFlux.fileName = actPartNode->Get("fileName")->As<std::string>();

      // loop over all components
      StdVector<Double> start(3), stop(3), inc(3);
      StdVector<UInt> numSamples(3);

      UInt totalPoints = 1;
      std::string comp;
      UInt compIndex;
      for( UInt iComp = 0; iComp < listNodes.GetSize(); iComp++ ) {
        PtrParamNode actCompNode = listNodes[iComp];
        actCompNode->GetValue("comp", comp);
        compIndex = domain->GetCoordSystem("default")->GetVecComponent(comp)-1;
        actCompNode->GetValue("start", start[compIndex]);
        actCompNode->GetValue("stop", stop[compIndex]);
        actCompNode->GetValue("inc", inc[compIndex]);
        numSamples[compIndex]  = 
            UInt(floor( (stop[compIndex]-start[compIndex]) / inc[compIndex] ) )+1;
        totalPoints *= numSamples[compIndex];
      }
      // create list

      actFlux.points.Resize(totalPoints);
      actFlux.flux.Resize(totalPoints);
      UInt pos = 0;
      for( UInt xSample = 0; xSample < numSamples[0]; xSample++ ) {
        Double actX = start[0] + xSample * inc[0];
        for( UInt ySample = 0; ySample < numSamples[1]; ySample++ ) {
          Double actY = start[1] + ySample * inc[1];
          for( UInt zSample = 0; zSample < numSamples[2]; zSample++ ) {
            Double actZ = start[2] + zSample * inc[2];
            Vector<Double> point(3);
            point[0] = actX;
            point[1] = actY;
            point[2] = actZ;
            LocPoint lp(point);
            actFlux.points[pos++] = lp;
          } // z
        } // y
      } // x
        

      
    } // loop over parts
  }


  // =======================================================================
  //   HELPER METHODS FOR CALCULATING AUXILIARY QUANTITIES
  // =======================================================================

  template<class TYPE>
  void MagEdgePDE::CalcFluxDensityAtIP( const Elem *el,
                                        LocPoint lp,
                                        Vector<TYPE>& field ) {
    Vector<TYPE> tempE, elemSol;
    NodeStoreSol<TYPE> & solhelp = dynamic_cast<NodeStoreSol<TYPE>&>(*sol_);
    field.Resize(dim_);
    LocPointMapped lpm;

    shared_ptr<ElemShapeMap> esm = ptgrid_->GetElemShapeMap( el, true );
    shared_ptr<BaseFeFunction> fct = GetFeFunction(MAG_POTENTIAL);

    CurlCurlEdgeInt * li = NULL;

    if ( regionNonLinType_[el->regionId] != NO_NONLINEARITY ) {
      li = nlinBilinForms_[el->regionId];
    } else {
      li = linBilinForms_[el->regionId];
    }
    BaseFE*  ptFe = fct->GetFeSpace()->GetFe( el->elemNum );
    lpm.Set(lp, esm );
    ElemList elist(ptgrid_);
    elist.SetElement(el);
    solhelp.GetElemSolution( elemSol, elist.GetIterator() );

    li->ApplyBMat( field, lpm, ptFe, elemSol);
  }
  
  template<class TYPE>
  void MagEdgePDE::CalcVecPotentialAtIP( EntityIterator it,
                                         LocPoint lp,
                                         Vector<TYPE>& field ) {
      Vector<TYPE> elemSol;
      NodeStoreSol<TYPE> & solhelp = dynamic_cast<NodeStoreSol<TYPE>&>(*sol_);
      field.Resize(dim_);
      LocPointMapped lpm;
      const Elem * el = it.GetElem();
      shared_ptr<ElemShapeMap> esm = ptgrid_->GetElemShapeMap( el, true );
      shared_ptr<BaseFeFunction> fct = GetFeFunction(MAG_POTENTIAL);
                                                     
      FeHCurl*  ptFe = dynamic_cast<FeHCurl*>(fct->GetFeSpace()->GetFe( it ));
      lpm.Set(Elem::shapes[el->type].midPointCoord, esm );
      solhelp.GetElemSolution( elemSol, it );
      
      Matrix<Double> shape;
      ptFe->GetShFnc(shape, lpm, el);
      field = shape * elemSol;
  }      
//  template<class TYPE>
//   void MagEdgePDE::CalcEddyCurrentAtIP( EntityIterator it,
//                                         UInt ip,
//                                         Vector<TYPE>& jEddy ) {
//
//     Vector<Double> lCoord;
//     Matrix<Double> shpFnc;
//     Vector<TYPE> magVecDeriv1Elem;
//     Double conductivity = 0.0;
//     
//      // Get the right material parameter for current element
//     RegionIdType actRegionId = it.GetElem()->regionId;
//     materials_[actRegionId]
//       ->GetScalar(conductivity,MAG_CONDUCTIVITY,Global::REAL);
//     BaseFE * ptEl = it.GetElem()->ptElem;
//
//     Matrix<Double> cornerCoords;
//     ptgrid_->GetElemNodesCoord( cornerCoords, it.GetElem()->connect, 
//                                 true );
//
//     jEddy.Resize(3);
//     
//     // case 1: dummy integration point
//     if( ip == 0 ) {
//       it.GetElem()->ptElem->GetCoordMidPoint( lCoord );
//       ptEl->CalcEdgeShapeFnc( shpFnc, lCoord, cornerCoords, it.GetElem() );
//     } 
//     else {
//       ptEl->CalcEdgeShapeFncAtIp( shpFnc, ip, cornerCoords, it.GetElem() );
//     }
//     
//     // 1) get part coming from vector potential
//     GetDerivSolVecOfElement(magVecDeriv1Elem,it,results_[0]);
//     jEddy.Init();
//     for( UInt iDof = 0; iDof < 3; iDof++ ) {
//       for( UInt i = 0; i < shpFnc.GetNumRows(); i++ ) {
//         jEddy[iDof] += shpFnc[i][iDof] * magVecDeriv1Elem[i];
//       }
//     }
//     jEddy *= -conductivity;   
//  }
//   

  template<class TYPE>
  void MagEdgePDE::CalcPermeability( shared_ptr<BaseResult> result ) {

    TYPE elemPerm;
    Vector<TYPE> elemFlux;
    Double bAbs, reluct;

    // fetch result object and convert data
    Result<TYPE> &  actSol = 
        dynamic_cast<Result<TYPE>&>(*result);
    Vector<TYPE> & actVal = actSol.GetVector();
    actVal.Resize( actSol.GetEntityList()->GetSize());

    // loop over elements
    EntityIterator it = actSol.GetEntityList()->GetIterator();
    for ( it.Begin(); !it.IsEnd(); it++ ) {

      // Determine regionId of element
      const Elem & actEl = *(it.GetElem());
      RegionIdType actRegion = actEl.regionId; 

      // Check, if region is nonlinear
      if ( regionNonLinType_[actRegion] == PERMEABILITY ) {

        // Obtain nonlinear approximation functional
        ApproxData * approx  = materials_[actRegion]->GetNonlinFncBH(MAG_PERMEABILITY);

        // Calculate flux density in element midpoint
        LocPoint lp;
        lp.coord = Elem::shapes[actEl.type].midPointCoord;
        CalcFluxDensityAtIP( it.GetElem(), lp, elemFlux );
        bAbs = elemFlux.NormL2();
        reluct = approx->EvaluateFuncNu(bAbs);
        elemPerm = 1.0 / reluct;
      } else {
        materials_[actRegion]->GetScalar( elemPerm, MAG_PERMEABILITY,Global::REAL); 
      }
      actVal[it.GetPos() ] = elemPerm;
    }
  }
//  template<class TYPE>
//  void MagEdgePDE::CalcMagPotentialDiv( shared_ptr<BaseResult> result ) {
//    // fetch result object and convert data
//    Result<TYPE> &  actSol = 
//      dynamic_cast<Result<TYPE>&>(*result);
//    EntityIterator it = actSol.GetEntityList()->GetIterator();
//    Vector<TYPE> & actVal = actSol.GetVector();
//    actVal.Resize( actSol.GetEntityList()->GetSize() );
//
//    StdVector<Matrix<Double> >  globDeriv;
//    Matrix<Double> cornerCoords;
//    Vector<Double> elemMidPoint;
//    Vector<TYPE> elemSol;
//    TYPE elemDiv;
//
//    for( it.Begin(); !it.IsEnd(); it++ ) {
//      const Elem * ptElem = it.GetElem();
//
//      // get coordinates of element
//      ptgrid_->GetElemNodesCoord( cornerCoords, ptElem->connect, true );
//
//      // fetch global gradient
//      ptElem->ptElem->GetCoordMidPoint( elemMidPoint );
//      ptElem->ptElem->GetEdgeGlobalDerivShapeFnc( globDeriv, elemMidPoint, 
//                                                  cornerCoords, ptElem );
//
//      // getch solution of element
//      sol_->GetElemSolution(elemSol, it);
//
//      // loop over shape functions
//      elemDiv = 0.0;
//      for( UInt i = 0; i < elemSol.GetSize(); i++ )  {
//
//        elemDiv += elemSol[i] * ( globDeriv[i][0][0] 
//                                + globDeriv[i][1][1] 
//                                + globDeriv[i][2][2] );
//      }
//      actVal[it.GetPos()] = elemDiv;
//    }
//  }

  
  // ======================================================
  // COUPLING SECTION
  // ======================================================


  void MagEdgePDE::InitCoupling(PDECoupling * Coupling) {

    isIterCoupled_ = true;
    ptCoupling_   = Coupling;

    // Enable update of geometry
    const UInt numCouplings = ptCoupling_->GetNumOutputCouplings();


    cplNodeNumPos_.Resize( numCouplings );

    for ( UInt actCoupling = 0; actCoupling < numCouplings; actCoupling++ ){

      if (ptCoupling_->GetOutputQuantity(actCoupling) == MAG_FORCE_LORENTZ) {

        // Intialize the memory of the coupling values
        ptCoupling_->CreateCouplingVector(actCoupling,isComplex_);

        //get the element-node to coupling node matching
        StdVector<std::string> couplRegions;
        StdVector<RegionIdType> regionIds;
        ptCoupling_->GetOutputRegions(actCoupling, couplRegions);
        ptgrid_->GetRegion().Parse( couplRegions , regionIds );

        // Check, that every coupling region is part of
        // the magnetic pde itself
        for( UInt iReg = 0; iReg < couplRegions.GetSize(); iReg++ ) {
          if( subdoms_.Find(regionIds[iReg]) == -1 ) {
            EXCEPTION( "Coupling region '" << couplRegions[iReg]
                       << "' is not contained in regions for"
                       << " magnetic PDE" );
          }
        }

        StdVector<UInt> * couplingnodes = NULL;
        ptCoupling_->GetOutputNodes(actCoupling, couplingnodes);

        if ( couplingnodes == NULL ) {
          EXCEPTION( "Pointer couplingnodes = NULL!" );
        }

        for( UInt iNode = 0; iNode < couplingnodes->GetSize(); iNode++ ) {
          UInt actNode = (*couplingnodes)[iNode];
          cplNodeNumPos_[actCoupling][actNode] = iNode;
        }
      }
    }
  }

  void MagEdgePDE::CalcOutputCoupling() {

    SolutionType quantity;
    StdVector<UInt> * couplingNodes = NULL;
    SingleVector * values = NULL;
    UInt forcesCount = 0;

    // at first, check if this PDE is iterative coupled
    if (isIterCoupled_ == false)
      return;

    // loop over all output coupling quantities
    for ( UInt actCoupling = 0;
          actCoupling < ptCoupling_->GetNumOutputCouplings();
          actCoupling++ ) {

      quantity = ptCoupling_->GetOutputQuantity(actCoupling);
      ptCoupling_->GetOutputValues(actCoupling, values);
      Vector<Double> *temp = dynamic_cast<Vector<Double> *>(values);

      switch(ptCoupling_->GetOutputType(actCoupling)) {

      case NODE:
        ptCoupling_->GetOutputNodes(actCoupling, couplingNodes);

        if (quantity == MAG_FORCE_LORENTZ) {
          CalcNodeForceLorentz( *temp, cplNodeNumPos_[forcesCount],
                                actCoupling, couplingNodes->GetSize());
          forcesCount++;
        }

        break;

      case ELEM:
        EXCEPTION( "No Element coupling output" );
      }
    }
  }

  void MagEdgePDE::
  CalcNodeForceLorentz( Vector<Double> & force, 
                        std::map<UInt, UInt> & cplNodeNumPos,
                        UInt actCoupling, UInt numCouplingNodes) {
    EXCEPTION("Not adjusted to new implementation");
    
//    //get the coupling regions
//    StdVector<std::string> couplRegions;
//    StdVector<RegionIdType> regionIds;
//    ptCoupling_->GetOutputRegions(actCoupling, couplRegions);
//    ptgrid_->GetRegion().Parse(couplRegions, regionIds);
//
//   
//    force.Init();
//    Vector<Double>  fluxAtIp(dim_), elemForce, fAtIp, jAtIp;
//    Matrix<Double> ptCoord;
//
//    for (UInt reg=0; reg<couplRegions.GetSize(); reg++) {
//
//      //find subdomain index
//      Integer sdIndex = subdoms_.Find( regionIds[reg] );
//      if( sdIndex == -1 ) {
//        EXCEPTION( "The region coupling region '" <<
//                   ptgrid_->GetRegion().ToString( regionIds[reg] )
//                   << "' was not found in magneticPDE" );
//      }
//      
//      RegionIdType actRegionId = subdoms_[sdIndex] ;
//      Double conductivity;
//      materials_[actRegionId]->GetScalar(conductivity,MAG_CONDUCTIVITY,Global::REAL);
//      
//      // Check if this region is a coil
//      Integer coilIndex = coilRegionId_.Find(actRegionId);
//      
//      ElemList actSDList(ptgrid_ );
//      actSDList.SetRegion(actRegionId );
//      
//       EntityIterator it = actSDList.GetIterator();
//       UInt actEl = 0;
//       // iterate over all elements of regions
//       for ( it.Begin(); !it.IsEnd(); it++, actEl++ ) {
//         
//         BaseFE * ptElem = it.GetElem()->ptElem;
//
//         //need info for ansatz functions of mechanical FE!
//         shared_ptr<AnsatzFct> fct(new LagrangeFct);
//         ptElem->SetAnsatzFct( fct );
//         UInt numFncs = ptElem->GetNumFncs(fct);
////          ptElem->SetAnsatzFct( results_[0]->fctType );
////          UInt numFncs = ptElem->GetNumFncs(results_[0]->fctType);
//
//         const UInt nrIntPts = ptElem->GetNumIntPoints();
//         const Vector<Double> & intWeights = ptElem->GetIntWeights();  
//         
//         Vector<Double> shpFncAtIp;         
//         ptgrid_->GetElemNodesCoord( ptCoord, 
//                                     it.GetElem()->connect,
//                                     true );
//         
//         // iterate over all integration points
//         fAtIp.Resize( dim_ * numFncs );
//         elemForce.Resize( dim_ * numFncs );
//         elemForce.Init();
//         for( UInt ip = 1; ip <= nrIntPts; ip++ ) {
//           ptElem->GetShFncAtIp(shpFncAtIp, ip, it.GetElem() );
//
//           // CHECK: If this region is a current coil, we simply take the
//           // prescribed current density value
//           if( coilIndex != -1 ) {
//             MathParser * mParser =  domain->GetMathParser();
//             std::string factor = coilDef_[coilIndex]->value_ + "/" 
//               + lexical_cast<std::string>(coilDef_[coilIndex]->windingCrossSection_ );
//             mParser->SetExpr( mHandle_, factor );
//             // TODO: Check if this is still needed
//             // Double currDens = mParser->Eval(mHandle_);
//           } 
//           else {
//             CalcEddyCurrentAtIP( it, 0, jAtIp );
//           }
//           
//           CalcFluxDensityAtIP( it, 0, fluxAtIp );
//           
//           // calculate cross product 
//           fAtIp.Init();
//           Vector<Double> tempCross;
//           jAtIp.CrossProduct( fluxAtIp, tempCross );
//           for (UInt iFnc=0; iFnc<numFncs; iFnc++) {
//             fAtIp[iFnc*3+0] = -tempCross[0] * shpFncAtIp[iFnc];
//             fAtIp[iFnc*3+1] = -tempCross[1] * shpFncAtIp[iFnc];
//             fAtIp[iFnc*3+2] = -tempCross[2] * shpFncAtIp[iFnc];
//           }
//
//           Double jacDet = ptElem->CalcJacobianDetAtIp(ip, ptCoord, 
//                                                       it.GetElem() );
//           elemForce += -fAtIp * (jacDet * intWeights[ip-1] );
//
//         }
//         StdVector<UInt> const & connecth = it.GetElem()->connect;
//
//         // Add the element force to the according coupling node
//         for (UInt ielemnode=0; ielemnode<connecth.GetSize(); ielemnode++) {
//           UInt pos = cplNodeNumPos[connecth[ielemnode]];
//           for( UInt idim=0; idim<dim_; idim++) {
//             force[pos*dim_+idim] += elemForce[ielemnode*dim_+idim];
//           }
//         }
//       }
//       
// 
//     }
  }
  
  std::map<SolutionType, shared_ptr<FeSpace> > MagEdgePDE::CreateFeSpaces(std::string formulation){
    //ok default case so we create grid based approximation H1 elements
    //and standard Gauss integration
    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    if(formulation == "default" || formulation == "H_CURL"){
      crSpaces[MAG_POTENTIAL] = FeSpace::CreateInstance(myParam_,FeSpace::HCURL);
      crSpaces[MAG_POTENTIAL]->Init();
    }else{
      EXCEPTION("The formulation " << formulation << "of magnetig edge PDE is not known!");
    }
    return crSpaces;
  }

  bool MagEdgePDE::HasOutput( SolutionType output ) {
    if (output == MAG_FORCE_LORENTZ) {
      return true;
    }
    return false;
  }

} // end of namespace

