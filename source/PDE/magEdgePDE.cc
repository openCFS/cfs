// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>

#include "magEdgePDE.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Driver/stdSolveStep.hh"
#include "Utils/Coil.hh"
#include "Utils/SmoothSpline.hh"
#include "Utils/LinInterpolate.hh"
#include "Forms/curlfieldop.hh"
#include "Forms/curlCurlEdgeInt.hh"
#include "Forms/massEdgeInt.hh"
#include "Forms/linearForm.hh"
#include "trapezoidal.hh"
#include "CoupledPDE/pdecoupling.hh"
#include "Domain/ansatzFct.hh"
#include "Driver/assemble.hh"
#include "Utils/coordSystem.hh"

#ifdef USE_SCRIPTING
#include "DataInOut/Scripting/cfsmessenger.hh"
#endif

namespace CoupledField {

  // **************
  //  Constructor
  // **************
  MagEdgePDE::MagEdgePDE( Grid * aptgrid, ParamNode* paramNode )
    :SinglePDE( aptgrid, paramNode ) {

    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_          = "magneticEdge";
    pdematerialclass_ = ELECTROMAGNETIC;
    maxTimeDerivOrder_ = 1;

    // check if we have a 3d setup
    bool is3d = param->Get("domain")->Get("geometryType")->AsString() == "3d";
    if ( !is3d )
      EXCEPTION("MagEdgePDE is just implemented for 3D setups!");
  }


  // *************
  //  Destructor
  // *************
  MagEdgePDE::~MagEdgePDE() {
  }


  // *********************************
  //  Read special boundary conitions
  // *********************************

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
    ParamNode * nonLinListNode = myParam_->Get("nonLinList", false );
    if( !nonLinListNode)
      return;

    // Get nonlinear types
    StdVector<ParamNode*> nonLinNodes = nonLinListNode->GetChildren();
    for( UInt i = 0; i < nonLinNodes.GetSize(); i++ ) {

      std::string actTypeString = nonLinNodes[i]->GetName();
      std::string actId = nonLinNodes[i]->Get("id")->AsString();

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
    StdVector<ParamNode*> regionNodes =
      myParam_->Get("regionList")->GetChildren();

    RegionIdType actRegionId;
    std::string actRegionName, actNonLinId;

    if( regionNodes.GetSize() > 0 ) {
      Info->PrintF( pdename_, "Non-linearity in following region(s)\n" );
    }
    for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {

      // get data
      regionNodes[i]->Get( "name", actRegionName );
      regionNodes[i]->Get( "nonLinId", actNonLinId );

      if( actNonLinId == "" )
        continue;

      actRegionId = ptgrid_->RegionNameToId( actRegionName );

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
    ParamNode * nonLinNode = myParam_->Get("nonLinear", false );
    nonLinMethod_ = "fixPoint";
    if( nonLinNode ) {
      nonLinNode->Get(  "method", nonLinMethod_, false );
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

   // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {

      // Set current region and material
       actRegion = it->first;
       actMat    = it->second;

      // Get current region node
      std::string regionName = ptgrid_->RegionIdToName( actRegion );

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( actRegion );

      if ( regionNonLinType_[actRegion] != NO_NONLINEARITY ) {

        Warning( "Nonlinearity for 3D magnetic edge formulation not yet ported from NACS" );
//        if ( regionNonLinType_[actRegion] == HYSTERESIS ) {
//          EXCEPTION("Magnetics with nonlineaity in 3D not supported");
//        }
//
//        nLinCurlCurlNode3DInt* curlcurlNL = 
//          new nLinCurlCurlNode3DInt( actMat, upLagrangeForm );
//
//        curlcurlNL->SetNonLinMethod( nonLinMethod_ );      
//        curlcurlNL->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
//        
//        BiLinFormContext * stiffContext = 
//          new BiLinFormContext( curlcurlNL, STIFFNESS );
//        stiffContext->SetPtPdes(this, this);   
//        stiffContext->SetResults( results_[0], results_[0],
//                                  actSDList, actSDList );     
//        assemble_->AddBiLinearForm( stiffContext);
//        
//        //save bilinearForm
//        pdeBilinearForms_[actRegion][curlcurlNL->GetName()] = curlcurlNL;
//        
//        if ( regionNonLinType_[actRegion] == PERMEABILITY ) {
//          // nonlinear RHS linearform!!
//          nLinMagNode3D_linFormInt* rhsSource 
//            = new nLinMagNode3D_linFormInt( actMat, upLagrangeForm);
//
//          rhsSource->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
//          LinearFormContext * rhsContext = 
//            new LinearFormContext( rhsSource );
//          rhsContext->SetPtPde( this );
//          rhsContext->SetResult( results_[0], actSDList );
//          assemble_->AddLinearForm( rhsContext );
//        }
      }
      else {
        CurlCurlEdgeInt* curlcurl =
          new CurlCurlEdgeInt( actMat, upLagrangeForm);

        BiLinFormContext * stiffContext =
          new BiLinFormContext(curlcurl, STIFFNESS );
        stiffContext->SetPtPdes(this, this);  
        stiffContext->SetResults( results_[0], results_[0],
                                  actSDList, actSDList );
        assemble_->AddBiLinearForm( stiffContext );

        //save bilinearForm
        pdeBilinearForms_[actRegion][curlcurl->GetName()] = curlcurl;


        if ( nonLin_ == true ) {
          // for nonlinear RHS linearform we need linear and nonlinear
          // subdomains
          nLinMagNode3D_linFormInt* rhsSource =
            new nLinMagNode3D_linFormInt( actMat, upLagrangeForm );

          rhsSource->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
          LinearFormContext * rhsContext =
            new LinearFormContext( rhsSource );
          rhsContext->SetPtPde( this );
          rhsContext->SetResult( results_[0], actSDList );
          assemble_->AddLinearForm( rhsContext );
        }
      }

      // Mass matrix
      Double conductivity = 0.0; // , maxPerm = 0.0; // TODO: Check if this is still needed
       bool scaleByEdgeSize = false;
       materials_[actRegion]->GetScalar(conductivity,MAG_CONDUCTIVITY,Global::REAL);
       
       if ( conductivity < 1e-10 || analysistype_ == STATIC ) {
         Double regularizationFactor, perm;
         // get tensor of permeability and determine max. value
           materials_[actRegion]->GetScalar( perm, MAG_PERMEABILITY, Global::REAL );
         scaleByEdgeSize = true;
         regularizationFactor = 1e-6;
         conductivity =  regularizationFactor / perm;
       } 
       
       MassEdgeInt *bilinear_mass = 
         new MassEdgeInt(conductivity, scaleByEdgeSize, upLagrangeForm );

       BiLinFormContext * massContext;
       if ( analysistype_ == STATIC) {
         // we have to guarantee, that we add some mass to 
         // the curl-curl-matrix!!
         massContext = 
           new BiLinFormContext(bilinear_mass, STIFFNESS );
       }
       else {
         massContext = 
           new BiLinFormContext(bilinear_mass, MASS );
       }

       massContext->SetPtPdes(this, this);
       massContext->SetResults( results_[0], results_[0],
                                actSDList, actSDList );        
       assemble_->AddBiLinearForm( massContext );
       
       //save bilinearForm
       pdeBilinearForms_[actRegion][bilinear_mass->GetName()] = bilinear_mass;

      // If this subdomain is a coil we have to do special things
      for ( UInt coil = 0; coil < coilDef_.GetSize(); coil++ ) {
        if ( actRegion == coilRegionId_[coil] ) {
          std::string factor = coilDef_[coil]->value_ + "/" +
            GenStr(coilDef_[coil]->windingCrossSection_);

          VolForceInt *coilSource3d =
            new LinearEdgeSrcInt ( 3, coilDef_[coil]->phase_,
                                   isaxi_ );

          StdVector<std::string> currDensity(3);
          currDensity[0] = factor + "*" + GenStr(coilDef_[coil]->locFlowDir_[0]);
          currDensity[1] = factor + "*" + GenStr(coilDef_[coil]->locFlowDir_[1]);
          currDensity[2] = factor + "*" + GenStr(coilDef_[coil]->locFlowDir_[2]);
          coilSource3d->SetVolForceVector( currDensity,
                                           coilDef_[coil]->flowCoordSys_,
                                           true, 1.0 );
          LinearFormContext * coilContext =
            new LinearFormContext( coilSource3d );
          coilContext->SetPtPde( this );
          coilContext->SetResult( results_[0], actSDList );
          assemble_->AddLinearForm( coilContext );
        }
      }

      // check, if this subdomain is a permanent magnet
      for ( UInt perm = 0; perm < magnetsDomain_.GetSize(); perm++ ) {
        if ( actRegion == magnetsDomain_[perm] ) {
          EXCEPTION("Currently magnetic 3D with edge elements do not support permanent magnets");

          Vector<Double> magnetization(dim_);
          magnetization[0] = magnetsOriX_[perm];
          magnetization[1] = magnetsOriY_[perm];
          magnetization[2] = magnetsOriZ_[perm];

          // Get reluctivity for this domain and perform consistency check
          Double reluctivity;
          actMat->GetScalar(reluctivity,MAG_RELUCTIVITY,Global::REAL);

          std::string fncname = "none";
          LinearForm *permSource =
            new MagPerm3DInt(magnetization, reluctivity,
                             isaxi_, upLagrangeForm );

          LinearFormContext * permContext =
            new LinearFormContext( permSource );
          permContext->SetPtPde( this );
          permContext->SetResult( results_[0], actSDList );
          assemble_->AddLinearForm( permContext );
        }
      }

      // Give result to equation numbering class
      eqnMap_->AddResult( *results_[0], actSDList );

    }
  }

  void MagEdgePDE::DefineSolveStep()
  {
    solveStep_ = new StdSolveStep(*this);
  }


  // ======================================================
  // TIME-STEPPING SECTION
  // ======================================================

  void MagEdgePDE::InitTimeStepping() {
    TS_alg_ = new Trapezoidal( algsys_ );
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
      
      case MAG_POTENTIAL_DIV:
        if( isComplex_ ) {
          CalcMagPotentialDiv<Complex>( res );
        } else {
          CalcMagPotentialDiv<Double>( res );
        }
        break;
      
      case MAG_EDDY_CURRENT:
         if( isComplex_ ) {
           CalcEddyCurrent<Complex>( res );
         } else {
           CalcEddyCurrent<Double>( res );
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
        Warning( "Resulttype not computable by magnetic PDE",
                 __FILE__, __LINE__ );
    }

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
      CalcFluxDensityAtIP( it, 0, elemFlux );
      for( UInt iDof = 0; iDof < dim_; iDof++ ) {
        actVal[it.GetPos()*dim_ + iDof] = elemFlux[iDof];
      }
    }
  }


  template<class TYPE>
  void MagEdgePDE::CalcEddyCurrent( shared_ptr<BaseResult> result ) {


    // fetch result object and convert data
    Result<TYPE> &  actSol =
      dynamic_cast<Result<TYPE>&>(*result);
    EntityIterator it = actSol.GetEntityList()->GetIterator();
    Vector<TYPE> & actVal = actSol.GetVector();
    UInt jEddyDofs = 3;
    actVal.Resize( actSol.GetEntityList()->GetSize() * jEddyDofs );

    Vector<TYPE> jEddyElem;
    for ( it.Begin(); !it.IsEnd(); it++ ) {
      CalcEddyCurrentAtIP( it, 0, jEddyElem );
      for( UInt iDof = 0; iDof < jEddyDofs; iDof++ ) {
        actVal[it.GetPos()*jEddyDofs + iDof] = jEddyElem[iDof];
      }
    }
  }
  
  template<class TYPE>
  void MagEdgePDE::CalcEnergy( shared_ptr<BaseResult> result ) {

    Matrix<Double> elemmat, ptCoord, massMat;

    StdVector<UInt> connecth;
    StdVector<Integer> Eqns;  
    Vector<TYPE> help;
    BaseForm *bilinear_stiff = NULL; 

    // get result
    Result<TYPE> &  actSol = 
      dynamic_cast<Result<TYPE>&>(*result);      
      EntityIterator regionIt = actSol.GetEntityList()->GetIterator();

      // resize vector
      Vector<TYPE> & actVal = actSol.GetVector();
      actVal.Resize( actSol.GetEntityList()->GetSize() );

      // Loop over regions
      for( regionIt.Begin(); !regionIt.IsEnd(); regionIt++ ) {


        // Create stiffness integrator
        if ( regionNonLinType_[regionIt.GetRegion()] != NO_NONLINEARITY ) {

//          bilinear_stiff = new nLinCurlCurlEdgeInt( materials_[regionIt.GetRegion()],
//                                                    true );
//
//          bilinear_stiff->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
//
//          // VERY IMPORTANT: Set nonlinear-method "fixpoint", as otherwise also
//          // the Frechet part of the stiffness is calculated!
//          bilinear_stiff->SetNonLinMethod( "fixPoint" );
        } else {
          bilinear_stiff = new CurlCurlEdgeInt( materials_[regionIt.GetRegion()], true);

        }
        ElemList actSDList(ptgrid_ );
        actSDList.SetRegion( regionIt.GetRegion() );
        EntityIterator elemIt = actSDList.GetIterator();

        TYPE energy = 0.0;
        Vector<TYPE> magvecpot;   
        for ( elemIt.Begin(); !elemIt.IsEnd(); elemIt++ ) {
          sol_->GetElemSolution(magvecpot, elemIt);
          bilinear_stiff->CalcElementMatrix(elemmat, elemIt, elemIt);
          help = elemmat * magvecpot;
          TYPE temp = 0.5 * (help * magvecpot);
          energy += temp;
        }
        
        actVal[regionIt.GetPos()] = energy;
        delete bilinear_stiff;   
      }
  }


  // ******************************************************
  //   Query parameter object for information about coils
  // ******************************************************
  void MagEdgePDE::ReadCoils() {

    // Check if the element "coils" is present at all.
    // Otherwise leave
    ParamNode * coilNode = myParam_->Get( "coils", false );
    if ( !coilNode )
      return;

    // Get single coil nodes
    StdVector<ParamNode*> coilNodes = coilNode->GetChildren();

    // Trigger reading in of definitions
    if( coilNodes.GetSize() > 0 ) {
      Info->PrintF( pdename_, "Using the following coils:\n" );
      for( UInt i = 0; i < coilNodes.GetSize(); i++ ) {

        // get region name of actual coil
        std::string regionName = coilNodes[i]->Get("name")->AsString();
        RegionIdType regionId = ptgrid_->RegionNameToId( regionName );

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
    ParamNode * magnetNode = myParam_->Get( "magnets", false );
    if ( !magnetNode )
      return;

    // Get single magnet nodes
    StdVector<ParamNode*> magnetNodes = magnetNode->GetChildren();

    // trigger definition of magnets
    if( magnetNodes.GetSize() > 0 ) {
      Info->PrintF( pdename_,
              "Found permanent magnets in the following regions:\n" );

      Double tmpDir;
      for( UInt i = 0; i < magnetNodes.GetSize(); i++ ) {

        // get region name of actual magnet
        std::string regionName = magnetNodes[i]->Get("name")->AsString();
        RegionIdType regionId = ptgrid_->RegionNameToId( regionName );

        magnetsDomain_.Push_back( regionId );

        // read orientation
        magnetNodes[i]->Get( "orientX", tmpDir );
        magnetsOriX_.Push_back( tmpDir );

        magnetNodes[i]->Get( "orientY", tmpDir );
        magnetsOriY_.Push_back( tmpDir );

        magnetNodes[i]->Get( "orientZ", tmpDir );
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
    shared_ptr<AnsatzFct> fct(new NedelecFct);
    res1->resultType = MAG_POTENTIAL;
    res1->dofNames = "";
    res1->unit = "Vs/m";
    res1->definedOn = ResultInfo::EDGE;
    res1->entryType = ResultInfo::SCALAR;
    res1->fctType = fct;
    
    results_.Push_back( res1 );
    availResults_.insert( res1 );
    
    // === MAGNETIC FLUX DENSITY ===
    shared_ptr<ResultInfo> flux(new ResultInfo);
    flux->resultType = MAG_FLUX_DENSITY;
    flux->dofNames = vecComponents;
    flux->unit = "Vs/m^2";
    flux->definedOn = ResultInfo::ELEMENT;
    flux->entryType = ResultInfo::VECTOR;
    flux->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( flux );

    // === EDDY CURRENT DENSITY ===
    shared_ptr<ResultInfo> eddy(new ResultInfo);
    eddy->resultType = MAG_EDDY_CURRENT;
    eddy->dofNames = vecComponents;
    eddy->unit = "A/m^2";
    eddy->definedOn = ResultInfo::ELEMENT;
    eddy->entryType = ResultInfo::VECTOR;
    eddy->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( eddy );
    
    // === DIVERGENCE OF MAGNETIC POTENTIAL ===
    shared_ptr<ResultInfo> div(new ResultInfo);
    div->resultType = MAG_POTENTIAL_DIV;
    div->dofNames = "";
    div->unit = "Vs/m^3";
    div->definedOn = ResultInfo::ELEMENT;
    div->entryType = ResultInfo::SCALAR;
    div->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( div );
    
    // === MAGNETIC ENERGY ===
     shared_ptr<ResultInfo> energy(new ResultInfo);
     energy->resultType = MAG_ENERGY;
     energy->dofNames = "";
     energy->unit = "Ws";
     energy->definedOn = ResultInfo::REGION;
     energy->entryType = ResultInfo::SCALAR;
     energy->fctType = shared_ptr<ConstFct>(new ConstFct() );
     availResults_.insert( energy );
  }

  void MagEdgePDE::ReadSpecialResults() {

  }


  // =======================================================================
  //   HELPER METHODS FOR CALCULATING AUXILIARY QUANTITIES
  // =======================================================================

  template<class TYPE>
  void MagEdgePDE::CalcFluxDensityAtIP( EntityIterator it,
                                    UInt ip,
                                    Vector<TYPE>& field ) {

    // get element solution
    Vector<TYPE> elSol;
    sol_->GetElemSolution(elSol, it);

    field.Resize(3);
    RegionIdType actRegionId = it.GetElem()->regionId;
    
    CurlCurlEdgeInt* curlOp = NULL;
    if ( regionNonLinType_[actRegionId] != NO_NONLINEARITY ) {
//      std::string bilinearName = "nLinCurlCurlEdgeInt";
//      curlOp = dynamic_cast<nLinCurlCurlEdgeInt*>(pdeBilinearForms_[actRegionId][bilinearName]);
    }
    else {
      std::string bilinearName = "CurlCurlEdgeInt";
      curlOp = dynamic_cast<CurlCurlEdgeInt*>(pdeBilinearForms_[actRegionId][bilinearName]);
    }

    //set element info
    curlOp->ExtractElemInfo( it );
    Matrix<Double> CornerCoords, bMatCurl;
    ptgrid_->GetElemNodesCoord( CornerCoords, it.GetElem()->connect,
                                true );

    StdVector< Matrix<Double> > xiDx;
    const UInt nrEdges = it.GetElem()->ptElem->GetNumEdges();
    xiDx.Resize(nrEdges);

    // case 1: element midpoint
    if( ip == 0 ) {
      Vector<Double> intPoint;
      it.GetElem()->ptElem->GetCoordMidPoint(intPoint);
      it.GetElem()->ptElem->GetEdgeGlobalDerivShapeFnc(xiDx, intPoint,  
                                                       CornerCoords, it.GetElem() );
      curlOp->CalcEdgeCurl( bMatCurl, xiDx );
    } else {
      // case2: real integration point
      it.GetElem()->ptElem->GetEdgeGlobDerivShFncAtIp(xiDx, ip,  CornerCoords, 
                                                      it.GetElem() );
      curlOp->CalcEdgeCurl( bMatCurl, xiDx );
    }
    field = bMatCurl * elSol;   
  }

  template<class TYPE>
   void MagEdgePDE::CalcEddyCurrentAtIP( EntityIterator it,
                                         UInt ip,
                                         Vector<TYPE>& jEddy ) {

     Vector<Double> lCoord;
     Matrix<Double> shpFnc;
     Vector<TYPE> magVecDeriv1Elem;
     Double conductivity = 0.0;
     
      // Get the right material parameter for current element
     RegionIdType actRegionId = it.GetElem()->regionId;
     materials_[actRegionId]
       ->GetScalar(conductivity,MAG_CONDUCTIVITY,Global::REAL);
     BaseFE * ptEl = it.GetElem()->ptElem;

     Matrix<Double> cornerCoords;
     ptgrid_->GetElemNodesCoord( cornerCoords, it.GetElem()->connect, 
                                 true );

     jEddy.Resize(3);
     
     // case 1: dummy integration point
     if( ip == 0 ) {
       it.GetElem()->ptElem->GetCoordMidPoint( lCoord );
       ptEl->CalcEdgeShapeFnc( shpFnc, lCoord, cornerCoords, it.GetElem() );
     } 
     else {
       ptEl->CalcEdgeShapeFncAtIp( shpFnc, ip, cornerCoords, it.GetElem() );
     }
     
     // 1) get part coming from vector potential
     GetDerivSolVecOfElement(magVecDeriv1Elem,it,results_[0]);
     jEddy.Init();
     for( UInt iDof = 0; iDof < 3; iDof++ ) {
       for( UInt i = 0; i < shpFnc.GetNumRows(); i++ ) {
         jEddy[iDof] += shpFnc[i][iDof] * magVecDeriv1Elem[i];
       }
     }
     jEddy *= -conductivity;   
  }
   
  template<class TYPE>
  void MagEdgePDE::CalcMagPotentialDiv( shared_ptr<BaseResult> result ) {
    // fetch result object and convert data
    Result<TYPE> &  actSol = 
      dynamic_cast<Result<TYPE>&>(*result);
    EntityIterator it = actSol.GetEntityList()->GetIterator();
    Vector<TYPE> & actVal = actSol.GetVector();
    actVal.Resize( actSol.GetEntityList()->GetSize() );

    StdVector<Matrix<Double> >  globDeriv;
    Matrix<Double> cornerCoords;
    Vector<Double> elemMidPoint;
    Vector<TYPE> elemSol;
    TYPE elemDiv;

    for( it.Begin(); !it.IsEnd(); it++ ) {
      const Elem * ptElem = it.GetElem();

      // get coordinates of element
      ptgrid_->GetElemNodesCoord( cornerCoords, ptElem->connect, true );

      // fetch global gradient
      ptElem->ptElem->GetCoordMidPoint( elemMidPoint );
      ptElem->ptElem->GetEdgeGlobalDerivShapeFnc( globDeriv, elemMidPoint, 
                                                  cornerCoords, ptElem );

      // getch solution of element
      sol_->GetElemSolution(elemSol, it);

      // loop over shape functions
      elemDiv = 0.0;
      for( UInt i = 0; i < elemSol.GetSize(); i++ )  {

        elemDiv += elemSol[i] * ( globDeriv[i][0][0] 
                                + globDeriv[i][1][1] 
                                + globDeriv[i][2][2] );
      }
      actVal[it.GetPos()] = elemDiv;
    }
  }

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
        ptgrid_->RegionNameToId( regionIds, couplRegions );

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
    
    //get the coupling regions
    StdVector<std::string> couplRegions;
    StdVector<RegionIdType> regionIds;
    ptCoupling_->GetOutputRegions(actCoupling, couplRegions);
    ptgrid_->RegionNameToId( regionIds, couplRegions );

   
    force.Init();
    Vector<Double>  fluxAtIp(dim_), elemForce, fAtIp, jAtIp;
    Matrix<Double> ptCoord;

    for (UInt reg=0; reg<couplRegions.GetSize(); reg++) {

      //find subdomain index
      Integer sdIndex = subdoms_.Find( regionIds[reg] );
      if( sdIndex == -1 ) {
        EXCEPTION( "The region coupling region '" <<
                   ptgrid_->RegionIdToName( regionIds[reg] )
                   << "' was not found in magneticPDE" );
      }
      
      RegionIdType actRegionId = subdoms_[sdIndex] ;
      Double conductivity;
      materials_[actRegionId]->GetScalar(conductivity,MAG_CONDUCTIVITY,Global::REAL);
      
      // Check if this region is a coil
      Integer coilIndex = coilRegionId_.Find(actRegionId);
      
      ElemList actSDList(ptgrid_ );
      actSDList.SetRegion(actRegionId );
      
       EntityIterator it = actSDList.GetIterator();
       UInt actEl = 0;
       // iterate over all elements of regions
       for ( it.Begin(); !it.IsEnd(); it++, actEl++ ) {
         
         BaseFE * ptElem = it.GetElem()->ptElem;

         //need info for ansatz functions of mechanical FE!
         shared_ptr<AnsatzFct> fct(new LagrangeFct);
         ptElem->SetAnsatzFct( fct );
         UInt numFncs = ptElem->GetNumFncs(fct);
//          ptElem->SetAnsatzFct( results_[0]->fctType );
//          UInt numFncs = ptElem->GetNumFncs(results_[0]->fctType);

         const UInt nrIntPts = ptElem->GetNumIntPoints();
         const Vector<Double> & intWeights = ptElem->GetIntWeights();  
         
         Vector<Double> shpFncAtIp;         
         ptgrid_->GetElemNodesCoord( ptCoord, 
                                     it.GetElem()->connect,
                                     true );
         
         // iterate over all integration points
         fAtIp.Resize( dim_ * numFncs );
         elemForce.Resize( dim_ * numFncs );
         elemForce.Init();
         for( UInt ip = 1; ip <= nrIntPts; ip++ ) {
           ptElem->GetShFncAtIp(shpFncAtIp, ip, it.GetElem() );

           // CHECK: If this region is a current coil, we simply take the
           // prescribed current density value
           if( coilIndex != -1 ) {
             MathParser * mParser =  domain->GetMathParser();
             std::string factor = coilDef_[coilIndex]->value_ + "/" 
               + GenStr(coilDef_[coilIndex]->windingCrossSection_ );
             mParser->SetExpr( mHandle_, factor );
             // TODO: Check if this is still needed
             // Double currDens = mParser->Eval(mHandle_);
           } 
           else {
             CalcEddyCurrentAtIP( it, 0, jAtIp );
           }
           
           CalcFluxDensityAtIP( it, 0, fluxAtIp );
           
           // calculate cross product 
           fAtIp.Init();
           Vector<Double> tempCross;
           jAtIp.CrossProduct( fluxAtIp, tempCross );
           for (UInt iFnc=0; iFnc<numFncs; iFnc++) {
             fAtIp[iFnc*3+0] = -tempCross[0] * shpFncAtIp[iFnc];
             fAtIp[iFnc*3+1] = -tempCross[1] * shpFncAtIp[iFnc];
             fAtIp[iFnc*3+2] = -tempCross[2] * shpFncAtIp[iFnc];
           }

           Double jacDet = ptElem->CalcJacobianDetAtIp(ip, ptCoord, 
                                                       it.GetElem() );
           elemForce += -fAtIp * (jacDet * intWeights[ip-1] );

         }
         StdVector<UInt> const & connecth = it.GetElem()->connect;

         // Add the element force to the according coupling node
         for (UInt ielemnode=0; ielemnode<connecth.GetSize(); ielemnode++) {
           UInt pos = cplNodeNumPos[connecth[ielemnode]];
           for( UInt idim=0; idim<dim_; idim++) {
             force[pos*dim_+idim] += elemForce[ielemnode*dim_+idim];
           }
         }
       }
       
 
     }
  }


  bool MagEdgePDE::HasOutput( SolutionType output ) {
    if (output == MAG_FORCE_LORENTZ) {
      return true;
    }
    return false;
  }

} // end of namespace

