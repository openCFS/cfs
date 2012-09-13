#include <fstream>

#include "MagEdgePDE.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Utils/Coil.hh"
#include "Utils/SmoothSpline.hh"
#include "Utils/LinInterpolate.hh"

#include "CoupledPDE/PDECoupling.hh"
#include "Driver/Assemble.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"
#include "FeBasis/HCurl/FeSpaceHCurlHi.hh"
#include "FeBasis/HCurl/HCurlElems.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

#include "Domain/CoefFunction/CoefFunctionExpression.hh"
#include "Domain/CoefFunction/CoefFunctionApprox.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefFunctionMulti.hh"
#include "Domain/CoefFunction/CoefXpr.hh"

// forms
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/BDUInt.hh"
#include "Forms/LinForms/KXInt.hh"
#include "Forms/Operators/CurlOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"

//time stepping
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"
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
    
    
    reluc_.reset(new CoefFunctionMulti());
    conduc_.reset(new CoefFunctionMulti());
  }


  // *************
  //  Destructor
  // *************
  MagEdgePDE::~MagEdgePDE() {
  }



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

    SinglePDE::InitNonLin();
  }


  // *****************************
  //   Definition of Integrators
  // *****************************
  void MagEdgePDE::DefineIntegrators() {

    RegionIdType actRegion;
    BaseMaterial * actMat = NULL;

    // flag for updatedLagrange formulation
    //bool upLagrangeForm = true;

    
    // initially, check for regularization factor
    Double regularizationFactor = 1e-6;
    myParam_->GetValue("penaltyFactor", regularizationFactor, ParamNode::PASS);
    
    //==============================================================
    //begin new implementation
    //==============================================================
    shared_ptr<FeSpace> mySpace = feFunctions_[MAG_POTENTIAL]->GetFeSpace();
    for(UInt iRegion = 0; iRegion < regions_.GetSize() ; iRegion ++){
      actRegion = regions_[iRegion];
      actMat    = materials_[actRegion];
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);

      shared_ptr<BaseFeFunction> feFunc = feFunctions_[MAG_POTENTIAL];
      shared_ptr<FeSpace> feSpace = feFunc->GetFeSpace();

      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();
      mySpace->SetRegionApproximation(actRegion, polyId,integId);
      
      //get possible nonlinearities defined in this region
      StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion]; 

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( actRegion );
      
      // pass entitylist ot fespace / fefunction
      feFunc->AddEntityList( actSDList );
      
      // Switch, if region is linear / nonlinear
      if ( nonLinTypes.GetSize() > 0 ) {
        
        
        // obtain BH-curve from nonlinear material and initialize it
        actMat->NeedApproxMatCurve( MAG_PERMEABILITY );
        actMat->InitApproxCurves();
        ApproxData * nLinFnc = actMat->GetNonlinFnc(MAG_PERMEABILITY);
        assert(nLinFnc);
        Double nuLin;
        actMat->GetScalar(nuLin , MAG_RELUCTIVITY, Global::REAL );
        // obtain initial value for material (cast to double)
        
        // ================================
        //  Nonlinear Stiffness Integrator 
        // =================================
        
        // create CoefFunctionApprox with nlinFnc, initial value and feFunction
        shared_ptr<CoefFunctionApprox<CurlOperator<FeHCurl,3, Double> > > nuNl (
            new CoefFunctionApprox<CurlOperator<FeHCurl,3, Double> >() );
        nuNl->Init( nuLin, nLinFnc, 
                    dynamic_pointer_cast<FeFunction<Double > > (feFunc) );
        // create stiffness integrator
        BaseBDBInt* stiff1 = NULL;
        stiff1 = new BBInt< CurlOperator<FeHCurl,3, Double> >(nuNl, 1.0) ;
        stiff1->SetName("CurlCurlIntegrator-NL");

       BiLinFormContext * stiffContext =
           new BiLinFormContext(stiff1, STIFFNESS );
       stiffContext->SetEntities( actSDList, actSDList );
       stiffContext->SetFeFunctions( feFunc, feFunc );
       assemble_->AddBiLinearForm( stiffContext );
       // Important: Add bdb-integrator to global list, as we need them later
       // for calculation of postprocessing results
       bdbInts_[actRegion] = stiff1;
       // add also material to global, distributed reluctivity coefficient function
       //reluc_->AddRegion(actRegion, nuNl);

       // ================================================
       //  Nonlinear Stiffness Integrator (only Newton )
       // ================================================
       // Note: currently we set the nonlinear method hard-coded to NEWTON for
       // testing purpose
       if( nonLinMethod_ == NEWTON ) {
         //create CoefFunctionApproxDeriv with nlinFnc and feFunction
         shared_ptr<CoefFunctionApproxDeriv<CurlOperator<FeHCurl,3, Double> > > nuDeriv (
             new CoefFunctionApproxDeriv<CurlOperator<FeHCurl,3, Double> >() );
         nuDeriv->Init( nLinFnc, dynamic_pointer_cast<FeFunction<Double > > (feFunc) );

         //create stiffness integrator
         BiLinearForm* stiff2 = NULL;
         stiff2 = new BDBInt< CurlOperator<FeHCurl,3, Double> >(nuDeriv, 1.0) ;
         stiff2->SetName("CurlCurlIntegrator-NL-Newton");

         BiLinFormContext * stiffContext2 =
             new BiLinFormContext(stiff2, STIFFNESS );
         stiffContext2->SetEntities( actSDList, actSDList );
         stiffContext2->SetFeFunctions( feFunc, feFunc );
         assemble_->AddBiLinearForm( stiffContext2 );
       }
        
       // =================================
       //  Nonlinear RHS-integrator
       // =================================
        LinearForm * rhsNlinForm = new KXIntegrator<Double>(stiff1, -1.0, feFunc );
        rhsNlinForm->SetName("RHSNonLinForm");
        LinearFormContext * rhsNlinContext =
            new LinearFormContext( rhsNlinForm );
        rhsNlinContext->SetEntities( actSDList );
        rhsNlinContext->SetFeFunction( feFunc );
        assemble_->AddLinearForm( rhsNlinContext );
      } else {

       // ***************************************
       // LINEAR PART
       // ***************************************

       // ===============================
       //  Standard Stiffness Integrator
       // ===============================
        PtrCoefFct curCoef = 
            actMat->GetTensorCoefFnc(MAG_RELUCTIVITY,FULL,Global::REAL );
        BaseBDBInt* curlcurl;
        curlcurl = new BDBInt< CurlOperator<FeHCurl,3, Double> >(curCoef,1.0) ;
        curlcurl->SetName("CurlCurlIntegrator");

       BiLinFormContext * stiffContext =
           new BiLinFormContext(curlcurl, STIFFNESS );
       stiffContext->SetEntities( actSDList, actSDList );
       stiffContext->SetFeFunctions( feFunc, feFunc );
       assemble_->AddBiLinearForm( stiffContext );
       // Important: Add bdb-integrator to global list, as we need them later
       // for calculation of postprocessing results
       bdbInts_[actRegion] = curlcurl;
       // add also material to global, distributed reluctivity coefficient function
       reluc_->AddRegion(actRegion, curCoef);

       // === Additional RHS integrator in case of Non-linearity ===
       if ( nonLin_ == true ) {
         REFACTOR;
         // =================================
         //  Nonlinear RHS-integrator
         // =================================
         LinearForm * rhsNlinForm = 
             new KXIntegrator<Double>(curlcurl, -1.0, feFunc );
         rhsNlinForm->SetName("RHSNonLinForm-Lin");
         LinearFormContext * rhsNlinContext =
             new LinearFormContext( rhsNlinForm );
         rhsNlinContext->SetEntities( actSDList );
         rhsNlinContext->SetFeFunction( feFunc );
         assemble_->AddLinearForm( rhsNlinContext );
       }
     } // END OF NONLIN / LIN PART

      // ============================
      // Standard Mass Matrix
      // ============================
      Double conductivity = 0.0; // 
      materials_[actRegion]->GetScalar(conductivity,MAG_CONDUCTIVITY,Global::REAL);
      bool scaleByEdgeSize = false;
      if ( conductivity < 1e-10 || analysistype_ == STATIC ) {
        Matrix<Double> reluc; 

        // get tensor of permeability and determine max. value
        materials_[actRegion]->GetTensor( reluc, MAG_RELUCTIVITY, Global::REAL );
        conductivity =  regularizationFactor * reluc[0][0];
        scaleByEdgeSize = true;
      }

      PtrCoefFct coeff =
          CoefFunction::Generate(Global::REAL, lexical_cast<std::string>(conductivity));
      // add also material to global, distributed reluctivity coefficient function
      conduc_->AddRegion(actRegion, coeff);
      BaseBDBInt *massInt;
     
      BiLinFormContext * massContext;
      if ( analysistype_ == STATIC) {
        // we have to guarantee, that we add some mass to curl-curl integrator.
        // Additionally, the integrator gets scaled by the edge size for a uniform
        // conditioning
        massInt = new BBIntMassEdge<ScaledByEdgeIdentityOperator<3,Double> >(coeff,1.0);
        massInt->SetName("MassIntegrator");
        massContext =  new BiLinFormContext(massInt, STIFFNESS );
      } else {
        // here we add the "normal" mass integrator, which gets not scaled by the 
        // edge size
        if( scaleByEdgeSize ) {
          massInt = new BBIntMassEdge<ScaledByEdgeIdentityOperator<3,Double> >(coeff,1.0);
        } else {
          massInt = new BBIntMassEdge<IdentityOperator<FeHCurl,3,1,Double> >(coeff,1.0);
        }
        massInt->SetName("MassIntegrator");
        massContext =
            new BiLinFormContext(massInt, MASS );
      }
      massContext->SetEntities( actSDList, actSDList );
      massContext->SetFeFunctions( feFunc, feFunc );
      assemble_->AddBiLinearForm( massContext );
      
      // insert mass integrator to list of defined mass integrators
      massInts_[actRegion] = massInt;

      // ============================
      // COIL INTEGRATORS
      // ============================
      Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
      // If this subdomain is a coil we have to do special things
      for ( UInt coil = 0; coil < coilDef_.GetSize(); coil++ ) {
        if ( actRegion == coilRegionId_[coil] ) {
         LinearForm* curInt;
          
          std::string factor = coilDef_[coil]->value_ + "/" +
            lexical_cast<std::string>(coilDef_[coil]->windingCrossSection_);

          StdVector<std::string> currDensity(3), imag(3);
          currDensity[0] = factor + "*" + lexical_cast<std::string>(coilDef_[coil]->locFlowDir_[0]);
          currDensity[1] = factor + "*" + lexical_cast<std::string>(coilDef_[coil]->locFlowDir_[1]);
          currDensity[2] = factor + "*" + lexical_cast<std::string>(coilDef_[coil]->locFlowDir_[2]);
          imag[0] = AmplPhaseToImag(currDensity[0], coilDef_[coil]->phase_ );
          imag[1] = AmplPhaseToImag(currDensity[1], coilDef_[coil]->phase_ );
          imag[2] = AmplPhaseToImag(currDensity[2], coilDef_[coil]->phase_ );
          
          PtrCoefFct coef(CoefFunction::Generate(part, currDensity, imag));
          coef->SetCoordinateSystem(coilDef_[coil]->flowCoordSys_);
          
          // remember coefficient for later use
          coilCoefs_[actRegion] = coef;
          
          if( isComplex_ ) {
            curInt = new BUIntegrator<IdentityOperator<FeHCurl,3,1,Complex>,Complex >(1.0, coef);
          }
          else {
            curInt = new BUIntegrator<IdentityOperator<FeHCurl,3,1,Double>,Double >(1.0, coef);
          }
          curInt->SetName("CoilIntegrator");
          LinearFormContext * coilContext =
              new LinearFormContext( curInt );
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
  
  
  void MagEdgePDE::DefineRhsLoadIntegrators() {
    LOG_TRACE(magEdgePde) << "Defining rhs load integrators for magEdgePDE";
       
    // Get FESpace and FeFunction of mechanical displacement
    shared_ptr<BaseFeFunction> myFct = feFunctions_[MAG_POTENTIAL];
    shared_ptr<FeSpace> mySpace = myFct->GetFeSpace();

    StdVector<shared_ptr<EntityList> > ent;
    StdVector<PtrCoefFct > coef;
    LinearForm * lin = NULL;
    StdVector<std::string> vecDofNames = myFct->GetResultInfo()->dofNames;
       
    // ==================
    //  FLUX DENSITY
    // ==================
    LOG_DBG(magEdgePde) << "Reading prescribed flux density";

    ReadRhsExcitation( "fluxDensity", vecDofNames, ResultInfo::VECTOR, isComplex_, 
                       ent, coef );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST ||
          ent[i]->GetType() == EntityList::SURF_ELEM_LIST ) {
        EXCEPTION("Prescribed magnetic flux density can only be defined im volume")
      }
      if(isComplex_) {
        lin = new BDUIntegrator<CurlOperator<FeHCurl,3, Double>, Complex>(Complex(1.0), coef[i],
                                                                          reluc_);
      } else {
        lin = new BDUIntegrator<CurlOperator<FeHCurl,3, Double>, Double>(1.0, coef[i],
                                                                         reluc_ );
      }
      lin->SetName("FluxIntegrator");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
    } // for
  }
  

  void MagEdgePDE::DefineSolveStep()
  {
    SolveStepMagEdge *magSolveStep = new SolveStepMagEdge(*this); 
    
    solveStep_ = magSolveStep; 
  }


  // ======================================================
  // TIME-STEPPING SECTION
  // ======================================================

  void MagEdgePDE::InitTimeStepping() {
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(TimeSchemeGLM::TRAPEZOIDAL, 0) );

    feFunctions_[MAG_POTENTIAL]->SetTimeScheme(myScheme);

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
      WARN("Adapt printing of coils to InfoNode");
      //Info->PrintF( pdename_, "Using the following coils:\n" );
      for( UInt i = 0; i < coilNodes.GetSize(); i++ ) {

        // get region name of actual coil
        std::string regionName = coilNodes[i]->Get("name")->As<std::string>();
        RegionIdType regionId = ptGrid_->GetRegion().Parse( regionName );

        coilRegionId_.Push_back( regionId );
        coilDef_.Push_back( shared_ptr<Coil>( new Coil( regionId,
                                                        coilNodes[i], ptGrid_) ) );
        //Info->PrintCoil( *coilDef_.Last(), analysistype_ );
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
      WARN("Adjust printing of permanent magnet definition to InfoNode");
//      Info->PrintF( pdename_,
//              "Found permanent magnets in the following regions:\n" );

      Double tmpDir;
      for( UInt i = 0; i < magnetNodes.GetSize(); i++ ) {

        // get region name of actual magnet
        std::string regionName = magnetNodes[i]->Get("name")->As<std::string>();
        RegionIdType regionId = ptGrid_->GetRegion().Parse( regionName );

        magnetsDomain_.Push_back( regionId );

        // read orientation
        magnetNodes[i]->GetValue( "orientX", tmpDir );
        magnetsOriX_.Push_back( tmpDir );

        magnetNodes[i]->GetValue( "orientY", tmpDir );
        magnetsOriY_.Push_back( tmpDir );

        magnetNodes[i]->GetValue( "orientZ", tmpDir );
        magnetsOriZ_.Push_back( tmpDir );

        // report name to logfile
        //Info->PrintF( pdename_, " %s\n", regionName.c_str());
      }
    }
  }


  void MagEdgePDE::DefinePrimaryResults() {

    StdVector<std::string> vecComponents;
    vecComponents = "x", "y", "z";

    // === MAGNETIC VECTOR POTENTIAL ===
    shared_ptr<ResultInfo> res1(new ResultInfo);
    res1->resultType = MAG_POTENTIAL;
    res1->dofNames = vecComponents;
    res1->unit = "Vs/m";
    res1->definedOn = ResultInfo::ELEMENT;
    res1->entryType = ResultInfo::VECTOR;
    
    results_.Push_back( res1 );
    availResults_.insert( res1 );
    
    feFunctions_[MAG_POTENTIAL]->SetResultInfo(res1);
    DefineFieldResult( feFunctions_[MAG_POTENTIAL], res1 );
   
    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    hdbcSolNameMap_[MAG_POTENTIAL] = "fluxParallel";
    idbcSolNameMap_[MAG_POTENTIAL] = "potential";
  }
  
  void MagEdgePDE::DefinePostProcResults() {

      StdVector<std::string> vecComponents;
      vecComponents = "x", "y", "z";
      
      Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
      shared_ptr<BaseFeFunction> feFct = feFunctions_[MAG_POTENTIAL];
      
      // === MAGNETIC VECTOR POTENTIAL - 1ST DERIVATIVE ===
      if( analysistype_ == TRANSIENT || analysistype_ == HARMONIC ) {
        shared_ptr<ResultInfo> aDot(new ResultInfo);
        aDot->resultType = MAG_POTENTIAL_DERIV1;
        aDot->dofNames = vecComponents;
        aDot->unit = "V/m";
        aDot->definedOn = ResultInfo::ELEMENT;
        aDot->entryType = ResultInfo::VECTOR;
        availResults_.insert( aDot );
        DefineTimeDerivResult( MAG_POTENTIAL_DERIV1, 1, MAG_POTENTIAL );
      }
      
      // === MAGNETIC FLUX DENSITY ===
      shared_ptr<ResultInfo> flux(new ResultInfo);
      flux->resultType = MAG_FLUX_DENSITY;
      flux->dofNames = vecComponents;
      flux->unit = "Vs/m^2";
      flux->definedOn = ResultInfo::ELEMENT;
      flux->entryType = ResultInfo::VECTOR;
      availResults_.insert( flux );
      shared_ptr<CoefFunctionFormBased> bFunc;
      if( isComplex_ ) {
        bFunc.reset(new CoefFunctionBOp<Complex>(feFct, flux));
      } else {
        bFunc.reset(new CoefFunctionBOp<Double>(feFct, flux));
      }
      DefineFieldResult( bFunc, flux );
      
      
      // === MECHANIC RHS ===
      shared_ptr<ResultInfo> rhs(new ResultInfo);
      rhs->resultType = MAG_RHS_LOAD;
      rhs->dofNames = vecComponents;
      rhs->unit = "-";
      rhs->entryType = ResultInfo::VECTOR;
      rhs->definedOn = ResultInfo::ELEMENT;
      availResults_.insert( rhs );
      rhsFeFunctions_[MAG_POTENTIAL]->SetResultInfo(rhs);
      DefineFieldResult( rhsFeFunctions_[MAG_POTENTIAL], rhs );
      //DefineFieldResult()
      
      
      // === EDDY CURRENT DENSITY ===
      shared_ptr<CoefFunctionFormBased> jFunc;
      if( analysistype_ != STATIC ) {
        shared_ptr<BaseFeFunction> aDotFct = 
            timeDerivFeFunctions_[MAG_POTENTIAL_DERIV1];
        shared_ptr<ResultInfo> eddy(new ResultInfo);
        eddy->resultType = MAG_EDDY_CURRENT_DENSITY;
        eddy->dofNames = vecComponents;
        eddy->unit = "A/m^2";
        eddy->definedOn = ResultInfo::ELEMENT;
        eddy->entryType = ResultInfo::VECTOR;
        availResults_.insert( eddy );

        if( isComplex_ ) {
          jFunc.reset(new CoefFunctionFlux<Complex>(aDotFct, eddy, -1.0));
        } else {
          jFunc.reset(new CoefFunctionFlux<Double>(aDotFct, eddy, -1.0));
        }
        DefineFieldResult( jFunc, eddy );
      }
      
      // === COIL CURRENT DENSITY ===
      shared_ptr<ResultInfo> ccd(new ResultInfo);
      ccd->resultType = MAG_COIL_CURRENT_DENSITY;
      ccd->dofNames = vecComponents;
      ccd->unit = "A/m^2";
      ccd->definedOn = ResultInfo::ELEMENT;
      ccd->entryType = ResultInfo::VECTOR;
      availResults_.insert( ccd );
      shared_ptr<CoefFunctionMulti> ccdCoef(new CoefFunctionMulti());
      // loop over all coil coefficients and add contribution to coef 
      std::map<RegionIdType, PtrCoefFct>::iterator coilIt = coilCoefs_.begin();
      for( ; coilIt != coilCoefs_.end(); ++coilIt ) {
        ccdCoef->AddRegion( coilIt->first, coilIt->second);
      }
      DefineFieldResult( ccdCoef, ccd );
      
      
      // === TOTAL CURRENT DENSITY ===
      shared_ptr<ResultInfo> tcd(new ResultInfo);
      tcd->resultType = MAG_TOTAL_CURRENT_DENSITY;
      tcd->dofNames = vecComponents;
      tcd->unit = "A/m^2";
      tcd->definedOn = ResultInfo::ELEMENT;
      tcd->entryType = ResultInfo::VECTOR;
      availResults_.insert( tcd );
      shared_ptr<CoefFunctionMulti> tcdCoef(new CoefFunctionMulti());
      // loop over all regions and assemble total current density:
      //  - if region is coil -> take coil current
      //  - if region is no coil and analyis is transient/harmonic -> eddy
      
      StdVector<RegionIdType>::iterator regIt = regions_.Begin();
      for( ; regIt != regions_.End(); ++regIt ) {
        RegionIdType actRegion = *regIt;
        if( coilCoefs_.find(actRegion) != coilCoefs_.end() ) {
          // region is a coil
          tcdCoef->AddRegion( actRegion, coilCoefs_[actRegion] );
        } else {
          // region is no coil
          if( analysistype_ == TRANSIENT || analysistype_ == HARMONIC ) {
            tcdCoef->AddRegion( actRegion, jFunc );
          }
        }
      }
      DefineFieldResult( tcdCoef, tcd );
      
      // === LORENTZ FORCE DENSITY ===
      if(analysistype_ == TRANSIENT ||
          analysistype_ == HARMONIC)  {
        shared_ptr<ResultInfo> lfd(new ResultInfo);
        lfd->resultType = MAG_FORCE_LORENTZ_DENSITY;
        lfd->dofNames = vecComponents;
        lfd->unit = "N/m^3";
        lfd->definedOn = ResultInfo::ELEMENT;
        lfd->entryType = ResultInfo::VECTOR;
        availResults_.insert( lfd );

        // assemble coefficient function F_L = J X B
        PtrCoefFct lfdFunc = CoefFunction::Generate( part, 
                 CoefXprBinOp( tcdCoef, bFunc, CoefXpr::OP_CROSS ) );
        DefineFieldResult( lfdFunc, lfd);

        // === LORENTZ FORCE (TOTAL) ===
        shared_ptr<ResultInfo> lf(new ResultInfo);
        lf->resultType = MAG_FORCE_LORENTZ;
        lf->dofNames = vecComponents;
        lf->unit = "N";
        lf->definedOn = ResultInfo::REGION;
        lf->entryType = ResultInfo::VECTOR;
        availResults_.insert( lf );

        // build result functor for integration
        shared_ptr<ResultFunctor> lfFunc;
        if( isComplex_ ) {
          lfFunc.reset(new ResultFunctorIntegrate<Complex>(lfdFunc, feFct, lf ) );
        } else {
          lfFunc.reset(new ResultFunctorIntegrate<Double>(lfdFunc, feFct, lf ) );
        }
        resultFunctors_[MAG_FORCE_LORENTZ] = lfFunc;
      }

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
      shared_ptr<ResultFunctor> energyFunc;
      if( isComplex_ ) {
        energyFunc.reset(new EnergyResultFunctor<Complex>(feFct, energy));
      } else {
        energyFunc.reset(new EnergyResultFunctor<Double>(feFct, energy));
      }
      resultFunctors_[MAG_ENERGY] = energyFunc;
      
      // ============================
      // Initialize result functors:
      // ============================
      // 1) Loop over all BDB-integrators
      std::map<RegionIdType, BaseBDBInt*>::iterator it = bdbInts_.begin();
      for( ; it != bdbInts_.end(); ++it ) {
        RegionIdType region = it->first;
        BaseBDBInt* bdb = it->second;

        // 2) pass integrators to functors
        bFunc->AddIntegrator(bdb, region);
        energyFunc->AddIntegrator(bdb, region);
      }
      
      // 2) Loop over all MASS-integrators
      it = massInts_.begin();
      for( ; it != massInts_.end(); ++it ) {
        RegionIdType region = it->first;
        BaseBDBInt* mass = it->second;

        // 2) pass integrators to functors
        if(jFunc )
          jFunc->AddIntegrator( mass, region );
      }
    }

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

      //get possible nonlinearities defined in this region
      StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion]; 

      // Check, if region is nonlinear
      if ( nonLinTypes[actRegion] == PERMEABILITY ) {

        // Obtain nonlinear approximation functional
        ApproxData * approx  = materials_[actRegion]->GetNonlinFnc(MAG_PERMEABILITY);

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
//      ptGrid_->GetElemNodesCoord( cornerCoords, ptElem->connect, true );
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
        ptGrid_->GetRegion().Parse( couplRegions , regionIds );

        // Check, that every coupling region is part of
        // the magnetic pde itself
        for( UInt iReg = 0; iReg < couplRegions.GetSize(); iReg++ ) {
          if( regions_.Find(regionIds[iReg]) == -1 ) {
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
EXCEPTION("Needs to be re-implemented");
//    SolutionType quantity;
//    StdVector<UInt> * couplingNodes = NULL;
//    SingleVector * values = NULL;
//    UInt forcesCount = 0;
//
//    // at first, check if this PDE is iterative coupled
//    if (isIterCoupled_ == false)
//      return;
//
//    // loop over all output coupling quantities
//    for ( UInt actCoupling = 0;
//          actCoupling < ptCoupling_->GetNumOutputCouplings();
//          actCoupling++ ) {
//
//      quantity = ptCoupling_->GetOutputQuantity(actCoupling);
//      ptCoupling_->GetOutputValues(actCoupling, values);
//      Vector<Double> *temp = dynamic_cast<Vector<Double> *>(values);
//
//      switch(ptCoupling_->GetOutputType(actCoupling)) {
//
//      case NODE:
//        ptCoupling_->GetOutputNodes(actCoupling, couplingNodes);
//
//        if (quantity == MAG_FORCE_LORENTZ) {
//          CalcNodeForceLorentz( *temp, cplNodeNumPos_[forcesCount],
//                                actCoupling, couplingNodes->GetSize());
//          forcesCount++;
//        }
//
//        break;
//
//      case ELEM:
//        EXCEPTION( "No Element coupling output" );
//        break;
//      }
//    }
  }

  std::map<SolutionType, shared_ptr<FeSpace> > 
  MagEdgePDE::CreateFeSpaces(const std::string& formulation,
                             PtrParamNode infoNode ) {
    //ok default case so we create grid based approximation H1 elements
    //and standard Gauss integration
    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    if(formulation == "default" || formulation == "H_CURL"){
      PtrParamNode potSpaceNode = infoNode->Get("magPotential");
      crSpaces[MAG_POTENTIAL] = 
          FeSpace::CreateInstance(myParam_, potSpaceNode, FeSpace::HCURL, ptGrid_ );
      crSpaces[MAG_POTENTIAL]->Init(solStrat_);
    }else{
      EXCEPTION("The formulation " << formulation 
                << "of magnetic edge PDE is not known!");
    }
    
    
    // in addition query, if special treatment of anisotropic elements
    // is activated
    if( myParam_->Has("thinElements") ) {
      Double aspectRatio = 0.0;
      aspectRatio = myParam_->Get("thinElements")->Get("maxAspectRatio")->As<Double>();
      dynamic_pointer_cast<FeSpaceHCurlHi>(crSpaces[MAG_POTENTIAL])
          ->TreatThinElements(aspectRatio);
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

