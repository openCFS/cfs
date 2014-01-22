#include <fstream>

#include "MagEdgePDE.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Utils/Coil.hh"
#include "Utils/SmoothSpline.hh"
#include "Utils/LinInterpolate.hh"

#include "Driver/Assemble.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"
#include "FeBasis/HCurl/FeSpaceHCurlHi.hh"
#include "FeBasis/HCurl/HCurlElems.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

#include "Domain/CoefFunction/CoefFunctionExpression.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefFunctionMulti.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
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
  MagEdgePDE::MagEdgePDE( Grid * aptgrid, PtrParamNode paramNode,
                          PtrParamNode infoNode,
                          shared_ptr<SimState> simState, Domain* domain )
    :SinglePDE( aptgrid, paramNode, infoNode, simState, domain ) {

    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_          = "magneticEdge";
    pdematerialclass_ = ELECTROMAGNETIC;
    
    //! Always use updated Lagrangian formulation 
    updatedGeo_        = true;
      
    // check if we have a 3d setup
    bool is3d = domain_->GetParamRoot()->Get("domain")->Get("geometryType")->As<std::string>() == "3d";
    if ( !is3d )
      EXCEPTION("MagEdgePDE is just implemented for 3D setups!");
    
    
    reluc_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, dim_, dim_, isComplex_));
    conduc_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1, isComplex_));
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
        
        // ================================
        //  Nonlinear Stiffness Integrator 
        // =================================
        
        // create stiffness integrator
        //BaseBOperator* bOp = new CurlOperator<FeHCurl,3, Double>();
        PtrCoefFct magFluxCoef = this->GetCoefFct(MAG_FLUX_DENSITY);
        PtrCoefFct nuNl = 
            actMat->GetScalCoefFncNonLin( MAG_RELUCTIVITY, Global::REAL, 
                                          magFluxCoef  );
        BaseBDBInt* stiff1 = NULL;
        stiff1 = new BBInt<>(new  CurlOperator<FeHCurl,3, Double>(), nuNl, 1.0, updatedGeo_) ;
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
       reluc_->AddRegion(actRegion, nuNl);

       // ================================================
       //  Nonlinear Stiffness Integrator (only Newton )
       // ================================================
       // Note: currently we set the nonlinear method hard-coded to NEWTON for
       // testing purpose
       if( nonLinMethod_ == NEWTON ) {

         //BaseBOperator* bOp = new CurlOperator<FeHCurl,3, Double>();
         PtrCoefFct magFluxCoef = this->GetCoefFct(MAG_FLUX_DENSITY);
         PtrCoefFct nuDeriv = actMat->GetTensorCoefFncNonLin( MAG_RELUCTIVITY_DERIV, FULL,  
                                                            Global::REAL, magFluxCoef );
         
         //create stiffness integrator
         BiLinearForm* stiff2 = NULL;
         stiff2 = new BDBInt<>(new CurlOperator<FeHCurl,3, Double>(), nuDeriv, 1.0, updatedGeo_) ;
         stiff2->SetName("CurlCurlIntegrator-NL-Newton");

         BiLinFormContext * stiffContext2 =
             new BiLinFormContext(stiff2, STIFFNESS );
         stiffContext2->SetEntities( actSDList, actSDList );
         stiffContext2->SetFeFunctions( feFunc, feFunc );
         assemble_->AddBiLinearForm( stiffContext2 );
       }
        
       if ( analysistype_ == STATIC ) {
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
       }
      } else {

       // ***************************************
       // LINEAR PART
       // ***************************************

       // ===============================
       //  Standard Stiffness Integrator
       // ===============================
        PtrCoefFct curCoef = 
            //actMat->GetTensorCoefFnc(MAG_RELUCTIVITY,FULL,Global::REAL );
            actMat->GetScalCoefFnc(MAG_RELUCTIVITY,Global::REAL );
        BaseBDBInt* curlcurl;
        //curlcurl = new BDBInt< CurlOperator<FeHCurl,3, Double> >(curCoef,1.0) ;
        curlcurl = new BBInt<>(new  CurlOperator<FeHCurl,3, Double>(), curCoef,1.0, updatedGeo_) ;
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
       if ( nonLin_ == true && analysistype_ == STATIC ) {
         //REFACTOR;
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
      // use gradient of shape functions?
      bool useGrad = true;
      if ( conductivity < 1e-10 || analysistype_ == STATIC ) {
        // do not use gradients for non-conductive regions (for regularization 
        // only the lowest order mass term is used)
        useGrad = false;
        
        Matrix<Double> reluc; 

        // get tensor of permeability and determine max. value
        materials_[actRegion]->GetTensor( reluc, MAG_RELUCTIVITY, Global::REAL );
        conductivity =  regularizationFactor * reluc[0][0];
        scaleByEdgeSize = true;
        
        // add region to set of "regularized" regions
        regularizedRegions_.insert(actRegion);
      }
      if(useGrad) {
        dynamic_pointer_cast<FeSpaceHCurlHi>(feSpace)->SetUseGradients(actRegion);
      }
      
      PtrCoefFct coeff =
          CoefFunction::Generate(mp_, Global::REAL, lexical_cast<std::string>(conductivity));
      // add also material to global, distributed reluctivity coefficient function
      conduc_->AddRegion(actRegion, coeff);
      BaseBDBInt *massInt;
     
      BiLinFormContext * massContext;
      if ( analysistype_ == STATIC) {
        // we have to guarantee, that we add some mass to curl-curl integrator.
        // Additionally, the integrator gets scaled by the edge size for a uniform
        // conditioning
        massInt = new BBIntMassEdge<>(new ScaledByEdgeIdentityOperator<3,Double>(),
                                      coeff,1.0);
        massInt->SetName("MassIntegrator");
        massContext =  new BiLinFormContext(massInt, STIFFNESS );
      } else {
        // here we add the "normal" mass integrator, which gets not scaled by the 
        // edge size
        if( scaleByEdgeSize ) {
          massInt = new BBIntMassEdge<>(new ScaledByEdgeIdentityOperator<3,Double>(),
                                        coeff,1.0, updatedGeo_);
        } else {
          massInt = new BBIntMassEdge<>(new IdentityOperator<FeHCurl,3,1,Double>(),
                                        coeff,1.0, updatedGeo_);
        }
        massInt->SetName("MassIntegrator");
        massContext =
            new BiLinFormContext(massInt, DAMPING );
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
          
          PtrCoefFct coef(CoefFunction::Generate(mp_, part, currDensity, imag));
          coef->SetCoordinateSystem(coilDef_[coil]->flowCoordSys_);
          
          // remember coefficient for later use
          coilCoefs_[actRegion] = coef;
          
          if( isComplex_ ) {
            curInt = new BUIntegrator<Complex>( new IdentityOperator<FeHCurl,3,1,Complex>(),
                                                1.0, coef, updatedGeo_);
          }
          else {
            curInt = new BUIntegrator<Double>( new IdentityOperator<FeHCurl,3,1,Double>(),
                                               1.0, coef, updatedGeo_);
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
    bool coefUpdateGeo = true;
    // ==================
    //  FLUX DENSITY
    // ==================
    LOG_DBG(magEdgePde) << "Reading prescribed flux density";

    ReadRhsExcitation( "fluxDensity", vecDofNames, ResultInfo::VECTOR, isComplex_, 
                       ent, coef, coefUpdateGeo );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST ||
          ent[i]->GetType() == EntityList::SURF_ELEM_LIST ) {
        EXCEPTION("Prescribed magnetic flux density can only be defined im volume")
      }
      Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
      PtrCoefFct factor = CoefFunction::Generate(mp_, part,
                                                 CoefXprBinOp(mp_, reluc_, coef[i] , CoefXpr::OP_MULT ) );
      
//      if(isComplex_) {
//        lin = new BDUIntegrator<CurlOperator<FeHCurl,3, Double>, Complex>(Complex(1.0), coef[i],
//                                                                          reluc_);
//      } else {
//        lin = new BDUIntegrator<CurlOperator<FeHCurl,3, Double>, Double>(1.0, coef[i],
//                                                                         reluc_ );
//      }
      EntityIterator it = ent[i]->GetIterator();
      it.Begin();
      if(isComplex_) {
             lin = new BUIntegrator<Complex>( new CurlOperator<FeHCurl,3, Complex>(),
                                              Complex(1.0), factor, coefUpdateGeo);
           } else {
             lin = new BUIntegrator<Double>( new CurlOperator<FeHCurl,3, Double>(),
                                             1.0, factor, coefUpdateGeo);
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
	// Use complete implicit scheme
    Double gamma = 1.0;
    GLMScheme * scheme = new Trapezoidal(gamma);
    TimeSchemeGLM::NonLinType nlType = (nonLin_)? TimeSchemeGLM::INCREMENTAL : TimeSchemeGLM::NONE;
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(scheme, 0, nlType) );

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
    shared_ptr<ResultInfo> fluxDens(new ResultInfo);
    fluxDens->resultType = MAG_FLUX_DENSITY;
    fluxDens->dofNames = vecComponents;
    fluxDens->unit = "Vs/m^2";
    fluxDens->definedOn = ResultInfo::ELEMENT;
    fluxDens->entryType = ResultInfo::VECTOR;
    availResults_.insert( fluxDens );
    shared_ptr<CoefFunctionFormBased> bFunc;
    if( isComplex_ ) {
      bFunc.reset(new CoefFunctionBOp<Complex>(feFct, fluxDens));
    } else {
      bFunc.reset(new CoefFunctionBOp<Double>(feFct, fluxDens));
    }
    DefineFieldResult( bFunc, fluxDens );
    stiffFormCoefs_.insert(bFunc);

    // === MAGNETIC NORMAL FLUX DENSITY ===
    shared_ptr<ResultInfo> normFlux(new ResultInfo);
    normFlux->resultType = MAG_NORMAL_FLUX_DENSITY;
    normFlux->dofNames = "";
    normFlux->unit = "Vs/m^2";
    normFlux->entryType = ResultInfo::SCALAR;
    normFlux->definedOn = ResultInfo::ELEMENT;
    shared_ptr<CoefFunctionSurf> sNormFDens;
    sNormFDens.reset(new CoefFunctionSurf(true, normFlux));
    DefineFieldResult( sNormFDens, normFlux );
    surfCoefFcts_[sNormFDens] = bFunc;

    // === MAGNETIC_FLUX ===
    shared_ptr<ResultInfo> flux(new ResultInfo);
    flux->resultType = MAG_FLUX;
    flux->dofNames = "";
    flux->unit = "Vs";
    flux->entryType = ResultInfo::SCALAR;
    flux->definedOn = ResultInfo::SURF_REGION;
    shared_ptr<ResultFunctor> fluxFct;
    if( isComplex_ ) {
      fluxFct.reset(new ResultFunctorIntegrate<Complex>(sNormFDens,
                                                          feFct, flux ) );
    } else {
      fluxFct.reset(new ResultFunctorIntegrate<Double>(sNormFDens,
                                                         feFct, flux ) );
    }
    resultFunctors_[MAG_FLUX] = fluxFct;
    availResults_.insert(flux);


    // === MAGNETIC RHS ===
    shared_ptr<ResultInfo> rhs(new ResultInfo);
    rhs->resultType = MAG_RHS_LOAD;
    rhs->dofNames = vecComponents;
    rhs->unit = "-";
    rhs->entryType = ResultInfo::VECTOR;
    rhs->definedOn = ResultInfo::ELEMENT;
    availResults_.insert( rhs );
    rhsFeFunctions_[MAG_POTENTIAL]->SetResultInfo(rhs);
    DefineFieldResult( rhsFeFunctions_[MAG_POTENTIAL], rhs );

    
    shared_ptr<CoefFunctionFormBased> jFunc;
    shared_ptr<CoefFunction> jPowerDensFunc;
    if( analysistype_ != STATIC ) {
      // === EDDY CURRENT DENSITY ===
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
      
      // === EDDY CURRENT (SURFACE RESULT) ===
      shared_ptr<ResultInfo> ec(new ResultInfo());
      ec->resultType = MAG_EDDY_CURRENT;
      ec->dofNames = "";
      ec->unit = "A";
      ec->definedOn = ResultInfo::SURF_REGION;
      ec->entryType = ResultInfo::SCALAR;
      availResults_.insert( ec );
      
      // first, create normal mapping
      shared_ptr<CoefFunctionSurf> ncd(new CoefFunctionSurf(true, ec));
      surfCoefFcts_[ncd] = jFunc;
          
      // then, integrate values
      shared_ptr<ResultFunctor> eddyCurrentFunc;
      if( isComplex_ ) {
        eddyCurrentFunc.reset(new ResultFunctorIntegrate<Complex>(ncd, 
                                                                  feFct, ec ) );
      } else {
        eddyCurrentFunc.reset(new ResultFunctorIntegrate<Double>(ncd, 
                                                           feFct, ec ) );
      }
      resultFunctors_[MAG_EDDY_CURRENT] = eddyCurrentFunc;
         
      // === EDDY POWER DENSITY ===
      shared_ptr<ResultInfo> epd(new ResultInfo());
      epd->resultType = MAG_EDDY_POWER_DENSITY;
      epd->dofNames = "";
      epd->unit = "W/m^3";
      epd->definedOn = ResultInfo::ELEMENT;
      epd->entryType = ResultInfo::SCALAR;
      shared_ptr<CoefFunctionFormBased> epdFunctor;
      if( isComplex_ ) { 
        epdFunctor.reset( new CoefFunctionBdBKernel<Complex>(aDotFct, 1.0));
      } else {
        epdFunctor.reset( new CoefFunctionBdBKernel<Double>(aDotFct, 1.0));
      }
      DefineFieldResult( epdFunctor, epd );
      massFormCoefs_.insert(epdFunctor);

      // === EDDY POWER ===
      shared_ptr<ResultInfo> ep(new ResultInfo());
      ep->resultType = MAG_EDDY_POWER;
      ep->dofNames = "";
      ep->unit = "W";
      ep->definedOn = ResultInfo::REGION;
      ep->entryType = ResultInfo::SCALAR;
      availResults_.insert( ep );
      shared_ptr<ResultFunctor> epFunctor;
      if( isComplex_ ) {
        epFunctor.reset(new EnergyResultFunctor<Complex>(aDotFct, ep, 1.0));
      } else {
        epFunctor.reset(new EnergyResultFunctor<Double>(aDotFct, ep, 1.0));
      }
      resultFunctors_[MAG_EDDY_POWER] = epFunctor;
      massFormFunctors_.insert(epFunctor);
    }
      
//      
//      // === EDDY POWER DENSITY ===
//      shared_ptr<ResultInfo> eddyPowerDens(new ResultInfo);
//      eddyPowerDens->resultType = MAG_EDDY_POWER_DENSITY;
//      eddyPowerDens->dofNames = "";
//      eddyPowerDens->unit = "W/m^3";
//      eddyPowerDens->definedOn = ResultInfo::ELEMENT;
//      eddyPowerDens->entryType = ResultInfo::SCALAR;
//      availResults_.insert( eddyPowerDens );
//      
//      jPowerDensFunc = CoefFunction::Generate(mp_, part, 
//                           CoefXprBinOp(mp_,
//                              CoefXprBinOp(mp_, jFunc, jFunc, CoefXpr::OP_MULT),
//                              conduc_, CoefXpr::OP_DIV) );
//      DefineFieldResult( jPowerDensFunc, eddyPowerDens);
//    }

    // === COIL CURRENT DENSITY ===
    shared_ptr<ResultInfo> ccd(new ResultInfo);
    ccd->resultType = MAG_COIL_CURRENT_DENSITY;
    ccd->dofNames = vecComponents;
    ccd->unit = "A/m^2";
    ccd->definedOn = ResultInfo::ELEMENT;
    ccd->entryType = ResultInfo::VECTOR;
    availResults_.insert( ccd );
    shared_ptr<CoefFunctionMulti> ccdCoef(new CoefFunctionMulti(CoefFunction::VECTOR,dim_,1, 
                                                                isComplex_));
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
    shared_ptr<CoefFunctionMulti> tcdCoef(new CoefFunctionMulti(CoefFunction::VECTOR, dim_,1,
                                                                isComplex_));
    DefineFieldResult( tcdCoef, tcd );


    // === LORENTZ FORCE DENSITY ===
    shared_ptr<ResultInfo> lfd(new ResultInfo);
    lfd->resultType = MAG_FORCE_LORENTZ_DENSITY;
    lfd->dofNames = vecComponents;
    lfd->unit = "N/m^3";
    lfd->definedOn = ResultInfo::ELEMENT;
    lfd->entryType = ResultInfo::VECTOR;
    availResults_.insert( lfd );

    // assemble coefficient function F_L = J X B
    PtrCoefFct lfdFunc = CoefFunction::Generate( mp_, part, 
                                                 CoefXprBinOp(mp_,  tcdCoef, bFunc, CoefXpr::OP_CROSS ) );
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

    // === PERMEABILITY  ===
    shared_ptr<ResultInfo> perm(new ResultInfo);
    perm->resultType = MAG_ELEM_PERMEABILITY;
    perm->dofNames = "";
    perm->unit = "Vs/Am";
    perm->definedOn = ResultInfo::ELEMENT;
    perm->entryType = ResultInfo::SCALAR;
    availResults_.insert( perm );

    // assemble coefficient function mu = 1 / nu
    Double oneOverMu0 = 1.0 / (4*PI*1e-7);
    PtrCoefFct muFunc = CoefFunction::Generate( mp_, part, 
                                                CoefXprBinOp( mp_, lexical_cast<std::string>(oneOverMu0), reluc_, CoefXpr::OP_DIV ) );
    DefineFieldResult( muFunc, perm);


    // === MAGNETIC FIELD INTENSITY ===
    shared_ptr<ResultInfo> magIntens(new ResultInfo);
    magIntens->resultType = MAG_FIELD_INTENSITY;
    magIntens->dofNames = vecComponents;
    magIntens->unit = "A/m";
    magIntens->definedOn = ResultInfo::ELEMENT;
    magIntens->entryType = ResultInfo::VECTOR;
    availResults_.insert( magIntens );

    // assemble coefficient function hFunc = reluctivity * flux density
    PtrCoefFct hFunc = CoefFunction::Generate( mp_, part, 
                                                CoefXprBinOp( mp_, reluc_, bFunc, CoefXpr::OP_MULT ) );
    DefineFieldResult( hFunc, magIntens);
    
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
      energyFunc.reset(new EnergyResultFunctor<Complex>(feFct, energy, 0.5));
    } else {
      energyFunc.reset(new EnergyResultFunctor<Double>(feFct, energy, 0.5));
    }
    resultFunctors_[MAG_ENERGY] = energyFunc;
    stiffFormFunctors_.insert(energyFunc);
  }
  
  void MagEdgePDE::FinalizePostProcResults() {

    // Initialize standard postprocessing results
    SinglePDE::FinalizePostProcResults();

    // Initialized results involving coils

    // === EDDY CURRENT DENSITY ===
    shared_ptr<CoefFunctionFormBased> jEddyCoef =
            dynamic_pointer_cast<CoefFunctionFormBased>(fieldCoefs_[MAG_EDDY_CURRENT_DENSITY]);
    std::map<RegionIdType, BaseBDBInt*>::iterator massIt = massInts_.begin();
    for( ; massIt != massInts_.end(); ++massIt ) {
      RegionIdType region = massIt->first;
      BaseBDBInt* massInt = massIt->second;
      
      // only assign region to jEddy, if it not a "regularized" region, i.e. 
      // only regions with "physical" conductivity get assigned
      if( regularizedRegions_.find(region) == regularizedRegions_.end()) {
        jEddyCoef->AddIntegrator( massInt, region);
      }
    }
    
    // === COIL CURRENT DENSITY ===
    shared_ptr<CoefFunctionMulti> ccdCoef =
        dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_COIL_CURRENT_DENSITY]);
    // loop over all coil coefficients and add contribution to coef 
    std::map<RegionIdType, PtrCoefFct>::iterator coilIt = coilCoefs_.begin();
    for( ; coilIt != coilCoefs_.end(); ++coilIt ) {
      ccdCoef->AddRegion( coilIt->first, coilIt->second);
    }

    // === TOTAL CURRENT DENSITY ===
    PtrCoefFct jEddy = GetCoefFct(MAG_EDDY_CURRENT_DENSITY);
    shared_ptr<CoefFunctionMulti> tcdCoef = 
        dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_TOTAL_CURRENT_DENSITY]);
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
          tcdCoef->AddRegion( actRegion, jEddy );
        }
      }
    }
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

} // end of namespace

