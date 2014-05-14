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
#include "Forms/BiLinForms/BiLinWrappredLinForm.hh"
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

    // initialize material coef functions covering all regions
    reluc_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, dim_, dim_, isComplex_));
    conduc_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1, isComplex_));

    // determine if there are coils excited by voltage
    hasVoltCoils_ = false;
    PtrParamNode coilNode = myParam_->Get( "coilList", ParamNode::PASS );
    if ( coilNode ){
      ParamNodeList coilNodes = coilNode->GetChildren();
      for( UInt k = 0; k < coilNodes.GetSize(); k++ ){
        if( coilNodes[k]->Has("source") ){
          std::string exType = coilNodes[k]->Get("source")->Get("type")->As<std::string>();
          if( exType == "voltage" ){
            hasVoltCoils_ = true;
            break;
          }
        }
      }
    }

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
    shared_ptr<BaseFeFunction> feFunc = feFunctions_[MAG_POTENTIAL];
    shared_ptr<FeSpace> feSpace = feFunc->GetFeSpace();
    
    
    for(UInt iRegion = 0; iRegion < regions_.GetSize() ; iRegion ++){
      actRegion = regions_[iRegion];
      actMat    = materials_[actRegion];
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);



      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();
      feSpace->SetRegionApproximation(actRegion,polyId,integId);
      
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
         stiff2->SetNewtonBilinearForm();

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
         
         // Set explicitly the solution dependency
         rhsNlinForm->SetSolDependent();
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
//      if(useGrad) {
//        dynamic_pointer_cast<FeSpaceHCurlHi>(feSpace)->SetUseGradients(actRegion);
//      }
      
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

    } // end for regions
    
    // ============================
    // COIL INTEGRATORS
    // ============================
    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
    std::map<Coil::IdType, shared_ptr<Coil> >::iterator it;
    it = coils_.begin();
    for( ; it != coils_.end(); it++ ) {
      Coil& actCoil = *(it->second);
      // run over all parts
      std::map<RegionIdType,shared_ptr<Coil::Part> >::iterator partIt;
      partIt = actCoil.parts_.begin();
      if( actCoil.sourceType_ == Coil::CURRENT ) {
        for( partIt = actCoil.parts_.begin(); 
            partIt != actCoil.parts_.end(); 
            partIt++ ) {
          Coil::Part & actPart = *(partIt->second);
          RegionIdType actRegion = partIt->first;
          shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
          actSDList->SetRegion( actRegion );

          // generate source current vector
          CoefXprVecScalOp iVec = CoefXprVecScalOp(mp_, actPart.jUnitVec, actCoil.srcVal_,
                                                   CoefXpr::OP_MULT);
          PtrCoefFct iFct = CoefFunction::Generate(mp_, part, iVec);

          CoefXprVecScalOp jVec = CoefXprVecScalOp(mp_, iFct, boost::lexical_cast<std::string>(actPart.wireCrossSect), 
                                                   CoefXpr::OP_DIV);
          PtrCoefFct jFct = CoefFunction::Generate(mp_, part, jVec);
          coilCurrentDens_[actRegion] = jFct;
          LinearForm* curInt;
          if( isComplex_ ) {
            curInt = new BUIntegrator<Complex>( new IdentityOperator<FeHCurl,3,1,Complex>(),
                                                1.0, jFct, updatedGeo_);
          }
          else {
            curInt = new BUIntegrator<Double>( new IdentityOperator<FeHCurl,3,1,Double>(),
                                               1.0, jFct, updatedGeo_);
          }
          curInt->SetName("CoilIntegrator");
          LinearFormContext * coilContext =
              new LinearFormContext( curInt );
          coilContext->SetEntities( actSDList );
          coilContext->SetFeFunction( feFunc );
          assemble_->AddLinearForm( coilContext );
          // obtain coefficient function
        } // loop: parts
      } else {
        EXCEPTION(" Implement voltage driven coils .. see the following code from NACS");
        
        
        // Note: Use the new BiLinWrappedLinForm to make a BiLinearForm out of 
        //       the LinearForms f_A and (f_A)^T in the following
   
        
//        // ===========================================
//        // 2) VOLTAGE driven coils
//        //
//        // Ref: M. Kaltenbacher, Numer. Sim. of. Mech.
//        //      Sens. and Act, 2nd edition, p. 211ff
//        // ===========================================
//        LOG_DBG(magpdebase) << "=> Type is VOLTAGE";
//
//        // The coupled equation system in this case looks like
//        //
//        //    ( M_A     0 ) ( A_dot ) + ( K_A -f_A  ) ( A ) = ( 0 )
//        //    ( (f_A)^T 0 ) ( i_dot )   ( 0     R   ) ( i )   ( u )
//
//
//        // Define new result for coil current (just once)
//        shared_ptr<ResultInfo> res = results_[resultIndex_[COIL_CURRENT]];
//
//        // Define entityList for coil
//        shared_ptr<CoilList> actCoilList (new CoilList(ptgrid_) );
//        actCoilList->AddCoil( it->second ); 
//
//        // calculate accumulated resistance
//        std::string totalR = "0.0 ";
//
//        // Run over single parts
//        for( partIt = actCoil.parts_.begin(); 
//            partIt != actCoil.parts_.end(); partIt++ ) {
//          Coil::Part & actPart = *(partIt->second);
//          RegionIdType actRegion = partIt->first;
//          shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
//          actSDList->SetRegion( actRegion );
//
//          // append resistance of this coil part
//          totalR += " + " + actPart.resistance;
//
//          LOG_DBG(magpdebase) << "Treating coil part on region '" << 
//              ptgrid_->RegionIdToName(actRegion) << "'";
//
//          // Define symmetric coupling integrators -f_A and and (f_A)^T with 
//          // 
//          //     f_A = 1 / wireCrossSection * int ... dOmega
//
//          // === Coupling integrator K_12 = -f_A ===
//          {
//            MechVolForceInt * cplInt = GetRHSInt( jDim, actCoil.phase_ );
//            StdVector<std::string> flowDir(1);
//            flowDir[0] = boost::lexical_cast<std::string>(actPart.orientFlag) 
//                                        + "* (-1.0)";
//            // check for solid / stranded conductor model
//            if( actPart.uniformCurrentDens_) {
//              flowDir[0] += " /" +  
//                  boost::lexical_cast<std::string>(actPart.wireCrossSect);
//            }
//            LOG_DBG3(magpdebase)  << "Vector for integrator -f_A: " << flowDir; 
//            cplInt->SetVolForceVector( flowDir, actPart.coordSys.get(),
//                                       true, 0.0 );
//
//            BiLinFormContext * cplContext = new BiLinFormContext( cplInt, STIFFNESS );
//            cplContext->SetCounterPart( false );
//            cplContext->SetPtPdes(this, this);
//            cplContext->SetResults( results_[0],
//                                    results_[resultIndex_[COIL_CURRENT]],
//                                    actSDList, actCoilList );
//            assemble_->AddBiLinearForm( cplContext );
//          }
//          // === Coupling integrator M_21 = (f_A)^T ===
//          {
//            MechVolForceInt * cplInt = GetRHSInt( jDim, actCoil.phase_ );
//            StdVector<std::string> flowDir(1);
//            flowDir[0] = boost::lexical_cast<std::string>(actPart.orientFlag);
//            // check for solid / stranded conductor model
//            if( actPart.uniformCurrentDens_) {
//              flowDir[0] += " / " + 
//                  boost::lexical_cast<std::string>(actPart.wireCrossSect);
//            }
//            LOG_DBG3(magpdebase)  << "Vector for integrator (f_A)^T: " << flowDir;
//            cplInt->SetVolForceVector( flowDir, actPart.coordSys.get(),
//                                       true, 0.0 );
//            cplInt->SetTransposed();
//            BiLinFormContext * cplContext = new BiLinFormContext( cplInt, MASS );
//            cplContext->SetCounterPart( false );
//            cplContext->SetPtPdes(this, this);
//            cplContext->SetResults( results_[resultIndex_[COIL_CURRENT]],
//                                    results_[0],
//                                    actCoilList, actSDList  );
//            assemble_->AddBiLinearForm( cplContext );
//          }
//        } // end loop over parts
//
//        LOG_DBG(magpdebase) << "Finished treating parts";
//
//        // === Diagonal Integrator for coil resistance K_22 = R ===
//        LOG_DBG(magpdebase) << "Resistance value is " << totalR;
//        SingleEntryInt * resistInt = 
//            new SingleEntryInt( totalR, 1, 1 );
//        BiLinFormContext * resistContext = 
//            new BiLinFormContext( resistInt, STIFFNESS );
//        resistContext->SetPtPdes(this, this);
//        resistContext->SetResults( results_[resultIndex_[COIL_CURRENT]],
//                                   results_[resultIndex_[COIL_CURRENT]],
//                                   actCoilList, actCoilList );
//        assemble_->AddBiLinearForm( resistContext );
//
//
//
//        // === RHS integrator (single entry) for coil voltageRHS_2 =
//        shared_ptr<LoadBc> actLoad( new LoadBc );
//
//        actLoad->entities = actCoilList;
//        actLoad->result = results_[resultIndex_[COIL_CURRENT]];
//        actLoad->value = actCoil.value_;
//        actLoad->phase = actCoil.phase_;
//        actLoad->eqnMap = eqnMap_;
//        actLoad->dof = 1;
//        LoadList loadList;
//        loadList.Push_back( actLoad );
//        assemble_->AddLoads( loadList );
//        LOG_DBG(magpdebase) << "Voltage value is " << actLoad->value;
//
//        // Pass result on this coils to the eqnMap class
//        eqnMap_->AddResult( *res, actCoilList );

        
      } //if: current / voltage driven
    } // loop: coils


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
    
//    solveStep_ = new StdSolveStep(*this);
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
    PtrParamNode coilNode = myParam_->Get( "coilList", ParamNode::PASS );
    PtrParamNode coilInfoNode = myInfo_->Get( "coilList", ParamNode::PASS );
    if ( !coilNode )
      return;

    // Get single coil nodes
    ParamNodeList coilNodes = coilNode->GetChildren();

    // Trigger reading in of definitions
    if( coilNodes.GetSize() > 0 ) {
      for( UInt i = 0; i < coilNodes.GetSize(); i++ ) {
        Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
        
        // get coil and id
        std::string coilId = coilNodes[i]->Get("id")->As<std::string>();

        // Check if coil with same ID already exists
        if( coils_.find(coilId) != coils_.end() ) {
          EXCEPTION("A coil with ID '" << coilId << "' was already defined.")
        }

        // Create new coil
        shared_ptr<Coil> actCoil( new Coil( coilNodes[i], coilInfoNode, 
                                            ptGrid_, mp_, part) );
        coils_[coilId] = actCoil;

        // Associate mapping of coil parts with regions
        std::map<RegionIdType, shared_ptr<Coil::Part> >::const_iterator it;
        for( it = actCoil->parts_.begin(); it != actCoil->parts_.end(); it++ ) {
          coilRegions_[it->first] = actCoil;
        }
      }

      // Adjust printing of coil information to info node
      WARN("Adapt printing of coils to InfoNode");
    }
  }


  void MagEdgePDE::DefinePrimaryResults() {

    StdVector<std::string> vecComponents;
    vecComponents = "x", "y", "z";

    // === MAGNETIC VECTOR POTENTIAL ===
    shared_ptr<ResultInfo> potInfo(new ResultInfo);
    potInfo->resultType = MAG_POTENTIAL;
    potInfo->dofNames = vecComponents;
    potInfo->unit = "Vs/m";
    potInfo->definedOn = ResultInfo::ELEMENT;
    potInfo->entryType = ResultInfo::VECTOR;

    feFunctions_[MAG_POTENTIAL]->SetResultInfo(potInfo);
    DefineFieldResult( feFunctions_[MAG_POTENTIAL], potInfo );

    // === COIL CURRENT ===
    if( hasVoltCoils_ ){
      shared_ptr<ResultInfo> currentInfo(new ResultInfo);
      currentInfo->resultType = COIL_CURRENT;
      currentInfo->dofNames = "";
      currentInfo->unit = "A";
      currentInfo->definedOn = ResultInfo::COIL;
      currentInfo->entryType = ResultInfo::SCALAR;

      feFunctions_[COIL_CURRENT]->SetResultInfo(currentInfo);
      DefineFieldResult( feFunctions_[COIL_CURRENT], currentInfo );
    }

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
    rhsFeFunctions_[MAG_POTENTIAL]->SetResultInfo(rhs);
    DefineFieldResult( rhsFeFunctions_[MAG_POTENTIAL], rhs );


    // === RESULTS RELATED TO EDDY CURRENTS ===
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
    std::map<RegionIdType, PtrCoefFct>::iterator coilIt = coilCurrentDens_.begin();
    for( ; coilIt != coilCurrentDens_.end(); ++coilIt ) {
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

    // assemble coefficient function mu = 1 / nu
    Double oneOverMu0 = 1.0 / (4e-7*PI);
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
    std::map<RegionIdType, PtrCoefFct>::iterator coilIt = coilCurrentDens_.begin();
    for( ; coilIt != coilCurrentDens_.end(); ++coilIt ) {
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
      if( coilCurrentDens_.find(actRegion) != coilCurrentDens_.end() ) {
        // region is a coil
        tcdCoef->AddRegion( actRegion, coilCurrentDens_[actRegion] );
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

    // ===================================================================
    // Create FeSpaceConst for coil currents (of coils excited by voltage)
    // ===================================================================
    if( hasVoltCoils_ ){
      PtrParamNode voltSpaceNode = infoNode->Get("coilCurrent");
      crSpaces[COIL_CURRENT] =
          FeSpace::CreateInstance(myParam_, voltSpaceNode, FeSpace::CONSTANT, ptGrid_);
      crSpaces[COIL_CURRENT]->Init(solStrat_);
    }

    return crSpaces;
  }

} // end of namespace

