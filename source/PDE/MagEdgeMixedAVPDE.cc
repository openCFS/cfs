#include <fstream>

#include "MagEdgeMixedAVPDE.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Utils/Coil.hh"
#include "Utils/SmoothSpline.hh"
#include "Utils/LinInterpolate.hh"

#include "Driver/Assemble.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"
#include "FeBasis/HCurl/FeSpaceHCurlHi.hh"
#include "FeBasis/HCurl/HCurlElems.hh"
#include "FeBasis/H1/FeSpaceH1Hi.hh"
#include "FeBasis/H1/H1Elems.hh"

#include "DataInOut/Logging/LogConfigurator.hh"

#include "Domain/CoefFunction/CoefFunctionExpression.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefFunctionMulti.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/CoefFunction/CoefXpr.hh"

// forms
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/BiLinForms/BiLinWrappedLinForm.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/BDUInt.hh"
#include "Forms/LinForms/KXInt.hh"
#include "Forms/LinForms/SingleEntryInt.hh"
#include "Forms/Operators/CurlOperator.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/DivOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"

//time stepping
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"


// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"
namespace CoupledField {

// declare class specific logging stream
DECLARE_LOG(magEdgeMixedAVPde)
DEFINE_LOG(magEdgeMixedAVPde, "magEdgeMixedAVPde")


  // **************
  //  Constructor
  // **************
  MagEdgeMixedAVPDE::MagEdgeMixedAVPDE( Grid * aptgrid, PtrParamNode paramNode,
                          PtrParamNode infoNode,
                          shared_ptr<SimState> simState, Domain* domain )
    :SinglePDE( aptgrid, paramNode, infoNode, simState, domain ) {

    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_          = "magneticEdge";
    pdematerialclass_ = ELECTROMAGNETIC;

    //! Always use updated Lagrangian formulation
    updatedGeo_        = true; //true;

    // check if we have a 3d setup
    bool is3d = domain_->GetParamRoot()->Get("domain")->Get("geometryType")->As<std::string>() == "3d";
    if ( !is3d )
      EXCEPTION("MagEdgeMixedAVPDE is just implemented for 3D setups!");

    // initialize material coef functions covering all regions
    reluc_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, dim_, dim_, isComplex_));
    conduc_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));

    // determine if there are coils excited by voltage
    hasVoltCoils_ = false;
    PtrParamNode coilNode = myParam_->Get( "coilList", ParamNode::PASS );



    if ( coilNode ){
      ParamNodeList coilNodes = coilNode->GetChildren();
      // there should only be one mixed flag and on formulation for all coils
      mixedCoil_ = coilNodes[0]->Get( "mixed" )->As<bool>();
      coilFormulation_ = coilNodes[0]->Get( "formulation" )->As<std::string>();
      for( UInt k = 0; k < coilNodes.GetSize(); k++ ){
        if( coilNodes[k]->Get( "formulation" )->As<std::string>() != coilFormulation_ ){
          EXCEPTION("The formulation of the coils must be the same for all coils!");
        }
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
  MagEdgeMixedAVPDE::~MagEdgeMixedAVPDE() {
  }

  shared_ptr<Coil> MagEdgeMixedAVPDE::GetCoilById(const Coil::IdType& id) {
    return coils_.at(id);
  }

  void MagEdgeMixedAVPDE::ReadSpecialBCs() {


    // --------------------------------------------------------------------
    //   Get information about coils and open files for measurement coils
    // --------------------------------------------------------------------
    ReadCoils();

  }


  // ****************************
  //  Initialize Nonlinearities
  // ****************************
  void MagEdgeMixedAVPDE::InitNonLin() {

    SinglePDE::InitNonLin();
  }


  // *****************************
  //  Definition of Integrators
  // *****************************
  void MagEdgeMixedAVPDE::DefineIntegrators() {

    RegionIdType actRegion;
    BaseMaterial * actMat = NULL;

    // initially, check for regularization factor
    Double regularizationFactor = 1e-6;
    myParam_->GetValue("penaltyFactor", regularizationFactor, ParamNode::PASS);

    shared_ptr<BaseFeFunction> magVecPotFeFunc = feFunctions_[MAG_POTENTIAL];
    shared_ptr<BaseFeFunction> elecScalPotFeFunc = feFunctions_[ELEC_POTENTIAL];
    shared_ptr<FeSpace> magVecPotFeSpace = magVecPotFeFunc->GetFeSpace();
    shared_ptr<FeSpace> elecScalPotFeSpace = elecScalPotFeFunc->GetFeSpace();


    PtrCoefFct magFluxCoef = this->GetCoefFct(MAG_FLUX_DENSITY);

    if ( (analysistype_ == HARMONIC) ) {
      EXCEPTION("MagEdgeMixedAVPDE currently only works for STATIC and TRANSIENT analysis")
    }

    for(UInt iRegion = 0; iRegion < regions_.GetSize() ; iRegion ++){
      actRegion = regions_[iRegion];
      actMat    = materials_[actRegion];
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());

      // Get polynomial and integration order for magnetic vector potential
      std::string magVecPolyId = curRegNode->Get("magVecPolyId")->As<std::string>();
      std::string magVecPolyIntegId = curRegNode->Get("magVecIntegId")->As<std::string>();
      magVecPotFeSpace->SetRegionApproximation(actRegion, magVecPolyId, magVecPolyIntegId);

      // Get polynomial and integration order for electric scalar potential
      std::string elecScalPolyId = curRegNode->Get("elecScalPolyId")->As<std::string>();
      std::string elecScalPolyIntegId = curRegNode->Get("elecScalPolyIntegId")->As<std::string>();
      elecScalPotFeSpace->SetRegionApproximation(actRegion, elecScalPolyId, elecScalPolyIntegId);

      // Get possible nonlinearities defined in this region
      StdVector<NonLinType> matDepenTypes = regionMatDepTypes_[actRegion]; // material dependency
      StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion];   // non-linearity

      // Create new entity list, based on the region
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( actRegion );

      // Pass entitylist to fespace / fefunction for magnetic vector and electric scalar potential
      magVecPotFeFunc->AddEntityList( actSDList );
      elecScalPotFeFunc->AddEntityList( actSDList );

      if(matDepenTypes.Find(NLELEC_CONDUCTIVITY) != -1){
        EXCEPTION("MagEdgeMixedAVPDE does not support nonlinear"
            "(temperatur dependent) electric conductivity yet!\n"
            "But the implementation is not big deal...I promise.\n"
            "Mostly c&p form MagEdgePDE.");
      }



      if ( nonLinTypes.GetSize() > 0 ){
        // =================================================================================
        //  NONLINEAR SECTION
        // =================================================================================
        EXCEPTION("MagEdgeMixedAVPDE does not support nonlinear reluctivity yet!";)
      } else {
        // =================================================================================
        //  LINEAR STIFFNESS SECTION
        // =================================================================================

        /* ==============================================
         * Handling of material parameters
           ============================================== */
        // Magnetic Reluctivity
        PtrCoefFct nuNl = NULL;
        nuNl = actMat->GetScalCoefFnc( MAG_RELUCTIVITY, Global::REAL);
        // Add material to global, distributed reluctivity coefficient function
        reluc_->AddRegion(actRegion, nuNl);

        // Magnetic Permeability
        PtrCoefFct constOne = CoefFunction::Generate( mp_, Global::REAL, "1.0");
        PtrCoefFct permeability = CoefFunction::Generate( mp_,  Global::REAL, CoefXprBinOp(mp_, constOne, nuNl, CoefXpr::OP_DIV ) );
        matCoefs_[MAG_ELEM_PERMEABILITY]->AddRegion(actRegion, permeability);

        // Electric Conductivity
        Double conductivity;
        materials_[actRegion]->GetScalar(conductivity,MAG_CONDUCTIVITY,Global::REAL);
        PtrCoefFct conducCoef = materials_[actRegion]->GetScalCoefFnc(MAG_CONDUCTIVITY,Global::REAL);


        /* ==============================================
         * Upper left STIFFNESS part:
         * curl(A) \cdot curl(A)
           ============================================== */
        BaseBDBInt* stiffUpperLeft = NULL;
        stiffUpperLeft = new BBInt<>(new  CurlOperator<FeHCurl,3, Double>(), nuNl, 1.0, updatedGeo_) ;
        stiffUpperLeft->SetName("CurlACurlAIntegratorUpperLeft");

        BiLinFormContext * stiffUpperLeftContext = new BiLinFormContext(stiffUpperLeft, STIFFNESS );
        stiffUpperLeftContext->SetEntities( actSDList, actSDList );
        stiffUpperLeftContext->SetFeFunctions( magVecPotFeFunc, magVecPotFeFunc );
        assemble_->AddBiLinearForm( stiffUpperLeftContext );


        /* ==============================================
         * Upper right STIFFNESS part:
         * \sigma grad(V) \cdot A
           ============================================== */
        BiLinearForm* stiffUpperRight = NULL;
        stiffUpperRight = new ABInt<>(new IdentityOperator<FeHCurl,3,1,Double>(),
                                      new GradientOperator<FeH1,3,1,Double>(),
                                      conducCoef, 1.0, updatedGeo_);
        stiffUpperRight->SetName("GradVIdentityAIntegratorUpperRight");

        BiLinFormContext * stiffUpperRightContext = new BiLinFormContext(stiffUpperRight, STIFFNESS );
        stiffUpperRightContext->SetEntities( actSDList, actSDList );
        stiffUpperRightContext->SetFeFunctions( magVecPotFeFunc, elecScalPotFeFunc );
        assemble_->AddBiLinearForm( stiffUpperRightContext );


        /* ==============================================
         * Lower right STIFFNESS part:
         * \sigma grad(V) \cdot grad(V)
           ============================================== */
        BiLinearForm* stiffLowerRight = NULL;
        stiffLowerRight = new BBInt<>(new  GradientOperator<FeH1,3,1,Double>(), conducCoef, 1.0, updatedGeo_) ;
        stiffLowerRight->SetName("GradVGradVIntegratorLowerRight");

        BiLinFormContext * stiffLowerRightContext = new BiLinFormContext(stiffLowerRight, STIFFNESS );
        stiffLowerRightContext->SetEntities( actSDList, actSDList );
        stiffLowerRightContext->SetFeFunctions( elecScalPotFeFunc, elecScalPotFeFunc );
        assemble_->AddBiLinearForm( stiffLowerRightContext );






        // =================================================================================
        //  LINEAR MASS SECTION
        // =================================================================================
        bool scaleByEdgeSize = false;
        if ( conductivity < 1e-10 || analysistype_ == STATIC ) {
          Matrix<Double> reluc;
          // Get tensor of permeability and determine max. value
          materials_[actRegion]->GetTensor( reluc, MAG_RELUCTIVITY, Global::REAL );
          conductivity =  regularizationFactor * reluc[0][0];
          scaleByEdgeSize = true;
          // Add region to set of "regularized" regions
          regularizedRegions_.insert(actRegion);
        }

        PtrCoefFct conducCoefReg = CoefFunction::Generate(mp_, Global::REAL, lexical_cast<std::string>(conductivity));
        // Also add material to global, distributed reluctivity coefficient function
        conduc_->AddRegion(actRegion, conducCoefReg);


        /* ==============================================
         * Upper left MASS part:
         * \sigma grad(V) \cdot grad(V)
           ============================================== */
        BaseBDBInt *massUpperLeftInt;
        BiLinFormContext * massUpperLeftContext;
        if ( analysistype_ == STATIC) {
          // We have to guarantee, that we add some mass to curl-curl integrator.
          // Additionally, the integrator gets scaled by the edge size for a uniform
          // conditioning
          massUpperLeftInt = new BBIntMassEdge<>(new ScaledByEdgeIdentityOperator<FeHCurl,3,Double>(), conducCoefReg,1.0);
          massUpperLeftInt->SetName("MassIntegratorUpperLeft");
          massUpperLeftContext =  new BiLinFormContext(massUpperLeftInt, STIFFNESS );
        } else {
          // here we add the "normal" mass integrator, which gets not scaled by the
          // edge size
          if( scaleByEdgeSize ) {
            massUpperLeftInt = new BBIntMassEdge<>(new ScaledByEdgeIdentityOperator<FeHCurl,3,Double>(), conducCoefReg,1.0, updatedGeo_);
            massUpperLeftContext = new BiLinFormContext(massUpperLeftInt, STIFFNESS );
          } else {
            massUpperLeftInt = new BBIntMassEdge<>(new IdentityOperator<FeHCurl,3,1,Double>(), conducCoefReg,1.0, updatedGeo_);
            massUpperLeftContext = new BiLinFormContext(massUpperLeftInt, DAMPING );
          }
          massUpperLeftInt->SetName("MassIntegratorUpperLeft");
        }
        massUpperLeftContext->SetEntities( actSDList, actSDList );
        massUpperLeftContext->SetFeFunctions( magVecPotFeFunc, magVecPotFeFunc );
        assemble_->AddBiLinearForm( massUpperLeftContext );



        /* ==============================================
         * Lower left MASS part:
         * \sigma div(A) V
           ============================================== */
        BaseBDBInt *massLowerLeftInt;
        BiLinFormContext * massLowerLeftContext;
        massLowerLeftInt = new ABInt<>(new DivOperator<FeHCurl,3,Double>(),
                           new IdentityOperator<FeH1,3,1,Double>(),
                           conducCoefReg, -1.0, updatedGeo_);
        massLowerLeftContext = new BiLinFormContext(massLowerLeftInt, DAMPING );
        massLowerLeftInt->SetName("MassIntegratorLowerLeft");

        massLowerLeftContext->SetEntities( actSDList, actSDList );
        massLowerLeftContext->SetFeFunctions( magVecPotFeFunc, elecScalPotFeFunc );
        assemble_->AddBiLinearForm( massLowerLeftContext );

     } // END OF NONLIN/LIN PART
    } // end for regions
  } // end DefineIntegrators



  void MagEdgeMixedAVPDE::DefineNcIntegrators() {
    EXCEPTION("Nonconforming interfaces not implemented for MagEdgeMixedAVPDE!!\n"
              "This will be a bit shitty");
  }


  void MagEdgeMixedAVPDE::DefineRhsLoadIntegrators() {
    LOG_TRACE(magEdgeMixedAVPde) << "Defining rhs load integrators for MagEdgeMixedAVPDE";

    shared_ptr<BaseFeFunction> magVecPotFeFunc = feFunctions_[MAG_POTENTIAL];
    shared_ptr<BaseFeFunction> elecScalPotFeFunc = feFunctions_[ELEC_POTENTIAL];


    StdVector<shared_ptr<EntityList> > ent;
    StdVector<PtrCoefFct > coef;
    LinearForm * upperRHSInt = NULL;

    StdVector<std::string> vecDofNames = magVecPotFeFunc->GetResultInfo()->dofNames;
    StdVector<std::string> scalDofNames = elecScalPotFeFunc->GetResultInfo()->dofNames;


    bool coefUpdateGeo = true;
    // ==================
    //  FLUX DENSITY
    // ==================
    LOG_DBG(magEdgeMixedAVPde) << "Reading prescribed flux density";

    ReadRhsExcitation( "fluxDensity", vecDofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST ||
          ent[i]->GetType() == EntityList::SURF_ELEM_LIST ) {
        EXCEPTION("Prescribed magnetic flux density can only be specified in a volume!")
      }
      Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
      PtrCoefFct factor = CoefFunction::Generate(mp_, part, CoefXprBinOp(mp_, reluc_, coef[i] , CoefXpr::OP_MULT ) );

      if(isComplex_) {
        upperRHSInt = new BUIntegrator<Complex>( new CurlOperator<FeHCurl,3, Complex>(), Complex(1.0), factor, coefUpdateGeo);
      } else {
        upperRHSInt = new BUIntegrator<Double>( new CurlOperator<FeHCurl,3, Double>(), 1.0, factor, coefUpdateGeo);
      }
      upperRHSInt->SetName("UpperRHSFluxInt");
      LinearFormContext *upperRHSContext = new LinearFormContext( upperRHSInt );
      upperRHSContext->SetEntities( ent[i] );
      upperRHSContext->SetFeFunction(magVecPotFeFunc);
      assemble_->AddLinearForm(upperRHSContext);
      magVecPotFeFunc->AddEntityList(ent[i]);

      bRHSRegions_[ent[i]->GetRegion()] = coef[i];
    } // for






/*



    // ==================
    //  FIELD INTENSITY
    // ==================
    LOG_DBG(magEdgeMixedAVPde) << "Reading prescribed field intensity";

    ReadRhsExcitation( "fieldIntensity", vecDofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST ||
          ent[i]->GetType() == EntityList::SURF_ELEM_LIST ) {
        EXCEPTION("Prescribed magnetic field intensity can only be specified in a volume!")
      }

      if(isComplex_) {
        lin = new BUIntegrator<Complex>( new CurlOperator<FeHCurl,3, Complex>(),
                                         Complex(1.0), coef[i], coefUpdateGeo);
      } else {
        lin = new BUIntegrator<Double>( new CurlOperator<FeHCurl,3, Double>(),
                                        1.0, coef[i], coefUpdateGeo);
      }
      lin->SetName("FieldIntensityIntegrator");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);

    } // for
*/
  }


  void MagEdgeMixedAVPDE::DefineSolveStep() {
    solveStep_ = new StdSolveStep(*this);
  }


  // ======================================================
  // TIME-STEPPING SECTION
  // ======================================================

  void MagEdgeMixedAVPDE::InitTimeStepping() {
	// Use complete implicit scheme
    Double gamma = 1.0;
    GLMScheme * scheme = new Trapezoidal(gamma);
    TimeSchemeGLM::NonLinType nlType = (nonLin_)? TimeSchemeGLM::INCREMENTAL : TimeSchemeGLM::NONE;
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(scheme, 0, nlType) );
    feFunctions_[MAG_POTENTIAL]->SetTimeScheme(myScheme);

    // Important: Create a new time scheme just for the elec potential unknowns, as otherwise the
    // size of the vectors does not match!
    GLMScheme * scheme2 = new Trapezoidal(gamma);
    shared_ptr<BaseTimeScheme> myScheme2(new TimeSchemeGLM(scheme2, 0, nlType) );
    feFunctions_[ELEC_POTENTIAL]->SetTimeScheme(myScheme2);
  }

  // ******************************************************
  //   Query parameter object for information about coils
  // ******************************************************
  void MagEdgeMixedAVPDE::ReadCoils() {
    PtrParamNode coilNode = myParam_->Get( "coilList", ParamNode::PASS );
    PtrParamNode coilInfoNode = myInfo_->Get( "coilList", ParamNode::PASS );
    if ( coilNode ){
      EXCEPTION("Currently no coils are supported for MagEdgeMixedAVPDE");
    }
  }


  void MagEdgeMixedAVPDE::DefinePrimaryResults() {

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

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    hdbcSolNameMap_[MAG_POTENTIAL] = "fluxParallel";
    idbcSolNameMap_[MAG_POTENTIAL] = "potential";


    // === ELECTRIC SCALAR POTENTIAL ===
    shared_ptr<BaseFeFunction> scalFct = feFunctions_[ELEC_POTENTIAL];
    shared_ptr<ResultInfo> res2(new ResultInfo);
    res2->resultType = ELEC_POTENTIAL;
    res2->dofNames = "";
    res2->unit = "V";
    res2->definedOn = ResultInfo::NODE;
    res2->entryType = ResultInfo::SCALAR;
    results_.Push_back( res2 );
    availResults_.insert( res2 );
    scalFct->SetResultInfo(res2);
    DefineFieldResult( scalFct, res2 );

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
//    hdbcSolNameMap_[ELEC_POTENTIAL] = "fluxParallel";
    idbcSolNameMap_[ELEC_POTENTIAL] = "elecPotential";


    // === PERMEABILITY ===
    shared_ptr<ResultInfo> permeability ( new ResultInfo );
    permeability->resultType = MAG_ELEM_PERMEABILITY;
    permeability->dofNames = "";
    permeability->unit = "Vs/Am";
    permeability->definedOn = ResultInfo::ELEMENT;
    permeability->entryType = ResultInfo::SCALAR;
    shared_ptr<CoefFunctionMulti> permFct(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1, false));
    matCoefs_[MAG_ELEM_PERMEABILITY] = permFct;
    DefineFieldResult(permFct, permeability);

  }

  void MagEdgeMixedAVPDE::DefinePostProcResults() {

    StdVector<std::string> vecComponents;
    vecComponents = "x", "y", "z";

//    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
    shared_ptr<BaseFeFunction> feFct = feFunctions_[MAG_POTENTIAL];

    // === TIME DERIVATIVES OF PRIMARY RESULTS ===
    if( analysistype_ == TRANSIENT || analysistype_ == HARMONIC ) {
      // === MAGNETIC VECTOR POTENTIAL - 1ST DERIVATIVE ===
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
    sNormFDens.reset(new CoefFunctionSurf(true, 1.0, normFlux));
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


  }

  void MagEdgeMixedAVPDE::FinalizePostProcResults() {

    // Initialize standard postprocessing results
    SinglePDE::FinalizePostProcResults();


  }

  std::map<SolutionType, shared_ptr<FeSpace> >
  MagEdgeMixedAVPDE::CreateFeSpaces(const std::string& formulation,
                             PtrParamNode infoNode ) {
    //ok default case so we create grid based approximation H1 elements
    //and standard Gauss integration
    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    PtrParamNode magVecPpotSpaceNode = infoNode->Get("magPotential");
    crSpaces[MAG_POTENTIAL] = FeSpace::CreateInstance(myParam_, magVecPpotSpaceNode, FeSpace::HCURL, ptGrid_ );
    crSpaces[MAG_POTENTIAL]->Init(solStrat_);

    PtrParamNode elecScalPotSpaceNode = infoNode->Get("elecPotential");
    crSpaces[ELEC_POTENTIAL] = FeSpace::CreateInstance(myParam_, elecScalPotSpaceNode, FeSpace::H1, ptGrid_);
    crSpaces[ELEC_POTENTIAL]->Init(solStrat_);


    return crSpaces;
  }

} // end of namespace

