#include <fstream>

#include "MagEdgeAdjPDE.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Utils/Coil.hh"
#include "Utils/SmoothSpline.hh"
#include "Utils/LinInterpolate.hh"

#include "Driver/Assemble.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"
#include "FeBasis/HCurl/FeSpaceHCurlHi.hh"
#include "FeBasis/HCurl/HCurlElems.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

#include "Domain/CoefFunction/CoefFunctionHarmBalance.hh"
#include "Domain/CoefFunction/CoefFunctionExpression.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefFunctionMulti.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/CoefFunction/CoefFunctionOpt.hh"



// forms
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/BiLinForms/BiLinWrappedLinForm.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/BDUInt.hh"
#include "Forms/LinForms/KXInt.hh"
#include "Forms/LinForms/SingleEntryInt.hh"
#include "Forms/Operators/CurlOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/ConvectiveOperator.hh"
#include "Forms/Operators/GradientOperator.hh"

//time stepping
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"

#include "Driver/MultiHarmonicDriver.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"
namespace CoupledField {

// declare class specific logging stream
DEFINE_LOG(magEdgeAdjPde, "magEdgeAdjPde")


  // **************
  //  Constructor
  // **************
  MagEdgeAdjPDE::MagEdgeAdjPDE( Grid * aptgrid, PtrParamNode paramNode,
                          PtrParamNode infoNode,
                          shared_ptr<SimState> simState, Domain* domain )
    :MagBasePDE( aptgrid, paramNode, infoNode, simState, domain ) {

    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_          = "magneticEdgeAdj";
    pdematerialclass_ = ELECTROMAGNETIC;
    formulation_ = MagBasePDE::EDGE;

    //! Always use updated Lagrangian formulation
    updatedGeo_        = true; //true;

    // default is false
    useGradFields_ = paramNode->Get("useGradientFields")->As<bool>();
    
    // check if we have a 3d setup
    //    bool is3d = domain_->GetParamRoot()->Get("domain")->Get("geometryType")->As<std::string>() == "3d";
    if ( !is3d_ )
      EXCEPTION("MagEdgeAdjPDE is just implemented for 3D setups!");

    nuDerivParam_ = false;
  }


  // *************
  //  Destructor
  // *************
  MagEdgeAdjPDE::~MagEdgeAdjPDE() {
  }

  // *****************************
  //  Definition of Integrators
  // *****************************
  void MagEdgeAdjPDE::DefineIntegrators() {
    this->DefineStandardIntegrators();
    // in MagBasePDE
    //DefineCoilIntegrators(1.0);
  }


  void MagEdgeAdjPDE::DefineStandardIntegrators(){
    RegionIdType actRegion;
    BaseMaterial * actMat = NULL;

    // initially, check for regularization factor
    Double regularizationFactor = 1e-6;
    myParam_->GetValue("penaltyFactor", regularizationFactor, ParamNode::PASS);

    shared_ptr<BaseFeFunction> feFunc = feFunctions_[MAG_POTENTIAL_ADJ];
    shared_ptr<FeSpace> feSpace = feFunc->GetFeSpace();

    PtrCoefFct magFluxCoef = this->GetCoefFct(MAG_FLUX_DENSITY);
    for(UInt iRegion = 0; iRegion < regions_.GetSize() ; iRegion ++){
      actRegion = regions_[iRegion];
      actMat    = materials_[actRegion];
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);

      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();
      feSpace->SetRegionApproximation(actRegion,polyId,integId);

      if(useGradFields_){
        feSpace->SetUseGradients(actRegion);
      }

      //get possible nonlinearities defined in this region
      StdVector<NonLinType> matDepenTypes = regionMatDepTypes_[actRegion]; // material dependency
      StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion];   // non-linearity

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( actRegion );

      // pass entitylist ot fespace / fefunction
      feFunc->AddEntityList( actSDList );


      PtrCoefFct curCoef;
      BaseBDBInt* curlcurl = NULL;
      // ===============================
      //  Standard Stiffness Integrator
      // ===============================
      curCoef = actMat->GetScalCoefFnc(MAG_RELUCTIVITY_SCALAR, Global::REAL);      
      
      //compute permeability
      PtrCoefFct constOne = CoefFunction::Generate( mp_, Global::REAL, "1.0");
      PtrCoefFct permeability = NULL;
      PtrCoefFct nuDeriv = NULL;

      //get reluctivity as coef-function
      if ( matDepenTypes.Find(RELUCTIVITY) != -1 ) {
        PtrCoefFct coef; 
        //dependency of magnetic reluctivity of magnetic flux density (computed by forward problem)            
        StdVector<std::string> dispDofNames;  
        shared_ptr<EntityList> ent = ptGrid_->GetEntityList( EntityList::ELEM_LIST, regionName );

        //get coeff-Fnc for magnetic flux density
        ReadMaterialDependency( "reluctivity", dispDofNames, ResultInfo::VECTOR, isComplex_,
                                ent, coef, updatedGeo_ );
        //std::cout << "ReadMaterialDependency OK, now define coefFunction nu" << std::endl;
        //coef-Fnc for magnetic reluctivity
        curCoef = actMat->GetScalCoefFncNonLin( MAG_RELUCTIVITY_SCALAR, Global::REAL, coef);
        //std::cout << "coef-Fnc for magnetic reluctivity OK" << std::endl;      

        //second part according to newton scheme (tangenial part)     
        nuDeriv = actMat->GetTensorCoefFncNonLin( MAG_RELUCTIVITY_DERIV, FULL, Global::REAL, coef );

        //coeff-Function for derivations of nu w.r.t. parameters
        PtrCoefFct derivParam1 = actMat->GetScalCoefFncNonLinDerivParam(MAG_RELUCTIVITY_DERIV_P1, Global::REAL, coef);
        PtrCoefFct derivParam2 = actMat->GetScalCoefFncNonLinDerivParam(MAG_RELUCTIVITY_DERIV_P2, Global::REAL, coef);
        PtrCoefFct derivParam3 = actMat->GetScalCoefFncNonLinDerivParam(MAG_RELUCTIVITY_DERIV_P3, Global::REAL, coef);
        PtrCoefFct derivParam4 = actMat->GetScalCoefFncNonLinDerivParam(MAG_RELUCTIVITY_DERIV_P4, Global::REAL, coef);

        nuDerivParamM1_[actRegion] = derivParam1;
        nuDerivParamM2_[actRegion] = derivParam2;
        nuDerivParamM3_[actRegion] = derivParam3;
        nuDerivParamM4_[actRegion] = derivParam4;
        bPostproc_[actRegion] = coef;
      } 
      else {
        curCoef = actMat->GetScalCoefFnc(MAG_RELUCTIVITY_SCALAR, Global::REAL ); 
      }

      permeability = CoefFunction::Generate( mp_,  Global::REAL, CoefXprBinOp(mp_, constOne, curCoef, CoefXpr::OP_DIV ) );
      matCoefs_[MAG_ELEM_PERMEABILITY]->AddRegion(actRegion, permeability);

      curlcurl = new BBInt<>(new  CurlOperator<FeHCurl,3, Double>(), curCoef,1.0, updatedGeo_) ;
      reluc_->AddRegion(actRegion, curCoef);

      curlcurl->SetName("CurlCurlIntegrator");
      BiLinFormContext * stiffContext = new BiLinFormContext(curlcurl, STIFFNESS );
      stiffContext->SetEntities( actSDList, actSDList );
      stiffContext->SetFeFunctions( feFunc, feFunc );
      assemble_->AddBiLinearForm( stiffContext );

      // Important: Add bdb-integrator to global list, as we need them later
      // for calculation of postprocessing results
      bdbInts_.insert( std::pair<RegionIdType, BaseBDBInt*>(actRegion,curlcurl) );

      if ( nuDeriv != NULL ) {
        //create tangential stiffness integrator according to Newton scheme
        BiLinearForm* stiff2 = NULL;
        stiff2 = new BDBInt<>(new CurlOperator<FeHCurl,3, Double>(), nuDeriv, 1.0, updatedGeo_) ;
        stiff2->SetName("CurlCurlIntegrator-NL-Newton");
        stiff2->SetNewtonBilinearForm();

        BiLinFormContext * stiffContext2 = new BiLinFormContext(stiff2, STIFFNESS );
        stiffContext2->SetEntities( actSDList, actSDList );
        stiffContext2->SetFeFunctions( feFunc, feFunc );
         assemble_->AddBiLinearForm( stiffContext2 );        
      }

      // ============================================================
      // Mass Matrix
      // ============================================================
      PtrCoefFct conductivityCoeff;
      Double conductivity = 0.0;
      if ( materials_[actRegion]->GetSymmetryType(MAG_CONDUCTIVITY_TENSOR) == BaseMaterial::GENERAL ){
        conductivityCoeff = materials_[actRegion]->GetTensorCoefFnc(MAG_CONDUCTIVITY_TENSOR, FULL, Global::REAL);
        Matrix<Double> conduc;
        materials_[actRegion]->GetTensor( conduc, MAG_CONDUCTIVITY_TENSOR, Global::REAL );
        conductivity = ( conduc[0][0]+conduc[1][1]+conduc[2][2] )/3;
      }else{
        materials_[actRegion]->GetScalar(conductivity,MAG_CONDUCTIVITY_SCALAR,Global::REAL);
        conductivityCoeff = CoefFunction::Generate(mp_, Global::REAL,lexical_cast<std::string>(conductivity));
      }
      // regularize
      bool scaleByEdgeSize = false;
      if ( conductivity < 1e-10 || analysistype_ == STATIC ) {
        Matrix<Double> reluc;
        // get tensor of permeability and determine max. value
        materials_[actRegion]->GetTensor( reluc, MAG_RELUCTIVITY_TENSOR, Global::REAL );
        conductivityCoeff = CoefFunction::Generate(mp_, Global::REAL,lexical_cast<std::string>(regularizationFactor * reluc[0][0]));
        scaleByEdgeSize = true;
        regularizedRegions_.insert(actRegion);
      }

      conduc_->AddRegion(actRegion, conductivityCoeff);
      BaseBDBInt *massInt;
      BiLinFormContext * massContext;
      if ( analysistype_ == STATIC) {
        // we have to guarantee, that we add some mass to curl-curl integrator.
        // Additionally, the integrator gets scaled by the edge size for a uniform
        // conditioning
        massInt = new BBIntMassEdge<>(new ScaledByEdgeIdentityOperator<FeHCurl,3,Double>(),
            conductivityCoeff,1.0);
        massInt->SetName("MassIntegrator");
        massContext =  new BiLinFormContext(massInt, STIFFNESS );
      } else {
        // here we add the "normal" mass integrator, which gets not scaled by the
        // edge size
        if( scaleByEdgeSize ) { // this is for non-conducting regions
          massInt = new BBIntMassEdge<>(new ScaledByEdgeIdentityOperator<FeHCurl,3,Double>(),
              conductivityCoeff,1.0, updatedGeo_);
          massContext = new BiLinFormContext(massInt, STIFFNESS );
        } else {
          if ( conductivityCoeff->GetDimType() == CoefFunction::TENSOR ) {
            massInt = new BDBInt<>(new IdentityOperator<FeHCurl,3,1,Double>(), conductivityCoeff, 1.0, updatedGeo_);
          } else {
            massInt = new BBIntMassEdge<>(new IdentityOperator<FeHCurl,3,1,Double>(), conductivityCoeff, 1.0, updatedGeo_);
          }
          massContext = new BiLinFormContext(massInt, DAMPING );
        }
      }
      
      massInt->SetName("MassIntegrator");
      massContext->SetEntities( actSDList, actSDList );
      massContext->SetFeFunctions( feFunc, feFunc );
      assemble_->AddBiLinearForm( massContext );
      // insert mass integrator to list of defined mass integrators
      massInts_[actRegion] = massInt;      
    }

    //check if reluctivity  depends on parameters
    if ( nuDerivParamM1_.size() > 0 )
      nuDerivParam_ = true;

  }


  BaseBDBInt* MagEdgeAdjPDE::GeHystStiffInt( Double factor, PtrCoefFct tensorReluctivity ){
    BaseBDBInt* stiffInt = NULL;

    // -------------------
    //  EDGE FORMULATION
    // -------------------
    // ===  3D CASE ===
    if( dim_ == 3 ) {
      stiffInt = new BDBInt<>(new CurlOperator<FeHCurl,3, Double>(), tensorReluctivity, factor, updatedGeo_) ;
    } else {
      // ===  2D / AXI CASE ===
      EXCEPTION("MagEdgeAdjPDE is just implemented for 3D setups!");
    }

    return stiffInt;
  }
  
  void MagEdgeAdjPDE::DefineNcIntegrators() {
    StdVector< NcInterfaceInfo >::iterator ncIt = ncInterfaces_.Begin(), endIt = ncInterfaces_.End();
    for ( ; ncIt != endIt; ++ncIt ) {
      switch (ncIt->type) {
        case NC_MORTAR:
          EXCEPTION("No Mortar nonconforming interface for magnetic PDE with edge elements.\n"
                    "Try using H1 nodal elements in MagneticPDE")
          break;
        case NC_NITSCHE:
          if (dim_ == 2)
            EXCEPTION("MagEdgeAdjPDE only works for 3D geometry!")
          else
            if(analysistype_ == MULTIHARMONIC){
              DefineNitscheCoupling<3,1>(MAG_POTENTIAL_ADJ, *ncIt, reluc_ );
            }else{
              DefineNitscheCoupling<3,1>(MAG_POTENTIAL_ADJ, *ncIt );
            }
          break;
        default:
          EXCEPTION("Unknown type of ncInterface");
          break;
      }
    }
  }


  void MagEdgeAdjPDE::DefineRhsLoadIntegrators() {
    LOG_TRACE(magEdgeAdjPde) << "Defining rhs load integrators for magEdgeAdjPDE";

    // Get FESpace and FeFunction 
    shared_ptr<BaseFeFunction> myFct = feFunctions_[MAG_POTENTIAL_ADJ];

    StdVector<shared_ptr<EntityList> > ent;
    StdVector<PtrCoefFct > coef;
    LinearForm * lin1 = NULL;    
    LinearForm * lin2 = NULL;  
    bool coefUpdateGeo = true;

    StdVector<std::string> dofNames(dim_);
    dofNames[0] = "Bx";
    dofNames[1] = "By";
    if ( dim_ == 3)
      dofNames[2] = "Bz";

    //read the measured flux density in the sensor positions and define RHS integrator
    ReadRhsExcitation("measuredFluxDensity", dofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo);

    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST ||
          ent[i]->GetType() == EntityList::SURF_ELEM_LIST ) {
        EXCEPTION("Measured magnetic flux density can only be specified in a volume!")
      }
      if(isComplex_) {
        lin1 = new BUIntegrator<Complex>( new CurlOperator<FeHCurl,3, Complex>(),
                                          Complex(1.0), coef[i], coefUpdateGeo);
      } else {
        lin1 = new BUIntegrator<Double>( new CurlOperator<FeHCurl,3, Double>(),
                                         1.0, coef[i], coefUpdateGeo);
      }
      lin1->SetName("MeasuredMagFluxIntensityInt");
      lin1->SetNormalizeToVol();
      LinearFormContext *ctx = new LinearFormContext( lin1 );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);

      bRHSRegions_[ent[i]->GetRegion()] = coef[i];      
    }

    // ==================
    //  FLUX DENSITY
    // ==================
    LOG_DBG(magEdgeAdjPde) << "Reading magnetic flux density from forward computation";

    //read magnetic field intensity from forward computation 
    //get coeff-Fnc for magnetic flux density
    StdVector<std::string> nameOfDofs;
    ReadRhsExcitation( "fluxDensity", nameOfDofs, ResultInfo::VECTOR, isComplex_,
                       ent, coef, coefUpdateGeo );

    //Please note: we have do adapt the vector B and set all components to zero., which 
    //where not measured by the sensors; e.g., when the sensor just measures the x-component 
    //then we have to set the y- and z-component to zero, which we do by the scalVec!
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST ||
          ent[i]->GetType() == EntityList::SURF_ELEM_LIST ) {
        EXCEPTION("Prescribed magnetic flux density can only be specified in a volume!")
      }

      Vector<Double> scalVec(dim_);
      scalVec.Init(0.0);
      std::istringstream iss(coef[i]->GetDofNames());
      do {
        std::string sub;
        iss >> sub;
        //sub = "all";
        if( sub == "all" ) {
          // add all dofs to the definedDofs
          for( UInt k = 0; k<dim_; k++ ) {
            scalVec[k] = 1.0;
          }
          break;
        } 
        else {
            if ( sub == "Bx" ) 
              scalVec[0] = 1.0;
            else if ( sub == "By")
              scalVec[1] = 1.0;
            else if ( sub == "Bz" && dim_ == 3 )
              scalVec[2] = 1.0;
        }
      } while (iss);

      //read magnetic field intensity from forward computation 
      //get coeff-Fnc for magnetic flux density
      //ReadMaterialDependency( "fluxDensity", dofNames, ResultInfo::VECTOR, isComplex_,
      //                        ent[i], coef[i], updatedGeo_ );

      // Here we store the B-field of the previous (forward) simulation
      Bmap_[ent[i]->GetRegion()] = coef[i];

      if(isComplex_) {
        lin2 = new BUIntegrator<Complex>( new CurlOperator<FeHCurl,3, Complex>(),
                                         Complex(-1.0), coef[i], coefUpdateGeo);
      } else {
        lin2 = new BUIntegrator<Double>( new CurlOperator<FeHCurl,3, Double>(),
                                        -1.0, coef[i], coefUpdateGeo);
      }
      lin2->SetName("MagFlux4ForwardIntegrator");
      lin2->SetScalVec(scalVec);
      //std::cout << "ScalVec: " << scalVec << std::endl;
      lin2->SetNormalizeToVol();
      LinearFormContext *ctx = new LinearFormContext( lin2 );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);

      bRHSRegions_[ent[i]->GetRegion()] = coef[i];
    } // for

  }

  // ======================================================
  // TIME-STEPPING SECTION > MagBasePDE
  // ======================================================

  void MagEdgeAdjPDE::DefinePrimaryResults() {

    StdVector<std::string> vecComponents;
    vecComponents = "x", "y", "z";

    // === MAGNETIC VECTOR POTENTIAL ===
    shared_ptr<ResultInfo> potInfo(new ResultInfo);
    potInfo->resultType = MAG_POTENTIAL_ADJ;
    potInfo->dofNames = vecComponents;
    potInfo->unit = "Vs/m";
    potInfo->definedOn = ResultInfo::ELEMENT;
    potInfo->entryType = ResultInfo::VECTOR;

    feFunctions_[MAG_POTENTIAL_ADJ]->SetResultInfo(potInfo);
    DefineFieldResult( feFunctions_[MAG_POTENTIAL_ADJ], potInfo );

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    hdbcSolNameMap_[MAG_POTENTIAL_ADJ] = "fluxParallel";
    idbcSolNameMap_[MAG_POTENTIAL_ADJ] = "potential";

    // // === PERMEABILITY ===
    shared_ptr<ResultInfo> permeability ( new ResultInfo );
    permeability->resultType = MAG_ELEM_PERMEABILITY;
    permeability->dofNames = "";
    permeability->unit = "Vs/Am";
    permeability->definedOn = ResultInfo::ELEMENT;
    permeability->entryType = ResultInfo::SCALAR;

    if(analysistype_ == MULTIHARMONIC){
      shared_ptr<CoefFunctionMulti> permFct(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1, true));
      matCoefs_[MAG_ELEM_PERMEABILITY] = permFct;
      DefineFieldResult(permFct, permeability);
    }else{
      //In multiharmonic analysis we have complex permeability
      shared_ptr<CoefFunctionMulti> permFct(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1, false));
      matCoefs_[MAG_ELEM_PERMEABILITY] = permFct;
      DefineFieldResult(permFct, permeability);
    }
  }

  void MagEdgeAdjPDE::DefinePostProcResults() {

    StdVector<std::string> vecComponents;
    vecComponents = "x", "y", "z";

    shared_ptr<BaseFeFunction> feFct = feFunctions_[MAG_POTENTIAL_ADJ];

    // === MAGNETIC RHS ===
    shared_ptr<ResultInfo> rhs(new ResultInfo);
    rhs->resultType = MAG_RHS_LOAD;
    rhs->dofNames = vecComponents;
    rhs->unit = "-";
    rhs->entryType = ResultInfo::VECTOR;
    rhs->definedOn = ResultInfo::ELEMENT;
    rhsFeFunctions_[MAG_POTENTIAL_ADJ]->SetResultInfo(rhs);
    DefineFieldResult( rhsFeFunctions_[MAG_POTENTIAL_ADJ], rhs );

    // === CURL OF ADJOINT SOLUTION ===
    shared_ptr<ResultInfo> curlAdj(new ResultInfo);
    curlAdj->resultType = MAG_CURL_ADJ;
    curlAdj->dofNames = vecComponents;
    curlAdj->unit = "Vs/m^2";
    curlAdj->definedOn = ResultInfo::ELEMENT;
    curlAdj->entryType = ResultInfo::VECTOR;
    availResults_.insert( curlAdj );
    shared_ptr<CoefFunctionFormBased> curlFunc;
    if( isComplex_ ) {
      curlFunc.reset(new CoefFunctionBOp<Complex>(feFct, curlAdj));
    } else {
      curlFunc.reset(new CoefFunctionBOp<Double>(feFct, curlAdj));
    }
    DefineFieldResult( curlFunc, curlAdj );
    stiffFormCoefs_.insert(curlFunc);

    // === Compute the volume
    shared_ptr<ResultInfo> vol(new ResultInfo);
    vol->resultType = VOLUME;
    vol->dofNames = "";
    vol->unit = "m^3";
    vol->entryType = ResultInfo::SCALAR;
    vol->definedOn = ResultInfo::REGION;
    shared_ptr<CoefFunction> coefVol = CoefFunction::Generate( mp_,  Global::REAL, "1.0");
    shared_ptr<ResultFunctor> volFunctor;
    volFunctor.reset(new ResultFunctorIntegrate<Double>(coefVol, feFct, vol));
    resultFunctors_[VOLUME] = volFunctor;
    availResults_.insert(vol);

    // === compute averaged magnetic flux density of forward simulation
    shared_ptr<ResultInfo> averagedB(new ResultInfo);
    averagedB->resultType = MAG_FLUX_DENSITY;
    averagedB->dofNames = "";
    averagedB->unit = "T";
    averagedB->entryType = ResultInfo::VECTOR;
    averagedB->definedOn = ResultInfo::ELEMENT;
    // The computation is defined in FinalizePostProcResults()
    shared_ptr<CoefFunctionMulti> averagedBfnc(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_));
    DefineFieldResult(averagedBfnc, averagedB);  

    //functor for averaged b-field of forward simulation
    shared_ptr<ResultInfo> resB1(new ResultInfo());
    resB1->resultType = MAG_AVERAGED_FLUX_DENSITY;
    resB1->dofNames = vecComponents;
    resB1->unit = "T";
    resB1->entryType = ResultInfo::VECTOR;
    resB1->definedOn = ResultInfo::REGION;   
    availResults_.insert(resB1);                             
    shared_ptr<ResultFunctor> averagedBfield;
    averagedBfield.reset(new ResultFunctorIntegrate<Double>(averagedBfnc, feFct, resB1));
    dynamic_pointer_cast< ResultFunctorIntegrate<Double> >(averagedBfield)->SetAveraged(true);
    resultFunctors_[MAG_AVERAGED_FLUX_DENSITY] = averagedBfield;   


  }

  void MagEdgeAdjPDE::FinalizePostProcResults() {

    // Initialize standard postprocessing results
    SinglePDE::FinalizePostProcResults();
    StdVector<std::string> vecComponents;
    vecComponents = "x", "y", "z";

    shared_ptr<BaseFeFunction> feFct = feFunctions_[MAG_POTENTIAL_ADJ];

    if (nuDerivParam_ == false) {
      // === Pure helper result ===
      shared_ptr<ResultInfo> ef(new ResultInfo);
      ef->resultType = OPT_RESULT_1; //RHS_PSEUDO_DENSITY;
      ef->dofNames = "";
      ef->unit = "";
      ef->definedOn = ResultInfo::ELEMENT;
      ef->entryType = ResultInfo::SCALAR;
      // The computation of the gradient of the design parameters is defined in FinalizePostProcResults()
      shared_ptr<CoefFunctionMulti> mult(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));
      DefineFieldResult(mult, ef);

      // === Gradient via adjoint ===
      shared_ptr<ResultInfo> gradAdjParam;
      gradAdjParam.reset(new ResultInfo);
      gradAdjParam->resultType = MAG_GRAD_ADJ_PARAM;
      gradAdjParam->dofNames = "";
      gradAdjParam->unit = "";
      gradAdjParam->entryType = ResultInfo::SCALAR;
      gradAdjParam->definedOn = ResultInfo::REGION;    
      availResults_.insert(gradAdjParam);                            
      shared_ptr<ResultFunctor> gradRegion;
      gradRegion.reset(new ResultFunctorIntegrate<Double>(mult, feFct, gradAdjParam));
      resultFunctors_[MAG_GRAD_ADJ_PARAM] = gradRegion;
    } 
    else {
      //the reluctivity depends on parameters
      
      //for parameter 1
      // === Pure helper result ===
      shared_ptr<ResultInfo> ef1(new ResultInfo);
      ef1->resultType = OPT_RESULT_1; 
      ef1->dofNames = "";
      ef1->unit = "";
      ef1->definedOn = ResultInfo::ELEMENT;
      ef1->entryType = ResultInfo::SCALAR;
      // The computation of the gradient of the design parameters is defined in FinalizePostProcResults()
      shared_ptr<CoefFunctionMulti> mult1(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));
      DefineFieldResult(mult1, ef1);

      // === Gradient via adjoint ===
      shared_ptr<ResultInfo> gradAdjParam1;
      gradAdjParam1.reset(new ResultInfo);
      gradAdjParam1->resultType = MAG_GRAD_ADJ_PARAM1;
      gradAdjParam1->dofNames = "";
      gradAdjParam1->unit = "";
      gradAdjParam1->entryType = ResultInfo::SCALAR;
      gradAdjParam1->definedOn = ResultInfo::REGION;    
      availResults_.insert(gradAdjParam1);                            
      shared_ptr<ResultFunctor> grad1Region;
      grad1Region.reset(new ResultFunctorIntegrate<Double>(mult1, feFct, gradAdjParam1));
      resultFunctors_[MAG_GRAD_ADJ_PARAM1] = grad1Region;

      //for parameter 2
      // === Pure helper result ===
      shared_ptr<ResultInfo> ef2(new ResultInfo);
      ef2->resultType = OPT_RESULT_2; 
      ef2->dofNames = "";
      ef2->unit = "";
      ef2->definedOn = ResultInfo::ELEMENT;
      ef2->entryType = ResultInfo::SCALAR;
      // The computation of the gradient of the design parameters is defined in FinalizePostProcResults()
      shared_ptr<CoefFunctionMulti> mult2(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));
      DefineFieldResult(mult2, ef2);

      // === Gradient via adjoint ===
      shared_ptr<ResultInfo> gradAdjParam2;
      gradAdjParam2.reset(new ResultInfo);
      gradAdjParam2->resultType = MAG_GRAD_ADJ_PARAM2;
      gradAdjParam2->dofNames = "";
      gradAdjParam2->unit = "";
      gradAdjParam2->entryType = ResultInfo::SCALAR;
      gradAdjParam2->definedOn = ResultInfo::REGION;    
      availResults_.insert(gradAdjParam2);                            
      shared_ptr<ResultFunctor> grad2Region;
      grad2Region.reset(new ResultFunctorIntegrate<Double>(mult2, feFct, gradAdjParam2));
      resultFunctors_[MAG_GRAD_ADJ_PARAM2] = grad2Region;

      //for parameter 3
      // === Pure helper result ===
      shared_ptr<ResultInfo> ef3(new ResultInfo);
      ef3->resultType = OPT_RESULT_3; 
      ef3->dofNames = "";
      ef3->unit = "";
      ef3->definedOn = ResultInfo::ELEMENT;
      ef3->entryType = ResultInfo::SCALAR;
      // The computation of the gradient of the design parameters is defined in FinalizePostProcResults()
      shared_ptr<CoefFunctionMulti> mult3(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));
      DefineFieldResult(mult3, ef3);

      // === Gradient via adjoint ===
      shared_ptr<ResultInfo> gradAdjParam3;
      gradAdjParam3.reset(new ResultInfo);
      gradAdjParam3->resultType = MAG_GRAD_ADJ_PARAM3;
      gradAdjParam3->dofNames = "";
      gradAdjParam3->unit = "";
      gradAdjParam3->entryType = ResultInfo::SCALAR;
      gradAdjParam3->definedOn = ResultInfo::REGION;    
      availResults_.insert(gradAdjParam3);                            
      shared_ptr<ResultFunctor> grad3Region;
      grad3Region.reset(new ResultFunctorIntegrate<Double>(mult3, feFct, gradAdjParam3));
      resultFunctors_[MAG_GRAD_ADJ_PARAM3] = grad3Region;

      //for parameter 4
      // === Pure helper result ===
      shared_ptr<ResultInfo> ef4(new ResultInfo);
      ef4->resultType = OPT_RESULT_4; 
      ef4->dofNames = "";
      ef4->unit = "";
      ef4->definedOn = ResultInfo::ELEMENT;
      ef4->entryType = ResultInfo::SCALAR;
      // The computation of the gradient of the design parameters is defined in FinalizePostProcResults()
      shared_ptr<CoefFunctionMulti> mult4(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));
      DefineFieldResult(mult4, ef4);

      // === Gradient via adjoint ===
      shared_ptr<ResultInfo> gradAdjParam4;
      gradAdjParam4.reset(new ResultInfo);
      gradAdjParam4->resultType = MAG_GRAD_ADJ_PARAM4;
      gradAdjParam4->dofNames = "";
      gradAdjParam4->unit = "";
      gradAdjParam4->entryType = ResultInfo::SCALAR;
      gradAdjParam4->definedOn = ResultInfo::REGION;    
      availResults_.insert(gradAdjParam4);                            
      shared_ptr<ResultFunctor> grad4Region;
      grad4Region.reset(new ResultFunctorIntegrate<Double>(mult4, feFct, gradAdjParam4));
      resultFunctors_[MAG_GRAD_ADJ_PARAM4] = grad4Region;
    }

    // ============ MAG_GRAD_ADJ_PARAM & avaeraged flux density from forward simulation ============
    // define functor to compute gradient part: 

    if ( nuDerivParam_ == false ) {
      shared_ptr<CoefFunctionMulti> scalMult = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[OPT_RESULT_1]);
      shared_ptr<CoefFunctionMulti> bField = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_FLUX_DENSITY]);
    
      StdVector<RegionIdType>::iterator regIt = regions_.Begin();
      regIt = regions_.Begin();
      for (; regIt != regions_.End(); ++regIt) {
        // ========= B field =============
        PtrCoefFct Bforward = NULL;
        PtrCoefFct mult = NULL;
        RegionIdType actRegion = *regIt;
        if(Bmap_.find(*regIt) != Bmap_.end()){
          // B from previous (forward) simulation
          Bforward = Bmap_[*regIt];
          // result MAG_CURL_ADJ is curl(p), where p is the adjoint solution!
          PtrCoefFct curlAdj = this->GetCoefFct(MAG_CURL_ADJ);
          mult = CoefFunction::Generate( mp_, Global::REAL,
                             CoefXprBinOp(mp_, Bforward, curlAdj, CoefXpr::OP_MULT) );        
          scalMult->AddRegion(actRegion, mult);

          // B from previous (forward) simulation
          bField->AddRegion(*regIt, Bforward);
        }
      } // loop over regions

    } 
    else {
      //the reluctivity depends on parameters
      shared_ptr<CoefFunctionMulti> scalMult1 = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[OPT_RESULT_1]);
      shared_ptr<CoefFunctionMulti> scalMult2 = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[OPT_RESULT_2]);
      shared_ptr<CoefFunctionMulti> scalMult3 = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[OPT_RESULT_3]);
      shared_ptr<CoefFunctionMulti> scalMult4 = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[OPT_RESULT_4]);            
      shared_ptr<CoefFunctionMulti> bField = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_FLUX_DENSITY]);
    
      StdVector<RegionIdType>::iterator regIt = regions_.Begin();
      regIt = regions_.Begin();
      for (; regIt != regions_.End(); ++regIt) {
        // ========= B field =============
        PtrCoefFct Bforward = NULL;
        PtrCoefFct multForAdj = NULL;
        PtrCoefFct multP1 = NULL;
        PtrCoefFct multP2 = NULL;
        PtrCoefFct multP3 = NULL;        
        PtrCoefFct multP4 = NULL;   
        PtrCoefFct relucParam1;
        PtrCoefFct relucParam2;        
        PtrCoefFct relucParam3;
        PtrCoefFct relucParam4;        
        RegionIdType actRegion = *regIt;

        if(Bmap_.find(*regIt) != Bmap_.end()){
          // B from previous (forward) simulation
          Bforward = Bmap_[*regIt];
          // B from previous (forward) simulation
          bField->AddRegion(*regIt, Bforward);
        }

        if (bPostproc_.find(*regIt) != bPostproc_.end()) {
          // B from previous (forward) simulation
          Bforward = bPostproc_[*regIt];
          bField->AddRegion(*regIt, Bforward);

          // result MAG_CURL_ADJ is curl(p), where p is the adjoint solution!
          PtrCoefFct curlAdj = this->GetCoefFct(MAG_CURL_ADJ);
          
          //multiply curl(Aforward) \cdot curl(adjoint)
          multForAdj = CoefFunction::Generate( mp_, Global::REAL,
                                  CoefXprBinOp(mp_, Bforward, curlAdj, CoefXpr::OP_MULT) );  

          //=======================
          //do all for parameter 1
          //=======================
          // get coeff function for derivation of reluctivity  w.r.t. parameter 1
          relucParam1 = nuDerivParamM1_[*regIt];      
          multP1 = CoefFunction::Generate( mp_, Global::REAL,
                             CoefXprBinOp(mp_, relucParam1, multForAdj, CoefXpr::OP_MULT) );                                 
          scalMult1->AddRegion(actRegion, multP1);

          //=======================
          //do all for parameter 2
          //=======================
          // get coeff function for derivation of reluctivity  w.r.t. parameter 1
          relucParam2 = nuDerivParamM2_[*regIt];       
          multP2 = CoefFunction::Generate( mp_, Global::REAL,
                             CoefXprBinOp(mp_, relucParam2, multForAdj, CoefXpr::OP_MULT) );                                 
          scalMult2->AddRegion(actRegion, multP2);          

          //=======================
          //do all for parameter 3
          //=======================
          // get coeff function for derivation of reluctivity  w.r.t. parameter 1
          relucParam3 = nuDerivParamM3_[*regIt];       
          multP3 = CoefFunction::Generate( mp_, Global::REAL,
                             CoefXprBinOp(mp_, relucParam3, multForAdj, CoefXpr::OP_MULT) );                                 
          scalMult3->AddRegion(actRegion, multP3);            

          //=======================
          //do all for parameter 4
          //=======================
          // get coeff function for derivation of reluctivity  w.r.t. parameter 1
          relucParam4 = nuDerivParamM4_[*regIt];       
          multP4 = CoefFunction::Generate( mp_, Global::REAL,
                             CoefXprBinOp(mp_, relucParam4, multForAdj, CoefXpr::OP_MULT) );                                 
          scalMult4->AddRegion(actRegion, multP4);            
        }
      } // loop over regions
    } 
  }

  std::map<SolutionType, shared_ptr<FeSpace> >
  MagEdgeAdjPDE::CreateFeSpaces(const std::string& formulation,
                             PtrParamNode infoNode ) {
    //ok default case so we create grid based approximation H1 elements
    //and standard Gauss integration
    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    if(formulation == "default" || formulation == "H_CURL"){
      PtrParamNode potSpaceNode = infoNode->Get("magPotentialAdj");
      crSpaces[MAG_POTENTIAL_ADJ] =
          FeSpace::CreateInstance(myParam_, potSpaceNode, FeSpace::HCURL, ptGrid_ );
      crSpaces[MAG_POTENTIAL_ADJ]->Init(solStrat_);
    }else{
      EXCEPTION("The formulation " << formulation
                << "of magnetic edge PDE is not known!");
    }

    // in addition query, if special treatment of anisotropic elements
    // is activated
    if( myParam_->Has("thinElements") ) {
      Double aspectRatio = 0.0;
      aspectRatio = myParam_->Get("thinElements")->Get("maxAspectRatio")->As<Double>();
      dynamic_pointer_cast<FeSpaceHCurlHi>(crSpaces[MAG_POTENTIAL_ADJ])
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
