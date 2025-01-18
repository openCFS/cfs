#include <fstream>
#include <string.h>
#include "FullWaveMaxwellEPDE.hh"

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
#include "Domain/CoefFunction/CoefFunctionPML.hh"

// forms
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/BDUInt.hh"
#include "Forms/LinForms/KXInt.hh"
#include "Forms/Operators/CurlOperator.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/DivOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"

// time stepping
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"
namespace CoupledField
{

  // declare class specific logging stream
  DEFINE_LOG(FullWaveMaxwellEPDE, "FullWaveMaxwellEPDE")

  // **************
  //  Constructor
  // **************
  FullWaveMaxwellEPDE::FullWaveMaxwellEPDE(Grid *aptgrid, PtrParamNode paramNode,
                                         PtrParamNode infoNode,
                                         shared_ptr<SimState> simState, Domain *domain)
      : SinglePDE(aptgrid, paramNode, infoNode, simState, domain)
  {

    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_ = "fullwave-E";
    // We reuse the Darwin material class (same as classical magnetic material plus permittivity)
    pdematerialclass_ = ELECTROMAGNETIC_DARWIN;

    //! Always use updated Lagrangian formulation
    updatedGeo_ = true; // true;

    formulationType_ = paramNode->Get("formulation")->As<std::string>();

    // check if we have a 3d setup
    bool is3d = domain_->GetParamRoot()->Get("domain")->Get("geometryType")->As<std::string>() == "3d";
    if (!is3d)
      EXCEPTION("FullWaveMaxwellEPDE is just implemented for 3D setups!");

    // initialize material coef functions covering all regions
    reluc_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, dim_, dim_, isComplex_));
    perm_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, dim_, dim_, isComplex_));
    conduc_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));
    eps_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));
  }

  // *************
  //  Destructor
  // *************
  FullWaveMaxwellEPDE::~FullWaveMaxwellEPDE()
  {
  }

  // ****************************
  //  Initialize Nonlinearities
  // ****************************
  void FullWaveMaxwellEPDE::InitNonLin()
  {

    SinglePDE::InitNonLin();
  }

  // ****************************
  //  Read damping information
  // ****************************
  void FullWaveMaxwellEPDE::ReadDampingInformation() {
    std::map<std::string, DampingType> idDampType;
    std::map<std::string, shared_ptr<RaylDampingData> > idRaylData;

    // try to get dampingList
    PtrParamNode dampListNode = myParam_->Get( "dampingList", ParamNode::PASS );
    if( dampListNode ) {

      // get specific damping nodes
      ParamNodeList dampNodes = dampListNode->GetChildren();

      for( UInt i = 0; i < dampNodes.GetSize(); i++ ) {

        std::string dampString = dampNodes[i]->GetName();
        std::string actId = dampNodes[i]->Get("id")->As<std::string>();

        // determine type of damping
        DampingType actType;
        String2Enum( dampString, actType );

        if( actType == RAYLEIGH ) {
          // set data for Rayleigh damping
          shared_ptr<RaylDampingData> actRaylDamp(new RaylDampingData());
          actRaylDamp->alpha = "0.0";
          actRaylDamp->beta = "0.0";
          actRaylDamp->adjustDamping = true;
          actRaylDamp->ratioDeltaF = 0.01;
          actRaylDamp->freq = 0.0;
          dampNodes[i]->GetValue( "freq", actRaylDamp->freq, ParamNode::PASS);
          dampNodes[i]->GetValue( "ratioDeltaF", actRaylDamp->ratioDeltaF,
                                  ParamNode::PASS );
          dampNodes[i]->GetValue( "adjustDamping",actRaylDamp->adjustDamping,  
                                  ParamNode::PASS);
          idRaylData[actId] = actRaylDamp;        
        }

        // store damping type string
        idDampType[actId] = actType;

      }
    }

    // Run over all region and set entry in "regionNonLinId"
    ParamNodeList regionNodes =
        myParam_->Get("regionList")->GetChildren();

    RegionIdType actRegionId;
    std::string actRegionName, actDampingId;

    for (UInt k = 0; k < regionNodes.GetSize(); k++) {
      regionNodes[k]->GetValue( "name", actRegionName );
      regionNodes[k]->GetValue( "dampingId", actDampingId );
      if( actDampingId == "" )
        continue;

      actRegionId = ptGrid_->GetRegion().Parse( actRegionName );

      // Check actDampingId was already registerd
      if( idDampType.count( actDampingId ) == 0 ) {
        EXCEPTION( "Damping with id '" << actDampingId
                   << "' was not defined in 'dampingList'" );
      }

      dampingList_[actRegionId] = idDampType[actDampingId];
      if ( dampingList_[actRegionId] == RAYLEIGH ){
        RaylDampingData actRayl = *(idRaylData[actDampingId]);
        Double dampFreq;

        if( actRayl.freq == 0.0 ) {
          materials_[actRegionId]->GetScalar(dampFreq,RAYLEIGH_FREQUENCY,Global::REAL);
        } else { 
          dampFreq = actRayl.freq;
        }
        // Compute Rayleigh damping parameters
        materials_[actRegionId]->
        ComputeRayleighDamping( actRayl.alpha, actRayl.beta,
                                dampFreq, actRayl.ratioDeltaF, 
                                actRayl.adjustDamping, isComplex_ );
        regionRaylDamping_[actRegionId] = actRayl;
      } else if(dampingList_[actRegionId] == PML &&
          analysistype_ == BasePDE::TRANSIENT ) {
        isTimeDomPML_ = true;
      }
      
    }
  }

  // *****************************
  //  Definition of Integrators
  // *****************************
  void FullWaveMaxwellEPDE::DefineIntegrators()
  {

    RegionIdType actRegion;
    BaseMaterial *actMat = NULL;

    // Double mu0 = 4*M_PI*1e-7;
    // Double eps0 = 8.854187817*1e-12;
    // Double k0_without_omega =  sqrt(mu0*eps0);
    // Double k0_without_omega_squared =  (mu0*eps0);
    // Double Z0 =  sqrt(mu0/eps0);
    // PtrCoefFct mu0coef = CoefFunction::Generate(mp_, Global::REAL, std::to_string(mu0));
    // PtrCoefFct eps0coef = CoefFunction::Generate(mp_, Global::REAL, std::to_string(eps0));



    shared_ptr<BaseFeFunction> eVecPotFeFunc = feFunctions_[ELEC_FIELD_INTENSITY];
    shared_ptr<FeSpace> eVecPotFeSpace = eVecPotFeFunc->GetFeSpace();

    for (UInt iRegion = 0; iRegion < regions_.GetSize(); iRegion++)
    {
      actRegion = regions_[iRegion];
      actMat = materials_[actRegion];
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region", "name", regionName.c_str());

      // Get polynomial and integration order for the E-field vector (primary unknown in the E-formulation)
      std::string eVecPolyId = curRegNode->Get("eVecPolyId")->As<std::string>();
      std::string eVecIntegId = curRegNode->Get("eVecIntegId")->As<std::string>();
      eVecPotFeSpace->SetRegionApproximation(actRegion, eVecPolyId, eVecIntegId);

      // Get possible nonlinearities defined in this region
      StdVector<NonLinType> matDepenTypes = regionMatDepTypes_[actRegion]; // material dependency
      StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion];   // non-linearity

      // Create new entity list, based on the region
      shared_ptr<ElemList> actSDList(new ElemList(ptGrid_));
      actSDList->SetRegion(actRegion);

      // Pass entitylist to fespace / fefunction for E-field (primary unknown in the E-formulation)
      eVecPotFeFunc->AddEntityList(actSDList);

      if (matDepenTypes.Find(NLELEC_CONDUCTIVITY) != -1)
      {
        EXCEPTION("FullWaveMaxwellEPDE does not support nonlinear"
                  "(temperatur dependent) electric conductivity yet!\n"
                  "But the implementation is not big deal...I promise.\n"
                  "Mostly c&p form MagEdgePDE.");
      }

      if (nonLinTypes.GetSize() > 0)
      {
        // =================================================================================
        //  NONLINEAR SECTION
        // =================================================================================
        EXCEPTION("FullWaveMaxwellEPDE does not support nonlinear reluctivity yet!";)
      }
      else
      {
        /* ==============================================
         * Handling of material parameters
           ============================================== */
        // Magnetic Reluctivity
        PtrCoefFct reluctivity = NULL;

        StdVector<NonLinType> matDepenTypes = regionMatDepTypes_[actRegion]; // material dependency
        if ( matDepenTypes.Find(MAT_RELUCTIVITY) != -1 ) {
          //we read the computed reluctvity from file
          shared_ptr<ResultInfo> resultInfo = GetResultInfo(MAG_ELEM_RELUCTIVITY); 
          std::string regionName = ptGrid_->GetRegion().ToString(actRegion);
          shared_ptr<EntityList> entity = ptGrid_->GetEntityList( EntityList::ELEM_LIST, regionName );
        
          //get coeff-Fnc for the magnetic permeability
          ReadMaterialDependency( "magElemReluctivity", resultInfo->dofNames, resultInfo->entryType, false,
                                    entity, reluctivity, updatedGeo_ );

        } else {
          reluctivity = actMat->GetScalCoefFnc(MAG_RELUCTIVITY_SCALAR, Global::REAL);
        }
        // Add material to global, distributed reluctivity coefficient function
        reluc_->AddRegion(actRegion, reluctivity);              
        matCoefs_[MAG_ELEM_RELUCTIVITY]->AddRegion(actRegion, reluctivity);

        // Magnetic Permeability
        PtrCoefFct constOne = CoefFunction::Generate(mp_, Global::REAL, "1.0");
        PtrCoefFct permeability = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, constOne, 
                                                         reluctivity, CoefXpr::OP_DIV));
        perm_->AddRegion(actRegion, permeability);  
        matCoefs_[MAG_ELEM_PERMEABILITY]->AddRegion(actRegion, permeability);


        // Electric / Magnetic Permittivity
        PtrCoefFct eps = NULL;        
        if ( matDepenTypes.Find(MAT_PERMITTIVITY) != -1 ) {
          //we read the computed permitivity from file
          shared_ptr<ResultInfo> resultInfo = GetResultInfo(MAG_ELEM_PERMITTIVITY); 
          std::string regionName = ptGrid_->GetRegion().ToString(actRegion);
          shared_ptr<EntityList> entity = ptGrid_->GetEntityList( EntityList::ELEM_LIST, regionName );
        
          //get coeff-Fnc for the magnetic permittivity
          ReadMaterialDependency( "magElemPermittivity", resultInfo->dofNames, resultInfo->entryType, false,
                                  entity, eps, updatedGeo_ );                                  
        } else {        
          eps = actMat->GetScalCoefFnc(MAG_PERMITTIVITY_SCALAR, Global::REAL);        
        }
        // Add material to global, distributed coefficient function
        eps_->AddRegion(actRegion, eps);                
        matCoefs_[MAG_ELEM_PERMITTIVITY]->AddRegion(actRegion, eps);

        // ====================================================================
        // Take account for pml (frequency domain only)
        // ====================================================================
        shared_ptr<CoefFunction> coeffPMLTensor, coeffPMLTensorInv;
        shared_ptr<CoefFunction> coeffPMLStiff;
        shared_ptr<CoefFunction> coeffPMLMass;
      
        bool harmonicPML = false;       
        if( dampingList_[actRegion] == PML ) {       
          //see book FINITE ELEMENT ANALYSIS OF ANTENNAS AND ARRAYS, Jian-Ming Jin, Douglas J. Riley
          std::string dampId;
          curRegNode->GetValue("dampingId",dampId);
          if (analysistype_ == HARMONIC) {
            // compute wave speed
            PtrCoefFct tmp = CoefFunction::Generate( mp_,  Global::REAL,
                                           CoefXprUnaryOp(mp_, CoefXprBinOp(mp_, permeability, eps,
                                           CoefXpr::OP_MULT), CoefXpr::OP_SQRT));                           
            PtrCoefFct c0  = CoefFunction::Generate( mp_, Global::REAL,
                                           CoefXprBinOp(mp_, constOne, tmp, CoefXpr::OP_DIV ) );

            PtrParamNode pmlNode = myParam_->Get("dampingList")->GetByVal("pml","id",dampId.c_str());            
            coeffPMLTensor.reset(new CoefFunctionPML<Complex>(pmlNode,c0,actSDList,regions_,
                                                              CoefFunction::CoefDimType::TENSOR));  
            coeffPMLTensor->SetPDEname(pdename_);  
            coeffPMLTensor->SetActive(false);                                                         
            coeffPMLTensorInv.reset(new CoefFunctionPML<Complex>(pmlNode,c0,actSDList,regions_,
                                                              CoefFunction::CoefDimType::TENSOR)); 
            coeffPMLTensorInv->SetPDEname(pdename_);  
            coeffPMLTensorInv->SetActive(true); //the inverse tensor will be applied!
          
            coeffPMLStiff  = CoefFunction::Generate( mp_, Global::COMPLEX,
                                                     CoefXprTensScalOp(mp_,coeffPMLTensorInv,reluctivity,
                                                    CoefXpr::OP_MULT));
            coeffPMLMass = CoefFunction::Generate( mp_,Global::COMPLEX,
                                                   CoefXprTensScalOp(mp_,coeffPMLTensor,eps,CoefXpr::OP_MULT));
            harmonicPML = true;
          }
        }
        else {
          harmonicPML = false;
        }       

        // =================================================================================
        // =================================================================================
        //  STIFFNESS SECTION
        // =================================================================================
        // =================================================================================

        /* ==============================================
         * K_E_E^(nu):
         * nu curl(E) \cdot curl(E)
           ============================================== */
        BaseBDBInt *K_E_E_nu = NULL;        
        if ( harmonicPML ) {
          K_E_E_nu = new BDBInt<Complex>(new CurlOperator<FeHCurl, 3, Complex>(), 
                                         coeffPMLStiff, 1.0, updatedGeo_);                               
        }
        else
          K_E_E_nu = new BBInt<>(new CurlOperator<FeHCurl, 3, Double>(), reluctivity, 1.0, updatedGeo_);

        K_E_E_nu->SetName("K_E_E_nu");
        BiLinFormContext *K_E_E_nuContext = new BiLinFormContext(K_E_E_nu, STIFFNESS);
        K_E_E_nuContext->SetEntities(actSDList, actSDList);
        K_E_E_nuContext->SetFeFunctions(eVecPotFeFunc, eVecPotFeFunc);
        assemble_->AddBiLinearForm(K_E_E_nuContext);
        // Add bdb-integrator to global list, needed for flux density evaluation
        bdbInts_.insert(std::pair<RegionIdType, BaseBDBInt *>(actRegion, K_E_E_nu));
        
        // =================================================================================
        // =================================================================================
        //  MASS SECTION
        // =================================================================================
        // =================================================================================
        BaseBDBInt *M_E_E_epsilon = NULL;       
         if ( harmonicPML ) {
          M_E_E_epsilon = new BDBInt<Complex>(new IdentityOperator<FeHCurl, 3, 1, Complex>(), 
                                              coeffPMLMass, 1.0, updatedGeo_);
         }
         else 
          M_E_E_epsilon = new BBInt<>(new IdentityOperator<FeHCurl, 3, 1, Double>(), eps, 
                                      1.0, updatedGeo_);
        M_E_E_epsilon->SetName("M_E_E_epsilon");

        BiLinFormContext *M_E_E_epsilonContext;
        M_E_E_epsilonContext = new BiLinFormContext(M_E_E_epsilon, MASS);
        M_E_E_epsilonContext->SetEntities(actSDList, actSDList);
        M_E_E_epsilonContext->SetFeFunctions(eVecPotFeFunc, eVecPotFeFunc);
        assemble_->AddBiLinearForm(M_E_E_epsilonContext);
      } // END OF NONLIN/LIN PART
    }   // end for regions
  }     // end DefineIntegrators

  void FullWaveMaxwellEPDE::DefineNcIntegrators()
  {
    StdVector<NcInterfaceInfo>::iterator ncIt = ncInterfaces_.Begin(), endIt = ncInterfaces_.End();
    for (; ncIt != endIt; ++ncIt)
    {
      switch (ncIt->type)
      {
      case NC_MORTAR:
        EXCEPTION("FullWaveMaxwellEPDE: Mortar NC interfaces not tested!");
      case NC_NITSCHE:
      {

        if (dim_ == 2){
            EXCEPTION("FullWaveMaxwellEPDE only works for 3D geometry!")
        } else {
            DefineNitscheCoupling<3,1>(ELEC_FIELD_INTENSITY, *ncIt );
        }
          break;
      }
      default:
        EXCEPTION("Unknown type of ncInterface");
        break;
      }
    }
  }

  void FullWaveMaxwellEPDE::DefineSurfaceIntegrators()
  {
  }

  void FullWaveMaxwellEPDE::DefineRhsLoadIntegrators()
  {
    //LOG_TRACE(FullWaveMaxwellEPDE) << "Defining rhs load integrators for FullWaveMaxwellEPDE";

    shared_ptr<BaseFeFunction> eVecPotFeFunc = feFunctions_[ELEC_FIELD_INTENSITY];

    StdVector<shared_ptr<EntityList>> ent;
    StdVector<PtrCoefFct> coef;
    LinearForm * lin = NULL;
    StdVector<std::string> vecDofNames = eVecPotFeFunc->GetResultInfo()->dofNames;
    
    bool coefUpdateGeo = true;
    // ==================
    //  IMPRESSED CURRENT DENSITY
    // ==================

    ReadRhsExcitation("impressedCurrentDensity", vecDofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo);
    for (UInt i = 0; i < ent.GetSize(); ++i)
    {
       // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST ||
          ent[i]->GetType() == EntityList::SURF_ELEM_LIST ) {
        EXCEPTION("Prescribed impressed current density can only be specified in a volume!")
      }
          
      if(isComplex_) {
        lin = new BUIntegrator<Complex>( new IdentityOperator<FeHCurl, 3, 1, Complex>(),
                                         Complex(-1.0), coef[i], coefUpdateGeo);
      } else {
        lin = new BUIntegrator<Double>( new IdentityOperator<FeHCurl, 3, 1, Double>(),
                                        1.0, coef[i], coefUpdateGeo);
      }

      lin->SetName("FluxIntegrator");
      
      // we need to multiply with j*omega and use the newly created possiblity to also
      // assign DAMPING, MASS, STIFFNESS to LinearForms (if nothing is defined, the default is STIFFNESS)
      lin->SetType(DAMPING);

      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(eVecPotFeFunc);
      assemble_->AddLinearForm(ctx);
      eVecPotFeFunc->AddEntityList(ent[i]);

      bRHSRegions_[ent[i]->GetRegion()] = coef[i];
    }
  }

  void FullWaveMaxwellEPDE::DefineSolveStep()
  {
    solveStep_ = new StdSolveStep(*this);
  }



  // ======================================================
  // TIME-STEPPING SECTION
  // ======================================================

  void FullWaveMaxwellEPDE::InitTimeStepping()
  {
    // Use complete implicit scheme
    Double gamma = 1.0;
    GLMScheme *scheme = new Trapezoidal(gamma);
    TimeSchemeGLM::NonLinType nlType = (nonLin_) ? TimeSchemeGLM::INCREMENTAL : TimeSchemeGLM::NONE;
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(scheme, 0, nlType));
    feFunctions_[ELEC_FIELD_INTENSITY]->SetTimeScheme(myScheme);
  }


  void FullWaveMaxwellEPDE::DefinePrimaryResults()
  {

    StdVector<std::string> vecComponents;
    vecComponents = "x", "y", "z";

    // === ELECTRIC FIELD INTENSITY ===
    shared_ptr<ResultInfo> elecIntens(new ResultInfo);
    elecIntens->resultType = ELEC_FIELD_INTENSITY;
    elecIntens->SetVectorDOFs(dim_, isaxi_);
    elecIntens->dofNames = vecComponents;
    elecIntens->unit = "V/m";
    elecIntens->definedOn = ResultInfo::ELEMENT;
    elecIntens->entryType = ResultInfo::VECTOR;
    feFunctions_[ELEC_FIELD_INTENSITY]->SetResultInfo(elecIntens);
    DefineFieldResult(feFunctions_[ELEC_FIELD_INTENSITY], elecIntens);

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    hdbcSolNameMap_[ELEC_FIELD_INTENSITY] = "fluxParallel";
    idbcSolNameMap_[ELEC_FIELD_INTENSITY] = "potential";

    // === PERMEABILITY ===
    shared_ptr<ResultInfo> permeability(new ResultInfo);
    permeability->resultType = MAG_ELEM_PERMEABILITY;
    permeability->dofNames = "";
    permeability->unit = "Vs/Am";
    permeability->definedOn = ResultInfo::ELEMENT;
    permeability->entryType = ResultInfo::SCALAR;
    shared_ptr<CoefFunctionMulti> permFct(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, false));
    matCoefs_[MAG_ELEM_PERMEABILITY] = permFct;
    DefineFieldResult(permFct, permeability);
  }

  void FullWaveMaxwellEPDE::DefinePostProcResults()
  {
    StdVector<std::string> vecComponents;
    vecComponents = "x", "y", "z";

    shared_ptr<BaseFeFunction> feFct = feFunctions_[ELEC_FIELD_INTENSITY];

    // === CURL OF ELECTRIC FIELD = ELECTRIC VORTICITY ===
    shared_ptr<ResultInfo> curlE(new ResultInfo);
    curlE->resultType = ELEC_FIELD_VORTICITY;
    curlE->dofNames = vecComponents;
    curlE->unit = "V/m^2";
    curlE->definedOn = ResultInfo::ELEMENT;
    curlE->entryType = ResultInfo::VECTOR;
    availResults_.insert( curlE );
    shared_ptr<CoefFunctionFormBased> curlFunc;
    if( isComplex_ ) {
      curlFunc.reset(new CoefFunctionBOp<Complex>(feFct, curlE));
    } else {
      curlFunc.reset(new CoefFunctionBOp<Double>(feFct, curlE));
    }
    DefineFieldResult( curlFunc, curlE );
    stiffFormCoefs_.insert(curlFunc);


    // === MAGNETIC ENERGY DENSITY INTEGRATED OVER PERIOD  (in the harmonic case)===
    shared_ptr<ResultInfo> jld(new ResultInfo);
    jld->resultType = MAG_ENERGY_DENSITY;
    jld->dofNames = "";
    jld->unit = "";
    jld->definedOn = ResultInfo::ELEMENT;
    jld->entryType = ResultInfo::SCALAR;
    shared_ptr<CoefFunctionMulti> jldCoef(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1, isComplex_));
    DefineFieldResult( jldCoef, jld );  

    // === MAGNETIC ENERGY ===
    shared_ptr<ResultInfo> jldRes(new ResultInfo());
    jldRes->resultType = MAG_ENERGY;
    jldRes->dofNames = "";
    jldRes->unit = "";
    jldRes->definedOn = ResultInfo::REGION;
    jldRes->entryType = ResultInfo::SCALAR;
    availResults_.insert( jldRes );
    shared_ptr<ResultFunctor> jldFunc;
    jldFunc.reset( new ResultFunctorIntegrate<Complex>(jldCoef, feFct, jldRes) );
    resultFunctors_[MAG_ENERGY] = jldFunc;

    // // === PERMEABILITY ===
    // shared_ptr<ResultInfo> permeability ( new ResultInfo );
    // permeability->resultType = MAG_ELEM_PERMEABILITY;
    // permeability->dofNames = "";
    // permeability->unit = "Vs/Am";
    // permeability->definedOn = ResultInfo::ELEMENT;
    // permeability->entryType = ResultInfo::SCALAR;
    // shared_ptr<CoefFunctionMulti> permFct(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1, false));
    // matCoefs_[MAG_ELEM_PERMEABILITY] = permFct;
    // DefineFieldResult(permFct, permeability);
    // availResults_.insert(permeability);

    // === RELUCTIVITY ===
    shared_ptr<ResultInfo> reluctivity ( new ResultInfo );
    reluctivity->resultType = MAG_ELEM_RELUCTIVITY;
    reluctivity->dofNames = "";
    reluctivity->unit = "Am/Vs";
    reluctivity->definedOn = ResultInfo::ELEMENT;
    reluctivity->entryType = ResultInfo::SCALAR;
    shared_ptr<CoefFunctionMulti> relFct(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1, false));
    matCoefs_[MAG_ELEM_RELUCTIVITY] = relFct;
    DefineFieldResult(relFct, reluctivity);
    availResults_.insert(reluctivity);

    // === PERMITTIVITY ===
    shared_ptr<ResultInfo> eps ( new ResultInfo );
    eps->resultType = MAG_ELEM_PERMITTIVITY;
    eps->dofNames = "";
    eps->unit = "As/Vm";
    eps->definedOn = ResultInfo::ELEMENT;
    eps->entryType = ResultInfo::SCALAR;
    shared_ptr<CoefFunctionMulti> epsFct(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1, false));
    matCoefs_[MAG_ELEM_PERMITTIVITY] = epsFct;
    DefineFieldResult(epsFct, eps);
    availResults_.insert(eps);
  }


  void FullWaveMaxwellEPDE::FinalizePostProcResults() {
    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
    StdVector<RegionIdType>::iterator regIt = regions_.Begin();
    
    // Initialize standard postprocessing results
    SinglePDE::FinalizePostProcResults();


    shared_ptr<CoefFunctionMulti> elecIntensCoef = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[ELEC_FIELD_INTENSITY]);

    // ============ E*E averaged over one period ============
    shared_ptr<CoefFunctionMulti> eCoef = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_ENERGY_DENSITY]);
    PtrCoefFct conjEinE = CoefFunction::Generate( mp_, part,
            CoefXprBinOp( mp_, GetCoefFct(ELEC_FIELD_INTENSITY),
                               GetCoefFct(ELEC_FIELD_INTENSITY), CoefXpr::OP_MULT_CONJ) );
    PtrCoefFct halfCoef = CoefFunction::Generate( mp_, part, "0.5");

    regIt = regions_.Begin();
    // for the sake of simplicity we should real with the total current density
    for( ; regIt != regions_.End(); ++regIt ) {
      RegionIdType actRegion = *regIt;
      eCoef->AddRegion(actRegion, CoefFunction::Generate( mp_, part,  CoefXprBinOp(mp_, conjEinE, halfCoef, CoefXpr::OP_MULT) ));
    }



  }



  std::map<SolutionType, shared_ptr<FeSpace>>
  FullWaveMaxwellEPDE::CreateFeSpaces(const std::string &formulation,
                                     PtrParamNode infoNode)
  {
    std::map<SolutionType, shared_ptr<FeSpace>> crSpaces;
    PtrParamNode eVecPotSpaceNode = infoNode->Get("elecFieldIntensity");
    crSpaces[ELEC_FIELD_INTENSITY] = FeSpace::CreateInstance(myParam_, eVecPotSpaceNode, FeSpace::HCURL, ptGrid_);
    crSpaces[ELEC_FIELD_INTENSITY]->Init(solStrat_);

    return crSpaces;
  }

} // end of namespace
