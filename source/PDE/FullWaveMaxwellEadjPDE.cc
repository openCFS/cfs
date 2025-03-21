#include <fstream>
#include <string.h>
#include "FullWaveMaxwellEadjPDE.hh"

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
#include "CoupledPDE/IterCoupledPDE.hh"

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
  DEFINE_LOG(FullWaveMaxwellEadjPDE, "FullWaveMaxwellEadjPDE")

  // **************
  //  Constructor
  // **************
  FullWaveMaxwellEadjPDE::FullWaveMaxwellEadjPDE(Grid *aptgrid, PtrParamNode paramNode,
                                         PtrParamNode infoNode,
                                         shared_ptr<SimState> simState, Domain *domain)
      : SinglePDE(aptgrid, paramNode, infoNode, simState, domain)
  {

    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_ = "fullwave-Eadj";
    std::cout << "In PDE: " << pdename_ << std::endl;
    
    // We reuse the Darwin material class (same as classical magnetic material plus permittivity)
    pdematerialclass_ = ELECTROMAGNETIC_DARWIN;

    //! Always use updated Lagrangian formulation
    updatedGeo_ = true; // true;

    // check if we have a 3d setup
    bool is3d = domain_->GetParamRoot()->Get("domain")->Get("geometryType")->As<std::string>() == "3d";
    if (!is3d)
      EXCEPTION("FullWaveMaxwellEadjPDE is just implemented for 3D setups!");

    // initialize material coef functions covering all regions
    reluc_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, dim_, dim_, isComplex_));
    //perm_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, dim_, dim_, isComplex_));
    conduc_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));
    eps_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));
  }

  // *************
  //  Destructor
  // *************
  FullWaveMaxwellEadjPDE::~FullWaveMaxwellEadjPDE()
  {
  }

  // ****************************
  //  Initialize Nonlinearities
  // ****************************
  void FullWaveMaxwellEadjPDE::InitNonLin()
  {

    SinglePDE::InitNonLin();
  }

  // ****************************
  //  Read damping information
  // ****************************
  void FullWaveMaxwellEadjPDE::ReadDampingInformation()
  {
    // get regions of current PDE
    ParamNodeList regionParamNodes = myParam_->Get("regionList")->GetChildren();
    // corresponding region id
    RegionIdType actRegionId;
    // corresponding region name and damping id
    std::string actRegionName, actDampingId;
    // try to get the dampingList and return if it is not specified
    PtrParamNode dampListNode = myParam_->Get("dampingList", ParamNode::PASS);
    if (dampListNode)
    {
      // get the individual damping nodes
      ParamNodeList dampNodes = dampListNode->GetChildren();
      // map of damping ids from the xml and corresponding damping types
      std::map<std::string, DampingType> idDampType;

      // Run over all region param nodes and assign the required damping information
      for (UInt iRegion = 0; iRegion < regionParamNodes.GetSize(); ++iRegion)
      {
        regionParamNodes[iRegion]->GetValue("name", actRegionName);
        regionParamNodes[iRegion]->GetValue("dampingId", actDampingId);

        // pass if no damping is defined
        if (actDampingId == "")
          continue;

        // parse region id from region name
        actRegionId = ptGrid_->GetRegion().Parse(actRegionName);

        // now, read the damping information from the dampingList
        for (UInt iChild = 0; iChild < dampNodes.GetSize(); ++iChild)
        {
          std::string dampString = dampNodes[iChild]->GetName();
          std::string actId = dampNodes[iChild]->Get("id")->As<std::string>();
          // only consider damping information for the current damping id
          if (actId != actDampingId)
            continue;

          // determine type of damping
          DampingType actType;
          String2Enum(dampString, actType);

          // store damping type string
          idDampType[actId] = actType;
          // break after the information is set as only one damping ID per region is possible
          break;
        }

        // check actDampingId was indeed registerd above
        if (idDampType.count(actDampingId) == 0)
          EXCEPTION("Damping with id '" << actDampingId << "' of region '" << actRegionName << "' was not found. Is it defined in the 'dampingList'?");

        // assign damping type to the region
        dampingList_[actRegionId] = idDampType[actDampingId];

        // if Rayleigh damping is specified, parse and store the additional damping information
        if (dampingList_[actRegionId] == RAYLEIGH)
        {
          RaylDampingData actRaylCoeffs;
          materials_[actRegionId]->GetRayleighCoeffStrings(actRaylCoeffs.alpha, actRaylCoeffs.beta);
          regionRaylDamping_[actRegionId] = actRaylCoeffs;
        }
        else if (dampingList_[actRegionId] == ADAPTED_LOSS_TANGENS_DELTA)
        {
          if (!(analysistype_ == BasePDE::HARMONIC))
            EXCEPTION("Adapted loss tangent delta damping is only allowed for harmonic analysis.");
          RaylDampingData actRaylCoeffs;
          materials_[actRegionId]->GetFreqAdaptedRayleighCoeffStrings(actRaylCoeffs.alpha, actRaylCoeffs.beta);
          regionRaylDamping_[actRegionId] = actRaylCoeffs;
        }
        else if (dampingList_[actRegionId] == GLOBAL_RAYLEIGH)
        {
          EXCEPTION("Global Rayleigh damping is not yet implemented.");
          if (dampNodes.GetSize() != 1)
            EXCEPTION("Global Rayleigh damping does not allow for additional damping nodes defined.");
          RaylDampingData actRaylCoeffs;
          actRaylCoeffs.alpha = dampNodes[0]->Get("alpha")->As<std::string>();
          actRaylCoeffs.beta = dampNodes[0]->Get("beta")->As<std::string>();
          regionRaylDamping_[actRegionId] = actRaylCoeffs;
        }

        // set flag to compute extra integrators for transient PML
        if (dampingList_[actRegionId] == PML && analysistype_ == BasePDE::TRANSIENT)
        {
          isTimeDomPML_ = true;
        }
      }
    }
  }

  // *****************************
  //  Definition of Integrators
  // *****************************
  void FullWaveMaxwellEadjPDE::DefineIntegrators()
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



    shared_ptr<BaseFeFunction> eVecPotFeFunc = feFunctions_[ELEC_FIELD_INTENSITY_ADJ];
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
        
          //get coeff-Fnc for the magnetic reluctivity
          ReadMaterialDependency( "matReluctivity", resultInfo->dofNames, resultInfo->entryType, false,
                                    entity, reluctivity, updatedGeo_ );

          // //get coeff-Fnc for the drivative of magnetic reluctivity w.r.t. parameter, and store it for postprocessing
          // PtrCoefFct reluctivityDeriv = NULL;                                    
          // ReadMaterialDependency( "magElemReluctivityDeriv", resultInfo->dofNames, resultInfo->entryType, false,
          //                           entity, reluctivityDeriv, updatedGeo_ );
          // if ( reluctivityDeriv == NULL )
          //   std::cout << "No magElemReluctivityDeriv!!" << std::endl;

          PtrCoefFct constOne = CoefFunction::Generate(mp_, Global::REAL, "1.0");                                    
          derivReluctivity_[actRegion] = constOne; //reluctivityDeriv;

          //get the curl of the electric field of forward simulation
          std::string pdeName = "fullwave-E";
          PtrCoefFct forwardEVortex = iterCplPde_->GetCouplingCoefFct(ELEC_FIELD_VORTICITY, actSDList, pdeName, updatedGeo_); 
          //forwardEVortex->IsComplex(true);
          EvorticityForward_[actRegion] = forwardEVortex; //iterCplPde_->GetCouplingCoefFct(ELEC_FIELD_VORTICITY, actSDList, pdeName, updatedGeo_);          
        } 
        else {
          reluctivity = actMat->GetScalCoefFnc(MAG_RELUCTIVITY_SCALAR, Global::REAL);
        }
        // Add material to global, distributed reluctivity coefficient function
        reluc_->AddRegion(actRegion, reluctivity);              
        matCoefs_[MAG_ELEM_RELUCTIVITY]->AddRegion(actRegion, reluctivity);

        // Magnetic Permeability
        PtrCoefFct constOne = CoefFunction::Generate(mp_, Global::REAL, "1.0");
        PtrCoefFct permeability = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, constOne, 
                                                         reluctivity, CoefXpr::OP_DIV));
        // perm_->AddRegion(actRegion, permeability);  
        // matCoefs_[MAG_ELEM_PERMEABILITY]->AddRegion(actRegion, permeability);


        // Electric / Magnetic Permittivity
        PtrCoefFct eps = NULL;        
        if ( matDepenTypes.Find(MAT_PERMITTIVITY) != -1 ) {
          //we read the computed permitivity from file
          shared_ptr<ResultInfo> resultInfo = GetResultInfo(MAG_ELEM_PERMITTIVITY); 
          std::string regionName = ptGrid_->GetRegion().ToString(actRegion);
          shared_ptr<EntityList> entity = ptGrid_->GetEntityList( EntityList::ELEM_LIST, regionName );
        
          //get coeff-Fnc for the magnetic permittivity
          ReadMaterialDependency( "matPermittivity", resultInfo->dofNames, resultInfo->entryType, false,
                                  entity, eps, updatedGeo_ );                                  

          //get coeff-Fnc for the drivative of magnetic permittivty w.r.t. parameter, and store it for postprocessing
          // PtrCoefFct permDeriv = NULL;                                    
          // ReadMaterialDependency( "matPermittivityDeriv", resultInfo->dofNames, resultInfo->entryType, false,
          //                         entity, permDeriv, updatedGeo_ );
          // if (permDeriv == NULL )
          //   std::cout << "No magElemPermittivityDeriv!!" << std::endl;

          PtrCoefFct constOne = CoefFunction::Generate(mp_, Global::REAL, "1.0");    
          derivPermittivity_[actRegion] = constOne;    

          //get the curl of the electric field of forward simulation
          std::string pdeName = "fullwave-E";
          EfieldForward_[actRegion] = iterCplPde_->GetCouplingCoefFct(ELEC_FIELD_INTENSITY, actSDList, pdeName, updatedGeo_);                                          
        } 
        else {        
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
        massInts_[actRegion] = M_E_E_epsilon;  
      } // END OF NONLIN/LIN PART
    }   // end for regions
  }     // end DefineIntegrators

  void FullWaveMaxwellEadjPDE::DefineNcIntegrators()
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
            DefineNitscheCoupling<3,1>(ELEC_FIELD_INTENSITY_ADJ, *ncIt );
        }
          break;
      }
      default:
        EXCEPTION("Unknown type of ncInterface");
        break;
      }
    }
  }

  void FullWaveMaxwellEadjPDE::DefineSurfaceIntegrators()
  {
  }

  void FullWaveMaxwellEadjPDE::DefineRhsLoadIntegrators()
  {
    //LOG_TRACE(FullWaveMaxwellEPDE) << "Defining rhs load integrators for FullWaveMaxwellEPDE";

    shared_ptr<BaseFeFunction> eVecPotFeFunc = feFunctions_[ELEC_FIELD_INTENSITY_ADJ];

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

    //read electric field intensity from forward simulation
    StdVector<std::string> nameOfDofs;
    ReadRhsExcitation("elecFieldIntensity", nameOfDofs, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo);    
    for (UInt i = 0; i < ent.GetSize(); ++i)
    {
       // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST ||
          ent[i]->GetType() == EntityList::SURF_ELEM_LIST ) {
        EXCEPTION("Electric field intensity from forward simulation can only be specified in a volume!")
      }
          
      if(isComplex_) {
        lin = new BUIntegrator<Complex>( new IdentityOperator<FeHCurl, 3, 1, Complex>(),
                                         Complex(1.0), coef[i], coefUpdateGeo);
      } else {
        lin = new BUIntegrator<Double>( new IdentityOperator<FeHCurl, 3, 1, Double>(),
                                        1.0, coef[i], coefUpdateGeo);
      }

      lin->SetName("RHSFieldintensityIntegrator");
      
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(eVecPotFeFunc);
      assemble_->AddLinearForm(ctx);
      eVecPotFeFunc->AddEntityList(ent[i]);

      bRHSRegions_[ent[i]->GetRegion()] = coef[i];
    }
  }

  void FullWaveMaxwellEadjPDE::DefineSolveStep()
  {
    solveStep_ = new StdSolveStep(*this);
  }



  // ======================================================
  // TIME-STEPPING SECTION
  // ======================================================

  void FullWaveMaxwellEadjPDE::InitTimeStepping()
  {
    // Use complete implicit scheme
    Double gamma = 1.0;
    GLMScheme *scheme = new Trapezoidal(gamma);
    TimeSchemeGLM::NonLinType nlType = (nonLin_) ? TimeSchemeGLM::INCREMENTAL : TimeSchemeGLM::NONE;
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(scheme, 0, nlType));
    feFunctions_[ELEC_FIELD_INTENSITY_ADJ]->SetTimeScheme(myScheme);
  }


  void FullWaveMaxwellEadjPDE::DefinePrimaryResults()
  {

    StdVector<std::string> vecComponents;
    vecComponents = "x", "y", "z";

    // === ELECTRIC FIELD INTENSITY ===
    shared_ptr<ResultInfo> elecIntens(new ResultInfo);
    elecIntens->resultType = ELEC_FIELD_INTENSITY_ADJ;
    elecIntens->SetVectorDOFs(dim_, isaxi_);
    elecIntens->dofNames = vecComponents;
    elecIntens->unit = "V/m";
    elecIntens->definedOn = ResultInfo::ELEMENT;
    elecIntens->entryType = ResultInfo::VECTOR;
    feFunctions_[ELEC_FIELD_INTENSITY_ADJ]->SetResultInfo(elecIntens);
    DefineFieldResult(feFunctions_[ELEC_FIELD_INTENSITY_ADJ], elecIntens);

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    hdbcSolNameMap_[ELEC_FIELD_INTENSITY_ADJ] = "fluxParallel";
    idbcSolNameMap_[ELEC_FIELD_INTENSITY_ADJ] = "potential";

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

  void FullWaveMaxwellEadjPDE::DefinePostProcResults()
  {
    StdVector<std::string> vecComponents;
    vecComponents = "x", "y", "z";

    shared_ptr<BaseFeFunction> feFct = feFunctions_[ELEC_FIELD_INTENSITY_ADJ];

    // === CURL OF ELECTRIC FIELD = ELECTRIC VORTICITY ===
    shared_ptr<ResultInfo> curlEadj(new ResultInfo);
    curlEadj->resultType = ELEC_FIELD_VORTICITY_ADJ;
    curlEadj->dofNames = vecComponents;
    curlEadj->unit = "V/m^2";
    curlEadj->definedOn = ResultInfo::ELEMENT;
    curlEadj->entryType = ResultInfo::VECTOR;
    availResults_.insert( curlEadj );
    shared_ptr<CoefFunctionFormBased> curlFunc;
    if( isComplex_ ) {
      curlFunc.reset(new CoefFunctionBOp<Complex>(feFct, curlEadj));
    } else {
      curlFunc.reset(new CoefFunctionBOp<Double>(feFct, curlEadj));
    }
    DefineFieldResult( curlFunc, curlEadj );
    stiffFormCoefs_.insert(curlFunc);

    // === Adjoint MAGNETIC RHS ===
    shared_ptr<ResultInfo> rhs(new ResultInfo);
    rhs->resultType = ELEC_RHS_LOAD_ADJ;
    rhs->dofNames = vecComponents;
    rhs->unit = "-";
    rhs->entryType = ResultInfo::VECTOR;
    rhs->definedOn = ResultInfo::ELEMENT;
    rhsFeFunctions_[ELEC_FIELD_INTENSITY_ADJ]->SetResultInfo(rhs);
    DefineFieldResult( rhsFeFunctions_[ELEC_FIELD_INTENSITY_ADJ], rhs );


    // === just define it: will be used in FinalizePostProcResult
    shared_ptr<ResultInfo> averagedEV(new ResultInfo);
    averagedEV->resultType = ELEC_FIELD_VORTICITY;
    averagedEV->dofNames = vecComponents;
    averagedEV->unit = "V/m^2";
    averagedEV->entryType = ResultInfo::VECTOR;
    averagedEV->definedOn = ResultInfo::ELEMENT;
    // The computation is defined in FinalizePostProcResults()
    shared_ptr<CoefFunctionMulti> averagedEVfnc(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_));
    DefineFieldResult(averagedEVfnc, averagedEV);  

    // === just define it: will be used in FinalizePostProcResult
    shared_ptr<ResultInfo> averagedEI(new ResultInfo);
    averagedEI->resultType = ELEC_FIELD_INTENSITY;
    averagedEI->dofNames = vecComponents;
    averagedEI->unit = "V/m";
    averagedEI->entryType = ResultInfo::VECTOR;
    averagedEI->definedOn = ResultInfo::ELEMENT;
    // The computation is defined in FinalizePostProcResults()
    shared_ptr<CoefFunctionMulti> averagedEIfnc(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_));
    DefineFieldResult(averagedEIfnc, averagedEI);  

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


  void FullWaveMaxwellEadjPDE::FinalizePostProcResults() {
    //Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
    StdVector<RegionIdType>::iterator regIt = regions_.Begin();
    
    // Initialize standard postprocessing results
    SinglePDE::FinalizePostProcResults();

    shared_ptr<CoefFunctionMulti> elecIntensCoef = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[ELEC_FIELD_INTENSITY_ADJ]);

    // ======================== Computation of gardients for SIMP =================
    StdVector<std::string> vecComponents;
    vecComponents = "x", "y", "z";

    shared_ptr<BaseFeFunction> feFct = feFunctions_[ELEC_FIELD_INTENSITY_ADJ];

      // === Pure helper results ===
      shared_ptr<ResultInfo> opt1(new ResultInfo);
      opt1->resultType = OPT_RESULT_1; 
      opt1->dofNames = "";
      opt1->unit = "";
      opt1->definedOn = ResultInfo::ELEMENT;
      opt1->entryType = ResultInfo::SCALAR;      
      shared_ptr<CoefFunctionMulti> optMult1(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));
      DefineFieldResult(optMult1, opt1);

      shared_ptr<ResultInfo> opt2(new ResultInfo);
      opt2->resultType = OPT_RESULT_2; 
      opt2->dofNames = "";
      opt2->unit = "";
      opt2->definedOn = ResultInfo::ELEMENT;
      opt2->entryType = ResultInfo::SCALAR;      
      shared_ptr<CoefFunctionMulti> optMult2(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));
      DefineFieldResult(optMult2, opt2);

      // === Gradient via adjoint ===
      // === Parameter 1 is reluctivity
      shared_ptr<ResultInfo> gradAdjParam1;
      gradAdjParam1.reset(new ResultInfo);
      gradAdjParam1->resultType = GRAD_PARAM1;
      gradAdjParam1->dofNames = "";
      gradAdjParam1->unit = "";
      gradAdjParam1->entryType = ResultInfo::SCALAR;
      gradAdjParam1->definedOn = ResultInfo::ELEMENT; //REGION;   
      availResults_.insert(gradAdjParam1);                         
      shared_ptr<ResultFunctor> gradElem1;
      gradElem1.reset(new ResultFunctorIntegrate<Complex>(optMult1, feFct, gradAdjParam1));
      resultFunctors_[GRAD_PARAM1] = gradElem1;

      // === Parameter 2 is permittivity
      shared_ptr<ResultInfo> gradAdjParam2;
      gradAdjParam2.reset(new ResultInfo);
      gradAdjParam2->resultType = GRAD_PARAM2;
      gradAdjParam2->dofNames = "";
      gradAdjParam2->unit = "";
      gradAdjParam2->entryType = ResultInfo::SCALAR;
      gradAdjParam2->definedOn = ResultInfo::ELEMENT; //REGION;   
      availResults_.insert(gradAdjParam2);                         
      shared_ptr<ResultFunctor> gradElem2;
      gradElem2.reset(new ResultFunctorIntegrate<Complex>(optMult2, feFct, gradAdjParam2));
      resultFunctors_[GRAD_PARAM2] = gradElem2;

      //compute gradient
      shared_ptr<CoefFunctionMulti> scalMult1 = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[OPT_RESULT_1]);    
      shared_ptr<CoefFunctionMulti> eVorticityField = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[ELEC_FIELD_VORTICITY]);      
      shared_ptr<CoefFunctionMulti> scalMult2 = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[OPT_RESULT_2]);    
      shared_ptr<CoefFunctionMulti> eField = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[ELEC_FIELD_INTENSITY]);  

      regIt = regions_.Begin();
      for (; regIt != regions_.End(); ++regIt) {
        // ========= Electric Vorticity Field =============
        PtrCoefFct curlForward = NULL;
        PtrCoefFct mult1 = NULL;
        RegionIdType actRegion = *regIt;
        if(EvorticityForward_.find(*regIt) != EvorticityForward_.end()){
          // Curl of E = ELEC_FIELD_VORTICITY from forward simulation
          curlForward = EvorticityForward_[*regIt];
          // get ELEC_FIELD_VORTICITY_ADJ is curl(Eadj), where Eadj is the adjoint solution!
          PtrCoefFct curlAdj = this->GetCoefFct(ELEC_FIELD_VORTICITY_ADJ);          
          mult1 = CoefFunction::Generate( mp_, Global::COMPLEX,
                             CoefXprBinOp(mp_, curlForward, curlAdj, CoefXpr::OP_MULT_CONJ) );        
          scalMult1->AddRegion(actRegion, mult1);

          //curl of E-field from forward simulation
          eVorticityField->AddRegion(*regIt, curlForward);    
        }

        PtrCoefFct eForward = NULL;
        PtrCoefFct mult2 = NULL;
        if(EfieldForward_.find(*regIt) != EfieldForward_.end()){
          // Electric field intensity = ELEC_FIELD_INTENSITY from forward simulation
          eForward = EfieldForward_[*regIt];
          // get coef of ELEC_FIELD_INTENSITY_ADJ: is adjoint solution!
          PtrCoefFct eFieldAdj = this->GetCoefFct(ELEC_FIELD_INTENSITY_ADJ);
          mult2 = CoefFunction::Generate( mp_, Global::COMPLEX,
                             CoefXprBinOp(mp_, eForward, eFieldAdj, CoefXpr::OP_MULT_CONJ) );        
          scalMult2->AddRegion(actRegion, mult2);

          //curl of E-field from forward simulation
          eField->AddRegion(*regIt, eForward);    
        }
      } // loop over regions

  }



  std::map<SolutionType, shared_ptr<FeSpace>>
  FullWaveMaxwellEadjPDE::CreateFeSpaces(const std::string &formulation,
                                     PtrParamNode infoNode)
  {
    std::map<SolutionType, shared_ptr<FeSpace>> crSpaces;
    PtrParamNode eVecPotSpaceNode = infoNode->Get("elecFieldIntensityAdj");
    crSpaces[ELEC_FIELD_INTENSITY_ADJ] = FeSpace::CreateInstance(myParam_, eVecPotSpaceNode, FeSpace::HCURL, ptGrid_);
    crSpaces[ELEC_FIELD_INTENSITY_ADJ]->Init(solStrat_);

    return crSpaces;
  }

} // end of namespace
