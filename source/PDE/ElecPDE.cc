#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <set>

#include "ElecPDE.hh"

#include "General/defs.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Domain/CoefFunction/CoefFunction.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/CoefFunction/CoefFunctionPML.hh"
#include "Domain/Mesh/NcInterfaces/MortarInterface.hh"
#include "Domain/CoefFunction/CoefFunctionHyst.hh"
#include "Domain/CoefFunction/CoefFunctionMapping.hh"
#include "Utils/StdVector.hh"
#include "Driver/SolveSteps/SolveStepElec.hh"
#include "Driver/SolveSteps/SolveStepHyst.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"
#include "Driver/Assemble.hh"
#include "Utils/ApproxData.hh"
#include "Utils/SmoothSpline.hh"
#include "Materials/Models/Hysteresis.hh"
#include "Materials/Models/Preisach.hh"
//#include "Materials/Models/VectorPreisach.hh"
#include "FeBasis/H1/FeSpaceH1Nodal.hh"
#include "FeBasis/FeFunctions.hh"

// new integrator concept
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/SingleEntryInt.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/DivOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/IdentityOperatorNormal.hh"
#include "Forms/Operators/SurfaceOperators.hh"
#include "Forms/Operators/SurfaceNormalFluxDensityOperator.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"

// multiharmonic stuff
#include "Driver/MultiHarmonicDriver.hh"
#include "Domain/CoefFunction/CoefFunctionHarmBalance.hh"
#include "Domain/CoefFunction/CoefFunctionMulti.hh"
#include <Domain/CoefFunction/CoefFunctionMaterialModel.hh>





namespace CoupledField {

  DEFINE_LOG(elecpde, "pde.electrostatic")
  
  // ***************
  //   Constructor
  // ***************
  ElecPDE::ElecPDE( Grid* aptgrid, PtrParamNode paramNode,
          PtrParamNode infoNode,
          shared_ptr<SimState> simState, Domain* domain )
  :SinglePDE( aptgrid, paramNode, infoNode, simState, domain ) {
    
    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_          = "electrostatic";
    pdematerialclass_ = ELECTROSTATIC;
    
    nonLin_         = false;
    nonLinMaterial_ = false;
    //    isAlwaysStatic_ = true;
    isPiezoCoupled_ = false;
    
    //! Always use updated Lagrangian formulation 
    updatedGeo_     = true;
    
    regionApproxSet_ = false;
    // Check the subtype of the problem
    paramNode->GetValue("subType", subType_);

    //this is the default value
    modelName_ = "nonlinearCurve";
    // for multiharmonic
    multiHarmCoef_.reset(new CoefFunctionHarmBalance<Complex>());
    // Be aware that the hysteresis which are build in via CoefFunctionMaterialModel does not use the isHysteresis_ flag!
    matModelCoef_.reset( new CoefFunctionMaterialModel<Complex>());
    //init the nlFluxCoef
    nlFluxCoef_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_, true));
  }
  
  void ElecPDE::InitNonLin() {
    
    SinglePDE::InitNonLin();
  }

  ElecPDE::~ElecPDE(){
    if(multiHarmCoef_){
      // Part of a bandaid solution to fix faulty memory management
      // TODO remove when proper memory management is implemented
      dynamic_cast<CoefFunctionHarmBalance<Complex>&>(*multiHarmCoef_).DeregisterMathParser();
      multiHarmCoef_.reset();
    }
  }
  
  void ElecPDE::ReadDampingInformation()
  {
    std::map<std::string, DampingType> idDampType;
    
    // try to get dampingList
    PtrParamNode dampListNode = myParam_->Get("dampingList", ParamNode::PASS);
    if (dampListNode)
    {
      // get specific damping nodes
      ParamNodeList dampNodes = dampListNode->GetChildren();
      
      for (UInt i = 0; i < dampNodes.GetSize(); i++)
      {
        std::string dampString = dampNodes[i]->GetName();
        std::string actId = dampNodes[i]->Get("id")->As<std::string>();
        
        // determine type of damping
        DampingType actType;
        String2Enum(dampString, actType);
        // store damping type string
        idDampType[actId] = actType;
      }
    }
    
    // Run over all region and set entry in "regionNonLinId"
    ParamNodeList regionNodes = myParam_->Get("regionList")->GetChildren();
    
    RegionIdType actRegionId;
    std::string actRegionName, actDampingId;
    
    for (UInt k = 0; k < regionNodes.GetSize(); k++)
    {
      regionNodes[k]->GetValue("name", actRegionName);
      regionNodes[k]->GetValue("dampingId", actDampingId);
      if (actDampingId == "")
        continue;
      
      actRegionId = ptGrid_->GetRegion().Parse(actRegionName);
      
      // Check actDampingId was already registerd
      if (idDampType.count(actDampingId) == 0)
      {
        EXCEPTION("Damping with id '" << actDampingId << "' was not defined in 'dampingList'");
      }
      
      dampingList_[actRegionId] = idDampType[actDampingId];
    }
    
    // Check, if all entries are identical
    for (UInt i = 1; i < dampingList_.size(); i++)
      if (dampingList_[regions_[i - 1]] != dampingList_[regions_[i]])
        break;
  }
  
  /*
   * Initialize hysteresis coeficient functions
   * 
   * Note to this function:
   *  as long as hysteresis is only used in elecPDE, it is sufficient to initialize
   *  the hysteresis coeficient function for each region directly before creating
   *  the integrators for stiffness matrix and rhs
   * 
   *  in the piezocoupled case, however, the coupled integrators are initialized
   *  prior to the single-field integrator; as those coupled integrators may depend
   *  on the hysteresis coef. function, it is important that these are already 
   *  initizalized; therefore, this extra function is called during the definition
   *  of the post-processing results (where the actual storage for the hyst-coef-functions
   *  gets obtained) which is prior to the creation of both coupled and single field
   *  integrators
   * 
   */
  void ElecPDE::InitHystCoefs() {
    //std::cout << "Init Hyst Coefs" << std::endl;
    RegionIdType actRegion;
    BaseMaterial * actSDMat = NULL;
    
    SubTensorType tensorType;
    
    if ( dim_ == 3  || subType_ == "2.5d" ) {
      tensorType = FULL;
    }
    else {
      if ( isaxi_ == true ) {
        tensorType = AXI;
      }
      else {
        // 2d: plane case
        tensorType = PLANE_STRAIN;
      }
    }
    
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    shared_ptr<FeSpace> mySpace = feFunctions_[ELEC_POTENTIAL]->GetFeSpace();
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      
      // Set current region and material
      actRegion = it->first;
      actSDMat = it->second;
      
      // Get current region name
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);
      
      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( actRegion );
      
      // ==========================================================
      // New implementation
      // ==========================================================
      // --- Set the approximation for the current region ---
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();
      mySpace->SetRegionApproximation(actRegion, polyId, integId);
      
      SDLists_[actRegion] = actSDList;
      
      const StdVector<NonLinType>& nonLinTypes = regionNonLinTypes_[actRegion];
      if ( nonLinTypes.Find(HYSTERESIS) != -1 ){
        shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
        actSDList->SetRegion( actRegion );
        
        std::string hystType;
        actSDMat->GetString(hystType, HYST_MODEL);
        
        if(hystType == "none"){
          std::string warnmsg = "Hysteresis set on region " + regionName + " but no hysteresis model was defined in mat file. Skip.";
          regionNonLinTypes_[actRegion] = NO_NONLINEARITY;
          WARN(warnmsg);
        } else {
          PtrCoefFct elecFieldCoef = this->GetCoefFct(ELEC_FIELD_INTENSITY);
          PtrCoefFct elecFieldCoefSurf = this->GetCoefFct(ELEC_FIELD_INTENSITY_SURF);
          PtrCoefFct hystPol(new CoefFunctionHyst( actSDMat, actSDList,
                  elecFieldCoef,elecFieldCoefSurf,tensorType,ELEC_PERMITTIVITY_TENSOR,mySpace));


          hysteresisCoefs_->AddRegion( actRegion, hystPol);

        }
      }
    }
    regionApproxSet_ = true;
    //std::cout << "Hyst Coefs set successfully" << std::endl;
  }
  
  void ElecPDE::DefineIntegrators() {
    
    LOG_DBG(elecpde) << "Define Integrators";
    RegionIdType actRegion;
    BaseMaterial * actSDMat = NULL;
    
    // CopyPasta from MagEdgePDE and hope this works:

    shared_ptr<BaseFeFunction> feFunc = feFunctions_[ELEC_POTENTIAL];
    shared_ptr<FeSpace> feSpace = feFunc->GetFeSpace();

    PtrCoefFct elecFieldCoef =  this->GetCoefFct(ELEC_FIELD_INTENSITY);

    UInt baseFreq=0, N, M, nFFT;

    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;

    // Init material model for transient analysis
    if (((analysistype_ == STATIC) || (analysistype_ == TRANSIENT)) && nonLin_ && (modelName_ != "nonlinearCurve")){
        matModelCoef_->Init( elecFieldCoef, modelName_, dim_);
    }
    // Create new harmonic balance coefficient function and register the regions and material
    if(analysistype_ == MULTIHARMONIC){
        if(modelName_ != "nonlinearCurve"){
          matModelCoef_->Init( elecFieldCoef, modelName_, dim_);
        }

        baseFreq = dynamic_cast<MultiHarmonicDriver*>(domain_->GetSingleDriver())->baseFreq_;
        N = dynamic_cast<MultiHarmonicDriver*>(domain_->GetSingleDriver())->numHarmonics_N_;
        M = dynamic_cast<MultiHarmonicDriver*>(domain_->GetSingleDriver())->numHarmonics_M_;
        nFFT = dynamic_cast<MultiHarmonicDriver*>(domain_->GetSingleDriver())->numFFT_;

        multiHarmCoef_->Init(feFunc, feSpace, regions_, materials_, ptGrid_, elecFieldCoef, N, M, baseFreq, nFFT, modelName_,matModelCoef_);

    }
    // now the handle for the mathparser should be in place, and ready to unlock!
    // bool upLagrangeForm = true;
    
    //transform the type
    SubTensorType tensorType;
    
    if ( dim_ == 3  || subType_ == "2.5d" ) {
      tensorType = FULL;
    }
    else {
      if ( isaxi_ == true ) {
        tensorType = AXI;
      }
      else {
        // 2d: plane case
        tensorType = PLANE_STRAIN;
      }
    }
    
    
    // if the pde is piezo-coupled, the electrostatic entries
    // have to be multiplied with -1
    std::string factor = "1.0";
    if ( isPiezoCoupled_ == true )
      factor = "-1.0";  // get applied to the stiffness matrix in function GetStiffIntegrator
    
    // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    
    shared_ptr<FeSpace> mySpace = feFunctions_[ELEC_POTENTIAL]->GetFeSpace();
    
    //flag indicating frequency PML formulation
    bool harmonicPML = false;

    // Changed iteration for the loop, like in magnetics, so i can use
    // copy-pasta code
    for(UInt iRegion = 0; iRegion < regions_.GetSize() ; iRegion ++){
      actRegion = regions_[iRegion];
      actSDMat    = materials_[actRegion];

      // Get current region name
      shared_ptr<BaseFeFunction> myFct = feFunctions_[ELEC_POTENTIAL];
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);
      
      shared_ptr<ElemList> actSDList;
      if(!regionApproxSet_){
        // create new entity list
        shared_ptr<ElemList> newSDList( new ElemList(ptGrid_ ) );
        actSDList = newSDList;
        actSDList->SetRegion( actRegion );
        SDLists_[actRegion] = actSDList;
      } else {
        actSDList = SDLists_[actRegion];
      }
      
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();
      
      if(!regionApproxSet_){
        // ==========================================================
        // New implementation
        // ==========================================================
        // --- Set the approximation for the current region ---
        mySpace->SetRegionApproximation(actRegion, polyId, integId);
      }
      
      myFct->AddEntityList( actSDList );
      
      // Take account for pml in frequency domain
      // 'coeffPMLScal' is the function used to scale the material tensor. If no PML is defined, it's unity
      PtrCoefFct coefPMLScal;
      PtrCoefFct coefPMLVec, speedOfSnd;
      PtrParamNode pmlNode;
      std::string pmlFormul;
      PtrCoefFct coefMAPScal, coefMAPVec;
      bool isMapping = false;
      
      const StdVector<NonLinType>& nonLinTypes = regionNonLinTypes_[actRegion];

      if (dampingList_[actRegion] == PML)
      {
        if (analysistype_ == HARMONIC)
        {
          harmonicPML = true;
          std::string dampId;
          curRegNode->GetValue("dampingId", dampId);
          pmlNode = myParam_->Get("dampingList")->GetByVal("pml", "id", dampId.c_str());
          pmlFormul = pmlNode->Get("formulation")->As<std::string>();
          
          // speed of sound is set to equal '1.0'
          speedOfSnd = CoefFunction::Generate(mp_, Global::REAL, "1.0");
          if (pmlFormul == "classic")
          {
            coefPMLScal.reset(new CoefFunctionPML<Complex>(pmlNode, speedOfSnd,
                    ptGrid_->GetEntityList(EntityList::ELEM_LIST, regionName), regions_, false));
            coefPMLVec.reset(new CoefFunctionPML<Complex>(pmlNode, speedOfSnd,
                    ptGrid_->GetEntityList(EntityList::ELEM_LIST, regionName), regions_, true));
          }
          else if (pmlFormul == "shifted")
          {
            coefPMLScal.reset(new CoefFunctionShiftedPML<Complex>(pmlNode, speedOfSnd,
                    ptGrid_->GetEntityList(EntityList::ELEM_LIST, regionName), regions_, false));
            coefPMLVec.reset(new CoefFunctionShiftedPML<Complex>(pmlNode, speedOfSnd,
                    ptGrid_->GetEntityList(EntityList::ELEM_LIST, regionName), regions_, true));
          }
          else
          {
            EXCEPTION("Unknown PML-formulation '" << pmlFormul << "'");
          }
        }
        else
          EXCEPTION("Not implemented yet");
      }
      else
      {
        if( dampingList_[actRegion] == MAPPING ) {
          // ====================================================================
          // Take account for mapping of an infinite domain
          // ====================================================================
          // Generate scalar valued coefficient function
          PtrCoefFct val1 = CoefFunction::Generate( mp_, Global::REAL, "1.0");
          std::string dampId;
          curRegNode->GetValue("dampingId",dampId);
          PtrParamNode mapNode = myParam_->Get("dampingList")->GetByVal("mapping","id",dampId.c_str());
          coefMAPVec.reset(new CoefFunctionMapping<Double>(mapNode,val1,actSDList,regions_,true));
          coefMAPScal.reset(new CoefFunctionMapping<Double>(mapNode,val1,actSDList,regions_,false));
          isMapping = true;
        }else{
          harmonicPML = false;
        }
      }

      // ----- standard real-valued stiffness integrator
      BaseBDBInt* stiffInt = NULL;
      PtrCoefFct epsilonNL = NULL;
      if (harmonicPML)
      {
        stiffInt = GetStiffIntegrator(actSDMat, tensorType, actRegion, coefPMLScal);
        stiffInt->SetBCoefFunctionOpA(coefPMLVec);
      }
      else{
        if(isMapping){
          // Infinte mapping case
          stiffInt = GetStiffIntegratorInfMap(actSDMat, tensorType, actRegion, coefMAPScal);
          stiffInt->SetBCoefFunctionOpA(coefMAPVec);

        } else if (nonLinTypes.Find(NLELEC_PERMITTIVITY) != -1){

          PtrCoefFct elecFieldCoef =  this->GetCoefFct(ELEC_FIELD_INTENSITY);

          if (analysistype_==MULTIHARMONIC){ //if nonlinear and Multiharmonic do this:
               bool nL = (nonLinTypes.GetSize() > 0)? true : false;
               epsilonNL = multiHarmCoef_->GenerateMatCoefFnc(iRegion, "Permittivity", nL, actSDList);
               if( dim_== 2){
                 EXCEPTION("Not implemented yet, but try simply CopyPasta the 3Dim Case.");
               }
               stiffInt = new BBInt<Complex>(new GradientOperator<FeH1,3>(),epsilonNL, (Complex)1.0, updatedGeo_);
          } else{ //transient - nonlinear or hysteretic
            PtrCoefFct elecFieldCoef =  this->GetCoefFct(ELEC_FIELD_INTENSITY);

            if(modelName_ == "JilesAthertonModel" ){
              std::map<std::string, double> ParameterMap;

              actSDMat->GetScalar(ParameterMap["Ps"], ELEC_PS_JILES, Global::REAL );
              actSDMat->GetScalar(ParameterMap["alpha"], ELEC_ALPHA_JILES, Global::REAL );
              actSDMat->GetScalar(ParameterMap["a"], ELEC_A_JILES, Global::REAL );
              actSDMat->GetScalar(ParameterMap["k"], ELEC_K_JILES, Global::REAL );
              actSDMat->GetScalar(ParameterMap["c"], ELEC_C_JILES, Global::REAL );

              ParameterMap["isMH"]=0;

              matModelCoef_->InitModel( ParameterMap,actSDList->GetSize());
              epsilonNL = matModelCoef_; //actSDMat->GetScalCoefFncModel( matModelCoef_ );

              PtrCoefFct hystFluxTmp = CoefFunction::Generate( mp_, part, CoefXprBinOp(mp_,epsilonNL,elecFieldCoef,CoefXpr::OP_MULT));
              nlFluxCoef_->AddRegion(actRegion, hystFluxTmp);

            }else{
              PtrCoefFct elecFieldCoef =  this->GetCoefFct(ELEC_FIELD_INTENSITY);
              epsilonNL = actSDMat->GetScalCoefFncNonLin( ELEC_PERMITTIVITY_SCALAR, Global::REAL, elecFieldCoef);

              PtrCoefFct hystFluxTmp = CoefFunction::Generate( mp_, part, CoefXprBinOp(mp_,epsilonNL,elecFieldCoef,CoefXpr::OP_MULT));
              nlFluxCoef_->AddRegion(actRegion, hystFluxTmp);
            }
            if( dim_ == 2 ) {
              stiffInt = new BBInt<>(new GradientOperator<FeH1,2>(), epsilonNL, 1.0, updatedGeo_);
            } else {
              stiffInt = new BBInt<>(new GradientOperator<FeH1,3>(), epsilonNL, 1.0, updatedGeo_);
            }
          }

        } else{
          stiffInt = GetStiffIntegrator(actSDMat, tensorType, actRegion);
        }
      }
      if (nonLinTypes.Find(NLELEC_PERMITTIVITY) != -1){
        stiffInt->SetName("NonLinElecIntegrator");

        BiLinFormContext* stiffIntDescr =
        new BiLinFormContext(stiffInt, STIFFNESS);

        stiffIntDescr->SetEntities( actSDList, actSDList );
        stiffIntDescr->SetFeFunctions( myFct, myFct);
        stiffInt->SetFeSpace( myFct->GetFeSpace());

        assemble_->AddBiLinearForm( stiffIntDescr );
        bdbInts_.insert( std::pair<RegionIdType, BaseBDBInt*>(actRegion,stiffInt) );

      } else{
        stiffInt->SetName("LinElecIntegrator");
        BiLinFormContext * stiffIntDescr =
                new BiLinFormContext(stiffInt, STIFFNESS );

        //myFct->AddEntityList( actSDList );

        //stiffIntDescr->SetPtPdes(this, this);
        stiffIntDescr->SetEntities( actSDList, actSDList );
        stiffIntDescr->SetFeFunctions( myFct, myFct);
        stiffInt->SetFeSpace( myFct->GetFeSpace());

        assemble_->AddBiLinearForm( stiffIntDescr );
        // Important: Add bdb-integrator to global list, as we need them later
        // for calculation of postprocessing results
        bdbInts_.insert( std::pair<RegionIdType, BaseBDBInt*>(actRegion,stiffInt) );
      
      }
    }
    
    // Define integrators for composite materials
    // (only for subType "flatShell")
    std::map<RegionIdType, Composite>::iterator compIt;
    for( compIt=compositeMaterials_.begin(); compIt!=compositeMaterials_.end();
            compIt++ ) {
      
      REFACTOR;
    }
    
    // define integrators for electric impedances
    DefineImpedanceIntegrators();
  }
  
  void ElecPDE::DefineSurfaceIntegrators()
  {
    PtrParamNode bcNode = myParam_->Get("bcsAndLoads", ParamNode::PASS);
    this->ptGrid_->MapEdges();
    this->ptGrid_->MapFaces();
    if (!bcNode)
      return;
    
    ParamNodeList blochNodesList = bcNode->GetList("blochPeriodic");
    for (UInt i = 0; i < blochNodesList.GetSize(); i++)
    {
      std::string str_value = blochNodesList[i]->Get("factor_value")->As<std::string>();
      std::string str_phase = blochNodesList[i]->Get("factor_phase")->As<std::string>();
      std::string formulation = blochNodesList[i]->Get("formulation")->As<std::string>();
      
      // propagation factor \gamma from xml-file
      std::string str_real, str_imag;
      str_real = AmplPhaseToReal(str_value, str_phase, true);
      str_imag = AmplPhaseToImag(str_value, str_phase, true);
      
      PtrCoefFct factor = CoefFunction::Generate(mp_, Global::COMPLEX, str_real, str_imag);
      PtrCoefFct one = CoefFunction::Generate(mp_, Global::REAL, "1.0", "0.0");
      
      ParamNodeList regionsList = blochNodesList[i]->GetList("region");
      for (UInt j = 0; j < regionsList.GetSize(); j++)
      {
        std::string ncRegionName = regionsList[j]->Get("name")->As<std::string>();
        shared_ptr<BaseNcInterface> ncIf = ptGrid_->GetNcInterface(ptGrid_->GetNcInterfaceId(ncRegionName));
        if (!ncIf)
        {
          EXCEPTION("No interface with the name '" << ncRegionName << "' found!");
        }
        shared_ptr<MortarInterface> mortarIf = dynamic_pointer_cast<MortarInterface>(ncIf);
        assert(mortarIf);
        
        // TODO: WHAT THE HELL? IT DOESN'T WORK WITH pz = -1...
        // To be checked later during the work with piezo-coupling
        Double pz = 1.0;
        if (isPiezoCoupled_ == true)
          pz = 1.0;
        
        if (formulation == "Nitsche")
        {
          PtrCoefFct matDataTensorMas, matDataTensorSla, matData;
          RegionIdType volMasterId = mortarIf->GetPrimaryVolRegion();
          RegionIdType volSlaveId = mortarIf->GetSecondaryVolRegion();
          
          matDataTensorMas = regionPermittivity_[volMasterId];
          matDataTensorSla = regionPermittivity_[volSlaveId];
          assert(matDataTensorMas);
          assert(matDataTensorSla);
          
          std::string nitFac = blochNodesList[i]->Get("nitscheFactor")->As<std::string>();
          Double nitscheFactor = boost::lexical_cast<Double>(nitFac);
          // master & slave penalty integrals
          BiLinearForm *pnlt_PhiM_PsiM = NULL;
          BiLinearForm *pnlt_PhiM_PsiS = NULL;
          BiLinearForm *pnlt_PhiS_PsiM = NULL;
          BiLinearForm *pnlt_PhiS_PsiS = NULL;
          // master & slave integrals with normal derivatives
          BiLinearForm *flux_DPhiM_PsiM = NULL;
          BiLinearForm *flux_DPhiM_PsiS = NULL;
          BiLinearForm *flux_PhiM_DPsiM = NULL;
          BiLinearForm *flux_PhiS_DPsiM = NULL;
          
          shared_ptr<ElemList> actSDList = ncIf->GetElemList();
          Double beta;
          PtrCoefFct factorSqr = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, factor, factor, CoefXpr::OP_MULT));
          
          // obtain a proper scaling for the penalty terms
          StdVector<Vector<Double> > points(1);
          Vector<Double> p1(dim_);
          p1.Init();
          points[0]= p1;
          
          if (matDataTensorMas->IsComplex())
          {
            Complex tmp(0.0, 0.0);
            StdVector<Matrix<Complex> > valsM, valsS;
            matDataTensorMas->GetTensorValuesAtCoords(points, valsM, this->ptGrid_);
            matDataTensorSla->GetTensorValuesAtCoords(points, valsS, this->ptGrid_);
            for (UInt k = 0; k < valsM[0].GetNumRows(); k++)
            {
              tmp += valsM[0][k][k]*conj(valsM[0][k][k]);
              tmp += valsS[0][k][k]*conj(valsS[0][k][k]);
            }
            beta = sqrt(0.5*tmp.real())*nitscheFactor;
          }
          else
          {
            Double tmp(0.0);
            StdVector<Matrix<Double> > valsM, valsS;
            matDataTensorMas->GetTensorValuesAtCoords(points, valsM, this->ptGrid_);
            matDataTensorSla->GetTensorValuesAtCoords(points, valsS, this->ptGrid_);
            for (UInt k = 0; k < valsM[0].GetNumRows(); k++)
            {
              tmp += valsM[0][k][k]*valsM[0][k][k];
              tmp += valsS[0][k][k]*valsS[0][k][k];
            }
            beta = sqrt(0.5*tmp)*nitscheFactor;
          }
          
          PtrCoefFct coefFuncPMLVec, coefFuncPMLScl;
          if (dampingList_[volMasterId] == PML)
          {
            if (analysistype_ == HARMONIC)
            {
              std::string regionName = ptGrid_->GetRegion().ToString(volMasterId);
              std::string dampId, pmlFormul;
              
              PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region", "name", regionName.c_str());
              curRegNode->GetValue("dampingId", dampId);
              PtrParamNode pmlNode = myParam_->Get("dampingList")->GetByVal("pml", "id", dampId.c_str());
              
              // speed of sound is set to equal '1.0'
              PtrCoefFct speedOfSnd = CoefFunction::Generate(mp_, Global::REAL, "1.0");
              pmlFormul = pmlNode->Get("formulation")->As<std::string>();
              
              if (pmlFormul == "classic")
              {
                coefFuncPMLVec.reset(new CoefFunctionPML<Complex>(pmlNode, speedOfSnd,
                        ptGrid_->GetEntityList(EntityList::ELEM_LIST, regionName), regions_, true));
                coefFuncPMLScl.reset(new CoefFunctionPML<Complex>(pmlNode, speedOfSnd,
                        ptGrid_->GetEntityList(EntityList::ELEM_LIST, regionName), regions_, false));
              }
              else if (pmlFormul == "shifted")
              {
                coefFuncPMLVec.reset(new CoefFunctionShiftedPML<Complex>(pmlNode, speedOfSnd,
                        ptGrid_->GetEntityList(EntityList::ELEM_LIST, regionName), regions_, true));
                coefFuncPMLScl.reset(new CoefFunctionShiftedPML<Complex>(pmlNode, speedOfSnd,
                        ptGrid_->GetEntityList(EntityList::ELEM_LIST, regionName), regions_, false));
              }
              else
              {
                EXCEPTION("Unknown PML-formulation '" << pmlFormul << "'");
              }
              
              matData = CoefFunction::Generate(mp_, Global::COMPLEX,
                      CoefXprTensScalOp(mp_, matDataTensorMas, coefFuncPMLScl, CoefXpr::OP_MULT));
            }
            else
              EXCEPTION("The analysis type '" << this->analysistype_ << "' is not supported!");
          }
          else
            matData = matDataTensorMas;
          
          // define bilinear forms for Nitsche coupling
          // penalty integrators
          pnlt_PhiM_PsiM = GetPenaltyIntegrator<Complex>(factor, beta, BiLinearForm::PRIM_PRIM);
          pnlt_PhiM_PsiS = GetPenaltyIntegrator<Complex>(factorSqr, -beta, BiLinearForm::PRIM_SEC);
          pnlt_PhiS_PsiM = GetPenaltyIntegrator<Double>(one, -beta, BiLinearForm::SEC_PRIM);
          pnlt_PhiS_PsiS = GetPenaltyIntegrator<Complex>(factor, beta, BiLinearForm::SEC_SEC);
          // flux integrators
          if (matData->IsComplex())
          {
            flux_DPhiM_PsiM = GetFluxIntegrator<Complex>(one, coefFuncPMLVec, -1.0*pz, BiLinearForm::PRIM_PRIM, true);
            flux_PhiM_DPsiM = GetFluxIntegrator<Complex>(factor, coefFuncPMLVec, -1.0*pz, BiLinearForm::PRIM_PRIM, false);
            flux_DPhiM_PsiS = GetFluxIntegrator<Complex>(factor, coefFuncPMLVec, 1.0*pz, BiLinearForm::PRIM_SEC, true);
            flux_PhiS_DPsiM = GetFluxIntegrator<Complex>(one, coefFuncPMLVec, 1.0*pz, BiLinearForm::SEC_PRIM, false);
          }
          else
          {
            flux_DPhiM_PsiM = GetFluxIntegrator<Double>(one, coefFuncPMLVec, -1.0*pz, BiLinearForm::PRIM_PRIM, true);
            flux_PhiM_DPsiM = GetFluxIntegrator<Complex>(factor, coefFuncPMLVec, -1.0*pz, BiLinearForm::PRIM_PRIM, false);
            flux_DPhiM_PsiS = GetFluxIntegrator<Complex>(factor, coefFuncPMLVec, 1.0*pz, BiLinearForm::PRIM_SEC, true);
            flux_PhiS_DPsiM = GetFluxIntegrator<Double>(one, coefFuncPMLVec, 1.0*pz, BiLinearForm::SEC_PRIM, false);
          }
          
          // pass material data to the flux operators
          flux_DPhiM_PsiM->SetBCoefFunctionOpA(matData);
          flux_PhiM_DPsiM->SetBCoefFunctionOpB(matData);
          flux_DPhiM_PsiS->SetBCoefFunctionOpA(matData);
          flux_PhiS_DPsiM->SetBCoefFunctionOpB(matData);
          
          // master-master
          pnlt_PhiM_PsiM->SetName("pnlt_PhiM_PsiM");
          flux_DPhiM_PsiM->SetName("flux_DPhiM_PsiM");
          flux_PhiM_DPsiM->SetName("flux_PhiM_DPsiM");
          //master-slave
          pnlt_PhiM_PsiS->SetName("pnlt_PhiM_PsiS");
          flux_DPhiM_PsiS->SetName("flux_DPhiM_PsiS");
          // slave-master
          pnlt_PhiS_PsiM->SetName("pnlt_PhiS_PsiM");
          flux_PhiS_DPsiM->SetName("flux_PhiS_DPsiM");
          //slave-slave
          pnlt_PhiS_PsiS->SetName("pnlt_PhiS_PsiS");
          
          // BiLinearForm::PRIM_PRIM
          SurfaceBiLinFormContext *pnlt_PhiM_PsiM_cont = new SurfaceBiLinFormContext(pnlt_PhiM_PsiM, STIFFNESS, BiLinearForm::PRIM_PRIM);
          SurfaceBiLinFormContext *flux_DPhiM_PsiM_cont = new SurfaceBiLinFormContext(flux_DPhiM_PsiM, STIFFNESS, BiLinearForm::PRIM_PRIM);
          SurfaceBiLinFormContext *flux_PhiM_DPsiM_cont = new SurfaceBiLinFormContext(flux_PhiM_DPsiM, STIFFNESS, BiLinearForm::PRIM_PRIM);
          // BiLinearForm::PRIM_SEC
          SurfaceBiLinFormContext *pnlt_PhiM_PsiS_cont = new SurfaceBiLinFormContext(pnlt_PhiM_PsiS, STIFFNESS, BiLinearForm::PRIM_SEC);
          SurfaceBiLinFormContext *flux_DPhiM_PsiS_cont = new SurfaceBiLinFormContext(flux_DPhiM_PsiS, STIFFNESS, BiLinearForm::PRIM_SEC);
          // BiLinearForm::SEC_PRIM
          SurfaceBiLinFormContext *pnlt_PhiS_PsiM_cont = new SurfaceBiLinFormContext(pnlt_PhiS_PsiM, STIFFNESS, BiLinearForm::SEC_PRIM);
          SurfaceBiLinFormContext *flux_PhiS_DPsiM_cont = new SurfaceBiLinFormContext(flux_PhiS_DPsiM, STIFFNESS, BiLinearForm::SEC_PRIM);
          // BiLinearForm::SEC_SEC
          SurfaceBiLinFormContext *pnlt_PhiS_PsiS_cont = new SurfaceBiLinFormContext(pnlt_PhiS_PsiS, STIFFNESS, BiLinearForm::SEC_SEC);
          
          pnlt_PhiM_PsiM_cont->SetEntities(actSDList, actSDList);
          flux_DPhiM_PsiM_cont->SetEntities(actSDList, actSDList);
          flux_PhiM_DPsiM_cont->SetEntities(actSDList, actSDList);
          pnlt_PhiM_PsiS_cont->SetEntities(actSDList, actSDList);
          flux_DPhiM_PsiS_cont->SetEntities(actSDList, actSDList);
          pnlt_PhiS_PsiM_cont->SetEntities(actSDList, actSDList);
          flux_PhiS_DPsiM_cont->SetEntities(actSDList, actSDList);
          pnlt_PhiS_PsiS_cont->SetEntities(actSDList, actSDList);
          
          pnlt_PhiM_PsiM_cont->SetFeFunctions(feFunctions_[ELEC_POTENTIAL], feFunctions_[ELEC_POTENTIAL]);
          flux_DPhiM_PsiM_cont->SetFeFunctions(feFunctions_[ELEC_POTENTIAL], feFunctions_[ELEC_POTENTIAL]);
          flux_PhiM_DPsiM_cont->SetFeFunctions(feFunctions_[ELEC_POTENTIAL], feFunctions_[ELEC_POTENTIAL]);
          pnlt_PhiM_PsiS_cont->SetFeFunctions(feFunctions_[ELEC_POTENTIAL], feFunctions_[ELEC_POTENTIAL]);
          flux_DPhiM_PsiS_cont->SetFeFunctions(feFunctions_[ELEC_POTENTIAL], feFunctions_[ELEC_POTENTIAL]);
          pnlt_PhiS_PsiM_cont->SetFeFunctions(feFunctions_[ELEC_POTENTIAL], feFunctions_[ELEC_POTENTIAL]);
          flux_PhiS_DPsiM_cont->SetFeFunctions(feFunctions_[ELEC_POTENTIAL], feFunctions_[ELEC_POTENTIAL]);
          pnlt_PhiS_PsiS_cont->SetFeFunctions(feFunctions_[ELEC_POTENTIAL], feFunctions_[ELEC_POTENTIAL]);
          
          assemble_->AddBiLinearForm(pnlt_PhiM_PsiM_cont);
          assemble_->AddBiLinearForm(flux_DPhiM_PsiM_cont);
          assemble_->AddBiLinearForm(flux_PhiM_DPsiM_cont);
          assemble_->AddBiLinearForm(pnlt_PhiM_PsiS_cont);
          assemble_->AddBiLinearForm(flux_DPhiM_PsiS_cont);
          assemble_->AddBiLinearForm(pnlt_PhiS_PsiM_cont);
          assemble_->AddBiLinearForm(flux_PhiS_DPsiM_cont);
          assemble_->AddBiLinearForm(pnlt_PhiS_PsiS_cont);
        } // end nitsche
        else if (formulation == "Mortar")
        {
          shared_ptr<SurfElemList> surfMasterGrid(new SurfElemList(ptGrid_)), surfSlaveGrid(new SurfElemList(ptGrid_));
          surfMasterGrid->SetRegion(mortarIf->GetPrimarySurfRegion());
          surfSlaveGrid->SetRegion(mortarIf->GetSecondarySurfRegion());
          
          // --- Set the approximation for Lagrange Multipliers for the current region ---
          RegionIdType regId = surfSlaveGrid->GetRegion();
          std::string polyId = regionsList[j]->Get("polyId")->As<std::string>();
          std::string integId = regionsList[j]->Get("integId")->As<std::string>();
          feFunctions_[LAGRANGE_MULT]->GetFeSpace()->SetRegionApproximation(regId, polyId, integId);
          
          BiLinearForm *intOne1 = 0;
          BiLinearForm *intFactor1 = 0;
          BiLinearForm *intOne2 = 0;
          BiLinearForm *intFactor2 = 0;
          
          if (dim_ == 2)
          {
            intOne1 = new SurfaceMortarABInt<Complex, Complex>(new IdentityOperator<FeH1, 2, 1, Complex>(),
                    new IdentityOperator<FeH1, 2, 1, Complex>(),
                    one, pz,
                    mortarIf->GetSecondaryVolRegion(), mortarIf->GetPrimaryVolRegion(),
                    mortarIf->IsCoplanar(), updatedGeo_, BiLinearForm::SEC_PRIM);
            intFactor1 = new BBInt<Complex, Complex>(new IdentityOperator<FeH1, 2, 1, Complex>(), factor, -pz, updatedGeo_);
            intOne2 = new BBInt<Complex, Complex>(new IdentityOperator<FeH1, 2, 1, Complex>(), one, pz, updatedGeo_);
            intFactor2 = new SurfaceMortarABInt<Complex, Complex>(new IdentityOperator<FeH1, 2, 1, Complex>(),
                    new IdentityOperator<FeH1, 2, 1, Complex>(),
                    factor, -pz,
                    mortarIf->GetPrimaryVolRegion(), mortarIf->GetSecondaryVolRegion(),
                    mortarIf->IsCoplanar(), updatedGeo_, BiLinearForm::PRIM_SEC);
          }
          else
          {
            intOne1 = new SurfaceMortarABInt<Complex, Complex>(new IdentityOperator<FeH1, 3, 1, Complex>(),
                    new IdentityOperator<FeH1, 3, 1, Complex>(),
                    one, pz,
                    mortarIf->GetSecondaryVolRegion(), mortarIf->GetPrimaryVolRegion(),
                    mortarIf->IsCoplanar(), updatedGeo_, BiLinearForm::SEC_PRIM);
            intFactor1 = new BBInt<Complex, Complex>(new IdentityOperator<FeH1, 3, 1, Complex>(), factor, -pz, updatedGeo_);
            intOne2 = new BBInt<Complex, Complex>(new IdentityOperator<FeH1, 3, 1, Complex>(), one, pz, updatedGeo_);
            intFactor2 = new SurfaceMortarABInt<Complex, Complex>(new IdentityOperator<FeH1, 3, 1, Complex>(),
                    new IdentityOperator<FeH1, 3, 1, Complex>(),
                    factor, -pz,
                    mortarIf->GetPrimaryVolRegion(), mortarIf->GetSecondaryVolRegion(),
                    mortarIf->IsCoplanar(), updatedGeo_, BiLinearForm::PRIM_SEC);
          }
          
          intOne1->SetName("master1Elec");
          intFactor1->SetName("slave1Elec");
          intOne2->SetName("master2Elec");
          intFactor2->SetName("slave2Elec");
          
          // (1) weak form of the periodic boundary conditions
          BiLinFormContext *lagPotContMaster = new BiLinFormContext(intFactor1, STIFFNESS);
          NcBiLinFormContext *lagPotContSlave = new NcBiLinFormContext(intOne1, STIFFNESS);
          
          lagPotContMaster->SetEntities(surfSlaveGrid, surfSlaveGrid);
          lagPotContSlave->SetEntities(mortarIf->GetElemList(), mortarIf->GetElemList());
          lagPotContMaster->SetFeFunctions(feFunctions_[LAGRANGE_MULT], feFunctions_[ELEC_POTENTIAL]);
          lagPotContSlave->SetFeFunctions(feFunctions_[LAGRANGE_MULT], feFunctions_[ELEC_POTENTIAL]);
          
          assemble_->AddBiLinearForm(lagPotContMaster);
          assemble_->AddBiLinearForm(lagPotContSlave);
          
          // (2) weak form of the boundary integrals of the PDE
          BiLinFormContext *potLagContMaster = new BiLinFormContext(intOne2, STIFFNESS);
          NcBiLinFormContext *potLagContSlave = new NcBiLinFormContext(intFactor2, STIFFNESS);
          
          potLagContMaster->SetEntities(surfSlaveGrid, surfSlaveGrid);
          potLagContSlave->SetEntities(mortarIf->GetElemList(), mortarIf->GetElemList());
          potLagContMaster->SetFeFunctions(feFunctions_[ELEC_POTENTIAL], feFunctions_[LAGRANGE_MULT]);
          potLagContSlave->SetFeFunctions(feFunctions_[ELEC_POTENTIAL], feFunctions_[LAGRANGE_MULT]);
          
          assemble_->AddBiLinearForm(potLagContMaster);
          assemble_->AddBiLinearForm(potLagContSlave);
          
          feFunctions_[LAGRANGE_MULT]->AddEntityList(surfSlaveGrid);
          feFunctions_[ELEC_POTENTIAL]->AddEntityList(surfMasterGrid);
          feFunctions_[ELEC_POTENTIAL]->AddEntityList(surfSlaveGrid);
        } // end mortar
        else
        {
          EXCEPTION("Unknown formulation: '" << formulation << "'!");
        }
      }
    }
    
    /*
     * Field-Parallel boundary
     * 
     * i.e. En = 0
     * 
     * for non-hysteretic materials 
     * field-parallel and flux-parallel (default-bc) coicide as
     *  D = eps*E
     * and thus
     *  En = 0 -> Dn = 0
     * 
     * for hysteretic materials
     * field-parallel and flux-parallel do not coincide as
     *  D = eps0*E + P
     * and thus
     *  En = 0 -> Dn = Pn
     * 
     * > in the later case, the additional surface integral
     * 
     *  -int_Gamma N Pn dGamma
     * 
     * has to be added where N is the standard shape functions
     *  
     * Note that this boundary conditions can also be applied to interfaces
     * between a hysteretic and a non-hysteretic materials
     * In that case, Dn = Pn will act like influences surface charges
     * Attention: In that case volumeTegion has to be the hysteretic material region
     *              > specify this in xml file
     */ 
    if( isHysteresis_ ){
      std::set<RegionIdType> volRegions (regions_.Begin(), regions_.End() );
      ParamNodeList fpNodeList = bcNode->GetList( "fieldParallel" );
      std::map<RegionIdType,PtrCoefFct > regionCoefs = hysteresisCoefs_->GetRegionCoefs();
      
      bool fullevaluation = true;
      const bool isSurface = true;
      bool coefUpdateGeo = true;
      
      shared_ptr<BaseFeFunction> myFct = feFunctions_[ELEC_POTENTIAL];
      // Note; in the piezoelectric case we have to multiply by -1
      Double factor = 1.0;
      if ( isPiezoCoupled_ ){
        factor = -1.0;
      }
      
      for( UInt i = 0; i < fpNodeList.GetSize(); i++ ) {
        std::string regionName = fpNodeList[i]->Get("name")->As<std::string>();
        shared_ptr<EntityList> actSDList =  ptGrid_->GetEntityList( EntityList::SURF_ELEM_LIST,regionName );
        std::string volRegName = fpNodeList[i]->Get("volumeRegion")->As<std::string>();
        
        RegionIdType volReg = ptGrid_->GetRegion().Parse(volRegName);
        
        // check if volReg has hyst material behaviour
        if(regionCoefs.find(volReg) == regionCoefs.end()){
          std::stringstream warnmsg;
          warnmsg << "Volume region " << volRegName << "has NO hysteretic material assigned." << std::endl;
          warnmsg << "Field parallel BC will thus act as default flux parallel BC." << std::endl;
//          WARN(warnmsg.str());
        } else {
          //std::cout << "Volume region " << volRegName << " has hysteretic material assigned. fieldParallel BC added" << std::endl;
          // create coef fnc delivering the boundary term (here just polarization)
          PtrCoefFct hystBC = regionCoefs[volReg]->GenerateRHSCoefFnc("ElecPolarization",isSurface);
          
          LinearForm * lin = NULL;
          
          if(isComplex_) {
            if( dim_ == 2 ) {
              //            lin = new BUIntegrator<Complex>( new GradientOperator<FeH1,2,1,Complex>(),
              //                                           Complex(factor),it->second,  coefUpdateGeo, fullevaluation);
              lin = new BUIntegrator<Complex,isSurface>( new IdentityOperatorNormal<FeH1,2>(),
                      Complex(factor),hystBC, volRegions, coefUpdateGeo, fullevaluation);
            } else {
              lin = new BUIntegrator<Complex,isSurface>( new IdentityOperatorNormal<FeH1,2>(),
                      Complex(factor),hystBC,volRegions,  coefUpdateGeo, fullevaluation);
            }
          } else  {
            if( dim_ == 2 ) {
              lin = new BUIntegrator<Double,isSurface>( new IdentityOperatorNormal<FeH1,2>(),
                      (factor),hystBC, volRegions, coefUpdateGeo, fullevaluation);
            } else {
              lin = new BUIntegrator<Double,isSurface>( new IdentityOperatorNormal<FeH1,2>(),
                      (factor),hystBC, volRegions, coefUpdateGeo, fullevaluation);
            }
          }
          //factor = factor*(-1.0);
          lin->SetName("fieldParallel_HystLoading");
          lin->SetSolDependent();
          LinearFormContext *ctx = new LinearFormContext( lin );
          ctx->SetEntities( actSDList );
          ctx->SetFeFunction(myFct);
          assemble_->AddLinearForm(ctx);
          // Add entity list will add nothing, if entities were already assigned
          myFct->AddEntityList(actSDList);
          
          // IMPORTANT: add surface elements to hyst operator such that it gets
          // storage space assigned
					bool onBoundary = true;
          regionCoefs[volReg]->AddAdditionalSDList(actSDList,volReg,onBoundary);
        }
      }
    }
    
  }
  
  void ElecPDE::DefineRhsLoadIntegrators() {
    LOG_TRACE(elecpde) << "Defining rhs load integrators for electrostatic PDE";
    
    // Get FESpace and FeFunction of electric potential
    shared_ptr<BaseFeFunction> myFct = feFunctions_[ELEC_POTENTIAL];
    shared_ptr<FeSpace> mySpace = myFct->GetFeSpace();
    
    StdVector<shared_ptr<EntityList> > ent;
    StdVector<PtrCoefFct > coef;
    LinearForm * lin = NULL;
    StdVector<std::string> dofNames;
    
    // Note; in the piezoelectric case we have to multiply by -1
    Double factor = 1.0;
    if ( isPiezoCoupled_ ){
      factor = -1.0;
    }
    
    // Flag, if coefficient function lives on updated geoemtry
    bool coefUpdateGeo = true;
    
    // =========================
    //  Charges (volume, nodal)
    // =========================
    LOG_DBG(elecpde) << "Reading charges";
    ReadRhsExcitation( "charge", dofNames, ResultInfo::VECTOR, 
            isComplex_, ent, coef,coefUpdateGeo );
 
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        
        // ---------------
        //  Nodal Charges 
        // ---------------
        // Nodal charge must be constant
        if(coef[i]->IsSpacialDependent())
          EXCEPTION("Nodal charges must not be spatial dependent");
        
        UInt numNodes = ent[i]->GetSize();
        Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;  
        coef[i] = CoefFunction::Generate(mp_, part, 
                CoefXprBinOp(mp_, coef[i], 
                boost::lexical_cast<std::string>(numNodes*factor), CoefXpr::OP_DIV) );
        
        lin = new SingleEntryInt(coef[i]);
        lin->SetName("NodalChargeInt");
        LinearFormContext *ctx = new LinearFormContext( lin );
        ctx->SetEntities( ent[i] );
        ctx->SetFeFunction(myFct);
        assemble_->AddLinearForm(ctx);
        myFct->AddEntityList(ent[i]);
      } else {
        // --------------------------
        //  Surface / Volume Charges 
        // --------------------------
        if( coef[i]->GetDependency() == CoefFunction::GENERAL || coef[i]->GetDependency() == CoefFunction::SPACE ) {
          EXCEPTION("The total prescribed charge must no be spatial dependent");
        }
        // "Divide" the total charge by the volume / surface of the current entity list
        Double volume = ptGrid_->CalcVolumeOfEntityList( ent[i], false );
        Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;  
        coef[i] = CoefFunction::Generate(mp_, part, 
                CoefXprVecScalOp(mp_, coef[i], 
                boost::lexical_cast<std::string>(volume), CoefXpr::OP_DIV) );
        if(isComplex_) {
          lin = new BUIntegrator<Complex>( new IdentityOperator<FeH1>(),
                  Complex(factor), coef[i],coefUpdateGeo);
        } else  {
          lin = new BUIntegrator<Double>( new IdentityOperator<FeH1>(),
                  factor, coef[i],coefUpdateGeo);
        }
        lin->SetName("ChargeDensityInt");
        LinearFormContext *ctx = new LinearFormContext( lin );
        ctx->SetEntities( ent[i] );
        ctx->SetFeFunction(myFct);
        assemble_->AddLinearForm(ctx);
        myFct->AddEntityList(ent[i]);
      }
    } // for
    
    // ================
    //  CHARGE DENSITY 
    // ================
    LOG_DBG(elecpde) << "Reading charge densities";
    ReadRhsExcitation( "chargeDensity", dofNames, 
            ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("Charge density must be defined on elements")
      }
      if(isComplex_) {
        lin = new BUIntegrator<Complex>( new IdentityOperator<FeH1>(),
                Complex(factor), coef[i], coefUpdateGeo);
      } else  {
        lin = new BUIntegrator<Double>( new IdentityOperator<FeH1>(),
                factor, coef[i], coefUpdateGeo);
      }
      lin->SetName("ChargeDensityInt");
      if(analysistype_ == MULTIHARMONIC){
        //ToDo Maybe can do it in a more generic way, like for the BC
        // So the framework is already there, if we want to implement for other pdes
        // but havent looked into this...
        ParamNodeList chargeDensitylist = myParam_->Get("bcsAndLoads", ParamNode::PASS)->GetList("chargeDensity");
        UInt harm = chargeDensitylist[i]->Get("harmonic")->As<UInt>();
        UInt harm_max=dynamic_cast<MultiHarmonicDriver*>(domain_->GetSingleDriver())->numHarmonics_N_;
        if(harm > harm_max ){
          EXCEPTION("RHS harmonic is bigger than the highest harmonic of the system.")
        }
        lin->SetHarm(harm);
      }
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
    } // for
    
    // =================================
    //  Polarisation -> from hysteresis (VOLUME)
    // =================================
    /*
     * NEW: we have to put P on the rhs also for deltaMaterial stepping!
     * -> see StdSolveStep::StepTransNonLinHysteresis for details
     */
    if ( isHysteresis_ ){
      // std::cout << "Putting polarisation to rhs" << std::endl;
      LOG_DBG(elecpde) << "Putting polarisation to rhs";
      
      std::map<RegionIdType,PtrCoefFct > regionCoefs = hysteresisCoefs_->GetRegionCoefs();
      std::map<RegionIdType, shared_ptr<CoefFunction> > ::iterator it;
      for( it = regionCoefs.begin(); it != regionCoefs.end(); it++) {
        
        // get regionIdType
        RegionIdType curReg = it->first;
        
        // get SDList
        shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
        actSDList->SetRegion( curReg );
        
        // set fullevaluation to trigger evaluation at each integration point
        // the nonlinear parameter "evaluation depth" determines if each
        // integration point gets mapped to midpoint (> fullevaluation = false)
        // or if hyst operator really is evaluated at the actual int. point
        bool fullevaluation = true;
        
        shared_ptr<CoefFunction> rhsPol = it->second->GenerateRHSCoefFnc("ElecPolarization");

        //factor = factor*(-1.0);
        if(isComplex_) {
          if( dim_ == 2 ) {
            //            lin = new BUIntegrator<Complex>( new GradientOperator<FeH1,2,1,Complex>(),
            //                                           Complex(factor),it->second,  coefUpdateGeo, fullevaluation);
            lin = new BUIntegrator<Complex>( new GradientOperator<FeH1,2,1,Complex>(),
                    Complex(factor),rhsPol,  coefUpdateGeo, fullevaluation);
          } else {
            lin = new BUIntegrator<Complex>( new GradientOperator<FeH1,3,1,Complex>(),
                    Complex(factor),rhsPol,  coefUpdateGeo, fullevaluation);
          }
        } else  {
          if( dim_ == 2 ) {
            lin = new BUIntegrator<Double>( new GradientOperator<FeH1,2> (),
                    (factor),rhsPol,  coefUpdateGeo, fullevaluation);
          } else {
            lin = new BUIntegrator<Double>( new GradientOperator<FeH1,3> (),
                    (factor),rhsPol,  coefUpdateGeo, fullevaluation);
          }
        }
        //factor = factor*(-1.0);
        lin->SetName("rhs_polarization");
        lin->SetSolDependent();
        LinearFormContext *ctx = new LinearFormContext( lin );
        ctx->SetEntities( actSDList );
        ctx->SetFeFunction(myFct);
        assemble_->AddLinearForm(ctx);
        // Add entity list will add nothing, if entities were already assigned
        myFct->AddEntityList(actSDList);
      }
    }
    
    // =================================
    //  FLUX DENSITY (VOLUME / SURFACE) 
    // =================================
    LOG_DBG(elecpde) << "Reading prescribed flux density";
    StdVector<std::string> vecDofNames;
    if(dim_ == 3 || subType_== "2.5d")
    {
      vecDofNames = "x", "y", "z";
    }
    else
    {
      if(dim_ == 2 && !isaxi_)
        vecDofNames = "x", "y";
      if(dim_ == 2 && isaxi_)
        vecDofNames = "r", "z";
    }
    ReadRhsExcitation( "fluxDensity", vecDofNames, 
            ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo );
    std::set<RegionIdType> volRegions (regions_.Begin(), regions_.End() );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("Charge density must be defined on elements")
      }
      
      // determine dimension
      EntityIterator it = ent[i]->GetIterator();
      UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
      if( elemDim == (dim_-1) ) {
        // === SURFACE ===
        if( dim_ == 2) {
          if(isComplex_) {
            lin = new BUIntegrator<Complex, true>( new IdentityOperatorNormal<FeH1,2>(),
                    Complex(factor), coef[i], volRegions, coefUpdateGeo);
          } else  {
            lin = new BUIntegrator<Double, true>( new IdentityOperatorNormal<FeH1,2>(),
                    factor, coef[i], volRegions, coefUpdateGeo);
          } 
        }else {
          if(isComplex_) {
            lin = new BUIntegrator<Complex, true>( new IdentityOperatorNormal<FeH1,3>(),
                    Complex(factor), coef[i], volRegions, coefUpdateGeo);
          } else  {
            lin = new BUIntegrator<Double, true>( new IdentityOperatorNormal<FeH1,3>(),
                    factor, coef[i], volRegions, coefUpdateGeo);
          } 
        }
        lin->SetName("SurfaceFluxDensityInt");
      } else {
        // === VOLUME ===
        EXCEPTION("Elec Flux density in the volume not yet implemented");
      }
      
      
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
    } // for
    
    // =========================================
    //  PERMANENT POLARIZATION (VOLUME) 
    // =========================================
    LOG_DBG(elecpde) << "Reading prescribed polarization";
    ReadRhsExcitation( "polarization", vecDofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) { EXCEPTION("Polarization must be defined on elements") }
      
      // determine elements dimension
      EntityIterator it = ent[i]->GetIterator();
      UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
      if (elemDim != dim_) { EXCEPTION("Polarization elements must have same dimension as geometry") }
      
      // === VOLUME ===
      if(isComplex_)
      {
        if( dim_ == 2 ) { lin = new BUIntegrator<Complex>( new GradientOperator<FeH1,2,1,Complex>(), Complex(factor), coef[i], coefUpdateGeo); }
        else { lin = new BUIntegrator<Complex>( new GradientOperator<FeH1,3,1,Complex>(), Complex(factor), coef[i], coefUpdateGeo ); }
      } 
      else
      {
        if( dim_ == 2 ) { lin = new BUIntegrator<Double>( new GradientOperator<FeH1,2> (), factor, coef[i], coefUpdateGeo ); }
        else { lin = new BUIntegrator<Double>( new GradientOperator<FeH1,3> (), factor, coef[i], coefUpdateGeo ); }
      }
      
      lin->SetName("PermanentPolarizationInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
    } // for
  }
  
  BaseBDBInt * ElecPDE::GetStiffIntegrator( BaseMaterial* actSDMat,
          SubTensorType tensorType,
          RegionIdType regionId ) {
    
    BaseBDBInt * integ = NULL;
    bool isComplex = complexMatData_[regionId];
    
    shared_ptr<CoefFunction > curCoef;
    //get possible nonlinearities defined in this region
    bool regionIsHyst = false;
    const StdVector<NonLinType>& nonLinTypes = regionNonLinTypes_[regionId];
    if ( nonLinTypes.Find(HYSTERESIS) != -1 ){
      //std::cout << "Found region with hysteresis" << std::endl;
      regionIsHyst = true;
      /* for both the delta material method as well as the std fixpoint method we have to know
       * which regions are affected by hysteresis
       */
      // NEW: coefFncHyst should already be created during DefinePostProcResults!      
      PtrCoefFct hystPol = hysteresisCoefs_->GetRegionCoef(regionId);
      
      //std::cout << "Generate Mat Coef Fnc" << std::endl;
      curCoef = hystPol->GenerateMatCoefFnc("Permittivity");
      
      //std::cout << "Generate Output Coef Fnc" << std::endl;
      PtrCoefFct hystOutput = hystPol->GenerateOutputCoefFnc("ElecPolarization");
      hysteresisPolarization_->AddRegion( regionId, hystOutput);
      
    }
    else {
      
      if ( isComplex ) {
        curCoef = actSDMat->GetTensorCoefFnc( ELEC_PERMITTIVITY_TENSOR,tensorType,
                Global::COMPLEX);
      }
      else {
        curCoef = actSDMat->GetTensorCoefFnc( ELEC_PERMITTIVITY_TENSOR,tensorType,
                Global::REAL);
      }
    }
    
    // store coefficient function for later use (e.g. in boundary integrators)
    regionPermittivity_[regionId] = curCoef;
    
    // Note; in the piezoelectric case we have to multiply by -1
    Double factor = 1.0;
    if ( isPiezoCoupled_ )
      factor = -1.0;
    
    
    if ( isComplex && (regionIsHyst == false)) {
      if( dim_ == 2 ) {
        if (subType_ == "2.5d")
          integ = new BDBInt<Complex, Complex>(new GradientOperator2p5D<FeH1, 2, 1, Complex>(), curCoef, factor, updatedGeo_);
        else
          integ = new BDBInt<Complex,Complex >(new GradientOperator<FeH1,2,1,Complex>(),
                  curCoef, factor, updatedGeo_ );
      } else {
        integ = new BDBInt<Complex,Complex >(new GradientOperator<FeH1,3,1,Complex>(),
                curCoef, factor, updatedGeo_ );
      }
    }
    else {
      if( dim_ == 2 ) {
        if (subType_ == "2.5d"){
          integ = new BDBInt<>(new GradientOperator2p5D<FeH1>(), curCoef, factor, updatedGeo_);
        }else{
          integ = new BDBInt<>(new GradientOperator<FeH1,2> (),
                  curCoef, factor, updatedGeo_ );
        }
      } else{
        integ = new BDBInt<>(new GradientOperator<FeH1,3>(),
                curCoef, factor, updatedGeo_ );
      }
    }
    return integ;
  }
  
  BaseBDBInt* ElecPDE::GetStiffIntegratorInfMap(BaseMaterial* actSDMat, SubTensorType tensorType,
          RegionIdType regionId, PtrCoefFct scalingFactor)
  {
    BaseBDBInt* integ = NULL;
    
    shared_ptr<CoefFunction > curCoef;
    curCoef = actSDMat->GetTensorCoefFnc(ELEC_PERMITTIVITY_TENSOR, tensorType, Global::REAL);

    // store coefficient function for later use (e.g. in boundary integrators)
    regionPermittivity_[regionId] = curCoef;
    PtrCoefFct curCoefScl = CoefFunction::Generate(mp_, Global::REAL, CoefXprTensScalOp(mp_, curCoef, scalingFactor, CoefXpr::OP_MULT_TENSOR));
    
    // Note; in the piezoelectric case we have to multiply by -1
    Double factor = 1.0;
    if (isPiezoCoupled_)
      factor = -1.0;
    
    BaseBOperator* bOp = NULL;
    //Infinite mapping case
    if (dim_ == 2)
    {
      if (subType_ == "2.5d"){
        EXCEPTION("2.5d version with infinite mapping not tested and therefore caught!");
      }else{
        bOp = new ScaledGradientOperator<FeH1, 2, Double>();
      }
    }
    else
      bOp = new ScaledGradientOperator<FeH1, 3, Double>();

    integ = new BDBInt<Double, Double>(bOp, curCoefScl, factor, updatedGeo_);

    
    return integ;
  }
  
  BaseBDBInt* ElecPDE::GetStiffIntegrator(BaseMaterial* actSDMat, SubTensorType tensorType,
          RegionIdType regionId, PtrCoefFct scalingFactor)
  {
    BaseBDBInt* integ = NULL;

    shared_ptr<CoefFunction > curCoef;
    curCoef = actSDMat->GetTensorCoefFnc(ELEC_PERMITTIVITY_TENSOR, tensorType, Global::REAL);


    // store coefficient function for later use (e.g. in boundary integrators)
    regionPermittivity_[regionId] = curCoef;
    PtrCoefFct curCoefScl = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprTensScalOp(mp_, curCoef, scalingFactor, CoefXpr::OP_MULT_TENSOR));

    // Note; in the piezoelectric case we have to multiply by -1
    Double factor = 1.0;
    if (isPiezoCoupled_)
      factor = -1.0;

    BaseBOperator* bOp = NULL;
    if (dim_ == 2)
    {
      if (subType_ == "2.5d")
        bOp = new ScaledGradientOperator2p5D<FeH1, 2, Complex>();
      else
        bOp = new ScaledGradientOperator<FeH1, 2, Complex>();
    }
    else
      bOp = new ScaledGradientOperator<FeH1, 3, Complex>();

    integ = new BDBInt<Complex, Complex>(bOp, curCoefScl, factor, updatedGeo_);

    return integ;
  }

  template<typename DATA_TYPE>
  BiLinearForm* ElecPDE::GetFluxIntegrator(PtrCoefFct scalCoefFunc, PtrCoefFct coefFnc, Double factor,
          BiLinearForm::CouplingDirection cplDir, bool fluxOpA)
  {
    BiLinearForm* integ = NULL;
    BaseBOperator *fluxOp = NULL, *idOp = NULL;
    
    // The flux operator will implement either a scaled gradient operator if a non-zero 'coefFnc' has bee passed,
    // or a gradient operator otherwise
    if (dim_ == 2)
    {
      idOp = new SurfaceIdentityOperator<FeH1, 2, 1>();
      if (coefFnc)
        fluxOp = new SurfaceNormalFluxDensityOperator<FeH1, 2, 1, Complex>(subType_, coefFnc);
      else
        fluxOp = new SurfaceNormalFluxDensityOperator<FeH1, 2, 1, Double>(subType_);
    }
    else
    {
      idOp = new SurfaceIdentityOperator<FeH1, 3, 1>();
      if (coefFnc)
        fluxOp = new SurfaceNormalFluxDensityOperator<FeH1, 3, 1, Complex>(subType_, coefFnc);
      else
        fluxOp = new SurfaceNormalFluxDensityOperator<FeH1, 3, 1, Double>(subType_);
    }
    
    // Check whether we have a du/dn*v or a u*dv/dn bilinear form
    if (fluxOpA)
      integ = new SurfaceNitscheABInt<DATA_TYPE, DATA_TYPE>(fluxOp, idOp, scalCoefFunc, factor, cplDir, updatedGeo_);
    else
      integ = new SurfaceNitscheABInt<DATA_TYPE, DATA_TYPE>(idOp, fluxOp, scalCoefFunc, factor, cplDir, updatedGeo_);
    
    return integ;
  }
  
  template<typename DATA_TYPE>
  BiLinearForm* ElecPDE::GetPenaltyIntegrator(PtrCoefFct scalCoefFunc, Double factor, BiLinearForm::CouplingDirection cplDir)
  {
    BiLinearForm* integ = NULL;
    
    if (dim_ == 2)
    {
      integ = new SurfaceNitscheABInt<DATA_TYPE, DATA_TYPE>(new SurfaceIdentityOperator<FeH1, 2, 1>(),
              new SurfaceIdentityOperator<FeH1, 2, 1>(),
              scalCoefFunc, factor, cplDir, updatedGeo_, false, true);
    }
    else
    {
      integ = new SurfaceNitscheABInt<DATA_TYPE, DATA_TYPE>(new SurfaceIdentityOperator<FeH1, 3, 1>(),
              new SurfaceIdentityOperator<FeH1, 3, 1>(),
              scalCoefFunc, factor, cplDir, updatedGeo_, false, true);
    }
    
    return integ;
  }
  
  void ElecPDE::DefineImpedanceIntegrators() {
    
    // The definition of the impedances is taken from:
    
    // Jian S. Wang and Dale F. Ostergaard:
    // A Finite Element-Electric Circuit Coupled Simulation Method for
    // Piezoelectric Transducer
    // IEEE Ultrasonics Symposium, 1999
    
    // =======================================================================
    // Integrators for electric impedances
    // =======================================================================
    
    for( UInt i = 0; i < impedances_.GetSize(); i++ ) {
      
      REFACTOR;
    }
  }
  
  void ElecPDE::DefineSolveStep() {
		if(isHysteresis_){
			solveStep_ = new SolveStepHyst(*this);
		} else {
			solveStep_ = new StdSolveStep(*this);
		}
  }
  
  void ElecPDE::InitTimeStepping()
  {
    // Check if time integration is defined in XML input
    PtrParamNode transientNode = myParam_->GetParent()->GetParent()->Get("analysis")->Get("transient", ParamNode::PASS);
    PtrParamNode integrationScheme = transientNode->Get("integrationScheme", ParamNode::PASS);

    // Create scheme from XML or use default Trapezoidal
    GLMScheme* scheme = nullptr;
    if (integrationScheme)
    {
      scheme = GetXmlDefinedScheme(integrationScheme);
    }
    else
    {
      scheme = new Trapezoidal(1.0);
    }

    TimeSchemeGLM::NonLinType nlType = (nonLin_ || isHysteresis_) ? TimeSchemeGLM::INCREMENTAL : TimeSchemeGLM::NONE;
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(scheme, 0, nlType));
    feFunctions_[ELEC_POTENTIAL]->SetTimeScheme(myScheme);
  }
  
  void ElecPDE::FinalizeAfterTimeStep() {
  }
  
  void ElecPDE::ReadSpecialBCs( ) 
  {
    //ReadImpedances();
  }
  
  void ElecPDE::ReadImpedances( ) 
  {
    REFACTOR;
  }
  
  void ElecPDE::CalcPolarizationField( shared_ptr<BaseResult> res )
  {
    REFACTOR;
  }
  
  
  
  template <class TYPE>
  void ElecPDE::CalcCharges( shared_ptr<BaseResult> res ) {
    REFACTOR;
    // This will be also moved to the resultFunctor class,
    // as soon as the surface B-operators are defined
  }
  
  void ElecPDE::SetPiezoCoupling()
  {
    isPiezoCoupled_ = true;
  }
  
  void ElecPDE::DefinePrimaryResults() {
    
    shared_ptr<BaseFeFunction> feFct = feFunctions_[ELEC_POTENTIAL];
    
    // Electric Potential
    shared_ptr<ResultInfo> res1( new ResultInfo);
    res1->resultType = ELEC_POTENTIAL;
    
    // check for special subtype 
    if( subType_  != "flatShell" ) {
      res1->dofNames = "";
      res1->unit = MapSolTypeToUnit(ELEC_POTENTIAL);
      res1->definedOn = ResultInfo::NODE;
      res1->entryType = ResultInfo::SCALAR;
    } else {
      
      // check number of composite materials:
      // Currently, we can handle just one electrostatic composite
      // material, as we would have a different number of electric
      // unknowns for different regions
      if( compositeMaterials_.size() > 1 ) {
        WARN("Currently the electrostatic flatShell PDE can only "
                "handle ONE type of composite material!");
      }
      
      Composite & actComp = compositeMaterials_.begin()->second;
      UInt numLaminas = actComp.thickness.GetSize();
      for( UInt i=0; i<numLaminas; i++ ) {
        std::string dofName = "ep";
        dofName += std::to_string(i+1);
        res1->dofNames.Push_back( dofName );
      }
      res1->unit = MapSolTypeToUnit(ELEC_POTENTIAL);
      res1->definedOn = ResultInfo::ELEMENT;
      res1->entryType = ResultInfo::SCALAR;
    }
    feFunctions_[ELEC_POTENTIAL]->SetResultInfo(res1);
    
    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    hdbcSolNameMap_[ELEC_POTENTIAL] = "ground";
    idbcSolNameMap_[ELEC_POTENTIAL] = "potential";
    
    res1->SetFeFunction(feFunctions_[ELEC_POTENTIAL]);
    results_.Push_back( res1 );
    DefineFieldResult( feFunctions_[ELEC_POTENTIAL], res1 );
    
    // Electric RHS Load
    // create new resultDof object
    shared_ptr<ResultInfo> rhs ( new ResultInfo );
    rhs->resultType = ELEC_RHS_LOAD;
    rhs->dofNames = "";
    rhs->unit = MapSolTypeToUnit(ELEC_RHS_LOAD);
    rhs->definedOn = results_[0]->definedOn;
    rhs->entryType = ResultInfo::SCALAR;
    rhsFeFunctions_[ELEC_POTENTIAL]->SetResultInfo(rhs);
    availResults_.insert( rhs );
    DefineFieldResult( rhsFeFunctions_[ELEC_POTENTIAL], rhs );
    
    // Check if an FE space for Lagrange multiplier has been created.
    // If any, define results for Lagrange multiplier in case of p.b.c.
    PtrParamNode lagSpaceNode = infoNode_->Get("feSpaces", ParamNode::PASS)->Get("lagrangeMultiplier", ParamNode::PASS);
    if (lagSpaceNode)
    {
      // <!> This is a hack, because there is no cross points handling
      hdbcSolNameMap_[LAGRANGE_MULT] = "ground";
      idbcSolNameMap_[LAGRANGE_MULT] = "potential";
      //
      shared_ptr<ResultInfo> lagMultElec(new ResultInfo);
      lagMultElec->resultType = LAGRANGE_MULT;
      lagMultElec->dofNames = "";
      lagMultElec->unit = "C/m^2";
      lagMultElec->entryType = ResultInfo::SCALAR;
      lagMultElec->SetFeFunction(feFunctions_[LAGRANGE_MULT]);
      lagMultElec->definedOn = ResultInfo::NODE;
      feFunctions_[LAGRANGE_MULT]->SetResultInfo(lagMultElec);
      
      results_.Push_back(lagMultElec);
      DefineFieldResult(feFunctions_[LAGRANGE_MULT], lagMultElec);
    }
  }
  
  void ElecPDE::DefineNcIntegrators() {
    StdVector< NcInterfaceInfo >::iterator ncIt = ncInterfaces_.Begin(),
            endIt = ncInterfaces_.End();
    for ( ; ncIt != endIt; ++ncIt ) {
      switch (ncIt->type) {
        case NC_MORTAR:
          if (dim_ == 2)
            DefineMortarCoupling<2,1>(ELEC_POTENTIAL, *ncIt);
          else
            DefineMortarCoupling<3,1>(ELEC_POTENTIAL, *ncIt);
          break;
        case NC_NITSCHE:
          EXCEPTION("ncInterface of Nitsche type is not implemented for ElecPDE");
          break;
        default:
          EXCEPTION("Unknown type of ncInterface");
          break;
      }
    }
  }
  
  void ElecPDE::DefinePostProcResults() {
    
    shared_ptr<BaseFeFunction> feFct = feFunctions_[ELEC_POTENTIAL];
    bool is2p5 = (subType_ == "2.5d");
    
    // === ELECTRIC FIELD INTENSITY ===
    // Electric field intensity \bm{E} = -\nabla \bm{\Phi}
    shared_ptr<ResultInfo> ef ( new ResultInfo );
    ef->resultType = ELEC_FIELD_INTENSITY;
    ef->SetVectorDOFs(dim_, isaxi_, is2p5);
    ef->unit = MapSolTypeToUnit(ELEC_FIELD_INTENSITY);
    ef->definedOn = ResultInfo::ELEMENT;
    ef->entryType = ResultInfo::VECTOR;
    shared_ptr<CoefFunctionFormBased> eFunc;

    if( isComplex_ ) {
      eFunc.reset(new CoefFunctionBOp<Complex>(feFct, ef, -1.0));

    } else {
      eFunc.reset(new CoefFunctionBOp<Double>(feFct, ef, -1.0));
    }
    DefineFieldResult( eFunc, ef );
    stiffFormCoefs_.insert(eFunc);
    
    // === SURF ELECTRIC FIELD INTENSITY ===
    // Background tangential component(s) of electric field on boundary
    //  required in case of hysteresis in combination with fieldParallel BC
    //  (in this case the hystOperator has to be evaluated with the tangential 
    //    component of E on the boundary)
    shared_ptr<ResultInfo> efs ( new ResultInfo );
    efs->resultType = ELEC_FIELD_INTENSITY_SURF;
    efs->SetVectorDOFs(dim_, isaxi_, is2p5);
    efs->unit = MapSolTypeToUnit(ELEC_FIELD_INTENSITY_SURF);
    efs->definedOn = ResultInfo::SURF_ELEM;
    efs->entryType = ResultInfo::VECTOR;
    shared_ptr<CoefFunctionSurf> eFuncS;
    
    eFuncS.reset(new CoefFunctionSurf(false, 1.0, efs));
    DefineFieldResult(eFuncS, efs);
    // original function required, not surface version
    surfCoefFcts_[eFuncS] = eFunc;
    
    // === ELECTRIC POLARIZATION ===
    // check if at least one region has hysteretic behavior
    if ( isHysteresis_){
      // this is the point where the actual storage is created
      hysteresisCoefs_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_,1,isComplex_));
      hysteresisPolarization_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_,1,isComplex_));
      shared_ptr<ResultInfo> elecP ( new ResultInfo );
      elecP->resultType = ELEC_POLARIZATION;
      elecP->SetVectorDOFs(dim_, isaxi_, is2p5);
      elecP->unit = MapSolTypeToUnit(ELEC_POLARIZATION);
      elecP->definedOn = ResultInfo::ELEMENT;
      elecP->entryType = ResultInfo::VECTOR;
      DefineFieldResult( hysteresisPolarization_, elecP );
      availResults_.insert( elecP );
      
      // for piezoCoupling, the hysteresis coef function has to be initialized
      // in order to assemble the required coupling tensors
      // however, the DefineIntegrators functions of the coupled pde is
      // called prior to the DefineIntegrators functions of the single pdes
      // (which is where the hystCoefFunctions were created so far
      // therefore, this helper functions is used to create the hystCoefFunctions
      //	before the call to DefineIntegrators
      InitHystCoefs();
    }
    //=== PERMITIVITY ==
    //TODO Maybe there is a more elegant way to do this. But for now it works
    shared_ptr<ResultInfo> perm (new ResultInfo);
    perm->resultType = ELEC_ELEM_PERMITTIVITY;
    perm->dofNames = "";
    perm->unit = MapSolTypeToUnit(ELEC_ELEM_PERMITTIVITY);
    perm->definedOn = ResultInfo::ELEMENT;
    perm->entryType = ResultInfo::SCALAR;
    shared_ptr<CoefFunctionFormBased> perm_coef;
    if(analysistype_== MULTIHARMONIC){
      perm_coef.reset(new CoefFunctionHomogenization<Complex, App::ELEC>(feFct));
    } else {
      perm_coef.reset(new CoefFunctionHomogenization<double, App::ELEC>(feFct));
    }
    DefineFieldResult(perm_coef , perm);
    stiffFormCoefs_.insert(perm_coef);
    
    // === ELECTRIC FLUX DENSITY ===
    shared_ptr<ResultInfo> flux ( new ResultInfo );
    flux->resultType = ELEC_FLUX_DENSITY;
    flux->SetVectorDOFs(dim_, isaxi_, is2p5);
    flux->unit = MapSolTypeToUnit(ELEC_FLUX_DENSITY);
    flux->definedOn = ResultInfo::ELEMENT;
    flux->entryType = ResultInfo::VECTOR;
    
    shared_ptr<CoefFunctionFormBased> fluxFunc;
    if( isComplex_ ) {
      fluxFunc.reset(new CoefFunctionFlux<Complex>(feFct, flux, -1.0));
    } else {
      fluxFunc.reset(new CoefFunctionFlux<Double>(feFct, flux, -1.0));
    }
    stiffFormCoefs_.insert(fluxFunc);
    
    PtrCoefFct combinedFlux;
    if ( isHysteresis_){
      /*
       * in case of hysteresis, fluxFunc will only contain the value eps0*E
       * -> we have to add polarization afterwards
       */
      if( isComplex_ ) {
        combinedFlux = CoefFunction::Generate(mp_,Global::COMPLEX,CoefXprBinOp(mp_,fluxFunc,hysteresisPolarization_,CoefXpr::OP_ADD));
      } else {
        combinedFlux = CoefFunction::Generate(mp_,Global::REAL,CoefXprBinOp(mp_,fluxFunc,hysteresisPolarization_,CoefXpr::OP_ADD));
      }
      DefineFieldResult( combinedFlux, flux );
    } else if(nonLin_){
      // for all other nonlinear cases
      // !!!
      // Be aware that the hysteresis which are build in via CoefFunctionMaterialModel does not use the isHysteresis_ flag!
      // !!!
      DefineFieldResult( nlFluxCoef_, flux );
    }else {
      DefineFieldResult( fluxFunc, flux );
    }
    
    // === ELECTRIC SURFACE CHARGE DENSITY ===
    // Electric surface charge density q_{s} = -\varepsilon \nabla \bm{\Phi} \cdot \bm{n}
    shared_ptr<ResultInfo> chargeD(new ResultInfo);
    chargeD->resultType = ELEC_SURFACE_CHARGE_DENSITY;
    chargeD->dofNames = "";
    chargeD->unit = MapSolTypeToUnit(ELEC_SURFACE_CHARGE_DENSITY);
    chargeD->definedOn = ResultInfo::SURF_ELEM;
    chargeD->entryType = ResultInfo::SCALAR;
    availResults_.insert( chargeD );
    // the coefficient function is defined later
    // Note: The positive normal direction in this case is defined as the
    //       inward facing one. 
    shared_ptr<CoefFunctionSurf> sChargeDens(new CoefFunctionSurf(true, -1.0, chargeD));
    DefineFieldResult( sChargeDens, chargeD);
    if(!isHysteresis_){
      surfCoefFcts_[sChargeDens] = fluxFunc;
    } else {
      surfCoefFcts_[sChargeDens] = combinedFlux;
    }
    
    // === TOTAL ELECTRIC CHARGE ===
    // Total electric charge Q = \int_{\Gamma} q_{s} \mathrm{d} \Gamma
    shared_ptr<ResultInfo> charge(new ResultInfo);
    charge->resultType = ELEC_CHARGE;
    charge->dofNames = "";
    charge->unit = MapSolTypeToUnit(ELEC_CHARGE);
    charge->definedOn = ResultInfo::SURF_REGION;
    charge->entryType = ResultInfo::SCALAR;
    availResults_.insert( charge );
    // build result functor for integration
    shared_ptr<ResultFunctor> chargeFunc;
    if( isComplex_ ) {
      chargeFunc.reset(new ResultFunctorIntegrate<Complex>(sChargeDens, feFct, charge ) );
    } else {
      chargeFunc.reset(new ResultFunctorIntegrate<Double>(sChargeDens, feFct, charge ) );
    }
    resultFunctors_[ELEC_CHARGE] = chargeFunc;
    
    // === ELECTRIC ENERGY DENSITY ===
    // Electric energy density for scalar permitiitvity e = \frac{1} {2} \varepsilon \bm{E} \cdo \bm{E}
    shared_ptr<ResultInfo> ed ( new ResultInfo );
    ed->resultType = ELEC_ENERGY_DENSITY;
    ed->dofNames = "";
    ed->unit = MapSolTypeToUnit(ELEC_ENERGY_DENSITY);
    ed->definedOn = ResultInfo::ELEMENT;
    ed->entryType = ResultInfo::SCALAR;
    shared_ptr<CoefFunctionFormBased> edFunc;

    // for both BdBKernel and EnergyResultFunctor, we need to apply the -1 factor
    // to get right sign in the results (even though the energy results are not really usable in the coupled case as they neglect the influnce of the coupled pde)
    Double factor = 1.0;
    if ( isPiezoCoupled_ ){
      factor = -1.0;
    }

    if( isComplex_ ) {
      edFunc.reset(new CoefFunctionBdBKernel<Complex>(feFct, factor*0.5));
    } else {
      edFunc.reset(new CoefFunctionBdBKernel<Double>(feFct, factor*0.5));
    }

    DefineFieldResult( edFunc, ed );
    stiffFormCoefs_.insert(edFunc);

    // Electric energy
    // Electric energy E = \frac{1} {2} \varepsilon \int_{\Omega} \bm{E} \cdot \bm{E} \mathrm{d} \Omega
    shared_ptr<ResultInfo> energy( new ResultInfo );
    energy->resultType = ELEC_ENERGY;
    energy->definedOn = ResultInfo::REGION;
    energy->entryType = ResultInfo::SCALAR;
    energy->dofNames = "";
    energy->unit = MapSolTypeToUnit(ELEC_ENERGY);
    availResults_.insert( energy );
    shared_ptr<ResultFunctor> energyFunc;
    if( isComplex_ ) {
      energyFunc.reset(new EnergyResultFunctor<Complex>(feFct, energy,factor*0.5));
    } else {
      energyFunc.reset(new EnergyResultFunctor<Double>(feFct, energy,factor*0.5));
    }
    resultFunctors_[ELEC_ENERGY] = energyFunc;
    stiffFormFunctors_.insert(energyFunc);


    StdVector<std::string> vecDofNames;
    if( ptGrid_->GetDim() == 3 ) {
      vecDofNames = "x", "y", "z";
    } else {
      if( ptGrid_->IsAxi() ) {
        vecDofNames = "r", "z";
      } else {
        vecDofNames = "x", "y";
      }
    }

    // === ELECTROSTATC FORCE DENSITY ===
    shared_ptr<ResultInfo> efd(new ResultInfo);
    efd->resultType = ELEC_FORCE_DENSITY;
    //efd->SetVectorDOFs(dim_, isaxi_, is2p5);
    efd->dofNames= vecDofNames;
    efd->unit = MapSolTypeToUnit(ELEC_FORCE_DENSITY);
    efd->definedOn = ResultInfo::SURF_ELEM;
    efd->entryType = ResultInfo::VECTOR;
    shared_ptr<CoefFunctionSurf> efdFuncS;
    //factor -1 because positive normal direction points into mechanic region
    efdFuncS.reset(new CoefFunctionSurf(true, -1.0, efd));
    DefineFieldResult(efdFuncS, efd);
    surfCoefFcts_[efdFuncS] = edFunc;

    // === ELECTROSTATC FORCE ===
    shared_ptr<ResultInfo> elForce(new ResultInfo);
    elForce->resultType = ELEC_FORCE;
    elForce->dofNames= vecDofNames;
    elForce->unit = MapSolTypeToUnit(ELEC_FORCE);
    elForce->definedOn = ResultInfo::SURF_REGION;
    elForce->entryType = ResultInfo::VECTOR;
    availResults_.insert( elForce );
    // build result functor for integration
    shared_ptr<ResultFunctor> elForceFunc;
    if( isComplex_ ) {
      elForceFunc.reset(new ResultFunctorIntegrate<Complex>(efdFuncS, feFct, elForce ) );
    } else {
      elForceFunc.reset(new ResultFunctorIntegrate<Double>(efdFuncS, feFct, elForce ) );
    }
    resultFunctors_[ELEC_FORCE] = elForceFunc;
  }
  
  
  std::map<SolutionType, shared_ptr<FeSpace> >
  ElecPDE::CreateFeSpaces(const std::string& formulation, PtrParamNode infoNode) {
    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    if(formulation == "default" || formulation == "H1"){
      PtrParamNode potSpaceNode = infoNode->Get("elecPotential");
      crSpaces[ELEC_POTENTIAL] =
              FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1, ptGrid_);
      crSpaces[ELEC_POTENTIAL]->Init(solStrat_);
      
      // Maybe we can use SingelPDE::CreateNcFeSpaces instead?
      CoupledField::ParamNodeList blochListMortar = myParam_->Get("bcsAndLoads")->GetListByVal("blochPeriodic", "formulation", "Mortar");
      if (blochListMortar.GetSize())
      {
        // Create FE-Space for Lagrange multiplier
        PtrParamNode lagSpaceNode = infoNode->Get("lagrangeMultiplier");
        crSpaces[LAGRANGE_MULT] = FeSpace::CreateInstance(myParam_, lagSpaceNode, FeSpace::H1, ptGrid_);
        crSpaces[LAGRANGE_MULT]->SetLagrSurfSpace();
        crSpaces[LAGRANGE_MULT]->Init(solStrat_);
      }
    }else{
      EXCEPTION("The formulation " << formulation << "of electric PDE is not known!");
    }
    return crSpaces;
  }
  
  template BiLinearForm* ElecPDE::GetPenaltyIntegrator<Double>(PtrCoefFct, Double, BiLinearForm::CouplingDirection);
  template BiLinearForm* ElecPDE::GetPenaltyIntegrator<Complex>(PtrCoefFct, Double, BiLinearForm::CouplingDirection);
  template BiLinearForm* ElecPDE::GetFluxIntegrator<Double>(PtrCoefFct, PtrCoefFct, Double, BiLinearForm::CouplingDirection, bool);
  template BiLinearForm* ElecPDE::GetFluxIntegrator<Complex>(PtrCoefFct, PtrCoefFct, Double, BiLinearForm::CouplingDirection, bool);
}
