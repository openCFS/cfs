#include "AcousticPDE.hh"

#include "General/defs.hh"

#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

//new integrator concept
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/BiLinForms/ADBInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/BDUInt.hh"
#include "Forms/LinForms/KXInt.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/IdentityOperatorNormal.hh"
#include "Forms/Operators/ConvectiveOperator.hh"
#include "Forms/Operators/SurfaceOperators.hh"
#include "Forms/Operators/DivOperator.hh"

#include "FeBasis/FeFunctions.hh"
#include "Utils/StdVector.hh"
#include "FeBasis/H1/FeSpaceH1Nodal.hh"

#include "Domain/Results/ResultFunctor.hh"

#include "Driver/Assemble.hh"
#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/CoefFunction/CoefFunctionCompound.hh"
#include "Domain/CoefFunction/CoefFunctionConst.hh"
#include "Domain/CoefFunction/CoefFunctionMulti.hh"
#include "Domain/CoefFunction/CoefFunctionPML.hh"
#include "Domain/CoefFunction/CoefFunctionCurvilinearPML.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/CoefFunction/CoefFunctionImpedanceModel.hh"
#include "Domain/CoefFunction/CoefFunctionOpt.hh"
#include "Domain/Mesh/NcInterfaces/BaseNcInterface.hh"
#include "Domain/Mesh/NcInterfaces/MortarInterface.hh"

#include <cmath>


#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"
#include "Materials/AcousticMaterial.hh"

namespace CoupledField{

   DEFINE_LOG(acousticpde, "pde.acoustic")


  AcousticPDE::AcousticPDE( Grid* aGrid, PtrParamNode paramNode,
                            PtrParamNode infoNode,
                            shared_ptr<SimState> simState, Domain* domain)
              : SinglePDE( aGrid, paramNode, infoNode, simState, domain ){

    pdename_           = "acoustic";
    pdematerialclass_  = ACOUSTIC;
    nonLin_            = false;
    isMechCoupled_     = false;

    //! Always use total Lagrangian formulation 
    updatedGeo_        = false;
    isTimeDomPML_      = false;
    isAPML_            = false;
    complexFluidFormulation_ = false;
    isOnlyOneMaterial_ = true;
    timeDomainEqFluidFormulation_ = false;

    //! set the PDE formulation
    std::string pdeFormulation = myParam_->Get("formulation")->As<std::string>();
    if (pdeFormulation == "acouPressureSOSatLaplace") {
      sosAtLaplace_ = true;
      formulation_ = ACOU_PRESSURE;
    } else if (pdeFormulation == "acouPressure") {
      sosAtLaplace_ = false;
      formulation_ = ACOU_PRESSURE;
    } else if (pdeFormulation == "acouPotential") {
      sosAtLaplace_ = false;
      formulation_ = ACOU_POTENTIAL;
    } else {
      EXCEPTION("Unknown PDE formulation '" << pdeFormulation << "' for AcousticPDE. " 
                << "Possible formulations are: 'acouPressure', 'acouPotential', 'acouPressureSOSatLaplace'.")
    }
  }

  std::map<SolutionType, shared_ptr<FeSpace> > AcousticPDE::CreateFeSpaces( const std::string&  formulation,
                  PtrParamNode infoNode )
  {
    if(this->analysistype_ == STATIC)
      EXCEPTION("No STATIC analysis in AcousticPDE");

    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    if(formulation == "default" || formulation == "H1"){
      std::string form = SolutionTypeEnum.ToString(formulation_);
      PtrParamNode potSpaceNode = infoNode->Get(form);
      crSpaces[formulation_] =
        FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1, ptGrid_);
      crSpaces[formulation_]->Init(solStrat_);
      
      // ===================================
      // blochPeriodic Mortar
      // ===================================
      CoupledField::PtrParamNode bcsList = myParam_->Get("bcsAndLoads", ParamNode::ActionType::PASS);
      if (bcsList != nullptr && bcsList->GetListByVal("blochPeriodic", "formulation", "Mortar").GetSize())
      {
        // Create FE-Space for Lagrange multiplier
        PtrParamNode lagSpaceNode = infoNode->Get("lagrangeMultiplier");
        crSpaces[LAGRANGE_MULT] = FeSpace::CreateInstance(myParam_, lagSpaceNode, FeSpace::H1, ptGrid_);
        crSpaces[LAGRANGE_MULT]->SetLagrSurfSpace();
        crSpaces[LAGRANGE_MULT]->Init(solStrat_);
      }
    }else{
      EXCEPTION("The formulation " << formulation << "of acoustic PDE is not known!");
    }

    // ===================================
    // Check for transient PML
    // ===================================
    if(this->analysistype_ == TRANSIENT && isTimeDomPML_){
      //now define the additional unknowns

      if(dim_==3 && !this->isAPML_){
        PtrParamNode scalarpml = infoNode->Get("TransientPMLScalarAuxVar");
        crSpaces[ACOU_PMLAUXSCALAR] =
            FeSpace::CreateInstance(myParam_,scalarpml,FeSpace::H1, ptGrid_);
        crSpaces[ACOU_PMLAUXSCALAR]->Init(solStrat_);
      }

      PtrParamNode vectorPML = infoNode->Get("TransientPMLVectorAuxVar");
      crSpaces[ACOU_PMLAUXVEC] =
          FeSpace::CreateInstance(myParam_,vectorPML,FeSpace::H1, ptGrid_);
      crSpaces[ACOU_PMLAUXVEC]->Init(solStrat_);
    }

    // ===================================
    // Check for complex fluid or TDEF regions
    // ===================================
    RegionIdType actRegion;
    bool evalRationalFunc = false;
    std::map<RegionIdType, BaseMaterial *>::iterator it;
    for (it = materials_.begin(); it != materials_.end(); it++) {
      // check for complex fluid formulation
      actRegion = it->first;
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);
      PtrParamNode curRegNode =
          myParam_->Get("regionList")->GetByVal("region", "name", regionName.c_str());
      if (curRegNode->Get("complexFluid")->As<std::string>() == "yes") {
        complexFluidFormulation_ = true;
        complexMatData_[actRegion] = true;  // this is important for optimization
        isMaterialComplex_ = true;
        if (this->analysistype_ != HARMONIC)
          EXCEPTION("Complex fluid region just allowed in harmonic analysis");
        // need an acoustic pressure formulation
        if (formulation_ != ACOU_PRESSURE)
          EXCEPTION("Complex fluid requires acoustic pressure formulation");
      }

      // check for time domain equivalent fluid (TDEF) formulation
      if (curRegNode->Get("timeDomainEqFluid")->As<std::string>() == "yes" && this->analysistype_ == TRANSIENT) {
        timeDomainEqFluidFormulation_ = true;
        // need an acoustic pressure formulation
        if (formulation_ != ACOU_PRESSURE)
          EXCEPTION("Time domain equivalent fluid formulation requires acoustic pressure formulation");
      }

      if (this->analysistype_ == HARMONIC && curRegNode->Get("timeDomainEqFluid")->As<std::string>() == "yes") {
        EXCEPTION("Time domain equivalent fluid formulation only possible for transient simulation. To use the rational function approximation for the complex fluid parameter, set useRationalMatApproximation to yes.");
      }

      // check if complex fluid parameters (for the complex fluid formulation) are provided using a rational function
      if (curRegNode->Get("useRationalMatApproximation")->As<std::string>() == "yes" && complexFluidFormulation_ == true) {
        std::cout << "Complex fluid parameters provided by rational function approximation. " << "\n";
        evalRationalFunc = true;
      }
    }

    // ===================================
    // Check for rational function approximation of the equivalent/complex fluid parameters and get no.
    // ===================================
    if ((this->analysistype_ == TRANSIENT && timeDomainEqFluidFormulation_) || (this->analysistype_ == HARMONIC && evalRationalFunc)) {
      // check for the actual size of the auxiliary variables
      TDEFpoleNumber_.realC.Resize(regions_.GetSize());
      TDEFpoleNumber_.complC.Resize(regions_.GetSize());
      TDEFpoleNumber_.realV.Resize(regions_.GetSize());
      TDEFpoleNumber_.complV.Resize(regions_.GetSize());
      isTDEFReg_.Resize(regions_.GetSize());

      RegionIdType actRegion;
      for (UInt iRegion = 0; iRegion < regions_.GetSize(); iRegion++) {
        actRegion = regions_[iRegion];
        std::string regionName = ptGrid_->GetRegion().ToString(actRegion);
        PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region", "name", regionName.c_str());

        if (curRegNode->Get("timeDomainEqFluid")->As<std::string>() == "yes" || curRegNode->Get("useRationalMatApproximation")->As<std::string>() == "yes") {
          LocPointMapped lpm;

          Vector<Double> vecAC;
          Vector<Double> vecBC;
          Vector<Double> vecAV;
          Vector<Double> vecBV;

          // get the vector of residues
          PtrCoefFct coefVecAC = materials_[actRegion]->GetVectorCoefFnc(ACOU_TDEF_INVBLK_A, Global::REAL);
          PtrCoefFct coefVecBC = materials_[actRegion]->GetVectorCoefFnc(ACOU_TDEF_INVBLK_B, Global::REAL);
          PtrCoefFct coefVecAV = materials_[actRegion]->GetVectorCoefFnc(ACOU_TDEF_INVDENS_A, Global::REAL);
          PtrCoefFct coefVecBV = materials_[actRegion]->GetVectorCoefFnc(ACOU_TDEF_INVDENS_B, Global::REAL);

          coefVecAC->GetVector(vecAC, lpm);
          coefVecBC->GetVector(vecBC, lpm);
          coefVecAV->GetVector(vecAV, lpm);
          coefVecBV->GetVector(vecBV, lpm);

          // assign the size of the parameter vector (no. real or complex poles)
          TDEFpoleNumber_.realC[iRegion] = vecAC.GetSize();
          TDEFpoleNumber_.complC[iRegion] = vecBC.GetSize();
          TDEFpoleNumber_.realV[iRegion] = vecAV.GetSize();
          TDEFpoleNumber_.complV[iRegion] = vecBV.GetSize();
        }
        else {
          // this region has no TDEF material defined, set size to 0
          TDEFpoleNumber_.realC[iRegion] = 0;
          TDEFpoleNumber_.complC[iRegion] = 0;
          TDEFpoleNumber_.realV[iRegion] = 0;
          TDEFpoleNumber_.complV[iRegion] = 0;
        }

        if (TDEFpoleNumber_.realC[iRegion] == 0 && TDEFpoleNumber_.complC[iRegion] == 0 && TDEFpoleNumber_.realV[iRegion] == 0 && TDEFpoleNumber_.complV[iRegion] == 0) {
          isTDEFReg_[iRegion] = false;
        }
        else {
          isTDEFReg_[iRegion] = true;
        }
      }
    }

    if (this->analysistype_ == TRANSIENT && timeDomainEqFluidFormulation_) {
      // define the additional unknowns (phiC, psiC, phiV, psiV)
      std::cout << std::endl
                << "Applying the TDEF formulation." << std::endl
                << std::endl;
      PtrParamNode spaceNode;
      for (unsigned int i = 0; i < TDEFpoleNumber_.realC.Max(); i++) {
        spaceNode = infoNode->Get(SolutionTypeEnum.ToString((SolutionType)(ACOU_TDEF_PHI_C_1 + i)));
        crSpaces[(SolutionType)(ACOU_TDEF_PHI_C_1 + i)] = FeSpace::CreateInstance(myParam_, spaceNode, FeSpace::H1, ptGrid_);
        crSpaces[(SolutionType)(ACOU_TDEF_PHI_C_1 + i)]->Init(solStrat_);
      }
      for (unsigned int i = 0; i < TDEFpoleNumber_.complC.Max(); i++) {
        spaceNode = infoNode->Get(SolutionTypeEnum.ToString((SolutionType)(ACOU_TDEF_PSI_C_1 + i)));
        crSpaces[(SolutionType)(ACOU_TDEF_PSI_C_1 + i)] = FeSpace::CreateInstance(myParam_, spaceNode, FeSpace::H1, ptGrid_);
        crSpaces[(SolutionType)(ACOU_TDEF_PSI_C_1 + i)]->Init(solStrat_);
      }
      for (unsigned int i = 0; i < TDEFpoleNumber_.realV.Max(); i++) {
        spaceNode = infoNode->Get(SolutionTypeEnum.ToString((SolutionType)(ACOU_TDEF_PHI_V_1 + i)));
        crSpaces[(SolutionType)(ACOU_TDEF_PHI_V_1 + i)] = FeSpace::CreateInstance(myParam_, spaceNode, FeSpace::H1, ptGrid_);
        crSpaces[(SolutionType)(ACOU_TDEF_PHI_V_1 + i)]->Init(solStrat_);
      }
      for (unsigned int i = 0; i < TDEFpoleNumber_.complV.Max(); i++) {
        spaceNode = infoNode->Get(SolutionTypeEnum.ToString((SolutionType)(ACOU_TDEF_PSI_V_1 + i)));
        crSpaces[(SolutionType)(ACOU_TDEF_PSI_V_1 + i)] = FeSpace::CreateInstance(myParam_, spaceNode, FeSpace::H1, ptGrid_);
        crSpaces[(SolutionType)(ACOU_TDEF_PSI_V_1 + i)]->Init(solStrat_);
      }
    }
    return crSpaces;
  }

  void AcousticPDE::ReadDampingInformation()
  {
    // get regions of current PDE
    ParamNodeList regionParamNodes = myParam_->Get("regionList")->GetChildren();
    // corresponding region id
    RegionIdType actRegionId;
    // corresponding region name and damping id
    std::string actRegionName, actDampingId;
    // try to get the dampingList and return if it is not specified
    PtrParamNode dampListNode = myParam_->Get("dampingList", ParamNode::PASS);
    if (dampListNode) {
      // get the individual damping nodes
      ParamNodeList dampNodes = dampListNode->GetChildren();
      // map of damping ids from the xml and corresponding damping types
      std::map<std::string, DampingType> idDampType;

      // Run over all region param nodes and assign the required damping information
      for (UInt iRegion = 0; iRegion < regionParamNodes.GetSize(); ++iRegion) {
        regionParamNodes[iRegion]->GetValue("name", actRegionName);
        regionParamNodes[iRegion]->GetValue("dampingId", actDampingId);

        // pass if no damping is defined
        if (actDampingId == "")
          continue;

        // parse region id from region name
        actRegionId = ptGrid_->GetRegion().Parse(actRegionName);

        // now, read the damping information from the dampingList
        for (UInt iChild = 0; iChild < dampNodes.GetSize(); ++iChild) {
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
        if (dampingList_[actRegionId] == RAYLEIGH) {
          RaylDampingData actRaylCoeffs;
          materials_[actRegionId]->GetRayleighCoeffStrings(actRaylCoeffs.alpha, actRaylCoeffs.beta);
          regionRaylDamping_[actRegionId] = actRaylCoeffs;
        }
        else if (dampingList_[actRegionId] == ADAPTED_LOSS_TANGENS_DELTA) {
          if (!(analysistype_ == BasePDE::HARMONIC))
            EXCEPTION("Adapted loss tangent delta damping is only allowed for harmonic analysis.");
          RaylDampingData actRaylCoeffs;
          materials_[actRegionId]->GetFreqAdaptedRayleighCoeffStrings(actRaylCoeffs.alpha, actRaylCoeffs.beta);
          regionRaylDamping_[actRegionId] = actRaylCoeffs;
        }
        else if (dampingList_[actRegionId] == GLOBAL_RAYLEIGH) {
          EXCEPTION("Global Rayleigh damping is not yet implemented.");
          if (dampNodes.GetSize() != 1)
            EXCEPTION("Global Rayleigh damping does not allow for additional damping nodes defined.");
          RaylDampingData actRaylCoeffs;
          actRaylCoeffs.alpha = dampNodes[0]->Get("alpha")->As<std::string>();
          actRaylCoeffs.beta = dampNodes[0]->Get("beta")->As<std::string>();
          regionRaylDamping_[actRegionId] = actRaylCoeffs;
        }

        // set flag to compute extra integrators for transient PML
        if (dampingList_[actRegionId] == PML && analysistype_ == BasePDE::TRANSIENT) {
          isTimeDomPML_ = true;
        }
      }
    }
  }

  void AcousticPDE::DefineIntegrators() {

    RegionIdType actRegion;
    // type of geometry
    isaxi_ = ptGrid_->IsAxi();

    // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial *>::iterator it;
    shared_ptr<FeSpace> mySpace = feFunctions_[formulation_]->GetFeSpace();

    LOG_DBG(acousticpde) << "DefineIntegrators BEGIN" << "\n";

    for (UInt iRegion = 0; iRegion < regions_.GetSize(); iRegion++) {
      actRegion = regions_[iRegion];

      // Get current region name
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);

      // create new entity list
      shared_ptr<ElemList> actSDList(new ElemList(ptGrid_));
      actSDList->SetRegion(actRegion);

      // --- Set the FE ansatz for the current region ---
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region", "name", regionName.c_str());
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();
      std::string useRationalAppr = curRegNode->Get("useRationalMatApproximation")->As<std::string>();
      mySpace->SetRegionApproximation(actRegion, polyId,integId);

      //=======================================================================
      // Generate coefficient functions 4 (Speed of Sound)
      //=======================================================================
      PtrCoefFct c0;
      PtrCoefFct dens;
      PtrCoefFct blk;
      PtrCoefFct freqCoef;
      PtrCoefFct constOne = CoefFunction::Generate(mp_, Global::REAL, "1.0");

      if (complexFluidFormulation_ && useRationalAppr == "no") {
        dens = materials_[actRegion]->GetScalCoefFnc(DENSITY, Global::COMPLEX);
        blk = materials_[actRegion]->GetScalCoefFnc(ACOU_BULK_MODULUS, Global::COMPLEX);
      }
      else if (complexFluidFormulation_ && useRationalAppr == "yes") {
        std::cout << "Evaluating the rational function approximation of the bulk modulus and density at each frequency." << std::endl;
         // get the current frequency
        freqCoef = CoefFunction::Generate(mp_, Global::REAL, "f");
        EvalTDEFRationalFncs(iRegion, freqCoef, TDEFcoeffs_);
        dens = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, constOne, invTDEFDens_, CoefXpr::OP_DIV));
        blk = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, constOne, invTDEFBlk_, CoefXpr::OP_DIV));
      }
      else if (timeDomainEqFluidFormulation_) {
        // Attention:
        // in case of the TDEF formulation, the high frequency limit is set for ACOU_BULK_MODULUS and DENSITY (for user input checks)
        // We also created ACOU_TDEF_INVDENS_CONST and ACOU_TDEF_INVBLK_CONST, which have the inverse values (we currently don't use them)
        PtrCoefFct densInvConst;
        PtrCoefFct blkInvConst;
        densInvConst = materials_[actRegion]->GetScalCoefFnc(DENSITY, Global::REAL);
        blkInvConst = materials_[actRegion]->GetScalCoefFnc(ACOU_BULK_MODULUS, Global::REAL);
        dens = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, constOne, densInvConst, CoefXpr::OP_DIV));
        blk = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, constOne, blkInvConst, CoefXpr::OP_DIV));
      }
      else {
        dens = materials_[actRegion]->GetScalCoefFnc(DENSITY, Global::REAL);
        blk = materials_[actRegion]->GetScalCoefFnc(ACOU_BULK_MODULUS, Global::REAL);
      }

      LOG_DBG(acousticpde) << "DefineIntegrators : dens   = " << dens->ToString() << "\n";
      LOG_DBG(acousticpde) << "DefineIntegrators : blk    = " << blk->ToString() << "\n";

      // ====================================================================
      // check for temperature field, which effects the speed of sound
      // ====================================================================
      std::string tempId = curRegNode->Get("temperatureId")->As<std::string>();
      PtrCoefFct regionTemp;

      if (tempId != "") {
        std::set<UInt> definedDofs;
        // just real valued temperature allowed
        shared_ptr<ResultInfo> tempInfo = GetResultInfo(HEAT_MEAN_TEMPERATURE);
        // Add the region information
        PtrParamNode tempNode = myParam_->Get("temperatureList")->GetByVal("temperature", "name", tempId.c_str());

        ReadUserFieldValues(actSDList, tempNode, tempInfo->dofNames, tempInfo->entryType, false, regionTemp, definedDofs, updatedGeo_);
        meanTemperatureCoef_->AddRegion(actRegion, regionTemp);
      }

      // compute speed of sound
      ComputeSOS(c0, dens, blk, regionTemp, tempId);
      // store coefficient functions
      matCoefs_[ELEM_DENSITY]->AddRegion(actRegion, dens);
      matCoefs_[ACOU_ELEM_SPEED_OF_SOUND]->AddRegion(actRegion, c0);

      // if pde couples with mechanic, we have to multiply the density by -1
      PtrCoefFct mechAcouFactor;
      CalcMechAcouFac(mechAcouFactor, dens);

      // basic coeff-functions for mass and stiffness matrix
      PtrCoefFct coeffM, coeffK;

      if (formulation_ == ACOU_PRESSURE && sosAtLaplace_ == true) {
        if (complexFluidFormulation_ || timeDomainEqFluidFormulation_)
          EXCEPTION("A complex fluid and sosAtLaplace-formulation not allowed!!");

        // pressure formulation with temperature depend speed of sound
        coeffM = CoefFunction::Generate(mp_, Global::REAL, "1.0");;
        coeffK = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, c0, c0, CoefXpr::OP_MULT));
      }
      else {
        if (complexFluidFormulation_) {
          // in this case c0 is actually 1/c0^2!!
          //! coeffM: 1/compressionModulus
          //! coeffK = 1/density
          coeffM = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, mechAcouFactor, blk, CoefXpr::OP_DIV));
          /// coeffM: 1/density
          coeffK = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, mechAcouFactor, dens, CoefXpr::OP_DIV));
        }
        else if (timeDomainEqFluidFormulation_) {
          // In the case of TDEF, dens and blk are already the inverse of the high-frequency limit of the rational function approximation.
          // They are assigned as coefficients to mass and stiffness matrix according to the TDEF formulation (see Diss Maurerlehner 2023).
          coeffK = dens;
          coeffM = blk;
        }
        else {
          coeffK = mechAcouFactor;
          // build coefficient for mass matrix as (factor / (c0*c0))
          PtrCoefFct constValOne = CoefFunction::Generate(mp_, Global::REAL, "1.0");
          PtrCoefFct coeffa, coeffb;
          coeffa = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, constValOne, c0, CoefXpr::OP_DIV));
          coeffb = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, coeffa, coeffa, CoefXpr::OP_MULT));
          coeffM = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, mechAcouFactor, coeffb, CoefXpr::OP_MULT));
        }
      }

      LOG_DBG(acousticpde) << "DefineIntegrators Fluid: coeffK = " << coeffK->ToString() << "\n";
      LOG_DBG(acousticpde) << "DefineIntegrators Fluid: coeffM = " << coeffM->ToString() << "\n";

      // Apply topology optimization factors
      if(domain->HasDesign())
      {
        if(!complexFluidFormulation_)
          throw Exception("Acoustic TopOpt needs the inhomogeneous Helmholtz equation (activate complexFluid)");
        // We wrap not the material coefficients directly but already their reciprocal!
        CoefFunctionOpt* tmpFnc1 = new CoefFunctionOpt(domain->GetDesign(), coeffM, ACOU_BULK_MODULUS, this); // 1/K
        coeffM.reset(tmpFnc1);
        CoefFunctionOpt* tmpFnc2 = new CoefFunctionOpt(domain->GetDesign(), coeffK, DENSITY, this); // 1/rho
        coeffK.reset(tmpFnc2);
      }

      BaseBDBInt *stiffInt = nullptr;
      BaseBDBInt *massInt = nullptr;

      // ====================================================================
      // PML integrators
      // ====================================================================
      if (dampingList_[actRegion] == PML) {
        // outsource integrator definitions for PML case
        if (dim_ == 2)
          DefinePMLIntegrators<2>(actRegion, actSDList, curRegNode, c0, coeffK, coeffM, tempId, stiffInt, massInt);
        else
          DefinePMLIntegrators<3>(actRegion, actSDList, curRegNode, c0, coeffK, coeffM, tempId, stiffInt, massInt);
      }
      else {
        // define complex or default integrators
        if (dim_ == 2) {
          if (complexFluidFormulation_) {
            stiffInt = new BBInt<Complex>(new GradientOperator<FeH1, 2>(), coeffK, 1.0, updatedGeo_);
            massInt = new BBInt<Complex>(new IdentityOperator<FeH1, 2, 1, Complex>(), coeffM, 1.0, updatedGeo_);
          }
          else {
            stiffInt = new BBInt<Double>(new GradientOperator<FeH1, 2>(), coeffK, 1.0, updatedGeo_);
            massInt = new BBInt<Double>(new IdentityOperator<FeH1, 2, 1, Double>(), coeffM, 1.0, updatedGeo_);
          }
        }
        else if (dim_ == 3) {
          if (complexFluidFormulation_) {
            stiffInt = new BBInt<Complex>(new GradientOperator<FeH1, 3>(), coeffK, 1.0, updatedGeo_);
            massInt = new BBInt<Complex>(new IdentityOperator<FeH1, 3, 1, Complex>(), coeffM, 1.0, updatedGeo_);
          }
          else {
            stiffInt = new BBInt<Double>(new GradientOperator<FeH1, 3>(), coeffK, 1.0, updatedGeo_);
            massInt = new BBInt<Double>(new IdentityOperator<FeH1, 3, 1, Double>(), coeffM, 1.0, updatedGeo_);
          }
        }
        else {
          EXCEPTION("Unsupported dimension for PDE");
        }
      }

      // the integrator has a coef function but for the optimization case the opt coef needs to know also the integrator
      if(domain->HasDesign()){
        dynamic_pointer_cast<CoefFunctionOpt>(coeffM)->SetForm(massInt);
        dynamic_pointer_cast<CoefFunctionOpt>(coeffK)->SetForm(stiffInt);
      }

      // set stiffness integrator context and mass integrator context
      SetIntegratorContext(stiffInt, massInt, actRegion, actSDList, coeffK, coeffM);

      // ====================================================================
      // flow integrators
      // ====================================================================
      if (dim_ == 2) {
        if (isComplex_) {
          DefineConvectiveIntegrators<2, true>(actRegion, curRegNode, actSDList, coeffM);
        }
        else {
          DefineConvectiveIntegrators<2, false>(actRegion, curRegNode, actSDList, coeffM);
        }
      }

      else { /* if (dim_ == 3) */
        if (isComplex_) {
          DefineConvectiveIntegrators<3, true>(actRegion, curRegNode, actSDList, coeffM);
        }
        else {
          DefineConvectiveIntegrators<3, false>(actRegion, curRegNode, actSDList, coeffM);
        }
      }

      // ====================================================================
      // TDEF formulation
      // ====================================================================
      if (timeDomainEqFluidFormulation_ && formulation_ == ACOU_PRESSURE &&
          curRegNode->Get("timeDomainEqFluid")->As<std::string>() == "yes") {
        // only available for the acou pressure formulation
        // additional terms in the wave equation
        ReadTDEFCoefficients(iRegion, TDEFcoeffs_);
        DefineTDEFIntegrators(iRegion, polyId, integId, actSDList, TDEFcoeffs_);
      }
    }
  }

  template <UInt DIM>
  void AcousticPDE::DefinePMLIntegrators(RegionIdType actRegion, shared_ptr<ElemList> &actSDList, PtrParamNode &curRegNode,
                                         PtrCoefFct &c0, PtrCoefFct &coeffK, PtrCoefFct &coeffM, std::string &tempId,
                                         BaseBDBInt *&stiffInt, BaseBDBInt *&massInt) {
    std::string pmlDampId = ""; // dampId of a PML region
    std::string pmlFormul = ""; // formulation of the PML region ("classic", "shifted" or "curvilinear")
    PtrCoefFct coeffPMLDeterminant, coeffPMLTensor, coeffPMLVector, coeffPMLStiff, coeffPMLMass;

    // get the dampingId of the current PML region
    curRegNode->GetValue("dampingId", pmlDampId);

    // check for PML formulation of current region: classic(=Cartesian) or curvilinear
    PtrParamNode pmlNode = myParam_->Get("dampingList")->GetByVal("pml", "id", pmlDampId.c_str());
    pmlFormul = pmlNode->Get("formulation")->As<std::string>();

    // Check if the analysis type is harmonic or inverse-source
    if (analysistype_ == HARMONIC || analysistype_ == INVERSESOURCE) {
      // For complex fluids, compute a real-valued speed of sound for PML definition
      shared_ptr<CoefFunction> densR, blkR, c0R;
      if (complexFluidFormulation_) {
        densR = materials_[actRegion]->GetScalCoefFnc(DENSITY, Global::REAL);
        blkR = materials_[actRegion]->GetScalCoefFnc(ACOU_BULK_MODULUS, Global::REAL);
        c0R = CoefFunction::Generate(mp_, Global::REAL, CoefXprUnaryOp(mp_, CoefXprBinOp(mp_, blkR, densR, CoefXpr::OP_DIV), CoefXpr::OP_SQRT));
      } else {
        c0R = c0;
      }

      if (pmlFormul == "classic") {
        // here we only have a diagonal Jacobi matrix, so we store values as vector
        coeffPMLVector.reset(new CoefFunctionPML<Complex>(pmlNode, c0R, actSDList, regions_, true));
        coeffPMLDeterminant.reset(new CoefFunctionPML<Complex>(pmlNode, c0R, actSDList, regions_, false));
        // store pml factor for the postprocessing result
        matCoefs_[PML_DAMP_FACTOR]->AddRegion(actRegion, coeffPMLVector);
        // compute factors for the integrators
        coeffPMLStiff = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, coeffPMLDeterminant, coeffK, CoefXpr::OP_MULT));
        coeffPMLMass = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, coeffPMLDeterminant, coeffM, CoefXpr::OP_MULT));

        // PML integrators (stiff + mass)
        stiffInt = new BBInt<Complex>(new ScaledGradientOperator<FeH1, DIM, Complex>(), coeffPMLStiff, 1.0, updatedGeo_);
        massInt = new BBInt<Complex>(new IdentityOperator<FeH1, DIM, 1, Complex>(), coeffPMLMass, 1.0, updatedGeo_);

        // set coordinate stretching coefFunction
        stiffInt->SetBCoefFunctionOpB(coeffPMLVector);
      } else if (pmlFormul == "curvilinear") {
        if (complexFluidFormulation_)
          WARN("Curvilinear PML with complex fluid is not tested. Check results and consider adding a testcase.");

        // pointer to object that handles the computation of the curvilinear PML damping tensor
        coeffPMLTensor.reset(new CoefFunctionCurvilinearPML<Complex>(pmlNode, c0R, actSDList, regions_, OutputType::TENSOR));
        coeffPMLDeterminant.reset(new CoefFunctionCurvilinearPML<Complex>(pmlNode, c0R, actSDList, regions_, OutputType::DETERMINANT));

        // create some more CoefFunctionCurvilinearPMLs to store info about the PML parameters
        PtrCoefFct coeffPMLDampFactor, coeffPMLDistance;
        coeffPMLDampFactor.reset(new CoefFunctionCurvilinearPML<Complex>(pmlNode, c0R, actSDList, regions_, OutputType::DAMP_FACTOR));
        coeffPMLDistance.reset(new CoefFunctionCurvilinearPML<Complex>(pmlNode, c0R, actSDList, regions_, OutputType::DISTANCE));

        // assign the coefFunctions to the matCoefs_ to make them available in the DefinePostProcResults()
        matCoefs_[PML_TENSOR]->AddRegion(actRegion, coeffPMLTensor);
        matCoefs_[PML_DETERMINANT]->AddRegion(actRegion, coeffPMLDeterminant);
        matCoefs_[PML_DAMP_FACTOR]->AddRegion(actRegion, coeffPMLDampFactor);
        matCoefs_[PML_DISTANCE]->AddRegion(actRegion, coeffPMLDistance);

        // the Jakobi determinant is used as scaling coefficient for the BBIntegrators
        coeffPMLStiff = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, coeffPMLDeterminant, coeffK, CoefXpr::OP_MULT));
        coeffPMLMass = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, coeffPMLDeterminant, coeffM, CoefXpr::OP_MULT));

        // define integrators for curvilinear PML
        stiffInt = new BBInt<Complex>(new TensorScaledGradientOperator<FeH1, DIM, Complex>(), coeffPMLStiff, 1.0, updatedGeo_);
        massInt = new BBInt<Complex>(new IdentityOperator<FeH1, DIM, 1, Complex>(), coeffPMLMass, 1.0, updatedGeo_);

        // set coordinate stretching coefFunction
        stiffInt->SetBCoefFunctionOpB(coeffPMLTensor);
      } else {
        EXCEPTION("Unknown PML formulation '" << pmlFormul << "' for AcousticPDE. Possible formulations: 'classic', 'curvilinear'.")
      }
    } else {
      // Define transient integrators for non-harmonic cases
      if (pmlFormul == "classic") {
        DefineTransientPMLInts<DIM>(actSDList, pmlDampId, actRegion, tempId);
        stiffInt = new BBInt<Double>(new GradientOperator<FeH1, DIM>(), coeffK, 1.0, updatedGeo_);
        massInt = new BBInt<Double>(new IdentityOperator<FeH1, DIM, 1, Double>(), coeffM, 1.0, updatedGeo_);
      } else {
        EXCEPTION("Transient PML is currently only implemented in 'classic' formulation.")
      }
    }
  }

  template <UInt DIM>
  void AcousticPDE::DefineTransientPMLInts(shared_ptr<ElemList> eList, std::string id, RegionIdType actRegion, std::string tempId)
  {
    // define some material coeffunction as above...
    PtrCoefFct dens = materials_[eList->GetRegion()]->GetScalCoefFnc(DENSITY, Global::REAL);
    PtrCoefFct blk = materials_[eList->GetRegion()]->GetScalCoefFnc(ACOU_BULK_MODULUS, Global::REAL);
    PtrCoefFct one = CoefFunction::Generate(mp_, Global::REAL, "1.0");

    PtrCoefFct mechAcouFactor;
    CalcMechAcouFac(mechAcouFactor,dens);

    if (timeDomainEqFluidFormulation_) {
      EXCEPTION("PML for TDEF formulation not implemented yet.")
      // Note: in the case of the TDEF formulation, multiplication by the inverse high-frequency limit of the density is required
    }

    PtrCoefFct regionTemp;
    std::set<UInt> definedDofs;
    if (tempId != "") {
      // just real valued temperature allowed
      shared_ptr<ResultInfo> tempInfo = GetResultInfo(HEAT_MEAN_TEMPERATURE);
      // Add the region information
      PtrParamNode tempNode = myParam_->Get("temperatureList")->GetByVal("temperature", "name", tempId.c_str());
      ReadUserFieldValues(eList, tempNode, tempInfo->dofNames, tempInfo->entryType, false, regionTemp, definedDofs, updatedGeo_);
    }

    // compute speed of sound
    PtrCoefFct c0;
    ComputeSOS(c0, dens, blk, regionTemp, tempId);

    // coeffc needed in PML
    shared_ptr<CoefFunction> coeffc;
    if (sosAtLaplace_) {
      // c0^2
      ComputeSOS_SQR(coeffc, dens, blk, regionTemp, tempId);
    }
    else {
      // 1/c0^2
      coeffc = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, dens, blk, CoefXpr::OP_DIV));
    }

    // check the PML node to retrieve information
    PtrParamNode pmlNode = myParam_->Get("dampingList")->GetByVal("pml", "id", id.c_str());
    // get the diagonal entries of the Jacobi matrix as vector
    shared_ptr<CoefFunction> coeffPMLVector;
    coeffPMLVector.reset(new CoefFunctionPML<Double>(pmlNode, c0, eList, regions_, true));

    // store pml factor
    matCoefs_[PML_DAMP_FACTOR]->AddRegion(eList->GetRegion(), coeffPMLVector);

    shared_ptr<CoefFunctionCompound<Double>> coefA(new CoefFunctionCompound<Double>(mp_));
    shared_ptr<CoefFunctionCompound<Double>> coefB(new CoefFunctionCompound<Double>(mp_));
    shared_ptr<CoefFunctionCompound<Double>> coefAlpha(new CoefFunctionCompound<Double>(mp_));
    shared_ptr<CoefFunctionCompound<Double>> coefBeta(new CoefFunctionCompound<Double>(mp_));

    // --- Set the FE ansatz for the current region ---
    PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region", "name", eList->GetName().c_str());
    std::string polyId = curRegNode->Get("polyId")->As<std::string>();
    std::string integId = curRegNode->Get("integId")->As<std::string>();

    // ===> DEFINE PML DAMPINGFUNCTIONS
    shared_ptr<FeSpace> vecSpace = feFunctions_[ACOU_PMLAUXVEC]->GetFeSpace();
    vecSpace->SetRegionApproximation(eList->GetRegion(), polyId, integId);

    // ===> DEFINE PML DAMPINGFUNCTIONS
    std::map<std::string, PtrCoefFct> vars;
    std::map<std::string, PtrCoefFct> var;
    vars["a"] = coeffPMLVector;
    vars["b"] = coeffc;
    var["a"] = coeffPMLVector;

    StdVector<std::string> matAReal;
    StdVector<std::string> matBReal;
    std::string alpha;
    std::string beta;

    if (DIM == 3) {
      if (sosAtLaplace_) {
        const std::string Amat[] = {"a_0_R", "0.0", "0.0", "0.0", "a_1_R", "0.0", "0.0", "0.0", "a_2_R"};
        const std::string Bmat[] = {"( (a_0_R - a_1_R - a_2_R) * b_R )", "0.0", "0.0", "0.0", "( (a_1_R - a_0_R - a_2_R) * b_R )", "0.0", "0.0", "0.0", "( (a_2_R - a_1_R - a_0_R) * b_R )"};
        alpha = "(a_0_R + a_1_R + a_2_R )";
        beta = "( ( a_0_R * a_1_R) + (a_0_R * a_2_R) + (a_1_R * a_2_R) )";
        matAReal.Import(Amat, 9);
        matBReal.Import(Bmat, 9);
        coefA->SetTensor(matAReal, DIM, DIM, var);
        coefB->SetTensor(matBReal, DIM, DIM, vars);
        coefAlpha->SetScalar(alpha, var);
        coefBeta->SetScalar(beta, var);
      }
      else { // dim = 2
        const std::string Amat[] = {"a_0_R", "0.0", "0.0", "0.0", "a_1_R", "0.0", "0.0", "0.0", "a_2_R"};
        const std::string Bmat[] = {"( a_0_R - a_1_R - a_2_R ) ", "0.0", "0.0", "0.0", "( a_1_R - a_0_R - a_2_R )", "0.0", "0.0", "0.0", "( a_2_R - a_1_R - a_0_R) "};
        alpha = "( ( a_0_R + a_1_R + a_2_R) * b_R )";
        beta = " ( ( (a_0_R * a_1_R) + (a_0_R * a_2_R) + (a_1_R * a_2_R) ) * b_R )";
        matAReal.Import(Amat, 9);
        matBReal.Import(Bmat, 9);
        coefA->SetTensor(matAReal, DIM, DIM, var);
        coefB->SetTensor(matBReal, DIM, DIM, var);
        coefAlpha->SetScalar(alpha, vars);
        coefBeta->SetScalar(beta, vars);
      }
    }
    else {
      if (sosAtLaplace_) {
        std::string Amat[] = {"a_0_R", "0.0", "0.0", "a_1_R"};
        std::string Bmat[] = {"( (a_0_R - a_1_R) * b_R )", "0.0", "0.0", "( (a_1_R - a_0_R) * b_R )"};
        alpha = "( a_0_R + a_1_R ) ";
        beta = " ( a_0_R * a_1_R ) ";
        matAReal.Import(Amat, 4);
        matBReal.Import(Bmat, 4);
        coefA->SetTensor(matAReal, DIM, DIM, var);
        coefB->SetTensor(matBReal, DIM, DIM, vars);
        coefAlpha->SetScalar(alpha, var);
        coefBeta->SetScalar(beta, var);
      }
      else {
        const std::string Amat[] = {"a_0_R", "0.0", "0.0", "a_1_R"};
        const std::string Bmat[] = {"( a_0_R - a_1_R)", "0.0", "0.0", "( a_1_R - a_0_R )"};
        alpha = "( (a_0_R + a_1_R ) * b_R )";
        beta = " ( (a_0_R * a_1_R ) * b_R )";
        matAReal.Import(Amat, 4);
        matBReal.Import(Bmat, 4);
        coefA->SetTensor(matAReal, DIM, DIM, var);
        coefB->SetTensor(matBReal, DIM, DIM, var);
        coefAlpha->SetScalar(alpha, vars);
        coefBeta->SetScalar(beta, vars);
      }
    }

    // now lets define some integrators
    // the vectorial auxiliary variable is called U...
    // in case of mechacoucoupling, we need to multiply all integrators which
    // couple into acouPDE with -density or density if eigenfrequency
    PtrCoefFct mAcouCorrect_CoefAlpha = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, mechAcouFactor, coefAlpha, CoefXpr::OP_MULT));
    PtrCoefFct mAcouCorrect_CoefBeta = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, mechAcouFactor, coefBeta, CoefXpr::OP_MULT));

    BaseBDBInt *dampdPdt = new BBInt<>(new IdentityOperator<FeH1, DIM>(), mAcouCorrect_CoefAlpha, 1.0, updatedGeo_);
    BaseBDBInt *dampP = new BBInt<>(new IdentityOperator<FeH1, DIM>(), mAcouCorrect_CoefBeta, 1.0, updatedGeo_);
    // this is already integrated by parts...
    BaseBDBInt *divU = new ABInt<>(new GradientOperator<FeH1, DIM>(), new IdentityOperator<FeH1, DIM, DIM>(), mechAcouFactor, 1.0, updatedGeo_);

    BaseBDBInt *dUdt = new BBInt<>(new IdentityOperator<FeH1, DIM, DIM>(), one, 1.0, updatedGeo_);
    BaseBDBInt *AU = new BDBInt<>(new IdentityOperator<FeH1, DIM, DIM>(), coefA, 1.0, updatedGeo_);
    BaseBDBInt *BgradP = new ADBInt<>(new IdentityOperator<FeH1, DIM, DIM>(), new GradientOperator<FeH1, DIM>(), coefB, 1.0, updatedGeo_);

    dampdPdt->SetName("acouPML_dampdPdt");
    dampP->SetName("acouPML_dampP");
    divU->SetName("acouPML_divU");
    dUdt->SetName("acouPML_dUdt");
    AU->SetName("acouPML_AU");
    BgradP->SetName("acouPML_BgradP");

    BiLinFormContext *Context_dampdPdt = new BiLinFormContext(dampdPdt, DAMPING);
    BiLinFormContext *Context_dampP = new BiLinFormContext(dampP, STIFFNESS);
    BiLinFormContext *Context_divU = new BiLinFormContext(divU, STIFFNESS);
    BiLinFormContext *Context_dUdt = new BiLinFormContext(dUdt, DAMPING);
    BiLinFormContext *Context_AU = new BiLinFormContext(AU, STIFFNESS);
    BiLinFormContext *Context_BgradP = new BiLinFormContext(BgradP, STIFFNESS);

    Context_dampdPdt->SetEntities(eList, eList);
    Context_dampdPdt->SetFeFunctions(feFunctions_[formulation_], feFunctions_[formulation_]);
    Context_dampP->SetEntities(eList, eList);
    Context_dampP->SetFeFunctions(feFunctions_[formulation_], feFunctions_[formulation_]);
    Context_divU->SetEntities(eList, eList);
    Context_divU->SetFeFunctions(feFunctions_[formulation_], feFunctions_[ACOU_PMLAUXVEC]);
    Context_dUdt->SetEntities(eList, eList);
    Context_dUdt->SetFeFunctions(feFunctions_[ACOU_PMLAUXVEC], feFunctions_[ACOU_PMLAUXVEC]);
    Context_AU->SetEntities(eList, eList);
    Context_AU->SetFeFunctions(feFunctions_[ACOU_PMLAUXVEC], feFunctions_[ACOU_PMLAUXVEC]);
    Context_BgradP->SetEntities(eList, eList);
    Context_BgradP->SetFeFunctions(feFunctions_[ACOU_PMLAUXVEC], feFunctions_[formulation_]);

    assemble_->AddBiLinearForm(Context_dampdPdt);
    assemble_->AddBiLinearForm(Context_dampP);
    assemble_->AddBiLinearForm(Context_divU);
    assemble_->AddBiLinearForm(Context_dUdt);
    assemble_->AddBiLinearForm(Context_AU);
    assemble_->AddBiLinearForm(Context_BgradP);

    feFunctions_[ACOU_PMLAUXVEC]->AddEntityList(eList);

    if (!this->isAPML_ && DIM == 3) {
      shared_ptr<FeSpace> scalSpace = feFunctions_[ACOU_PMLAUXSCALAR]->GetFeSpace();
      scalSpace->SetRegionApproximation(eList->GetRegion(), polyId, integId);
      feFunctions_[ACOU_PMLAUXSCALAR]->AddEntityList(eList);

      shared_ptr<CoefFunctionCompound<Double>> coefGamma(new CoefFunctionCompound<Double>(mp_));
      shared_ptr<CoefFunctionCompound<Double>> coefC(new CoefFunctionCompound<Double>(mp_));

      StdVector<std::string> matCReal;

      if (sosAtLaplace_) {
        std::string gamma = "( a_0_R * a_1_R * a_2_R ) ";
        const std::string Cmat[] = {" (a_1_R * a_2_R * b_R) ", "0.0", "0.0", "0.0", " (a_0_R * a_2_R * b_R) ", "0.0", "0.0", "0.0", " (a_0_R * a_1_R * b_R) "};
        matCReal.Import(Cmat, 9);
        coefC->SetTensor(matCReal, 3, 3, vars);
        coefGamma->SetScalar(gamma, var);
      }
      else {
        std::string gamma = "( (a_0_R * a_1_R * a_2_R) * b_R )";
        const std::string Cmat[] = {" (a_1_R * a_2_R) ", "0.0", "0.0", "0.0", " (a_0_R * a_2_R) ", "0.0", "0.0", "0.0", " (a_0_R * a_1_R) "};
        matCReal.Import(Cmat, 9);
        coefC->SetTensor(matCReal, 3, 3, var);
        coefGamma->SetScalar(gamma, vars);
      }

      PtrCoefFct mAcouCorrect_CoefGamma = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, mechAcouFactor, coefGamma, CoefXpr::OP_MULT));

      BaseBDBInt *dampNu = new BBInt<>(new IdentityOperator<FeH1, DIM>(), mAcouCorrect_CoefGamma, 1.0, updatedGeo_);
      BaseBDBInt *CgradNu = new ADBInt<>(new IdentityOperator<FeH1, DIM, DIM>(), new GradientOperator<FeH1, DIM>(), coefC, -1.0, updatedGeo_);
      BaseBDBInt *dNudt = new BBInt<>(new IdentityOperator<FeH1, DIM>(), one, 1.0, updatedGeo_);
      BaseBDBInt *P = new BBInt<>(new IdentityOperator<FeH1, DIM>(), one, -1.0, updatedGeo_);

      // the scalar auxiliary variable is called Nu...
      dampNu->SetName("acouPML_dampNu");
      CgradNu->SetName("acouPML_CgIdentityOperaradNu");
      dNudt->SetName("acouPML_dNudt");
      P->SetName("acouPML_P");

      BiLinFormContext *Context_dampNu = new BiLinFormContext(dampNu, STIFFNESS);
      BiLinFormContext *Context_CgradNu = new BiLinFormContext(CgradNu, STIFFNESS);
      BiLinFormContext *Context_dNudt = new BiLinFormContext(dNudt, DAMPING);
      BiLinFormContext *Context_P = new BiLinFormContext(P, STIFFNESS);

      Context_dampNu->SetEntities(eList, eList);
      Context_dampNu->SetFeFunctions(feFunctions_[formulation_], feFunctions_[ACOU_PMLAUXSCALAR]);
      Context_CgradNu->SetEntities(eList, eList);
      Context_CgradNu->SetFeFunctions(feFunctions_[ACOU_PMLAUXVEC], feFunctions_[ACOU_PMLAUXSCALAR]);
      Context_dNudt->SetEntities(eList, eList);
      Context_dNudt->SetFeFunctions(feFunctions_[ACOU_PMLAUXSCALAR], feFunctions_[ACOU_PMLAUXSCALAR]);
      Context_P->SetEntities(eList, eList);
      Context_P->SetFeFunctions(feFunctions_[ACOU_PMLAUXSCALAR], feFunctions_[formulation_]);

      assemble_->AddBiLinearForm(Context_dampNu);
      assemble_->AddBiLinearForm(Context_CgradNu);
      assemble_->AddBiLinearForm(Context_dNudt);
      assemble_->AddBiLinearForm(Context_P);
    }
  }

  void AcousticPDE::DefineTDEFIntegrators(UInt iRegion, string polyId, string integId, shared_ptr<ElemList> actSDList, rationalCoeffsTDEF &TDEFcoeffs)
  {
    PtrCoefFct constOne = CoefFunction::Generate(mp_, Global::REAL, "1.0");
    RegionIdType actRegion = regions_[iRegion];

    // ====================================================================
    // Loop over real poles of inverse bulk modulus (compressibility)
    for (UInt iPole = 0; iPole < TDEFpoleNumber_.realC[iRegion]; iPole++) {

      // ====================================================================
      // K_PPHI1 (TDEF): stiffness integrator, TDEF (A_j^C term)
      // \int_{Omega_1} A_j^C p^\prime \phi_j^C d\Omega
      // ====================================================================
      shared_ptr<FeSpace> spacePhiC = feFunctions_[(SolutionType)(ACOU_TDEF_PHI_C_1 + iPole)]->GetFeSpace();
      spacePhiC->SetRegionApproximation(actRegion, polyId, integId);

      BiLinearForm *stiffIntTDEFPPHI1 = nullptr;
      if (dim_ == 2) {
        stiffIntTDEFPPHI1 = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), TDEFcoeffs.AC[iPole], 1.0, updatedGeo_);
      }
      else {
        stiffIntTDEFPPHI1 = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), TDEFcoeffs.AC[iPole], 1.0, updatedGeo_);
      }

      stiffIntTDEFPPHI1->SetName("AcousticStiffIntTDEFPPHI1_" + std::to_string(iPole));
      BiLinFormContext *stiffIntTDEFPPHI1Context = nullptr;
      stiffIntTDEFPPHI1Context = new BiLinFormContext(stiffIntTDEFPPHI1, STIFFNESS);
      stiffIntTDEFPPHI1Context->SetEntities(actSDList, actSDList);
      stiffIntTDEFPPHI1Context->SetFeFunctions(feFunctions_[formulation_], feFunctions_[(SolutionType)(ACOU_TDEF_PHI_C_1 + iPole)]);
      assemble_->AddBiLinearForm(stiffIntTDEFPPHI1Context);
      feFunctions_[(SolutionType)(ACOU_TDEF_PHI_C_1 + iPole)]->AddEntityList(actSDList);
    } // end loop stiffIntTDEFPPHI1

    // ====================================================================
    // Loop over complex-conjugated poles of inverse bulk modulus (compressibility)
    for (UInt iPole = 0; iPole < TDEFpoleNumber_.complC[iRegion]; iPole++) {

      // ====================================================================
      // D_PPSI1 (TDEF): damping integrator, TDEF (\frac{B_k^C,\gamma_k^C} term)
      // \int_{Omega_1} \frac{B_k^C,\gamma_k^C} p^\prime \dot{\psi}_k^C d\Omega
      // ====================================================================
      shared_ptr<FeSpace> spacePsiC = feFunctions_[(SolutionType)(ACOU_TDEF_PSI_C_1 + iPole)]->GetFeSpace();
      spacePsiC->SetRegionApproximation(actRegion, polyId, integId);

      BiLinearForm *dampIntTDEFPPSI1 = nullptr;
      PtrCoefFct fncBCgammaC;
      fncBCgammaC = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, TDEFcoeffs.BC[iPole], TDEFcoeffs.GammaC[iPole], CoefXpr::OP_DIV));
      if (dim_ == 2) {
        dampIntTDEFPPSI1 = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), fncBCgammaC, 1.0, updatedGeo_);
      }
      else {
        dampIntTDEFPPSI1 = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), fncBCgammaC, 1.0, updatedGeo_);
      }

      dampIntTDEFPPSI1->SetName("AcousticDampIntTDEFPPSI1_" + std::to_string(iPole));
      BiLinFormContext *dampIntTDEFPPSI1Context = nullptr;
      dampIntTDEFPPSI1Context = new BiLinFormContext(dampIntTDEFPPSI1, DAMPING);
      dampIntTDEFPPSI1Context->SetEntities(actSDList, actSDList);
      dampIntTDEFPPSI1Context->SetFeFunctions(feFunctions_[formulation_], feFunctions_[(SolutionType)(ACOU_TDEF_PSI_C_1 + iPole)]);
      assemble_->AddBiLinearForm(dampIntTDEFPPSI1Context);
      feFunctions_[(SolutionType)(ACOU_TDEF_PSI_C_1 + iPole)]->AddEntityList(actSDList);

      // ====================================================================
      // K_PPSI1 (TDEF): stiffness integrator, TDEF (D_k^C term)
      // \int_{Omega_1} D_k^C p^\prime \psi_k^C d\Omega
      // ====================================================================
      BiLinearForm *stiffIntTDEFPPSI1 = nullptr;
      PtrCoefFct fncQuotC;
      fncQuotC = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, TDEFcoeffs.BetaC[iPole], TDEFcoeffs.GammaC[iPole], CoefXpr::OP_DIV));
      PtrCoefFct fncTermC;
      fncTermC = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, TDEFcoeffs.BC[iPole], fncQuotC, CoefXpr::OP_MULT));
      PtrCoefFct fncDC;
      fncDC = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, TDEFcoeffs.CC[iPole], fncTermC, CoefXpr::OP_ADD));

      if (dim_ == 2) {
        stiffIntTDEFPPSI1 = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), fncDC, 1.0, updatedGeo_);
      }
      else {
        stiffIntTDEFPPSI1 = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), fncDC, 1.0, updatedGeo_);
      }

      stiffIntTDEFPPSI1->SetName("AcousticStiffIntTDEFPPSI1_" + std::to_string(iPole));
      BiLinFormContext *stiffIntTDEFPPSI1Context = nullptr;
      stiffIntTDEFPPSI1Context = new BiLinFormContext(stiffIntTDEFPPSI1, STIFFNESS);
      stiffIntTDEFPPSI1Context->SetEntities(actSDList, actSDList);
      stiffIntTDEFPPSI1Context->SetFeFunctions(feFunctions_[formulation_], feFunctions_[(SolutionType)(ACOU_TDEF_PSI_C_1 + iPole)]);
      assemble_->AddBiLinearForm(stiffIntTDEFPPSI1Context);
    } // end loop stiffIntTDEFPPSI1 and dampIntTDEFPPSI1

    // ====================================================================
    // Loop over real poles of inverse density (specific volume)
    for (UInt iPole = 0; iPole < TDEFpoleNumber_.realV[iRegion]; iPole++) {

      // ====================================================================
      // K_PPHI2 (TDEF): stiffness integrator, TDEF (A_l^V term)
      // \int_{Omega_1} A_l^V \nabla p^\prime \nabla \phi_l^V d\Omega
      // ====================================================================
      shared_ptr<FeSpace> spacePhiV = feFunctions_[(SolutionType)(ACOU_TDEF_PHI_V_1 + iPole)]->GetFeSpace();
      spacePhiV->SetRegionApproximation(actRegion, polyId, integId);
      BiLinearForm *stiffIntTDEFPPHI2 = nullptr;

      if (dim_ == 2) {
        stiffIntTDEFPPHI2 = new BBInt<Double>(new GradientOperator<FeH1, 2>(), TDEFcoeffs.AV[iPole], 1.0, updatedGeo_);
      }
      else {
        stiffIntTDEFPPHI2 = new BBInt<Double>(new GradientOperator<FeH1, 3>(), TDEFcoeffs.AV[iPole], 1.0, updatedGeo_);
      }

      stiffIntTDEFPPHI2->SetName("AcousticStiffIntTDEFPPHI2_" + std::to_string(iPole));
      BiLinFormContext *stiffIntTDEFPPHI2Context = nullptr;
      stiffIntTDEFPPHI2Context = new BiLinFormContext(stiffIntTDEFPPHI2, STIFFNESS);
      stiffIntTDEFPPHI2Context->SetEntities(actSDList, actSDList);
      stiffIntTDEFPPHI2Context->SetFeFunctions(feFunctions_[formulation_], feFunctions_[(SolutionType)(ACOU_TDEF_PHI_V_1 + iPole)]);
      assemble_->AddBiLinearForm(stiffIntTDEFPPHI2Context);
      feFunctions_[(SolutionType)(ACOU_TDEF_PHI_V_1 + iPole)]->AddEntityList(actSDList);
    } // end loop stiffIntTDEFPPHI2

    // ====================================================================
    // Loop over complex-conjugated poles of inverse density (specific volume)
    for (UInt iPole = 0; iPole < TDEFpoleNumber_.complV[iRegion]; iPole++) {

      // ====================================================================
      // D_PPSI2 (TDEF): damping integrator, TDEF (\frac{B_k^V,\gamma_k^V} term)
      // \int_{Omega_1} \frac{B_k^V,\gamma_k^V} p^\prime \dot{\psi}_m^C d\Omega
      // ====================================================================
      shared_ptr<FeSpace> spacePsiV = feFunctions_[(SolutionType)(ACOU_TDEF_PSI_V_1 + iPole)]->GetFeSpace();
      spacePsiV->SetRegionApproximation(actRegion, polyId, integId);
      BiLinearForm *dampIntTDEFPPSI2 = nullptr;
      PtrCoefFct fncBVgammaV;
      fncBVgammaV = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, TDEFcoeffs.BV[iPole], TDEFcoeffs.GammaV[iPole], CoefXpr::OP_DIV));
      if (dim_ == 2) {
        dampIntTDEFPPSI2 = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), fncBVgammaV, 1.0, updatedGeo_);
      }
      else {
        dampIntTDEFPPSI2 = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), fncBVgammaV, 1.0, updatedGeo_);
      }

      dampIntTDEFPPSI2->SetName("AcousticDampIntTDEFPPSI2_" + std::to_string(iPole));
      BiLinFormContext *dampIntTDEFPPSI2Context = nullptr;
      dampIntTDEFPPSI2Context = new BiLinFormContext(dampIntTDEFPPSI2, DAMPING);
      dampIntTDEFPPSI2Context->SetEntities(actSDList, actSDList);
      dampIntTDEFPPSI2Context->SetFeFunctions(feFunctions_[formulation_], feFunctions_[(SolutionType)(ACOU_TDEF_PSI_V_1 + iPole)]);
      assemble_->AddBiLinearForm(dampIntTDEFPPSI2Context);
      feFunctions_[(SolutionType)(ACOU_TDEF_PSI_V_1 + iPole)]->AddEntityList(actSDList);

      // ====================================================================
      // K_PPSI2 (TDEF): stiffness integrator, TDEF (D_k^V term)
      // \int_{Omega_1} D_m^V \nabla p^\prime \nabla \psi_m^V d\Omega
      // ====================================================================
      BiLinearForm *stiffIntTDEFPPSI2 = nullptr;
      PtrCoefFct fncQuotV;
      fncQuotV = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, TDEFcoeffs.BetaV[iPole], TDEFcoeffs.GammaV[iPole], CoefXpr::OP_DIV));
      PtrCoefFct fncTermV;
      fncTermV = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, TDEFcoeffs.BV[iPole], fncQuotV, CoefXpr::OP_MULT));
      PtrCoefFct fncDV;
      fncDV = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, TDEFcoeffs.CV[iPole], fncTermV, CoefXpr::OP_ADD));

      if (dim_ == 2) {
        stiffIntTDEFPPSI2 = new BBInt<Double>(new GradientOperator<FeH1, 2>(), fncDV, 1.0, updatedGeo_);
      }
      else {
        stiffIntTDEFPPSI2 = new BBInt<Double>(new GradientOperator<FeH1, 3>(), fncDV, 1.0, updatedGeo_);
      }

      stiffIntTDEFPPSI2->SetName("AcousticStiffIntTDEFPPSI2_" + std::to_string(iPole));
      BiLinFormContext *stiffIntTDEFPPSI2Context = nullptr;
      stiffIntTDEFPPSI2Context = new BiLinFormContext(stiffIntTDEFPPSI2, STIFFNESS);
      stiffIntTDEFPPSI2Context->SetEntities(actSDList, actSDList);
      stiffIntTDEFPPSI2Context->SetFeFunctions(feFunctions_[formulation_], feFunctions_[(SolutionType)(ACOU_TDEF_PSI_V_1 + iPole)]);
      assemble_->AddBiLinearForm(stiffIntTDEFPPSI2Context);
    } // end loop stiffIntTDEFPPSI2 and dampIntTDEFPPSI2

    // ######################################################
    // Auxiliar differential equations (ADE) section
    // ######################################################

    // ====================================================================
    // Loop over real poles of inverse compression modulus (compressibility)
    for (UInt iPole = 0; iPole < TDEFcoeffs.AlphaC.GetSize(); iPole++) {

      // ====================================================================
      // D_PHI1PHI1 (TDEF): damping integrator, TDEF
      // \int_{Omega_1} \phi_j^{C,\prime} \prime{\phi}_j^C d\Omega
      // ====================================================================
      BiLinearForm *dampIntTDEFPHI1PHI1 = nullptr;
      if (dim_ == 2) {
        dampIntTDEFPHI1PHI1 = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), constOne, 1.0, updatedGeo_);
      }
      else {
        dampIntTDEFPHI1PHI1 = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), constOne, 1.0, updatedGeo_);
      }

      dampIntTDEFPHI1PHI1->SetName("AcousticDampIntTDEFPHI1PHI1_" + std::to_string(iPole));
      BiLinFormContext *dampIntTDEFPHI1PHI1Context = nullptr;
      dampIntTDEFPHI1PHI1Context = new BiLinFormContext(dampIntTDEFPHI1PHI1, DAMPING);
      dampIntTDEFPHI1PHI1Context->SetEntities(actSDList, actSDList);
      dampIntTDEFPHI1PHI1Context->SetFeFunctions(feFunctions_[(SolutionType)(ACOU_TDEF_PHI_C_1 + iPole)], feFunctions_[(SolutionType)(ACOU_TDEF_PHI_C_1 + iPole)]);
      assemble_->AddBiLinearForm(dampIntTDEFPHI1PHI1Context);

      // ====================================================================
      // K_PHI1PHI1 (TDEF): stiffness integrator, TDEF (\alpha_j^C term)
      // \int_{Omega_1} \alpha_j^C \phi_j^{C,\prime} \phi_j^C d\Omega
      // ====================================================================
      BiLinearForm *stiffIntTDEFPHI1PHI1 = nullptr;
      if (dim_ == 2) {
        stiffIntTDEFPHI1PHI1 = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), TDEFcoeffs.AlphaC[iPole], 1.0, updatedGeo_);
      }
      else {
        stiffIntTDEFPHI1PHI1 = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), TDEFcoeffs.AlphaC[iPole], 1.0, updatedGeo_);
      }

      stiffIntTDEFPHI1PHI1->SetName("AcousticStiffIntTDEFPHI1PHI1_" + std::to_string(iPole));
      BiLinFormContext *stiffIntTDEFPHI1PHI1Context = nullptr;
      stiffIntTDEFPHI1PHI1Context = new BiLinFormContext(stiffIntTDEFPHI1PHI1, STIFFNESS);
      stiffIntTDEFPHI1PHI1Context->SetEntities(actSDList, actSDList);
      stiffIntTDEFPHI1PHI1Context->SetFeFunctions(feFunctions_[(SolutionType)(ACOU_TDEF_PHI_C_1 + iPole)], feFunctions_[(SolutionType)(ACOU_TDEF_PHI_C_1 + iPole)]);
      assemble_->AddBiLinearForm(stiffIntTDEFPHI1PHI1Context);

      // ====================================================================
      // M_PHI1P (TDEF): mass integrator, TDEF
      // -\int_{Omega_1} \phi_j^{C,\prime} \ddot{p}_a d\Omega
      // ====================================================================
      BiLinearForm *massIntTDEFPHI1P = nullptr;
      if (dim_ == 2) {
        massIntTDEFPHI1P = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), constOne, -1.0, updatedGeo_);
      }
      else {
        massIntTDEFPHI1P = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), constOne, -1.0, updatedGeo_);
      }

      massIntTDEFPHI1P->SetName("AcousticMassIntTDEFPHI1P_" + std::to_string(iPole));
      BiLinFormContext *massIntTDEFPHI1PContext = nullptr;
      massIntTDEFPHI1PContext = new BiLinFormContext(massIntTDEFPHI1P, MASS);
      massIntTDEFPHI1PContext->SetEntities(actSDList, actSDList);
      massIntTDEFPHI1PContext->SetFeFunctions(feFunctions_[(SolutionType)(ACOU_TDEF_PHI_C_1 + iPole)], feFunctions_[formulation_]);
      assemble_->AddBiLinearForm(massIntTDEFPHI1PContext);
    } // end loop dampIntTDEFPHI1PHI1, stiffIntTDEFPHI1PHI1 and massInteTDEFPHI1P

    // ====================================================================
    // Loop over complex poles of inverse compression modulus (compressibility)
    for (UInt iPole = 0; iPole < TDEFcoeffs.GammaC.GetSize(); iPole++) {

      // ====================================================================
      // M_PSI1PSI1 (TDEF): mass integrator, TDEF (\frac{1,\gamma_k^C} term)
      // \int_{Omega_1} \frac{1,\gamma_k^C} \psi_k^{C,\prime} \prime{\psi}_k^C d\Omega
      // ====================================================================
      BiLinearForm *massIntTDEFPSI1PSI1 = nullptr;
      PtrCoefFct coefOneOverGammaC;
      coefOneOverGammaC = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, constOne, TDEFcoeffs.GammaC[iPole], CoefXpr::OP_DIV));
      if (dim_ == 2) {
        massIntTDEFPSI1PSI1 = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), coefOneOverGammaC, 1.0, updatedGeo_);
      }
      else {
        massIntTDEFPSI1PSI1 = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), coefOneOverGammaC, 1.0, updatedGeo_);
      }

      massIntTDEFPSI1PSI1->SetName("AcousticMassIntTDEFPSI1PSI1_" + std::to_string(iPole));
      BiLinFormContext *massIntTDEFPSI1PSI1Context = nullptr;
      massIntTDEFPSI1PSI1Context = new BiLinFormContext(massIntTDEFPSI1PSI1, MASS);
      massIntTDEFPSI1PSI1Context->SetEntities(actSDList, actSDList);
      massIntTDEFPSI1PSI1Context->SetFeFunctions(feFunctions_[(SolutionType)(ACOU_TDEF_PSI_C_1 + iPole)], feFunctions_[(SolutionType)(ACOU_TDEF_PSI_C_1 + iPole)]);
      assemble_->AddBiLinearForm(massIntTDEFPSI1PSI1Context);

      // ====================================================================
      // D_PSI1PSI1 (TDEF): damping integrator, TDEF (\frac{2 \beta_k^C,\gamma_k^C} term)
      // \int_{Omega_1} \frac{2 \beta_k^C,\gamma_k^C} \psi_j^{C,\prime} \dot{\psi}_j^C d\Omega
      // ====================================================================
      BiLinearForm *dampIntTDEFPSI1PSI1 = nullptr;
      PtrCoefFct coefTwoBetaC;
      coefTwoBetaC = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, CoefFunction::Generate(mp_, Global::REAL, "2.0"), TDEFcoeffs.BetaC[iPole], CoefXpr::OP_MULT));
      PtrCoefFct coefTwoBetaCOverGammaC;
      coefTwoBetaCOverGammaC = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, coefTwoBetaC, TDEFcoeffs.GammaC[iPole], CoefXpr::OP_DIV));
      if (dim_ == 2) {
        dampIntTDEFPSI1PSI1 = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), coefTwoBetaCOverGammaC, 1.0, updatedGeo_);
      }
      else {
        dampIntTDEFPSI1PSI1 = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), coefTwoBetaCOverGammaC, 1.0, updatedGeo_);
      }

      dampIntTDEFPSI1PSI1->SetName("AcousticDampIntTDEFPSI1PSI1_" + std::to_string(iPole));
      BiLinFormContext *dampIntTDEFPSI1PSI1Context = nullptr;
      dampIntTDEFPSI1PSI1Context = new BiLinFormContext(dampIntTDEFPSI1PSI1, DAMPING);
      dampIntTDEFPSI1PSI1Context->SetEntities(actSDList, actSDList);
      dampIntTDEFPSI1PSI1Context->SetFeFunctions(feFunctions_[(SolutionType)(ACOU_TDEF_PSI_C_1 + iPole)], feFunctions_[(SolutionType)(ACOU_TDEF_PSI_C_1 + iPole)]);
      assemble_->AddBiLinearForm(dampIntTDEFPSI1PSI1Context);

      // ====================================================================
      // K_PSI1PSI1 (TDEF): stiffness integrator, TDEF (\delta_j^C term)
      // \int_{Omega_1} \delta_j^C \psi_j^{C,\prime} \psi_j^C d\Omega
      // ====================================================================
      BiLinearForm *stiffIntTDEFPSI1PSI1 = nullptr;
      PtrCoefFct fncQuotC;
      fncQuotC = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, TDEFcoeffs.BetaC[iPole], TDEFcoeffs.GammaC[iPole], CoefXpr::OP_DIV));
      PtrCoefFct fncTermC;
      fncTermC = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, TDEFcoeffs.BetaC[iPole], fncQuotC, CoefXpr::OP_MULT));
      PtrCoefFct fncDeltaC;
      fncDeltaC = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, TDEFcoeffs.GammaC[iPole], fncTermC, CoefXpr::OP_ADD));
      if (dim_ == 2) {
        stiffIntTDEFPSI1PSI1 = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), fncDeltaC, 1.0, updatedGeo_);
      }
      else {
        stiffIntTDEFPSI1PSI1 = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), fncDeltaC, 1.0, updatedGeo_);
      }

      stiffIntTDEFPSI1PSI1->SetName("AcousticStiffIntTDEFPSI1PSI1_" + std::to_string(iPole));
      BiLinFormContext *stiffIntTDEFPSI1PSI1Context = nullptr;
      stiffIntTDEFPSI1PSI1Context = new BiLinFormContext(stiffIntTDEFPSI1PSI1, STIFFNESS);
      stiffIntTDEFPSI1PSI1Context->SetEntities(actSDList, actSDList);
      stiffIntTDEFPSI1PSI1Context->SetFeFunctions(feFunctions_[(SolutionType)(ACOU_TDEF_PSI_C_1 + iPole)], feFunctions_[(SolutionType)(ACOU_TDEF_PSI_C_1 + iPole)]);
      assemble_->AddBiLinearForm(stiffIntTDEFPSI1PSI1Context);

      // ====================================================================
      // M_PSI1P (TDEF): mass integrator, TDEF
      // -\int_{Omega_1} \psi_j^{C,\prime} \ddot{p}_a d\Omega
      // ====================================================================
      BiLinearForm *massIntTDEFPSI1P = nullptr;
      if (dim_ == 2) {
        massIntTDEFPSI1P = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), constOne, -1.0, updatedGeo_);
      }
      else {
        massIntTDEFPSI1P = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), constOne, -1.0, updatedGeo_);
      }

      massIntTDEFPSI1P->SetName("AcousticMassIntTDEFPSI1P_" + std::to_string(iPole));
      BiLinFormContext *massIntTDEFPSI1PContext = nullptr;
      massIntTDEFPSI1PContext = new BiLinFormContext(massIntTDEFPSI1P, MASS);
      massIntTDEFPSI1PContext->SetEntities(actSDList, actSDList);
      massIntTDEFPSI1PContext->SetFeFunctions(feFunctions_[(SolutionType)(ACOU_TDEF_PSI_C_1 + iPole)], feFunctions_[formulation_]);
      assemble_->AddBiLinearForm(massIntTDEFPSI1PContext);
    } // end loop massIntTDEFPSI1PSI1, dampIntTDEFPSI1PSI1, stiffIntTDEFPSI1PSI1 and massInteTDEFPSI1P

    // ====================================================================
    // Loop over real poles of inverse density (specific volume)
    for (UInt iPole = 0; iPole < TDEFcoeffs.AlphaV.GetSize(); iPole++) {

      // ====================================================================
      // D_PHI2PHI2 (TDEF): damping integrator, TDEF
      // \int_{Omega_1} \phi_l^{V,\prime} \prime{\phi}_l^V d\Omega
      // ====================================================================
      BiLinearForm *dampIntTDEFPHI2PHI2 = nullptr;
      if (dim_ == 2) {
        dampIntTDEFPHI2PHI2 = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), constOne, 1.0, updatedGeo_);
      }
      else {
        dampIntTDEFPHI2PHI2 = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), constOne, 1.0, updatedGeo_);
      }

      dampIntTDEFPHI2PHI2->SetName("AcousticDampIntTDEFPHI2PHI2_" + std::to_string(iPole));
      BiLinFormContext *dampIntTDEFPHI2PHI2Context = nullptr;
      dampIntTDEFPHI2PHI2Context = new BiLinFormContext(dampIntTDEFPHI2PHI2, DAMPING);
      dampIntTDEFPHI2PHI2Context->SetEntities(actSDList, actSDList);
      dampIntTDEFPHI2PHI2Context->SetFeFunctions(feFunctions_[(SolutionType)(ACOU_TDEF_PHI_V_1 + iPole)], feFunctions_[(SolutionType)(ACOU_TDEF_PHI_V_1 + iPole)]);
      assemble_->AddBiLinearForm(dampIntTDEFPHI2PHI2Context);

      // ====================================================================
      // K_PHI2PHI2 (TDEF): stiffness integrator, TDEF (\alpha_j^C term)
      // \int_{Omega_1} \alpha_l^V \phi_l^{V,\prime} \phi_l^V d\Omega
      // ====================================================================
      BiLinearForm *stiffIntTDEFPHI2PHI2 = nullptr;
      if (dim_ == 2) {
        stiffIntTDEFPHI2PHI2 = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), TDEFcoeffs.AlphaV[iPole], 1.0, updatedGeo_);
      }
      else {
        stiffIntTDEFPHI2PHI2 = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), TDEFcoeffs.AlphaV[iPole], 1.0, updatedGeo_);
      }

      stiffIntTDEFPHI2PHI2->SetName("AcousticStiffIntTDEFPHI2PHI2_" + std::to_string(iPole));
      BiLinFormContext *stiffIntTDEFPHI2PHI2Context = nullptr;
      stiffIntTDEFPHI2PHI2Context = new BiLinFormContext(stiffIntTDEFPHI2PHI2, STIFFNESS);
      stiffIntTDEFPHI2PHI2Context->SetEntities(actSDList, actSDList);
      stiffIntTDEFPHI2PHI2Context->SetFeFunctions(feFunctions_[(SolutionType)(ACOU_TDEF_PHI_V_1 + iPole)], feFunctions_[(SolutionType)(ACOU_TDEF_PHI_V_1 + iPole)]);
      assemble_->AddBiLinearForm(stiffIntTDEFPHI2PHI2Context);

      // ====================================================================
      // K_PHI2P (TDEF): stiffness integrator, TDEF
      // -\int_{Omega_1} \phi_l^{V,\prime} \ddot{p}_a d\Omega
      // ====================================================================
      BiLinearForm *stiffIntTDEFPHI2P = nullptr;
      if (dim_ == 2) {
        stiffIntTDEFPHI2P = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), constOne, -1.0, updatedGeo_);
      }
      else {
        stiffIntTDEFPHI2P = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), constOne, -1.0, updatedGeo_);
      }

      stiffIntTDEFPHI2P->SetName("AcousticStiffIntTDEFPHI2P_" + std::to_string(iPole));
      BiLinFormContext *stiffIntTDEFPHI2PContext = nullptr;
      stiffIntTDEFPHI2PContext = new BiLinFormContext(stiffIntTDEFPHI2P, STIFFNESS);
      stiffIntTDEFPHI2PContext->SetEntities(actSDList, actSDList);
      stiffIntTDEFPHI2PContext->SetFeFunctions(feFunctions_[(SolutionType)(ACOU_TDEF_PHI_V_1 + iPole)], feFunctions_[formulation_]);
      assemble_->AddBiLinearForm(stiffIntTDEFPHI2PContext);
    } // end loop dampIntTDEFPHI2PHI2, stiffIntTDEFPHI2PHI2 and stiffIntTDEFPHI2P

    // ====================================================================
    // Loop over complex-conjugated poles of inverse density (specific volume)
    for (UInt iPole = 0; iPole < TDEFcoeffs.GammaV.GetSize(); iPole++) {

      // ====================================================================
      // M_PSI2PSI2 (TDEF): mass integrator, TDEF (\frac{1,\gamma_k^V} term)
      // \int_{Omega_1} \frac{1,\gamma_m^V} \psi_m^{V,\prime} \prime{\psi}_m^V d\Omega
      // ====================================================================
      BiLinearForm *massIntTDEFPSI2PSI2 = nullptr;
      PtrCoefFct coefOneOverGammaV;
      coefOneOverGammaV = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, constOne, TDEFcoeffs.GammaV[iPole], CoefXpr::OP_DIV));
      if (dim_ == 2) {
        massIntTDEFPSI2PSI2 = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), coefOneOverGammaV, 1.0, updatedGeo_);
      }
      else {
        massIntTDEFPSI2PSI2 = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), coefOneOverGammaV, 1.0, updatedGeo_);
      }

      massIntTDEFPSI2PSI2->SetName("AcousticMassIntTDEFPSI2PSI2_" + std::to_string(iPole));
      BiLinFormContext *massIntTDEFPSI2PSI2Context = nullptr;
      massIntTDEFPSI2PSI2Context = new BiLinFormContext(massIntTDEFPSI2PSI2, MASS);
      massIntTDEFPSI2PSI2Context->SetEntities(actSDList, actSDList);
      massIntTDEFPSI2PSI2Context->SetFeFunctions(feFunctions_[(SolutionType)(ACOU_TDEF_PSI_V_1 + iPole)], feFunctions_[(SolutionType)(ACOU_TDEF_PSI_V_1 + iPole)]);
      assemble_->AddBiLinearForm(massIntTDEFPSI2PSI2Context);

      // ====================================================================
      // D_PSI2PSI2 (TDEF): damping integrator, TDEF (\frac{2 \beta_m^V,\gamma_m^V} term)
      // \int_{Omega_1} \frac{2 \beta_m^V,\gamma_m^V} \psi_m^{V,\prime} \dot{\psi}_m^V d\Omega
      // ====================================================================
      BiLinearForm *dampIntTDEFPSI2PSI2 = nullptr;
      PtrCoefFct coefTwoBetaV;
      coefTwoBetaV = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, CoefFunction::Generate(mp_, Global::REAL, "2.0"), TDEFcoeffs.BetaV[iPole], CoefXpr::OP_MULT));
      PtrCoefFct coefTwoBetaVOverGammaV;
      coefTwoBetaVOverGammaV = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, coefTwoBetaV, TDEFcoeffs.GammaV[iPole], CoefXpr::OP_DIV));
      if (dim_ == 2) {
        dampIntTDEFPSI2PSI2 = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), coefTwoBetaVOverGammaV, 1.0, updatedGeo_);
      }
      else {
        dampIntTDEFPSI2PSI2 = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), coefTwoBetaVOverGammaV, 1.0, updatedGeo_);
      }

      dampIntTDEFPSI2PSI2->SetName("AcousticDampIntTDEFPSI2PSI2_" + std::to_string(iPole));
      BiLinFormContext *dampIntTDEFPSI2PSI2Context = nullptr;
      dampIntTDEFPSI2PSI2Context = new BiLinFormContext(dampIntTDEFPSI2PSI2, DAMPING);
      dampIntTDEFPSI2PSI2Context->SetEntities(actSDList, actSDList);
      dampIntTDEFPSI2PSI2Context->SetFeFunctions(feFunctions_[(SolutionType)(ACOU_TDEF_PSI_V_1 + iPole)], feFunctions_[(SolutionType)(ACOU_TDEF_PSI_V_1 + iPole)]);
      assemble_->AddBiLinearForm(dampIntTDEFPSI2PSI2Context);

      // ====================================================================
      // K_PSI2PSI2 (TDEF): stiffness integrator, TDEF (\delta_m^V term)
      // \int_{Omega_1} \delta_m^V \psi_m^{V,\prime} \psi_m^V d\Omega
      // ====================================================================
      BiLinearForm *stiffIntTDEFPSI2PSI2 = nullptr;
      PtrCoefFct fncQuotV;
      fncQuotV = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, TDEFcoeffs.BetaV[iPole], TDEFcoeffs.GammaV[iPole], CoefXpr::OP_DIV));
      PtrCoefFct fncTermV;
      fncTermV = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, TDEFcoeffs.BetaV[iPole], fncQuotV, CoefXpr::OP_MULT));
      PtrCoefFct fncDeltaV;
      fncDeltaV = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, TDEFcoeffs.GammaV[iPole], fncTermV, CoefXpr::OP_ADD));
      if (dim_ == 2) {
        stiffIntTDEFPSI2PSI2 = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), fncDeltaV, 1.0, updatedGeo_);
      }
      else {
        stiffIntTDEFPSI2PSI2 = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), fncDeltaV, 1.0, updatedGeo_);
      }

      stiffIntTDEFPSI2PSI2->SetName("AcousticStiffIntTDEFPSI2PSI12_" + std::to_string(iPole));
      BiLinFormContext *stiffIntTDEFPSI2PSI2Context = nullptr;
      stiffIntTDEFPSI2PSI2Context = new BiLinFormContext(stiffIntTDEFPSI2PSI2, STIFFNESS);
      stiffIntTDEFPSI2PSI2Context->SetEntities(actSDList, actSDList);
      stiffIntTDEFPSI2PSI2Context->SetFeFunctions(feFunctions_[(SolutionType)(ACOU_TDEF_PSI_V_1 + iPole)], feFunctions_[(SolutionType)(ACOU_TDEF_PSI_V_1 + iPole)]);
      assemble_->AddBiLinearForm(stiffIntTDEFPSI2PSI2Context);

      // ====================================================================
      // M_PSI2P (TDEF): mass integrator, TDEF
      // -\int_{Omega_1} \psi_m^{V,\prime} \ddot{p}_a d\Omega
      // ====================================================================
      BiLinearForm *massIntTDEFPSI2P = nullptr;
      if (dim_ == 2) {
        massIntTDEFPSI2P = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), constOne, -1.0, updatedGeo_);
      }
      else {
        massIntTDEFPSI2P = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), constOne, -1.0, updatedGeo_);
      }

      massIntTDEFPSI2P->SetName("AcousticMassIntTDEFPSI2P_" + std::to_string(iPole));
      BiLinFormContext *massIntTDEFPSI2PContext = nullptr;
      massIntTDEFPSI2PContext = new BiLinFormContext(massIntTDEFPSI2P, MASS);
      massIntTDEFPSI2PContext->SetEntities(actSDList, actSDList);
      massIntTDEFPSI2PContext->SetFeFunctions(feFunctions_[(SolutionType)(ACOU_TDEF_PSI_V_1 + iPole)], feFunctions_[formulation_]);
      assemble_->AddBiLinearForm(massIntTDEFPSI2PContext);
    } // end loop massIntTDEFPSI2PSI2, dampIntTDEFPSI2PSI2, stiffIntTDEFPSI2PSI2 and massInteTDEFPSI2P
  }

  template <UInt DIM, bool IS_COMPLEX>
  void AcousticPDE::DefineConvectiveIntegrators(RegionIdType actRegion, PtrParamNode curRegNode,
                                                shared_ptr<ElemList> actSDList, PtrCoefFct coeffM) {
    std::string flowId = curRegNode->Get("flowId")->As<std::string>();
    if (flowId != "") {
      std::cout << "Assigning convective integrators for AcousticPDE" << std::endl;
      if (complexFluidFormulation_) {
        EXCEPTION("Complex fluid and flow currently not allowed");
      }

      // Get result info object for flow
      shared_ptr<ResultInfo> flowInfo = GetResultInfo(MEAN_FLUIDMECH_VELOCITY);

      // Add the region information
      PtrParamNode flowNode = myParam_->Get("flowList")->GetByVal("flow", "name", flowId.c_str());

      bool fullForm = false;

      if (myParam_->Get("flowFormulation")->As<std::string>() == "withDivergence") {
        fullForm = true;
      }

      // Read coefficient flow coefficient function for this region and add it to flow functor
      PtrCoefFct regionFlow;
      std::set<UInt> definedDofs;
      bool coefUpdateGeo;
      ReadUserFieldValues(actSDList, flowNode, flowInfo->dofNames, flowInfo->entryType,
                          IS_COMPLEX, regionFlow, definedDofs, coefUpdateGeo);
      meanFlowCoef_->AddRegion(actRegion, regionFlow);

      PtrCoefFct divRegionFlow;
      PtrCoefFct divUFactors;

      if (fullForm) {
        ReadUserFieldValues(actSDList, flowNode, flowInfo->dofNames, flowInfo->entryType,
                            IS_COMPLEX, divRegionFlow, definedDofs, coefUpdateGeo);
        divRegionFlow->SetDerivativeOperation(CoefFunction::VECTOR_DIVERGENCE);
        divMeanFlowCoef_->AddRegion(actRegion, divRegionFlow);
      }

      // Create integrators
      BaseBDBInt *convectiveStiff = nullptr;
      BiLinearForm *convectiveStiffDivU = nullptr;
      BiLinearForm *convectiveDamp = nullptr;
      BiLinearForm *convectiveDampT = nullptr;
      BiLinearForm *convectiveDampDivU = nullptr;

      // Define aliases depending on IS_COMPLEX (template parameter)
      using ScalarType = typename std::conditional<IS_COMPLEX, Complex, Double>::type;
      using ConvectiveOperatorType = typename std::conditional<IS_COMPLEX, ConvectiveOperator<FeH1, DIM, 1, Complex>, ConvectiveOperator<FeH1, DIM, 1>>::type;
      using IdentityOperatorType = typename std::conditional<IS_COMPLEX, IdentityOperator<FeH1, DIM, 1, Complex>, IdentityOperator<FeH1, DIM, 1>>::type;

      convectiveDamp = new ABInt<ScalarType>(new IdentityOperatorType(), new ConvectiveOperatorType(),
                                             coeffM, 1.0, coefUpdateGeo);
      convectiveDampT = new ABInt<ScalarType>(new ConvectiveOperatorType(), new IdentityOperatorType(),
                                              coeffM, -1.0, coefUpdateGeo);
      convectiveStiff = new BBInt<ScalarType>(new ConvectiveOperatorType(), coeffM, -1.0, coefUpdateGeo);

      if (fullForm) {
        divUFactors = CoefFunction::Generate(mp_, IS_COMPLEX ? Global::COMPLEX : Global::REAL,
                                             CoefXprBinOp(mp_, coeffM, divRegionFlow, CoefXpr::OP_MULT));
        convectiveDampDivU = new BBInt<ScalarType>(new IdentityOperatorType(), divUFactors, 1.0, coefUpdateGeo);
        convectiveStiffDivU = new ABInt<ScalarType>(new IdentityOperatorType(), new ConvectiveOperatorType(),
                                                    divUFactors, 1.0, coefUpdateGeo);
      }

      convectiveStiff->SetBCoefFunctionOpB(meanFlowCoef_);
      convectiveStiff->SetName("convectiveStiffPierce");
      convectiveDamp->SetBCoefFunctionOpB(meanFlowCoef_);
      convectiveDamp->SetName("convectiveDampPierce");
      convectiveDampT->SetBCoefFunctionOpA(meanFlowCoef_);
      convectiveDampT->SetName("convectiveDampPierceTransposed");

      convectiveInts_[actRegion] = convectiveStiff;

      BiLinFormContext *convectiveContextStiff = new BiLinFormContext(convectiveStiff, STIFFNESS);
      BiLinFormContext *convectiveContextDamp = new BiLinFormContext(convectiveDamp, DAMPING);
      BiLinFormContext *convectiveContextDampT = new BiLinFormContext(convectiveDampT, DAMPING);

      convectiveContextDamp->SetEntities(actSDList, actSDList);
      convectiveContextDamp->SetFeFunctions(feFunctions_[formulation_], feFunctions_[formulation_]);
      convectiveContextDampT->SetEntities(actSDList, actSDList);
      convectiveContextDampT->SetFeFunctions(feFunctions_[formulation_], feFunctions_[formulation_]);
      convectiveContextStiff->SetEntities(actSDList, actSDList);
      convectiveContextStiff->SetFeFunctions(feFunctions_[formulation_], feFunctions_[formulation_]);
      assemble_->AddBiLinearForm(convectiveContextDamp);
      assemble_->AddBiLinearForm(convectiveContextDampT);
      assemble_->AddBiLinearForm(convectiveContextStiff);

      if (fullForm) {
        convectiveStiffDivU->SetBCoefFunctionOpB(meanFlowCoef_);
        convectiveStiffDivU->SetName("convectiveStiffPierceDivU");
        convectiveDampDivU->SetName("convectiveDampPierceDivU");
        BiLinFormContext *convectiveContextStiffDivU = new BiLinFormContext(convectiveStiffDivU, STIFFNESS);
        BiLinFormContext *convectiveContextDampDivU = new BiLinFormContext(convectiveDampDivU, DAMPING);

        convectiveContextDampDivU->SetEntities(actSDList, actSDList);
        convectiveContextDampDivU->SetFeFunctions(feFunctions_[formulation_], feFunctions_[formulation_]);
        convectiveContextStiffDivU->SetEntities(actSDList, actSDList);
        convectiveContextStiffDivU->SetFeFunctions(feFunctions_[formulation_], feFunctions_[formulation_]);

        assemble_->AddBiLinearForm(convectiveContextDampDivU);
        assemble_->AddBiLinearForm(convectiveContextStiffDivU);
      }
    }
  }

  void AcousticPDE::DefineNcIntegrators() {
    StdVector<NcInterfaceInfo>::iterator ncIt = ncInterfaces_.Begin(),
                                         endIt = ncInterfaces_.End();
    for (; ncIt != endIt; ++ncIt) {
      switch (ncIt->type) {
      case NC_MORTAR:
        if (dim_ == 2) {
          DefineMortarCoupling<2, 1>(formulation_, *ncIt);
        } else {
          DefineMortarCoupling<3, 1>(formulation_, *ncIt);
        }
        break;
      case NC_NITSCHE:
        if (dim_ == 2) {
          if (isComplex_)
            DefineNitscheCoupling<2, true>(*ncIt);
          else
            DefineNitscheCoupling<2, false>(*ncIt);
        } else { /* if (dim_ == 3) */
          if (isComplex_)
            DefineNitscheCoupling<3, true>(*ncIt);
          else
            DefineNitscheCoupling<3, false>(*ncIt);
        }
        break;
      default:
        EXCEPTION("Unknown type of ncInterface");
        break;
      }
    }
  }

  template <UInt DIM, bool IS_COMPLEX>
  void AcousticPDE::DefineNitscheCoupling(NcInterfaceInfo &iface) {

    // Define aliases depending on IS_COMPLEX (template parameter)
    auto complType = IS_COMPLEX ? Global::COMPLEX : Global::REAL;
    using ScalarType = typename std::conditional<IS_COMPLEX, Complex, Double>::type;

    shared_ptr<BaseNcInterface> ncIf = ptGrid_->GetNcInterface(iface.interfaceId);
    MortarInterface *nitscheIf = dynamic_cast<MortarInterface *>(ncIf.get());
    assert(nitscheIf);

    // in case of Nitsche coupling edge/face information is required
    this->ptGrid_->MapEdges();
    this->ptGrid_->MapFaces();

    // currently we have a moving formulation only for acoustics
    updatedGeo_ = updatedGeo_ || ncIf->IsMoving();

    // create new entity list
    shared_ptr<ElemList> actSDList = ncIf->GetElemList();

    // we set here the penalty factor
    Double beta = iface.nitscheFactor;

    // get material density
    PtrCoefFct dens = materials_[nitscheIf->GetPrimaryVolRegion()]->GetScalCoefFnc(DENSITY, complType);

    // here we consider the case of mechAcou coupling in the potential formulation
    PtrCoefFct factor = CoefFunction::Generate(mp_, complType, "1.0"); // *1, default
    if (isMechCoupled_ == true && formulation_ == ACOU_POTENTIAL) {    // *rho or *-rho
      // Important: In case of a general / quadratic EV problem, we must
      // ensure to have a "positive definite" matrix, i.e. we are not allowed
      // to multiply all matrices by -1! (but still multiply by density)
      std::string stringFac = (analysistype_ != EIGENFREQUENCY) ? "-1.0" : "1.0";
      factor = CoefFunction::Generate(mp_, complType, CoefXprBinOp(mp_, dens, stringFac, CoefXpr::OP_MULT));
    } else {
      if (IS_COMPLEX) // *1/rho
        factor = CoefFunction::Generate(mp_, complType, CoefXprBinOp(mp_, factor, dens, CoefXpr::OP_DIV));
    }

    // notation> assume the test function is called v
    BiLinearForm *penalty_v1_u1 = nullptr;
    BiLinearForm *penalty_v1_u2 = nullptr; // penalty_v2_u1 is covered by SetCounterPart(true)
    BiLinearForm *penalty_v2_u2 = nullptr;
    // now bilinear forms related to the normal derivatives
    // du1 refers to the normal derivative directing from 1 to 2
    BiLinearForm *flux_dv1_u1 = nullptr;
    BiLinearForm *flux_dv1_u2 = nullptr;
    BiLinearForm *flux_v1_du1 = nullptr;

    // NOTE: the algebraic system sets the system matrix to
    // nonSym  if any bilinear form with the same fctID1 and fctID2 is nonSym.
    // We set here the symmetric flag to true in the constructor
    // of the SurfaceNitscheABInt even though the bilinear form itself is
    // not symmetric. Nitsche formulation is basically sym due to the
    // set counterpart directive for the context.

    // first, assign BOperators that include solution u1 and and test funciton v1 of the primary side...
    // primary part of penalty BOperator (u1*v1)
    penalty_v1_u1 = new SurfaceNitscheABInt<ScalarType, ScalarType>(new SurfaceIdentityOperator<FeH1, DIM, 1>(),
                                                                    new SurfaceIdentityOperator<FeH1, DIM, 1>(),
                                                                    factor, beta, BiLinearForm::PRIM_PRIM,
                                                                    updatedGeo_, true, true);

    // primary part of flux BOperator
    flux_dv1_u1 = new SurfaceNitscheABInt<ScalarType, ScalarType>(new SurfaceNormalDerivOperator<FeH1, DIM, 1>(),
                                                                  new SurfaceIdentityOperator<FeH1, DIM, 1>(),
                                                                  factor, -1.0, BiLinearForm::PRIM_PRIM,
                                                                  updatedGeo_, true, false);

    // symmetrization term
    flux_v1_du1 = new SurfaceNitscheABInt<ScalarType, ScalarType>(new SurfaceIdentityOperator<FeH1, DIM, 1>(),
                                                                  new SurfaceNormalDerivOperator<FeH1, DIM, 1>(),
                                                                  factor, -1.0, BiLinearForm::PRIM_PRIM,
                                                                  updatedGeo_, true, false);

    // now, assign BOperators that include solution u1 of the primary and and test funciton v2 of the secondary side...
    // mixed part of penalty term (-u1*v2)
    penalty_v1_u2 = new SurfaceNitscheABInt<ScalarType, ScalarType>(new SurfaceIdentityOperator<FeH1, DIM, 1>(),
                                                                    new SurfaceIdentityOperator<FeH1, DIM, 1>(),
                                                                    factor, beta * -1.0, BiLinearForm::PRIM_SEC,
                                                                    updatedGeo_, true, true);

    // mixed part of flux term
    flux_dv1_u2 = new SurfaceNitscheABInt<ScalarType, ScalarType>(new SurfaceNormalDerivOperator<FeH1, DIM, 1>(),
                                                                  new SurfaceIdentityOperator<FeH1, DIM, 1>(),
                                                                  factor, 1.0, BiLinearForm::PRIM_SEC,
                                                                  updatedGeo_, true, false);

    // finally, the terms living purely on the secondary domain...
    // here we need to define only the penalty term as for the flux we enforce that (du1 = du2), so only the mixed flux term is required
    penalty_v2_u2 = new SurfaceNitscheABInt<ScalarType, ScalarType>(new SurfaceIdentityOperator<FeH1, DIM, 1>(),
                                                                    new SurfaceIdentityOperator<FeH1, DIM, 1>(),
                                                                    factor, beta, BiLinearForm::SEC_SEC,
                                                                    updatedGeo_, true, true);

    // Nitsche coupling matrix is a type of stiffness matrix
    FEMatrixType targetMatrix;
    if (ncIf->IsMoving()) {
      targetMatrix = STIFFNESS_UPDATE;
    } else {
      targetMatrix = STIFFNESS;
    }

    // now the BOperators are set, so define the contexts...
    SurfaceBiLinFormContext *penalty_v1_u1_Context = new SurfaceBiLinFormContext(penalty_v1_u1, targetMatrix, BiLinearForm::PRIM_PRIM);
    SurfaceBiLinFormContext *penalty_v2_u2_Context = new SurfaceBiLinFormContext(penalty_v2_u2, targetMatrix, BiLinearForm::SEC_SEC);
    SurfaceBiLinFormContext *penalty_v1_u2_Context = new SurfaceBiLinFormContext(penalty_v1_u2, targetMatrix, BiLinearForm::PRIM_SEC);
    SurfaceBiLinFormContext *flux_dv1_u1_Context = new SurfaceBiLinFormContext(flux_dv1_u1, targetMatrix, BiLinearForm::PRIM_PRIM);
    SurfaceBiLinFormContext *flux_v1_du1_Context = new SurfaceBiLinFormContext(flux_v1_du1, targetMatrix, BiLinearForm::PRIM_PRIM);
    SurfaceBiLinFormContext *flux_dv1_u2_Context = new SurfaceBiLinFormContext(flux_dv1_u2, targetMatrix, BiLinearForm::PRIM_SEC);
    // assign motion to the contexts
    penalty_v1_u1_Context->SetMotion(ncIf->IsMoving());
    penalty_v2_u2_Context->SetMotion(ncIf->IsMoving());
    penalty_v1_u2_Context->SetMotion(ncIf->IsMoving());
    flux_dv1_u1_Context->SetMotion(ncIf->IsMoving());
    flux_v1_du1_Context->SetMotion(ncIf->IsMoving());
    flux_dv1_u2_Context->SetMotion(ncIf->IsMoving());
    // assign names to the operators
    penalty_v1_u1->SetName("penalty_v1_u1");
    penalty_v2_u2->SetName("penalty_v2_u2");
    penalty_v1_u2->SetName("penalty_v1_u2");
    flux_dv1_u1->SetName("flux_dv1_u1");
    flux_v1_du1->SetName("flux_v1_du1");
    flux_dv1_u2->SetName("flux_dv1_u2");
    // assign the entity list to the operators
    penalty_v1_u1_Context->SetEntities(actSDList, actSDList);
    penalty_v2_u2_Context->SetEntities(actSDList, actSDList);
    penalty_v1_u2_Context->SetEntities(actSDList, actSDList);
    flux_dv1_u1_Context->SetEntities(actSDList, actSDList);
    flux_v1_du1_Context->SetEntities(actSDList, actSDList);
    flux_dv1_u2_Context->SetEntities(actSDList, actSDList);
    // assign fe functions to the operators
    penalty_v1_u1_Context->SetFeFunctions(feFunctions_[formulation_], feFunctions_[formulation_]);
    penalty_v2_u2_Context->SetFeFunctions(feFunctions_[formulation_], feFunctions_[formulation_]);
    penalty_v1_u2_Context->SetFeFunctions(feFunctions_[formulation_], feFunctions_[formulation_]);
    flux_dv1_u1_Context->SetFeFunctions(feFunctions_[formulation_], feFunctions_[formulation_]);
    flux_v1_du1_Context->SetFeFunctions(feFunctions_[formulation_], feFunctions_[formulation_]);
    flux_dv1_u2_Context->SetFeFunctions(feFunctions_[formulation_], feFunctions_[formulation_]);
    // assign symmetry
    penalty_v1_u2_Context->SetCounterPart(true);
    flux_dv1_u2_Context->SetCounterPart(true);
    // passt the contexts to the assembler
    assemble_->AddBiLinearForm(penalty_v1_u1_Context);
    assemble_->AddBiLinearForm(penalty_v2_u2_Context);
    assemble_->AddBiLinearForm(penalty_v1_u2_Context);
    assemble_->AddBiLinearForm(flux_dv1_u1_Context);
    assemble_->AddBiLinearForm(flux_v1_du1_Context);
    assemble_->AddBiLinearForm(flux_dv1_u2_Context);
    // register integrators
    ncIf->RegisterIntegrator(penalty_v1_u1_Context);
    ncIf->RegisterIntegrator(penalty_v2_u2_Context);
    ncIf->RegisterIntegrator(penalty_v1_u2_Context);
    ncIf->RegisterIntegrator(flux_dv1_u1_Context);
    ncIf->RegisterIntegrator(flux_v1_du1_Context);
    ncIf->RegisterIntegrator(flux_dv1_u2_Context);

    // check for eulerian formulation of moving grid
    DefineEulerianSystem<DIM, 1>(formulation_, iface);
  }

  void AcousticPDE::SetIntegratorContext(BaseBDBInt *&stiffInt, BaseBDBInt *&massInt, RegionIdType actRegion, shared_ptr<ElemList> &actSDList, PtrCoefFct &coeffK, PtrCoefFct &coeffM) {
    // stiffness context
    stiffInt->SetName("LaplaceIntegrator");
    BiLinFormContext *stiffIntDescr = new BiLinFormContext(stiffInt, STIFFNESS);
    feFunctions_[formulation_]->AddEntityList(actSDList);
    stiffIntDescr->SetEntities(actSDList, actSDList);
    stiffIntDescr->SetFeFunctions(feFunctions_[formulation_], feFunctions_[formulation_]);
    stiffInt->SetFeSpace(feFunctions_[formulation_]->GetFeSpace());

    // mass context
    massInt->SetName("MassIntegrator");
    massInt->SetFeSpace(feFunctions_[formulation_]->GetFeSpace());
    BiLinFormContext *massContext = new BiLinFormContext(massInt, MASS);
    massContext->SetEntities(actSDList, actSDList);
    massContext->SetFeFunctions(feFunctions_[formulation_], feFunctions_[formulation_]);

    // Check for Rayleigh damping
    if (dampingList_[actRegion] == RAYLEIGH || dampingList_[actRegion] == ADAPTED_LOSS_TANGENS_DELTA || dampingList_[actRegion] == GLOBAL_RAYLEIGH) {
      if (complexFluidFormulation_)
      EXCEPTION("Complex fluid region and Rayleigh damping not allowed!!");
      RaylDampingData &actDamp = regionRaylDamping_[actRegion];
      massContext->SetSecDestMat(DAMPING, actDamp.alpha);
      stiffIntDescr->SetSecDestMat(DAMPING, actDamp.beta);
    }

    // add bilinear form descriptions to the assembler
    assemble_->AddBiLinearForm(stiffIntDescr);
    assemble_->AddBiLinearForm(massContext);
    // Add bdb-integrator to global list, as we need them later for calculation of postprocessing results.
    bdbInts_.insert(std::pair<RegionIdType, BaseBDBInt *>(actRegion, stiffInt));
    // Add mass-integrator to global list, as we need them later for calculation of postprocessing results
    massInts_[actRegion] = massInt;
  }

  void AcousticPDE::DefineSurfaceIntegrators( ){
    //========================================================================================
    // ABC boundaries
    //========================================================================================
    PtrParamNode bcNode = myParam_->Get( "bcsAndLoads", ParamNode::PASS );
    LOG_DBG(acousticpde) << "Define Surface Integrator BEGIN"<< "\n";

    if( bcNode ) {
      //========================================================================================
      // Absorbing boundary condition (ABC)
      //========================================================================================
      ParamNodeList abcNodes = bcNode->GetList( "absorbingBCs" );
      LOG_DBG(acousticpde) << "ABCs count :  " << abcNodes.GetSize() <<  "\n" ;

      std::set<RegionIdType> adjVolRegion;
      for( UInt i = 0; i < abcNodes.GetSize(); i++ ) {
        std::string regionName = abcNodes[i]->Get("name")->As<std::string>();
        shared_ptr<EntityList> actSDList =  ptGrid_->GetEntityList( EntityList::SURF_ELEM_LIST,regionName );
        std::string volRegName = abcNodes[i]->Get("volumeRegion")->As<std::string>();
        LOG_DBG(acousticpde) << "ABCs volRegName :  " << volRegName <<  "\n" ;

        RegionIdType aRegion = ptGrid_->GetRegion().Parse(volRegName);
        UInt iRegion = regions_.Find(aRegion);

        //check, if region has complex fluid
        PtrParamNode curRegNode =
            myParam_->Get("regionList")->GetByVal("region","name",volRegName.c_str());
        std::string useRationalAppr = curRegNode->Get("useRationalMatApproximation")->As<std::string>();

        // c0 = sqrt(bulk_modulus / density)
        PtrCoefFct dens;
        PtrCoefFct blk;
        PtrCoefFct constOne = CoefFunction::Generate(mp_, Global::COMPLEX, "1.0", "0.0");
        PtrCoefFct omegaTrg;

        if (complexFluidFormulation_ && useRationalAppr == "no") {
          dens = materials_[aRegion]->GetScalCoefFnc(DENSITY, Global::COMPLEX);
          blk = materials_[aRegion]->GetScalCoefFnc(ACOU_BULK_MODULUS, Global::COMPLEX);
        }
        else if (complexFluidFormulation_ && useRationalAppr == "yes") {
          // get the current frequency
          PtrCoefFct freqCoef;
          freqCoef = CoefFunction::Generate(mp_, Global::REAL, "f");
          EvalTDEFRationalFncs(iRegion, freqCoef, TDEFcoeffs_);
          dens = CoefFunction::Generate(mp_, Global::COMPLEX,
                                        CoefXprBinOp(mp_, constOne, invTDEFDens_, CoefXpr::OP_DIV));
          blk = CoefFunction::Generate(mp_, Global::COMPLEX,
                                       CoefXprBinOp(mp_, constOne, invTDEFBlk_, CoefXpr::OP_DIV));
        }
        else if (timeDomainEqFluidFormulation_ && isTDEFReg_[iRegion]) {
          // ABC adjacent to TDEF region
          LocPointMapped lpm;

          std::cout << "Establishing narrowband ABC for TDEF region for the defined target frequency." << std::endl;
          std::cout << "ATTENTION: for dispersive media, the frequency content depends on the distance from the source! \n"
                    << std::endl
                    << std::endl;
          Double ftrg = abcNodes[i]->Get("targetFrequency")->As<UInt>();

          PtrCoefFct targetFreg = CoefFunction::Generate(mp_, Global::REAL, std::to_string(ftrg));
          omegaTrg = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, targetFreg, CoefFunction::Generate(mp_, Global::REAL, "2*pi"), CoefXpr::OP_MULT));
          // get density and bulk modulus at target frequency
          EvalTDEFRationalFncs(iRegion, targetFreg, TDEFcoeffs_);
          dens = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, constOne, invTDEFDens_, CoefXpr::OP_DIV));
          blk = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, constOne, invTDEFBlk_, CoefXpr::OP_DIV));
          LOG_DBG(acousticpde) << "Narrowband-NRBC: Bulk modulus at f=" << targetFreg->ToString() << "Hz: K=" << blk->ToString() << "Pa \n";
          LOG_DBG(acousticpde) << "Narrowband-NRBC: Density at f=" << targetFreg->ToString() << "Hz: rho=" << dens->ToString() << "kg/m3 \n";
        }
        else {
          dens = materials_[aRegion]->GetScalCoefFnc( DENSITY, Global::REAL );
          blk = materials_[aRegion]->GetScalCoefFnc( ACOU_BULK_MODULUS, Global::REAL );
        }

        LOG_DBG(acousticpde) << "ABC: dens = " << dens->ToString() << "\n";
        LOG_DBG(acousticpde) << "ABC: blk  = " << blk->ToString() << "\n";

        PtrCoefFct c0;
        PtrCoefFct coeffStiffTDEF;

        //check for temperature dependency
        curRegNode = myParam_->Get("regionList")->GetByVal("region","name",volRegName.c_str());
        std::string tempId = curRegNode->Get("temperatureId")->As<std::string>();
        if (tempId != "") {
          // Get result info object for flow
          shared_ptr<ResultInfo> tempInfo = GetResultInfo(HEAT_MEAN_TEMPERATURE);

          //Add the region information
          PtrParamNode tempNode = myParam_->Get("temperatureList")->GetByVal("temperature","name",tempId.c_str());

          // Read flow coefficient function for the region and add it to flow functor
          PtrCoefFct regionTemp;
          std::set<UInt> definedDofs;
          //we assume that the temperature values are real (not complex!)
          ReadUserFieldValues( actSDList, tempNode, tempInfo->dofNames, tempInfo->entryType,
                             false, regionTemp, definedDofs, updatedGeo_ );
          // gasR=287.058 J/kg K   ... universal gas constant
          // kappa=1.402, adabatic exponent for air
          PtrCoefFct constVal = CoefFunction::Generate( mp_, Global::REAL, "402.4553160");
          c0 = CoefFunction::Generate( mp_,  Global::REAL,
                           CoefXprUnaryOp( mp_, CoefXprBinOp(mp_, constVal, regionTemp, CoefXpr::OP_MULT),
                           CoefXpr::OP_SQRT) );
          LOG_DBG(acousticpde) << "Define Surface Integrators standard c0 =" << c0->ToString() << "\n";
        }
        else {
        	if (complexFluidFormulation_ ) {
        		// correct factor is 1/(density*c0) = 1/sqrt(dens*compressModulus)
        		// so we modify c0 to : sqrt(density*compressModulus)
        		// division is done later
        		// note that c0 in this case does not represent the speed of sound
        		c0 = CoefFunction::Generate( mp_,  Global::COMPLEX, CoefXprUnaryOp(mp_, CoefXprBinOp(mp_, blk, dens,
                                            CoefXpr::OP_MULT), CoefXpr::OP_SQRT) );
        	}
          else if (timeDomainEqFluidFormulation_ && isTDEFReg_[iRegion]) {
            // Narrow-band ABC for TDEF region according to Abdulkareem et al. 2018
            c0 = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprUnaryOp(mp_, CoefXprBinOp(mp_, blk, dens, CoefXpr::OP_DIV), CoefXpr::OP_SQRT));
            LOG_DBG(acousticpde) << "Narrowband-NRBC: Complex wave speed: c = " << c0->ToString() << "m/s \n";
          }
          else
            c0 = CoefFunction::Generate(mp_, Global::REAL, CoefXprUnaryOp(mp_, CoefXprBinOp(mp_, blk, dens, CoefXpr::OP_DIV), CoefXpr::OP_SQRT));
        }
        LOG_DBG(acousticpde) << "Def Surface Integrator:  c0 =" << c0->ToString() << "\n";


        // the following part was missing which is why abc did not function for acouPotential + mechanic
        // if pde couples with mechanic, we have to multiply the density by -1
        PtrCoefFct mechAcouFactor;
        CalcMechAcouFac(mechAcouFactor,dens);
  
        LOG_DBG(acousticpde) << "Def Surface Integrator: mechAcouFactor =" << mechAcouFactor->ToString() << "\n";

        PtrCoefFct coeffDamp;
        if ( sosAtLaplace_ ) {
          // factor for damping matrix: mechAcouFactor * c0
          coeffDamp = CoefFunction::Generate( mp_, Global::REAL,
                         			CoefXprBinOp(mp_, mechAcouFactor, c0, CoefXpr::OP_MULT ) );
        }
        else if (timeDomainEqFluidFormulation_) {
          CalcCoefsTDEFABC(iRegion, mechAcouFactor, c0, dens, omegaTrg, coeffDamp, coeffStiffTDEF);
        }
        else if (complexFluidFormulation_) {
          // factor for damping matrix: factor / c0
          coeffDamp = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, mechAcouFactor, c0, CoefXpr::OP_DIV));
        }
        else {
          // factor for damping matrix: mechAcouFactor / c0
         if (complexFluidFormulation_ )
           coeffDamp = CoefFunction::Generate( mp_, Global::COMPLEX, CoefXprBinOp(mp_, mechAcouFactor, c0, CoefXpr::OP_DIV ) );
         else
           coeffDamp = CoefFunction::Generate( mp_, Global::REAL, CoefXprBinOp(mp_, mechAcouFactor, c0, CoefXpr::OP_DIV ) );

        }
        LOG_DBG(acousticpde) << "Define Surface Integrator: coeffDamp =" << coeffDamp->ToString() << "\n";
        
        BiLinearForm * abcInt = NULL;
        if (complexFluidFormulation_) {
          if( dim_ == 2 )
            abcInt = new BBInt<Complex>(new IdentityOperator<FeH1,2,1, Complex>(), coeffDamp, 1.0, updatedGeo_ );
          else
            abcInt = new BBInt<Complex>(new IdentityOperator<FeH1,3,1, Complex>(), coeffDamp, 1.0, updatedGeo_ );
        }

        else {
          if( dim_ == 2 )
            abcInt = new BBInt<>(new IdentityOperator<FeH1,2,1>(), coeffDamp, 1.0, updatedGeo_ );
          else
        	abcInt = new BBInt<>(new IdentityOperator<FeH1,3,1>(), coeffDamp, 1.0, updatedGeo_ );

        }

        FEMatrixType targetMatrix = DAMPING;
        if(updatedGeo_){
          targetMatrix = DAMPING_UPDATE;
        } else if(this->analysistype_ == HARMONIC25D) {
          targetMatrix = DAMPING_AUX;
        }

        abcInt->SetName("abcIntegrator");
        BiLinFormContext *abcContext = new BiLinFormContext(abcInt, targetMatrix );

        abcContext->SetEntities( actSDList, actSDList );
        abcContext->SetFeFunctions( feFunctions_[formulation_] , feFunctions_[formulation_]);
        feFunctions_[formulation_]->AddEntityList( actSDList );
        assemble_->AddBiLinearForm( abcContext );

        // ABC on a TDEF-volume region requires an additional stiffness term
        // rho dp/dn = -  rho/c dp/dt - rho beta p; with c = real(c_complex) , beta = - imag(c_complex)
        if (timeDomainEqFluidFormulation_ && isTDEFReg_[iRegion]) {
          DefineTDEFABCSurfaceIntegrators(coeffStiffTDEF, aRegion, adjVolRegion, actSDList, TDEFcoeffs_);
        }
      } // end abc

      //========================================================================================
      // Impedance boundaries
      // TODO: implement impedance BC
      //========================================================================================
      ParamNodeList impedNodes = bcNode->GetList( "impedance" );

      for( UInt i = 0; i < impedNodes.GetSize(); ++i ) {
        BiLinearForm * impedInt = NULL;
        std::string regionName = impedNodes[i]->Get("name")->As<std::string>();
        shared_ptr<EntityList> actSDList =  ptGrid_->GetEntityList( EntityList::SURF_ELEM_LIST, regionName );
        std::string volRegName = impedNodes[i]->Get("volumeRegion")->As<std::string>();
        RegionIdType aRegion = ptGrid_->GetRegion().Parse(volRegName);
        bool isNormalised = impedNodes[i]->Get("isNormalised")->As<bool>();
        std::string impId = impedNodes[i]->Get("impedanceId")->As<std::string>();

        shared_ptr<CoefFunctionImpedanceModel<Complex> >
                   Z_impMod(new CoefFunctionImpedanceModel<Complex>(mp_, materials_[aRegion], isNormalised));

        PtrParamNode impListNodes = myParam_->Get( "impedanceList", ParamNode::PASS );
        if ( !impListNodes ) {
          EXCEPTION("No impedance list in *.xml");
        }
        ParamNodeList impChildren = impListNodes->GetChildren();
        UInt i_imp = 0;
        for ( ; i_imp < impChildren.GetSize(); ++i_imp) {
          std::string impType = impChildren[i_imp]->GetName();
          std::string actId = impChildren[i_imp]->Get("id")->As<std::string>();

          if (actId != impId) {
            continue;
          }
          if (impType == "mpp") {
            std::string impSubType = impChildren[i_imp]->Get("subtype")->As<std::string>();
            if (impSubType == "slit") {
              Z_impMod->GenerateSlitMpp(impChildren[i_imp]);
            } else if (impSubType == "circ") {
              Z_impMod->GenerateCircMpp(impChildren[i_imp]);
            } else {
              EXCEPTION("No such impedance type: " << impType << " with subtype: " << impSubType);
            }
          } else if (impType == "interpolImpedance"){
            Z_impMod->GenerateInterpolImpedance(impChildren[i_imp]);
          } else if (impType == "fctImpedance") {
            Z_impMod->GenerateImpedanceFct(impChildren[i_imp]);
          } else {
            EXCEPTION("No such impedance type: " << impType);
          }
          break;

        }
        if (i_imp == impChildren.GetSize()) {
          EXCEPTION("No such impedance id found: " << impId);
        }

        // for complexFluid (inhomogeneous) PDE, we need to multiply the impedance with rho for alignment with the different PDE structure
        // Z_impMode normal = rho**2 / Z
        // Z_impMode complex = rho / Z
        PtrCoefFct Z_impScaled = Z_impMod;
        if ( complexFluidFormulation_ ){
          PtrCoefFct dens = materials_[aRegion]->GetScalCoefFnc( DENSITY, Global::COMPLEX );
          Z_impScaled = CoefFunction::Generate( mp_, Global::COMPLEX, CoefXprBinOp(mp_, Z_impScaled, dens, CoefXpr::OP_DIV ));
        }

        if( dim_ == 2 ) {
          impedInt = new BBInt<Complex>(new IdentityOperator<FeH1,2,1, Complex>(), Z_impScaled, 1.0, updatedGeo_ );
        } else {
          impedInt = new BBInt<Complex>(new IdentityOperator<FeH1,3,1, Complex>(), Z_impScaled, 1.0, updatedGeo_ );
        }

        FEMatrixType targetMatrix = DAMPING;
        if(updatedGeo_){
          targetMatrix = DAMPING_UPDATE;
        }

        impedInt->SetName("impedIntegrator");
        BiLinFormContext *impedContext = new BiLinFormContext(impedInt, targetMatrix );

        impedContext->SetEntities( actSDList, actSDList );
        impedContext->SetFeFunctions( feFunctions_[formulation_] , feFunctions_[formulation_]);
        feFunctions_[formulation_]->AddEntityList( actSDList );
        assemble_->AddBiLinearForm( impedContext );
      } // end impedance bc

      //========================================================================================
      // boundary Layers
      //========================================================================================
      ParamNodeList blNodes = bcNode->GetList( "boundaryLayer" );

      for( UInt i = 0; i < blNodes.GetSize(); ++i ) {
        if ( !(formulation_ == ACOU_PRESSURE) || !(this->analysistype_ == HARMONIC) ) EXCEPTION("bounadryLayer only available for harmonic pressure formulation");

        std::string regionName = blNodes[i]->Get("name")->As<std::string>();
        shared_ptr<EntityList> actSDList =  ptGrid_->GetEntityList( EntityList::SURF_ELEM_LIST, regionName );
        std::string volRegName = blNodes[i]->Get("volumeRegion")->As<std::string>();
        RegionIdType volRegion = ptGrid_->GetRegion().Parse(volRegName);

        PtrCoefFct rho0 = materials_[volRegion]->GetScalCoefFnc( DENSITY, Global::REAL );
        PtrCoefFct K = materials_[volRegion]->GetScalCoefFnc( ACOU_BULK_MODULUS, Global::REAL );
        PtrCoefFct cp = CoefFunction::Generate( mp_,  Global::REAL, blNodes[i]->Get("cp")->As<std::string>() );
        PtrCoefFct cv = CoefFunction::Generate( mp_,  Global::REAL, blNodes[i]->Get("cv")->As<std::string>() );
        PtrCoefFct nu = CoefFunction::Generate( mp_,  Global::REAL, blNodes[i]->Get("nu")->As<std::string>() );
        PtrCoefFct k = CoefFunction::Generate( mp_,  Global::REAL, blNodes[i]->Get("k")->As<std::string>() );

        LOG_DBG(acousticpde) << "Define Surface Integrator: regionName =" << regionName << "\n";
        LOG_DBG(acousticpde) << "Define Surface Integrator: volRegName =" << volRegName << "\n";
        LOG_DBG(acousticpde) << "Define Surface Integrator: rho0       =" << rho0->ToString() << "\n";
        LOG_DBG(acousticpde) << "Define Surface Integrator: K          =" << K->ToString() << "\n";
        LOG_DBG(acousticpde) << "Define Surface Integrator: k          =" << k->ToString() << "\n";

        PtrCoefFct omegaHalv = CoefFunction::Generate( mp_,  Global::REAL, "pi*f");//
        // deltaV = sqrt( 2*nu/omega )
        PtrCoefFct deltaV = CoefFunction::Generate(mp_,Global::REAL, CoefXprUnaryOp(mp_, CoefXprBinOp(mp_, nu, omegaHalv, CoefXpr::OP_DIV ) , CoefXpr::OP_SQRT));
        // deltaT = sqrt( 2*k/(omega*rho0*cp) )
        PtrCoefFct deltaT = CoefFunction::Generate(mp_,Global::COMPLEX,
            CoefXprUnaryOp(mp_,
                CoefXprBinOp(mp_,
                    k, CoefXprBinOp(mp_,CoefXprBinOp(mp_, cp, rho0, CoefXpr::OP_MULT ),omegaHalv, CoefXpr::OP_MULT ),
                    CoefXpr::OP_DIV ),
                    CoefXpr::OP_SQRT));
        // 1-i ... common factor for both terms
        PtrCoefFct oneMinusI =  CoefFunction::Generate( mp_,  Global::COMPLEX, "1", "-1");
        // gamma = cp/cv -> gamma-1 = (cp-cv)/cv
        PtrCoefFct gammaMinusOne = CoefFunction::Generate(mp_,Global::COMPLEX, CoefXprBinOp(mp_, CoefXprBinOp(mp_,cp,cv,CoefXpr::OP_SUB), cv, CoefXpr::OP_DIV));

        // Mass matrix
        PtrCoefFct coefM =  CoefFunction::Generate(mp_,Global::COMPLEX,
            CoefXprBinOp(mp_,
                CoefXprBinOp(mp_,
                    CoefXprBinOp(mp_, deltaT, CoefXprBinOp(mp_, rho0, K, CoefXpr::OP_DIV), CoefXpr::OP_MULT),
                    oneMinusI,
                    CoefXpr::OP_MULT),
                    gammaMinusOne,
                    CoefXpr::OP_MULT)
        );
        BiLinearForm * blmInt = NULL;
        if( dim_ == 2 ) {
          blmInt = new BBInt<Complex>(new IdentityOperator<FeH1,2,1, Complex>(), coefM, 0.5, updatedGeo_ );
        } else {
          blmInt = new BBInt<Complex>(new IdentityOperator<FeH1,3,1, Complex>(), coefM, 0.5, updatedGeo_ );
        }
        blmInt->SetName("blMassIntegrator");
        BiLinFormContext *blmContext = new BiLinFormContext(blmInt, MASS);
        blmContext->SetEntities(actSDList, actSDList);
        blmContext->SetFeFunctions(feFunctions_[formulation_],feFunctions_[formulation_]);
        feFunctions_[formulation_]->AddEntityList(actSDList);
        assemble_->AddBiLinearForm(blmContext);

        LOG_DBG(acousticpde) << "Define Surface Integrator boundary layer: coefM =" << coefM->ToString() << "\n";

        // Stiffness matrix
        PtrCoefFct coefK = CoefFunction::Generate(mp_,Global::COMPLEX, CoefXprBinOp(mp_, deltaV, oneMinusI, CoefXpr::OP_MULT) );
        BiLinearForm * blkInt = NULL;
        std::set<RegionIdType> volRegions;
        volRegions.insert(volRegion);
        if( dim_ == 2 ) {
          blkInt = new BBInt<Complex>(new GradientOperator<FeH1,2,1, Complex>(), coefK, -0.5, updatedGeo_ );
        } else {
          blkInt = new BBInt<Complex>(new GradientOperator<FeH1,3,1, Complex>(), coefK, -0.5, updatedGeo_ );
        }
        blkInt->SetName("blStiffnessIntegrator");
        BiLinFormContext *blkContext = new BiLinFormContext(blkInt, STIFFNESS);
        blkContext->SetEntities(actSDList, actSDList);
        blkContext->SetFeFunctions(feFunctions_[formulation_], feFunctions_[formulation_]);
        feFunctions_[formulation_]->AddEntityList(actSDList);
        assemble_->AddBiLinearForm(blkContext);

        LOG_DBG(acousticpde) << "Define Surface Integrator boundary layer: coefK =" << coefK->ToString() << "\n";
      } // end boundary Layers

      //========================================================================================
      // blochPeriodic
      //========================================================================================
      ParamNodeList blochNodesList = bcNode->GetList("blochPeriodic");
      for (PtrParamNode& node : blochNodesList)
      {
        if(formulation_ != ACOU_PRESSURE)
          throw Exception("blochPeriodic only tested for acouPressure!");

        // probably also works for normal acoustics
        if(!complexFluidFormulation_)
          throw Exception("blochPeriodic only tested for complexFluid!");

        std::string str_value = node->Get("factor_value")->As<std::string>();
        std::string str_phase = node->Get("factor_phase")->As<std::string>();
        std::string formulation = node->Get("formulation")->As<std::string>();
        
        // propagation factor \gamma from xml-file
        std::string str_real, str_imag;
        str_real = AmplPhaseToReal(str_value, str_phase, true);
        str_imag = AmplPhaseToImag(str_value, str_phase, true);
        
        PtrCoefFct factor = CoefFunction::Generate(mp_, Global::COMPLEX, str_real, str_imag);
        PtrCoefFct one = CoefFunction::Generate(mp_, Global::REAL, "1.0", "0.0");
        
        ParamNodeList regionsList = node->GetList("region");
        for (PtrParamNode& regnode : regionsList)
        {
          std::string ncRegionName = regnode->Get("name")->As<std::string>();
          shared_ptr<BaseNcInterface> ncIf = ptGrid_->GetNcInterface(ptGrid_->GetNcInterfaceId(ncRegionName));
          if (!ncIf)
          {
            EXCEPTION("No interface with the name '" << ncRegionName << "' found!");
          }
          shared_ptr<MortarInterface> mortarIf = dynamic_pointer_cast<MortarInterface>(ncIf);
          assert(mortarIf);
          
          if (formulation == "Mortar")
          {
            shared_ptr<SurfElemList> surfPrmGrid(new SurfElemList(ptGrid_)), surfSndGrid(new SurfElemList(ptGrid_));
            surfPrmGrid->SetRegion(mortarIf->GetPrimarySurfRegion());
            surfSndGrid->SetRegion(mortarIf->GetSecondarySurfRegion());
            
            // --- Set the approximation for Lagrange Multipliers for the current region ---
            RegionIdType regId = surfSndGrid->GetRegion();
            std::string polyId = regnode->Get("polyId")->As<std::string>();
            std::string integId = regnode->Get("integId")->As<std::string>();
            feFunctions_[LAGRANGE_MULT]->GetFeSpace()->SetRegionApproximation(regId, polyId, integId);
            
            BiLinearForm* intOne1 = nullptr;
            BiLinearForm* intFactor1 = nullptr;
            BiLinearForm* intOne2 = nullptr;
            BiLinearForm* intFactor2 = nullptr;

            if (dim_ == 2)
            {
              intOne1 = new SurfaceMortarABInt<Complex, Complex>(new IdentityOperator<FeH1, 2, 1, Complex>(),
                    new IdentityOperator<FeH1, 2, 1, Complex>(),
                    one, 1.0,
                    mortarIf->GetSecondaryVolRegion(), mortarIf->GetPrimaryVolRegion(),
                    mortarIf->IsCoplanar(), updatedGeo_, BiLinearForm::SEC_PRIM);
              intFactor1 = new BBInt<Complex, Complex>(new IdentityOperator<FeH1, 2, 1, Complex>(), factor, -1.0, updatedGeo_);
              intOne2 = new BBInt<Complex, Complex>(new IdentityOperator<FeH1, 2, 1, Complex>(), one, 1.0, updatedGeo_);
              intFactor2 = new SurfaceMortarABInt<Complex, Complex>(new IdentityOperator<FeH1, 2, 1, Complex>(),
                      new IdentityOperator<FeH1, 2, 1, Complex>(),
                      factor, -1.0,
                      mortarIf->GetPrimaryVolRegion(), mortarIf->GetSecondaryVolRegion(),
                      mortarIf->IsCoplanar(), updatedGeo_, BiLinearForm::PRIM_SEC);
            }
            else
            {
              intOne1 = new SurfaceMortarABInt<Complex, Complex>(new IdentityOperator<FeH1, 3, 1, Complex>(),
                      new IdentityOperator<FeH1, 3, 1, Complex>(),
                      one, 1.0,
                      mortarIf->GetSecondaryVolRegion(), mortarIf->GetPrimaryVolRegion(),
                      mortarIf->IsCoplanar(), updatedGeo_, BiLinearForm::SEC_PRIM);
              intFactor1 = new BBInt<Complex, Complex>(new IdentityOperator<FeH1, 3, 1, Complex>(), factor, -1.0, updatedGeo_);
              intOne2 = new BBInt<Complex, Complex>(new IdentityOperator<FeH1, 3, 1, Complex>(), one, 1.0, updatedGeo_);
              intFactor2 = new SurfaceMortarABInt<Complex, Complex>(new IdentityOperator<FeH1, 3, 1, Complex>(),
                      new IdentityOperator<FeH1, 3, 1, Complex>(),
                      factor, -1.0,
                      mortarIf->GetPrimaryVolRegion(), mortarIf->GetSecondaryVolRegion(),
                      mortarIf->IsCoplanar(), updatedGeo_, BiLinearForm::PRIM_SEC);
            }
            
            intOne1->SetName("primary1Acou");
            intFactor1->SetName("secondary1Acou");
            intOne2->SetName("primary2Acou");
            intFactor2->SetName("secondary2Acou");
            
            // (1) weak form of the periodic boundary conditions
            BiLinFormContext* lagPotContPrm = new BiLinFormContext(intFactor1, STIFFNESS);
            NcBiLinFormContext* lagPotContSnd = new NcBiLinFormContext(intOne1, STIFFNESS);
            
            // I dont understand why this is different to the ElecPDE, but it works...
            lagPotContPrm->SetEntities(surfSndGrid, surfSndGrid);
            lagPotContSnd->SetEntities(mortarIf->GetElemList(), mortarIf->GetElemList());
            lagPotContPrm->SetFeFunctions(feFunctions_[LAGRANGE_MULT], feFunctions_[ACOU_PRESSURE]);
            lagPotContSnd->SetFeFunctions(feFunctions_[LAGRANGE_MULT], feFunctions_[ACOU_PRESSURE]);
            
            assemble_->AddBiLinearForm(lagPotContPrm);
            assemble_->AddBiLinearForm(lagPotContSnd);
            
            // (2) weak form of the boundary integrals of the PDE
            BiLinFormContext* potLagContPrm = new BiLinFormContext(intOne2, STIFFNESS);
            NcBiLinFormContext* potLagContSnd = new NcBiLinFormContext(intFactor2, STIFFNESS);
            
            potLagContPrm->SetEntities(surfSndGrid, surfSndGrid);
            potLagContSnd->SetEntities(mortarIf->GetElemList(), mortarIf->GetElemList());
            potLagContPrm->SetFeFunctions(feFunctions_[ACOU_PRESSURE], feFunctions_[LAGRANGE_MULT]);
            potLagContSnd->SetFeFunctions(feFunctions_[ACOU_PRESSURE], feFunctions_[LAGRANGE_MULT]);
            
            assemble_->AddBiLinearForm(potLagContPrm);
            assemble_->AddBiLinearForm(potLagContSnd);
            
            feFunctions_[LAGRANGE_MULT]->AddEntityList(surfSndGrid);
            feFunctions_[ACOU_PRESSURE]->AddEntityList(surfPrmGrid);
            feFunctions_[ACOU_PRESSURE]->AddEntityList(surfSndGrid);
          } // end mortar
          else if (formulation == "Nitsche")
          {
            EXCEPTION("Nitsche not implemented, use Mortar instead!");
          }
          else 
          {
            EXCEPTION("Unknown formulation: '" << formulation << "'!");
          }
        }
      } // end blochPeriodic
    } // end if ( bcNode )
    LOG_DBG(acousticpde) << "Define Surface Integrator END"<< "\n";
  } // DefineSurfaceIntegrators

  void AcousticPDE::DefineTDEFABCSurfaceIntegrators(PtrCoefFct &coeffStiffTDEF, RegionIdType &aRegion, std::set<RegionIdType> &adjVolRegion, shared_ptr<EntityList> &actSDList, rationalCoeffsTDEF &TDEFcoeffs)
  {
    UInt iRegion = regions_.Find(aRegion);
    BiLinearForm *abcIntSiffTDEF = NULL;

    if (dim_ == 2)
      abcIntSiffTDEF = new BBInt<>(new IdentityOperator<FeH1, 2, 1>(), coeffStiffTDEF, 1.0, updatedGeo_);
    else
      abcIntSiffTDEF = new BBInt<>(new IdentityOperator<FeH1, 3, 1>(), coeffStiffTDEF, 1.0, updatedGeo_);

    FEMatrixType targetMatrix = STIFFNESS;
    if (updatedGeo_) {
      targetMatrix = STIFFNESS_UPDATE;
    }
    abcIntSiffTDEF->SetName("abcIntegratorStiffnTermTDEF");
    BiLinFormContext *abcContextStiffBeta = new BiLinFormContext(abcIntSiffTDEF, targetMatrix);

    abcContextStiffBeta->SetEntities(actSDList, actSDList);
    abcContextStiffBeta->SetFeFunctions(feFunctions_[formulation_], feFunctions_[formulation_]);
    feFunctions_[formulation_]->AddEntityList(actSDList);
    assemble_->AddBiLinearForm(abcContextStiffBeta);

    // Loop over REAL poles of inverse density (specific volume)
    // (We do not have to care for the aux variables related to the inverse blk modulus)
    for (UInt iPole = 0; iPole < TDEFpoleNumber_.realV[iRegion]; iPole++) {
      adjVolRegion.insert(aRegion);

      BiLinearForm *stiffIntTDEFABCPDPHI1 = NULL;
      if (dim_ == 2) {
        stiffIntTDEFABCPDPHI1 = new SurfaceABInt<>(new IdentityOperator<FeH1, 2, 1>(),
                                                   new SurfaceNormalDerivOperator<FeH1, 2, 1>(),
                                                   TDEFcoeffs.AV[iPole], -1.0, adjVolRegion, updatedGeo_);
      }
      else {
        stiffIntTDEFABCPDPHI1 = new SurfaceABInt<>(new IdentityOperator<FeH1, 3, 1>(),
                                                   new SurfaceNormalDerivOperator<FeH1, 3, 1>(),
                                                   TDEFcoeffs.AV[iPole], -1.0, adjVolRegion, updatedGeo_);
      }
      FEMatrixType targetMatrix = STIFFNESS;
      if (updatedGeo_) {
        targetMatrix = STIFFNESS_UPDATE;
      }
      stiffIntTDEFABCPDPHI1->SetName("abcIntegratorStiffnTermTDEFAUX" + std::to_string(iPole));
      stiffIntTDEFABCPDPHI1->SetUseVolEqnB(true);
      BiLinFormContext *abcContextStiffAux = new BiLinFormContext(stiffIntTDEFABCPDPHI1, targetMatrix);

      abcContextStiffAux->SetEntities(actSDList, actSDList);
      abcContextStiffAux->SetFeFunctions(feFunctions_[formulation_], feFunctions_[(SolutionType)(ACOU_TDEF_PHI_V_1 + iPole)]);
      feFunctions_[formulation_]->AddEntityList(actSDList);
      assemble_->AddBiLinearForm(abcContextStiffAux);
    }
    if (TDEFpoleNumber_.complC[iRegion] > 0) {
      EXCEPTION("ABC for TDEF with complex-conjugated poles not implemented yet.");
    }
  }

  void AcousticPDE::CalcCoefsTDEFABC(UInt iRegion, PtrCoefFct &factor, PtrCoefFct &c0, PtrCoefFct &dens, PtrCoefFct &omegaTrg, PtrCoefFct &coeffDamp, PtrCoefFct &coeffStiffTDEF)
  {
    PtrCoefFct abcCoefTDEFre;
    PtrCoefFct abcCoefTDEFim;
    PtrCoefFct abcCoefTDEFreSQ;
    PtrCoefFct abcCoefTDEFimSQ;
    PtrCoefFct abcCoefTDEFdenom;
    abcCoefTDEFre = CoefFunction::Generate(mp_, Global::REAL, CoefXprUnaryOp(mp_, c0, CoefXpr::OP_RE));
    RegionIdType aRegion = regions_[iRegion];

    if (TDEFpoleNumber_.realC[iRegion] == 0 && TDEFpoleNumber_.complC[iRegion] == 0 && TDEFpoleNumber_.realV[iRegion] == 0 && TDEFpoleNumber_.complV[iRegion] == 0) {
      // in the case that no pole are provided, c0 is real-valued
      abcCoefTDEFim = CoefFunction::Generate(mp_, Global::REAL, "0.0");
    }
    else {
      abcCoefTDEFim = CoefFunction::Generate(mp_, Global::REAL, CoefXprUnaryOp(mp_, c0, CoefXpr::OP_IM));
    }
    abcCoefTDEFreSQ = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, abcCoefTDEFre, abcCoefTDEFre, CoefXpr::OP_MULT));
    abcCoefTDEFimSQ = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, abcCoefTDEFim, abcCoefTDEFim, CoefXpr::OP_MULT));
    abcCoefTDEFdenom = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, abcCoefTDEFreSQ, abcCoefTDEFimSQ, CoefXpr::OP_ADD));
    PtrCoefFct coeffDampTmp = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, abcCoefTDEFre, abcCoefTDEFdenom, CoefXpr::OP_DIV));

    if (isTDEFReg_[iRegion]) {
      // scale by high frequency limit of inverse density
      PtrCoefFct invDensInf;
      invDensInf = materials_[aRegion]->GetScalCoefFnc(ACOU_TDEF_INVDENS_CONST, Global::REAL);
      coeffDamp = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, coeffDampTmp, invDensInf, CoefXpr::OP_MULT));
      PtrCoefFct alpha = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, coeffDampTmp, omegaTrg, CoefXpr::OP_MULT));

      // In a TDEF region, an additional stiffness term is required for the ABC
      PtrCoefFct betaTmp = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, abcCoefTDEFim, abcCoefTDEFdenom, CoefXpr::OP_DIV));
      PtrCoefFct beta = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, betaTmp, omegaTrg, CoefXpr::OP_MULT));
      coeffStiffTDEF = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, beta, invDensInf, CoefXpr::OP_MULT));
    }
    else {
      // factor for damping matrix: 1 / density for non-dissipative fluid with real-valued material parameters (e.g. air)
      coeffDampTmp = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, factor, dens, CoefXpr::OP_DIV));
      coeffDamp = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, coeffDampTmp, c0, CoefXpr::OP_DIV));
    }
  }

  void AcousticPDE::DefineRhsLoadIntegrators() {
    LOG_TRACE(acousticpde) << "Defining rhs load integrators for acoustic PDE";
    
    // Get FESpace and FeFunction
    shared_ptr<BaseFeFunction> myFct = feFunctions_[formulation_];
    shared_ptr<FeSpace> mySpace = myFct->GetFeSpace();
    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
    StdVector<shared_ptr<EntityList> > ent;
    StdVector<PtrCoefFct > coef;
    StdVector<std::string> empty;
    LinearForm * lin = NULL;

    // obtain density
    shared_ptr<CoefFunctionMulti> densFct  = matCoefs_[ELEM_DENSITY];
    shared_ptr<CoefFunctionSurf> surfDens(new CoefFunctionSurf(false));
    surfDens->SetVolumeCoefs( densFct->GetRegionCoefs() );
    

    // In the case of acou-mech coupling we have to multiply the
    // integrators by -density
    Double scalFactor = 1.0;
    std::set<RegionIdType> volRegions (regions_.Begin(), regions_.End() );

    bool coefUpdateGeo;
    // ===========================
    //  NORMAL VELOCITY (surface)
    // ===========================
    LOG_DBG(acousticpde) << "Reading normal velocity"; 
    
    ReadRhsExcitation( "normalVelocity", empty, ResultInfo::SCALAR, isComplex_,
                          ent, coef, coefUpdateGeo );

    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      if ( sosAtLaplace_)
        EXCEPTION("NormalVelocity-BC and speed of sound at Laplace operator currently not possible");

      // ensure that list contains only surface elements
      EntityIterator it = ent[i]->GetIterator();
      UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
      if( elemDim != (dim_-1) ) {
        EXCEPTION("Normal velocity can only be defined on surface elements");
      }
      PtrCoefFct exValue;
      CalcMechAcouFacWithCoef(exValue,coef[i],surfDens,scalFactor,part);
      
      
      if( formulation_ == ACOU_POTENTIAL ) {
        if( dim_ == 2) {
          if(isComplex_) {
            lin = new BUIntegrator<Complex,true>( new IdentityOperator<FeH1,2,1>(),
                                                  scalFactor, exValue, volRegions, coefUpdateGeo);
          } else {
            lin = new BUIntegrator<Double,true>( new IdentityOperator<FeH1,2,1>(),
                                                 scalFactor, exValue, volRegions, coefUpdateGeo);
          }
        } else  {
          if(isComplex_) {
            lin = new BUIntegrator<Complex,true>( new IdentityOperator<FeH1,3,1>(),
                                                  1.0, exValue, volRegions, coefUpdateGeo);
          } else {
            lin = new BUIntegrator<Double,true>( new IdentityOperator<FeH1,3,1>(),
                                                 1.0, exValue , volRegions, coefUpdateGeo);
          }
        }
      } else if( formulation_ == ACOU_PRESSURE && isComplex_) {

        // in this case the pressure can be related to the 
        // normal velocity as 
        // p_n = - j*Omega*v_n*rho
        PtrCoefFct tmp = CoefFunction::Generate( mp_, part,"0.0",
                                                        "-2*pi*f");
        if ( complexFluidFormulation_ ) {
          //do not need multiplication with density
          exValue = CoefFunction::Generate( mp_, part,
                                            CoefXprBinOp(mp_, tmp, exValue, CoefXpr::OP_MULT) );
        }
        else {
          PtrCoefFct tmp2 = CoefFunction::Generate( mp_, part,
                                                    CoefXprBinOp(mp_, tmp, exValue, CoefXpr::OP_MULT) );
          exValue = CoefFunction::Generate( mp_, part, CoefXprBinOp(mp_, tmp2, surfDens,
                                            CoefXpr::OP_MULT) );
        }

          if( dim_ == 2) {
            lin = new BUIntegrator<Complex,true>( new IdentityOperator<FeH1,2,1>(),
                                                  1.0, exValue, volRegions, coefUpdateGeo);
          } else  {
            lin = new BUIntegrator<Complex,true>( new IdentityOperator<FeH1,3,1>(),
                                                  1.0, exValue, volRegions, coefUpdateGeo);
          }

      } else {
        EXCEPTION( "Normal velocity can only be prescribed for potential "
            << "formulation or for harmonic pressure formulation" )
      }
      lin->SetName("NormalVelocityIntegrator");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
    }

    
    // ===========================
    //  TOTAL VELOCITY (surface)
    // ===========================
    LOG_DBG(acousticpde) << "Reading total velocity";
    StdVector<std::string> vecDofNames;
    if(dim_ == 3)
      vecDofNames = "x", "y", "z";
    if(dim_ == 2 && !isaxi_)
      vecDofNames = "x", "y";
    if(dim_ == 2 && isaxi_)
      vecDofNames = "r", "z";
    ReadRhsExcitation( "velocity", vecDofNames, ResultInfo::VECTOR, isComplex_,
                       ent, coef,coefUpdateGeo );
    
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      if ( sosAtLaplace_)
        EXCEPTION("RHS and speed of sound at Laplace operator currently not possible");

      // ensure that list contains only surface elements
      EntityIterator it = ent[i]->GetIterator();
      UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
      if( elemDim != (dim_-1) ) {
        EXCEPTION("Normal velocity can only be defined on surface elements");
      }

      PtrCoefFct exValue;
      CalcMechAcouFacWithCoef( exValue,coef[i],surfDens,scalFactor,part);
  
      if( formulation_ == ACOU_POTENTIAL ) {
        if( dim_ == 2) {
          if(isComplex_) {
            lin = new BUIntegrator<Complex,true>( new IdentityOperatorNormal<FeH1,2>(),
                                                  scalFactor, exValue, volRegions, coefUpdateGeo);
          } else {
            lin = new BUIntegrator<Double,true>( new IdentityOperatorNormal<FeH1,2>(),
                                                 scalFactor, exValue, volRegions, coefUpdateGeo);
          }
        } else  {
          if(isComplex_) {
            lin = new BUIntegrator<Complex,true>( new IdentityOperatorNormal<FeH1,3>(),
                                                  scalFactor, exValue, volRegions, coefUpdateGeo);
          } else {
            lin = new BUIntegrator<Double,true>( new IdentityOperatorNormal<FeH1,3>(),
                                                 scalFactor, exValue, volRegions, coefUpdateGeo);
          }
        }
      } else if( formulation_ == ACOU_PRESSURE && isComplex_) {

        // in this case the pressure can be related to the 
        // normal velocity as 
        // p_n = - j*Omega*v_n*rho
        PtrCoefFct tmp = CoefFunction::Generate( mp_, part,"0.0",
                                                "-2*pi*f");
        if ( complexFluidFormulation_ ) {
        	//do not need multiplication with density
        	exValue = CoefFunction::Generate( mp_, part,
        			CoefXprBinOp(mp_, tmp, exValue, CoefXpr::OP_MULT) );
        }
        else {
          PtrCoefFct tmp2 = CoefFunction::Generate( mp_, part, CoefXprBinOp(mp_, tmp, exValue, CoefXpr::OP_MULT) );
          exValue = CoefFunction::Generate( mp_, part, CoefXprBinOp(mp_, surfDens, tmp2, CoefXpr::OP_MULT) );
        }

        if( dim_ == 2) {
          lin = new BUIntegrator<Complex,true>( new IdentityOperatorNormal<FeH1,2>(),
                                                scalFactor, exValue, volRegions, coefUpdateGeo);
        } else  {
          lin = new BUIntegrator<Complex,true>( new IdentityOperatorNormal<FeH1,3>(),
                                                scalFactor, exValue, volRegions, coefUpdateGeo);
        }
      } else {
        EXCEPTION( "Normal velocity can only be prescribed for potential "
            << "formulation or for harmonic pressure formulation" )
      }
      lin->SetName("VelocityIntegrator");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
    }
    

    // ===========================
    //  NORMAL ACCELERATION (surface)
    // ===========================
    LOG_DBG(acousticpde) << "Reading normal acceleration";

    ReadRhsExcitation( "normalAcceleration", empty, ResultInfo::SCALAR, isComplex_,
                          ent, coef, coefUpdateGeo );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      if ( sosAtLaplace_)
        EXCEPTION("RHS and speed of sound at Laplace operator currently not possible");

      // ensure that list contains only surface elements
      EntityIterator it = ent[i]->GetIterator();
      UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
      if( elemDim != (dim_-1) ) {
        EXCEPTION("Normal acceleration can only be defined on surface elements");
      }
      scalFactor = 1.0;
      PtrCoefFct exValue;
      if (isMechCoupled_ == true && formulation_ !=  ACOU_PRESSURE) {
        EXCEPTION("Normal acceleration can only be prescribed for pressure formulation");
      } else {
        exValue = coef[i];
      }

      if( formulation_ == ACOU_PRESSURE ) {
    	  exValue =
    			  CoefFunction::Generate( mp_, part,
    					  CoefXprBinOp(mp_, coef[i],surfDens, CoefXpr::OP_MULT) );

        if( dim_ == 2) {
          if(isComplex_) {
            lin = new BUIntegrator<Complex,true>( new IdentityOperator<FeH1,2,1>(),
            		scalFactor, exValue, volRegions, coefUpdateGeo);
          } else {
            lin = new BUIntegrator<Double,true>( new IdentityOperator<FeH1,2,1>(),
            		scalFactor, exValue, volRegions, coefUpdateGeo);
          }
        } else  {
          if(isComplex_) {
            lin = new BUIntegrator<Complex,true>( new IdentityOperator<FeH1,3,1>(),
            			scalFactor, exValue, volRegions, coefUpdateGeo);
          } else {
            lin = new BUIntegrator<Double,true>( new IdentityOperator<FeH1,3,1>(),
            		scalFactor, exValue , volRegions, coefUpdateGeo);
          }
        }
      } else if( formulation_ == ACOU_POTENTIAL) {
        EXCEPTION( "Normal acceleration can only be prescribed for pressure "
            << "formulation" )
      }
      lin->SetName("NormalAccelerationIntegrator");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
    }


    // ===========================
    //  PRESSURE (surface) 
    // ===========================
    // Only possible for harmonic potential formulation
    ReadRhsExcitation( "pressure", empty, ResultInfo::SCALAR, isComplex_,
                       ent, coef, coefUpdateGeo );
    if( formulation_ == ACOU_POTENTIAL && isComplex_ ) { 
      for( UInt i = 0; i < ent.GetSize(); ++i ) {
        if ( sosAtLaplace_)
          EXCEPTION("RHS and speed of sound at Laplace operator currently not possible");

        // ensure that list contains only surface elements
        EntityIterator it = ent[i]->GetIterator();
        UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
        if( elemDim != (dim_-1) ) {
          EXCEPTION("Pressure can only be defined on surface elements");
        }
        PtrCoefFct exValue;
        CalcMechAcouFacWithCoef( exValue,coef[i],surfDens,scalFactor,part);

        // psi_n = j * 1 / (omega*rho) * p
        PtrCoefFct tmp = CoefFunction::Generate( mp_, Global::COMPLEX, "0","1/(2*pi*f)");
        PtrCoefFct tmp2 = CoefFunction::Generate( mp_, Global::COMPLEX,
                                                  CoefXprBinOp(mp_, tmp, surfDens, 
                                                  CoefXpr::OP_DIV ) );
        exValue = CoefFunction::Generate( mp_, part, 
                                          CoefXprBinOp(mp_, exValue, tmp2, 
                                          CoefXpr::OP_MULT) );
        shared_ptr<InhomDirichletBc> actBc ( new InhomDirichletBc );
        actBc->entities = ent[i];
        actBc->result = myFct->GetResultInfo();
        actBc->value = exValue;
        actBc->dofs.insert(0);
        actBc->updatedGeo = coefUpdateGeo;
        myFct->AddInhomDirichletBc(actBc);
      } // loop entities
    } // if

    // =====================
    //  ACOUSTIC SOURCE DENSITY
    // =====================
    //LOG_DBG(heatcondpde) << "Reading acoustic source density";

    ReadRhsExcitation( "acouSourceDensity", empty,
                       ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("Acoustic source density must be defined on elements")
      }
      EntityIterator it = ent[i]->GetIterator();
      it.Begin();

      if(isComplex_) {
        lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1>(),
                                          Complex(1.0), coef[i], coefUpdateGeo);
      } else  {
        lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1>(),
                                         1.0, coef[i], coefUpdateGeo);
      }
      lin->SetName("AcousticSourceDensityInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);

      RegionIdType regionId = ent[i]->GetRegion();
      acousticSourceDensityCoef_->AddRegion(  regionId, coef[i] );
    } // for

    // =====================================
    //  rhsValues for e.g. for aeroacoustics
    // =====================================
    ReadRhsExcitation( "rhsValues", empty, ResultInfo::SCALAR, isComplex_,
                          ent, coef, coefUpdateGeo );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      if ( sosAtLaplace_)
        EXCEPTION("RHS and speed of sound at Laplace operator currently not possible");

      coef[i]->SetConservative(true);
      //std::cout << "ADD" << std::endl;
      coef[i]->SetFeFunction( feFunctions_[formulation_], formulation_);

      if ( !approxSourceWithDeltaFnc_ && coef[i]->GetInverseType() == CoefFunction::INVSOURCE )  {
    	  std::cout << "Add BUIntegrator SRC" << std::endl;
    	  // check type of entitylist
          if (ent[i]->GetType() == EntityList::NODE_LIST) {
            EXCEPTION("Acoustic source density must be defined on elements")
          }
          EntityIterator it = ent[i]->GetIterator();
          it.Begin();

          if(isComplex_) {
            lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1>(),
                                              Complex(1.0), coef[i], coefUpdateGeo);
          } else  {
            lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1>(),
                                             1.0, coef[i], coefUpdateGeo);
          }
          lin->SetName("AcousticSourceDensityInverseInt");
          LinearFormContext *ctx = new LinearFormContext( lin );
          ctx->SetEntities( ent[i] );
          ctx->SetFeFunction(myFct);
          assemble_->AddLinearForm(ctx);

          RegionIdType regionId = ent[i]->GetRegion();
          acousticSourceDensityCoef_->AddRegion(  regionId, coef[i] );
      }

      //if ( coef[i]->GetInverseType() == CoefFunction::INVSOURCE )
      this->rhsFeFunctions_[formulation_]->AddLoadCoefFunction(coef[i], ent[i]);
    }
    
       // ================
       //  RHS DENSITY Vector
       //	For reading vectoriel source terms of the lightill analogy (div(T)) 
       //      Diss Freidhager 2022
       // ================
       LOG_DBG(acousticpde) << "Reading rhsDensityVector";
       ReadRhsExcitation( "rhsDensityVector", vecDofNames, ResultInfo::VECTOR,
    		   isComplex_, ent, coef, updatedGeo_ );

       //perform checks



       for( UInt i = 0; i < ent.GetSize(); ++i ) {
    	   	   	   	   	   	   	   	   std::cout << "ausgabe     " << i << std::endl;
         // check type of entitylist
         //if (ent[i]->GetType() == EntityList::NODE_LIST) {
         //  EXCEPTION("Rhs density Vector must be defined on elements, because BUIntegrator requires a density")
         //}
         EntityIterator it = ent[i]->GetIterator();
         it.Begin();

         if( formulation_ == ACOU_POTENTIAL){
      	   EXCEPTION("Using rhsDensityVector based on Lighthill's equation -> AcouPressure Formulation!")
         } //closes if AcouPotential
         else{
        	 if(isComplex_) {
           		 std::cout << "\t \t \t Using complexe vector formulation of Lighthill div(T) in grad(phi)" << std::endl;
                 scalFactor = 1.0;
        		 if( dim_ == 2) {

        		 	 lin = new BUIntegrator<Complex>( new GradientOperator<FeH1,2,1,Complex()>(),Complex(scalFactor), coef[i], updatedGeo_);
        		 }else{

            		 lin = new BUIntegrator<Complex>( new GradientOperator<FeH1,3,1,Complex()>(),Complex(scalFactor), coef[i], updatedGeo_);
        		 }

				} else  { //Transiente Simulation

				std::cout << "\t \t \t Using vector formulation of Lighthill div(T) in grad(phi)" << std::endl;
				scalFactor = 1.0;
				if( dim_ == 2) {
					lin = new BUIntegrator<Double,false>( new GradientOperator<FeH1,2,1>(),scalFactor, coef[i], updatedGeo_);
				}else{
					lin = new BUIntegrator<Double,false>( new GradientOperator<FeH1,3,1>(),scalFactor, coef[i], updatedGeo_);
				}
				}
         }

         lin->SetName("RhsDensityVectorInt");
         LinearFormContext *ctx = new LinearFormContext( lin );
         ctx->SetEntities( ent[i] );
         ctx->SetFeFunction(myFct);
         assemble_->AddLinearForm(ctx);
     }


       // ================
       //  BC Vector
       //  Requires a vector v for which v in n is computed on BC surfaces.
       //  Therefore the vector v has to be interpolated onto BC surfaces first
       //  v has to be on elements, since it is a numeric density
       //  For Lighthill thsi means computing div(T)+grad(c_0^2 rho) as BC Vector
       //  Diss Freidhager 2022
       // ================
       LOG_DBG(acousticpde) << "Reading BCVector";
             ReadRhsExcitation( "BCVector", vecDofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo );

             for( UInt i = 0; i < ent.GetSize(); ++i ) {
               EntityIterator it = ent[i]->GetIterator();
               it.Begin();
               // check type of entitylist
		// if (ent[i]->GetType() == EntityList::NODE_LIST) {
		// EXCEPTION("BCVector must be defined on elements, because BUIntegrator requires a density")
		// }

               // ensure that list contains only surface elements
               UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
               if( elemDim != (dim_-1) ) {
                 EXCEPTION("BCVector can only be defined on surface elements");
               }


               if ( sosAtLaplace_) {
            	   EXCEPTION("BCVector not defined for sosAtLaplace");
             	}//closes sosAtLaplace
              else{
            	  if(isComplex_) {
            		  std::cout << "\t \t \t Surface Term of vector formulation of Lighthill's analogy" << std::endl;
            		  scalFactor = 1.0;
		          if( dim_ == 2) {
			      lin = new BUIntegrator<Complex,true>( new IdentityOperatorNormal<FeH1,2>(),
								Complex(scalFactor), coef[i],volRegions, coefUpdateGeo); //IdentityOperatorNormal

			   } else  {
			      lin = new BUIntegrator<Complex,true>( new IdentityOperatorNormal<FeH1,3>(),
								Complex(scalFactor), coef[i],volRegions, coefUpdateGeo);
			   }
            	  }
            	  else {

            	          std::cout << "\t \t \t Surface Term of vector formulation of Lighthill's analogy" << std::endl;
            	          scalFactor = 1.0;
                         if( dim_ == 2) {
            	         	lin = new BUIntegrator<Double,true>( new IdentityOperatorNormal<FeH1,2>(),
            	         			scalFactor, coef[i],volRegions, coefUpdateGeo);
                         } else  {
                              lin = new BUIntegrator<Double,true>( new IdentityOperatorNormal<FeH1,3>(),
            	         			scalFactor, coef[i],volRegions, coefUpdateGeo);
            	         }
            	  }
               }

               lin->SetName("BCVector");
               LinearFormContext *ctx = new LinearFormContext( lin );
               ctx->SetEntities( ent[i] );
               ctx->SetFeFunction(myFct);
               assemble_->AddLinearForm(ctx);
               myFct->AddEntityList(ent[i]);

//               RegionIdType regionId = ent[i]->GetRegion();
//               BCSurfaceVectorCoef_->AddRegion(  regionId, coef[i] );

             }
                      

    // ================
    //  RHS DENSITY
    // ================
    ReadRhsExcitation( "rhsDensity", empty,
                       ResultInfo::SCALAR, isComplex_, ent, coef, updatedGeo_ );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("Rhs density must be defined on elements")
      }
      if ( sosAtLaplace_) {
    	//get region id
    	std::string volRegName = ent[i]->GetName();
    	RegionIdType aRegion = ptGrid_->GetRegion().Parse(volRegName);

    	//get adiabatic exponent
    	PtrCoefFct density = materials_[aRegion]->GetScalCoefFnc( DENSITY, Global::REAL );

    	//get speed of sound and square it
    	std::map<RegionIdType,PtrCoefFct > c0Fcts = matCoefs_[ACOU_ELEM_SPEED_OF_SOUND]->GetRegionCoefs();
    	PtrCoefFct regionC0 = c0Fcts[aRegion];
    	PtrCoefFct regionC02 = CoefFunction::Generate( mp_, Global::REAL,
    			CoefXprBinOp(mp_, regionC0, regionC0, CoefXpr::OP_MULT ) );

    	PtrCoefFct tmpCoef = CoefFunction::Generate( mp_, Global::REAL,
              		  CoefXprBinOp(mp_, regionC02, density, CoefXpr::OP_MULT ) );

    	if(isComplex_) {
    	  PtrCoefFct coefRHS = CoefFunction::Generate( mp_, Global::COMPLEX,
    	      	    		  CoefXprBinOp(mp_, tmpCoef, coef[i], CoefXpr::OP_MULT ) );
    	  lin = new BUIntegrator<Complex>( new IdentityOperator<FeH1>(),
    				  Complex(scalFactor), coefRHS, updatedGeo_);
    	}
    	else {
    	  PtrCoefFct coefRHS = CoefFunction::Generate( mp_, Global::REAL,
    		         		  CoefXprBinOp(mp_, tmpCoef, coef[i], CoefXpr::OP_MULT ) );
    	  lin = new BUIntegrator<Double>( new IdentityOperator<FeH1>(),
    			                          scalFactor, coefRHS, updatedGeo_);
    	}
      }
      else {
        if(isComplex_) {
        	  lin = new BUIntegrator<Complex>( new IdentityOperator<FeH1>(),
        				  Complex(scalFactor), coef[i], updatedGeo_);
        }
        else {
        	lin = new BUIntegrator<Double>( new IdentityOperator<FeH1>(),
        			                          scalFactor, coef[i], updatedGeo_);
        }
      }
      lin->SetName("RhsDensityInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
    } // for


    // ================
    //  RHS Mass time derivative
    // ================
    ReadRhsExcitation( "rhsMassTimeDeriv", empty,
                       ResultInfo::SCALAR, isComplex_, ent, coef, updatedGeo_ );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("rhsMassTimeDeriv must be defined on elements");
      }
      if ( !sosAtLaplace_ )
    	  EXCEPTION("rhsMassTimeDeriv currently just available for speed of sound at Laplace operator");

      //source term is partial time derivative of speed of sound squared
      //multiplied by mass, which has been evaluated by "ReadRhsExcitation"
      //and stored in coef[i]
//      if ( this->analysistype_ == HARMONIC )
//    	  EXCEPTION("rhsMassTimeDeriv in harmonic case for speed of sound at Laplace operator\
//    		         currently not available");

      //get region id
      std::string volRegName = ent[i]->GetName();
      RegionIdType aRegion = ptGrid_->GetRegion().Parse(volRegName);

      //get adiabatic exponent
      PtrCoefFct aExp = materials_[aRegion]->GetScalCoefFnc( FLUID_ADIABATIC_EXPONENT, Global::REAL );

      //get speed of sound and square it
      std::map<RegionIdType,PtrCoefFct > c0Fcts = matCoefs_[ACOU_ELEM_SPEED_OF_SOUND]->GetRegionCoefs();
      PtrCoefFct regionC0 = c0Fcts[aRegion];
      PtrCoefFct regionC02 = CoefFunction::Generate( mp_, Global::REAL,
    		  CoefXprBinOp(mp_, regionC0, regionC0, CoefXpr::OP_MULT ) );

      PtrCoefFct tmpCoef = CoefFunction::Generate( mp_, Global::REAL,
    		  CoefXprBinOp(mp_, regionC02, aExp, CoefXpr::OP_DIV ) );

      if (isComplex_) {
    	  PtrCoefFct coefRHS = CoefFunction::Generate( mp_, Global::COMPLEX,
    		  CoefXprBinOp(mp_, tmpCoef, coef[i], CoefXpr::OP_MULT ) );

    	  lin = new BUIntegrator<Complex>( new IdentityOperator<FeH1>(),
        		                          scalFactor, coefRHS, updatedGeo_);
      }
      else {
          PtrCoefFct coefRHS = CoefFunction::Generate( mp_, Global::REAL,
        		  CoefXprBinOp(mp_, tmpCoef, coef[i], CoefXpr::OP_MULT ) );

          lin = new BUIntegrator<Double>( new IdentityOperator<FeH1>(),
            		                          scalFactor, coefRHS, updatedGeo_);
      }
      lin->SetName("RhsMassTimeDerivInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
    } // for


    // ================
    //  RHS Mass convective term
    // ================
    ReadRhsExcitation( "rhsMassConvective", empty,
                       ResultInfo::SCALAR, isComplex_, ent, coef, updatedGeo_ );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("rhsMassConvective must be defined on elements");
      }

      //source term is convective derivative of speed of sound squared
      //multiplied by mass, which has been evaluated by "ReadRhsExcitation"
      //and stored in coef[i]

      //get region id
      std::string volRegName = ent[i]->GetName();
      RegionIdType aRegion = ptGrid_->GetRegion().Parse(volRegName);

      //get adiabatic exponent
      PtrCoefFct aExp = materials_[aRegion]->GetScalCoefFnc( FLUID_ADIABATIC_EXPONENT, Global::REAL );

      //get speed of sound and square it
      std::map<RegionIdType,PtrCoefFct > c0Fcts = matCoefs_[ACOU_ELEM_SPEED_OF_SOUND]->GetRegionCoefs();
      PtrCoefFct regionC0 = c0Fcts[aRegion];
      PtrCoefFct regionC02 = CoefFunction::Generate( mp_, Global::REAL,
    		  CoefXprBinOp(mp_, regionC0, regionC0, CoefXpr::OP_MULT ) );

      PtrCoefFct tmpCoef = CoefFunction::Generate( mp_, Global::REAL,
          		  CoefXprBinOp(mp_, regionC02, aExp, CoefXpr::OP_DIV ) );

      BaseBOperator* bOp;
      if ( dim_ == 2)
    	  bOp = new ConvectiveOperator<FeH1,2,1>();
      else
    	  bOp = new ConvectiveOperator<FeH1,3,1>();
      bOp->SetCoefFunction( meanFlowCoef_ ); //meanFlow )

      if (isComplex_) {
    	  PtrCoefFct coefRHS = CoefFunction::Generate( mp_, Global::COMPLEX,
    	    		  CoefXprBinOp(mp_, tmpCoef, coef[i], CoefXpr::OP_MULT ) );
    	  lin = new BUIntegrator<Complex>( bOp, scalFactor, coefRHS, updatedGeo_);
      }
      else {
    	  PtrCoefFct coefRHS = CoefFunction::Generate( mp_, Global::REAL,
    	    		  CoefXprBinOp(mp_, tmpCoef, coef[i], CoefXpr::OP_MULT ) );
    	  lin = new BUIntegrator<Double>( bOp, scalFactor, coefRHS, updatedGeo_);
      }

      lin->SetName("RhsMassConvectiveInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
    } // for

  }

  void AcousticPDE::DefineSolveStep(){
    solveStep_ = new StdSolveStep(*this);
  }

  void AcousticPDE::DefinePrimaryResults(){

    // === Primary result according to definition ===
    shared_ptr<ResultInfo> res1( new ResultInfo);
    if ( formulation_ ==  ACOU_PRESSURE) {
      res1->resultType = ACOU_PRESSURE;
      res1->dofNames = "";
      res1->unit = MapSolTypeToUnit(ACOU_PRESSURE);
    } else {
      res1->resultType = ACOU_POTENTIAL;
      res1->dofNames = "";
      res1->unit = MapSolTypeToUnit(ACOU_POTENTIAL);
    }
    res1->definedOn = ResultInfo::NODE;
    res1->entryType = ResultInfo::SCALAR;
    feFunctions_[formulation_]->SetResultInfo(res1);
    results_.Push_back( res1 );
    res1->SetFeFunction(feFunctions_[formulation_]);
    DefineFieldResult( feFunctions_[formulation_], res1 );
    

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    if ( formulation_ ==  ACOU_PRESSURE) {
      hdbcSolNameMap_[ACOU_PRESSURE] = "soundSoft";
      idbcSolNameMap_[ACOU_PRESSURE] = "pressure";
    } else {
      hdbcSolNameMap_[ACOU_POTENTIAL] = "soundSoft";
      idbcSolNameMap_[ACOU_POTENTIAL] = "potential";
    }
    
    // === ACOUSTIC RHS ===
    shared_ptr<ResultInfo> rhs ( new ResultInfo );
    rhs->resultType = ACOU_RHS_LOAD;
    rhs->dofNames = "";
    rhs->unit = "?";
    rhs->definedOn = results_[0]->definedOn;
    rhs->entryType = ResultInfo::SCALAR;
    this->rhsFeFunctions_[formulation_]->SetResultInfo(rhs);
    DefineFieldResult( this->rhsFeFunctions_[formulation_], rhs );
    results_.Push_back( rhs );
    availResults_.insert( rhs );

    //creates the mean flow
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

    //create result info for mean flow velocity
    CreateMeanFlowFunction(vecDofNames);
    
    //create result info for mean temperature field
    //(mean in the sense of time averaged)
    shared_ptr<ResultInfo> temperature( new ResultInfo);
    temperature->resultType = HEAT_MEAN_TEMPERATURE;
    temperature->dofNames = "";
    temperature->unit = MapSolTypeToUnit(HEAT_MEAN_TEMPERATURE);

    temperature->definedOn = ResultInfo::NODE;
    temperature->entryType = ResultInfo::SCALAR;

   	meanTemperatureCoef_.reset(new CoefFunctionMulti(CoefFunction::SCALAR,1,1,false));
    DefineFieldResult( meanTemperatureCoef_, temperature );

    // ==============================================
    //  Define CoefFunctions for material parameters
    // ==============================================
    
    // === DENSITY ===
    // Density \rho_{\mathrm{a}} = {\frac{p_{\mathrm{a}}} {c^{2}}}
    shared_ptr<ResultInfo> density ( new ResultInfo );
    density->resultType = ELEM_DENSITY;
    density->dofNames = "";
    density->unit = MapSolTypeToUnit(ELEM_DENSITY);
    density->definedOn = ResultInfo::ELEMENT;
    density->entryType = ResultInfo::SCALAR;
    shared_ptr<CoefFunctionMulti> densFct(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1,
    		                              complexFluidFormulation_));
    matCoefs_[ELEM_DENSITY] = densFct;
    DefineFieldResult(densFct, density);
    
    // === SPEED OF SOUND ===
    // Speed of sound c = \sqrt{\kappa \frac{p_0} {\rho_0}} = \sqrt{\kappa R T_{0}}
    shared_ptr<ResultInfo> sos ( new ResultInfo );
    sos->resultType = ACOU_ELEM_SPEED_OF_SOUND;
    sos->dofNames = "";
    sos->unit = MapSolTypeToUnit(ACOU_ELEM_SPEED_OF_SOUND);
    sos->definedOn = ResultInfo::ELEMENT;
    sos->entryType = ResultInfo::SCALAR;
    shared_ptr<CoefFunctionMulti> sosFct(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1,
    		                             complexFluidFormulation_));
    matCoefs_[ACOU_ELEM_SPEED_OF_SOUND] = sosFct;
    DefineFieldResult(sosFct, sos);

    // Check if an FE space for Lagrange multiplier has been created.
    // If any, define results for Lagrange multiplier in case of p.b.c.
    PtrParamNode lagSpaceNode = infoNode_->Get("feSpaces", ParamNode::PASS)->Get("lagrangeMultiplier", ParamNode::PASS);
    if (lagSpaceNode)
    {
      // <!> This is a hack, because there is no cross points handling
      hdbcSolNameMap_[LAGRANGE_MULT] = "ground";
      idbcSolNameMap_[LAGRANGE_MULT] = "potential";
      //
      shared_ptr<ResultInfo> lagMultAcou(new ResultInfo);
      lagMultAcou->resultType = LAGRANGE_MULT;
      lagMultAcou->dofNames = "";
      lagMultAcou->unit = "m/s";
      lagMultAcou->entryType = ResultInfo::SCALAR;
      lagMultAcou->SetFeFunction(feFunctions_[LAGRANGE_MULT]);
      lagMultAcou->definedOn = ResultInfo::NODE;
      feFunctions_[LAGRANGE_MULT]->SetResultInfo(lagMultAcou);
      
      results_.Push_back(lagMultAcou);
      DefineFieldResult(feFunctions_[LAGRANGE_MULT], lagMultAcou);
    }
  }
  
  void AcousticPDE::FinalizePostProcResults(){
    //first call base class method
    SinglePDE::FinalizePostProcResults();
  }

  void AcousticPDE::CheckIfIsOnlyOneMaterial(){

    // isOnlyOneMaterial_ initialized to true (assume the best case, if we find one counter example, set it to false)

    // check if there is only one region defined (one region, one material)
    if( regions_.size() <= 1 ) {
      return;
    }

    // get list of parameter nodes for region definitions
    UInt numRegions = 0; // this if fine because it skips everything below if it is not set otherwise
    ParamNodeList regionNodes;
    PtrParamNode regionListNode = domain_->GetParamRoot()->Get("domain")->Get("regionList",ParamNode::PASS );

    if( regionListNode ) {
      regionNodes = regionListNode->GetList("region");
      numRegions = regionNodes.GetSize();
    }

    // get the reference material for the first region
    RegionIdType actRegion = regions_[0];
    std::string regionName = ptGrid_->GetRegion().ToString(actRegion);

    std::string regionNameDomain, refMaterial;
    // iterate over all regions
    for( UInt i = 0; i < numRegions; ++i ) {
      // get data from node
      regionNodes[i]->GetValue("name", regionNameDomain);

      if( regionNameDomain==regionName ) {
        regionNodes[i]->GetValue("material", refMaterial);
      }
    }

    std::string actMaterial;
    // start by 1 and not 0 because 0 is the reference region
    for ( UInt iRegion = 1; iRegion < regions_.GetSize(); iRegion++ ) {
      actRegion = regions_[iRegion];
      regionName = ptGrid_->GetRegion().ToString(actRegion);

      
      for( UInt i = 0; i < numRegions; ++i ) {
        // get data from node
        regionNodes[i]->GetValue("name", regionNameDomain);

        if( regionNameDomain == regionName) {
          regionNodes[i]->GetValue("material", actMaterial);
          
          if ( refMaterial != actMaterial ){
            isOnlyOneMaterial_ = false;
          }
        }
      }
    }
  }

  void AcousticPDE::DefinePostProcResults(){
    
    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
    shared_ptr<BaseFeFunction> feFct = feFunctions_[formulation_];
    shared_ptr<ResultInfo> res1 = feFct->GetResultInfo();
    
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

    StdVector<std::string> tensorDofNames;
    if( ptGrid_->GetDim() == 3 ) {
      tensorDofNames.Resize(9);
      tensorDofNames[0] = "xx";
      tensorDofNames[1] = "xy";
      tensorDofNames[2] = "xz";
      tensorDofNames[3] = "yx";
      tensorDofNames[4] = "yy";
      tensorDofNames[5] = "yz";
      tensorDofNames[6] = "zx";
      tensorDofNames[7] = "zy";
      tensorDofNames[8] = "zz";
    } else {
      if( ptGrid_->IsAxi() ) {
        // Tensor DOF names in are not yet used in axisymmetric problems.
        tensorDofNames.Resize(4);
        tensorDofNames[0] = "rr";
        tensorDofNames[1] = "rz";
        tensorDofNames[2] = "zr";
        tensorDofNames[3] = "zz";
      } else { // 2D case
        tensorDofNames.Resize(4);
        tensorDofNames[0] = "xx";
        tensorDofNames[1] = "xy";
        tensorDofNames[2] = "yx";
        tensorDofNames[3] = "yy";
      }
    }

    // === PRESSURE / POTENTIAL - 1.DERIVATIVE ===
    // first time derivative of potential 
    // \frac{\partial \psi_mathrm{a}} {\partial t} or j \omega \psi_mathrm{a}
    // or of pressure 
    // \frac{\partial p_{\mathrm{a}}} {\partial t} or j \omega p_{\mathrm{a}}
    shared_ptr<ResultInfo> deriv1(new ResultInfo);
    if( formulation_ == ACOU_POTENTIAL ) {
      deriv1->resultType = ACOU_POTENTIAL_DERIV_1;
      deriv1->dofNames = "";
      deriv1->unit = MapSolTypeToUnit(ACOU_POTENTIAL_DERIV_1);
    } else {
      deriv1->resultType = ACOU_PRESSURE_DERIV_1;
      deriv1->dofNames = "";
      deriv1->unit = MapSolTypeToUnit(ACOU_PRESSURE_DERIV_1);
    }
    deriv1->entryType = res1->entryType;
    deriv1->definedOn = res1->definedOn;
    availResults_.insert( deriv1 );
    DefineTimeDerivResult( deriv1->resultType, 1, formulation_ );

    // === PRESSURE / POTENTIAL - 2.DERIVATIVE ===
    // second time derivative of potential 
    // \frac{\partial^{2} \psi_mathrm{a}} {\partial t^{2}} \mathrm or -\omega^{2} \psi_mathrm{a}
    // or of pressure 
    // \frac{\partial^{2} p_{\mathrm{a}}} {\partial t^{2}} or -\omega^{2} p_{\mathrm{a}}
    shared_ptr<ResultInfo> deriv2(new ResultInfo);
    if( formulation_ == ACOU_POTENTIAL ) {
      deriv2->resultType = ACOU_POTENTIAL_DERIV_2;
      deriv2->dofNames = "";
      deriv2->unit = MapSolTypeToUnit(ACOU_POTENTIAL_DERIV_2);
    } else {
      deriv2->resultType = ACOU_PRESSURE_DERIV_2;
      deriv2->dofNames = "";
      deriv2->unit = MapSolTypeToUnit(ACOU_PRESSURE_DERIV_2);
    }
    deriv2->entryType = res1->entryType;
    deriv2->definedOn = res1->definedOn;
    availResults_.insert( deriv2 );
    DefineTimeDerivResult( deriv2->resultType, 2, formulation_ );

    
    // define acoustic pressure, depending on formulation
    PtrCoefFct presFct;
    if( formulation_ == ACOU_POTENTIAL ) {
      shared_ptr<ResultInfo> pres(new ResultInfo);
      pres->resultType = ACOU_PRESSURE;
      pres->dofNames = "";
      pres->unit = MapSolTypeToUnit(ACOU_PRESSURE);
      CheckIfIsOnlyOneMaterial();
      pres->entryType = ResultInfo::SCALAR;
      if(isOnlyOneMaterial_){ // only one material is defined in the whole computational domain
        pres->definedOn = ResultInfo::NODE; //can be defined as nodal quantity - no material jump
      } else{ // if more than one material (i.e. material jump occuring)
        pres->definedOn = ResultInfo::ELEMENT; // result downgraded from nodal to element
        WARN("More than one material is defined in the whole computational domain:\n"
        "-> nodeResult for acouPressure is downgraded to elemResult")
      }
      // Define pressure as p = rho * dPsi/dt
      PtrCoefFct potD1Fct= this->GetCoefFct( ACOU_POTENTIAL_DERIV_1 );
      PtrCoefFct densFct = this->GetCoefFct( ELEM_DENSITY);
      presFct  =
               CoefFunction::Generate( mp_, part,
                                      CoefXprBinOp(mp_, potD1Fct, densFct, CoefXpr::OP_MULT ) );

     //in case of flow, we need the substantial derivative
      // THIS WOULD WORK IF ONLY CefFunctionBOp WOULD BE DEFINED IN A CONSISTENT MANNER!
      // USAGE OF THIS COEFFUNCTION in MechPDE and FluidPDE PROHIBIT THIS AT THE MOMENT!
     /*if( isComplex_ ) {
       convectiveCoef_.reset(new CoefFunctionBOp<Complex>(feFct, res1));
     } else {
       convectiveCoef_.reset(new CoefFunctionBOp<Double>(feFct, res1));
     }


     PtrCoefFct p2Fct  =
              CoefFunction::Generate( mp_, part,
                                     CoefXprBinOp(mp_, densFct,convectiveCoef_, CoefXpr::OP_MULT ) );

      presFct  =
          CoefFunction::Generate( mp_, part, CoefXprBinOp(mp_, p1Fct, p2Fct, CoefXpr::OP_ADD ) );
*/
      DefineFieldResult(presFct, pres);
    } else{
      presFct = this->GetCoefFct(ACOU_PRESSURE);
    }
    
    // === ACOU_SURFACE_PRESSURE ===
    shared_ptr<ResultInfo> surfPres(new ResultInfo);
    surfPres->resultType = ACOU_SURFPRESSURE;
    surfPres->dofNames = "";
    surfPres->unit = MapSolTypeToUnit(ACOU_SURFPRESSURE);
    surfPres->entryType = ResultInfo::SCALAR;
    surfPres->definedOn = ResultInfo::SURF_ELEM;
    DefineFieldResult(presFct, surfPres);

    // optimization results are provided in DesignSpace::ExtractResults()
    // === PSEUDO_DENISTY ===
    shared_ptr<ResultInfo> mpd(new ResultInfo);
    mpd->resultType = PSEUDO_DENSITY;
    mpd->entryType = ResultInfo::SCALAR;
    mpd->definedOn = ResultInfo::ELEMENT;
    mpd->dofNames = MapSolTypeToUnit(PSEUDO_DENSITY);
    mpd->fromOptimization = true;
    DefineFieldResult(shared_ptr<FeFunction<double> >(new FeFunction<double>(NULL)), mpd); // the fe-function is only a dummy

    // === PHYSICAL_PSEUDO_DENISTY ===
    shared_ptr<ResultInfo> ppd(new ResultInfo);
    ppd->resultType = PHYSICAL_PSEUDO_DENSITY;
    ppd->entryType = ResultInfo::SCALAR;
    ppd->definedOn = ResultInfo::ELEMENT;
    ppd->dofNames = MapSolTypeToUnit(PHYSICAL_PSEUDO_DENSITY);
    ppd->fromOptimization = true;
    DefineFieldResult(shared_ptr<FeFunction<double> >(new FeFunction<double>(NULL)), ppd);
  
    // === ACOU_FORCE ===
    // Force F_{\mathrm{a}} = \int_{\Gamma} p_{\mathrm{a}} \ \mathrm{d} \Gamma
    shared_ptr<ResultInfo> force(new ResultInfo);
    force->resultType = ACOU_FORCE;
    force->dofNames = "";
    force->unit = MapSolTypeToUnit(ACOU_FORCE);
    force->entryType = ResultInfo::SCALAR;
    force->definedOn = ResultInfo::SURF_REGION;
    // Force F = \int_Gamma p *n dGamma
    // Integrate pressure over surface
    shared_ptr<ResultFunctor> forceFct;
    if( isComplex_ ) {
      forceFct.reset(new ResultFunctorIntegrate<Complex>(presFct, feFct, force ) );
    } else {
      forceFct.reset(new ResultFunctorIntegrate<Double>(presFct, feFct, force ) );
    }
    resultFunctors_[ACOU_FORCE] = forceFct;
    availResults_.insert(force);


    shared_ptr<CoefFunctionFormBased> presGradFct;
    
    // some results are only available in potential formulation (both transient and harmonic) or
    // pressure formulation (harmonic only)
    if( formulation_ == ACOU_POTENTIAL || isComplex_ ){
      // === ACOU_VELOCITY ===
      // for potential formulation
      // Velocity \bm{v}_{\mathrm{a}} = -\nabla \psi_\mathrm{a}
      // or for pressure formulation
      // Velocity \bm{v}_\mathrm{a} = j \frac{1} {\omega* \rho} \nabla p_\mathrm{a}
      shared_ptr<ResultInfo> vel(new ResultInfo);
      vel->resultType = ACOU_VELOCITY;
      vel->dofNames = vecDofNames;
      vel->unit = MapSolTypeToUnit(ACOU_VELOCITY);
      vel->entryType = ResultInfo::VECTOR;
      vel->definedOn = ResultInfo::ELEMENT;
      PtrCoefFct velFct;
      shared_ptr<CoefFunctionFormBased> velFctPot;
      if( formulation_ == ACOU_POTENTIAL) {
        // Velocity v = - grad Psi
        if( isComplex_ ) {
          velFctPot.reset(new CoefFunctionBOp<Complex>(feFct, vel, -1.0));
        } else {
          velFctPot.reset(new CoefFunctionBOp<Double>(feFct, vel, -1.0));
        }
        velFct = velFctPot;
        stiffFormCoefs_.insert(velFctPot);
        DefineFieldResult( velFct, vel );
      }
      else if( formulation_ == ACOU_PRESSURE && isComplex_ ) {
        // Velocity v = j* 1/(omega*rho) * grad(p)
        // a) define gradient of pressure
        presGradFct.reset(new CoefFunctionBOp<Complex>(feFct, vel, 1.0));  
        stiffFormCoefs_.insert(presGradFct);
        
        // b) multiply by factor
        PtrCoefFct densFct = this->GetCoefFct( ELEM_DENSITY);
        PtrCoefFct factor = 
            CoefFunction::Generate( mp_, Global::COMPLEX, "0","1/(2*pi*f)");
        PtrCoefFct factor2 = 
            CoefFunction::Generate( mp_, Global::COMPLEX,
                                  CoefXprBinOp(mp_,factor, densFct, CoefXpr::OP_DIV ) );
        velFct = 
            CoefFunction::Generate( mp_,  part,
                                    CoefXprBinOp( mp_, factor2, presGradFct, CoefXpr::OP_MULT ) );
        DefineFieldResult( velFct, vel );
      }

      // === ACOU_NORMAL_VELOCITY ===
      // Normal velocity v_{\mathrm{a,n}} = -\nabla \psi_mathrm{a} \cdot \bm{n}
      shared_ptr<ResultInfo> velNormal(new ResultInfo);
      velNormal->resultType = ACOU_NORMAL_VELOCITY;
      velNormal->dofNames = "";
      velNormal->unit = MapSolTypeToUnit(ACOU_NORMAL_VELOCITY);
      velNormal->entryType = ResultInfo::SCALAR;
      velNormal->definedOn = ResultInfo::SURF_ELEM;
      shared_ptr<CoefFunctionSurf> velFctNormal;
      velFctNormal.reset(new CoefFunctionSurf(true, 1.0, velNormal));
      std::string quantity = "acouNormalVelocity";
      this->SetSurfVolNeighborRegion(velFctNormal, quantity);
      DefineFieldResult(velFctNormal, velNormal);
      if( formulation_ == ACOU_POTENTIAL ) {
        surfCoefFcts_[velFctNormal] = velFctPot;
      }
      else if ( formulation_ == ACOU_PRESSURE && isComplex_ ) {
        surfCoefFcts_[velFctNormal] = velFct;
      }
      
      // === ACOU_SURFIMPEDANCE ===
      shared_ptr<ResultInfo> surfImpedance(new ResultInfo);
      surfImpedance->resultType = ACOU_SURFIMPEDANCE;
      surfImpedance->dofNames = "";
      surfImpedance->unit = MapSolTypeToUnit(ACOU_SURFIMPEDANCE);
      surfImpedance->entryType = ResultInfo::SCALAR;
      surfImpedance->definedOn = ResultInfo::SURF_ELEM;
      PtrCoefFct surfImpFct = CoefFunction::Generate( mp_, part, CoefXprBinOp(mp_, presFct, velFctNormal, CoefXpr::OP_DIV ) );
      DefineFieldResult(surfImpFct, surfImpedance);
      
      // === ACOU_SURF_AVG_IMPEDANCE ===
      shared_ptr<ResultInfo> impedance(new ResultInfo);
      impedance->resultType = ACOU_SURF_AVG_IMPEDANCE;
      impedance->dofNames = "";
      impedance->unit = MapSolTypeToUnit(ACOU_SURF_AVG_IMPEDANCE);
      impedance->entryType = ResultInfo::SCALAR;
      impedance->definedOn = ResultInfo::SURF_REGION;
      // area weighted average
      shared_ptr<ResultFunctor> impedanceFct;
      if( isComplex_ ) {
        impedanceFct.reset(new ResultFunctorIntegrate<Complex>(surfImpFct, 
                                                            feFct, impedance ) );
        dynamic_pointer_cast< ResultFunctorIntegrate<Complex> >(impedanceFct)->SetAveraged(true);
      } else {
        impedanceFct.reset(new ResultFunctorIntegrate<Double>(surfImpFct, 
                                                          feFct, impedance ) );
        dynamic_pointer_cast< ResultFunctorIntegrate<Double> >(impedanceFct)->SetAveraged(true);                                                  
      }
      resultFunctors_[ACOU_SURF_AVG_IMPEDANCE] = impedanceFct;
      availResults_.insert(impedance);

      // === ACOU_INTENSITY ===
      // Intensity \bm{\mathrm{I}}_{\mathrm{a}} = p_{\mathrm{a}} \overline{\bm{v}_{\mathrm{a}}}
      shared_ptr<ResultInfo> intensity(new ResultInfo);
      intensity->resultType = ACOU_INTENSITY;
      intensity->dofNames = vecDofNames;
      intensity->unit = MapSolTypeToUnit(ACOU_INTENSITY);
      intensity->entryType = ResultInfo::VECTOR;
      intensity->definedOn = ResultInfo::ELEMENT;
      PtrCoefFct intensFct;
      // Intensity I = p * conj(v)
      intensFct = 
          CoefFunction::Generate( mp_, part,
                                CoefXprBinOp(mp_, presFct, velFct, CoefXpr::OP_MULT_CONJ ) );
      DefineFieldResult(intensFct, intensity);
      
      // === ACOU_SURFINTENSITY ===
      // Surface intensity \bm{I}_{\mathrm{a,surf}} =  p_{\mathrm{a}} \bm{v_{\mathrm{a}}} \cdot \bm{n}
      shared_ptr<ResultInfo> surfIntensity(new ResultInfo);
      surfIntensity->resultType = ACOU_SURFINTENSITY;
      surfIntensity->dofNames = vecDofNames;
      surfIntensity->unit = MapSolTypeToUnit(ACOU_SURFINTENSITY);
      surfIntensity->entryType = ResultInfo::VECTOR;
      surfIntensity->definedOn = ResultInfo::SURF_ELEM;
      shared_ptr<CoefFunctionSurf> sIntens;
      sIntens.reset(new CoefFunctionSurf(false, 1.0, surfIntensity));
      DefineFieldResult(sIntens, surfIntensity);
      surfCoefFcts_[sIntens] = intensFct;

      // === ACOU_NORMAL_INTENSITY ===
      // Normal intensity I_{\mathrm{a}} = p_{\mathrm{a}} \bm{v_{\mathrm{a}}} \cdot \bm{n}
      shared_ptr<ResultInfo> intensNormal(new ResultInfo);
      intensNormal->resultType = ACOU_NORMAL_INTENSITY;
      intensNormal->dofNames = "";
      intensNormal->unit = MapSolTypeToUnit(ACOU_NORMAL_INTENSITY);
      intensNormal->entryType = ResultInfo::SCALAR;
      intensNormal->definedOn = ResultInfo::SURF_ELEM;
      shared_ptr<CoefFunctionSurf> sNormIntens;
      sNormIntens.reset(new CoefFunctionSurf(true, 1.0, intensNormal));
      DefineFieldResult( sNormIntens, intensNormal );
      surfCoefFcts_[sNormIntens] = intensFct;

      // === ACOU_POWER ===
      // Power P_{\mathrm{a}} = \int_{\Gamma} I_{\mathrm{a}} \ \mathrm{d} \Gamma
      shared_ptr<ResultInfo> power(new ResultInfo);
      power->resultType = ACOU_POWER;
      power->dofNames = "";
      power->unit = MapSolTypeToUnit(ACOU_POWER);
      power->entryType = ResultInfo::SCALAR;
      power->definedOn = ResultInfo::SURF_REGION;
      // Power p = \int_Gamma I *n dGamma
      // Integrate normal intensity
      shared_ptr<ResultFunctor> powerFct;
      if( isComplex_ ) {
        powerFct.reset(new ResultFunctorIntegrate<Complex>(sNormIntens, 
                                                            feFct, power ) );
      } else {
        powerFct.reset(new ResultFunctorIntegrate<Double>(sNormIntens, 
                                                          feFct, power ) );
      }
      resultFunctors_[ACOU_POWER] = powerFct;
      availResults_.insert(power);
    }

    // some results are only available in pressure formulation with harmonic simulation
    if ( formulation_ == ACOU_PRESSURE && isComplex_ ) {
      // --------------------------------
      //  COMPLEX & PRESSURE FORMULATION
      // --------------------------------
      
      // === ACOU_POSITION ===
      // Particle Position u_\mathrm{a} = \frac{1} {\rho* \omega^2} \nabla p_\mathrm{a}
      shared_ptr<ResultInfo> pos(new ResultInfo);
      pos->resultType = ACOU_POSITION;
      pos->dofNames = vecDofNames;
      pos->unit = MapSolTypeToUnit(ACOU_POSITION);
      pos->entryType = ResultInfo::VECTOR;
      pos->definedOn = ResultInfo::ELEMENT;
      // u = 1/(rho*omega^2) * grad(p)
      PtrCoefFct one = CoefFunction::Generate( mp_, Global::REAL, "1.0");
      PtrCoefFct densFct = this->GetCoefFct( ELEM_DENSITY);
      PtrCoefFct oneOverOmega2rho = CoefFunction::Generate( mp_, part,
          CoefXprBinOp( mp_, one,
            CoefXprBinOp(mp_,CoefFunction::Generate( mp_, Global::REAL, "4*pi*pi*f*f"), densFct, CoefXpr::OP_MULT ),
          CoefXpr::OP_DIV ));
      PtrCoefFct posFct = CoefFunction::Generate( mp_,  part, CoefXprBinOp( mp_, oneOverOmega2rho, presGradFct, CoefXpr::OP_MULT ) );
      DefineFieldResult( posFct, pos );

      // === ACOU_NORMAL_INTENSITY ===
      // Normal intensity I_{\mathrm{a}} = p_{\mathrm{a}} \bm{v_{\mathrm{a}}} \cdot \bm{n}
      shared_ptr<ResultInfo> intensNormalPW;
      intensNormalPW.reset(new ResultInfo);
      intensNormalPW->resultType = ACOU_NORMAL_INTENSITY_PLANEWAVE;
      intensNormalPW->dofNames = "";
      intensNormalPW->unit = MapSolTypeToUnit(ACOU_NORMAL_INTENSITY_PLANEWAVE);
      intensNormalPW->entryType = ResultInfo::SCALAR;
      intensNormalPW->definedOn = ResultInfo::SURF_ELEM;

      //compute 1/(2*rho*c0)
      //shared_ptr<CoefFunctionMulti> c0Fct  = matCoefs_[ACOU_ELEM_SPEED_OF_SOUND];
      //shared_ptr<CoefFunctionSurf> surfC0(new CoefFunctionSurf(false));
      //surfC0->SetVolumeCoefs( c0Fct->GetRegionCoefs() );

      PtrCoefFct c0Fct = this->GetCoefFct( ACOU_ELEM_SPEED_OF_SOUND);
      PtrCoefFct constVal = CoefFunction::Generate( mp_, Global::REAL, "0.5");
      PtrCoefFct val1 =
                CoefFunction::Generate( mp_, Global::COMPLEX,
                                      CoefXprBinOp(mp_,c0Fct, densFct, CoefXpr::OP_MULT ) );
      PtrCoefFct val2 =
                      CoefFunction::Generate( mp_, Global::COMPLEX,
                                            CoefXprBinOp(mp_,constVal, val1, CoefXpr::OP_DIV ) );

      PtrCoefFct velFctPW =
          CoefFunction::Generate( mp_,  part,
                                  CoefXprBinOp( mp_, val2, presFct, CoefXpr::OP_MULT ) );

      PtrCoefFct intensPWfnc =
                      CoefFunction::Generate( mp_, Global::COMPLEX,
                                            CoefXprBinOp(mp_, velFctPW, presFct, CoefXpr::OP_MULT_CONJ ) );

      // normal acoustic intensity for plane waves
      shared_ptr<CoefFunctionSurf> sNormIntensPW;
      sNormIntensPW.reset(new CoefFunctionSurf(false, 1.0, intensNormalPW));
      DefineFieldResult( sNormIntensPW, intensNormalPW );
      surfCoefFcts_[sNormIntensPW] = intensPWfnc;

      // === ACOU_POWER WITH PLANE WAVE ASSUMPTION: p/vn = \rho c ===
      // Power P = \int_{\Gamma} \frac{p_{\mathrm{a}}^2} {2 \rho c} \ \mathrm{d} \Gamma
      shared_ptr<ResultInfo> powerPW;
      powerPW.reset(new ResultInfo);
      powerPW->resultType = ACOU_POWER_PLANEWAVE;
      powerPW->dofNames = "";
      powerPW->unit = MapSolTypeToUnit(ACOU_POWER_PLANEWAVE);
      powerPW->entryType = ResultInfo::SCALAR;
      powerPW->definedOn = ResultInfo::SURF_REGION;
      shared_ptr<ResultFunctor> powerFctPW;
      powerFctPW.reset(new ResultFunctorIntegrate<Complex>(sNormIntensPW,
                                                        feFct, powerPW ) );
      resultFunctors_[ACOU_POWER_PLANEWAVE] = powerFctPW;
      availResults_.insert(powerPW);
    }
    
    // === ACOUSTIC KINETIC ENERGY ===
    // kinetic Energy w_{\mathrm{a}}^{\mathrm{kin}} = \frac{1} {2}\rho_{0} \bm{v}_{\mathrm{a}} \cdot \bm{v}_{\mathrm{a}}
    shared_ptr<ResultFunctor> keFunc;
    shared_ptr<ResultInfo> kinEnergy(new ResultInfo);
    kinEnergy->resultType = ACOU_KIN_ENERGY;
    kinEnergy->dofNames = "";
    kinEnergy->unit = MapSolTypeToUnit(ACOU_KIN_ENERGY);
    kinEnergy->entryType = ResultInfo::SCALAR;
    kinEnergy->definedOn = ResultInfo::REGION;
    availResults_.insert ( kinEnergy );

    shared_ptr<BaseFeFunction> deriv1vFct;
    if ( formulation_ == ACOU_PRESSURE )
      deriv1vFct = timeDerivFeFunctions_[ACOU_PRESSURE_DERIV_1];
    else
      deriv1vFct = timeDerivFeFunctions_[ACOU_POTENTIAL_DERIV_1];

    if( isComplex_ ) {
      keFunc.reset(new EnergyResultFunctor<Complex>(deriv1vFct, kinEnergy, 0.5));
    } else {
      keFunc.reset(new EnergyResultFunctor<Double>(deriv1vFct, kinEnergy, 0.5));
    }
    resultFunctors_[ACOU_KIN_ENERGY] = keFunc;
    massFormFunctors_.insert(keFunc);

    // === ACOUSTIC POTENTIAL ENERGY ===
    // potential Energy w_{\mathrm{a}}^{\mathrm{pot}} = \frac{p_{\mathrm{a}}^{2}} {2 \rho_{0} c^{2}}
    shared_ptr<ResultFunctor> keFuncPot;
    shared_ptr<ResultInfo> potEnergy(new ResultInfo);
    potEnergy->resultType = ACOU_POT_ENERGY;
    potEnergy->dofNames = "";
    potEnergy->unit = MapSolTypeToUnit(ACOU_POT_ENERGY);
    potEnergy->entryType = ResultInfo::SCALAR;
    potEnergy->definedOn = ResultInfo::REGION;
    availResults_.insert ( potEnergy );
    if( isComplex_ ) {
      keFuncPot.reset(new EnergyResultFunctor<Complex>(feFct, potEnergy, 0.5));
    } else {
      keFuncPot.reset(new EnergyResultFunctor<Double>(feFct, potEnergy, 0.5));
    }
    resultFunctors_[ACOU_POT_ENERGY] = keFuncPot;
    stiffFormFunctors_.insert(keFuncPot);

    //== ACOUSTIC LOAD DDENSITY  ====
    shared_ptr<ResultInfo> loadDensity( new ResultInfo);
    loadDensity->resultType = ACOU_RHS_LOAD_DENSITY;
    loadDensity->dofNames = "";
    loadDensity->unit = "";

    loadDensity->definedOn = ResultInfo::NODE;
    loadDensity->entryType = ResultInfo::SCALAR;

    acousticSourceDensityCoef_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1,isComplex_));
    DefineFieldResult( acousticSourceDensityCoef_,loadDensity );

    // === PML Output Parameters ===
    // Vector holding the eigenvalues of the PML matrix
    // this are the diagonal elements for CLASSIC formulation, 
    // or the diagonal elements of the inner (stretching) matrix 
    // for the CURVILINEAR PML formulation 
    shared_ptr<ResultInfo> pmlDampFactor ( new ResultInfo );
    pmlDampFactor->resultType = PML_DAMP_FACTOR;
    pmlDampFactor->dofNames = vecDofNames;
    pmlDampFactor->unit = "";
    pmlDampFactor->definedOn = ResultInfo::ELEMENT;
    pmlDampFactor->entryType = ResultInfo::VECTOR;
    shared_ptr<CoefFunctionMulti> pmlDampFactorCoefFct(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_));
    matCoefs_[PML_DAMP_FACTOR] = pmlDampFactorCoefFct;
    DefineFieldResult(pmlDampFactorCoefFct, pmlDampFactor);

    // Curvilinear PML: currently only available for harmonic simulation
    if (this->isComplex_) {
      // Tensor holding the complete the PML matrix. 
      // only set for CURVILINEAR PML formulation
      shared_ptr<ResultInfo> pmlTensor ( new ResultInfo );
      pmlTensor->resultType = PML_TENSOR;
      pmlTensor->dofNames = tensorDofNames;
      pmlTensor->unit = "";
      pmlTensor->definedOn = ResultInfo::ELEMENT;
      pmlTensor->entryType = ResultInfo::TENSOR;
      shared_ptr<CoefFunctionMulti> pmlTensorCoefFct(new CoefFunctionMulti(CoefFunction::TENSOR, dim_, dim_, isComplex_));
      matCoefs_[PML_TENSOR] = pmlTensorCoefFct;
      DefineFieldResult(pmlTensorCoefFct, pmlTensor);

      // Scalar holding the determinant of the PML matrix. 
      // only set for CURVILINEAR PML formulation
      shared_ptr<ResultInfo> pmlDeterminant ( new ResultInfo );
      pmlDeterminant->resultType = PML_DETERMINANT;
      pmlDeterminant->dofNames = "";
      pmlDeterminant->unit = "";
      pmlDeterminant->definedOn = ResultInfo::ELEMENT;
      pmlDeterminant->entryType = ResultInfo::SCALAR;
      shared_ptr<CoefFunctionMulti> pmlDetCoefFct(new CoefFunctionMulti(CoefFunction::SCALAR, dim_, dim_, isComplex_));
      matCoefs_[PML_DETERMINANT] = pmlDetCoefFct;
      DefineFieldResult(pmlDetCoefFct, pmlDeterminant);

      // Scalar holding the distance between points and the PML interface. 
      // only set for CURVILINEAR PML formulation
      shared_ptr<ResultInfo> pmlDistance ( new ResultInfo );
      pmlDistance->resultType = PML_DISTANCE;
      pmlDistance->dofNames = "";
      pmlDistance->unit = "";
      pmlDistance->definedOn = ResultInfo::ELEMENT;
      pmlDistance->entryType = ResultInfo::SCALAR;
      shared_ptr<CoefFunctionMulti> pmlDistanceCoefFct(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));
      matCoefs_[PML_DISTANCE] = pmlDistanceCoefFct;
      DefineFieldResult(pmlDistanceCoefFct, pmlDistance);
    }

    // === AUX Variables for transient PML===
    if(this->isTimeDomPML_){
      if(!this->isAPML_ && dim_ == 3){
        shared_ptr<ResultInfo> pmlScal ( new ResultInfo );
        pmlScal->resultType = ACOU_PMLAUXSCALAR;
        pmlScal->dofNames = "";
        pmlScal->unit = "-";
        pmlScal->definedOn = ResultInfo::NODE;
        pmlScal->entryType = ResultInfo::SCALAR;
        feFunctions_[ACOU_PMLAUXSCALAR]->SetResultInfo(pmlScal);
        results_.Push_back( pmlScal );
        pmlScal->SetFeFunction(feFunctions_[ACOU_PMLAUXSCALAR]);
        DefineFieldResult( feFunctions_[ACOU_PMLAUXSCALAR], pmlScal );
      }

      shared_ptr<ResultInfo> pmlVec ( new ResultInfo );
      pmlVec->resultType = ACOU_PMLAUXVEC;
      pmlVec->dofNames = vecDofNames;
      pmlVec->unit = "-";
      pmlVec->definedOn = ResultInfo::NODE;
      pmlVec->entryType = ResultInfo::VECTOR;
      feFunctions_[ACOU_PMLAUXVEC]->SetResultInfo(pmlVec);
      results_.Push_back( pmlVec );
      pmlVec->SetFeFunction(feFunctions_[ACOU_PMLAUXVEC]);
      DefineFieldResult( feFunctions_[ACOU_PMLAUXVEC], pmlVec );
    }

    // === TDEF AUX Variables ===
    if (this->timeDomainEqFluidFormulation_) {
      for (unsigned int i = 0; i < TDEFpoleNumber_.realC.Max(); i++) {
        shared_ptr<ResultInfo> resTDEFPhiC(new ResultInfo);
        resTDEFPhiC->resultType = (SolutionType)(ACOU_TDEF_PHI_C_1 + i);
        resTDEFPhiC->dofNames = "";
        resTDEFPhiC->unit = "-";
        resTDEFPhiC->definedOn = ResultInfo::NODE;
        resTDEFPhiC->entryType = ResultInfo::SCALAR;
        feFunctions_[(SolutionType)(ACOU_TDEF_PHI_C_1 + i)]->SetResultInfo(resTDEFPhiC);
        results_.Push_back(resTDEFPhiC);
        resTDEFPhiC->SetFeFunction(feFunctions_[(SolutionType)(ACOU_TDEF_PHI_C_1 + i)]);
        DefineFieldResult(feFunctions_[(SolutionType)(ACOU_TDEF_PHI_C_1 + i)], resTDEFPhiC);
      }
      for (unsigned int i = 0; i < TDEFpoleNumber_.complC.Max(); i++) {
        shared_ptr<ResultInfo> resTDEFPsiC(new ResultInfo);
        resTDEFPsiC->resultType = (SolutionType)(ACOU_TDEF_PSI_C_1 + i);
        resTDEFPsiC->dofNames = "";
        resTDEFPsiC->unit = "-";
        resTDEFPsiC->definedOn = ResultInfo::NODE;
        resTDEFPsiC->entryType = ResultInfo::SCALAR;
        feFunctions_[(SolutionType)(ACOU_TDEF_PSI_C_1 + i)]->SetResultInfo(resTDEFPsiC);
        results_.Push_back(resTDEFPsiC);
        resTDEFPsiC->SetFeFunction(feFunctions_[(SolutionType)(ACOU_TDEF_PSI_C_1 + i)]);
        DefineFieldResult(feFunctions_[(SolutionType)(ACOU_TDEF_PSI_C_1 + i)], resTDEFPsiC);
      }
      for (unsigned int i = 0; i < TDEFpoleNumber_.realV.Max(); i++) {
        shared_ptr<ResultInfo> resTDEFPhiV(new ResultInfo);
        resTDEFPhiV->resultType = (SolutionType)(ACOU_TDEF_PHI_V_1 + i);
        resTDEFPhiV->dofNames = "";
        resTDEFPhiV->unit = "-";
        resTDEFPhiV->definedOn = ResultInfo::NODE;
        resTDEFPhiV->entryType = ResultInfo::SCALAR;
        feFunctions_[(SolutionType)(ACOU_TDEF_PHI_V_1 + i)]->SetResultInfo(resTDEFPhiV);
        results_.Push_back(resTDEFPhiV);
        resTDEFPhiV->SetFeFunction(feFunctions_[(SolutionType)(ACOU_TDEF_PHI_V_1 + i)]);
        DefineFieldResult(feFunctions_[(SolutionType)(ACOU_TDEF_PHI_V_1 + i)], resTDEFPhiV);
      }
      for (unsigned int i = 0; i < TDEFpoleNumber_.complV.Max(); i++) {
        shared_ptr<ResultInfo> resTDEFPsiV(new ResultInfo);
        resTDEFPsiV->resultType = (SolutionType)(ACOU_TDEF_PSI_V_1 + i);
        resTDEFPsiV->dofNames = "";
        resTDEFPsiV->unit = "-";
        resTDEFPsiV->definedOn = ResultInfo::NODE;
        resTDEFPsiV->entryType = ResultInfo::SCALAR;
        feFunctions_[(SolutionType)(ACOU_TDEF_PSI_V_1 + i)]->SetResultInfo(resTDEFPsiV);
        results_.Push_back(resTDEFPsiV);
        resTDEFPsiV->SetFeFunction(feFunctions_[(SolutionType)(ACOU_TDEF_PSI_V_1 + i)]);
        DefineFieldResult(feFunctions_[(SolutionType)(ACOU_TDEF_PSI_V_1 + i)], resTDEFPsiV);
      }
    }
  }


   void AcousticPDE::CreateMeanFlowFunction(StdVector<std::string> dofNames){
     //// === MEAN FLUIDMECH VELOCITY ===
     shared_ptr<ResultInfo> flowvelocity( new ResultInfo);
     flowvelocity->resultType = MEAN_FLUIDMECH_VELOCITY;
     flowvelocity->dofNames = dofNames;
     flowvelocity->unit = MapSolTypeToUnit(MEAN_FLUIDMECH_VELOCITY);

     flowvelocity->definedOn = ResultInfo::NODE;
     flowvelocity->entryType = ResultInfo::VECTOR;
     
     meanFlowCoef_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_,1,isComplex_));
     DefineFieldResult( meanFlowCoef_, flowvelocity );

     results_.Push_back( flowvelocity );
     availResults_.insert( flowvelocity );
     if(myParam_->Get("flowFormulation")->As<std::string>() == "withDivergence"){
       //// === DIVERGENCE OF MEAN FLOW ===
       shared_ptr<ResultInfo> divflowvelocity( new ResultInfo);
       divflowvelocity->resultType = DIV_MEAN_FLUIDMECH_VELOCITY;
       divflowvelocity->dofNames = "";
       divflowvelocity->unit = MapSolTypeToUnit(DIV_MEAN_FLUIDMECH_VELOCITY);

       divflowvelocity->definedOn = ResultInfo::ELEMENT;
       divflowvelocity->entryType = ResultInfo::SCALAR;

       divMeanFlowCoef_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1,isComplex_));
       DefineFieldResult( divMeanFlowCoef_, divflowvelocity );

       results_.Push_back( divflowvelocity );
       availResults_.insert( divflowvelocity );
     }
   }

  //! Init the time stepping
  void AcousticPDE::InitTimeStepping(){
      
    // Check if time integration is defined in XML input
    PtrParamNode transientNode = myParam_->GetParent()->GetParent()->Get("analysis")->Get("transient", ParamNode::PASS);
    PtrParamNode integrationScheme = transientNode->Get("integrationScheme", ParamNode::PASS);

    PtrParamNode timeStepAlphaNode = this->myParam_->Get("timeStepAlpha", ParamNode::PASS);
    if (integrationScheme && timeStepAlphaNode)
      throw Exception("Both 'integrationScheme' and 'timeStepAlpha' are specified for the acoustic PDE. "
                      "Please use 'integrationScheme' only, as it provides more flexibility and "
                      "supersedes the legacy 'timeStepAlpha' parameter.");

    // Helper lambda to create the appropriate scheme
    auto makeScheme = [&]() -> GLMScheme* {
      if (integrationScheme)
        return GetXmlDefinedScheme(integrationScheme);
      else {
        Double alpha = this->myParam_->Get("timeStepAlpha")->As<Double>();
        return new Newmark(0.5, 0.25, alpha);
      }
    };
  
    // Main formulation scheme
    shared_ptr<BaseTimeScheme> acouScheme(new TimeSchemeGLM(makeScheme(), 0));
    feFunctions_[formulation_]->SetTimeScheme(acouScheme);
  
    // PML auxiliary vector scheme
    if (this->isTimeDomPML_)
    {
      shared_ptr<BaseTimeScheme> vecScheme(new TimeSchemeGLM(makeScheme(), 0));
      feFunctions_[ACOU_PMLAUXVEC]->SetTimeScheme(vecScheme);
    
      if (!this->isAPML_ && dim_ == 3)
      {
        shared_ptr<BaseTimeScheme> scalScheme(new TimeSchemeGLM(makeScheme(), 0));
        feFunctions_[ACOU_PMLAUXSCALAR]->SetTimeScheme(scalScheme);
      }
    }
  
    // TDEF pole schemes
    if (this->timeDomainEqFluidFormulation_)
    {
      for (unsigned int i = 0; i < TDEFpoleNumber_.realC.Max(); i++)
        feFunctions_[(SolutionType)(ACOU_TDEF_PHI_C_1 + i)]->SetTimeScheme(
          shared_ptr<BaseTimeScheme>(new TimeSchemeGLM(makeScheme(), 0)));
        
      for (unsigned int i = 0; i < TDEFpoleNumber_.complC.Max(); i++)
        feFunctions_[(SolutionType)(ACOU_TDEF_PSI_C_1 + i)]->SetTimeScheme(
          shared_ptr<BaseTimeScheme>(new TimeSchemeGLM(makeScheme(), 0)));
        
      for (unsigned int i = 0; i < TDEFpoleNumber_.realV.Max(); i++)
        feFunctions_[(SolutionType)(ACOU_TDEF_PHI_V_1 + i)]->SetTimeScheme(
          shared_ptr<BaseTimeScheme>(new TimeSchemeGLM(makeScheme(), 0)));
        
      for (unsigned int i = 0; i < TDEFpoleNumber_.complV.Max(); i++)
        feFunctions_[(SolutionType)(ACOU_TDEF_PSI_V_1 + i)]->SetTimeScheme(
          shared_ptr<BaseTimeScheme>(new TimeSchemeGLM(makeScheme(), 0)));
    }
  }

  void AcousticPDE::ComputeSOS(PtrCoefFct& c0, PtrCoefFct dens, PtrCoefFct blk,
		  PtrCoefFct regionTemp, std::string tempId) {
//TODO 1 erase inconsistent coding style 1/c0^2
	  if ( complexFluidFormulation_ ) {
		  //here c0 is actually 1/c0^2!!
		  c0 = CoefFunction::Generate( mp_,  Global::COMPLEX,
				  CoefXprBinOp(mp_, dens, blk, CoefXpr::OP_DIV ) );
	  }
	  else {
		  if (tempId != "") {
			  if( formulation_ == ACOU_POTENTIAL )
				  EXCEPTION("We need for temperature dependent speed of sound a pressure formulation");

			  // gasR=287.058 J/kg K   ... universal gas constant
			  // kappa=1.402, adabatic exponent for air
        shared_ptr< CoefFunctionConst<Double> > constVal(new CoefFunctionConst<Double>());
        constVal->SetScalar(402.4553160);

			  c0 = CoefFunction::Generate( mp_,  Global::REAL,
					  CoefXprUnaryOp( mp_, CoefXprBinOp(mp_, constVal, regionTemp, CoefXpr::OP_MULT),
							  CoefXpr::OP_SQRT) );
		  }
		  else {
			  // c0 = sqrt(bulk_modulus / density)
			  c0 = CoefFunction::Generate( mp_,  Global::REAL,
					  CoefXprUnaryOp( mp_, CoefXprBinOp(mp_, blk, dens, CoefXpr::OP_DIV),
							  CoefXpr::OP_SQRT) );
		  }
	  }
  }

  void AcousticPDE::ComputeSOS_SQR(PtrCoefFct& cSQR, PtrCoefFct dens, PtrCoefFct blk,
		                           PtrCoefFct regionTemp, std::string tempId) {
    if (tempId != "") {
      if( formulation_ == ACOU_POTENTIAL )
          EXCEPTION("We need for temperature dependent speed of sound a pressure formulation");

      // gasR=287.058 J/kg K   ... universal gas constant
      // kappa=1.402, adabatic exponent for air
      shared_ptr< CoefFunctionConst<Double> > constVal(new CoefFunctionConst<Double>());
      constVal->SetScalar(402.4553160);

      cSQR = CoefFunction::Generate( mp_,  Global::REAL,
                                    CoefXprBinOp(mp_, constVal, regionTemp, CoefXpr::OP_MULT) );
    }
    else {
      // c^2 = bulk_modulus / density
      cSQR = CoefFunction::Generate( mp_, Global::REAL, CoefXprBinOp(mp_, blk, dens, CoefXpr::OP_DIV) );
    }
  }

  void AcousticPDE::CalcMechAcouFac( PtrCoefFct& mechAcouFactor, PtrCoefFct dens){

    if (isMechCoupled_ == true && formulation_ == ACOU_POTENTIAL) {
      // Important: In case of a general / quadratic EV problem, we must
      // ensure to have a "positive definite" matrix, i.e. we are not allowed
      // to multiply all matrices by -1!
      std::string stringFac = (analysistype_ != EIGENFREQUENCY) ? "-1.0" : "1.0";
      mechAcouFactor = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, dens, stringFac, CoefXpr::OP_MULT));
    }
    else {
      mechAcouFactor = CoefFunction::Generate(mp_, Global::REAL, "1.0");
    }
  }

  void AcousticPDE::CalcMechAcouFacWithCoef( PtrCoefFct& exValue, PtrCoefFct coef, PtrCoefFct surfDens, Double& scalFactor, Global::ComplexPart part){
    scalFactor = 1.0;
    if ( isMechCoupled_ == true && formulation_ ==  ACOU_POTENTIAL ) {
      scalFactor = -1.0;
      exValue = CoefFunction::Generate( mp_, part, CoefXprBinOp(mp_, coef,surfDens, CoefXpr::OP_MULT) );
    } else {
      exValue = coef;
    }
  }

  } // namespace CoupledField
  void AcousticPDE::ReadTDEFCoefficients(UInt iRegion, rationalCoeffsTDEF &TDEFcoeffs)
  {
    unsigned int actRegion = regions_[iRegion];
    LocPointMapped lpm;

    Vector<Double> vecAC; // auxiliary vectors to loop over them
    Vector<Double> vecBC;
    Vector<Double> vecCC;
    Vector<Double> vecAlphaC;
    Vector<Double> vecBetaC;
    Vector<Double> vecGammaC;
    Vector<Double> vecAV;
    Vector<Double> vecBV;
    Vector<Double> vecCV;
    Vector<Double> vecAlphaV;
    Vector<Double> vecBetaV;
    Vector<Double> vecGammaV;

    if (TDEFpoleNumber_.realC[iRegion] > 0) {
      // Equivalent compressibility: real poles
      PtrCoefFct coefVecAC = materials_[actRegion]->GetVectorCoefFnc(ACOU_TDEF_INVBLK_A, Global::REAL);
      PtrCoefFct coefVecAlphaC = materials_[actRegion]->GetVectorCoefFnc(ACOU_TDEF_INVBLK_ALPHA, Global::REAL);
      coefVecAC->GetVector(vecAC, lpm);
      coefVecAlphaC->GetVector(vecAlphaC, lpm);
    }
    else {
      vecAC.Resize(0);
      vecAlphaC.Resize(0);
    }

    if (TDEFpoleNumber_.complC[iRegion] > 0) {
      // Equivalent compressibility: complex-conjugated poles
      PtrCoefFct coefVecBC = materials_[actRegion]->GetVectorCoefFnc(ACOU_TDEF_INVBLK_B, Global::REAL);
      PtrCoefFct coefVecCC = materials_[actRegion]->GetVectorCoefFnc(ACOU_TDEF_INVBLK_C, Global::REAL);
      PtrCoefFct coefVecBetaC = materials_[actRegion]->GetVectorCoefFnc(ACOU_TDEF_INVBLK_BETA, Global::REAL);
      PtrCoefFct coefVecGammaC = materials_[actRegion]->GetVectorCoefFnc(ACOU_TDEF_INVBLK_GAMMA, Global::REAL);
      coefVecBC->GetVector(vecBC, lpm);
      coefVecCC->GetVector(vecCC, lpm);
      coefVecBetaC->GetVector(vecBetaC, lpm);
      coefVecGammaC->GetVector(vecGammaC, lpm);
    }
    else {
      vecBC.Resize(0);
      vecCC.Resize(0);
      vecBetaC.Resize(0);
      vecGammaC.Resize(0);
    }

    if (TDEFpoleNumber_.realV[iRegion] > 0) {
      // Equivalent specific volume (inv. density): real poles
      PtrCoefFct coefVecAV = materials_[actRegion]->GetVectorCoefFnc(ACOU_TDEF_INVDENS_A, Global::REAL);
      PtrCoefFct coefVecAlphaV = materials_[actRegion]->GetVectorCoefFnc(ACOU_TDEF_INVDENS_ALPHA, Global::REAL);
      coefVecAV->GetVector(vecAV, lpm);
      coefVecAlphaV->GetVector(vecAlphaV, lpm);
    }
    else {
      vecAV.Resize(0);
      vecAlphaV.Resize(0);
    }

    if (TDEFpoleNumber_.complV[iRegion] > 0) {
      // Equivalent specific volume (inv. density): complex-conjugated poles
      PtrCoefFct coefVecBV = materials_[actRegion]->GetVectorCoefFnc(ACOU_TDEF_INVDENS_B, Global::REAL);
      PtrCoefFct coefVecCV = materials_[actRegion]->GetVectorCoefFnc(ACOU_TDEF_INVDENS_C, Global::REAL);
      PtrCoefFct coefVecBetaV = materials_[actRegion]->GetVectorCoefFnc(ACOU_TDEF_INVDENS_BETA, Global::REAL);
      PtrCoefFct coefVecGammaV = materials_[actRegion]->GetVectorCoefFnc(ACOU_TDEF_INVDENS_GAMMA, Global::REAL);
      coefVecBV->GetVector(vecBV, lpm);
      coefVecCV->GetVector(vecCV, lpm);
      coefVecBetaV->GetVector(vecBetaV, lpm);
      coefVecGammaV->GetVector(vecGammaV, lpm);
    }
    else {
      vecBV.Resize(0);
      vecCV.Resize(0);
      vecBetaV.Resize(0);
      vecGammaV.Resize(0);
    }

    TDEFcoeffs.AC.Resize(vecAC.GetSize());
    TDEFcoeffs.BC.Resize(vecBC.GetSize());
    TDEFcoeffs.CC.Resize(vecCC.GetSize());
    TDEFcoeffs.AlphaC.Resize(vecAlphaC.GetSize());
    TDEFcoeffs.BetaC.Resize(vecBetaC.GetSize());
    TDEFcoeffs.GammaC.Resize(vecGammaC.GetSize());

    TDEFcoeffs.AV.Resize(vecAV.GetSize());
    TDEFcoeffs.BV.Resize(vecBV.GetSize());
    TDEFcoeffs.CV.Resize(vecCV.GetSize());
    TDEFcoeffs.AlphaV.Resize(vecAlphaV.GetSize());
    TDEFcoeffs.BetaV.Resize(vecBetaV.GetSize());
    TDEFcoeffs.GammaV.Resize(vecGammaV.GetSize());

    // Get real poles: eq. compressibility
    for (UInt ii = 0; ii < vecAC.GetSize(); ii++) {
      TDEFcoeffs.AC[ii] = CoefFunction::Generate(mp_, Global::REAL, std::to_string(vecAC[ii]));
    }
    for (UInt ii = 0; ii < vecAlphaC.GetSize(); ii++) {
      TDEFcoeffs.AlphaC[ii] = CoefFunction::Generate(mp_, Global::REAL, std::to_string(vecAlphaC[ii]));
    }

    // Get real poles: eq. specific volume
    for (UInt ii = 0; ii < vecAV.GetSize(); ii++) {
      TDEFcoeffs.AV[ii] = CoefFunction::Generate(mp_, Global::REAL, std::to_string(vecAV[ii]));
    }
    for (UInt ii = 0; ii < vecAlphaV.GetSize(); ii++) {
      TDEFcoeffs.AlphaV[ii] = CoefFunction::Generate(mp_, Global::REAL, std::to_string(vecAlphaV[ii]));
    }

    // Get complex-conjugated poles: compressibility
    for (UInt ii = 0; ii < vecBC.GetSize(); ii++) {
      TDEFcoeffs.BC[ii] = CoefFunction::Generate(mp_, Global::REAL, std::to_string(vecBC[ii]));
    }
    for (UInt ii = 0; ii < vecCC.GetSize(); ii++) {
      TDEFcoeffs.CC[ii] = CoefFunction::Generate(mp_, Global::REAL, std::to_string(vecCC[ii]));
    }
    for (UInt ii = 0; ii < vecBetaC.GetSize(); ii++) {
      TDEFcoeffs.BetaC[ii] = CoefFunction::Generate(mp_, Global::REAL, std::to_string(vecBetaC[ii]));
    }
    for (UInt ii = 0; ii < vecGammaC.GetSize(); ii++) {
      TDEFcoeffs.GammaC[ii] = CoefFunction::Generate(mp_, Global::REAL, std::to_string(vecGammaC[ii]));
    }

    // Get complex-conjugated poles: specific volume
    for (UInt ii = 0; ii < vecBV.GetSize(); ii++) {
      TDEFcoeffs.BV[ii] = CoefFunction::Generate(mp_, Global::REAL, std::to_string(vecBV[ii]));
    }
    for (UInt ii = 0; ii < vecCV.GetSize(); ii++) {
      TDEFcoeffs.CV[ii] = CoefFunction::Generate(mp_, Global::REAL, std::to_string(vecCV[ii]));
    }
    for (UInt ii = 0; ii < vecBetaV.GetSize(); ii++) {
      TDEFcoeffs.BetaV[ii] = CoefFunction::Generate(mp_, Global::REAL, std::to_string(vecBetaV[ii]));
    }
    for (UInt ii = 0; ii < vecGammaV.GetSize(); ii++) {
      TDEFcoeffs.GammaV[ii] = CoefFunction::Generate(mp_, Global::REAL, std::to_string(vecGammaV[ii]));
    }

    // check if all coefFunctions vectors are smaller than 15 (FE unknowns are limited to max. 15)
    if (TDEFcoeffs.AC.GetSize() > 15 || TDEFcoeffs.BC.GetSize() > 15 || TDEFcoeffs.CC.GetSize() > 15 || TDEFcoeffs.AlphaC.GetSize() > 15 || TDEFcoeffs.BetaC.GetSize() > 15 || TDEFcoeffs.GammaC.GetSize() > 15 ||
        TDEFcoeffs.AV.GetSize() > 15 || TDEFcoeffs.BV.GetSize() > 15 || TDEFcoeffs.CV.GetSize() > 15 || TDEFcoeffs.AlphaV.GetSize() > 15 || TDEFcoeffs.BetaV.GetSize() > 15 || TDEFcoeffs.GammaV.GetSize() > 15) {
      EXCEPTION("TDEF: Only 15 poles are allowed, please reduce the number of poles!");
    }
  }

  void AcousticPDE::EvalTDEFRationalFncs(UInt iRegion, PtrCoefFct actFreq, rationalCoeffsTDEF &TDEFcoeffs)
  {
    ReadTDEFCoefficients(iRegion, TDEFcoeffs);

    PtrCoefFct imagUnit = CoefFunction::Generate(mp_, Global::COMPLEX, "0.0", "1.0");
    PtrCoefFct constOne = CoefFunction::Generate(mp_, Global::COMPLEX, "1.0", "0.0");
    PtrCoefFct omega = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, actFreq, CoefFunction::Generate(mp_, Global::REAL, "2*pi"), CoefXpr::OP_MULT));
    PtrCoefFct omegaIm = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, imagUnit, omega, CoefXpr::OP_MULT));

    PtrCoefFct fncTerm;
    PtrCoefFct fncDenom;
    PtrCoefFct fncImag;
    PtrCoefFct fncNumerator;

    // initialize the sums with the high frequency limits
    unsigned int aRegion = regions_[iRegion];
    invTDEFBlk_ = materials_[aRegion]->GetScalCoefFnc(ACOU_TDEF_INVBLK_CONST, Global::COMPLEX);
    invTDEFDens_ = materials_[aRegion]->GetScalCoefFnc(ACOU_TDEF_INVDENS_CONST, Global::COMPLEX);

    // sum over real poles of inverse bulk modulus
    for (UInt ii = 0; ii < TDEFpoleNumber_.realC[iRegion]; ii++) {
      fncDenom = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, TDEFcoeffs.AlphaC[ii], omegaIm, CoefXpr::OP_ADD));
      fncTerm = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, TDEFcoeffs.AC[ii], fncDenom, CoefXpr::OP_DIV));
      invTDEFBlk_ = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, invTDEFBlk_, fncTerm, CoefXpr::OP_ADD));
    }

    // sum over complex-conjugated poles of inverse bulk modulus
    for (UInt ii = 0; ii < TDEFpoleNumber_.complC[iRegion]; ii++) {

      // pole 1 of complex conjugated pole pair
      fncImag = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, TDEFcoeffs.CC[ii], imagUnit, CoefXpr::OP_MULT));
      fncNumerator = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, TDEFcoeffs.BC[ii], fncImag, CoefXpr::OP_ADD));
      fncImag = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, TDEFcoeffs.GammaC[ii], imagUnit, CoefXpr::OP_MULT));
      fncDenom = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, TDEFcoeffs.BetaC[ii], fncImag, CoefXpr::OP_ADD));
      fncDenom = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, fncDenom, omegaIm, CoefXpr::OP_ADD));
      fncTerm = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, fncNumerator, fncDenom, CoefXpr::OP_DIV));
      invTDEFBlk_ = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, invTDEFBlk_, fncTerm, CoefXpr::OP_ADD));

      // pole 2 of complex conjugated pole pair
      fncImag = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, TDEFcoeffs.CC[ii], imagUnit, CoefXpr::OP_MULT));
      fncNumerator = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, TDEFcoeffs.BC[ii], fncImag, CoefXpr::OP_SUB));
      fncImag = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, TDEFcoeffs.GammaC[ii], imagUnit, CoefXpr::OP_MULT));
      fncDenom = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, TDEFcoeffs.BetaC[ii], fncImag, CoefXpr::OP_SUB));
      fncDenom = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, fncDenom, omegaIm, CoefXpr::OP_ADD));
      fncTerm = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, fncNumerator, fncDenom, CoefXpr::OP_DIV));
      invTDEFBlk_ = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, invTDEFBlk_, fncTerm, CoefXpr::OP_ADD));
    }

    // sum over real poles of inverse density
    for (UInt ii = 0; ii < TDEFpoleNumber_.realV[iRegion]; ii++) {
      fncDenom = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, TDEFcoeffs.AlphaV[ii], omegaIm, CoefXpr::OP_ADD));
      fncTerm = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, TDEFcoeffs.AV[ii], fncDenom, CoefXpr::OP_DIV));
      invTDEFDens_ = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, invTDEFDens_, fncTerm, CoefXpr::OP_ADD));
    }

    // sum over complex poles of inverse density
    for (UInt ii = 0; ii < TDEFpoleNumber_.complV[iRegion]; ii++) {

      // pole 1 of complex conjugated pole pair
      fncImag = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, TDEFcoeffs.CV[ii], imagUnit, CoefXpr::OP_MULT));
      fncNumerator = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, TDEFcoeffs.BV[ii], fncImag, CoefXpr::OP_ADD));
      fncImag = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, TDEFcoeffs.GammaV[ii], imagUnit, CoefXpr::OP_MULT));
      fncDenom = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, TDEFcoeffs.BetaV[ii], fncImag, CoefXpr::OP_ADD));
      fncDenom = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, fncDenom, omegaIm, CoefXpr::OP_ADD));
      fncTerm = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, fncNumerator, fncDenom, CoefXpr::OP_DIV));
      invTDEFDens_ = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, invTDEFDens_, fncTerm, CoefXpr::OP_ADD));

      // pole 2 of complex conjugated pole pair
      fncImag = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, TDEFcoeffs.CV[ii], imagUnit, CoefXpr::OP_MULT));
      fncNumerator = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, TDEFcoeffs.BV[ii], fncImag, CoefXpr::OP_SUB));
      fncImag = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, TDEFcoeffs.GammaV[ii], imagUnit, CoefXpr::OP_MULT));
      fncDenom = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, TDEFcoeffs.BetaV[ii], fncImag, CoefXpr::OP_SUB));
      fncDenom = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, fncDenom, omegaIm, CoefXpr::OP_ADD));
      fncTerm = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, fncNumerator, fncDenom, CoefXpr::OP_DIV));
      invTDEFDens_ = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, invTDEFDens_, fncTerm, CoefXpr::OP_ADD));
    }
  }

template void AcousticPDE::DefineTransientPMLInts<2>(shared_ptr<ElemList>, std::string,
    RegionIdType actRegion, std::string tempId);
template void AcousticPDE::DefineTransientPMLInts<3>(shared_ptr<ElemList>, std::string,
    RegionIdType actRegion, std::string tempId);
template void AcousticPDE::DefineConvectiveIntegrators<2, true>(RegionIdType actRegion, PtrParamNode curRegNode, 
    shared_ptr<ElemList> actSDList, PtrCoefFct coeffM);
template void AcousticPDE::DefineConvectiveIntegrators<2, false>(RegionIdType actRegion, PtrParamNode curRegNode, 
    shared_ptr<ElemList> actSDList, PtrCoefFct coeffM);
template void AcousticPDE::DefineConvectiveIntegrators<3, true>(RegionIdType actRegion, PtrParamNode curRegNode, 
    shared_ptr<ElemList> actSDList, PtrCoefFct coeffM);
template void AcousticPDE::DefineConvectiveIntegrators<3, false>(RegionIdType actRegion, PtrParamNode curRegNode, 
    shared_ptr<ElemList> actSDList, PtrCoefFct coeffM);
template void AcousticPDE::DefinePMLIntegrators<2>(RegionIdType actRegion, shared_ptr<ElemList>& actSDList,
    PtrParamNode& curRegNode, PtrCoefFct& c0, PtrCoefFct& coeffK, PtrCoefFct& coeffM, std::string& tempId,
    BaseBDBInt*& stiffInt, BaseBDBInt*& massInt);
template void AcousticPDE::DefinePMLIntegrators<3>(RegionIdType actRegion, shared_ptr<ElemList>& actSDList,
    PtrParamNode& curRegNode, PtrCoefFct& c0, PtrCoefFct& coeffK, PtrCoefFct& coeffM, std::string& tempId,
    BaseBDBInt*& stiffInt, BaseBDBInt*& massInt);
