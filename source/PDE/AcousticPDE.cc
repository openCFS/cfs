#include "AcousticPDE.hh"

#include "General/defs.hh"

#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

// new integrator concept
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
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/CoefFunction/CoefFunctionImpedanceModel.hh"
#include "Domain/Mesh/NcInterfaces/BaseNcInterface.hh"

#include <boost/lexical_cast.hpp>
#include <cmath>

#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"
#include "Materials/AcousticMaterial.hh"

namespace CoupledField
{

  DEFINE_LOG(acousticpde, "pde.acoustic")

  AcousticPDE::AcousticPDE(Grid *aGrid, PtrParamNode paramNode,
                           PtrParamNode infoNode,
                           shared_ptr<SimState> simState, Domain *domain)
      : SinglePDE(aGrid, paramNode, infoNode, simState, domain)
  {

    pdename_ = "acoustic";
    pdematerialclass_ = ACOUSTIC;
    nonLin_ = false;
    isMechCoupled_ = false;
    formulation_ = ACOU_PRESSURE;

    //! Always use total Lagrangian formulation
    updatedGeo_ = false;
    isTimeDomPML_ = false;
    isAPML_ = false;
    complexFluidFormulation_ = false;
    timeDomainEqFluidFormulation_ = false;

    std::string pdeFormulation = myParam_->Get("formulation")->As<std::string>();
    // check for pressure or potential formulation
    sosAtLaplace_ = (pdeFormulation == "acouPressureSOSatLaplace") ? true : false;

    if (pdeFormulation != "default" && pdeFormulation != "acouPressureSOSatLaplace")
    {
      formulation_ = SolutionTypeEnum.Parse(pdeFormulation);
    }
  }

  std::map<SolutionType, shared_ptr<FeSpace>> AcousticPDE::CreateFeSpaces(const std::string &formulation,
                                                                          PtrParamNode infoNode)
  {

    if (this->analysistype_ == STATIC)
      EXCEPTION("No STATIC analysis in AcousticPDE");

    std::map<SolutionType, shared_ptr<FeSpace>> crSpaces;
    if (formulation == "default" || formulation == "H1")
    {
      std::string form = SolutionTypeEnum.ToString(formulation_);
      PtrParamNode potSpaceNode = infoNode->Get(form);
      crSpaces[formulation_] =
          FeSpace::CreateInstance(myParam_, potSpaceNode, FeSpace::H1, ptGrid_);
      crSpaces[formulation_]->Init(solStrat_);
    }
    else
    {
      EXCEPTION("The formulation " << formulation << "of acoustic PDE is not known!");
    }

    // This check has been moved from DefinePrimaryResults since it is needed earlier for the TDEF formulation
    RegionIdType actRegion;
    std::map<RegionIdType, BaseMaterial *>::iterator it;
    for (it = materials_.begin(); it != materials_.end(); it++)
    {
      // check for complex fluid formulation
      actRegion = it->first;
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);
      PtrParamNode curRegNode =
          myParam_->Get("regionList")->GetByVal("region", "name", regionName.c_str());
      if (curRegNode->Get("complexFluid")->As<std::string>() == "yes")
      {
        complexFluidFormulation_ = true;
        isMaterialComplex_ = true;
        if (this->analysistype_ != HARMONIC)
          EXCEPTION("Complex fluid region just allowed in harmonic analysis");
        // need an acoustic pressure formulation
        if (formulation_ != ACOU_PRESSURE)
          EXCEPTION("Complex fluid needs acoustic pressure formulation");
      }

      // check for time domain equivalent fluid (TDEF) formulation
      if (curRegNode->Get("timeDomainEqFluid")->As<std::string>() == "yes")
      {
        timeDomainEqFluidFormulation_ = true;
        if (this->analysistype_ != TRANSIENT)
          EXCEPTION("Time domain equivalent fluid formulation only possible in transient analysis");
        // need an acoustic pressure formulation
        if (formulation_ != ACOU_PRESSURE)
          EXCEPTION("Time domain equivalent fluid formulation needs acoustic pressure formulation");
      }
    }

    // ===================================
    // Check for transient PML
    // ===================================
    if (this->analysistype_ == TRANSIENT && isTimeDomPML_)
    {
      // now define the additional unknowns

      if (dim_ == 3 && !this->isAPML_)
      {
        PtrParamNode scalarpml = infoNode->Get("TransientPMLScalarAuxVar");
        crSpaces[ACOU_PMLAUXSCALAR] =
            FeSpace::CreateInstance(myParam_, scalarpml, FeSpace::H1, ptGrid_);
        crSpaces[ACOU_PMLAUXSCALAR]->Init(solStrat_);
      }

      PtrParamNode vectorPML = infoNode->Get("TransientPMLVectorAuxVar");
      crSpaces[ACOU_PMLAUXVEC] =
          FeSpace::CreateInstance(myParam_, vectorPML, FeSpace::H1, ptGrid_);
      crSpaces[ACOU_PMLAUXVEC]->Init(solStrat_);
    }

    // ===================================
    // Check for transient TDEF
    // ===================================
    if (this->analysistype_ == TRANSIENT && timeDomainEqFluidFormulation_ && this->formulation_ == ACOU_PRESSURE)
    {
      // check for the actual size of the auxiliary variables
      nAuxFncAC_.Resize(regions_.GetSize());
      nAuxFncBC_.Resize(regions_.GetSize());
      nAuxFncAV_.Resize(regions_.GetSize());
      nAuxFncBV_.Resize(regions_.GetSize());
      isTDEFReg_.Resize(regions_.GetSize());

      RegionIdType actRegion;
      for (UInt iRegion = 0; iRegion < regions_.GetSize(); iRegion++)
      {
        actRegion = regions_[iRegion];
        std::string regionName = ptGrid_->GetRegion().ToString(actRegion);
        PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region", "name", regionName.c_str());

        if (curRegNode->Get("timeDomainEqFluid")->As<std::string>() == "yes")
        {
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
          nAuxFncAC_[iRegion] = vecAC.GetSize();
          nAuxFncBC_[iRegion] = vecBC.GetSize();
          nAuxFncAV_[iRegion] = vecAV.GetSize();
          nAuxFncBV_[iRegion] = vecBV.GetSize();

        } else {
          // this region has no TDEF material defined, set size to 0
          nAuxFncAC_[iRegion] = 0;
          nAuxFncBC_[iRegion] = 0;
          nAuxFncAV_[iRegion] = 0;
          nAuxFncBV_[iRegion] = 0;
        }

        if (nAuxFncAC_[iRegion]==0 && nAuxFncBC_[iRegion]==0 && nAuxFncAV_[iRegion]==0 && nAuxFncBV_[iRegion]==0){
          isTDEFReg_[iRegion] = false;
        }
        else{
          isTDEFReg_[iRegion] = true;
        }

      }
      
      // define the additional unknowns (phiC, psiC, phiV, psiV)
      PtrParamNode spaceNode;
      for (unsigned int i = 0; i < nAuxFncAC_.Max(); i++)
      {
        spaceNode = infoNode->Get(SolutionTypeEnum.ToString((SolutionType)(ACOU_TDEF_PHI_C_1 + i)));
        crSpaces[(SolutionType)(ACOU_TDEF_PHI_C_1 + i)] = FeSpace::CreateInstance(myParam_, spaceNode, FeSpace::H1, ptGrid_);
        crSpaces[(SolutionType)(ACOU_TDEF_PHI_C_1 + i)]->Init(solStrat_);
      }
      for (unsigned int i = 0; i < nAuxFncBC_.Max(); i++)
      {
        spaceNode = infoNode->Get(SolutionTypeEnum.ToString((SolutionType)(ACOU_TDEF_PSI_C_1 + i)));
        crSpaces[(SolutionType)(ACOU_TDEF_PSI_C_1 + i)] = FeSpace::CreateInstance(myParam_, spaceNode, FeSpace::H1, ptGrid_);
        crSpaces[(SolutionType)(ACOU_TDEF_PSI_C_1 + i)]->Init(solStrat_);
      }
      for (unsigned int i = 0; i < nAuxFncAV_.Max(); i++)
      {
        spaceNode = infoNode->Get(SolutionTypeEnum.ToString((SolutionType)(ACOU_TDEF_PHI_V_1 + i)));
        crSpaces[(SolutionType)(ACOU_TDEF_PHI_V_1 + i)] = FeSpace::CreateInstance(myParam_, spaceNode, FeSpace::H1, ptGrid_);
        crSpaces[(SolutionType)(ACOU_TDEF_PHI_V_1 + i)]->Init(solStrat_);
      }
      for (unsigned int i = 0; i < nAuxFncBV_.Max(); i++)
      {
        spaceNode = infoNode->Get(SolutionTypeEnum.ToString((SolutionType)(ACOU_TDEF_PSI_V_1 + i)));
        crSpaces[(SolutionType)(ACOU_TDEF_PSI_V_1 + i)] = FeSpace::CreateInstance(myParam_, spaceNode, FeSpace::H1, ptGrid_);
        crSpaces[(SolutionType)(ACOU_TDEF_PSI_V_1 + i)]->Init(solStrat_);
      }
    }

    return crSpaces;
  }

  void AcousticPDE::ReadDampingInformation()
  {
    std::map<std::string, DampingType> idDampType;
    std::map<std::string, shared_ptr<RaylDampingData>> idRaylData;

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

        if (actType == RAYLEIGH)
        {
          // set data for Rayleigh damping
          shared_ptr<RaylDampingData> actRaylDamp(new RaylDampingData());
          actRaylDamp->alpha = "0.0";
          actRaylDamp->beta = "0.0";
          actRaylDamp->adjustDamping = true;
          actRaylDamp->ratioDeltaF = 0.01;
          actRaylDamp->freq = 0.0;
          dampNodes[i]->GetValue("freq", actRaylDamp->freq, ParamNode::PASS);
          dampNodes[i]->GetValue("ratioDeltaF", actRaylDamp->ratioDeltaF,
                                 ParamNode::PASS);
          dampNodes[i]->GetValue("adjustDamping", actRaylDamp->adjustDamping,
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
        EXCEPTION("Damping with id '" << actDampingId
                                      << "' was not defined in 'dampingList'");
      }

      dampingList_[actRegionId] = idDampType[actDampingId];
      if (dampingList_[actRegionId] == RAYLEIGH)
      {
        RaylDampingData actRayl = *(idRaylData[actDampingId]);
        Double dampFreq;

        if (actRayl.freq == 0.0)
        {
          materials_[actRegionId]->GetScalar(dampFreq, RAYLEIGH_FREQUENCY, Global::REAL);
        }
        else
        {
          dampFreq = actRayl.freq;
        }
        // Compute Rayleigh damping parameters
        materials_[actRegionId]->ComputeRayleighDamping(actRayl.alpha, actRayl.beta,
                                                        dampFreq, actRayl.ratioDeltaF,
                                                        actRayl.adjustDamping, isComplex_);
        regionRaylDamping_[actRegionId] = actRayl;
      }
      else if (dampingList_[actRegionId] == PML &&
               analysistype_ == BasePDE::TRANSIENT)
      {
        isTimeDomPML_ = true;
      }
    }
  }

  /// @brief 
  void AcousticPDE::DefineIntegrators()
  {

    RegionIdType actRegion;
    ;

    // flag that states if we are in a PML region while performing a harmonic analysis
    bool harmonicPML = false;

    // type of geometry
    isaxi_ = ptGrid_->IsAxi();

    // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial *>::iterator it;
    shared_ptr<FeSpace> mySpace = feFunctions_[formulation_]->GetFeSpace();

    LOG_DBG(acousticpde) << "DefineIntegrators BEGIN"
                         << "\n";

    for (UInt iRegion = 0; iRegion < regions_.GetSize(); iRegion++)
    {
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
      mySpace->SetRegionApproximation(actRegion, polyId, integId);

      //=======================================================================
      // Generate coefficient functions 4 (Speed of Sound)
      //=======================================================================
      PtrCoefFct c0;
      PtrCoefFct dens;
      PtrCoefFct blk;
      PtrCoefFct constOne = CoefFunction::Generate(mp_, Global::REAL, "1.0");

      if (complexFluidFormulation_ && useRationalAppr == "false")
      {
        dens = materials_[actRegion]->GetScalCoefFnc(DENSITY, Global::COMPLEX);
        blk = materials_[actRegion]->GetScalCoefFnc(ACOU_BULK_MODULUS, Global::COMPLEX);
      }
      else if (complexFluidFormulation_ && useRationalAppr == "true"){
        std::cout << "Using the rational function approximation prvovided by material-file for bulk modulus and density." << std::endl;


//########################################################################

         //TODO: get current frequency


          // the following passage is taken from CoefFunctionPML.cc



        //this is just to be up to date with the desired frequency!
        // obtain handle from internal variable coefficient function


        // mp_ = domain->GetMathParser();
        // mHandle_ = mp_->GetNewHandle(true);

        // mp_->SetExpr(mHandle_,"f");

        // // register callback mechanism if expression changes
        // mp_->AddExpChangeCallBack(
        //     boost::bind(&AcousticPDE::UpdateFreq, this ),
        //     mHandle_ );
        // // important: Trigger first-time calculation
        // UpdateFreq();


        // //It might be necessary to just disconnect the callback instead of releasing the handle
        // mp_->ReleaseHandle( mHandle_ );


//########################################################################


        EvalRationalFncs(iRegion, freq_);
        dens = CoefFunction::Generate(mp_, Global::COMPLEX,
                                      CoefXprBinOp(mp_, constOne, invTDEFDens_, CoefXpr::OP_DIV));
        blk = CoefFunction::Generate(mp_, Global::COMPLEX,
                                      CoefXprBinOp(mp_, constOne, invTDEFBlk_, CoefXpr::OP_DIV));

      }

      else if (timeDomainEqFluidFormulation_)
      {
        // Attention:
        // in case of the TDEF formulation, the high frequency limit is set for ACOU_BULK_MODULUS and DENSITY (for user input checks)
        // We also created ACOU_TDEF_INVDENS_CONST and ACOU_TDEF_INVBLK_CONST, which have the inverse values (we currently don't use them)
        PtrCoefFct densInvConst;
        PtrCoefFct blkInvConst;
        densInvConst = materials_[actRegion]->GetScalCoefFnc(DENSITY, Global::REAL);
        blkInvConst = materials_[actRegion]->GetScalCoefFnc(ACOU_BULK_MODULUS, Global::REAL);
        dens = CoefFunction::Generate(mp_, Global::REAL,
                                      CoefXprBinOp(mp_, constOne, densInvConst, CoefXpr::OP_DIV));
        blk = CoefFunction::Generate(mp_, Global::REAL,
                                     CoefXprBinOp(mp_, constOne, blkInvConst, CoefXpr::OP_DIV));                         
      }
      else
      {
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

      if (tempId != "")
      {
        std::set<UInt> definedDofs;
        // just real valued temperature allowed
        shared_ptr<ResultInfo> tempInfo = GetResultInfo(HEAT_MEAN_TEMPERATURE);
        // Add the region information
        PtrParamNode tempNode = myParam_->Get("temperatureList")->GetByVal("temperature", "name", tempId.c_str());

        ReadUserFieldValues(actSDList, tempNode, tempInfo->dofNames, tempInfo->entryType,
                            false, regionTemp, definedDofs, updatedGeo_);
        meanTemperatureCoef_->AddRegion(actRegion, regionTemp);
      }

      // compute speed of sound
      ComputeSOS(c0, dens, blk, regionTemp, tempId);
      // store coefficient functions
      matCoefs_[ELEM_DENSITY]->AddRegion(actRegion, dens);
      matCoefs_[ACOU_ELEM_SPEED_OF_SOUND]->AddRegion(actRegion, c0);

      // if pde couples with mechanic, we have to multiply the density by -1
      PtrCoefFct factor;
      if (isMechCoupled_ == true && formulation_ != ACOU_PRESSURE)
      {
        if (complexFluidFormulation_ || timeDomainEqFluidFormulation_)
          EXCEPTION("Complex fluid and coupled mechanical-acoustic simulation not allowed");

        if (timeDomainEqFluidFormulation_)
          EXCEPTION("TDEF formulation and coupled mechanical-acoustic simulation currently not allowed");
          // ATTENTION: in the case of mechanical excitation of a porous material, the assumption of a rigid
          // frame (skeleton) of the material is not valid anymore in a strict sense. Instead of for example 
          // JCAL, Biot's theory must be applied accounting for coupling of the solid and fluid phase.

        // Important: In case of a general / quadratic EV problem, we must
        // ensure to have a "positive definite" matrix, i.e. we are not allowed
        // to multiply all matrices by -1!
        std::string stringFac = (analysistype_ != EIGENFREQUENCY) ? "-1.0" : "1.0";

        factor = CoefFunction::Generate(mp_, Global::REAL,
                                        CoefXprBinOp(mp_, dens, stringFac, CoefXpr::OP_MULT));
      }
      else
      {
        factor = constOne;
      }

      // basic coeff-functions for mass and stiffness matrix
      PtrCoefFct coeffM, coeffK;

      if (formulation_ == ACOU_PRESSURE && sosAtLaplace_ == true)
      {
        if (complexFluidFormulation_ || timeDomainEqFluidFormulation_)
          EXCEPTION("A complex/equivalent fluid and sosAtLaplace-formulation not allowed!!");

        // pressure formulation with temperature depend speed of sound
        coeffM = constOne;
        coeffK = CoefFunction::Generate(mp_, Global::REAL,
                                        CoefXprBinOp(mp_, c0, c0, CoefXpr::OP_MULT));
      }
      else
      {
        if (complexFluidFormulation_)
        {
          // in this case c0 is actually 1/c0^2!!
          //! coeffM: 1/compressionModulus
          //! coeffK = 1/density
          coeffM = CoefFunction::Generate(mp_, Global::COMPLEX,
                                          CoefXprBinOp(mp_, factor, blk, CoefXpr::OP_DIV));
          /// coeffM: 1/density
          coeffK = CoefFunction::Generate(mp_, Global::COMPLEX,
                                          CoefXprBinOp(mp_, factor, dens, CoefXpr::OP_DIV));
        }
        else if (timeDomainEqFluidFormulation_)
        {
          // Here, dens and blk are already the inverse of the high-frequency limit of the 
          // rational function approximation.
          // As density and blk modulus they occur in their inverse form, the inverse density
          // comes to the stiffness matrix and the inverse blk modulus to the mass matrix
          coeffK = dens;   
          coeffM = blk;
        }
        
        else
        {
          coeffK = factor;
          // build coefficient for mass matrix as (factor / (c0*c0))
          PtrCoefFct constValOne = CoefFunction::Generate(mp_, Global::REAL, "1.0");
          PtrCoefFct coeffa, coeffb;
          coeffa = CoefFunction::Generate(mp_, Global::REAL,
                                          CoefXprBinOp(mp_, constValOne, c0, CoefXpr::OP_DIV));
          coeffb = CoefFunction::Generate(mp_, Global::REAL,
                                          CoefXprBinOp(mp_, coeffa, coeffa, CoefXpr::OP_MULT));
          coeffM = CoefFunction::Generate(mp_, Global::REAL,
                                          CoefXprBinOp(mp_, factor, coeffb, CoefXpr::OP_MULT));
        }
      }

      LOG_DBG(acousticpde) << "DefineIntegrators Fluid: coeffK = " << coeffK->ToString() << "\n";
      LOG_DBG(acousticpde) << "DefineIntegrators Fluid: coeffM = " << coeffM->ToString() << "\n";

      // ====================================================================
      // Take account for pml
      // ====================================================================
      shared_ptr<CoefFunction> coeffPMLScal, coeffPMLVec;
      shared_ptr<CoefFunction> coeffPMLStiff;
      shared_ptr<CoefFunction> coeffPMLMass;

      if (dampingList_[actRegion] == PML)
      {
        std::string dampId;
        curRegNode->GetValue("dampingId", dampId);

        // check for PML formulation: classic(=Cartesian) or curvilinear
        std::string pmlFormul;
        PtrParamNode pmlNode = myParam_->Get("dampingList")->GetByVal("pml", "id", dampId.c_str());
        pmlFormul = pmlNode->Get("formulation")->As<std::string>();

        // check if harmonic or inverse-source analysis is performed
        if (analysistype_ == HARMONIC || analysistype_ == BasePDE::INVERSESOURCE)
        {
          harmonicPML = true;
          shared_ptr<CoefFunction> densR, blkR, c0R;
          if (complexFluidFormulation_ == true) // in case of complexFluid, speed of sound is complex; so we need to compute a real-valued one for the PML definition
          {
            densR = materials_[actRegion]->GetScalCoefFnc(DENSITY, Global::REAL);
            blkR = materials_[actRegion]->GetScalCoefFnc(ACOU_BULK_MODULUS, Global::REAL);
            c0R = CoefFunction::Generate(mp_, Global::REAL,
                                         CoefXprUnaryOp(mp_, CoefXprBinOp(mp_, blkR, densR, CoefXpr::OP_DIV),
                                                        CoefXpr::OP_SQRT));
          }
          else
            c0R = c0;

          if (pmlFormul == "classic")
          {
            coeffPMLVec.reset(new CoefFunctionPML<Complex>(pmlNode, c0R, actSDList, regions_, true));
            coeffPMLScal.reset(new CoefFunctionPML<Complex>(pmlNode, c0R, actSDList, regions_, false));
            // store pml factor
            matCoefs_[PML_DAMP_FACTOR]->AddRegion(actRegion, coeffPMLVec);
            coeffPMLStiff = CoefFunction::Generate(mp_, Global::COMPLEX,
                                                   CoefXprBinOp(mp_, coeffPMLScal, coeffK, CoefXpr::OP_MULT));

            coeffPMLMass = CoefFunction::Generate(mp_, Global::COMPLEX,
                                                  CoefXprBinOp(mp_, coeffPMLScal, coeffM, CoefXpr::OP_MULT));
          }
          else if (pmlFormul == "curvilinear")
          {
            EXCEPTION("Formulation '" << pmlFormul << "' for AcousticPDE "
                                      << "is not implemented yet!")
            // here, I need to define my PML coefficients....

            // at first, I need to determine the geometry...

            // somehow, I need to take account for the new Jakobi matrix, as this one will not be diagonal..

            // the Jakobi determinant will likely have the same properties as the Cartesian one (but is different in its value)
          }
          else // when pmlFormul is invalid...
          {
            EXCEPTION("Unknown PML-formulation '" << pmlFormul << "' for AcousticPDE. "
                                                  << "Possible formulations are: 'classic' (default), 'curvilinear';")
          }
        }
        else // if not harmonic, define the transient integrators
        {
          harmonicPML = false;
          if (pmlFormul == "classic")
          {
            if (dim_ == 2)
              DefineTransientPMLInts<2>(actSDList, dampId, actRegion, tempId);
            else
              DefineTransientPMLInts<3>(actSDList, dampId, actRegion, tempId);
          }
          else
            EXCEPTION("Transient PML is currently only implemented in 'classic' formulation.")
        }
      }
      else // (re-)set harmonicPML to false if we are not currently in a PML region
        harmonicPML = false;

      // ====================================================================
      // standard stiffness integrator
      // ====================================================================
      BaseBDBInt *stiffInt = NULL;
      if (dim_ == 2)
      {
        if (harmonicPML)
        {
          stiffInt = new BBInt<Complex>(new ScaledGradientOperator<FeH1, 2, Complex>(),
                                        coeffPMLStiff, 1.0, updatedGeo_);
          stiffInt->SetBCoefFunctionOpB(coeffPMLVec);
        }
        else if (complexFluidFormulation_)
        {
          stiffInt = new BBInt<Complex>(new GradientOperator<FeH1, 2>(), coeffK,
                                        1.0, updatedGeo_);
        }
        else
        {
          stiffInt = new BBInt<Double>(new GradientOperator<FeH1, 2>(), coeffK,
                                       1.0, updatedGeo_);
        }
      }
      else
      {
        if (harmonicPML)
        {
          stiffInt = new BBInt<Complex>(new ScaledGradientOperator<FeH1, 3, Complex>(),
                                        coeffPMLStiff, 1.0, updatedGeo_);
          stiffInt->SetBCoefFunctionOpB(coeffPMLVec);
        }
        else if (complexFluidFormulation_)
        {
          stiffInt = new BBInt<Complex>(new GradientOperator<FeH1, 3>(), coeffK,
                                        1.0, updatedGeo_);
        }
        else
        {
          stiffInt = new BBInt<Double>(new GradientOperator<FeH1, 3>(), coeffK,
                                       1.0, updatedGeo_);
        }
      }

      stiffInt->SetName("LaplaceIntegrator");

      BiLinFormContext *stiffIntDescr =
          new BiLinFormContext(stiffInt, STIFFNESS);

      // check for damping
      if (dampingList_[actRegion] == RAYLEIGH)
      {
        if (complexFluidFormulation_ || timeDomainEqFluidFormulation_)
          EXCEPTION("Complex fluid region and Rayleigh damping not allowed!!");

        RaylDampingData &actDamp = (regionRaylDamping_[actRegion]);
        stiffIntDescr->SetSecDestMat(DAMPING, actDamp.beta);
      }

      feFunctions_[formulation_]->AddEntityList(actSDList);

      stiffIntDescr->SetEntities(actSDList, actSDList);
      stiffIntDescr->SetFeFunctions(feFunctions_[formulation_], feFunctions_[formulation_]);
      stiffInt->SetFeSpace(feFunctions_[formulation_]->GetFeSpace());

      assemble_->AddBiLinearForm(stiffIntDescr);
      // Important: Add bdb-integrator to global list, as we need them later
      // for calculation of postprocessing results
      bdbInts_.insert(std::pair<RegionIdType, BaseBDBInt *>(actRegion, stiffInt));

      // ====================================================================
      // standard mass integrator
      // ====================================================================

      BaseBDBInt *massInt = NULL;

      if (dim_ == 2)
      {
        if (harmonicPML)
          massInt = new BBInt<Complex>(new IdentityOperator<FeH1, 2, 1, Complex>(),
                                       coeffPMLMass, 1.0, updatedGeo_);
        else if (complexFluidFormulation_)
          massInt = new BBInt<Complex>(new IdentityOperator<FeH1, 2, 1, Double>, coeffM,
                                       1.0, updatedGeo_);
        else
          massInt = new BBInt<Double>(new IdentityOperator<FeH1, 2, 1, Double>, coeffM,
                                      1.0, updatedGeo_);
      }
      else
      {
        if (harmonicPML)
          massInt = new BBInt<Complex>(new IdentityOperator<FeH1, 3, 1, Complex>(), coeffPMLMass,
                                       1.0, updatedGeo_);
        else if (complexFluidFormulation_)
          massInt = new BBInt<Complex>(new IdentityOperator<FeH1, 3, 1, Double>, coeffM,
                                       1.0, updatedGeo_);
        else
          massInt = new BBInt<Double>(new IdentityOperator<FeH1, 3, 1, Double>, coeffM,
                                      1.0, updatedGeo_);
      }

      massInt->SetName("MassIntegrator");
      massInt->SetFeSpace(feFunctions_[formulation_]->GetFeSpace());

      BiLinFormContext *massContext = new BiLinFormContext(massInt, MASS);

      // Check for damping (mass part)
      if (dampingList_[actRegion] == RAYLEIGH)
      {
        if (complexFluidFormulation_ || timeDomainEqFluidFormulation_)
          EXCEPTION("Complex fluid region and Rayleigh damping not allowed!!");

        RaylDampingData &actDamp = regionRaylDamping_[actRegion];
        massContext->SetSecDestMat(DAMPING, actDamp.alpha);
      }

      massContext->SetEntities(actSDList, actSDList);
      massContext->SetFeFunctions(feFunctions_[formulation_], feFunctions_[formulation_]);
      assemble_->AddBiLinearForm(massContext);
      // Important: Add mass-integrator to global list, as we need them later
      // for calculation of postprocessing results
      massInts_[actRegion] = massInt;

      // ====================================================================
      // check for flow (Pierce equation)
      // ====================================================================
      std::string flowId = curRegNode->Get("flowId")->As<std::string>();
      if (flowId != "")
      {
        std::cout << "DO Convective Wave Operator!!!!!!!!!!!!!!!!!" << std::endl
                  << std::endl;
        if (complexFluidFormulation_ || timeDomainEqFluidFormulation_)
          EXCEPTION("Complex fluid and flow currently not allowed");

        // Get result info object for flow
        shared_ptr<ResultInfo> flowInfo = GetResultInfo(MEAN_FLUIDMECH_VELOCITY);

        // Add the region information
        PtrParamNode flowNode = myParam_->Get("flowList")->GetByVal("flow", "name", flowId.c_str());

        bool fullForm = false;
        if (myParam_->Get("flowFormulation")->As<std::string>() == "withDivergence")
          fullForm = true;

        // Read coefficient flow coefficient function for this region and add it to flow functor
        PtrCoefFct regionFlow;
        std::set<UInt> definedDofs;
        bool coefUpdateGeo;
        // we assume that mean flow is pure real
        ReadUserFieldValues(actSDList, flowNode, flowInfo->dofNames, flowInfo->entryType,
                            isComplex_, regionFlow, definedDofs, coefUpdateGeo);
        meanFlowCoef_->AddRegion(actRegion, regionFlow);

        PtrCoefFct divRegionFlow;
        PtrCoefFct divUFactors;
        if (fullForm)
        {
          ReadUserFieldValues(actSDList, flowNode, flowInfo->dofNames, flowInfo->entryType,
                              isComplex_, divRegionFlow, definedDofs, coefUpdateGeo);
          divRegionFlow->SetDerivativeOperation(CoefFunction::VECTOR_DIVERGENCE);
          divMeanFlowCoef_->AddRegion(actRegion, divRegionFlow);
        }

        // now create the integrators
        BaseBDBInt *convectiveStiff = NULL;
        BiLinearForm *convectiveStiffDivU = NULL;
        BiLinearForm *convectiveDamp = NULL;
        BiLinearForm *convectiveDampT = NULL;
        BiLinearForm *convectiveDampDivU = NULL;
        if (dim_ == 2)
        {
          if (isComplex_)
          {

            convectiveDamp = new ABInt<Complex>(new IdentityOperator<FeH1, 2, 1>(),
                                                new ConvectiveOperator<FeH1, 2, 1, Complex>(),
                                                coeffM, 1.0, coefUpdateGeo);
            convectiveDampT = new ABInt<Complex>(new ConvectiveOperator<FeH1, 2, 1, Complex>(),
                                                 new IdentityOperator<FeH1, 2, 1>(),
                                                 coeffM, -1.0, coefUpdateGeo);

            convectiveStiff = new BBInt<Complex>(new ConvectiveOperator<FeH1, 2, 1, Complex>(),
                                                 coeffM, -1.0, coefUpdateGeo);
            if (fullForm)
            {
              divUFactors = CoefFunction::Generate(mp_, Global::COMPLEX,
                                                   CoefXprBinOp(mp_, coeffM, divRegionFlow, CoefXpr::OP_MULT));
              convectiveDampDivU = new BBInt<Complex>(new IdentityOperator<FeH1, 2, 1>(),
                                                      divUFactors, 1.0, coefUpdateGeo);

              convectiveStiffDivU = new ABInt<Complex>(new IdentityOperator<FeH1, 2, 1>(),
                                                       new ConvectiveOperator<FeH1, 2, 1, Complex>(),
                                                       divUFactors, 1.0, coefUpdateGeo);
            }
          }
          else
          {

            convectiveDamp = new ABInt<>(new IdentityOperator<FeH1, 2, 1>(),
                                         new ConvectiveOperator<FeH1, 2, 1>(),
                                         coeffM, 1.0, coefUpdateGeo);
            convectiveStiff = new BBInt<>(new ConvectiveOperator<FeH1, 2, 1>(),
                                          coeffM, -1.0, coefUpdateGeo);

            convectiveDampT = new ABInt<>(new ConvectiveOperator<FeH1, 2, 1>(),
                                          new IdentityOperator<FeH1, 2, 1>(),
                                          coeffM, -1.0, coefUpdateGeo);
            if (fullForm)
            {
              divUFactors = CoefFunction::Generate(mp_, Global::REAL,
                                                   CoefXprBinOp(mp_, divRegionFlow, coeffM, CoefXpr::OP_MULT));
              convectiveDampDivU = new BBInt<>(new IdentityOperator<FeH1, 2, 1>(),
                                               divUFactors, 1.0, coefUpdateGeo);

              convectiveStiffDivU = new ABInt<>(new IdentityOperator<FeH1, 2, 1>(),
                                                new ConvectiveOperator<FeH1, 2, 1>(),
                                                divUFactors, 1.0, coefUpdateGeo);
            }
          }
        }
        else
        {
          if (isComplex_)
          {

            convectiveDamp = new ABInt<Complex>(new IdentityOperator<FeH1, 3, 1>(),
                                                new ConvectiveOperator<FeH1, 3, 1, Complex>(),
                                                coeffM, 1.0, coefUpdateGeo);
            convectiveDampT = new ABInt<Complex>(new ConvectiveOperator<FeH1, 3, 1, Complex>(),
                                                 new IdentityOperator<FeH1, 3, 1>(),
                                                 coeffM, -1.0, coefUpdateGeo);

            convectiveStiff = new BBInt<Complex>(new ConvectiveOperator<FeH1, 3, 1, Complex>(),
                                                 coeffM, -1.0, coefUpdateGeo);
            if (fullForm)
            {
              divUFactors = CoefFunction::Generate(mp_, Global::COMPLEX,
                                                   CoefXprBinOp(mp_, coeffM, divRegionFlow, CoefXpr::OP_MULT));

              convectiveDampDivU = new BBInt<Complex>(new IdentityOperator<FeH1, 3, 1>(),
                                                      divUFactors, 1.0, coefUpdateGeo);

              convectiveStiffDivU = new ABInt<Complex>(new IdentityOperator<FeH1, 3, 1>(),
                                                       new ConvectiveOperator<FeH1, 3, 1, Complex>(),
                                                       divUFactors, 1.0, coefUpdateGeo);
            }
          }
          else
          {

            convectiveDamp = new ABInt<>(new IdentityOperator<FeH1, 3, 1>(),
                                         new ConvectiveOperator<FeH1, 3, 1>(),
                                         coeffM, 1.0, coefUpdateGeo);
            convectiveDampT = new ABInt<>(new ConvectiveOperator<FeH1, 3, 1>(),
                                          new IdentityOperator<FeH1, 3, 1>(),
                                          coeffM, -1.0, coefUpdateGeo);

            convectiveStiff = new BBInt<>(new ConvectiveOperator<FeH1, 3, 1>(),
                                          coeffM, -1.0, coefUpdateGeo);

            if (fullForm)
            {
              divUFactors = CoefFunction::Generate(mp_, Global::REAL,
                                                   CoefXprBinOp(mp_, coeffM, divRegionFlow, CoefXpr::OP_MULT));

              convectiveDampDivU = new BBInt<>(new IdentityOperator<FeH1, 3, 1>(),
                                               divUFactors, 1.0, coefUpdateGeo);

              convectiveStiffDivU = new ABInt<>(new IdentityOperator<FeH1, 3, 1>(),
                                                new ConvectiveOperator<FeH1, 3, 1>(),
                                                divUFactors, 1.0, coefUpdateGeo);
            }
          }
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

        if (fullForm)
        {
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

      } // end flow (convective terms)

      // ====================================================================
      // check for TDEF formulation
      // ====================================================================
      if (timeDomainEqFluidFormulation_ && formulation_ == ACOU_PRESSURE &&
          curRegNode->Get("timeDomainEqFluid")->As<std::string>() == "yes" )
      { // only available for the acou pressure formulation
        // additional terms in the wave equation
        ReadTDEFCoefficients(iRegion);
        DefineTDEFIntegrators(iRegion, actRegion, polyId, integId, actSDList);
      }
    }
  }

  template <UInt DIM>
  void AcousticPDE::DefineTransientPMLInts(shared_ptr<ElemList> eList, std::string id,
                                           RegionIdType actRegion, std::string tempId)
  {

    // define some material coeffunction as above...
    PtrCoefFct dens = materials_[eList->GetRegion()]->GetScalCoefFnc(DENSITY, Global::REAL);
    PtrCoefFct blk = materials_[eList->GetRegion()]->GetScalCoefFnc(ACOU_BULK_MODULUS, Global::REAL);
    PtrCoefFct one = CoefFunction::Generate(mp_, Global::REAL, "1.0");

    PtrCoefFct mechAcouFactor;
    PtrCoefFct corrFactor;
    if (isMechCoupled_ == true && formulation_ != ACOU_PRESSURE)
    {
      // Important: In case of a general / quadratic EV problem, we must
      // ensure to have a "positive definite" matrix, i.e. we are not allowed
      // to multiply all matrices by -1!
      std::string stringFac = (analysistype_ != EIGENFREQUENCY) ? "-1.0" : "1.0";

      mechAcouFactor = CoefFunction::Generate(mp_, Global::REAL,
                                              CoefXprBinOp(mp_, dens, stringFac, CoefXpr::OP_MULT));
    }
    else
    {
      mechAcouFactor = CoefFunction::Generate(mp_, Global::REAL, "1.0");
    }

    // In the case of the TDEF formulation, multiplicatoin by the inverse high-frequency limit of the density is required
    PtrCoefFct TDEFFactor;
    if (timeDomainEqFluidFormulation_){
      TDEFFactor = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, one, dens, CoefXpr::OP_DIV));
      std::cout << "Attention: PML adjacent to TDEF region currently implemented. However, it works for weakly damped materials not accurate (for strongly damped materials)." << std::endl;
    }
    else{
      TDEFFactor = one;
    }

    // calculate the resulting correction factor
    corrFactor = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, TDEFFactor, mechAcouFactor, CoefXpr::OP_MULT));

    PtrCoefFct regionTemp;
    std::set<UInt> definedDofs;
    if (tempId != "")
    {
      // just real valued temperature allowed
      shared_ptr<ResultInfo> tempInfo = GetResultInfo(HEAT_MEAN_TEMPERATURE);
      // Add the region information
      PtrParamNode tempNode = myParam_->Get("temperatureList")->GetByVal("temperature", "name", tempId.c_str());

      ReadUserFieldValues(eList, tempNode, tempInfo->dofNames, tempInfo->entryType,
                          false, regionTemp, definedDofs, updatedGeo_);
    }

    // compute speed of sound
    PtrCoefFct c0;
    ComputeSOS(c0, dens, blk, regionTemp, tempId);

    // coeffc needed in PML
    shared_ptr<CoefFunction> coeffc;
    if (sosAtLaplace_)
      // c0^2
      ComputeSOS_SQR(coeffc, dens, blk, regionTemp, tempId);
    else
    {
      // 1/c0^2
      coeffc = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, dens, blk, CoefXpr::OP_DIV));
    }

    PtrParamNode pmlNode = myParam_->Get("dampingList")->GetByVal("pml", "id", id.c_str());
    shared_ptr<CoefFunction> coeffPMLVec;
    coeffPMLVec.reset(new CoefFunctionPML<Double>(pmlNode, c0, eList, regions_, true));

    // store pml factor
    matCoefs_[PML_DAMP_FACTOR]->AddRegion(eList->GetRegion(), coeffPMLVec);

    shared_ptr<CoefFunctionCompound<Double>> coefA(new CoefFunctionCompound<Double>(mp_));
    shared_ptr<CoefFunctionCompound<Double>> coefB(new CoefFunctionCompound<Double>(mp_));
    shared_ptr<CoefFunctionCompound<Double>> coefAlpha(new CoefFunctionCompound<Double>(mp_));
    shared_ptr<CoefFunctionCompound<Double>> coefBeta(new CoefFunctionCompound<Double>(mp_));

    // --- Set the FE ansatz for the current region ---
    PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region", "name", eList->GetName().c_str());
    std::string polyId = curRegNode->Get("polyId")->As<std::string>();
    std::string integId = curRegNode->Get("integId")->As<std::string>();

    shared_ptr<FeSpace> vecSpace = feFunctions_[ACOU_PMLAUXVEC]->GetFeSpace();
    vecSpace->SetRegionApproximation(eList->GetRegion(), polyId, integId);

    // ===> DEFINE PML DAMPINGFUNCTIONS
    std::map<std::string, PtrCoefFct> vars;
    std::map<std::string, PtrCoefFct> var;
    vars["a"] = coeffPMLVec;
    vars["b"] = coeffc;
    var["a"] = coeffPMLVec;

    StdVector<std::string> matAReal;
    StdVector<std::string> matBReal;
    std::string alpha;
    std::string beta;

    if (DIM == 3)
    {
      if (sosAtLaplace_)
      {
        const std::string Amat[] = {"a_0_R", "0.0", "0.0", "0.0", "a_1_R", "0.0", "0.0", "0.0", "a_2_R"};
        const std::string Bmat[] = {"( (a_0_R - a_1_R - a_2_R) * b_R )",
                                    "0.0", "0.0", "0.0",
                                    "( (a_1_R - a_0_R - a_2_R) * b_R )",
                                    "0.0", "0.0", "0.0",
                                    "( (a_2_R - a_1_R - a_0_R) * b_R )"};
        alpha = "(a_0_R + a_1_R + a_2_R )";
        beta = "( (a_0_R * a_1_R) + (a_0_R * a_2_R) + (a_1_R * a_2_R) )";
        matAReal.Import(Amat, 9);
        matBReal.Import(Bmat, 9);
        coefA->SetTensor(matAReal, DIM, DIM, var);
        coefB->SetTensor(matBReal, DIM, DIM, vars);
        coefAlpha->SetScalar(alpha, var);
        coefBeta->SetScalar(beta, var);
      }
      else
      {
        const std::string Amat[] = {"a_0_R", "0.0", "0.0", "0.0", "a_1_R", "0.0", "0.0", "0.0", "a_2_R"};
        const std::string Bmat[] = {"( a_0_R - a_1_R - a_2_R ) ",
                                    "0.0", "0.0", "0.0",
                                    "( a_1_R - a_0_R - a_2_R )",
                                    "0.0", "0.0", "0.0",
                                    "( a_2_R - a_1_R - a_0_R) "};
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
    else
    {
      // dim = 2
      if (sosAtLaplace_)
      {
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
      else
      {
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

    /// now lets define some integrators
    // the vectorial auxiliary variable is called U...
    // in case of mechacoucoupling, we need to multiply all integrators which couple into acouPDE with -density or density
    // if eigenfrequency

    PtrCoefFct mAcouCorrect_CoefAlpha =
        CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, corrFactor, coefAlpha, CoefXpr::OP_MULT));
    PtrCoefFct mAcouCorrect_CoefBeta =
        CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, corrFactor, coefBeta, CoefXpr::OP_MULT));

    BaseBDBInt *dampdPdt = new BBInt<>(new IdentityOperator<FeH1, DIM>(), mAcouCorrect_CoefAlpha, 1.0, updatedGeo_);
    BaseBDBInt *dampP = new BBInt<>(new IdentityOperator<FeH1, DIM>(), mAcouCorrect_CoefBeta, 1.0, updatedGeo_);
    // this is already integrated by parts...
    BaseBDBInt *divU = new ABInt<>(new GradientOperator<FeH1, DIM>(), new IdentityOperator<FeH1, DIM, DIM>(),
                                   corrFactor, 1.0, updatedGeo_);

    BaseBDBInt *dUdt = new BBInt<>(new IdentityOperator<FeH1, DIM, DIM>(), one, 1.0, updatedGeo_);
    BaseBDBInt *AU = new BDBInt<>(new IdentityOperator<FeH1, DIM, DIM>(), coefA, 1.0, updatedGeo_);
    BaseBDBInt *BgradP = new ADBInt<>(new IdentityOperator<FeH1, DIM, DIM>(),
                                      new GradientOperator<FeH1, DIM>(), coefB, 1.0, updatedGeo_);

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

    if (!this->isAPML_ && DIM == 3)
    {
      shared_ptr<FeSpace> scalSpace = feFunctions_[ACOU_PMLAUXSCALAR]->GetFeSpace();
      scalSpace->SetRegionApproximation(eList->GetRegion(), polyId, integId);
      feFunctions_[ACOU_PMLAUXSCALAR]->AddEntityList(eList);

      shared_ptr<CoefFunctionCompound<Double>> coefGamma(new CoefFunctionCompound<Double>(mp_));
      shared_ptr<CoefFunctionCompound<Double>> coefC(new CoefFunctionCompound<Double>(mp_));

      StdVector<std::string> matCReal;

      if (sosAtLaplace_)
      {
        std::string gamma = "( a_0_R * a_1_R * a_2_R ) ";
        const std::string Cmat[] = {" (a_1_R * a_2_R * b_R) ",
                                    "0.0", "0.0", "0.0",
                                    " (a_0_R * a_2_R * b_R) ",
                                    "0.0", "0.0", "0.0",
                                    " (a_0_R * a_1_R * b_R) "};
        matCReal.Import(Cmat, 9);
        coefC->SetTensor(matCReal, 3, 3, vars);
        coefGamma->SetScalar(gamma, var);
      }
      else
      {
        std::string gamma = "( (a_0_R * a_1_R * a_2_R) * b_R )";
        const std::string Cmat[] = {" (a_1_R * a_2_R) ", "0.0", "0.0", "0.0", " (a_0_R * a_2_R) ", "0.0", "0.0", "0.0", " (a_0_R * a_1_R) "};
        matCReal.Import(Cmat, 9);
        coefC->SetTensor(matCReal, 3, 3, var);
        coefGamma->SetScalar(gamma, vars);
      }

      PtrCoefFct mAcouCorrect_CoefGamma =
          CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, corrFactor, coefGamma, CoefXpr::OP_MULT));

      // the scalar auxiliary variable is called Nu...
      BaseBDBInt *dampNu = new BBInt<>(new IdentityOperator<FeH1, DIM>(), mAcouCorrect_CoefGamma, 1.0, updatedGeo_);
      BaseBDBInt *CgradNu = new ADBInt<>(new IdentityOperator<FeH1, DIM, DIM>(),
                                         new GradientOperator<FeH1, DIM>(), coefC, -1.0, updatedGeo_);
      BaseBDBInt *dNudt = new BBInt<>(new IdentityOperator<FeH1, DIM>(), one, 1.0, updatedGeo_);
      BaseBDBInt *P = new BBInt<>(new IdentityOperator<FeH1, DIM>(), one, -1.0, updatedGeo_);

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

  void AcousticPDE::DefineTDEFIntegrators(UInt iRegion, RegionIdType actRegion, string polyId, string integId, shared_ptr<ElemList> actSDList)
  {
    std::cout << "Establishing the TDEF integrators:" << std::endl
              << std::endl;
    PtrCoefFct constOne = CoefFunction::Generate(mp_, Global::REAL, "1.0");

    // ====================================================================
    // Loop over real poles of inverse blk modulus (compressibility)
    for (UInt ii = 0; ii < fncAC_.GetSize(); ii++)
    {
      std::cout << "re Poles compressibility: " << ii + 1 << " from " << fncAC_.GetSize() << std::endl;
      // ====================================================================
      // K_PPHI1 (TDEF): stiffness integrator, TDEF (A_j^C term)
      // \int_{Omega_1} A_j^C p^\prime \phi_j^C d\Omega
      // ====================================================================
      shared_ptr<FeSpace> spacePhiC = feFunctions_[(SolutionType)(ACOU_TDEF_PHI_C_1 + ii)]->GetFeSpace();
      spacePhiC->SetRegionApproximation(actRegion, polyId, integId);

      BiLinearForm *stiffIntTDEFPPHI1 = NULL;

      if (dim_ == 2)
      {
        stiffIntTDEFPPHI1 = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), fncAC_[ii], 1.0, updatedGeo_);
      }
      else
      {
        stiffIntTDEFPPHI1 = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), fncAC_[ii], 1.0, updatedGeo_);
      }

      stiffIntTDEFPPHI1->SetName("AcousticStiffIntTDEFPPHI1_" + std::to_string(ii));

      BiLinFormContext *stiffIntTDEFPPHI1Context = NULL;
      stiffIntTDEFPPHI1Context = new BiLinFormContext(stiffIntTDEFPPHI1, STIFFNESS);

      stiffIntTDEFPPHI1Context->SetEntities(actSDList, actSDList);
      stiffIntTDEFPPHI1Context->SetFeFunctions(feFunctions_[formulation_], feFunctions_[(SolutionType)(ACOU_TDEF_PHI_C_1 + ii)]);
      assemble_->AddBiLinearForm(stiffIntTDEFPPHI1Context);

      feFunctions_[(SolutionType)(ACOU_TDEF_PHI_C_1 + ii)]->AddEntityList(actSDList);
    } // end loop stiffIntTDEFPPHI1

    // ====================================================================
    // Loop over complex poles  of inverse blk modulus (compressibility)
    for (UInt ii = 0; ii < fncBC_.GetSize(); ii++)
    {
      std::cout << "compl. poles compressibility: " << ii + 1 << " from " << fncBC_.GetSize() << std::endl;

      // ====================================================================
      // D_PPSI1 (TDEF): damping integrator, TDEF (\frac{B_k^C,\gamma_k^C} term)
      // \int_{Omega_1} \frac{B_k^C,\gamma_k^C} p^\prime \dot{\psi}_k^C d\Omega
      // ====================================================================
      shared_ptr<FeSpace> spacePsiC = feFunctions_[(SolutionType)(ACOU_TDEF_PSI_C_1 + ii)]->GetFeSpace();
      spacePsiC->SetRegionApproximation(actRegion, polyId, integId);

      BiLinearForm *dampIntTDEFPPSI1 = NULL;
      PtrCoefFct fncBCgammaC;
      fncBCgammaC = CoefFunction::Generate(mp_, Global::REAL,
                                           CoefXprBinOp(mp_, fncBC_[ii], fncGammaC_[ii], CoefXpr::OP_DIV));
      if (dim_ == 2)
      {
        dampIntTDEFPPSI1 = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), fncBCgammaC, 1.0, updatedGeo_);
      }
      else
      {
        dampIntTDEFPPSI1 = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), fncBCgammaC, 1.0, updatedGeo_);
      }

      dampIntTDEFPPSI1->SetName("AcousticDampIntTDEFPPSI1_" + std::to_string(ii));

      BiLinFormContext *dampIntTDEFPPSI1Context = NULL;
      dampIntTDEFPPSI1Context = new BiLinFormContext(dampIntTDEFPPSI1, DAMPING);

      dampIntTDEFPPSI1Context->SetEntities(actSDList, actSDList);
      dampIntTDEFPPSI1Context->SetFeFunctions(feFunctions_[formulation_], feFunctions_[(SolutionType)(ACOU_TDEF_PSI_C_1 + ii)]);
      assemble_->AddBiLinearForm(dampIntTDEFPPSI1Context);

      feFunctions_[(SolutionType)(ACOU_TDEF_PSI_C_1 + ii)]->AddEntityList(actSDList);

      // ====================================================================
      // K_PPSI1 (TDEF): stiffness integrator, TDEF (D_k^C term)
      // \int_{Omega_1} D_k^C p^\prime \psi_k^C d\Omega
      // ====================================================================

      BiLinearForm *stiffIntTDEFPPSI1 = NULL;

      PtrCoefFct fncQuotC;
      fncQuotC = CoefFunction::Generate(mp_, Global::REAL,
                                        CoefXprBinOp(mp_, fncBetaC_[ii], fncGammaC_[ii], CoefXpr::OP_DIV));
      PtrCoefFct fncTermC;
      fncTermC = CoefFunction::Generate(mp_, Global::REAL,
                                        CoefXprBinOp(mp_, fncBC_[ii], fncQuotC, CoefXpr::OP_MULT));
      PtrCoefFct fncDC;
      fncDC = CoefFunction::Generate(mp_, Global::REAL,
                                     CoefXprBinOp(mp_, fncCC_[ii], fncTermC, CoefXpr::OP_ADD));

      if (dim_ == 2)
      {
        stiffIntTDEFPPSI1 = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), fncDC, 1.0, updatedGeo_);
      }
      else
      {
        stiffIntTDEFPPSI1 = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), fncDC, 1.0, updatedGeo_);
      }

      stiffIntTDEFPPSI1->SetName("AcousticStiffIntTDEFPPSI1_" + std::to_string(ii));

      BiLinFormContext *stiffIntTDEFPPSI1Context = NULL;
      stiffIntTDEFPPSI1Context = new BiLinFormContext(stiffIntTDEFPPSI1, STIFFNESS);

      stiffIntTDEFPPSI1Context->SetEntities(actSDList, actSDList);
      stiffIntTDEFPPSI1Context->SetFeFunctions(feFunctions_[formulation_], feFunctions_[(SolutionType)(ACOU_TDEF_PSI_C_1 + ii)]);
      assemble_->AddBiLinearForm(stiffIntTDEFPPSI1Context);
    } // end loop stiffIntTDEFPPSI1 and dampIntTDEFPPSI1

    // ====================================================================
    // Loop over real poles of inverse density (specific volume)
    for (UInt ii = 0; ii < fncAV_.GetSize(); ii++)
    {
      std::cout << "re poles density: " << ii + 1 << " from " << fncAV_.GetSize() << std::endl;
      // ====================================================================
      // K_PPHI2 (TDEF): stiffness integrator, TDEF (A_l^V term)
      // \int_{Omega_1} A_l^V \nabla p^\prime \nabla \phi_l^V d\Omega
      // ====================================================================
      shared_ptr<FeSpace> spacePhiV = feFunctions_[(SolutionType)(ACOU_TDEF_PHI_V_1 + ii)]->GetFeSpace();
      spacePhiV->SetRegionApproximation(actRegion, polyId, integId);

      BiLinearForm *stiffIntTDEFPPHI2 = NULL;

      if (dim_ == 2)
      {
        stiffIntTDEFPPHI2 = new BBInt<Double>(new GradientOperator<FeH1, 2>(), fncAV_[ii], 1.0, updatedGeo_);
      }
      else
      {
        stiffIntTDEFPPHI2 = new BBInt<Double>(new GradientOperator<FeH1, 3>(), fncAV_[ii], 1.0, updatedGeo_);
      }

      stiffIntTDEFPPHI2->SetName("AcousticStiffIntTDEFPPHI2_" + std::to_string(ii));

      BiLinFormContext *stiffIntTDEFPPHI2Context = NULL;
      stiffIntTDEFPPHI2Context = new BiLinFormContext(stiffIntTDEFPPHI2, STIFFNESS);

      stiffIntTDEFPPHI2Context->SetEntities(actSDList, actSDList);
      stiffIntTDEFPPHI2Context->SetFeFunctions(feFunctions_[formulation_], feFunctions_[(SolutionType)(ACOU_TDEF_PHI_V_1 + ii)]);
      assemble_->AddBiLinearForm(stiffIntTDEFPPHI2Context);

      feFunctions_[(SolutionType)(ACOU_TDEF_PHI_V_1 + ii)]->AddEntityList(actSDList);
    } // end loop stiffIntTDEFPPHI2

    // ====================================================================
    // Loop over complex poles of inverse density (specific volume)
    for (UInt ii = 0; ii < fncBV_.GetSize(); ii++)
    {
      std::cout << "compl. poles density: " << ii + 1 << " from " << fncBV_.GetSize() << std::endl;
      // ====================================================================
      // D_PPSI2 (TDEF): damping integrator, TDEF (\frac{B_k^V,\gamma_k^V} term)
      // \int_{Omega_1} \frac{B_k^V,\gamma_k^V} p^\prime \dot{\psi}_m^C d\Omega
      // ====================================================================
      shared_ptr<FeSpace> spacePsiV = feFunctions_[(SolutionType)(ACOU_TDEF_PSI_V_1 + ii)]->GetFeSpace();
      spacePsiV->SetRegionApproximation(actRegion, polyId, integId);

      BiLinearForm *dampIntTDEFPPSI2 = NULL;
      PtrCoefFct fncBVgammaV;
      fncBVgammaV = CoefFunction::Generate(mp_, Global::REAL,
                                           CoefXprBinOp(mp_, fncBV_[ii], fncGammaV_[ii], CoefXpr::OP_DIV));
      if (dim_ == 2)
      {
        dampIntTDEFPPSI2 = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), fncBVgammaV, 1.0, updatedGeo_);
      }
      else
      {
        dampIntTDEFPPSI2 = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), fncBVgammaV, 1.0, updatedGeo_);
      }

      dampIntTDEFPPSI2->SetName("AcousticDampIntTDEFPPSI2_" + std::to_string(ii));

      BiLinFormContext *dampIntTDEFPPSI2Context = NULL;
      dampIntTDEFPPSI2Context = new BiLinFormContext(dampIntTDEFPPSI2, DAMPING);

      dampIntTDEFPPSI2Context->SetEntities(actSDList, actSDList);
      dampIntTDEFPPSI2Context->SetFeFunctions(feFunctions_[formulation_], feFunctions_[(SolutionType)(ACOU_TDEF_PSI_V_1 + ii)]);
      assemble_->AddBiLinearForm(dampIntTDEFPPSI2Context);

      feFunctions_[(SolutionType)(ACOU_TDEF_PSI_V_1 + ii)]->AddEntityList(actSDList);

      // ====================================================================
      // K_PPSI2 (TDEF): stiffness integrator, TDEF (D_k^V term)
      // \int_{Omega_1} D_m^V \nabla p^\prime \nabla \psi_m^V d\Omega
      // ====================================================================
      BiLinearForm *stiffIntTDEFPPSI2 = NULL;

      PtrCoefFct fncQuotV;
      fncQuotV = CoefFunction::Generate(mp_, Global::REAL,
                                        CoefXprBinOp(mp_, fncBetaV_[ii], fncGammaV_[ii], CoefXpr::OP_DIV));
      PtrCoefFct fncTermV;
      fncTermV = CoefFunction::Generate(mp_, Global::REAL,
                                        CoefXprBinOp(mp_, fncBV_[ii], fncQuotV, CoefXpr::OP_MULT));
      PtrCoefFct fncDV;
      fncDV = CoefFunction::Generate(mp_, Global::REAL,
                                     CoefXprBinOp(mp_, fncCV_[ii], fncTermV, CoefXpr::OP_ADD));
      if (dim_ == 2)
      {
        stiffIntTDEFPPSI2 = new BBInt<Double>(new GradientOperator<FeH1, 2>(), fncDV, 1.0, updatedGeo_);
      }
      else
      {
        stiffIntTDEFPPSI2 = new BBInt<Double>(new GradientOperator<FeH1, 3>(), fncDV, 1.0, updatedGeo_);
      }

      stiffIntTDEFPPSI2->SetName("AcousticStiffIntTDEFPPSI2_" + std::to_string(ii));

      BiLinFormContext *stiffIntTDEFPPSI2Context = NULL;
      stiffIntTDEFPPSI2Context = new BiLinFormContext(stiffIntTDEFPPSI2, STIFFNESS);

      stiffIntTDEFPPSI2Context->SetEntities(actSDList, actSDList);
      stiffIntTDEFPPSI2Context->SetFeFunctions(feFunctions_[formulation_], feFunctions_[(SolutionType)(ACOU_TDEF_PSI_V_1 + ii)]);
      assemble_->AddBiLinearForm(stiffIntTDEFPPSI2Context);

    } // end loop stiffIntTDEFPPSI2 and dampIntTDEFPPSI2

    // ######################################################
    // ADE section
    // ######################################################

    // ====================================================================
    // Loop over real poles of inverse compression modulus (compressibility)
    for (UInt ii = 0; ii < fncAlphaC_.GetSize(); ii++)
    {
      std::cout << "re. ADE compress.: " << ii + 1 << " from " << fncAlphaC_.GetSize() << std::endl;
      // ====================================================================
      // D_PHI1PHI1 (TDEF): damping integrator, TDEF
      // \int_{Omega_1} \phi_j^{C,\prime} \prime{\phi}_j^C d\Omega
      // ====================================================================
      BiLinearForm *dampIntTDEFPHI1PHI1 = NULL;

      if (dim_ == 2)
      {
        dampIntTDEFPHI1PHI1 = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), constOne, 1.0, updatedGeo_);
      }
      else
      {
        dampIntTDEFPHI1PHI1 = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), constOne, 1.0, updatedGeo_);
      }

      dampIntTDEFPHI1PHI1->SetName("AcousticDampIntTDEFPHI1PHI1_" + std::to_string(ii));

      BiLinFormContext *dampIntTDEFPHI1PHI1Context = NULL;
      dampIntTDEFPHI1PHI1Context = new BiLinFormContext(dampIntTDEFPHI1PHI1, DAMPING);

      dampIntTDEFPHI1PHI1Context->SetEntities(actSDList, actSDList);
      dampIntTDEFPHI1PHI1Context->SetFeFunctions(feFunctions_[(SolutionType)(ACOU_TDEF_PHI_C_1 + ii)], feFunctions_[(SolutionType)(ACOU_TDEF_PHI_C_1 + ii)]);
      assemble_->AddBiLinearForm(dampIntTDEFPHI1PHI1Context);

      // ====================================================================
      // K_PHI1PHI1 (TDEF): stiffness integrator, TDEF (\alpha_j^C term)
      // \int_{Omega_1} \alpha_j^C \phi_j^{C,\prime} \phi_j^C d\Omega
      // ====================================================================
      BiLinearForm *stiffIntTDEFPHI1PHI1 = NULL;

      if (dim_ == 2)
      {
        stiffIntTDEFPHI1PHI1 = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), fncAlphaC_[ii], 1.0, updatedGeo_);
      }
      else
      {
        stiffIntTDEFPHI1PHI1 = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), fncAlphaC_[ii], 1.0, updatedGeo_);
      }

      stiffIntTDEFPHI1PHI1->SetName("AcousticStiffIntTDEFPHI1PHI1_" + std::to_string(ii));

      BiLinFormContext *stiffIntTDEFPHI1PHI1Context = NULL;
      stiffIntTDEFPHI1PHI1Context = new BiLinFormContext(stiffIntTDEFPHI1PHI1, STIFFNESS);

      stiffIntTDEFPHI1PHI1Context->SetEntities(actSDList, actSDList);
      stiffIntTDEFPHI1PHI1Context->SetFeFunctions(feFunctions_[(SolutionType)(ACOU_TDEF_PHI_C_1 + ii)], feFunctions_[(SolutionType)(ACOU_TDEF_PHI_C_1 + ii)]);
      assemble_->AddBiLinearForm(stiffIntTDEFPHI1PHI1Context);

      // ====================================================================
      // M_PHI1P (TDEF): mass integrator, TDEF
      // -\int_{Omega_1} \phi_j^{C,\prime} \ddot{p}_a d\Omega
      // ====================================================================
      BiLinearForm *massIntTDEFPHI1P = NULL;

      if (dim_ == 2)
      {
        massIntTDEFPHI1P = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), constOne, -1.0, updatedGeo_);
      }
      else
      {
        massIntTDEFPHI1P = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), constOne, -1.0, updatedGeo_);
      }

      massIntTDEFPHI1P->SetName("AcousticMassIntTDEFPHI1P_" + std::to_string(ii));

      BiLinFormContext *massIntTDEFPHI1PContext = NULL;
      massIntTDEFPHI1PContext = new BiLinFormContext(massIntTDEFPHI1P, MASS);

      massIntTDEFPHI1PContext->SetEntities(actSDList, actSDList);
      massIntTDEFPHI1PContext->SetFeFunctions(feFunctions_[(SolutionType)(ACOU_TDEF_PHI_C_1 + ii)], feFunctions_[formulation_]);
      assemble_->AddBiLinearForm(massIntTDEFPHI1PContext);
    } // end loop dampIntTDEFPHI1PHI1, stiffIntTDEFPHI1PHI1 and massInteTDEFPHI1P

    // ====================================================================
    // Loop over complex poles of inverse compression modulus (compressibility)
    for (UInt ii = 0; ii < fncGammaC_.GetSize(); ii++)
    {
      std::cout << "compl. ADE compress.: " << ii + 1 << " from " << fncGammaC_.GetSize() << std::endl;
      // ====================================================================
      // M_PSI1PSI1 (TDEF): mass integrator, TDEF (\frac{1,\gamma_k^C} term)
      // \int_{Omega_1} \frac{1,\gamma_k^C} \psi_k^{C,\prime} \prime{\psi}_k^C d\Omega
      // ====================================================================
      BiLinearForm *massIntTDEFPSI1PSI1 = NULL;

      PtrCoefFct coefOneOverGammaC;
      coefOneOverGammaC = CoefFunction::Generate(mp_, Global::REAL,
                                                 CoefXprBinOp(mp_, constOne, fncGammaC_[ii], CoefXpr::OP_DIV));
      if (dim_ == 2)
      {
        massIntTDEFPSI1PSI1 = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), coefOneOverGammaC, 1.0, updatedGeo_);
      }
      else
      {
        massIntTDEFPSI1PSI1 = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), coefOneOverGammaC, 1.0, updatedGeo_);
      }

      massIntTDEFPSI1PSI1->SetName("AcousticMassIntTDEFPSI1PSI1_" + std::to_string(ii));

      BiLinFormContext *massIntTDEFPSI1PSI1Context = NULL;
      massIntTDEFPSI1PSI1Context = new BiLinFormContext(massIntTDEFPSI1PSI1, MASS);

      massIntTDEFPSI1PSI1Context->SetEntities(actSDList, actSDList);
      massIntTDEFPSI1PSI1Context->SetFeFunctions(feFunctions_[(SolutionType)(ACOU_TDEF_PSI_C_1 + ii)], feFunctions_[(SolutionType)(ACOU_TDEF_PSI_C_1 + ii)]);
      assemble_->AddBiLinearForm(massIntTDEFPSI1PSI1Context);

      // ====================================================================
      // D_PSI1PSI1 (TDEF): damping integrator, TDEF (\frac{2 \beta_k^C,\gamma_k^C} term)
      // \int_{Omega_1} \frac{2 \beta_k^C,\gamma_k^C} \psi_j^{C,\prime} \dot{\psi}_j^C d\Omega
      // ====================================================================
      BiLinearForm *dampIntTDEFPSI1PSI1 = NULL;

      PtrCoefFct coefTwoBetaC;
      coefTwoBetaC = CoefFunction::Generate(mp_, Global::REAL,
                                            CoefXprBinOp(mp_, CoefFunction::Generate(mp_, Global::REAL, "2.0"), fncBetaC_[ii], CoefXpr::OP_MULT));
      PtrCoefFct coefTwoBetaCOverGammaC;
      coefTwoBetaCOverGammaC = CoefFunction::Generate(mp_, Global::REAL,
                                                      CoefXprBinOp(mp_, coefTwoBetaC, fncGammaC_[ii], CoefXpr::OP_DIV));
      if (dim_ == 2)
      {
        dampIntTDEFPSI1PSI1 = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), coefTwoBetaCOverGammaC, 1.0, updatedGeo_);
      }
      else
      {
        dampIntTDEFPSI1PSI1 = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), coefTwoBetaCOverGammaC, 1.0, updatedGeo_);
      }

      dampIntTDEFPSI1PSI1->SetName("AcousticDampIntTDEFPSI1PSI1_" + std::to_string(ii));

      BiLinFormContext *dampIntTDEFPSI1PSI1Context = NULL;
      dampIntTDEFPSI1PSI1Context = new BiLinFormContext(dampIntTDEFPSI1PSI1, DAMPING);

      dampIntTDEFPSI1PSI1Context->SetEntities(actSDList, actSDList);
      dampIntTDEFPSI1PSI1Context->SetFeFunctions(feFunctions_[(SolutionType)(ACOU_TDEF_PSI_C_1 + ii)], feFunctions_[(SolutionType)(ACOU_TDEF_PSI_C_1 + ii)]);
      assemble_->AddBiLinearForm(dampIntTDEFPSI1PSI1Context);

      // ====================================================================
      // K_PSI1PSI1 (TDEF): stiffness integrator, TDEF (\delta_j^C term)
      // \int_{Omega_1} \delta_j^C \psi_j^{C,\prime} \psi_j^C d\Omega
      // ====================================================================
      BiLinearForm *stiffIntTDEFPSI1PSI1 = NULL;
      PtrCoefFct fncQuotC;
      fncQuotC = CoefFunction::Generate(mp_, Global::REAL,
                                        CoefXprBinOp(mp_, fncBetaC_[ii], fncGammaC_[ii], CoefXpr::OP_DIV));
      PtrCoefFct fncTermC;
      fncTermC = CoefFunction::Generate(mp_, Global::REAL,
                                        CoefXprBinOp(mp_, fncBetaC_[ii], fncQuotC, CoefXpr::OP_MULT));
      PtrCoefFct fncDeltaC;
      fncDeltaC = CoefFunction::Generate(mp_, Global::REAL,
                                         CoefXprBinOp(mp_, fncGammaC_[ii], fncTermC, CoefXpr::OP_ADD));
      if (dim_ == 2)
      {
        stiffIntTDEFPSI1PSI1 = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), fncDeltaC, 1.0, updatedGeo_);
      }
      else
      {
        stiffIntTDEFPSI1PSI1 = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), fncDeltaC, 1.0, updatedGeo_);
      }

      stiffIntTDEFPSI1PSI1->SetName("AcousticStiffIntTDEFPSI1PSI1_" + std::to_string(ii));

      BiLinFormContext *stiffIntTDEFPSI1PSI1Context = NULL;
      stiffIntTDEFPSI1PSI1Context = new BiLinFormContext(stiffIntTDEFPSI1PSI1, STIFFNESS);

      stiffIntTDEFPSI1PSI1Context->SetEntities(actSDList, actSDList);
      stiffIntTDEFPSI1PSI1Context->SetFeFunctions(feFunctions_[(SolutionType)(ACOU_TDEF_PSI_C_1 + ii)], feFunctions_[(SolutionType)(ACOU_TDEF_PSI_C_1 + ii)]);
      assemble_->AddBiLinearForm(stiffIntTDEFPSI1PSI1Context);

      // ====================================================================
      // M_PSI1P (TDEF): mass integrator, TDEF
      // -\int_{Omega_1} \psi_j^{C,\prime} \ddot{p}_a d\Omega
      // ====================================================================
      BiLinearForm *massIntTDEFPSI1P = NULL;

      if (dim_ == 2)
      {
        massIntTDEFPSI1P = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), constOne, -1.0, updatedGeo_);
      }
      else
      {
        massIntTDEFPSI1P = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), constOne, -1.0, updatedGeo_);
      }

      massIntTDEFPSI1P->SetName("AcousticMassIntTDEFPSI1P_" + std::to_string(ii));

      BiLinFormContext *massIntTDEFPSI1PContext = NULL;
      massIntTDEFPSI1PContext = new BiLinFormContext(massIntTDEFPSI1P, MASS);

      massIntTDEFPSI1PContext->SetEntities(actSDList, actSDList);
      massIntTDEFPSI1PContext->SetFeFunctions(feFunctions_[(SolutionType)(ACOU_TDEF_PSI_C_1 + ii)], feFunctions_[formulation_]);
      assemble_->AddBiLinearForm(massIntTDEFPSI1PContext);
    } // end loop massIntTDEFPSI1PSI1, dampIntTDEFPSI1PSI1, stiffIntTDEFPSI1PSI1 and massInteTDEFPSI1P

    // ====================================================================
    // Loop over real poles of inverse density (specific volume)
    for (UInt ii = 0; ii < fncAlphaV_.GetSize(); ii++)
    {
      std::cout << "re. ADE density.: " << ii + 1 << " from " << fncAlphaV_.GetSize() << std::endl
                << std::endl;
      // ====================================================================
      // D_PHI2PHI2 (TDEF): damping integrator, TDEF
      // \int_{Omega_1} \phi_l^{V,\prime} \prime{\phi}_l^V d\Omega
      // ====================================================================
      BiLinearForm *dampIntTDEFPHI2PHI2 = NULL;

      if (dim_ == 2)
      {
        dampIntTDEFPHI2PHI2 = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), constOne, 1.0, updatedGeo_);
      }
      else
      {
        dampIntTDEFPHI2PHI2 = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), constOne, 1.0, updatedGeo_);
      }

      dampIntTDEFPHI2PHI2->SetName("AcousticDampIntTDEFPHI2PHI2_" + std::to_string(ii));

      BiLinFormContext *dampIntTDEFPHI2PHI2Context = NULL;
      dampIntTDEFPHI2PHI2Context = new BiLinFormContext(dampIntTDEFPHI2PHI2, DAMPING);

      dampIntTDEFPHI2PHI2Context->SetEntities(actSDList, actSDList);
      dampIntTDEFPHI2PHI2Context->SetFeFunctions(feFunctions_[(SolutionType)(ACOU_TDEF_PHI_V_1 + ii)], feFunctions_[(SolutionType)(ACOU_TDEF_PHI_V_1 + ii)]);
      assemble_->AddBiLinearForm(dampIntTDEFPHI2PHI2Context);

      // ====================================================================
      // K_PHI2PHI2 (TDEF): stiffness integrator, TDEF (\alpha_j^C term)
      // \int_{Omega_1} \alpha_l^V \phi_l^{V,\prime} \phi_l^V d\Omega
      // ====================================================================
      BiLinearForm *stiffIntTDEFPHI2PHI2 = NULL;

      if (dim_ == 2)
      {
        stiffIntTDEFPHI2PHI2 = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), fncAlphaV_[ii], 1.0, updatedGeo_);
      }
      else
      {
        stiffIntTDEFPHI2PHI2 = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), fncAlphaV_[ii], 1.0, updatedGeo_);
      }

      stiffIntTDEFPHI2PHI2->SetName("AcousticStiffIntTDEFPHI2PHI2_" + std::to_string(ii));

      BiLinFormContext *stiffIntTDEFPHI2PHI2Context = NULL;
      stiffIntTDEFPHI2PHI2Context = new BiLinFormContext(stiffIntTDEFPHI2PHI2, STIFFNESS);

      stiffIntTDEFPHI2PHI2Context->SetEntities(actSDList, actSDList);
      stiffIntTDEFPHI2PHI2Context->SetFeFunctions(feFunctions_[(SolutionType)(ACOU_TDEF_PHI_V_1 + ii)], feFunctions_[(SolutionType)(ACOU_TDEF_PHI_V_1 + ii)]);
      assemble_->AddBiLinearForm(stiffIntTDEFPHI2PHI2Context);

      // ====================================================================
      // K_PHI2P (TDEF): stiffness integrator, TDEF
      // -\int_{Omega_1} \phi_l^{V,\prime} \ddot{p}_a d\Omega
      // ====================================================================
      BiLinearForm *stiffIntTDEFPHI2P = NULL;

      if (dim_ == 2)
      {
        stiffIntTDEFPHI2P = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), constOne, -1.0, updatedGeo_);
      }
      else
      {
        stiffIntTDEFPHI2P = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), constOne, -1.0, updatedGeo_);
      }

      stiffIntTDEFPHI2P->SetName("AcousticStiffIntTDEFPHI2P_" + std::to_string(ii));

      BiLinFormContext *stiffIntTDEFPHI2PContext = NULL;
      stiffIntTDEFPHI2PContext = new BiLinFormContext(stiffIntTDEFPHI2P, STIFFNESS);

      stiffIntTDEFPHI2PContext->SetEntities(actSDList, actSDList);
      stiffIntTDEFPHI2PContext->SetFeFunctions(feFunctions_[(SolutionType)(ACOU_TDEF_PHI_V_1 + ii)], feFunctions_[formulation_]);
      assemble_->AddBiLinearForm(stiffIntTDEFPHI2PContext);
    } // end loop dampIntTDEFPHI2PHI2, stiffIntTDEFPHI2PHI2 and stiffIntTDEFPHI2P

    // ====================================================================
    // Loop over complex poles of inverse density (specific volume)
    for (UInt ii = 0; ii < fncGammaV_.GetSize(); ii++)
    {
      std::cout << "compl. ADE density.: " << ii + 1 << " from " << fncGammaV_.GetSize() << std::endl;
      // ====================================================================
      // M_PSI2PSI2 (TDEF): mass integrator, TDEF (\frac{1,\gamma_k^V} term)
      // \int_{Omega_1} \frac{1,\gamma_m^V} \psi_m^{V,\prime} \prime{\psi}_m^V d\Omega
      // ====================================================================
      BiLinearForm *massIntTDEFPSI2PSI2 = NULL;

      PtrCoefFct coefOneOverGammaV;
      coefOneOverGammaV = CoefFunction::Generate(mp_, Global::REAL,
                                                 CoefXprBinOp(mp_, constOne, fncGammaV_[ii], CoefXpr::OP_DIV));
      if (dim_ == 2)
      {
        massIntTDEFPSI2PSI2 = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), coefOneOverGammaV, 1.0, updatedGeo_);
      }
      else
      {
        massIntTDEFPSI2PSI2 = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), coefOneOverGammaV, 1.0, updatedGeo_);
      }

      massIntTDEFPSI2PSI2->SetName("AcousticMassIntTDEFPSI2PSI2_" + std::to_string(ii));

      BiLinFormContext *massIntTDEFPSI2PSI2Context = NULL;
      massIntTDEFPSI2PSI2Context = new BiLinFormContext(massIntTDEFPSI2PSI2, MASS);

      massIntTDEFPSI2PSI2Context->SetEntities(actSDList, actSDList);
      massIntTDEFPSI2PSI2Context->SetFeFunctions(feFunctions_[(SolutionType)(ACOU_TDEF_PSI_V_1 + ii)], feFunctions_[(SolutionType)(ACOU_TDEF_PSI_V_1 + ii)]);
      assemble_->AddBiLinearForm(massIntTDEFPSI2PSI2Context);

      // ====================================================================
      // D_PSI2PSI2 (TDEF): damping integrator, TDEF (\frac{2 \beta_m^V,\gamma_m^V} term)
      // \int_{Omega_1} \frac{2 \beta_m^V,\gamma_m^V} \psi_m^{V,\prime} \dot{\psi}_m^V d\Omega
      // ====================================================================
      BiLinearForm *dampIntTDEFPSI2PSI2 = NULL;

      PtrCoefFct coefTwoBetaV;
      coefTwoBetaV = CoefFunction::Generate(mp_, Global::REAL,
                                            CoefXprBinOp(mp_, CoefFunction::Generate(mp_, Global::REAL, "2.0"), fncBetaV_[ii], CoefXpr::OP_MULT));
      PtrCoefFct coefTwoBetaVOverGammaV;
      coefTwoBetaVOverGammaV = CoefFunction::Generate(mp_, Global::REAL,
                                                      CoefXprBinOp(mp_, coefTwoBetaV, fncGammaV_[ii], CoefXpr::OP_DIV));
      if (dim_ == 2)
      {
        dampIntTDEFPSI2PSI2 = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), coefTwoBetaVOverGammaV, 1.0, updatedGeo_);
      }
      else
      {
        dampIntTDEFPSI2PSI2 = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), coefTwoBetaVOverGammaV, 1.0, updatedGeo_);
      }

      dampIntTDEFPSI2PSI2->SetName("AcousticDampIntTDEFPSI2PSI2_" + std::to_string(ii));

      BiLinFormContext *dampIntTDEFPSI2PSI2Context = NULL;
      dampIntTDEFPSI2PSI2Context = new BiLinFormContext(dampIntTDEFPSI2PSI2, DAMPING);

      dampIntTDEFPSI2PSI2Context->SetEntities(actSDList, actSDList);
      dampIntTDEFPSI2PSI2Context->SetFeFunctions(feFunctions_[(SolutionType)(ACOU_TDEF_PSI_V_1 + ii)], feFunctions_[(SolutionType)(ACOU_TDEF_PSI_V_1 + ii)]);
      assemble_->AddBiLinearForm(dampIntTDEFPSI2PSI2Context);

      // ====================================================================
      // K_PSI2PSI2 (TDEF): stiffness integrator, TDEF (\delta_m^V term)
      // \int_{Omega_1} \delta_m^V \psi_m^{V,\prime} \psi_m^V d\Omega
      // ====================================================================
      BiLinearForm *stiffIntTDEFPSI2PSI2 = NULL;
      PtrCoefFct fncQuotV;
      fncQuotV = CoefFunction::Generate(mp_, Global::REAL,
                                        CoefXprBinOp(mp_, fncBetaV_[ii], fncGammaV_[ii], CoefXpr::OP_DIV));
      PtrCoefFct fncTermV;
      fncTermV = CoefFunction::Generate(mp_, Global::REAL,
                                        CoefXprBinOp(mp_, fncBetaV_[ii], fncQuotV, CoefXpr::OP_MULT));
      PtrCoefFct fncDeltaV;
      fncDeltaV = CoefFunction::Generate(mp_, Global::REAL,
                                         CoefXprBinOp(mp_, fncGammaV_[ii], fncTermV, CoefXpr::OP_ADD));
      if (dim_ == 2)
      {
        stiffIntTDEFPSI2PSI2 = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), fncDeltaV, 1.0, updatedGeo_);
      }
      else
      {
        stiffIntTDEFPSI2PSI2 = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), fncDeltaV, 1.0, updatedGeo_);
      }

      stiffIntTDEFPSI2PSI2->SetName("AcousticStiffIntTDEFPSI2PSI12_" + std::to_string(ii));

      BiLinFormContext *stiffIntTDEFPSI2PSI2Context = NULL;
      stiffIntTDEFPSI2PSI2Context = new BiLinFormContext(stiffIntTDEFPSI2PSI2, STIFFNESS);

      stiffIntTDEFPSI2PSI2Context->SetEntities(actSDList, actSDList);
      stiffIntTDEFPSI2PSI2Context->SetFeFunctions(feFunctions_[(SolutionType)(ACOU_TDEF_PSI_V_1 + ii)], feFunctions_[(SolutionType)(ACOU_TDEF_PSI_V_1 + ii)]);
      assemble_->AddBiLinearForm(stiffIntTDEFPSI2PSI2Context);

      // ====================================================================
      // M_PSI2P (TDEF): mass integrator, TDEF
      // -\int_{Omega_1} \psi_m^{V,\prime} \ddot{p}_a d\Omega
      // ====================================================================
      BiLinearForm *massIntTDEFPSI2P = NULL;

      if (dim_ == 2)
      {
        massIntTDEFPSI2P = new BBInt<Double>(new IdentityOperator<FeH1, 2>(), constOne, -1.0, updatedGeo_);
      }
      else
      {
        massIntTDEFPSI2P = new BBInt<Double>(new IdentityOperator<FeH1, 3>(), constOne, -1.0, updatedGeo_);
      }

      massIntTDEFPSI2P->SetName("AcousticMassIntTDEFPSI2P_" + std::to_string(ii));

      BiLinFormContext *massIntTDEFPSI2PContext = NULL;
      massIntTDEFPSI2PContext = new BiLinFormContext(massIntTDEFPSI2P, MASS);

      massIntTDEFPSI2PContext->SetEntities(actSDList, actSDList);
      massIntTDEFPSI2PContext->SetFeFunctions(feFunctions_[(SolutionType)(ACOU_TDEF_PSI_V_1 + ii)], feFunctions_[formulation_]);
      assemble_->AddBiLinearForm(massIntTDEFPSI2PContext);
    } // end loop massIntTDEFPSI2PSI2, dampIntTDEFPSI2PSI2, stiffIntTDEFPSI2PSI2 and massInteTDEFPSI2P
  }

  void AcousticPDE::DefineNcIntegrators()
  {
    //	if ( complexFluidFormulation_ )
    //		EXCEPTION("Complex fluid and NC-interfaces currently not allowed");

    StdVector<NcInterfaceInfo>::iterator ncIt = ncInterfaces_.Begin(),
                                         endIt = ncInterfaces_.End();
    for (; ncIt != endIt; ++ncIt)
    {
      switch (ncIt->type)
      {
      case NC_MORTAR:
        if (dim_ == 2)
          DefineMortarCoupling<2, 1>(formulation_, *ncIt);
        else
          DefineMortarCoupling<3, 1>(formulation_, *ncIt);
        break;
      case NC_NITSCHE:
        if (dim_ == 2)
          DefineNitscheCoupling<2, 1>(formulation_, *ncIt);
        else
          DefineNitscheCoupling<3, 1>(formulation_, *ncIt);
        break;
      default:
        EXCEPTION("Unknown type of ncInterface");
        break;
      }
    }
  }

  void AcousticPDE::DefineSurfaceIntegrators()
  {
    //========================================================================================
    // ABC boundaries
    //========================================================================================
    PtrParamNode bcNode = myParam_->Get("bcsAndLoads", ParamNode::PASS);
    LOG_DBG(acousticpde) << "Define Surface Integrator BEGIN"
                         << "\n";

    if (bcNode)
    {
      ParamNodeList abcNodes = bcNode->GetList("absorbingBCs");
      LOG_DBG(acousticpde) << "ABCs count :  " << abcNodes.GetSize() << "\n";

      for (UInt i = 0; i < abcNodes.GetSize(); i++)
      {
        std::string regionName = abcNodes[i]->Get("name")->As<std::string>();
        shared_ptr<EntityList> actSDList = ptGrid_->GetEntityList(EntityList::SURF_ELEM_LIST, regionName);
        std::string volRegName = abcNodes[i]->Get("volumeRegion")->As<std::string>();
        LOG_DBG(acousticpde) << "ABCs volRegName :  " << volRegName << "\n";

        RegionIdType aRegion = ptGrid_->GetRegion().Parse(volRegName);
        UInt iRegion = regions_.Find(aRegion);

        // check, if region has complex fluid
        PtrParamNode curRegNode =
            myParam_->Get("regionList")->GetByVal("region", "name", volRegName.c_str());

        // c0 = sqrt(bulk_modulus / density)
        PtrCoefFct dens;
        PtrCoefFct blk;
        PtrCoefFct constOne = CoefFunction::Generate(mp_, Global::REAL, "1.0");
        PtrCoefFct omegaTrg;

        if (complexFluidFormulation_)
        {
          dens = materials_[aRegion]->GetScalCoefFnc(DENSITY, Global::COMPLEX);
          blk = materials_[aRegion]->GetScalCoefFnc(ACOU_BULK_MODULUS, Global::COMPLEX);
        }
        else if (timeDomainEqFluidFormulation_ && isTDEFReg_[iRegion]) // volume region is TDEF region
        {
          LocPointMapped lpm;

          std::cout << "TDEF ABC" << std::endl;
          Double ftrg = abcNodes[i]->Get("targetFrequency")->As<UInt>();

          PtrCoefFct targetFreg = CoefFunction::Generate(mp_, Global::REAL, std::to_string(ftrg));
          omegaTrg = CoefFunction::Generate(mp_, Global::REAL,
                                                      CoefXprBinOp(mp_, targetFreg, CoefFunction::Generate(mp_, Global::REAL, "2*pi"), CoefXpr::OP_MULT));

          // TODO get dens blk at trg freq ftrg
          EvalRationalFncs(iRegion, ftrg);

          dens = CoefFunction::Generate(mp_, Global::COMPLEX,
                                        CoefXprBinOp(mp_, constOne, invTDEFDens_, CoefXpr::OP_DIV));
          blk = CoefFunction::Generate(mp_, Global::COMPLEX,
                                       CoefXprBinOp(mp_, constOne, invTDEFBlk_, CoefXpr::OP_DIV));
        }
        else
        {
          dens = materials_[aRegion]->GetScalCoefFnc(DENSITY, Global::REAL);
          blk = materials_[aRegion]->GetScalCoefFnc(ACOU_BULK_MODULUS, Global::REAL);
        }

        LOG_DBG(acousticpde) << "ABC: dens = " << dens->ToString() << "\n";
        LOG_DBG(acousticpde) << "ABC: blk  = " << blk->ToString() << "\n";

        PtrCoefFct c0;
        PtrCoefFct coeffStiffTDEF;

        // check for temperature dependency
        curRegNode = myParam_->Get("regionList")->GetByVal("region", "name", volRegName.c_str());
        std::string tempId = curRegNode->Get("temperatureId")->As<std::string>();
        if (tempId != "")
        {
          // Get result info object for flow
          shared_ptr<ResultInfo> tempInfo = GetResultInfo(HEAT_MEAN_TEMPERATURE);

          // Add the region information
          PtrParamNode tempNode = myParam_->Get("temperatureList")->GetByVal("temperature", "name", tempId.c_str());

          // Read coefficient flow coefficient function for this region and add it to flow functor
          PtrCoefFct regionTemp;
          std::set<UInt> definedDofs;
          // we assume that the temperature values are real (not complex!)
          ReadUserFieldValues(actSDList, tempNode, tempInfo->dofNames, tempInfo->entryType,
                              false, regionTemp, definedDofs, updatedGeo_);
          // gasR=287.058 J/kg K   ... universal gas constant
          // kappa=1.402, adabatic exponent for air
          PtrCoefFct constVal = CoefFunction::Generate(mp_, Global::REAL, "402.4553160");
          c0 = CoefFunction::Generate(mp_, Global::REAL,
                                      CoefXprUnaryOp(mp_, CoefXprBinOp(mp_, constVal, regionTemp, CoefXpr::OP_MULT),
                                                     CoefXpr::OP_SQRT));
          LOG_DBG(acousticpde) << "Define Surface Integrators standard c0 =" << c0->ToString() << "\n";
        }

        else
        {
          if (complexFluidFormulation_)
          {
            // correct factor is 1/(density*c0) = 1/sqrt(dens*compressModulus)
            // so we modify c0 to : sqrt(density*compressModulus)
            // division is done later
            // note that c0 in this case does not represent the speed of sound
            c0 = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprUnaryOp(mp_, CoefXprBinOp(mp_, blk, dens, CoefXpr::OP_MULT), CoefXpr::OP_SQRT));
          }
          else if (timeDomainEqFluidFormulation_ && isTDEFReg_[iRegion])
          {
            std::cout << "establishing narrow-band ABC for TDEF formulation" << std::endl;
            // narrow-band ABC for TDEF region according to Abdulkareem et al. 2018
            c0 = CoefFunction::Generate(mp_, Global::COMPLEX,
                                        CoefXprUnaryOp(mp_, CoefXprBinOp(mp_, blk, dens, CoefXpr::OP_DIV), CoefXpr::OP_SQRT));

            std::cout << "complex sos: c = " << c0->ToString() << "\n";
          }
          else
            c0 = CoefFunction::Generate(mp_, Global::REAL,
                                        CoefXprUnaryOp(mp_, CoefXprBinOp(mp_, blk, dens, CoefXpr::OP_DIV), CoefXpr::OP_SQRT));
        }
        LOG_DBG(acousticpde) << "Def Surface Integrator:  c0 =" << c0->ToString() << "\n";

        // the following part was missing which is why abc did not function for acouPotential + mechanic
        // if pde couples with mechanic, we have to multiply the density by -1
        PtrCoefFct factor;
        if (isMechCoupled_ == true && formulation_ != ACOU_PRESSURE)
        {
          // Important: In case of a general / quadratic EV problem, we must
          // ensure to have a "positive definite" matrix, i.e. we are not allowed
          // to multiply all matrices by -1!
          std::string stringFac = (analysistype_ != EIGENFREQUENCY) ? "-1.0" : "1.0";

          factor = CoefFunction::Generate(mp_, Global::REAL,
                                          CoefXprBinOp(mp_, dens, stringFac, CoefXpr::OP_MULT));
        }
        else
        {
          factor = CoefFunction::Generate(mp_, Global::REAL, "1.0");
        }
        LOG_DBG(acousticpde) << "Def Surface Integrator: factor =" << factor->ToString() << "\n";

        PtrCoefFct coeffDamp;
        if (sosAtLaplace_)
        {
          // factor for damping matrix: factor * c0
          coeffDamp = CoefFunction::Generate(mp_, Global::REAL,
                                             CoefXprBinOp(mp_, factor, c0, CoefXpr::OP_MULT));
        }
        else if(timeDomainEqFluidFormulation_){
          // multiplication by denstity required in the case of the TDEF formulation

          PtrCoefFct abcCoefTDEFre;
          PtrCoefFct abcCoefTDEFim;
          PtrCoefFct abcCoefTDEFreSQ;
          PtrCoefFct abcCoefTDEFimSQ;
          PtrCoefFct abcCoefTDEFdenom;
          abcCoefTDEFre = CoefFunction::Generate(mp_, Global::REAL, CoefXprUnaryOp(mp_, c0, CoefXpr::OP_RE));

          if(nAuxFncAC_[iRegion]==0 && nAuxFncBC_[iRegion]==0 && nAuxFncAV_[iRegion]==0 && nAuxFncBV_[iRegion]==0){
            // in the case that no pole are provided, c0 is real-valued
            abcCoefTDEFim = CoefFunction::Generate(mp_, Global::REAL, "0.0");
          }
          else{
            abcCoefTDEFim = CoefFunction::Generate(mp_, Global::REAL, CoefXprUnaryOp(mp_, c0, CoefXpr::OP_IM));
          }
          abcCoefTDEFreSQ = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, abcCoefTDEFre, abcCoefTDEFre, CoefXpr::OP_MULT));
          abcCoefTDEFimSQ = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, abcCoefTDEFim, abcCoefTDEFim, CoefXpr::OP_MULT));
          abcCoefTDEFdenom = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, abcCoefTDEFreSQ, abcCoefTDEFimSQ, CoefXpr::OP_ADD));
          
          PtrCoefFct coeffDampTmp = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, abcCoefTDEFre, abcCoefTDEFdenom, CoefXpr::OP_DIV));


          // scale by high frequency limit of inverse density
          PtrCoefFct invDensInf;
          invDensInf = materials_[aRegion]->GetScalCoefFnc(ACOU_TDEF_INVDENS_CONST, Global::REAL);
          coeffDamp = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, coeffDampTmp, invDensInf, CoefXpr::OP_MULT));
          std::cout << "coeffDamp = " << coeffDamp->ToString() << "\n";


          if(isTDEFReg_[iRegion]){

            std::cout << "omegaTrg = " << omegaTrg->ToString() << "\n";
            PtrCoefFct alpha = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, coeffDampTmp, omegaTrg, CoefXpr::OP_MULT));
            std::cout << "alpha = " << alpha->ToString() << "\n";

            // for a TDEF region case, an additional stiffnes term is required for the (narrow-band for the giben target frequency) ABC
            
            PtrCoefFct betaTmp = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, abcCoefTDEFim, abcCoefTDEFdenom, CoefXpr::OP_DIV));
            PtrCoefFct beta = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, betaTmp, omegaTrg, CoefXpr::OP_MULT));
            std::cout << "beta = " << beta->ToString() << "\n";

            coeffStiffTDEF = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, beta, invDensInf, CoefXpr::OP_MULT));            
            std::cout << "coeffStiffTDEF = " << coeffStiffTDEF->ToString() << "\n";
          }
        }
        else
        {
          // factor for damping matrix: factor / c0
          if (complexFluidFormulation_)
            coeffDamp = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, factor, c0, CoefXpr::OP_DIV));
          else
            coeffDamp = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, factor, c0, CoefXpr::OP_DIV));
        }
        LOG_DBG(acousticpde) << "Define Surface Integrator: coeffDamp =" << coeffDamp->ToString() << "\n";

        BiLinearForm *abcInt = NULL;
        if (complexFluidFormulation_)
        {
          if (dim_ == 2)
            abcInt = new BBInt<Complex>(new IdentityOperator<FeH1, 2, 1, Complex>(), coeffDamp, 1.0, updatedGeo_);
          else
            abcInt = new BBInt<Complex>(new IdentityOperator<FeH1, 3, 1, Complex>(), coeffDamp, 1.0, updatedGeo_);
        }

        else
        {
          if (dim_ == 2)
            abcInt = new BBInt<>(new IdentityOperator<FeH1, 2, 1>(), coeffDamp, 1.0, updatedGeo_);
          else
            abcInt = new BBInt<>(new IdentityOperator<FeH1, 3, 1>(), coeffDamp, 1.0, updatedGeo_);
        }

        FEMatrixType targetMatrix = DAMPING;
        if (updatedGeo_)
        {
          targetMatrix = DAMPING_UPDATE;
        }

        abcInt->SetName("abcIntegrator");
        BiLinFormContext *abcContext = new BiLinFormContext(abcInt, targetMatrix);

        abcContext->SetEntities(actSDList, actSDList);
        abcContext->SetFeFunctions(feFunctions_[formulation_], feFunctions_[formulation_]);
        feFunctions_[formulation_]->AddEntityList(actSDList);
        assemble_->AddBiLinearForm(abcContext);

        // ABC adjacent to TDEF-volume region requires an additional stiffness term
        // rho dp/dn = -  rho/c dp/dt - rho beta p; with c = real(c_complex) , beta = - imag(c_complex)
        if (timeDomainEqFluidFormulation_ && isTDEFReg_[iRegion])
        {
          BiLinearForm *abcIntSiffTDEF = NULL;
          if (dim_ == 2)
            abcIntSiffTDEF = new BBInt<>(new IdentityOperator<FeH1, 2, 1>(), coeffStiffTDEF, 1.0, updatedGeo_);
          else
            abcIntSiffTDEF = new BBInt<>(new IdentityOperator<FeH1, 3, 1>(), coeffStiffTDEF, 1.0, updatedGeo_);

          FEMatrixType targetMatrix = STIFFNESS;
          if (updatedGeo_)
          {
            targetMatrix = STIFFNESS_UPDATE;
          }

          abcIntSiffTDEF->SetName("abcIntegratorStiffnTermTDEF");
          BiLinFormContext *abcContextStiffBeta = new BiLinFormContext(abcIntSiffTDEF, targetMatrix);

          abcContextStiffBeta->SetEntities(actSDList, actSDList);
          abcContextStiffBeta->SetFeFunctions(feFunctions_[formulation_], feFunctions_[formulation_]);
          feFunctions_[formulation_]->AddEntityList(actSDList);
          assemble_->AddBiLinearForm(abcContextStiffBeta);
        }
      }

      //========================================================================================
      // Impedance boundaries
      // TODO: implement impedance BC
      //========================================================================================
      ParamNodeList impedNodes = bcNode->GetList("impedance");

      for (UInt i = 0; i < impedNodes.GetSize(); ++i)
      {
        if (complexFluidFormulation_ || timeDomainEqFluidFormulation_)
          EXCEPTION("Complex fluid and impedance currently not working");

        BiLinearForm *impedInt = NULL;
        std::string regionName = impedNodes[i]->Get("name")->As<std::string>();
        shared_ptr<EntityList> actSDList = ptGrid_->GetEntityList(EntityList::SURF_ELEM_LIST, regionName);
        std::string volRegName = impedNodes[i]->Get("volumeRegion")->As<std::string>();
        RegionIdType aRegion = ptGrid_->GetRegion().Parse(volRegName);
        bool isNormalised = impedNodes[i]->Get("isNormalised")->As<bool>();
        std::string impId = impedNodes[i]->Get("impedanceId")->As<std::string>();

        shared_ptr<CoefFunctionImpedanceModel<Complex>>
            Z_impMod(new CoefFunctionImpedanceModel<Complex>(mp_, materials_[aRegion], isNormalised));

        PtrParamNode impListNodes = myParam_->Get("impedanceList", ParamNode::PASS);
        if (!impListNodes)
        {
          EXCEPTION("No impedance list in *.xml");
        }
        ParamNodeList impChildren = impListNodes->GetChildren();
        UInt i_imp = 0;
        for (; i_imp < impChildren.GetSize(); ++i_imp)
        {
          std::string impType = impChildren[i_imp]->GetName();
          std::string actId = impChildren[i_imp]->Get("id")->As<std::string>();

          if (actId != impId)
          {
            continue;
          }
          if (impType == "mpp")
          {
            std::string impSubType = impChildren[i_imp]->Get("subtype")->As<std::string>();
            if (impSubType == "slit")
            {
              Z_impMod->GenerateSlitMpp(impChildren[i_imp]);
            }
            else if (impSubType == "circ")
            {
              Z_impMod->GenerateCircMpp(impChildren[i_imp]);
            }
            else
            {
              EXCEPTION("No such impedance type: " << impType << " with subtype: " << impSubType);
            }
          }
          else if (impType == "interpolImpedance")
          {
            Z_impMod->GenerateInterpolImpedance(impChildren[i_imp]);
          }
          else if (impType == "fctImpedance")
          {
            Z_impMod->GenerateImpedanceFct(impChildren[i_imp]);
          }
          else
          {
            EXCEPTION("No such impedance type: " << impType);
          }
          break;
        }
        if (i_imp == impChildren.GetSize())
        {
          EXCEPTION("No such impedance id found: " << impId);
        }

        if (dim_ == 2)
        {
          impedInt = new BBInt<Complex>(new IdentityOperator<FeH1, 2, 1, Complex>(), Z_impMod, 1.0, updatedGeo_);
        }
        else
        {
          impedInt = new BBInt<Complex>(new IdentityOperator<FeH1, 3, 1, Complex>(), Z_impMod, 1.0, updatedGeo_);
        }

        FEMatrixType targetMatrix = DAMPING;
        if (updatedGeo_)
        {
          targetMatrix = DAMPING_UPDATE;
        }

        impedInt->SetName("impedIntegrator");
        BiLinFormContext *impedContext = new BiLinFormContext(impedInt, targetMatrix);

        impedContext->SetEntities(actSDList, actSDList);
        impedContext->SetFeFunctions(feFunctions_[formulation_], feFunctions_[formulation_]);
        feFunctions_[formulation_]->AddEntityList(actSDList);
        assemble_->AddBiLinearForm(impedContext);
      }

      //========================================================================================
      // boundary Layers
      //========================================================================================
      ParamNodeList blNodes = bcNode->GetList("boundaryLayer");

      for (UInt i = 0; i < blNodes.GetSize(); ++i)
      {
        if (!(formulation_ == ACOU_PRESSURE) || !(this->analysistype_ == HARMONIC))
          EXCEPTION("bounadryLayer only available for harmonic pressure formulation");

        std::string regionName = blNodes[i]->Get("name")->As<std::string>();
        shared_ptr<EntityList> actSDList = ptGrid_->GetEntityList(EntityList::SURF_ELEM_LIST, regionName);
        std::string volRegName = blNodes[i]->Get("volumeRegion")->As<std::string>();
        RegionIdType volRegion = ptGrid_->GetRegion().Parse(volRegName);

        PtrCoefFct rho0 = materials_[volRegion]->GetScalCoefFnc(DENSITY, Global::REAL);
        PtrCoefFct K = materials_[volRegion]->GetScalCoefFnc(ACOU_BULK_MODULUS, Global::REAL);
        PtrCoefFct cp = CoefFunction::Generate(mp_, Global::REAL, blNodes[i]->Get("cp")->As<std::string>());
        PtrCoefFct cv = CoefFunction::Generate(mp_, Global::REAL, blNodes[i]->Get("cv")->As<std::string>());
        PtrCoefFct nu = CoefFunction::Generate(mp_, Global::REAL, blNodes[i]->Get("nu")->As<std::string>());
        PtrCoefFct k = CoefFunction::Generate(mp_, Global::REAL, blNodes[i]->Get("k")->As<std::string>());

        LOG_DBG(acousticpde) << "Define Surface Integrator: regionName =" << regionName << "\n";
        LOG_DBG(acousticpde) << "Define Surface Integrator: volRegName =" << volRegName << "\n";
        LOG_DBG(acousticpde) << "Define Surface Integrator: rho0       =" << rho0->ToString() << "\n";
        LOG_DBG(acousticpde) << "Define Surface Integrator: K          =" << K->ToString() << "\n";
        LOG_DBG(acousticpde) << "Define Surface Integrator: k          =" << k->ToString() << "\n";

        PtrCoefFct omegaHalv = CoefFunction::Generate(mp_, Global::REAL, "pi*f"); //
        // deltaV = sqrt( 2*nu/omega )
        PtrCoefFct deltaV = CoefFunction::Generate(mp_, Global::REAL, CoefXprUnaryOp(mp_, CoefXprBinOp(mp_, nu, omegaHalv, CoefXpr::OP_DIV), CoefXpr::OP_SQRT));
        // deltaT = sqrt( 2*k/(omega*rho0*cp) )
        PtrCoefFct deltaT = CoefFunction::Generate(mp_, Global::COMPLEX,
                                                   CoefXprUnaryOp(mp_,
                                                                  CoefXprBinOp(mp_,
                                                                               k, CoefXprBinOp(mp_, CoefXprBinOp(mp_, cp, rho0, CoefXpr::OP_MULT), omegaHalv, CoefXpr::OP_MULT),
                                                                               CoefXpr::OP_DIV),
                                                                  CoefXpr::OP_SQRT));
        // 1-i ... common factor for both terms
        PtrCoefFct oneMinusI = CoefFunction::Generate(mp_, Global::COMPLEX, "1", "-1");
        // gamma = cp/cv -> gamma-1 = (cp-cv)/cv
        PtrCoefFct gammaMinusOne = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, CoefXprBinOp(mp_, cp, cv, CoefXpr::OP_SUB), cv, CoefXpr::OP_DIV));

        // Mass matrix
        PtrCoefFct coefM = CoefFunction::Generate(mp_, Global::COMPLEX,
                                                  CoefXprBinOp(mp_,
                                                               CoefXprBinOp(mp_,
                                                                            CoefXprBinOp(mp_, deltaT, CoefXprBinOp(mp_, rho0, K, CoefXpr::OP_DIV), CoefXpr::OP_MULT),
                                                                            oneMinusI,
                                                                            CoefXpr::OP_MULT),
                                                               gammaMinusOne,
                                                               CoefXpr::OP_MULT));
        BiLinearForm *blmInt = NULL;
        if (dim_ == 2)
        {
          blmInt = new BBInt<Complex>(new IdentityOperator<FeH1, 2, 1, Complex>(), coefM, 0.5, updatedGeo_);
        }
        else
        {
          blmInt = new BBInt<Complex>(new IdentityOperator<FeH1, 3, 1, Complex>(), coefM, 0.5, updatedGeo_);
        }
        blmInt->SetName("blMassIntegrator");
        BiLinFormContext *blmContext = new BiLinFormContext(blmInt, MASS);
        blmContext->SetEntities(actSDList, actSDList);
        blmContext->SetFeFunctions(feFunctions_[formulation_], feFunctions_[formulation_]);
        feFunctions_[formulation_]->AddEntityList(actSDList);
        assemble_->AddBiLinearForm(blmContext);

        LOG_DBG(acousticpde) << "Define Surface Integrator boundary layer: coefM =" << coefM->ToString() << "\n";

        // Stiffness matrix
        PtrCoefFct coefK = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, deltaV, oneMinusI, CoefXpr::OP_MULT));
        BiLinearForm *blkInt = NULL;
        std::set<RegionIdType> volRegions;
        volRegions.insert(volRegion);
        if (dim_ == 2)
        {
          blkInt = new BBInt<Complex>(new GradientOperator<FeH1, 2, 1, Complex>(), coefK, -0.5, updatedGeo_);
        }
        else
        {
          blkInt = new BBInt<Complex>(new GradientOperator<FeH1, 3, 1, Complex>(), coefK, -0.5, updatedGeo_);
        }
        blkInt->SetName("blStiffnessIntegrator");
        BiLinFormContext *blkContext = new BiLinFormContext(blkInt, STIFFNESS);
        blkContext->SetEntities(actSDList, actSDList);
        blkContext->SetFeFunctions(feFunctions_[formulation_], feFunctions_[formulation_]);
        feFunctions_[formulation_]->AddEntityList(actSDList);
        assemble_->AddBiLinearForm(blkContext);

        LOG_DBG(acousticpde) << "Define Surface Integrator boundary layer: coefK =" << coefK->ToString() << "\n";
      } // boundary Layers
    }   // end if ( bcNode )
    LOG_DBG(acousticpde) << "Define Surface Integrator END"
                         << "\n";
  } // DefineSurfaceIntegrators

  void AcousticPDE::DefineRhsLoadIntegrators()
  {
    LOG_TRACE(acousticpde) << "Defining rhs load integrators for acoustic PDE";

    // Get FESpace and FeFunction
    shared_ptr<BaseFeFunction> myFct = feFunctions_[formulation_];
    shared_ptr<FeSpace> mySpace = myFct->GetFeSpace();
    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
    StdVector<shared_ptr<EntityList>> ent;
    StdVector<PtrCoefFct> coef;
    StdVector<std::string> empty;
    LinearForm *lin = NULL;

    // obtain density
    shared_ptr<CoefFunctionMulti> densFct = matCoefs_[ELEM_DENSITY];
    shared_ptr<CoefFunctionSurf> surfDens(new CoefFunctionSurf(false));
    surfDens->SetVolumeCoefs(densFct->GetRegionCoefs());

    // In the case of acou-mech coupling we have to multiply the
    // integrators by -density
    Double scalFactor = 1.0;
    std::set<RegionIdType> volRegions(regions_.Begin(), regions_.End());

    bool coefUpdateGeo;
    // ===========================
    //  NORMAL VELOCITY (surface)
    // ===========================
    LOG_DBG(acousticpde) << "Reading normal velocity";

    ReadRhsExcitation("normalVelocity", empty, ResultInfo::SCALAR, isComplex_,
                      ent, coef, coefUpdateGeo);

    for (UInt i = 0; i < ent.GetSize(); ++i)
    {
      if (sosAtLaplace_)
        EXCEPTION("NormalVelocity-BC and speed of sound at Laplace operator currently not possible");

      // ensure that list contains only surface elements
      EntityIterator it = ent[i]->GetIterator();
      UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
      if (elemDim != (dim_ - 1))
      {
        EXCEPTION("Normal velocity can only be defined on surface elements");
      }
      PtrCoefFct exValue;
      if (isMechCoupled_ == true && formulation_ != ACOU_PRESSURE)
      {
        scalFactor = -1.0;
        exValue =
            CoefFunction::Generate(mp_, part,
                                   CoefXprBinOp(mp_, coef[i], surfDens, CoefXpr::OP_MULT));
      }
      else
      {
        exValue = coef[i];
      }

      if (formulation_ == ACOU_POTENTIAL)
      {
        if (dim_ == 2)
        {
          if (isComplex_)
          {
            lin = new BUIntegrator<Complex, true>(new IdentityOperator<FeH1, 2, 1>(),
                                                  scalFactor, exValue, volRegions, coefUpdateGeo);
          }
          else
          {
            lin = new BUIntegrator<Double, true>(new IdentityOperator<FeH1, 2, 1>(),
                                                 scalFactor, exValue, volRegions, coefUpdateGeo);
          }
        }
        else
        {
          if (isComplex_)
          {
            lin = new BUIntegrator<Complex, true>(new IdentityOperator<FeH1, 3, 1>(),
                                                  1.0, exValue, volRegions, coefUpdateGeo);
          }
          else
          {
            lin = new BUIntegrator<Double, true>(new IdentityOperator<FeH1, 3, 1>(),
                                                 1.0, exValue, volRegions, coefUpdateGeo);
          }
        }
      }
      else if (formulation_ == ACOU_PRESSURE && isComplex_)
      {

        // in this case the pressure can be related to the
        // normal velocity as
        // p_n = - j*Omega*v_n*rho
        PtrCoefFct tmp = CoefFunction::Generate(mp_, part, "0.0",
                                                "-2*pi*f");
        if (complexFluidFormulation_)
        {
          // do not need multiplication with density
          exValue = CoefFunction::Generate(mp_, part,
                                           CoefXprBinOp(mp_, tmp, exValue, CoefXpr::OP_MULT));
        }
        else
        {
          PtrCoefFct tmp2 = CoefFunction::Generate(mp_, part,
                                                   CoefXprBinOp(mp_, tmp, exValue, CoefXpr::OP_MULT));
          exValue = CoefFunction::Generate(mp_, part, CoefXprBinOp(mp_, tmp2, surfDens, CoefXpr::OP_MULT));
        }

        if (dim_ == 2)
        {
          lin = new BUIntegrator<Complex, true>(new IdentityOperator<FeH1, 2, 1>(),
                                                1.0, exValue, volRegions, coefUpdateGeo);
        }
        else
        {
          lin = new BUIntegrator<Complex, true>(new IdentityOperator<FeH1, 3, 1>(),
                                                1.0, exValue, volRegions, coefUpdateGeo);
        }
      }
      else
      {
        EXCEPTION("Normal velocity can only be prescribed for potential "
                  << "formulation or for harmonic pressure formulation")
      }
      lin->SetName("NormalVelocityIntegrator");
      LinearFormContext *ctx = new LinearFormContext(lin);
      ctx->SetEntities(ent[i]);
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
    }

    // ===========================
    //  TOTAL VELOCITY (surface)
    // ===========================
    LOG_DBG(acousticpde) << "Reading total velocity";
    StdVector<std::string> vecDofNames;
    if (dim_ == 3)
      vecDofNames = "x", "y", "z";
    if (dim_ == 2 && !isaxi_)
      vecDofNames = "x", "y";
    if (dim_ == 2 && isaxi_)
      vecDofNames = "r", "z";
    ReadRhsExcitation("velocity", vecDofNames, ResultInfo::VECTOR, isComplex_,
                      ent, coef, coefUpdateGeo);

    for (UInt i = 0; i < ent.GetSize(); ++i)
    {
      if (sosAtLaplace_)
        EXCEPTION("RHS and speed of sound at Laplace operator currently not possible");

      // ensure that list contains only surface elements
      EntityIterator it = ent[i]->GetIterator();
      UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
      if (elemDim != (dim_ - 1))
      {
        EXCEPTION("Normal velocity can only be defined on surface elements");
      }

      PtrCoefFct exValue;
      if (isMechCoupled_ == true && formulation_ != ACOU_PRESSURE)
      {
        scalFactor = -1.0;
        exValue =
            CoefFunction::Generate(mp_, part,
                                   CoefXprBinOp(mp_, coef[i], surfDens, CoefXpr::OP_MULT));
      }
      else
      {
        exValue = coef[i];
      }
      if (formulation_ == ACOU_POTENTIAL)
      {
        if (dim_ == 2)
        {
          if (isComplex_)
          {
            lin = new BUIntegrator<Complex, true>(new IdentityOperatorNormal<FeH1, 2>(),
                                                  scalFactor, exValue, volRegions, coefUpdateGeo);
          }
          else
          {
            lin = new BUIntegrator<Double, true>(new IdentityOperatorNormal<FeH1, 2>(),
                                                 scalFactor, exValue, volRegions, coefUpdateGeo);
          }
        }
        else
        {
          if (isComplex_)
          {
            lin = new BUIntegrator<Complex, true>(new IdentityOperatorNormal<FeH1, 3>(),
                                                  scalFactor, exValue, volRegions, coefUpdateGeo);
          }
          else
          {
            lin = new BUIntegrator<Double, true>(new IdentityOperatorNormal<FeH1, 3>(),
                                                 scalFactor, exValue, volRegions, coefUpdateGeo);
          }
        }
      }
      else if (formulation_ == ACOU_PRESSURE && isComplex_)
      {

        // in this case the pressure can be related to the
        // normal velocity as
        // p_n = - j*Omega*v_n*rho
        PtrCoefFct tmp = CoefFunction::Generate(mp_, part, "0.0",
                                                "-2*pi*f");
        if (complexFluidFormulation_)
        {
          // do not need multiplication with density
          exValue = CoefFunction::Generate(mp_, part,
                                           CoefXprBinOp(mp_, tmp, exValue, CoefXpr::OP_MULT));
        }
        else
        {
          PtrCoefFct tmp2 = CoefFunction::Generate(mp_, part, CoefXprBinOp(mp_, tmp, exValue, CoefXpr::OP_MULT));
          exValue = CoefFunction::Generate(mp_, part, CoefXprBinOp(mp_, surfDens, tmp2, CoefXpr::OP_MULT));
        }

        if (dim_ == 2)
        {
          lin = new BUIntegrator<Complex, true>(new IdentityOperatorNormal<FeH1, 2>(),
                                                scalFactor, exValue, volRegions, coefUpdateGeo);
        }
        else
        {
          lin = new BUIntegrator<Complex, true>(new IdentityOperatorNormal<FeH1, 3>(),
                                                scalFactor, exValue, volRegions, coefUpdateGeo);
        }
      }
      else
      {
        EXCEPTION("Normal velocity can only be prescribed for potential "
                  << "formulation or for harmonic pressure formulation")
      }
      lin->SetName("VelocityIntegrator");
      LinearFormContext *ctx = new LinearFormContext(lin);
      ctx->SetEntities(ent[i]);
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
    }

    // ===========================
    //  NORMAL ACCELERATION (surface)
    // ===========================
    LOG_DBG(acousticpde) << "Reading normal acceleration";

    ReadRhsExcitation("normalAcceleration", empty, ResultInfo::SCALAR, isComplex_,
                      ent, coef, coefUpdateGeo);
    for (UInt i = 0; i < ent.GetSize(); ++i)
    {
      if (sosAtLaplace_)
        EXCEPTION("RHS and speed of sound at Laplace operator currently not possible");

      // ensure that list contains only surface elements
      EntityIterator it = ent[i]->GetIterator();
      UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
      if (elemDim != (dim_ - 1))
      {
        EXCEPTION("Normal acceleration can only be defined on surface elements");
      }
      scalFactor = 1.0;
      PtrCoefFct exValue;
      if (isMechCoupled_ == true && formulation_ != ACOU_PRESSURE)
      {
        EXCEPTION("Normal acceleration can only be prescribed for pressure"
                  << "formulation");
      }
      else
      {
        exValue = coef[i];
      }

      if (formulation_ == ACOU_PRESSURE)
      {
        exValue =
            CoefFunction::Generate(mp_, part,
                                   CoefXprBinOp(mp_, coef[i], surfDens, CoefXpr::OP_MULT));

        if (dim_ == 2)
        {
          if (isComplex_)
          {
            lin = new BUIntegrator<Complex, true>(new IdentityOperator<FeH1, 2, 1>(),
                                                  scalFactor, exValue, volRegions, coefUpdateGeo);
          }
          else
          {
            lin = new BUIntegrator<Double, true>(new IdentityOperator<FeH1, 2, 1>(),
                                                 scalFactor, exValue, volRegions, coefUpdateGeo);
          }
        }
        else
        {
          if (isComplex_)
          {
            lin = new BUIntegrator<Complex, true>(new IdentityOperator<FeH1, 3, 1>(),
                                                  scalFactor, exValue, volRegions, coefUpdateGeo);
          }
          else
          {
            lin = new BUIntegrator<Double, true>(new IdentityOperator<FeH1, 3, 1>(),
                                                 scalFactor, exValue, volRegions, coefUpdateGeo);
          }
        }
      }
      else if (formulation_ == ACOU_POTENTIAL)
      {
        EXCEPTION("Normal acceleration can only be prescribed for pressure "
                  << "formulation")
      }
      lin->SetName("NormalAccelerationIntegrator");
      LinearFormContext *ctx = new LinearFormContext(lin);
      ctx->SetEntities(ent[i]);
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
    }

    // ===========================
    //  PRESSURE (surface)
    // ===========================
    // Only possible for harmonic potential formulation
    ReadRhsExcitation("pressure", empty, ResultInfo::SCALAR, isComplex_,
                      ent, coef, coefUpdateGeo);
    if (formulation_ == ACOU_POTENTIAL && isComplex_)
    {
      for (UInt i = 0; i < ent.GetSize(); ++i)
      {
        if (sosAtLaplace_)
          EXCEPTION("RHS and speed of sound at Laplace operator currently not possible");

        // ensure that list contains only surface elements
        EntityIterator it = ent[i]->GetIterator();
        UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
        if (elemDim != (dim_ - 1))
        {
          EXCEPTION("Pressure can only be defined on surface elements");
        }
        PtrCoefFct exValue;
        if (isMechCoupled_ == true && formulation_ != ACOU_PRESSURE)
        {
          scalFactor = -1.0;
          exValue =
              CoefFunction::Generate(mp_, part,
                                     CoefXprBinOp(mp_, coef[i], surfDens, CoefXpr::OP_MULT));
        }
        else
        {
          exValue = coef[i];
        }
        // psi_n = j * 1 / (omega*rho) * p
        PtrCoefFct tmp =
            CoefFunction::Generate(mp_, Global::COMPLEX, "0", "1/(2*pi*f)");
        PtrCoefFct tmp2 =
            CoefFunction::Generate(mp_, Global::COMPLEX,
                                   CoefXprBinOp(mp_, tmp, surfDens, CoefXpr::OP_DIV));
        exValue = CoefFunction::Generate(mp_, part,
                                         CoefXprBinOp(mp_, exValue, tmp2,
                                                      CoefXpr::OP_MULT));
        shared_ptr<InhomDirichletBc> actBc(new InhomDirichletBc);
        actBc->entities = ent[i];
        actBc->result = myFct->GetResultInfo();
        actBc->value = exValue;
        actBc->dofs.insert(0);
        actBc->updatedGeo = coefUpdateGeo;
        myFct->AddInhomDirichletBc(actBc);
      } // loop entities
    }   // if

    // =====================
    //  ACOUSTIC SOURCE DENSITY
    // =====================
    // LOG_DBG(heatcondpde) << "Reading acoustic source density";

    ReadRhsExcitation("acouSourceDensity", empty,
                      ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo);
    for (UInt i = 0; i < ent.GetSize(); ++i)
    {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST)
      {
        EXCEPTION("Acoustic source density must be defined on elements")
      }
      EntityIterator it = ent[i]->GetIterator();
      it.Begin();

      if (isComplex_)
      {
        lin = new BUIntegrator<Complex>(new IdentityOperator<FeH1>(),
                                        Complex(1.0), coef[i], coefUpdateGeo);
      }
      else
      {
        lin = new BUIntegrator<Double>(new IdentityOperator<FeH1>(),
                                       1.0, coef[i], coefUpdateGeo);
      }
      lin->SetName("AcousticSourceDensityInt");
      LinearFormContext *ctx = new LinearFormContext(lin);
      ctx->SetEntities(ent[i]);
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);

      RegionIdType regionId = ent[i]->GetRegion();
      acousticSourceDensityCoef_->AddRegion(regionId, coef[i]);
    } // for

    // =====================================
    //  rhsValues for e.g. for aeroacoustics
    // =====================================
    ReadRhsExcitation("rhsValues", empty, ResultInfo::SCALAR, isComplex_,
                      ent, coef, coefUpdateGeo);
    for (UInt i = 0; i < ent.GetSize(); ++i)
    {
      if (sosAtLaplace_)
        EXCEPTION("RHS and speed of sound at Laplace operator currently not possible");

      coef[i]->SetConservative(true);
      // std::cout << "ADD" << std::endl;
      coef[i]->SetFeFunction(feFunctions_[formulation_], formulation_);

      if (!approxSourceWithDeltaFnc_ && coef[i]->GetInverseType() == CoefFunction::INVSOURCE)
      {
        std::cout << "Add BUIntegrator SRC" << std::endl;
        // check type of entitylist
        if (ent[i]->GetType() == EntityList::NODE_LIST)
        {
          EXCEPTION("Acoustic source density must be defined on elements")
        }
        EntityIterator it = ent[i]->GetIterator();
        it.Begin();

        if (isComplex_)
        {
          lin = new BUIntegrator<Complex>(new IdentityOperator<FeH1>(),
                                          Complex(1.0), coef[i], coefUpdateGeo);
        }
        else
        {
          lin = new BUIntegrator<Double>(new IdentityOperator<FeH1>(),
                                         1.0, coef[i], coefUpdateGeo);
        }
        lin->SetName("AcousticSourceDensityInverseInt");
        LinearFormContext *ctx = new LinearFormContext(lin);
        ctx->SetEntities(ent[i]);
        ctx->SetFeFunction(myFct);
        assemble_->AddLinearForm(ctx);

        RegionIdType regionId = ent[i]->GetRegion();
        acousticSourceDensityCoef_->AddRegion(regionId, coef[i]);
      }

      // if ( coef[i]->GetInverseType() == CoefFunction::INVSOURCE )
      this->rhsFeFunctions_[formulation_]->AddLoadCoefFunction(coef[i], ent[i]);
    }

    // ================
    //  RHS DENSITY Vector
    //	For reading vectoriel source terms of the lightill analogy (div(T))
    //      Diss Freidhager 2022
    // ================
    LOG_DBG(acousticpde) << "Reading rhsDensityVector";
    ReadRhsExcitation("rhsDensityVector", vecDofNames, ResultInfo::VECTOR,
                      isComplex_, ent, coef, updatedGeo_);

    // perform checks

    for (UInt i = 0; i < ent.GetSize(); ++i)
    {
      std::cout << "ausgabe     " << i << std::endl;
      // check type of entitylist
      // if (ent[i]->GetType() == EntityList::NODE_LIST) {
      //  EXCEPTION("Rhs density Vector must be defined on elements, because BUIntegrator requires a density")
      //}
      EntityIterator it = ent[i]->GetIterator();
      it.Begin();

      if (formulation_ == ACOU_POTENTIAL)
      {
        EXCEPTION("Using rhsDensityVector based on Lighthill's equation -> AcouPressure Formulation!")
      } // closes if AcouPotential
      else
      {
        if (isComplex_)
        {
          std::cout << "\t \t \t Using complexe vector formulation of Lighthill div(T) in grad(phi)" << std::endl;
          scalFactor = 1.0;
          if (dim_ == 2)
          {

            lin = new BUIntegrator<Complex>(new GradientOperator<FeH1, 2, 1, Complex()>(), Complex(scalFactor), coef[i], updatedGeo_);
          }
          else
          {

            lin = new BUIntegrator<Complex>(new GradientOperator<FeH1, 3, 1, Complex()>(), Complex(scalFactor), coef[i], updatedGeo_);
          }
        }
        else
        { // Transiente Simulation

          std::cout << "\t \t \t Using vector formulation of Lighthill div(T) in grad(phi)" << std::endl;
          scalFactor = 1.0;
          if (dim_ == 2)
          {
            lin = new BUIntegrator<Double, false>(new GradientOperator<FeH1, 2, 1>(), scalFactor, coef[i], updatedGeo_);
          }
          else
          {
            lin = new BUIntegrator<Double, false>(new GradientOperator<FeH1, 3, 1>(), scalFactor, coef[i], updatedGeo_);
          }
        }
      }

      lin->SetName("RhsDensityVectorInt");
      LinearFormContext *ctx = new LinearFormContext(lin);
      ctx->SetEntities(ent[i]);
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
    ReadRhsExcitation("BCVector", vecDofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo);

    for (UInt i = 0; i < ent.GetSize(); ++i)
    {
      EntityIterator it = ent[i]->GetIterator();
      it.Begin();
      // check type of entitylist
      // if (ent[i]->GetType() == EntityList::NODE_LIST) {
      // EXCEPTION("BCVector must be defined on elements, because BUIntegrator requires a density")
      // }

      // ensure that list contains only surface elements
      UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
      if (elemDim != (dim_ - 1))
      {
        EXCEPTION("BCVector can only be defined on surface elements");
      }

      if (sosAtLaplace_)
      {
        EXCEPTION("BCVector not defined for sosAtLaplace");
      } // closes sosAtLaplace
      else
      {
        if (isComplex_)
        {
          std::cout << "\t \t \t Surface Term of vector formulation of Lighthill's analogy" << std::endl;
          scalFactor = 1.0;
          if (dim_ == 2)
          {
            lin = new BUIntegrator<Complex, true>(new IdentityOperatorNormal<FeH1, 2>(),
                                                  Complex(scalFactor), coef[i], volRegions, coefUpdateGeo); // IdentityOperatorNormal
          }
          else
          {
            lin = new BUIntegrator<Complex, true>(new IdentityOperatorNormal<FeH1, 3>(),
                                                  Complex(scalFactor), coef[i], volRegions, coefUpdateGeo);
          }
        }
        else
        {

          std::cout << "\t \t \t Surface Term of vector formulation of Lighthill's analogy" << std::endl;
          scalFactor = 1.0;
          if (dim_ == 2)
          {
            lin = new BUIntegrator<Double, true>(new IdentityOperatorNormal<FeH1, 2>(),
                                                 scalFactor, coef[i], volRegions, coefUpdateGeo);
          }
          else
          {
            lin = new BUIntegrator<Double, true>(new IdentityOperatorNormal<FeH1, 3>(),
                                                 scalFactor, coef[i], volRegions, coefUpdateGeo);
          }
        }
      }

      lin->SetName("BCVector");
      LinearFormContext *ctx = new LinearFormContext(lin);
      ctx->SetEntities(ent[i]);
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);

      //               RegionIdType regionId = ent[i]->GetRegion();
      //               BCSurfaceVectorCoef_->AddRegion(  regionId, coef[i] );
    }

    // ================
    //  RHS DENSITY
    // ================
    ReadRhsExcitation("rhsDensity", empty,
                      ResultInfo::SCALAR, isComplex_, ent, coef, updatedGeo_);
    for (UInt i = 0; i < ent.GetSize(); ++i)
    {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST)
      {
        EXCEPTION("Rhs density must be defined on elements")
      }
      if (sosAtLaplace_)
      {
        // get region id
        std::string volRegName = ent[i]->GetName();
        RegionIdType aRegion = ptGrid_->GetRegion().Parse(volRegName);

        // get adiabatic exponent
        PtrCoefFct density = materials_[aRegion]->GetScalCoefFnc(DENSITY, Global::REAL);

        // get speed of sound and square it
        std::map<RegionIdType, PtrCoefFct> c0Fcts = matCoefs_[ACOU_ELEM_SPEED_OF_SOUND]->GetRegionCoefs();
        PtrCoefFct regionC0 = c0Fcts[aRegion];
        PtrCoefFct regionC02 = CoefFunction::Generate(mp_, Global::REAL,
                                                      CoefXprBinOp(mp_, regionC0, regionC0, CoefXpr::OP_MULT));

        PtrCoefFct tmpCoef = CoefFunction::Generate(mp_, Global::REAL,
                                                    CoefXprBinOp(mp_, regionC02, density, CoefXpr::OP_MULT));

        if (isComplex_)
        {
          PtrCoefFct coefRHS = CoefFunction::Generate(mp_, Global::COMPLEX,
                                                      CoefXprBinOp(mp_, tmpCoef, coef[i], CoefXpr::OP_MULT));
          lin = new BUIntegrator<Complex>(new IdentityOperator<FeH1>(),
                                          Complex(scalFactor), coefRHS, updatedGeo_);
        }
        else
        {
          PtrCoefFct coefRHS = CoefFunction::Generate(mp_, Global::REAL,
                                                      CoefXprBinOp(mp_, tmpCoef, coef[i], CoefXpr::OP_MULT));
          lin = new BUIntegrator<Double>(new IdentityOperator<FeH1>(),
                                         scalFactor, coefRHS, updatedGeo_);
        }
      }
      else
      {
        if (isComplex_)
        {
          lin = new BUIntegrator<Complex>(new IdentityOperator<FeH1>(),
                                          Complex(scalFactor), coef[i], updatedGeo_);
        }
        else
        {
          lin = new BUIntegrator<Double>(new IdentityOperator<FeH1>(),
                                         scalFactor, coef[i], updatedGeo_);
        }
      }
      lin->SetName("RhsDensityInt");
      LinearFormContext *ctx = new LinearFormContext(lin);
      ctx->SetEntities(ent[i]);
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
    } // for

    // ================
    //  RHS Mass time derivative
    // ================
    ReadRhsExcitation("rhsMassTimeDeriv", empty,
                      ResultInfo::SCALAR, isComplex_, ent, coef, updatedGeo_);
    for (UInt i = 0; i < ent.GetSize(); ++i)
    {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST)
      {
        EXCEPTION("rhsMassTimeDeriv must be defined on elements");
      }
      if (!sosAtLaplace_)
        EXCEPTION("rhsMassTimeDeriv currently just available for speed of sound at Laplace operator");

      // source term is partial time derivative of speed of sound squared
      // multiplied by mass, which has been evaluated by "ReadRhsExcitation"
      // and stored in coef[i]
      //      if ( this->analysistype_ == HARMONIC )
      //    	  EXCEPTION("rhsMassTimeDeriv in harmonic case for speed of sound at Laplace operator\
//    		         currently not available");

      // get region id
      std::string volRegName = ent[i]->GetName();
      RegionIdType aRegion = ptGrid_->GetRegion().Parse(volRegName);

      // get adiabatic exponent
      PtrCoefFct aExp = materials_[aRegion]->GetScalCoefFnc(FLUID_ADIABATIC_EXPONENT, Global::REAL);

      // get speed of sound and square it
      std::map<RegionIdType, PtrCoefFct> c0Fcts = matCoefs_[ACOU_ELEM_SPEED_OF_SOUND]->GetRegionCoefs();
      PtrCoefFct regionC0 = c0Fcts[aRegion];
      PtrCoefFct regionC02 = CoefFunction::Generate(mp_, Global::REAL,
                                                    CoefXprBinOp(mp_, regionC0, regionC0, CoefXpr::OP_MULT));

      PtrCoefFct tmpCoef = CoefFunction::Generate(mp_, Global::REAL,
                                                  CoefXprBinOp(mp_, regionC02, aExp, CoefXpr::OP_DIV));

      if (isComplex_)
      {
        PtrCoefFct coefRHS = CoefFunction::Generate(mp_, Global::COMPLEX,
                                                    CoefXprBinOp(mp_, tmpCoef, coef[i], CoefXpr::OP_MULT));

        lin = new BUIntegrator<Complex>(new IdentityOperator<FeH1>(),
                                        scalFactor, coefRHS, updatedGeo_);
      }
      else
      {
        PtrCoefFct coefRHS = CoefFunction::Generate(mp_, Global::REAL,
                                                    CoefXprBinOp(mp_, tmpCoef, coef[i], CoefXpr::OP_MULT));

        lin = new BUIntegrator<Double>(new IdentityOperator<FeH1>(),
                                       scalFactor, coefRHS, updatedGeo_);
      }
      lin->SetName("RhsMassTimeDerivInt");
      LinearFormContext *ctx = new LinearFormContext(lin);
      ctx->SetEntities(ent[i]);
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
    } // for

    // ================
    //  RHS Mass convective term
    // ================
    ReadRhsExcitation("rhsMassConvective", empty,
                      ResultInfo::SCALAR, isComplex_, ent, coef, updatedGeo_);
    for (UInt i = 0; i < ent.GetSize(); ++i)
    {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST)
      {
        EXCEPTION("rhsMassConvective must be defined on elements");
      }

      // source term is convective derivative of speed of sound squared
      // multiplied by mass, which has been evaluated by "ReadRhsExcitation"
      // and stored in coef[i]

      // get region id
      std::string volRegName = ent[i]->GetName();
      RegionIdType aRegion = ptGrid_->GetRegion().Parse(volRegName);

      // get adiabatic exponent
      PtrCoefFct aExp = materials_[aRegion]->GetScalCoefFnc(FLUID_ADIABATIC_EXPONENT, Global::REAL);

      // get speed of sound and square it
      std::map<RegionIdType, PtrCoefFct> c0Fcts = matCoefs_[ACOU_ELEM_SPEED_OF_SOUND]->GetRegionCoefs();
      PtrCoefFct regionC0 = c0Fcts[aRegion];
      PtrCoefFct regionC02 = CoefFunction::Generate(mp_, Global::REAL,
                                                    CoefXprBinOp(mp_, regionC0, regionC0, CoefXpr::OP_MULT));

      PtrCoefFct tmpCoef = CoefFunction::Generate(mp_, Global::REAL,
                                                  CoefXprBinOp(mp_, regionC02, aExp, CoefXpr::OP_DIV));

      BaseBOperator *bOp;
      if (dim_ == 2)
        bOp = new ConvectiveOperator<FeH1, 2, 1>();
      else
        bOp = new ConvectiveOperator<FeH1, 3, 1>();
      bOp->SetCoefFunction(meanFlowCoef_); // meanFlow )

      if (isComplex_)
      {
        PtrCoefFct coefRHS = CoefFunction::Generate(mp_, Global::COMPLEX,
                                                    CoefXprBinOp(mp_, tmpCoef, coef[i], CoefXpr::OP_MULT));
        lin = new BUIntegrator<Complex>(bOp, scalFactor, coefRHS, updatedGeo_);
      }
      else
      {
        PtrCoefFct coefRHS = CoefFunction::Generate(mp_, Global::REAL,
                                                    CoefXprBinOp(mp_, tmpCoef, coef[i], CoefXpr::OP_MULT));
        lin = new BUIntegrator<Double>(bOp, scalFactor, coefRHS, updatedGeo_);
      }

      lin->SetName("RhsMassConvectiveInt");
      LinearFormContext *ctx = new LinearFormContext(lin);
      ctx->SetEntities(ent[i]);
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
    } // for
  }

  void AcousticPDE::DefineSolveStep()
  {
    solveStep_ = new StdSolveStep(*this);
  }

  void AcousticPDE::DefinePrimaryResults()
  {

    // === Primary result according to definition ===
    shared_ptr<ResultInfo> res1(new ResultInfo);
    if (formulation_ == ACOU_PRESSURE)
    {
      res1->resultType = ACOU_PRESSURE;
      res1->dofNames = "";
      res1->unit = "Pa";
    }
    else
    {
      res1->resultType = ACOU_POTENTIAL;
      res1->dofNames = "";
      res1->unit = "m^2/s";
    }
    res1->definedOn = ResultInfo::NODE;
    res1->entryType = ResultInfo::SCALAR;
    feFunctions_[formulation_]->SetResultInfo(res1);
    results_.Push_back(res1);
    res1->SetFeFunction(feFunctions_[formulation_]);
    DefineFieldResult(feFunctions_[formulation_], res1);

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    if (formulation_ == ACOU_PRESSURE)
    {
      hdbcSolNameMap_[ACOU_PRESSURE] = "soundSoft";
      idbcSolNameMap_[ACOU_PRESSURE] = "pressure";
    }
    else
    {
      hdbcSolNameMap_[ACOU_POTENTIAL] = "soundSoft";
      idbcSolNameMap_[ACOU_POTENTIAL] = "potential";
    }

    // === ACOUSTIC RHS ===
    shared_ptr<ResultInfo> rhs(new ResultInfo);
    rhs->resultType = ACOU_RHS_LOAD;
    rhs->dofNames = "";
    rhs->unit = "?";
    rhs->definedOn = results_[0]->definedOn;
    rhs->entryType = ResultInfo::SCALAR;
    this->rhsFeFunctions_[formulation_]->SetResultInfo(rhs);
    DefineFieldResult(this->rhsFeFunctions_[formulation_], rhs);
    results_.Push_back(rhs);
    availResults_.insert(rhs);

    // creates the mean flow
    StdVector<std::string> vecDofNames;
    if (ptGrid_->GetDim() == 3)
    {
      vecDofNames = "x", "y", "z";
    }
    else
    {
      if (ptGrid_->IsAxi())
      {
        vecDofNames = "r", "z";
      }
      else
      {
        vecDofNames = "x", "y";
      }
    }

    // create result info for mean flow velocity
    CreateMeanFlowFunction(vecDofNames);

    // create result info for mean temperature field
    //(mean in the sense of time averaged)
    shared_ptr<ResultInfo> temperature(new ResultInfo);
    temperature->resultType = HEAT_MEAN_TEMPERATURE;
    temperature->dofNames = "";
    temperature->unit = "K";

    temperature->definedOn = ResultInfo::NODE;
    temperature->entryType = ResultInfo::SCALAR;

    meanTemperatureCoef_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, false));
    DefineFieldResult(meanTemperatureCoef_, temperature);

    // ==============================================
    //  Define CoefFunctions for material parameters
    // ==============================================

    // === DENSITY ===
    shared_ptr<ResultInfo> density(new ResultInfo);
    density->resultType = ELEM_DENSITY;
    density->dofNames = "";
    density->unit = "kg/m^3";
    density->definedOn = ResultInfo::ELEMENT;
    density->entryType = ResultInfo::SCALAR;
    shared_ptr<CoefFunctionMulti> densFct(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1,
                                                                complexFluidFormulation_));
    matCoefs_[ELEM_DENSITY] = densFct;
    DefineFieldResult(densFct, density);

    // === SPEED OF SOUND ===
    shared_ptr<ResultInfo> sos(new ResultInfo);
    sos->resultType = ACOU_ELEM_SPEED_OF_SOUND;
    sos->dofNames = "";
    sos->unit = "m/s";
    sos->definedOn = ResultInfo::ELEMENT;
    sos->entryType = ResultInfo::SCALAR;
    shared_ptr<CoefFunctionMulti> sosFct(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1,
                                                               complexFluidFormulation_));
    matCoefs_[ACOU_ELEM_SPEED_OF_SOUND] = sosFct;
    DefineFieldResult(sosFct, sos);

    // === PML DAMPING FACTORS ===
    // if( matCoefs_.find(PML_DAMP_FACTOR) != matCoefs_.end() ) {
    shared_ptr<ResultInfo> pml(new ResultInfo);
    pml->resultType = PML_DAMP_FACTOR;
    pml->dofNames = vecDofNames;
    // pml->dofNames = "";
    pml->unit = "";
    pml->definedOn = ResultInfo::ELEMENT;
    pml->entryType = ResultInfo::VECTOR;
    shared_ptr<CoefFunctionMulti> pmlFct(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1,
                                                               isComplex_));
    matCoefs_[PML_DAMP_FACTOR] = pmlFct;
    DefineFieldResult(pmlFct, pml);
    //}

    // === PML AUX Variables ===
    if (this->isTimeDomPML_)
    {
      if (!this->isAPML_ && dim_ == 3)
      {
        shared_ptr<ResultInfo> pmlScal(new ResultInfo);
        pmlScal->resultType = ACOU_PMLAUXSCALAR;
        pmlScal->dofNames = "";
        pmlScal->unit = "-";
        pmlScal->definedOn = ResultInfo::NODE;
        pmlScal->entryType = ResultInfo::SCALAR;
        feFunctions_[ACOU_PMLAUXSCALAR]->SetResultInfo(pmlScal);
        results_.Push_back(pmlScal);
        pmlScal->SetFeFunction(feFunctions_[ACOU_PMLAUXSCALAR]);
        DefineFieldResult(feFunctions_[ACOU_PMLAUXSCALAR], pmlScal);
      }

      shared_ptr<ResultInfo> pmlVec(new ResultInfo);
      pmlVec->resultType = ACOU_PMLAUXVEC;
      pmlVec->dofNames = vecDofNames;
      pmlVec->unit = "-";
      pmlVec->definedOn = ResultInfo::NODE;
      pmlVec->entryType = ResultInfo::VECTOR;
      feFunctions_[ACOU_PMLAUXVEC]->SetResultInfo(pmlVec);
      results_.Push_back(pmlVec);
      pmlVec->SetFeFunction(feFunctions_[ACOU_PMLAUXVEC]);
      DefineFieldResult(feFunctions_[ACOU_PMLAUXVEC], pmlVec);
    }

    // === TDEF AUX Variables ===

    // loop over all materials and use largest number of poles (if multiple can be defined, which I think is the case atm)

    if (this->timeDomainEqFluidFormulation_)
    {
      for (unsigned int i = 0; i < nAuxFncAC_.Max(); i++)
      {
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
      for (unsigned int i = 0; i < nAuxFncBC_.Max(); i++)
      {
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
      for (unsigned int i = 0; i < nAuxFncAV_.Max(); i++)
      {
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
      for (unsigned int i = 0; i < nAuxFncBV_.Max(); i++)
      {
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

  void AcousticPDE::FinalizePostProcResults()
  {
    // first call base class method
    SinglePDE::FinalizePostProcResults();
  }

  void AcousticPDE::DefinePostProcResults()
  {

    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
    shared_ptr<BaseFeFunction> feFct = feFunctions_[formulation_];
    shared_ptr<ResultInfo> res1 = feFct->GetResultInfo();

    StdVector<std::string> vecDofNames;
    if (ptGrid_->GetDim() == 3)
    {
      vecDofNames = "x", "y", "z";
    }
    else
    {
      if (ptGrid_->IsAxi())
      {
        vecDofNames = "r", "z";
      }
      else
      {
        vecDofNames = "x", "y";
      }
    }

    // === PRESSURE / POTENTIAL - 1.DERIVATIVE ===
    shared_ptr<ResultInfo> deriv1(new ResultInfo);
    if (formulation_ == ACOU_POTENTIAL)
    {
      deriv1->resultType = ACOU_POTENTIAL_DERIV_1;
      deriv1->dofNames = "";
      deriv1->unit = "m^2/s^2";
    }
    else
    {
      deriv1->resultType = ACOU_PRESSURE_DERIV_1;
      deriv1->dofNames = "";
      deriv1->unit = "Pa/s";
    }
    deriv1->entryType = res1->entryType;
    deriv1->definedOn = res1->definedOn;
    availResults_.insert(deriv1);
    DefineTimeDerivResult(deriv1->resultType, 1, formulation_);

    // === PRESSURE / POTENTIAL - 2.DERIVATIVE ===
    shared_ptr<ResultInfo> deriv2(new ResultInfo);
    if (formulation_ == ACOU_POTENTIAL)
    {
      deriv2->resultType = ACOU_POTENTIAL_DERIV_2;
      deriv2->dofNames = "";
      deriv2->unit = "m^2/s^3";
    }
    else
    {
      deriv2->resultType = ACOU_PRESSURE_DERIV_2;
      deriv2->dofNames = "";
      deriv2->unit = "Pa/s^2";
    }
    deriv2->entryType = res1->entryType;
    deriv2->definedOn = res1->definedOn;
    availResults_.insert(deriv2);
    DefineTimeDerivResult(deriv2->resultType, 2, formulation_);

    // define acoustic pressure, depending on formulation
    PtrCoefFct presFct;
    if (formulation_ == ACOU_POTENTIAL)
    {
      shared_ptr<ResultInfo> pres(new ResultInfo);
      pres->resultType = ACOU_PRESSURE;
      pres->dofNames = "";
      pres->unit = "Pa";
      pres->entryType = ResultInfo::SCALAR;
      pres->definedOn = ResultInfo::ELEMENT;
      // Define pressure as p = rho * dPsi/dt
      PtrCoefFct potD1Fct = this->GetCoefFct(ACOU_POTENTIAL_DERIV_1);
      PtrCoefFct densFct = this->GetCoefFct(ELEM_DENSITY);
      presFct =
          CoefFunction::Generate(mp_, part,
                                 CoefXprBinOp(mp_, potD1Fct, densFct, CoefXpr::OP_MULT));

      // in case of flow, we need the substantial derivative
      //  THIS WOULD WORK IF ONLY CefFunctionBOp WOULD BE DEFINED IN A CONSISTENT MANNER!
      //  USAGE OF THIS COEFFUNCTION in MechPDE and FluidPDE PROHIBIT THIS AT THE MOMENT!
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
    }
    else
    {
      presFct = this->GetCoefFct(ACOU_PRESSURE);
    }

    shared_ptr<ResultInfo> pos, vel, velNormal, intensity, surfIntensity, intensNormal, power, pres;
    PtrCoefFct intensFct, velFct, velFctPW, posFct;
    shared_ptr<CoefFunctionSurf> sNormIntens, sIntens, velFctNormal, sNormIntensPW;
    shared_ptr<CoefFunctionFormBased> presGradFct, velFctPot;
    shared_ptr<ResultFunctor> powerFct;

    PtrCoefFct one = CoefFunction::Generate(mp_, Global::REAL, "1.0");

    // some results are only available in potential and / or
    // harmonic simulation
    if (formulation_ == ACOU_POTENTIAL || isComplex_)
    {
      // === ACOU_PRESSURE ===
      pres.reset(new ResultInfo);
      pres->resultType = ACOU_PRESSURE;
      pres->dofNames = "";
      pres->unit = "Pa";
      pres->entryType = ResultInfo::SCALAR;
      pres->definedOn = ResultInfo::ELEMENT;

      // === ACOU_VELOCITY ===
      vel.reset(new ResultInfo);
      vel->resultType = ACOU_VELOCITY;
      vel->dofNames = vecDofNames;
      vel->unit = "m/s";
      vel->entryType = ResultInfo::VECTOR;
      vel->definedOn = ResultInfo::ELEMENT;

      // === ACOU_POSITION ===
      pos.reset(new ResultInfo);
      pos->resultType = ACOU_POSITION;
      pos->dofNames = vecDofNames;
      pos->unit = "m";
      pos->entryType = ResultInfo::VECTOR;
      pos->definedOn = ResultInfo::ELEMENT;

      // === ACOU_NORMAL_VELOCITY ===
      velNormal.reset(new ResultInfo);
      velNormal->resultType = ACOU_NORMAL_VELOCITY;
      velNormal->dofNames = "";
      velNormal->unit = "m/s";
      velNormal->entryType = ResultInfo::SCALAR;
      velNormal->definedOn = ResultInfo::SURF_ELEM;

      // === ACOU_INTENSITY ===
      intensity.reset(new ResultInfo);
      intensity->resultType = ACOU_INTENSITY;
      intensity->dofNames = vecDofNames;
      intensity->unit = "W/m^2";
      intensity->entryType = ResultInfo::VECTOR;
      intensity->definedOn = ResultInfo::ELEMENT;

      // === ACOU_SURFINTENSITY ===
      surfIntensity.reset(new ResultInfo);
      surfIntensity->resultType = ACOU_SURFINTENSITY;
      surfIntensity->dofNames = vecDofNames;
      surfIntensity->unit = "W/m^2";
      surfIntensity->entryType = ResultInfo::VECTOR;
      surfIntensity->definedOn = ResultInfo::SURF_ELEM;

      // === ACOU_NORMAL_INTENSITY ===
      intensNormal.reset(new ResultInfo);
      intensNormal->resultType = ACOU_NORMAL_INTENSITY;
      intensNormal->dofNames = "";
      intensNormal->unit = "W/m^2";
      intensNormal->entryType = ResultInfo::SCALAR;
      intensNormal->definedOn = ResultInfo::SURF_ELEM;

      // === ACOU_POWER ===
      power.reset(new ResultInfo);
      power->resultType = ACOU_POWER;
      power->dofNames = "";
      power->unit = "W";
      power->entryType = ResultInfo::SCALAR;
      power->definedOn = ResultInfo::SURF_REGION;
    }

    if (formulation_ == ACOU_POTENTIAL)
    {
      // --------------------------
      //  POTENTENTIAL FORMULATION
      // --------------------------
      // === ACOU_VELOCITY ===
      // Velocity v = - grad Psi
      if (isComplex_)
      {
        velFctPot.reset(new CoefFunctionBOp<Complex>(feFct, vel, -1.0));
      }
      else
      {
        velFctPot.reset(new CoefFunctionBOp<Double>(feFct, vel, -1.0));
      }
      velFct = velFctPot;
      stiffFormCoefs_.insert(velFctPot);
      DefineFieldResult(velFct, vel);

      // === ACOU_NORMAL_VELOCITY ===
      velFctNormal.reset(new CoefFunctionSurf(true, 1.0, velNormal));
      DefineFieldResult(velFctNormal, velNormal);
      surfCoefFcts_[velFctNormal] = velFctPot;

      // === ACOU_INTENSITY ===
      // Intensity I = p * conj(v)
      intensFct =
          CoefFunction::Generate(mp_, part,
                                 CoefXprBinOp(mp_, presFct, velFct, CoefXpr::OP_MULT_CONJ));
      DefineFieldResult(intensFct, intensity);

      sIntens.reset(new CoefFunctionSurf(false, 1.0, surfIntensity));
      DefineFieldResult(sIntens, surfIntensity);
      surfCoefFcts_[sIntens] = intensFct;

      // === ACOU_NORMAL_INTENSITY ===
      sNormIntens.reset(new CoefFunctionSurf(true, 1.0, intensNormal));
      DefineFieldResult(sNormIntens, intensNormal);
      surfCoefFcts_[sNormIntens] = intensFct;

      // === ACOU_POWER ===
      // Power p = \int_Gamma I *n dGamma
      // Integrate normal intensity
      if (isComplex_)
      {
        powerFct.reset(new ResultFunctorIntegrate<Complex>(sNormIntens,
                                                           feFct, power));
      }
      else
      {
        powerFct.reset(new ResultFunctorIntegrate<Double>(sNormIntens,
                                                          feFct, power));
      }
      resultFunctors_[ACOU_POWER] = powerFct;
      availResults_.insert(power);
    }

    if (formulation_ == ACOU_PRESSURE && isComplex_)
    {
      // --------------------------------
      //  COMPLEX & PRESSURE FORMULATION
      // --------------------------------

      // === ACOU_VELOCITY ===

      // Velocity v = j* 1/(omega*rho) * grad(p)
      // a) define gradient of pressure
      if (isComplex_)
      {
        presGradFct.reset(new CoefFunctionBOp<Complex>(feFct, vel, 1.0));
      }
      else
      {
        presGradFct.reset(new CoefFunctionBOp<Double>(feFct, vel, 1.0));
      }
      stiffFormCoefs_.insert(presGradFct);

      // b) multiply by factor
      PtrCoefFct densFct = this->GetCoefFct(ELEM_DENSITY);
      PtrCoefFct factor =
          CoefFunction::Generate(mp_, Global::COMPLEX, "0", "1/(2*pi*f)");
      PtrCoefFct factor2 =
          CoefFunction::Generate(mp_, Global::COMPLEX,
                                 CoefXprBinOp(mp_, factor, densFct, CoefXpr::OP_DIV));
      velFct =
          CoefFunction::Generate(mp_, part,
                                 CoefXprBinOp(mp_, factor2, presGradFct, CoefXpr::OP_MULT));
      DefineFieldResult(velFct, vel);

      // === ACOU Particle Position ===
      // u = 1/(rho*omega^2) * grad(p)
      PtrCoefFct oneOverOmega2rho = CoefFunction::Generate(mp_, part,
                                                           CoefXprBinOp(mp_, one,
                                                                        CoefXprBinOp(mp_, CoefFunction::Generate(mp_, Global::REAL, "4*pi*pi*f*f"), densFct, CoefXpr::OP_MULT),
                                                                        CoefXpr::OP_DIV));
      posFct = CoefFunction::Generate(mp_, part, CoefXprBinOp(mp_, oneOverOmega2rho, presGradFct, CoefXpr::OP_MULT));
      DefineFieldResult(posFct, pos);

      // === ACOU_NORMAL_VELOCITY ===
      velFctNormal.reset(new CoefFunctionSurf(true, 1.0, velNormal));
      DefineFieldResult(velFctNormal, velNormal);
      surfCoefFcts_[velFctNormal] = velFct;

      // === ACOU_INTENSITY ===
      // Intensity I = p * v
      intensFct =
          CoefFunction::Generate(mp_, part,
                                 CoefXprBinOp(mp_, presFct, velFct, CoefXpr::OP_MULT_CONJ));

      DefineFieldResult(intensFct, intensity);

      sIntens.reset(new CoefFunctionSurf(false, 1.0, surfIntensity));
      DefineFieldResult(sIntens, surfIntensity);
      surfCoefFcts_[sIntens] = intensFct;

      // === ACOU_NORMAL_INTENSITY ===
      sNormIntens.reset(new CoefFunctionSurf(true, 1.0, intensNormal));
      DefineFieldResult(sNormIntens, intensNormal);
      surfCoefFcts_[sNormIntens] = intensFct;

      // === ACOU_POWER ===
      // Power p = \int_Gamma I *n dGammac
      //  Integrate normal intensity
      powerFct.reset(new ResultFunctorIntegrate<Complex>(sNormIntens,
                                                         feFct, power));
      resultFunctors_[ACOU_POWER] = powerFct;
      availResults_.insert(power);

      // === ACOU_POWER WITH PLANE WAVE ASSUMPTION: p/vn = \rho c ===
      // Power P = \int_Gamma (p^2)/(2*rho*c) dGamma
      // === ACOU_POWER ===
      shared_ptr<ResultInfo> powerPW;
      powerPW.reset(new ResultInfo);
      powerPW->resultType = ACOU_POWER_PLANEWAVE;
      powerPW->dofNames = "";
      powerPW->unit = "W";
      powerPW->entryType = ResultInfo::SCALAR;
      powerPW->definedOn = ResultInfo::SURF_REGION;

      // === ACOU_NORMAL_INTENSITY ===
      shared_ptr<ResultInfo> intensNormalPW;
      intensNormalPW.reset(new ResultInfo);
      intensNormalPW->resultType = ACOU_NORMAL_INTENSITY_PLANEWAVE;
      intensNormalPW->dofNames = "";
      intensNormalPW->unit = "W/m^2";
      intensNormalPW->entryType = ResultInfo::SCALAR;
      intensNormalPW->definedOn = ResultInfo::SURF_ELEM;

      // compute 1/(2*rho*c0)
      // shared_ptr<CoefFunctionMulti> c0Fct  = matCoefs_[ACOU_ELEM_SPEED_OF_SOUND];
      // shared_ptr<CoefFunctionSurf> surfC0(new CoefFunctionSurf(false));
      // surfC0->SetVolumeCoefs( c0Fct->GetRegionCoefs() );

      PtrCoefFct c0Fct = this->GetCoefFct(ACOU_ELEM_SPEED_OF_SOUND);
      PtrCoefFct constVal = CoefFunction::Generate(mp_, Global::REAL, "0.5");
      PtrCoefFct val1 =
          CoefFunction::Generate(mp_, Global::COMPLEX,
                                 CoefXprBinOp(mp_, c0Fct, densFct, CoefXpr::OP_MULT));
      PtrCoefFct val2 =
          CoefFunction::Generate(mp_, Global::COMPLEX,
                                 CoefXprBinOp(mp_, constVal, val1, CoefXpr::OP_DIV));

      velFctPW =
          CoefFunction::Generate(mp_, part,
                                 CoefXprBinOp(mp_, val2, presFct, CoefXpr::OP_MULT));

      PtrCoefFct intensPWfnc =
          CoefFunction::Generate(mp_, Global::COMPLEX,
                                 CoefXprBinOp(mp_, velFctPW, presFct, CoefXpr::OP_MULT_CONJ));

      // normal acoustic intensity for plane waves
      sNormIntensPW.reset(new CoefFunctionSurf(false, 1.0, intensNormalPW));
      DefineFieldResult(sNormIntensPW, intensNormalPW);
      surfCoefFcts_[sNormIntensPW] = intensPWfnc;

      shared_ptr<ResultFunctor> powerFctPW;
      powerFctPW.reset(new ResultFunctorIntegrate<Complex>(sNormIntensPW,
                                                           feFct, powerPW));
      resultFunctors_[ACOU_POWER_PLANEWAVE] = powerFctPW;
      availResults_.insert(powerPW);
    }

    // === ACOUSTIC KINETIC ENERGY ===
    shared_ptr<ResultFunctor> keFunc;
    shared_ptr<ResultInfo> kinEnergy(new ResultInfo);
    kinEnergy->resultType = ACOU_KIN_ENERGY;
    kinEnergy->dofNames = "";
    kinEnergy->unit = "Ws";
    kinEnergy->entryType = ResultInfo::SCALAR;
    kinEnergy->definedOn = ResultInfo::REGION;
    availResults_.insert(kinEnergy);

    shared_ptr<BaseFeFunction> deriv1vFct;
    if (formulation_ == ACOU_PRESSURE)
      deriv1vFct = timeDerivFeFunctions_[ACOU_PRESSURE_DERIV_1];
    else
      deriv1vFct = timeDerivFeFunctions_[ACOU_POTENTIAL_DERIV_1];

    if (isComplex_)
    {
      keFunc.reset(new EnergyResultFunctor<Complex>(deriv1vFct, kinEnergy, 0.5));
    }
    else
    {
      keFunc.reset(new EnergyResultFunctor<Double>(deriv1vFct, kinEnergy, 0.5));
    }
    resultFunctors_[ACOU_KIN_ENERGY] = keFunc;
    massFormFunctors_.insert(keFunc);

    // === ACOUSTIC POTENTIAL ENERGY ===
    shared_ptr<ResultFunctor> keFuncPot;
    shared_ptr<ResultInfo> potEnergy(new ResultInfo);
    potEnergy->resultType = ACOU_POT_ENERGY;
    potEnergy->dofNames = "";
    potEnergy->unit = "Ws";
    potEnergy->entryType = ResultInfo::SCALAR;
    potEnergy->definedOn = ResultInfo::REGION;
    availResults_.insert(potEnergy);
    if (isComplex_)
    {
      keFuncPot.reset(new EnergyResultFunctor<Complex>(feFct, potEnergy, 0.5));
    }
    else
    {
      keFuncPot.reset(new EnergyResultFunctor<Double>(feFct, potEnergy, 0.5));
    }
    resultFunctors_[ACOU_POT_ENERGY] = keFuncPot;
    stiffFormFunctors_.insert(keFuncPot);

    //== ACOUSTIC LOAD DDENSITY  ====
    shared_ptr<ResultInfo> loadDensity(new ResultInfo);
    loadDensity->resultType = ACOU_RHS_LOAD_DENSITY;
    loadDensity->dofNames = "";
    loadDensity->unit = "";

    loadDensity->definedOn = ResultInfo::NODE;
    loadDensity->entryType = ResultInfo::SCALAR;

    acousticSourceDensityCoef_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));
    DefineFieldResult(acousticSourceDensityCoef_, loadDensity);
  }

  void AcousticPDE::CreateMeanFlowFunction(StdVector<std::string> dofNames)
  {
    //// === MEAN FLUIDMECH VELOCITY ===
    shared_ptr<ResultInfo> flowvelocity(new ResultInfo);
    flowvelocity->resultType = MEAN_FLUIDMECH_VELOCITY;
    flowvelocity->dofNames = dofNames;
    flowvelocity->unit = "m/s";

    flowvelocity->definedOn = ResultInfo::NODE;
    flowvelocity->entryType = ResultInfo::VECTOR;

    meanFlowCoef_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_));
    DefineFieldResult(meanFlowCoef_, flowvelocity);

    results_.Push_back(flowvelocity);
    availResults_.insert(flowvelocity);
    if (myParam_->Get("flowFormulation")->As<std::string>() == "withDivergence")
    {
      //// === DIVERGENCE OF MEAN FLOW ===
      shared_ptr<ResultInfo> divflowvelocity(new ResultInfo);
      divflowvelocity->resultType = DIV_MEAN_FLUIDMECH_VELOCITY;
      divflowvelocity->dofNames = "";
      divflowvelocity->unit = "1/s";

      divflowvelocity->definedOn = ResultInfo::ELEMENT;
      divflowvelocity->entryType = ResultInfo::SCALAR;

      divMeanFlowCoef_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));
      DefineFieldResult(divMeanFlowCoef_, divflowvelocity);

      results_.Push_back(divflowvelocity);
      availResults_.insert(divflowvelocity);
    }
  }

  //! Init the time stepping
  void AcousticPDE::InitTimeStepping()
  {

    Double alpha = this->myParam_->Get("timeStepAlpha")->As<Double>();
    GLMScheme *scheme1 = new Newmark(0.5, 0.25, alpha);
    shared_ptr<BaseTimeScheme> acouScheme(new TimeSchemeGLM(scheme1, 0));
    feFunctions_[formulation_]->SetTimeScheme(acouScheme);

    if (this->isTimeDomPML_)
    {
      // the choice of alpha schemes depends on the damping matrix (e.g. PML)
      GLMScheme *scheme2 = new Newmark(0.5, 0.25, alpha);
      shared_ptr<BaseTimeScheme> vecScheme(new TimeSchemeGLM(scheme2, 0));
      feFunctions_[ACOU_PMLAUXVEC]->SetTimeScheme(vecScheme);

      if (!this->isAPML_ && dim_ == 3)
      {
        GLMScheme *scheme3 = new Newmark(0.5, 0.25, alpha);
        shared_ptr<BaseTimeScheme> scalScheme(new TimeSchemeGLM(scheme3, 0));
        feFunctions_[ACOU_PMLAUXSCALAR]->SetTimeScheme(scalScheme);
      }
    }

    if (this->timeDomainEqFluidFormulation_)
    {
      for (unsigned int i = 0; i < nAuxFncAC_.Max(); i++)
      {
        shared_ptr<BaseTimeScheme> schemePhiC(new TimeSchemeGLM(new Newmark(0.5, 0.25, alpha), 0));
        feFunctions_[(SolutionType)(ACOU_TDEF_PHI_C_1 + i)]->SetTimeScheme(schemePhiC);
      }
      for (unsigned int i = 0; i < nAuxFncBC_.Max(); i++)
      {
        shared_ptr<BaseTimeScheme> schemePsiC(new TimeSchemeGLM(new Newmark(0.5, 0.25, alpha), 0));
        feFunctions_[(SolutionType)(ACOU_TDEF_PSI_C_1 + i)]->SetTimeScheme(schemePsiC);
      }
      for (unsigned int i = 0; i < nAuxFncAV_.Max(); i++)
      {
        shared_ptr<BaseTimeScheme> schemePhiV(new TimeSchemeGLM(new Newmark(0.5, 0.25, alpha), 0));
        feFunctions_[(SolutionType)(ACOU_TDEF_PHI_V_1 + i)]->SetTimeScheme(schemePhiV);
      }
      for (unsigned int i = 0; i < nAuxFncBV_.Max(); i++)
      {
        shared_ptr<BaseTimeScheme> schemePsiV(new TimeSchemeGLM(new Newmark(0.5, 0.25, alpha), 0));
        feFunctions_[(SolutionType)(ACOU_TDEF_PSI_V_1 + i)]->SetTimeScheme(schemePsiV);
      }
    }
  }

  void AcousticPDE::ComputeSOS(PtrCoefFct &c0, PtrCoefFct dens, PtrCoefFct blk,
                               PtrCoefFct regionTemp, std::string tempId)
  {
    // TODO 1 erase inconsistent coding style 1/c0^2
    if (complexFluidFormulation_)
    {
      // here c0 is actually 1/c0^2!!
      c0 = CoefFunction::Generate(mp_, Global::COMPLEX,
                                  CoefXprBinOp(mp_, dens, blk, CoefXpr::OP_DIV));
    }
    else
    {
      if (tempId != "")
      {
        if (formulation_ == ACOU_POTENTIAL)
          EXCEPTION("We need for temperature dependent speed of sound a pressure formulation");

        // gasR=287.058 J/kg K   ... universal gas constant
        // kappa=1.402, adabatic exponent for air
        shared_ptr<CoefFunctionConst<Double>> constVal(new CoefFunctionConst<Double>());
        constVal->SetScalar(402.4553160);

        c0 = CoefFunction::Generate(mp_, Global::REAL,
                                    CoefXprUnaryOp(mp_, CoefXprBinOp(mp_, constVal, regionTemp, CoefXpr::OP_MULT),
                                                   CoefXpr::OP_SQRT));
      }
      else
      {
        // c0 = sqrt(bulk_modulus / density)
        c0 = CoefFunction::Generate(mp_, Global::REAL,
                                    CoefXprUnaryOp(mp_, CoefXprBinOp(mp_, blk, dens, CoefXpr::OP_DIV),
                                                   CoefXpr::OP_SQRT));
      }
    }
  }

  void AcousticPDE::ComputeSOS_SQR(PtrCoefFct &cSQR, PtrCoefFct dens, PtrCoefFct blk,
                                   PtrCoefFct regionTemp, std::string tempId)
  {

    if (tempId != "")
    {
      if (formulation_ == ACOU_POTENTIAL)
        EXCEPTION("We need for temperature dependent speed of sound a pressure formulation");

      // gasR=287.058 J/kg K   ... universal gas constant
      // kappa=1.402, adabatic exponent for air
      shared_ptr<CoefFunctionConst<Double>> constVal(new CoefFunctionConst<Double>());
      constVal->SetScalar(402.4553160);

      cSQR = CoefFunction::Generate(mp_, Global::REAL,
                                    CoefXprBinOp(mp_, constVal, regionTemp, CoefXpr::OP_MULT));
    }
    else
    {
      // c^2 = bulk_modulus / density
      cSQR = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, blk, dens, CoefXpr::OP_DIV));
    }
  }


  void AcousticPDE::ReadTDEFCoefficients(UInt iRegion)
  {
    unsigned int actRegion = regions_[iRegion];

    StdVector<PtrCoefFct> fncAC;
    StdVector<PtrCoefFct> fncBC;
    StdVector<PtrCoefFct> fncCC;
    StdVector<PtrCoefFct> fncAlphaC;
    StdVector<PtrCoefFct> fncBetaC;
    StdVector<PtrCoefFct> fncGammaC;

    StdVector<PtrCoefFct> fncAV;
    StdVector<PtrCoefFct> fncBV;
    StdVector<PtrCoefFct> fncCV;
    StdVector<PtrCoefFct> fncAlphaV;
    StdVector<PtrCoefFct> fncBetaV;
    StdVector<PtrCoefFct> fncGammaV;

    LocPointMapped lpm;

    Vector<Double> vecAC;
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

    if (nAuxFncAC_[iRegion] > 0)
    {
      PtrCoefFct coefVecAC = materials_[actRegion]->GetVectorCoefFnc(ACOU_TDEF_INVBLK_A, Global::REAL);
      PtrCoefFct coefVecAlphaC = materials_[actRegion]->GetVectorCoefFnc(ACOU_TDEF_INVBLK_ALPHA, Global::REAL);

      coefVecAC->GetVector(vecAC, lpm);
      coefVecAlphaC->GetVector(vecAlphaC, lpm);
    }
    else
    {
      vecAC.Resize(0);
      vecAlphaC.Resize(0);
    }

    if (nAuxFncBC_[iRegion] > 0)
    {
      PtrCoefFct coefVecBC = materials_[actRegion]->GetVectorCoefFnc(ACOU_TDEF_INVBLK_B, Global::REAL);
      PtrCoefFct coefVecCC = materials_[actRegion]->GetVectorCoefFnc(ACOU_TDEF_INVBLK_C, Global::REAL);
      PtrCoefFct coefVecBetaC = materials_[actRegion]->GetVectorCoefFnc(ACOU_TDEF_INVBLK_BETA, Global::REAL);
      PtrCoefFct coefVecGammaC = materials_[actRegion]->GetVectorCoefFnc(ACOU_TDEF_INVBLK_GAMMA, Global::REAL);

      coefVecBC->GetVector(vecBC, lpm);
      coefVecCC->GetVector(vecCC, lpm);
      coefVecBetaC->GetVector(vecBetaC, lpm);
      coefVecGammaC->GetVector(vecGammaC, lpm);
    }
    else
    {
      vecBC.Resize(0);
      vecCC.Resize(0);
      vecBetaC.Resize(0);
      vecGammaC.Resize(0);
    }

    if (nAuxFncAV_[iRegion] > 0)
    {
      PtrCoefFct coefVecAV = materials_[actRegion]->GetVectorCoefFnc(ACOU_TDEF_INVDENS_A, Global::REAL);
      PtrCoefFct coefVecAlphaV = materials_[actRegion]->GetVectorCoefFnc(ACOU_TDEF_INVDENS_ALPHA, Global::REAL);

      coefVecAV->GetVector(vecAV, lpm);
      coefVecAlphaV->GetVector(vecAlphaV, lpm);
    }
    else
    {
      vecAV.Resize(0);
      vecAlphaV.Resize(0);
    }

    if (nAuxFncBV_[iRegion] > 0)
    {
      PtrCoefFct coefVecBV = materials_[actRegion]->GetVectorCoefFnc(ACOU_TDEF_INVDENS_B, Global::REAL);
      PtrCoefFct coefVecCV = materials_[actRegion]->GetVectorCoefFnc(ACOU_TDEF_INVDENS_C, Global::REAL);
      PtrCoefFct coefVecBetaV = materials_[actRegion]->GetVectorCoefFnc(ACOU_TDEF_INVDENS_BETA, Global::REAL);
      PtrCoefFct coefVecGammaV = materials_[actRegion]->GetVectorCoefFnc(ACOU_TDEF_INVDENS_GAMMA, Global::REAL);

      coefVecBV->GetVector(vecBV, lpm);
      coefVecCV->GetVector(vecCV, lpm);
      coefVecBetaV->GetVector(vecBetaV, lpm);
      coefVecGammaV->GetVector(vecGammaV, lpm);
    }
    else
    {
      vecBV.Resize(0);
      vecCV.Resize(0);
      vecBetaV.Resize(0);
      vecGammaV.Resize(0);
    }

    // Test creat matrix as vector of vector of coef functions
    StdVector<StdVector<PtrCoefFct>> coefTDEFC;
    coefTDEFC = StdVector<StdVector<PtrCoefFct>>();

    // same for coefTDEFV

    fncAC_.Resize(vecAC.GetSize());
    fncBC_.Resize(vecBC.GetSize());
    fncCC_.Resize(vecCC.GetSize());
    fncAlphaC_.Resize(vecAlphaC.GetSize());
    fncBetaC_.Resize(vecBetaC.GetSize());
    fncGammaC_.Resize(vecGammaC.GetSize());

    fncAV_.Resize(vecAV.GetSize());
    fncBV_.Resize(vecBV.GetSize());
    fncCV_.Resize(vecCV.GetSize());
    fncAlphaV_.Resize(vecAlphaV.GetSize());
    fncBetaV_.Resize(vecBetaV.GetSize());
    fncGammaV_.Resize(vecGammaV.GetSize());

    // real poles
    for (UInt ii = 0; ii < vecAC.GetSize(); ii++)
    {
      fncAC_[ii] = CoefFunction::Generate(mp_, Global::REAL, std::to_string(vecAC[ii]));
      std::cout << "Coef AC: " << std::to_string(vecAC[ii]) << std::endl;
    }

    for (UInt ii = 0; ii < vecAlphaC.GetSize(); ii++)
    {
      fncAlphaC_[ii] = CoefFunction::Generate(mp_, Global::REAL, std::to_string(vecAlphaC[ii]));
      std::cout << "Coef AlphaC: " << std::to_string(vecAlphaC[ii]) << std::endl;
    }

    for (UInt ii = 0; ii < vecAV.GetSize(); ii++)
    {
      fncAV_[ii] = CoefFunction::Generate(mp_, Global::REAL, std::to_string(vecAV[ii]));
      std::cout << "Coef AV: " << std::to_string(vecAV[ii]) << std::endl;
    }

    for (UInt ii = 0; ii < vecAlphaV.GetSize(); ii++)
    {
      fncAlphaV_[ii] = CoefFunction::Generate(mp_, Global::REAL, std::to_string(vecAlphaV[ii]));
      std::cout << "Coef AlphaV: " << std::to_string(vecAlphaV[ii]) << std::endl;
    }

    // Complex Poles
    for (UInt ii = 0; ii < vecBC.GetSize(); ii++)
    {
      fncBC_[ii] = CoefFunction::Generate(mp_, Global::REAL, std::to_string(vecBC[ii]));
      std::cout << "Coef BC: " << std::to_string(vecBC[ii]) << std::endl;
    }
    for (UInt ii = 0; ii < vecCC.GetSize(); ii++)
    {
      fncCC_[ii] = CoefFunction::Generate(mp_, Global::REAL, std::to_string(vecCC[ii]));
      std::cout << "Coef CC: " << std::to_string(vecCC[ii]) << std::endl;
    }
    for (UInt ii = 0; ii < vecBetaC.GetSize(); ii++)
    {
      fncBetaC_[ii] = CoefFunction::Generate(mp_, Global::REAL, std::to_string(vecBetaC[ii]));
      std::cout << "Coef BetaC: " << std::to_string(vecBetaC[ii]) << std::endl;
    }
    for (UInt ii = 0; ii < vecGammaC.GetSize(); ii++)
    {
      fncGammaC_[ii] = CoefFunction::Generate(mp_, Global::REAL, std::to_string(vecGammaC[ii]));
      std::cout << "Coef GammaC: " << std::to_string(vecGammaC[ii]) << std::endl;
    }

    for (UInt ii = 0; ii < vecBV.GetSize(); ii++)
    {
      fncBV_[ii] = CoefFunction::Generate(mp_, Global::REAL, std::to_string(vecBV[ii]));
      std::cout << "Coef BV: " << std::to_string(vecBV[ii]) << std::endl;
    }
    for (UInt ii = 0; ii < vecCV.GetSize(); ii++)
    {
      fncCV_[ii] = CoefFunction::Generate(mp_, Global::REAL, std::to_string(vecCV[ii]));
      std::cout << "Coef CV: " << std::to_string(vecCV[ii]) << std::endl;
    }
    for (UInt ii = 0; ii < vecBetaV.GetSize(); ii++)
    {
      fncBetaV_[ii] = CoefFunction::Generate(mp_, Global::REAL, std::to_string(vecBetaV[ii]));
      std::cout << "Coef BetaV: " << std::to_string(vecBetaV[ii]) << std::endl;
    }
    for (UInt ii = 0; ii < vecGammaV.GetSize(); ii++)
    {
      fncGammaV_[ii] = CoefFunction::Generate(mp_, Global::REAL, std::to_string(vecGammaV[ii]));
      std::cout << "Coef GammaV: " << std::to_string(vecGammaV[ii]) << std::endl;
    }

    // check if all coefFunctions vectors are smaller than 15
    if (fncAC_.GetSize() > 15 || fncBC_.GetSize() > 15 || fncCC_.GetSize() > 15 || fncAlphaC_.GetSize() > 15 || fncBetaC_.GetSize() > 15 || fncGammaC_.GetSize() > 15 ||
        fncAV_.GetSize() > 15 || fncBV_.GetSize() > 15 || fncCV_.GetSize() > 15 || fncAlphaV_.GetSize() > 15 || fncBetaV_.GetSize() > 15 || fncGammaV_.GetSize() > 15)
    {
      EXCEPTION("TDEF: Only 15 ples are allowed, please reduce the number of poles!");
    }

    std::cout << std::endl;
    std::cout << "TDEF parameters: " << std::endl;
    std::cout << "Rational model of compressibility (inverse bulk modulus): " << std::endl;
    std::cout << "size fncAC: " << fncAC_.GetSize() << std::endl;
    std::cout << "size fncAlphaC: " << fncAlphaC_.GetSize() << std::endl;
    std::cout << "size fncBC: " << fncBC_.GetSize() << std::endl;
    std::cout << "size fncBetaC: " << fncBetaC_.GetSize() << std::endl;
    std::cout << "size fncCC: " << fncCC_.GetSize() << std::endl;
    std::cout << "size fncGammaC: " << fncGammaC_.GetSize() << std::endl;

    std::cout << "Rational model of specific volume (inverse density): " << std::endl;
    std::cout << "size fncAv: " << fncAV_.GetSize() << std::endl;
    std::cout << "size fncAlphav: " << fncAlphaV_.GetSize() << std::endl;
    std::cout << "size fncBv: " << fncBV_.GetSize() << std::endl;
    std::cout << "size fncBetav: " << fncBetaV_.GetSize() << std::endl;
    std::cout << "size fncCv: " << fncCV_.GetSize() << std::endl;
    std::cout << "size fncGammav: " << fncGammaV_.GetSize() << std::endl;

    std::cout << "Finished reading TDEF coefficients." << std::endl
              << std::endl;
  }


  // void AcousticPDE::UpdateFreq(){
  // freq_ = this->mp_->Eval(mHandle_);
  // }


  void AcousticPDE::EvalRationalFncs(UInt iRegion, Double ftrg)
  {
    ReadTDEFCoefficients(iRegion);

    PtrCoefFct imagUnit = CoefFunction::Generate(mp_, Global::COMPLEX, "0.0", "1.0");
    PtrCoefFct constOne = CoefFunction::Generate(mp_, Global::REAL, "1.0");

    PtrCoefFct targetFreg = CoefFunction::Generate(mp_, Global::REAL, std::to_string(ftrg));
    PtrCoefFct omegaTrg = CoefFunction::Generate(mp_, Global::REAL,
                                                 CoefXprBinOp(mp_, targetFreg, CoefFunction::Generate(mp_, Global::REAL, "2*pi"), CoefXpr::OP_MULT));
    PtrCoefFct omegaIm = CoefFunction::Generate(mp_, Global::COMPLEX,
                                                CoefXprBinOp(mp_, omegaTrg, imagUnit, CoefXpr::OP_MULT));

    PtrCoefFct fncTerm;
    PtrCoefFct fncDenom;

    PtrCoefFct fncImag;
    PtrCoefFct fncNumerator;

    // initialize the sums with the high frequency limits
    unsigned int aRegion = regions_[iRegion];
    invTDEFBlk_ = materials_[aRegion]->GetScalCoefFnc(ACOU_TDEF_INVBLK_CONST, Global::REAL);
    invTDEFDens_ = materials_[aRegion]->GetScalCoefFnc(ACOU_TDEF_INVDENS_CONST, Global::REAL);

    // TODO: We have to compute alpha + i omega instead alpha - i omega
    // to get the correct sign for the imaginary part

    // sum over real poles of inverse bulk modulus
    for (UInt ii = 0; ii < fncAC_.GetSize(); ii++)
    {
      fncDenom = CoefFunction::Generate(mp_, Global::COMPLEX,
                                        CoefXprBinOp(mp_, fncAlphaC_[ii], omegaIm, CoefXpr::OP_ADD));

      fncTerm = CoefFunction::Generate(mp_, Global::COMPLEX,
                                       CoefXprBinOp(mp_, fncAC_[ii], fncDenom, CoefXpr::OP_DIV));

      invTDEFBlk_ = CoefFunction::Generate(mp_, Global::COMPLEX,
                                           CoefXprBinOp(mp_, invTDEFBlk_, fncTerm, CoefXpr::OP_ADD));

      std::cout << "added real pole: invTDEFBlk_ = " << invTDEFBlk_->ToString() << "\n";
    }

    // sum over complex poles of inverse bulk modulus
    for (UInt ii = 0; ii < fncBC_.GetSize(); ii++)
    {

      // pole 1 of complex conjugated pole pair
      fncImag = CoefFunction::Generate(mp_, Global::COMPLEX,
                                       CoefXprBinOp(mp_, fncCC_[ii], imagUnit, CoefXpr::OP_MULT));

      fncNumerator = CoefFunction::Generate(mp_, Global::COMPLEX,
                                            CoefXprBinOp(mp_, fncBC_[ii], fncImag, CoefXpr::OP_ADD));

      fncImag = CoefFunction::Generate(mp_, Global::COMPLEX,
                                       CoefXprBinOp(mp_, fncGammaC_[ii], imagUnit, CoefXpr::OP_MULT));

      fncDenom = CoefFunction::Generate(mp_, Global::COMPLEX,
                                        CoefXprBinOp(mp_, fncBetaC_[ii], fncImag, CoefXpr::OP_ADD));

      fncDenom = CoefFunction::Generate(mp_, Global::COMPLEX,
                                        CoefXprBinOp(mp_, fncDenom, omegaIm, CoefXpr::OP_ADD));

      fncTerm = CoefFunction::Generate(mp_, Global::COMPLEX,
                                       CoefXprBinOp(mp_, fncNumerator, fncDenom, CoefXpr::OP_DIV));

      invTDEFBlk_ = CoefFunction::Generate(mp_, Global::COMPLEX,
                                           CoefXprBinOp(mp_, invTDEFBlk_, fncTerm, CoefXpr::OP_ADD));

      // pole 2 of complex conjugated pole pair
      fncImag = CoefFunction::Generate(mp_, Global::COMPLEX,
                                       CoefXprBinOp(mp_, fncCC_[ii], imagUnit, CoefXpr::OP_MULT));

      fncNumerator = CoefFunction::Generate(mp_, Global::COMPLEX,
                                            CoefXprBinOp(mp_, fncBC_[ii], fncImag, CoefXpr::OP_SUB));

      fncImag = CoefFunction::Generate(mp_, Global::COMPLEX,
                                       CoefXprBinOp(mp_, fncGammaC_[ii], imagUnit, CoefXpr::OP_MULT));

      fncDenom = CoefFunction::Generate(mp_, Global::COMPLEX,
                                        CoefXprBinOp(mp_, fncBetaC_[ii], fncImag, CoefXpr::OP_SUB));

      fncDenom = CoefFunction::Generate(mp_, Global::COMPLEX,
                                        CoefXprBinOp(mp_, fncDenom, omegaIm, CoefXpr::OP_ADD));

      fncTerm = CoefFunction::Generate(mp_, Global::COMPLEX,
                                       CoefXprBinOp(mp_, fncNumerator, fncDenom, CoefXpr::OP_DIV));

      invTDEFBlk_ = CoefFunction::Generate(mp_, Global::COMPLEX,
                                           CoefXprBinOp(mp_, invTDEFBlk_, fncTerm, CoefXpr::OP_ADD));

      std::cout << "added complex pole: invTDEFBlk_ = " << invTDEFBlk_->ToString() << "\n";
    }

    // sum over real poles of inverse density
    for (UInt ii = 0; ii < fncAV_.GetSize(); ii++)
    {
      fncDenom = CoefFunction::Generate(mp_, Global::COMPLEX,
                                        CoefXprBinOp(mp_, fncAlphaV_[ii], omegaIm, CoefXpr::OP_ADD));

      fncTerm = CoefFunction::Generate(mp_, Global::COMPLEX,
                                       CoefXprBinOp(mp_, fncAV_[ii], fncDenom, CoefXpr::OP_DIV));

      invTDEFDens_ = CoefFunction::Generate(mp_, Global::COMPLEX,
                                            CoefXprBinOp(mp_, invTDEFDens_, fncTerm, CoefXpr::OP_ADD));

      std::cout << "added real pole: invTDEFDens_ = " << invTDEFDens_->ToString() << "\n";
    }

    // sum over complex poles of inverse density
    for (UInt ii = 0; ii < fncBV_.GetSize(); ii++)
    {

      // pole 1 of complex conjugated pole pair
      fncImag = CoefFunction::Generate(mp_, Global::COMPLEX,
                                       CoefXprBinOp(mp_, fncCV_[ii], imagUnit, CoefXpr::OP_MULT));

      fncNumerator = CoefFunction::Generate(mp_, Global::COMPLEX,
                                            CoefXprBinOp(mp_, fncBV_[ii], fncImag, CoefXpr::OP_ADD));

      fncImag = CoefFunction::Generate(mp_, Global::COMPLEX,
                                       CoefXprBinOp(mp_, fncGammaV_[ii], imagUnit, CoefXpr::OP_MULT));

      fncDenom = CoefFunction::Generate(mp_, Global::COMPLEX,
                                        CoefXprBinOp(mp_, fncBetaV_[ii], fncImag, CoefXpr::OP_ADD));

      fncDenom = CoefFunction::Generate(mp_, Global::COMPLEX,
                                        CoefXprBinOp(mp_, fncDenom, omegaIm, CoefXpr::OP_ADD));

      fncTerm = CoefFunction::Generate(mp_, Global::COMPLEX,
                                       CoefXprBinOp(mp_, fncNumerator, fncDenom, CoefXpr::OP_DIV));

      invTDEFDens_ = CoefFunction::Generate(mp_, Global::COMPLEX,
                                            CoefXprBinOp(mp_, invTDEFDens_, fncTerm, CoefXpr::OP_ADD));

      // pole 2 of complex conjugated pole pair
      fncImag = CoefFunction::Generate(mp_, Global::COMPLEX,
                                       CoefXprBinOp(mp_, fncCV_[ii], imagUnit, CoefXpr::OP_MULT));

      fncNumerator = CoefFunction::Generate(mp_, Global::COMPLEX,
                                            CoefXprBinOp(mp_, fncBV_[ii], fncImag, CoefXpr::OP_SUB));

      fncImag = CoefFunction::Generate(mp_, Global::COMPLEX,
                                       CoefXprBinOp(mp_, fncGammaV_[ii], imagUnit, CoefXpr::OP_MULT));

      fncDenom = CoefFunction::Generate(mp_, Global::COMPLEX,
                                        CoefXprBinOp(mp_, fncBetaV_[ii], fncImag, CoefXpr::OP_SUB));

      fncDenom = CoefFunction::Generate(mp_, Global::COMPLEX,
                                        CoefXprBinOp(mp_, fncDenom, omegaIm, CoefXpr::OP_ADD));

      fncTerm = CoefFunction::Generate(mp_, Global::COMPLEX,
                                       CoefXprBinOp(mp_, fncNumerator, fncDenom, CoefXpr::OP_DIV));

      invTDEFDens_ = CoefFunction::Generate(mp_, Global::COMPLEX,
                                            CoefXprBinOp(mp_, invTDEFDens_, fncTerm, CoefXpr::OP_ADD));

      std::cout << "added complex pole: invTDEFDens_ = " << invTDEFDens_->ToString() << "\n\n";
    }

    PtrCoefFct dens = CoefFunction::Generate(mp_, Global::COMPLEX,
                                             CoefXprBinOp(mp_, constOne, invTDEFDens_, CoefXpr::OP_DIV));
    PtrCoefFct Blk = CoefFunction::Generate(mp_, Global::COMPLEX,
                                            CoefXprBinOp(mp_, constOne, invTDEFBlk_, CoefXpr::OP_DIV));

    std::cout << "final: density = " << dens->ToString() << "\n";
    std::cout << "final: blk mod = " << Blk->ToString() << "\n\n";
  }

}

template void AcousticPDE::DefineTransientPMLInts<2>(shared_ptr<ElemList>, std::string,
                                                     RegionIdType actRegion, std::string tempId);
template void AcousticPDE::DefineTransientPMLInts<3>(shared_ptr<ElemList>, std::string,
                                                     RegionIdType actRegion, std::string tempId);
