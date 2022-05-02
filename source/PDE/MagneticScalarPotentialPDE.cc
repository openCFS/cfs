#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <set>

#include "MagneticScalarPotentialPDE.hh"

#include "General/defs.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Domain/CoefFunction/CoefFunction.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/CoefFunction/CoefFunctionMapping.hh"
#include "Utils/StdVector.hh"
#include "Driver/Assemble.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"
#include "FeBasis/H1/FeSpaceH1Nodal.hh"
#include "FeBasis/FeFunctions.hh"

// new integrator concept
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/SingleEntryInt.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/IdentityOperatorNormal.hh"
#include "Driver/SolveSteps/SolveStepElec.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"
#include "Domain/CoefFunction/CoefFunctionMulti.hh"

namespace CoupledField
{

  DEFINE_LOG(magscalpde, "magneticScalarPotentialPDE")

  // ***************
  //   Constructor
  // ***************
  MagneticScalarPotentialPDE::MagneticScalarPotentialPDE(Grid *aptgrid, PtrParamNode paramNode,
                                                         PtrParamNode infoNode,
                                                         shared_ptr<SimState> simState, Domain *domain)
      : SinglePDE(aptgrid, paramNode, infoNode, simState, domain)
  {
    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_ = "magnetic";
    pdematerialclass_ = ELECTROMAGNETIC;

    nonLin_ = false;
    nonLinMaterial_ = false;

    //! Always use updated Lagrangian formulation
    updatedGeo_ = true;

    // this is the default value
    modelName_ = "nonlinearCurve";

    // Be aware that the hysteresis which are build in via CoefFunctionMaterialModel does not use the isHysteresis_ flag!
    matModelCoef_.reset(new CoefFunctionMaterialModel<Complex>());
    // init the nlFluxCoef
    nlFluxCoef_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_, true));

    //perm_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, dim_, dim_, isComplex_));
  }

  void MagneticScalarPotentialPDE::DefineIntegrators()
  {

    RegionIdType actRegion;
    BaseMaterial *actSDMat = NULL;

    // get FEFunctions and space
    shared_ptr<BaseFeFunction> feFunc_reduced = feFunctions_[MAG_POTENTIAL];
    shared_ptr<FeSpace> feSpace_reduced = feFunc_reduced->GetFeSpace();

    // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial *>::iterator it;

    
    // Init material model for transient analysis
    // if (((analysistype_ == STATIC) || (analysistype_ == TRANSIENT)) && nonLin_ && (modelName_ != "nonlinearCurve"))
    // {
    //   matModelCoef_->Init(magFieldCoef, modelName_);
    // }

    // iterate over the region (or materials)
    for (UInt iRegion = 0; iRegion < regions_.GetSize(); iRegion++)
    {
      Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;

      // set current region and material
      actRegion = regions_[iRegion];
      actSDMat = materials_[actRegion];
      StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion];

      // Get current region name
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);

      // create new entity list
      shared_ptr<ElemList> actSDList(new ElemList(ptGrid_)); // changed newSDList to actSDList
      /* actSDList = newSDList; */
      actSDList->SetRegion(actRegion);

      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region", "name", regionName.c_str());
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();

      // --- Set the approximation for the current region ---
      feSpace_reduced->SetRegionApproximation(actRegion, polyId, integId);
      feFunc_reduced->AddEntityList(actSDList);

      PtrCoefFct magFieldCoef = this->GetCoefFct(MAG_FIELD_INTENSITY);

      // ====================================================================
      // stiffness integrator - NONLINEAR
      // ====================================================================
      BaseBDBInt *stiffInt = NULL;
      PtrCoefFct muNL = NULL;
      if (nonLinTypes.Find(PERMEABILITY) != -1)
      {
        // ====================================================================
        // stiffness integrator - NONLINEAR
        // ====================================================================

        // hysteretic - Energy Based
        if (modelName_ == "JilesAthertonModel")
        {
          EXCEPTION("Jiles-Atherton model not implemented for MagneticScalarPotentialPDE");
        }
        else if (modelName_ == "EBHysteresis")
        {
          EXCEPTION("EB Hysteresis model not yet finished for MagneticScalarPotentialPDE");
          std::map<std::string, double> ParameterMap;
          actSDMat->GetScalar(ParameterMap["Ps"], MAG_PS_EB, Global::REAL);
          actSDMat->GetScalar(ParameterMap["A"], MAG_A_EB, Global::REAL);
          actSDMat->GetScalar(ParameterMap["mu0"], MAG_MU0_EB, Global::REAL);
          actSDMat->GetScalar(ParameterMap["numS"], MAG_NUMS_EB, Global::REAL);
          actSDMat->GetScalar(ParameterMap["chi_factor"], MAG_CHI_FACTOR_EB, Global::REAL);
          ParameterMap["isMH"] = 0;
          matModelCoef_->InitModel(ParameterMap, actSDList->GetSize());

          muNL = actSDMat->GetTensorCoefFncModel(matModelCoef_);

          // TODO kroppert: create functionality to get the B_previous in the same way as muNL_...the values are there
          //                in the hysteresis model, we only need to pipe it through here
          PtrCoefFct hystFluxTmp = CoefFunction::Generate(mp_, part, CoefXprBinOp(mp_, muNL, magFieldCoef, CoefXpr::OP_MULT));
          nlFluxCoef_->AddRegion(actRegion, hystFluxTmp);

          if (dim_ == 2)
          {
            stiffInt = new BDBInt<>(new GradientOperator<FeH1, 2>(), muNL, 1.0, updatedGeo_);
          }
          else
          {
            stiffInt = new BDBInt<>(new GradientOperator<FeH1, 3>(), muNL, 1.0, updatedGeo_);
          }
        }
        else
        {
          // The problem is, that we need the nonlinearity as a function of H and not of B, like in the usual
          // magnetic nonlinear cases!!

          // classical nonlinearity with analytic prescription
          PtrCoefFct muNl = actSDMat->GetScalCoefFncNonLin( MAG_PERMEABILITY_SCALAR, Global::REAL, magFieldCoef);
          
          if (dim_ == 2)
          {
            stiffInt = new BBInt<>(new GradientOperator<FeH1, 2>(), muNl, 1.0, updatedGeo_);
          }
          else
          {
            stiffInt = new BBInt<>(new GradientOperator<FeH1, 3>(), muNl, 1.0, updatedGeo_);
          }
          matCoefs_[MAG_ELEM_PERMEABILITY]->AddRegion(actRegion, muNl);
//          perm_->AddRegion(actRegion, muNl);
        }

        stiffInt->SetName("StiffnessIntegratorHysteresis");
      }
      else
      {
        // ====================================================================
        // stiffness integrator - LINEAR
        // ====================================================================
        PtrCoefFct mu = actSDMat->GetScalCoefFnc(MAG_PERMEABILITY_SCALAR, Global::REAL);
        //matCoefs_[MAG_ELEM_PERMEABILITY]->AddRegion(actRegion, mu);
        if (dim_ == 2)
        {
          stiffInt = new BBInt<>(new GradientOperator<FeH1, 2>(), mu, 1.0, updatedGeo_);
        }
        else
        {
          stiffInt = new BBInt<>(new GradientOperator<FeH1, 3>(), mu, 1.0, updatedGeo_);
        }
        stiffInt->SetName("StiffnessIntegrator");
        matCoefs_[MAG_ELEM_PERMEABILITY]->AddRegion(actRegion, mu);
      }

      BiLinFormContext *stiffIntDescr = new BiLinFormContext(stiffInt, STIFFNESS);

      stiffIntDescr->SetEntities(actSDList, actSDList);
      stiffIntDescr->SetFeFunctions(feFunc_reduced, feFunc_reduced);
      stiffInt->SetFeSpace(feSpace_reduced);

      assemble_->AddBiLinearForm(stiffIntDescr);
      bdbInts_.insert(std::pair<RegionIdType, BaseBDBInt *>(actRegion, stiffInt));
    }
  }

  void MagneticScalarPotentialPDE::DefineRhsLoadIntegrators()
  {

    // get FEFunctions and space
    shared_ptr<BaseFeFunction> feFunc_reduced = feFunctions_[MAG_POTENTIAL];
    shared_ptr<FeSpace> feSpace_reduced = feFunc_reduced->GetFeSpace();

    StdVector<shared_ptr<EntityList>> ent;
    StdVector<PtrCoefFct> coef;
    LinearForm *lin = NULL;
    StdVector<std::string> dofNames;

    std::set<RegionIdType> volRegions (regions_.Begin(), regions_.End() );

    bool coefUpdateGeo = true;
    // =====================
    //  Current Vector Potential (T or Hs)
    // =====================
    ReadRhsExcitation("fieldIntensity", dofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo);
    for (UInt i = 0; i < ent.GetSize(); ++i)
    {
      UInt entDim = ent[i]->GetGrid()->GetEntityDim(ptGrid_->GetRegionName(ent[i]->GetRegion()));
      std::string regName = ptGrid_->GetRegionName(ent[i]->GetRegion());
      std::string entName = ent[i]->GetName();
      LOG_DBG(magscalpde) << "Entity with name " << entName
                          << " is in region with name "
                          << entName
                          << "\n and has the spatial dimension "
                          << entDim;

      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST)
      {
        EXCEPTION("Magnetic field intensity must be defined on elements")
      }
      if (regName.compare(entName) != 0)
      {
        EXCEPTION("There seems to be an error with region and entity names")
      }
      

      // Here we store the Hs field for every region to have it ready for postprocessing
      Hsmap_[ent[i]->GetRegion()] = coef[i];

      lin = new BUIntegrator<Double>(new GradientOperator<FeH1, 3>(), 1.0, coef[i], volRegions, coefUpdateGeo);
      
      lin->SetName("SourceMagFieldIntensityInt");
      LinearFormContext *ctx = new LinearFormContext(lin);
      ctx->SetEntities(ent[i]);
      ctx->SetFeFunction(feFunc_reduced);
      assemble_->AddLinearForm(ctx);

    } // end loop over entities
  }


  // ****************************
  //  Initialize Nonlinearities
  // ****************************
  void MagneticScalarPotentialPDE::InitNonLin() {

    SinglePDE::InitNonLin();
 //  nonLinTotalFormulation_ = true;

  }

  void MagneticScalarPotentialPDE::DefineSolveStep()
  {
    solveStep_ = new StdSolveStep(*this);
  }

  void MagneticScalarPotentialPDE::InitTimeStepping()
  {
    Double gamma = 1.0;
    GLMScheme *scheme = new Trapezoidal(gamma);

    TimeSchemeGLM::NonLinType nlType = (nonLin_ || isHysteresis_) ? TimeSchemeGLM::INCREMENTAL : TimeSchemeGLM::NONE;
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(scheme, 0, nlType));
    feFunctions_[MAG_POTENTIAL]->SetTimeScheme(myScheme);
  }

  void MagneticScalarPotentialPDE::DefinePrimaryResults()
  {

    shared_ptr<BaseFeFunction> feFct = feFunctions_[MAG_POTENTIAL];

    // Magnetic Potential
    shared_ptr<ResultInfo> res1(new ResultInfo);
    res1->resultType = MAG_POTENTIAL;
    res1->dofNames = "";
    res1->unit = "A";
    res1->definedOn = ResultInfo::NODE;
    res1->entryType = ResultInfo::SCALAR;
    feFunctions_[MAG_POTENTIAL]->SetResultInfo(res1);
    res1->SetFeFunction(feFunctions_[MAG_POTENTIAL]);
    results_.Push_back(res1);
    availResults_.insert(res1);
    DefineFieldResult(feFunctions_[MAG_POTENTIAL], res1);

    // === PERMEABILITY  ===
    shared_ptr<ResultInfo> perm(new ResultInfo);
    perm->resultType = MAG_ELEM_PERMEABILITY;
    perm->dofNames = "";
    perm->unit = "Vs/Am";
    perm->definedOn = ResultInfo::ELEMENT;
    perm->entryType = ResultInfo::SCALAR;
    shared_ptr<CoefFunctionMulti> permFct(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1, false));
    matCoefs_[MAG_ELEM_PERMEABILITY] = permFct;
    DefineFieldResult( permFct, perm );


    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    hdbcSolNameMap_[MAG_POTENTIAL] = "ground";
    idbcSolNameMap_[MAG_POTENTIAL] = "potential";
  }

  void MagneticScalarPotentialPDE::DefinePostProcResults()
  {

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

    shared_ptr<BaseFeFunction> feFct = feFunctions_[MAG_POTENTIAL];

    // === MAGNETIC FIELD INTENSITY ===
    shared_ptr<ResultInfo> ef(new ResultInfo);
    ef->resultType = MAG_FIELD_INTENSITY;
    ef->dofNames = vecDofNames;
    ef->unit = "A/m";
    ef->definedOn = ResultInfo::ELEMENT;
    ef->entryType = ResultInfo::VECTOR;
    // The recipe about how to actually evaluate H, is defined in FinalizePostProcResults()
    shared_ptr<CoefFunctionMulti> hIntensFunc(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_));
    DefineFieldResult(hIntensFunc, ef);

    // === MAGNETIC FLUX DENSITY ===
    shared_ptr<ResultInfo> bf(new ResultInfo);
    bf->resultType = MAG_FLUX_DENSITY;
    bf->dofNames = vecDofNames;
    bf->unit = "T";
    bf->definedOn = ResultInfo::ELEMENT;
    bf->entryType = ResultInfo::VECTOR;
    // The recipe about how to actually evaluate B, is defined in FinalizePostProcResults()
    shared_ptr<CoefFunctionMulti> bFunc(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_));
    DefineFieldResult(bFunc, bf);

    // === - GRADIENT PHI (helper result)===
    shared_ptr<ResultInfo> hf(new ResultInfo);
    hf->resultType = MAG_POTENTIAL_DIV;
    hf->dofNames = vecDofNames;
    hf->unit = "A/m";
    hf->definedOn = ResultInfo::ELEMENT;
    hf->entryType = ResultInfo::VECTOR;
    shared_ptr<CoefFunctionFormBased> hFunc;
    if (isComplex_)
    {
      hFunc.reset(new CoefFunctionBOp<Complex>(feFct, hf, -1.0));
    }
    else
    {
      hFunc.reset(new CoefFunctionBOp<Double>(feFct, hf, -1.0));
    }
    DefineFieldResult(hFunc, hf);
    stiffFormCoefs_.insert(hFunc);
  }

  void MagneticScalarPotentialPDE::FinalizePostProcResults()
  {
    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;


    // Initialize standard postprocessing results
    SinglePDE::FinalizePostProcResults();

    // ============ MAGNETIC FIELD INTENSITY ============
    // Assemble coefficient function for
    // H = Hs - \grad \Phi
    // whereas Hs is stored in a CoefFunction and typically read in from a previous sequence
    // step as the solution of a curl curl type problem (magnetostatics)

    // Hs
    shared_ptr<CoefFunctionMulti> hIntensCoef = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_FIELD_INTENSITY]);
    shared_ptr<CoefFunctionMulti> bCoef = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_FLUX_DENSITY]);
    shared_ptr<CoefFunctionMulti> mu = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_ELEM_PERMEABILITY]);

    StdVector<RegionIdType>::iterator regIt = regions_.Begin();
    regIt = regions_.Begin();
    for (; regIt != regions_.End(); ++regIt)
    {
      // ========= H field =============
      PtrCoefFct hs, magIntens;
      if(Hsmap_.find(*regIt) != Hsmap_.end()){
        // Hs is prescribed in the region
        hs = Hsmap_[*regIt];
        // result MAG_POTENTIAL_DIV is (-grad(Phi))!!
        magIntens = CoefFunction::Generate(mp_, part, CoefXprBinOp(mp_, hs, GetCoefFct(MAG_POTENTIAL_DIV), CoefXpr::OP_ADD));        
        hIntensCoef->AddRegion(*regIt, magIntens);
      }else{
        // Hs is NOT prescribed in the region
        hIntensCoef->AddRegion(*regIt, GetCoefFct(MAG_POTENTIAL_DIV));
      }

      // ========= B field =============
      PtrCoefFct b = CoefFunction::Generate( mp_, part, CoefXprVecScalOp(mp_, hIntensCoef, mu, CoefXpr::OP_MULT));
      bCoef->AddRegion(*regIt, b);

    } // loop over regions
  }

  // create the FEspace (H1 in this case)
  std::map<SolutionType, shared_ptr<FeSpace>>
  MagneticScalarPotentialPDE::CreateFeSpaces(const std::string &formulation, PtrParamNode infoNode)
  {
    std::map<SolutionType, shared_ptr<FeSpace>> crSpaces;
    if (formulation == "default" || formulation == "H1")
    {
      PtrParamNode potSpaceNode_reduced = infoNode->Get("magPotential");
      crSpaces[MAG_POTENTIAL] = FeSpace::CreateInstance(myParam_, potSpaceNode_reduced, FeSpace::H1, ptGrid_);
      crSpaces[MAG_POTENTIAL]->Init(solStrat_);
    }
    else
    {
      EXCEPTION("The formulation " << formulation << "of magnetic PDE is not known!");
    }
    return crSpaces;
  }
}