#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <set>

#include "ElecCurrentPDE.hh"

#include "General/defs.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Domain/CoefFunction/CoefFunction.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Utils/StdVector.hh"
#include "Driver/Assemble.hh"
#include "Utils/ApproxData.hh"
#include "Utils/SmoothSpline.hh"
#include "Materials/Models/Hysteresis.hh"
#include "Materials/Models/Preisach.hh"
#include "FeBasis/H1/FeSpaceH1Nodal.hh"
#include "FeBasis/FeFunctions.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"


// new integrator concept
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/SingleEntryInt.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/LinForms/KXInt.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/IdentityOperatorNormal.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"


namespace CoupledField {

  DEFINE_LOG(eleccurrentpde, "pde.electricCurrent")

  // ***************
  //   Constructor
  // ***************
  ElecCurrentPDE::ElecCurrentPDE( Grid* aptgrid, PtrParamNode paramNode,
                    PtrParamNode infoNode,
                    shared_ptr<SimState> simState, Domain* domain )
    :SinglePDE( aptgrid, paramNode, infoNode, simState, domain ) {


    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_          = "elecConduction";
    pdematerialclass_ = ELECTRICCONDUCTION;
 
    nonLin_         = false;
    nonLinMaterial_ = false;
    
    //! Always use updated Lagrangian formulation 
    updatedGeo_     = true;
  }
  
  void ElecCurrentPDE::InitNonLin() {

    SinglePDE::InitNonLin();

  }

  void ElecCurrentPDE::GetPoleRegionIds(NonLinType nl, StdVector<RegionIdType> & regIds) {

    // TODO: also check pole id
    if ( !myParam_->Has("poleList") ) 
      EXCEPTION("B/-Tripole was defined but no xml section 'poleList' is present.");

    StdVector<std::string> poleNames, regNames;
    StdVector<RegionIdType> volRegions;
    PtrParamNode xml;
    if ( nl == NLELEC_BIPOLE || nl == NLELEC_BIPOLE_TEMP_DEP) {

      poleNames.Push_back ("anode");
      poleNames.Push_back ("cathode");

      regNames.Resize(2);
      volRegions.Resize(2);
      regIds.Resize(2);
      
      xml = myParam_->Get("poleList")->Get("Bipole"); 
    }
    else if ( nl == NLELEC_TRIPOLE || nl == NLELEC_TRIPOLE_TEMP_DEP) {
      poleNames.Push_back ("gate");
      poleNames.Push_back ("drain");
      poleNames.Push_back ("source");

      regNames.Resize(3);
      volRegions.Resize(3);
      regIds.Resize(3);
      
      xml = myParam_->Get("poleList")->Get("Tripole"); 
    }
    else  {
      EXCEPTION("NonLinType not understood.");
    }

        
    for (UInt i=0; i< poleNames.GetSize(); i++) {
      if(! xml->Has(poleNames[i]) ) {
        EXCEPTION("poleList must define " << poleNames[i] );
      }
      xml->Get(poleNames[i])->GetValue("region", regNames[i]);

      regIds[i] = ptGrid_->GetRegion().Parse( regNames[i] );

      StdVector<RegionIdType> volRegs;
      ptGrid_->GetListOfVolumeRegions(regIds[i], volRegs);

      for (UInt reg = 0; reg < volRegs.GetSize(); reg++){
        if( regions_.Find(volRegs[reg]) < 0)
          EXCEPTION("Pole regions '" << regNames[i] << "' must refer to (surf) regions inside the ELEC CURRENT pde");
      }
    }
  }

  void ElecCurrentPDE::DefineIntegrators() {

    RegionIdType actRegion;
    BaseMaterial * actSDMat = NULL;
    
    // bool upLagrangeForm = true;
    
    //transform the type
    SubTensorType tensorType;

    if ( dim_ == 3 ) {
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

    // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    shared_ptr<FeSpace> mySpace = feFunctions_[ELEC_POTENTIAL]->GetFeSpace();

    for(UInt iRegion = 0; iRegion < regions_.GetSize() ; iRegion++){
      actRegion = regions_[iRegion];
      actSDMat    = materials_[actRegion];

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
      mySpace->SetRegionApproximation(actRegion, polyId,integId);


      PtrCoefFct coef; // if only a single coef function

      StdVector<NonLinType> matDepenTypes = regionMatDepTypes_[actRegion]; // material dependency
      StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion];   // non-linearity

      if ( matDepenTypes.Find(NLELEC_CONDUCTIVITY) != -1 && nonLinTypes.GetSize() == 0) {
        // purely temperature dependent (no further non linearity)
        shared_ptr<BaseFeFunction> myFct = feFunctions_[ELEC_POTENTIAL];
        StdVector<std::string> dispDofNames = myFct->GetResultInfo()->dofNames;
        shared_ptr<EntityList> ent = ptGrid_->GetEntityList( EntityList::ELEM_LIST, regionName );

        //get coeff-Fnc to evaluate the temperature
        ReadMaterialDependency( "elecConductivity", dispDofNames, ResultInfo::SCALAR, isComplex_,
                                ent, coef, updatedGeo_ );
        //coef-Fnc for electric conductivity
        PtrCoefFct condNL =
                  actSDMat->GetScalCoefFncNonLin( ELEC_CONDUCTIVITY_SCALAR, Global::REAL, coef);

        // create stiffness integrator
        BaseBDBInt* stiffInt = NULL;
        if( dim_ == 2 ) {
          stiffInt = new BBInt<>(new GradientOperator<FeH1,2>(), condNL,
                                 1.0, updatedGeo_ );
        } else {
          stiffInt = new BBInt<>(new GradientOperator<FeH1,3>(), condNL,
                                 1.0, updatedGeo_ );
        }
        stiffInt->SetName("StiffnessIntegrator-Temperatur-Depend");

        BiLinFormContext * stiffContext =
          new BiLinFormContext(stiffInt, STIFFNESS );
        stiffContext->SetEntities( actSDList, actSDList );
        stiffContext->SetFeFunctions( feFunctions_[ELEC_POTENTIAL], feFunctions_[ELEC_POTENTIAL] );

        feFunctions_[ELEC_POTENTIAL]->AddEntityList( actSDList );
        assemble_->AddBiLinearForm( stiffContext );
        bdbInts_.insert( std::pair<RegionIdType, BaseBDBInt*>(actRegion,stiffInt) );
      }
      else if ( nonLinTypes.Find(NLELEC_BIPOLE) != -1 || 
                nonLinTypes.Find(NLELEC_BIPOLE_TEMP_DEP) != -1 ||
                nonLinTypes.Find(NLELEC_TRIPOLE) != -1 ||
                nonLinTypes.Find(NLELEC_TRIPOLE_TEMP_DEP) != -1 
              )  {
        shared_ptr<BaseFeFunction> myFct = feFunctions_[ELEC_POTENTIAL];
        StdVector<std::string> dispDofNames = myFct->GetResultInfo()->dofNames;
        shared_ptr<EntityList> ent = ptGrid_->GetEntityList( EntityList::ELEM_LIST, regionName );

        PtrCoefFct heatCoef;

        // get elec Pot coef
        PtrCoefFct elPotCoef = this->GetCoefFct(ELEC_POTENTIAL);
        PtrCoefFct condNLNew;

        StdVector<RegionIdType> depRegionIds;
        StdVector<PtrCoefFct> dep;
        std::string intName = "";
        if (nonLinTypes.Find(NLELEC_BIPOLE) != -1 ) {
          GetPoleRegionIds(NLELEC_BIPOLE, depRegionIds);
          dep.Push_back(elPotCoef);
          condNLNew = 
            actSDMat->GetScalCoefFncMultivariateNonLin( ELEC_CONDUCTIVITY_SCALAR, NLELEC_BIPOLE, Global::REAL, dep, depRegionIds);
            /*  ELEC_CONDUCTIVITY means gamma(T) */ 
          intName = "ElecStiffnessIntegrator-Bipole-Voltage-Depend";
        }
        else if ( nonLinTypes.Find(NLELEC_BIPOLE_TEMP_DEP) != -1)  {
          // both nonlinearity and mat dependency
          //get coeff-Fnc to evaluate the temperature
          ReadMaterialDependency( "elecConductivity", dispDofNames, ResultInfo::SCALAR, isComplex_,
                                  ent, heatCoef, updatedGeo_ );
          
          GetPoleRegionIds(NLELEC_BIPOLE_TEMP_DEP, depRegionIds);

          dep.Push_back(elPotCoef);
          dep.Push_back(heatCoef);
          condNLNew = 
            actSDMat->GetScalCoefFncMultivariateNonLin( ELEC_CONDUCTIVITY_SCALAR, NLELEC_BIPOLE_TEMP_DEP,
                         Global::REAL, dep, depRegionIds);
          intName = "ElecStiffnessIntegrator-Bipole-Voltage-Temperature-Depend";
        }
        else if ( nonLinTypes.Find(NLELEC_TRIPOLE) != -1 ) {
          GetPoleRegionIds(NLELEC_TRIPOLE, depRegionIds);
          dep.Push_back(elPotCoef);
          condNLNew = 
            actSDMat->GetScalCoefFncMultivariateNonLin( ELEC_CONDUCTIVITY_SCALAR, NLELEC_TRIPOLE, Global::REAL, dep, depRegionIds);
          intName = "ElecStiffnessIntegrator-Tripole-Voltage-Depend";
        }
        else if ( nonLinTypes.Find(NLELEC_TRIPOLE_TEMP_DEP) != -1)  {
          // both nonlinearity and mat dependency
          //get coeff-Fnc to evaluate the temperature
          ReadMaterialDependency( "elecConductivity", dispDofNames, ResultInfo::SCALAR, isComplex_,
                                  ent, heatCoef, updatedGeo_ );
          
          GetPoleRegionIds(NLELEC_TRIPOLE_TEMP_DEP, depRegionIds);

          dep.Push_back(elPotCoef);
          dep.Push_back(heatCoef);
          condNLNew = 
            actSDMat->GetScalCoefFncMultivariateNonLin( ELEC_CONDUCTIVITY_SCALAR, NLELEC_TRIPOLE_TEMP_DEP,
                         Global::REAL, dep, depRegionIds);
          intName = "ElecStiffnessIntegrator-TRipole-Voltage-Temperature-Depend";
        }
        else 
          EXCEPTION("Nonlinarity not defined.");

        // create stiffness integrator
        BaseBDBInt* stiffInt = NULL;
        if( dim_ == 2 ) {
          stiffInt = new BBInt<>(new GradientOperator<FeH1,2>(), condNLNew,
                                 1.0, updatedGeo_ );
        } else {
          stiffInt = new BBInt<>(new GradientOperator<FeH1,3>(), condNLNew,
                                 1.0, updatedGeo_ );
        }
        stiffInt->SetName(intName);

        BiLinFormContext * stiffContext =
          new BiLinFormContext(stiffInt, STIFFNESS );
        stiffContext->SetEntities( actSDList, actSDList );
        stiffContext->SetFeFunctions( feFunctions_[ELEC_POTENTIAL], feFunctions_[ELEC_POTENTIAL] );

        feFunctions_[ELEC_POTENTIAL]->AddEntityList( actSDList );
        assemble_->AddBiLinearForm( stiffContext );
        bdbInts_.insert( std::pair<RegionIdType, BaseBDBInt*>(actRegion,stiffInt) );
      }
      else {
        // ----- standard real-valued stiffness integrator with constant conductivity
        BaseBDBInt * stiffInt = GetStiffIntegrator(actSDMat, tensorType, actRegion);
        stiffInt->SetName("LinElecCurrentIntegrator");
        BiLinFormContext * stiffIntDescr =
            new BiLinFormContext(stiffInt, STIFFNESS );
        feFunctions_[ELEC_POTENTIAL]->AddEntityList( actSDList );

        //stiffIntDescr->SetPtPdes(this, this);
        stiffIntDescr->SetEntities( actSDList, actSDList );
        stiffIntDescr->SetFeFunctions( feFunctions_[ELEC_POTENTIAL],feFunctions_[ELEC_POTENTIAL]);
        stiffInt->SetFeSpace( feFunctions_[ELEC_POTENTIAL]->GetFeSpace());

        assemble_->AddBiLinearForm( stiffIntDescr );
        // Important: Add bdb-integrator to global list, as we need them later
        // for calculation of postprocessing results
        bdbInts_.insert( std::pair<RegionIdType, BaseBDBInt*>(actRegion,stiffInt) );
      }
    }

  }
  
  
  BaseBDBInt * ElecCurrentPDE::GetStiffIntegrator( BaseMaterial* actSDMat,
                                            SubTensorType tensorType,
                                            RegionIdType regionId ) {

    BaseBDBInt * integ = NULL;
    bool isComplex = complexMatData_[regionId];

    shared_ptr<CoefFunction > curCoef;
    if ( isComplex ) {
      curCoef = 
        actSDMat->GetTensorCoefFnc( ELEC_CONDUCTIVITY_TENSOR,tensorType,
                                   Global::COMPLEX);
    }
    else {
      curCoef = 
        actSDMat->GetTensorCoefFnc( ELEC_CONDUCTIVITY_TENSOR,tensorType,
                                   Global::REAL);
    }
    
    if ( isComplex ) {
      if( dim_ == 2 ) {
        integ = new BDBInt<Complex,Complex >(new GradientOperator<FeH1,2,1,Complex>(),
                                             curCoef, 1.0, updatedGeo_ );
      } else {
        integ = new BDBInt<Complex,Complex >(new GradientOperator<FeH1,3,1,Complex>(),
                                             curCoef, 1.0, updatedGeo_ );
      }
    }
    else {
      if( dim_ == 2 ) {
        integ = new BDBInt<>(new GradientOperator<FeH1,2> (), 
            curCoef, 1.0, updatedGeo_ );
      } else {
        integ = new BDBInt<>(new GradientOperator<FeH1,3>(),
            curCoef, 1.0, updatedGeo_ );
      }
    }

    return integ;
  }
  
  void ElecCurrentPDE::DefineRhsLoadIntegrators() {
    shared_ptr<BaseFeFunction> myFct = feFunctions_[ELEC_POTENTIAL];
    LinearForm * lin = NULL;
    StdVector<std::string> vecDofNames;
    if(dim_ == 3) {
      vecDofNames = "x", "y", "z";
    } else {
      if(dim_ == 2 && !isaxi_)
        vecDofNames = "x", "y";
      if(dim_ == 2 && isaxi_)
        vecDofNames = "r", "z";
    }

    StdVector<shared_ptr<EntityList> > ent;
    StdVector<PtrCoefFct > coef;
    bool coefUpdateGeo;

    //=======================================================
    // NORMAL CURRENT DENSITY
    // sign is negative if brought to the RHS because surface
    // integral is negative on LHS but J = -grad(phi)*gamma
    //=======================================================
    ReadRhsExcitation( "normalCurrentDensity", vecDofNames,
                       ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo );

    std::set<RegionIdType> volRegions (regions_.Begin(), regions_.End() );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("Normal current density must be defined on elements.")
      }

      // determine dimension
      EntityIterator it = ent[i]->GetIterator();
      UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
      if( elemDim == (dim_-1) ) {
        // === SURFACE ===
        if( dim_ == 2 ) {
          if(isComplex_) {
            lin = new BUIntegrator<Complex, true>( new IdentityOperatorNormal<FeH1,2>(),
                                                   Complex(-1.0), coef[i], volRegions, coefUpdateGeo);
          } else {
            lin = new BUIntegrator<Double, true>( new IdentityOperatorNormal<FeH1,2>(),
                                                  -1.0, coef[i], volRegions, coefUpdateGeo);
          }
        } else {
          if(isComplex_) {
            lin = new BUIntegrator<Complex, true>( new IdentityOperatorNormal<FeH1,3>(),
                                                   Complex(-1.0), coef[i], volRegions, coefUpdateGeo);
          } else {
            lin = new BUIntegrator<Double, true>( new IdentityOperatorNormal<FeH1,3>(),
                                                  -1.0, coef[i], volRegions, coefUpdateGeo);
          }
        }
        lin->SetName("NormalCurrentDensityInt");
      } else {
        // === VOLUME ===
        EXCEPTION("Specifying current density in a volume is not supported.");
      }

      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
    }
  }


  void ElecCurrentPDE::DefineSolveStep() {
    solveStep_ = new StdSolveStep(*this);
  }

  void ElecCurrentPDE::InitTimeStepping() {

    // Until now no effective mass formulation in the trapezoidal
    //  integration scheme is implemented!
    //TS_alg_ = new Trapezoidal( algsys_, olasNode_ );
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(GLMScheme::TRAPEZOIDAL, 0) );

    feFunctions_[ELEC_POTENTIAL]->SetTimeScheme(myScheme);
  }
  
  void ElecCurrentPDE::DefinePrimaryResults() {
    
    shared_ptr<BaseFeFunction> feFct = feFunctions_[ELEC_POTENTIAL];
    
    // Electric Potential
    shared_ptr<ResultInfo> res1( new ResultInfo);
    res1->resultType = ELEC_POTENTIAL;
    
    res1->dofNames = "";
    res1->unit = MapSolTypeToUnit(ELEC_POTENTIAL);
    res1->definedOn = ResultInfo::MapSolTypeToDefinedOn(ELEC_POTENTIAL);
    res1->entryType = ResultInfo::SCALAR;
    feFunctions_[ELEC_POTENTIAL]->SetResultInfo(res1);

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    hdbcSolNameMap_[ELEC_POTENTIAL] = "ground";
    idbcSolNameMap_[ELEC_POTENTIAL] = "potential";
    
    res1->SetFeFunction(feFunctions_[ELEC_POTENTIAL]);
    results_.Push_back( res1 );
    DefineFieldResult( feFunctions_[ELEC_POTENTIAL], res1 );
    
  }
  
  void ElecCurrentPDE::DefineNcIntegrators() {
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
        if (dim_ == 2)
          DefineNitscheCoupling<2,1>(ELEC_POTENTIAL, *ncIt );
        else
          DefineNitscheCoupling<3,1>(ELEC_POTENTIAL, *ncIt );
        break;
      default:
        EXCEPTION("Unknown type of ncInterface");
        break;
      }
    }
  }
  
  void ElecCurrentPDE::DefinePostProcResults() {

    shared_ptr<BaseFeFunction> feFct = feFunctions_[ELEC_POTENTIAL];
    
    // === ELECTRIC FIELD INTENSITY ===
    shared_ptr<ResultInfo> ef ( new ResultInfo );
    ef->resultType = ELEC_FIELD_INTENSITY;
    ef->SetVectorDOFs(dim_, isaxi_);
    ef->unit = MapSolTypeToUnit(ELEC_FIELD_INTENSITY);
    ef->definedOn = ResultInfo::MapSolTypeToDefinedOn(ELEC_FIELD_INTENSITY);
    ef->entryType = ResultInfo::VECTOR;
    shared_ptr<CoefFunctionFormBased> eFunc;
    if( isComplex_ ) {
      eFunc.reset(new CoefFunctionBOp<Complex>(feFct, ef, -1.0));
    } else {
      eFunc.reset(new CoefFunctionBOp<Double>(feFct, ef, -1.0));
    }
    DefineFieldResult( eFunc, ef );
    stiffFormCoefs_.insert(eFunc);
    
    // === ELECTRIC CURRENT DENSITY ===
    shared_ptr<ResultInfo> flux ( new ResultInfo );
    flux->resultType = ELEC_CURRENT_DENSITY;
    flux->SetVectorDOFs(dim_, isaxi_);
    flux->unit = MapSolTypeToUnit(ELEC_CURRENT_DENSITY);
    flux->definedOn = ResultInfo::MapSolTypeToDefinedOn(ELEC_CURRENT_DENSITY);
    flux->entryType = ResultInfo::VECTOR;
    shared_ptr<CoefFunctionFormBased> fluxFunc;
    if( isComplex_ ) {
      fluxFunc.reset(new CoefFunctionFlux<Complex>(feFct, flux, Complex(-1.0)));
    } else {
      fluxFunc.reset(new CoefFunctionFlux<Double>(feFct, flux, -1.0));
    }
    DefineFieldResult( fluxFunc, flux );
    stiffFormCoefs_.insert(fluxFunc);


    // \int_\Omega_C grad(V) \cdot grad(V)
    shared_ptr<ResultInfo> elecGradVInt( new ResultInfo );
    elecGradVInt->resultType = ELEC_GRAD_V_INT;
    elecGradVInt->definedOn = ResultInfo::MapSolTypeToDefinedOn(ELEC_GRAD_V_INT);
    elecGradVInt->entryType = ResultInfo::SCALAR;
    elecGradVInt->dofNames = "";
    elecGradVInt->unit = MapSolTypeToUnit(ELEC_GRAD_V_INT);
    shared_ptr<CoefFunctionFormBased> sigmaElecGradVFunc;
    if( isComplex_ ) {
      sigmaElecGradVFunc.reset(new CoefFunctionBdBKernel<Complex>(feFct, 1.0));
    } else {
      sigmaElecGradVFunc.reset(new CoefFunctionBdBKernel<Double>(feFct, 1.0));
    }
    stiffFormCoefs_.insert(sigmaElecGradVFunc);

    shared_ptr<ResultFunctor> ElecGradVIntFunc;
    if( isComplex_ ) {
      ElecGradVIntFunc.reset(new ResultFunctorIntegrate<Complex>(sigmaElecGradVFunc, feFct, elecGradVInt));
    } else {
      ElecGradVIntFunc.reset(new ResultFunctorIntegrate<Double>(sigmaElecGradVFunc, feFct, elecGradVInt));
    }
    resultFunctors_[ELEC_GRAD_V_INT] = ElecGradVIntFunc;
    availResults_.insert( elecGradVInt );


    // == ELECTRIC_NORMAL_CURRENT_DENSITY ==
    shared_ptr<ResultInfo> fluxNormal ( new ResultInfo );
    fluxNormal->resultType = ELEC_NORMAL_CURRENT_DENSITY;
    fluxNormal->dofNames = "";
    fluxNormal->unit = MapSolTypeToUnit(ELEC_NORMAL_CURRENT_DENSITY);
    fluxNormal->entryType = ResultInfo::SCALAR;
    fluxNormal->definedOn = ResultInfo::MapSolTypeToDefinedOn(ELEC_NORMAL_CURRENT_DENSITY);

    shared_ptr<CoefFunctionSurf> fluxFctNormal;
    fluxFctNormal.reset(new CoefFunctionSurf(true, 1.0, fluxNormal));
    DefineFieldResult( fluxFctNormal, fluxNormal );
    surfCoefFcts_[fluxFctNormal] = fluxFunc;

    // === ELEC_CURRENT ===
    shared_ptr<ResultInfo> current (new ResultInfo );
    current.reset(new ResultInfo);
    current->resultType = ELEC_CURRENT;
    current->dofNames = "";
    current->unit = MapSolTypeToUnit(ELEC_CURRENT);
    current->entryType = ResultInfo::SCALAR;
    current->definedOn = ResultInfo::MapSolTypeToDefinedOn(ELEC_CURRENT);
    // Current = \int \vec j \dot \vec n dA
    shared_ptr<ResultFunctor> currentFct;
    if (isComplex_) {
      currentFct.reset( new ResultFunctorIntegrate<Complex>(fluxFctNormal, 
                                                           feFct, current) );
    } else {
      currentFct.reset( new ResultFunctorIntegrate<Double>(fluxFctNormal, 
                                                           feFct, current) );
    }
    resultFunctors_[ELEC_CURRENT] = currentFct;
    availResults_.insert(current);


    // === ELECTRIC POWER DENSITY ===
    shared_ptr<ResultInfo> ed ( new ResultInfo );
    ed->resultType = ELEC_POWER_DENSITY;
    ed->dofNames = "";
    ed->unit = MapSolTypeToUnit(ELEC_POWER_DENSITY);
    ed->definedOn = ResultInfo::MapSolTypeToDefinedOn(ELEC_POWER_DENSITY);
    ed->entryType = ResultInfo::SCALAR;
    availResults_.insert( ed );
    shared_ptr<CoefFunctionFormBased> edFunc;
    if( isComplex_ ) {
      /* for harmonic analysis, the effective value is 1/sqrt(2) \hat I \times 1/sqrt(2) \hat U */
      edFunc.reset(new CoefFunctionBdBKernel<Complex>(feFct, 0.5));
    } else {
      edFunc.reset(new CoefFunctionBdBKernel<Double>(feFct, 1));
    }
    DefineFieldResult( edFunc, ed );
    stiffFormCoefs_.insert(edFunc);
    
    // Electric power
    shared_ptr<ResultInfo> energy( new ResultInfo );
    energy->resultType = ELEC_POWER;
    energy->definedOn = ResultInfo::MapSolTypeToDefinedOn(ELEC_POWER);
    energy->entryType = ResultInfo::SCALAR;
    energy->dofNames = "";
    energy->unit = MapSolTypeToUnit(ELEC_POWER);
    availResults_.insert( energy );
    shared_ptr<ResultFunctor> energyFunc;
    if( isComplex_ ) {
      /* for harmonic analysis, the effective value is 1/sqrt(2) \hat I \times 1/sqrt(2) \hat U */
      energyFunc.reset(new EnergyResultFunctor<Complex>(feFct, energy,0.5));
    } else {
      energyFunc.reset(new EnergyResultFunctor<Double>(feFct, energy, 1));
    }
    resultFunctors_[ELEC_POWER] = energyFunc;
    stiffFormFunctors_.insert(energyFunc);
    
  }

  
  std::map<SolutionType, shared_ptr<FeSpace> >
  ElecCurrentPDE::CreateFeSpaces(const std::string& formulation, PtrParamNode infoNode) {
    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    if(formulation == "default" || formulation == "H1"){
      PtrParamNode potSpaceNode = infoNode->Get("elecPotential");
      crSpaces[ELEC_POTENTIAL] =
        FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1, ptGrid_);
      crSpaces[ELEC_POTENTIAL]->Init(solStrat_);
    }else{
      EXCEPTION("The formulation " << formulation << "of electricCurrent PDE is not known!");
    }
    return crSpaces;
  }
}
