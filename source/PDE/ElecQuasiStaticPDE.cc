#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <set>

#include "ElecQuasiStaticPDE.hh"

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

  DEFINE_LOG(elecquasistaticpde, "pde.elecquasistatic")

  // ***************
  //   Constructor
  // ***************
  ElecQuasistaticPDE::ElecQuasistaticPDE( Grid* aptgrid, PtrParamNode paramNode,
                    PtrParamNode infoNode,
                    shared_ptr<SimState> simState, Domain* domain )
    :SinglePDE( aptgrid, paramNode, infoNode, simState, domain ) {


    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_          = "elecQuasistatic";
    pdematerialclass_ = ELECQUASISTATIC;
 
    nonLin_         = false;
    nonLinMaterial_ = false;
    
    // Set PDE subtype
    subType_ = "3d";
    if ( ptGrid_->GetDim() == 2 ) {
      if ( ptGrid_->IsAxi() ) {
        subType_ = "axi";
        isaxi_ = true;
      } 
      else {
        subType_ = "plane";
        isaxi_ = false;
      }
    } 

    //! Always use updated Lagrangian formulation 
    updatedGeo_     = true;

    //! currently we just support harmonic analysis
    isMaterialComplex_ = true;
  }

  
  void ElecQuasistaticPDE::InitNonLin() {

    SinglePDE::InitNonLin();

  }


  void ElecQuasistaticPDE::DefineIntegrators() {

    RegionIdType actRegion;
    BaseMaterial * actSDMat = NULL;

    //check analysis type:
    if ( analysistype_ != HARMONIC )
      EXCEPTION("ElecQuasistaticPDE: Currently we just support HARMONIC analysis!");

//  //set tensor type
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
      
      // --- Set the approximation for the current region ---
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();
      mySpace->SetRegionApproximation(actRegion, polyId,integId);
       
      //       
      // stiffness integrator: int (electric conductivity * grad N_i grad N_j) dx
      //
      shared_ptr<CoefFunction > condCoef;
      const StdVector<NonLinType>& matDepenTypes = regionMatDepTypes_[actRegion]; // material dependency
      if ( matDepenTypes.Find(NLELEC_CONDUCTIVITY) != -1 ) {
        //we read the computed conductivity tensor for each element from file
        //shared_ptr<BaseFeFunction> myFct = feFunctions_[ELEC_POTENTIAL];
        //StdVector<std::string> dispDofNames = myFct->GetResultInfo()->dofNames;
        shared_ptr<ResultInfo> resultInfo = GetResultInfo(ELEC_COND_TENSOR); 
        shared_ptr<EntityList> ent = ptGrid_->GetEntityList( EntityList::ELEM_LIST, regionName );

        //get coeff-Fnc for the electric conductivity
        ReadMaterialDependency( "elecConductivity", resultInfo->dofNames, resultInfo->entryType, false,
                                ent, condCoef, updatedGeo_ );
      }
      else {
        condCoef = actSDMat->GetTensorCoefFnc(ELEC_CONDUCTIVITY_TENSOR, tensorType, Global::REAL);
      }

      BaseBDBInt *stiffInt;
      if( dim_ == 2 ) {
        stiffInt = new BDBInt<>(new GradientOperator<FeH1,2>(), condCoef, 1.0, updatedGeo_ );
      } else {
        stiffInt = new BDBInt<>(new GradientOperator<FeH1,3>(), condCoef, 1.0, updatedGeo_ );
      }
      stiffInt->SetName("LinElecQuasistaticStiffIntegrator");
      BiLinFormContext * stiffIntDescr = new BiLinFormContext(stiffInt, STIFFNESS);
      feFunctions_[ELEC_POTENTIAL]->AddEntityList( actSDList );

      stiffIntDescr->SetEntities( actSDList, actSDList );
      stiffIntDescr->SetFeFunctions(feFunctions_[ELEC_POTENTIAL],feFunctions_[ELEC_POTENTIAL]);
      stiffInt->SetFeSpace(feFunctions_[ELEC_POTENTIAL]->GetFeSpace());

      assemble_->AddBiLinearForm( stiffIntDescr );
      bdbInts_.insert( std::pair<RegionIdType, BaseBDBInt*>(actRegion,stiffInt) ); //for postprocessing results

      //
      // ----- mass integrator: int (permability * grad N_i grad N_j) dx
      //
      shared_ptr<CoefFunction > permCoef;
      permCoef = actSDMat->GetTensorCoefFnc(ELEC_PERMITTIVITY_TENSOR, tensorType, Global::REAL);

      BaseBDBInt *massInt;
      if( dim_ == 2 ) {
        massInt = new BDBInt<Complex,Complex>(new GradientOperator<FeH1,2>(), permCoef, 1.0, updatedGeo_ );
      } else {
        massInt = new BDBInt<Complex,Complex>(new GradientOperator<FeH1,3>(), permCoef, 1.0, updatedGeo_ );
      }

      massInt->SetName("LinElecQuasistaticMassIntegrator");
      massInt->SetFeSpace(feFunctions_[ELEC_POTENTIAL]->GetFeSpace());

      BiLinFormContext *massIntDescr = new BiLinFormContext(massInt, DAMPING);
      massIntDescr->SetEntities( actSDList, actSDList );
      massIntDescr->SetFeFunctions(feFunctions_[ELEC_POTENTIAL],feFunctions_[ELEC_POTENTIAL]);
      assemble_->AddBiLinearForm( massIntDescr );

      //for calculation of postprocessing results
      //bdbInts_.insert( std::pair<RegionIdType, BaseBDBInt*>(actRegion,massInt) );
      massInts_.insert( std::pair<RegionIdType, BaseBDBInt*>(actRegion,massInt) ); //for postprocessing results 
    }

  }
  
  
  void ElecQuasistaticPDE::DefineRhsLoadIntegrators() {
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
                       ResultInfo::SCALAR, isComplex_, ent, coef, coefUpdateGeo );

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
            lin = new BUIntegrator<Complex, true>( new IdentityOperator<FeH1,2>(),
                                                   Complex(-1.0), coef[i], volRegions, coefUpdateGeo);
          } else {
            lin = new BUIntegrator<Double, true>( new IdentityOperator<FeH1,2>(),
                                                  -1.0, coef[i], volRegions, coefUpdateGeo);
          }
        } else {
          if(isComplex_) {
            lin = new BUIntegrator<Complex, true>( new IdentityOperator<FeH1,3>(),
                                                   Complex(-1.0), coef[i], volRegions, coefUpdateGeo);
          } else {
            lin = new BUIntegrator<Double, true>( new IdentityOperator<FeH1,3>(),
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


  void ElecQuasistaticPDE::DefineSolveStep() {
    solveStep_ = new StdSolveStep(*this);
  }

  void ElecQuasistaticPDE::InitTimeStepping()
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

    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(scheme, 1.0));
    feFunctions_[ELEC_POTENTIAL]->SetTimeScheme(myScheme);
  }
  
  void ElecQuasistaticPDE::DefinePrimaryResults() {
    
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
    
    // === define conductivity tensor ================
    // Check for subType
    StdVector<std::string> condDofNames;
    
    if(subType_ == "3d")
      condDofNames = "xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz";
    else if(subType_ == "plane")
      condDofNames = "xx", "xy", "yx" "yy";
    else if(subType_ == "axi")
      EXCEPTION("Conductivity Tensor for axisymmetric case not implemented");
    
    shared_ptr<ResultInfo> condTensor(new ResultInfo);
    condTensor->resultType = ELEC_COND_TENSOR;
    condTensor->dofNames = condDofNames;
    condTensor->unit = MapSolTypeToUnit(ELEC_COND_TENSOR);
    condTensor->entryType = ResultInfo::TENSOR;
    condTensor->SetFeFunction(feFunctions_[ELEC_POTENTIAL]);
    condTensor->definedOn = ResultInfo::MapSolTypeToDefinedOn(ELEC_COND_TENSOR);
    availResults_.insert(condTensor);
  }
  
  void ElecQuasistaticPDE::DefineNcIntegrators() {
    StdVector< NcInterfaceInfo >::iterator ncIt = ncInterfaces_.Begin(),
                                           endIt = ncInterfaces_.End();
    for ( ; ncIt != endIt; ++ncIt ) {
      switch (ncIt->type) {
      case NC_MORTAR:
        EXCEPTION("ncInterface of Mortar type is not implemented for ElecQuasistaticPDE")
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
  
  void ElecQuasistaticPDE::DefinePostProcResults() {

    shared_ptr<BaseFeFunction> feFct = feFunctions_[ELEC_POTENTIAL];
    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
    StdVector<std::string> vecComponents;
    if(dim_ == 3) {
      vecComponents = "x", "y", "z";
    } else {
      if(dim_ == 2 && !isaxi_)
        vecComponents = "x", "y";
      if(dim_ == 2 && isaxi_)
        vecComponents = "r", "z";
    }
 
    // === ELECTRIC POTENTIAL - 1ST DERIVATIVE ===
    shared_ptr<ResultInfo> potDot(new ResultInfo);
    potDot->resultType = ELEC_POTENTIAL_DERIV_1;
    potDot->dofNames = "";
    potDot->unit = MapSolTypeToUnit(ELEC_POTENTIAL_DERIV_1);
    potDot->definedOn = ResultInfo::MapSolTypeToDefinedOn(ELEC_POTENTIAL_DERIV_1);
    potDot->entryType = ResultInfo::SCALAR;
    availResults_.insert( potDot );    
    DefineTimeDerivResult( ELEC_POTENTIAL_DERIV_1, 1, ELEC_POTENTIAL );      

    // === ELECTRIC FIELD INTENSITY ===
    shared_ptr<ResultInfo> ef ( new ResultInfo );
    ef->resultType = ELEC_FIELD_INTENSITY;
    ef->SetVectorDOFs(dim_, isaxi_);
    ef->dofNames = vecComponents;
    ef->unit = MapSolTypeToUnit(ELEC_FIELD_INTENSITY);
    ef->definedOn = ResultInfo::MapSolTypeToDefinedOn(ELEC_FIELD_INTENSITY);
    ef->entryType = ResultInfo::VECTOR;
    availResults_.insert(ef);
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
    availResults_.insert(flux);
    shared_ptr<CoefFunctionFormBased> fluxFunc;
    if( isComplex_ ) {
      fluxFunc.reset(new CoefFunctionFlux<Complex>(feFct, flux, Complex(-1.0)));
    } else {
      fluxFunc.reset(new CoefFunctionFlux<Double>(feFct, flux, -1.0));
    }
    DefineFieldResult( fluxFunc, flux );
    stiffFormCoefs_.insert(fluxFunc);

    // == ELECTRIC_NORMAL_CURRENT_DENSITY ==
    shared_ptr<ResultInfo> fluxNormal ( new ResultInfo );
    fluxNormal->resultType = ELEC_NORMAL_CURRENT_DENSITY;
    fluxNormal->dofNames = "";
    fluxNormal->unit = MapSolTypeToUnit(ELEC_NORMAL_CURRENT_DENSITY);
    fluxNormal->entryType = ResultInfo::SCALAR;
    fluxNormal->definedOn = ResultInfo::MapSolTypeToDefinedOn(ELEC_NORMAL_CURRENT_DENSITY);
    availResults_.insert(fluxNormal);
    shared_ptr<CoefFunctionSurf> fluxFctNormal;
    fluxFctNormal.reset(new CoefFunctionSurf(true, -1.0, fluxNormal));
    DefineFieldResult( fluxFctNormal, fluxNormal );
    surfCoefFcts_[fluxFctNormal] = fluxFunc;

    // === ELEC_CURRENT ===
    shared_ptr<ResultInfo> elecCurrent (new ResultInfo );
    elecCurrent.reset(new ResultInfo);
    elecCurrent->resultType = ELEC_CURRENT;
    elecCurrent->dofNames = "";
    elecCurrent->unit = MapSolTypeToUnit(ELEC_CURRENT);
    elecCurrent->entryType = ResultInfo::SCALAR;
    elecCurrent->definedOn = ResultInfo::MapSolTypeToDefinedOn(ELEC_CURRENT);
    // Current = \int \vec j \dot \vec n dA
    shared_ptr<ResultFunctor> elecCurrentFct;
    if (isComplex_) {
      elecCurrentFct.reset( new ResultFunctorIntegrate<Complex>(fluxFctNormal, 
                                                           feFct, elecCurrent) );
    } else {
      elecCurrentFct.reset( new ResultFunctorIntegrate<Double>(fluxFctNormal, 
                                                           feFct, elecCurrent) );
    }
    resultFunctors_[ELEC_CURRENT] = elecCurrentFct;
    availResults_.insert(elecCurrent);


    // === ELECTRIC DISPLACEMENT CURRENT DENSITY ===
    shared_ptr<ResultInfo> dflux ( new ResultInfo );
    dflux->resultType = DISPLACEMENT_CURRENT_FIELD_INTENSITY;
    dflux->SetVectorDOFs(dim_, isaxi_);
    dflux->unit = MapSolTypeToUnit(DISPLACEMENT_CURRENT_FIELD_INTENSITY);
    dflux->definedOn = ResultInfo::MapSolTypeToDefinedOn(DISPLACEMENT_CURRENT_FIELD_INTENSITY);
    dflux->entryType = ResultInfo::VECTOR;
    availResults_.insert(dflux);
    shared_ptr<BaseFeFunction> potDotFct =
          timeDerivFeFunctions_[ELEC_POTENTIAL_DERIV_1];
    shared_ptr<CoefFunctionFormBased> dfluxFunc;
    if( isComplex_ ) {
      dfluxFunc.reset(new CoefFunctionFlux<Complex>(potDotFct, dflux, Complex(-1.0)));
    } else {
      dfluxFunc.reset(new CoefFunctionFlux<Double>(potDotFct, dflux, -1.0));
    }
    DefineFieldResult( dfluxFunc, dflux );
    massFormCoefs_.insert(dfluxFunc);

    // == DISPLACEMENT_NORMAL_CURRENT_DENSITY ==
    shared_ptr<ResultInfo> dfluxNormal ( new ResultInfo );
    dfluxNormal->resultType = DISPLACEMENT_NORMAL_CURRENT_DENSITY;
    dfluxNormal->dofNames = "";
    dfluxNormal->unit = MapSolTypeToUnit(DISPLACEMENT_NORMAL_CURRENT_DENSITY);
    dfluxNormal->entryType = ResultInfo::SCALAR;
    dfluxNormal->definedOn = ResultInfo::MapSolTypeToDefinedOn(DISPLACEMENT_NORMAL_CURRENT_DENSITY);
    availResults_.insert(dfluxNormal);
    shared_ptr<CoefFunctionSurf> dfluxFctNormal;
    dfluxFctNormal.reset(new CoefFunctionSurf(true, -1.0, dfluxNormal));
    DefineFieldResult( dfluxFctNormal, dfluxNormal );
    surfCoefFcts_[dfluxFctNormal] = dfluxFunc;

    // === DISPLACEMENT_CURRENT ===
    shared_ptr<ResultInfo> dispCurrent (new ResultInfo );
    dispCurrent.reset(new ResultInfo);
    dispCurrent->resultType = DISPLACEMENT_CURRENT_SURF;
    dispCurrent->dofNames = "";
    dispCurrent->unit = MapSolTypeToUnit(DISPLACEMENT_CURRENT_SURF);
    dispCurrent->entryType = ResultInfo::SCALAR;
    dispCurrent->definedOn = ResultInfo::MapSolTypeToDefinedOn(DISPLACEMENT_CURRENT_SURF);
    // Current = \int \vec j \dot \vec n dA
    shared_ptr<ResultFunctor> dispCurrentFct;
    if (isComplex_) {
      dispCurrentFct.reset( new ResultFunctorIntegrate<Complex>(dfluxFctNormal, 
                                                           feFct, dispCurrent) );
    } else {
      dispCurrentFct.reset( new ResultFunctorIntegrate<Double>(dfluxFctNormal, 
                                                           feFct, dispCurrent) );
    }
    resultFunctors_[DISPLACEMENT_CURRENT_SURF] = dispCurrentFct;
    availResults_.insert(dispCurrent);


 // === ELECTRIC CURRENT PLUS ELECTRIC DISPLACEMENT CURRENT DENSITY ===
    shared_ptr<ResultInfo> edflux ( new ResultInfo );
    edflux->resultType = ELECTRIC_AND_DISPLACEMENT_CURRENT_DENSITY;
    edflux->SetVectorDOFs(dim_, isaxi_);
    edflux->unit = MapSolTypeToUnit(ELECTRIC_AND_DISPLACEMENT_CURRENT_DENSITY);
    edflux->definedOn = ResultInfo::MapSolTypeToDefinedOn(ELECTRIC_AND_DISPLACEMENT_CURRENT_DENSITY);
    edflux->entryType = ResultInfo::VECTOR;
    availResults_.insert( edflux );
    PtrCoefFct totalCurrDensity;
    totalCurrDensity =
                CoefFunction::Generate( mp_, part,
                CoefXprBinOp(mp_,fluxFunc,dfluxFunc,CoefXpr::OP_ADD));
    DefineFieldResult( totalCurrDensity, edflux );

 // === ELECTRIC CURRENT PLUS ELECTRIC DISPLACEMENT NORMAL CURRENT DENSITY ===
    shared_ptr<ResultInfo> edfluxNormal ( new ResultInfo );
    edfluxNormal->resultType = ELECTRIC_AND_DISPLACEMENT_NORMAL_CURRENT_DENSITY;
    edfluxNormal->dofNames = "";
    edfluxNormal->unit = MapSolTypeToUnit(ELECTRIC_AND_DISPLACEMENT_NORMAL_CURRENT_DENSITY);
    edfluxNormal->entryType = ResultInfo::SCALAR;
    edfluxNormal->definedOn = ResultInfo::MapSolTypeToDefinedOn(ELECTRIC_AND_DISPLACEMENT_NORMAL_CURRENT_DENSITY);
    availResults_.insert( edfluxNormal );
    PtrCoefFct totalCurrNormalDensity;
    totalCurrNormalDensity =
                CoefFunction::Generate( mp_, part,
                CoefXprBinOp(mp_,fluxFctNormal,dfluxFctNormal,CoefXpr::OP_ADD));
    DefineFieldResult( totalCurrNormalDensity, edfluxNormal );

    // === TOTAL CURRENT ===
    shared_ptr<ResultInfo> current (new ResultInfo );
    current.reset(new ResultInfo);
    current->resultType = ELEC_AND_DISPLACEMENT_CURRENT;
    current->dofNames = "";
    current->unit = MapSolTypeToUnit(ELEC_AND_DISPLACEMENT_CURRENT);
    current->entryType = ResultInfo::SCALAR;
    current->definedOn = ResultInfo::MapSolTypeToDefinedOn(ELEC_AND_DISPLACEMENT_CURRENT);
    availResults_.insert(current);
    // Current = \int \vec j \dot \vec n dA
    shared_ptr<ResultFunctor> currentFct;
    if (isComplex_) {
      currentFct.reset( new ResultFunctorIntegrate<Complex>(totalCurrNormalDensity, 
                                                           feFct, current) );
    } else {
      currentFct.reset( new ResultFunctorIntegrate<Double>(totalCurrNormalDensity, 
                                                           feFct, current) );
    }
    resultFunctors_[ELEC_AND_DISPLACEMENT_CURRENT] = currentFct;
   
  }


  std::map<SolutionType, shared_ptr<FeSpace> >
  ElecQuasistaticPDE::CreateFeSpaces(const std::string& formulation, PtrParamNode infoNode) {
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
