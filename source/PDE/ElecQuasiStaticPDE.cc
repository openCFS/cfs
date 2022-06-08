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

//    //transform the type
//    SubTensorType tensorType;
//
//    if ( dim_ == 3 ) {
//      tensorType = FULL;
//    }
//    else {
//      if ( isaxi_ == true ) {
//        tensorType = AXI;
//      }
//      else {
//        // 2d: plane case
//        tensorType = PLANE_STRAIN;
//      }
//    }

    // Define integrators for "standard" materials
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
      
      // --- Set the approximation for the current region ---
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();
      mySpace->SetRegionApproximation(actRegion, polyId,integId);

      // ----- standard complex valued stiffness integrator with commplex material parameter
      // ----- alpha = electricConductivity + i*omega*permittivity
      PtrCoefFct omega = CoefFunction::Generate( mp_,  Global::COMPLEX, "0", "2*pi*f");
      shared_ptr<CoefFunction > condCoef;
      condCoef = actSDMat->GetScalCoefFnc(ELEC_CONDUCTIVITY_SCALAR, Global::REAL);
      shared_ptr<CoefFunction > permCoef;
      permCoef = actSDMat->GetScalCoefFnc(ELEC_PERMITTIVITY_SCALAR, Global::REAL);
      shared_ptr<CoefFunction > matCoef;
      shared_ptr<CoefFunction > omegaPermCoef;
      omegaPermCoef = CoefFunction::Generate( mp_, Global::COMPLEX,
                                    CoefXprBinOp(mp_, permCoef, omega, CoefXpr::OP_MULT));
      matCoef = CoefFunction::Generate( mp_, Global::COMPLEX,
                              CoefXprBinOp(mp_, condCoef, omegaPermCoef, CoefXpr::OP_ADD));

      BaseBDBInt *stiffInt;
      if( dim_ == 2 ) {
        stiffInt = new BBInt<Complex,Complex >(new GradientOperator<FeH1,2,1,Complex>(),
                                                matCoef, 1.0, updatedGeo_ );
      } else {
        stiffInt = new BBInt<Complex,Complex >(new GradientOperator<FeH1,3,1,Complex>(),
                                                matCoef, 1.0, updatedGeo_ );
      }
      stiffInt->SetName("LinElecQuasistaticStiffIntegrator");
      stiffInt->SetFeSpace(feFunctions_[ELEC_POTENTIAL]->GetFeSpace());
      BiLinFormContext * stiffIntDescr = new BiLinFormContext(stiffInt, STIFFNESS);
      feFunctions_[ELEC_POTENTIAL]->AddEntityList( actSDList );

      stiffIntDescr->SetEntities( actSDList, actSDList );
      stiffIntDescr->SetFeFunctions(feFunctions_[ELEC_POTENTIAL],feFunctions_[ELEC_POTENTIAL]);

      assemble_->AddBiLinearForm( stiffIntDescr );
      bdbInts_.insert( std::pair<RegionIdType, BaseBDBInt*>(actRegion,stiffInt) );  //for calculation of postprocessing results


//      shared_ptr<CoefFunction > stiffCoef;
//      stiffCoef = actSDMat->GetTensorCoefFnc(ELEC_CONDUCTIVITY_TENSOR,tensorType, Global::REAL);
//      BaseBDBInt *stiffInt;
//      if( dim_ == 2 ) {
//        stiffInt = new BDBInt<>(new GradientOperator<FeH1,2>(), stiffCoef, 1.0, updatedGeo_ );
//      } else {
//        stiffInt = new BDBInt<>(new GradientOperator<FeH1,3>(), stiffCoef, 1.0, updatedGeo_ );
//      }
//      stiffInt->SetName("LinElecQuasistaticStiffIntegrator");
//      BiLinFormContext * stiffIntDescr = new BiLinFormContext(stiffInt, STIFFNESS);
//      feFunctions_[ELEC_POTENTIAL]->AddEntityList( actSDList );
//
//      stiffIntDescr->SetEntities( actSDList, actSDList );
//      stiffIntDescr->SetFeFunctions(feFunctions_[ELEC_POTENTIAL],feFunctions_[ELEC_POTENTIAL]);
//      stiffInt->SetFeSpace(feFunctions_[ELEC_POTENTIAL]->GetFeSpace());
//
//      assemble_->AddBiLinearForm( stiffIntDescr );
//      bdbInts_[actRegion] = stiffInt; // for calculation of postprocessing results
//
//      // ----- standard real-valued mass integrator with constant permittivity
//      shared_ptr<CoefFunction > massCoef;
//      massCoef = actSDMat->GetTensorCoefFnc(ELEC_PERMITTIVITY_TENSOR,tensorType, Global::REAL);
//      BaseBDBInt *massInt;
//      if( dim_ == 2 ) {
//        massInt = new BDBInt<>(new GradientOperator<FeH1,2>(), massCoef, 1.0, updatedGeo_ );
//      } else {
//        massInt = new BDBInt<>(new GradientOperator<FeH1,3>(), massCoef, 1.0, updatedGeo_ );
//      }
//
//      massInt->SetName("LinElecQuasistaticMassIntegrator");
//      massInt->SetFeSpace(feFunctions_[ELEC_POTENTIAL]->GetFeSpace());
//
//      BiLinFormContext *massIntDescr = new BiLinFormContext(massInt, DAMPING);
//      massIntDescr->SetEntities( actSDList, actSDList );
//      massIntDescr->SetFeFunctions(feFunctions_[ELEC_POTENTIAL],feFunctions_[ELEC_POTENTIAL]);
//
//      assemble_->AddBiLinearForm( massIntDescr );
//      massInts_[actRegion] = massInt; // for calculation of postprocessing results
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

  void ElecQuasistaticPDE::InitTimeStepping() {

    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(GLMScheme::TRAPEZOIDAL, 1.0) );

    feFunctions_[ELEC_POTENTIAL]->SetTimeScheme(myScheme);
  }
  
  void ElecQuasistaticPDE::DefinePrimaryResults() {
    
    shared_ptr<BaseFeFunction> feFct = feFunctions_[ELEC_POTENTIAL];
    
    // Electric Potential
    shared_ptr<ResultInfo> res1( new ResultInfo);
    res1->resultType = ELEC_POTENTIAL;
    
    res1->dofNames = "";
    res1->unit = "V";
    res1->definedOn = ResultInfo::NODE;
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
    
    // === ELECTRIC FIELD INTENSITY ===
    shared_ptr<ResultInfo> ef ( new ResultInfo );
    ef->resultType = ELEC_FIELD_INTENSITY;
    ef->SetVectorDOFs(dim_, isaxi_);
    ef->unit = "V/m";
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
    
    // === ELECTRIC CURRENT DENSITY ===
    shared_ptr<ResultInfo> flux ( new ResultInfo );
    flux->resultType = ELEC_CURRENT_DENSITY;
    flux->SetVectorDOFs(dim_, isaxi_);
    flux->unit = "A/m^2";
    flux->definedOn = ResultInfo::ELEMENT;
    flux->entryType = ResultInfo::VECTOR;
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
    fluxNormal->unit = "A/m^2";
    fluxNormal->entryType = ResultInfo::SCALAR;
    fluxNormal->definedOn = ResultInfo::SURF_ELEM;

    shared_ptr<CoefFunctionSurf> fluxFctNormal;
    fluxFctNormal.reset(new CoefFunctionSurf(true, -1.0, fluxNormal));
    DefineFieldResult( fluxFctNormal, fluxNormal );
    surfCoefFcts_[fluxFctNormal] = fluxFunc;

    // === ELEC_CURRENT ===
    shared_ptr<ResultInfo> current (new ResultInfo );
    current.reset(new ResultInfo);
    current->resultType = ELEC_CURRENT;
    current->dofNames = "";
    current->unit = "A";
    current->entryType = ResultInfo::SCALAR;
    current->definedOn = ResultInfo::SURF_REGION; /* cannot output to gmsh */
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


//    // === ELECTRIC POWER DENSITY ===
//    shared_ptr<ResultInfo> ed ( new ResultInfo );
//    ed->resultType = ELEC_POWER_DENSITY;
//    ed->dofNames = "";
//    ed->unit = "W/m^3";
//    ed->definedOn = ResultInfo::ELEMENT;
//    ed->entryType = ResultInfo::SCALAR;
//    availResults_.insert( ed );
//    shared_ptr<CoefFunctionFormBased> edFunc;
//    if( isComplex_ ) {
//      /* for harmonic analysis, the effective value is 1/sqrt(2) \hat I \times 1/sqrt(2) \hat U */
//      edFunc.reset(new CoefFunctionBdBKernel<Complex>(feFct, 0.5));
//    } else {
//      edFunc.reset(new CoefFunctionBdBKernel<Double>(feFct, 1.));
//    }
//    DefineFieldResult( edFunc, ed );
//    stiffFormCoefs_.insert(edFunc);
//
//
//    // Electric power
//    shared_ptr<ResultInfo> energy( new ResultInfo );
//    energy->resultType = ELEC_POWER;
//    energy->definedOn = ResultInfo::REGION;
//    energy->entryType = ResultInfo::SCALAR;
//    energy->dofNames = "";
//    energy->unit = "W";
//    availResults_.insert( energy );
//    shared_ptr<ResultFunctor> energyFunc;
//    if( isComplex_ ) {
//      /* for harmonic analysis, the effective value is 1/sqrt(2) \hat I \times 1/sqrt(2) \hat U */
//      energyFunc.reset(new EnergyResultFunctor<Complex>(feFct, energy,0.5));
//    } else {
//      energyFunc.reset(new EnergyResultFunctor<Double>(feFct, energy,1.));
//    }
//    resultFunctors_[ELEC_POWER] = energyFunc;
//    stiffFormFunctors_.insert(energyFunc);
    
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
