#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>
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
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/IdentityOperatorNormal.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"


namespace CoupledField {

  DECLARE_LOG(eleccurrentpde)
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
      mySpace->SetRegionApproximation(actRegion, polyId,integId);

      
      
      // ----- standard real-valued stiffness integrator
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
      bdbInts_[actRegion] = stiffInt;
      
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
  
  void ElecCurrentPDE::DefineNcIntegrators() {
    StdVector< NcInterfaceInfo >::iterator ncIt = ncInterfaces_.Begin(),
                                           endIt = ncInterfaces_.End();
    for ( ; ncIt != endIt; ++ncIt ) {
      switch (ncIt->type) {
      case NC_MORTAR:
        DefineMortarCoupling(ELEC_POTENTIAL, *ncIt);
        break;
      case NC_NITSCHE:
        EXCEPTION("ncInterface of Nitsche type is not implemented for ElecCurrentPDE");
        break;
      default:
        EXCEPTION("Unknown type of ncInterface");
        break;
      }
    }
  }
  
  void ElecCurrentPDE::DefinePostProcResults() {

    shared_ptr<BaseFeFunction> feFct = feFunctions_[ELEC_POTENTIAL];
    
    // === ELECTRIC CURRENT INTENSITY ===
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
    
    // === ELECTRIC FLUX DENSITY ===
    shared_ptr<ResultInfo> flux ( new ResultInfo );
    flux->resultType = ELEC_CURRENT_DENSITY;
    flux->SetVectorDOFs(dim_, isaxi_);
    flux->unit = "A/m^2";
    flux->definedOn = ResultInfo::ELEMENT;
    flux->entryType = ResultInfo::VECTOR;
    shared_ptr<CoefFunctionFormBased> fluxFunc;
    if( isComplex_ ) {
      fluxFunc.reset(new CoefFunctionFlux<Complex>(feFct, flux));
    } else {
      fluxFunc.reset(new CoefFunctionFlux<Double>(feFct, flux));
    }
    DefineFieldResult( fluxFunc, flux );
    stiffFormCoefs_.insert(fluxFunc);


    // === ELECTRIC POWER DENSITY ===
    shared_ptr<ResultInfo> ed ( new ResultInfo );
    ed->resultType = ELEC_POWER_DENSITY;
    ed->dofNames = "";
    ed->unit = "W/m^3";
    ed->definedOn = ResultInfo::ELEMENT;
    ed->entryType = ResultInfo::SCALAR;
    availResults_.insert( ed );
    shared_ptr<CoefFunctionFormBased> edFunc;
    if( isComplex_ ) {
      edFunc.reset(new CoefFunctionBdBKernel<Complex>(feFct, 0.5));
    } else {
      edFunc.reset(new CoefFunctionBdBKernel<Double>(feFct, 0.5));
    }
    DefineFieldResult( edFunc, ed );
    stiffFormCoefs_.insert(edFunc);
    
    // Electric power
    shared_ptr<ResultInfo> energy( new ResultInfo );
    energy->resultType = ELEC_POWER;
    energy->definedOn = ResultInfo::REGION;
    energy->entryType = ResultInfo::SCALAR;
    energy->dofNames = "";
    energy->unit = "W";
    availResults_.insert( energy );
    shared_ptr<ResultFunctor> energyFunc;
    if( isComplex_ ) {
      energyFunc.reset(new EnergyResultFunctor<Complex>(feFct, energy,0.5));
    } else {
      energyFunc.reset(new EnergyResultFunctor<Double>(feFct, energy,0.5));
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
