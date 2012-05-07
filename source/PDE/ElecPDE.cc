#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>
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
#include "Utils/StdVector.hh"
#include "Driver/SolveSteps/SolveStepElec.hh"
#include "CoupledPDE/PDECoupling.hh"
#include "Driver/Assemble.hh"
#include "Utils/ApproxData.hh"
#include "Utils/SmoothSpline.hh"
#include "Materials/Models/Hysteresis.hh"
#include "Materials/Models/Preisach.hh"
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

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"


namespace CoupledField {

  DECLARE_LOG(elecpde)
  DEFINE_LOG(elecpde, "pde.electrostatic")

  // ***************
  //   Constructor
  // ***************
  ElecPDE::ElecPDE( Grid* aptgrid, PtrParamNode paramNode )
    :SinglePDE( aptgrid, paramNode ) {


    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_          = "electrostatic";
    pdematerialclass_ = ELECTROSTATIC;
    maxTimeDerivOrder_ = 0;
 
    nonLin_    = false;
    nonLinMaterial_ = false;
    isAlwaysStatic_ = true;
    isPiezoCoupled_ = false;

    needSolPrev_ = true;
 
    // Check the subtype of the problem
    paramNode->GetValue("subType", subType_);
  }
  
  void ElecPDE::InitNonLin() {

    SinglePDE::InitNonLin();

    //now do PDE specifics
    std::map<std::string, NonLinType>::iterator it;
    for ( it=nonLinTypes_.begin() ; it != nonLinTypes_.end(); it++ ) {
      if ( (*it).second == HYSTERESIS ) {
        isHysteresis_ = true;
      }
    }

  }

  void ElecPDE::DefineIntegrators() {

    RegionIdType actRegion;
    BaseMaterial * actSDMat = NULL;
    
    // flag for updatedLagrange formulation
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

    // if the pde is piezo-coupled, the electrostatic entries
    // have to be multiplied with -1
    std::string factor = "1.0";
    if ( isPiezoCoupled_ == true )
      factor = "-1.0";  

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
      stiffInt->SetName("LinElecIntegrator");
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
    
    // Define integrators for composite materials
    // (only for subType "flatShell")
    std::map<RegionIdType, Composite>::iterator compIt;
    for( compIt=compositeMaterials_.begin(); compIt!=compositeMaterials_.end();
         compIt++ ) {
      
      REFACTOR;
    }

    // =======================================================================
    // Integrators for NonConforming Interfaces
    // =======================================================================
    ///OUT FOR REFACTOR
    // define integrators for electric impedances
    DefineImpedanceIntegrators();
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
    // Double factor = 1.0;
    // if ( isPiezoCoupled_ )
    //   factor = -1.0;
    
    // =========================
    //  Charges (volume, nodal)
    // =========================
    LOG_DBG(elecpde) << "Reading charges";
    ReadRhsExcitation( "charge", dofNames, ResultInfo::VECTOR, 
                       isComplex_, ent, coef );

    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {

        // ---------------
        //  Nodal Charges 
        // ---------------
        // Nodal charge must be constant
        if( coef[i]->GetDependency() == CoefFunction::GENERAL ) {
          EXCEPTION("Nodal charges must not be spatial dependent");
        }

        UInt numNodes = ent[i]->GetSize();
        if( numNodes > 1 ) {
          // Here we would divide the nodal force by the number of nodes
          // in the list, in order to ensure that the whole force corresponds
          // to the prescribed value. However, this requires modification of 
          // the expressions of the coefficient functions, which depends on real/harm
          // and the specific type (const, timefreq, variable).
          WARN("The charge value will not be divided by the number of nodes and thus "
              << "depends on the number of nodes" );
        }

        lin = new SingleEntryInt(coef[i]);
        lin->SetName("NodalChargeInt");
        LinearFormContext *ctx = new LinearFormContext( lin );
        ctx->SetEntities( ent[i] );
        ctx->SetFeFunction(myFct);
        assemble_->AddLinearForm(ctx);
      } else {
        // --------------------------
        //  Surface / Volume Charges 
        // --------------------------
        if( coef[i]->GetDependency() == CoefFunction::GENERAL ) {
          EXCEPTION("The total prescribed charge must no be spatial dependend");
        }
        // "Divide" the total charge by the volume / surface of the current entity list
        Double volume = ptGrid_->CalcVolumeOfEntityList( ent[i], false );
        Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;  
        coef[i] = CoefFunction::Generate(part, 
                   CoefXprVecScalOp(coef[i], 
                   boost::lexical_cast<std::string>(volume), CoefXpr::OP_DIV) );
        if(isComplex_) {
          lin = new BUIntegrator<IdentityOperator<FeH1>, Complex>(Complex(1.0), coef[i]);
        } else  {
          lin = new BUIntegrator<IdentityOperator<FeH1>, Double>(1.0, coef[i]);
        }
        lin->SetName("ChargeDensityInt");
        LinearFormContext *ctx = new LinearFormContext( lin );
        ctx->SetEntities( ent[i] );
        ctx->SetFeFunction(myFct);
        assemble_->AddLinearForm(ctx);
      }
    } // for

    // ================
    //  CHARGE DENSITY 
    // ================
    LOG_DBG(elecpde) << "Reading charge densities";
    ReadRhsExcitation( "chargeDensity", dofNames, 
                       ResultInfo::VECTOR, isComplex_, ent, coef );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("Charge density must be defined on elements")
      }
      if(isComplex_) {
        lin = new BUIntegrator<IdentityOperator<FeH1>, Complex>(Complex(1.0), coef[i]);
      } else  {
        lin = new BUIntegrator<IdentityOperator<FeH1>, Double>(1.0, coef[i]);
      }
      lin->SetName("ChargeDensityInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
    } // for
    
    // =================================
    //  FLUX DENSITY (VOLUME / SURFACE) 
    // =================================
    LOG_DBG(elecpde) << "Reading prescribed flux density";
    StdVector<std::string> vecDofNames;
    if(dim_ == 3)
      vecDofNames = "x", "y", "z";
    if(dim_ == 2 && !isaxi_)
      vecDofNames = "x", "y";
    if(dim_ == 2 && isaxi_)
      vecDofNames = "r", "z";
    ReadRhsExcitation( "fluxDensity", vecDofNames, 
                       ResultInfo::VECTOR, isComplex_, ent, coef );
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
            lin = new BUIntegrator<IdentityOperatorNormal<FeH1,2>, 
                Complex, true>(Complex(1.0), coef[i], volRegions);
          } else  {
            lin = new BUIntegrator<IdentityOperatorNormal<FeH1,2>, 
                Double, true>(1.0, coef[i], volRegions);
          } 
        }else {
          if(isComplex_) {
            lin = new BUIntegrator<IdentityOperatorNormal<FeH1,3>, 
                Complex, true>(Complex(1.0), coef[i], volRegions);
          } else  {
            lin = new BUIntegrator<IdentityOperatorNormal<FeH1,3>, 
                Double, true>(1.0, coef[i], volRegions);
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
  }

  
  BaseBDBInt * ElecPDE::GetStiffIntegrator( BaseMaterial* actSDMat,
                                            SubTensorType tensorType,
                                            RegionIdType regionId ) {

    BaseBDBInt * integ = NULL;
    shared_ptr<CoefFunction > curCoef = 
        actSDMat->GetTensorCoefFnc( ELEC_PERMITTIVITY,tensorType,
                                   Global::REAL);

    // Note; in the piezoelectric case we have to multiply by -1
    Double factor = 1.0;
    if ( isPiezoCoupled_ )
      factor = -1.0;
    
    if( dim_ == 2 ) {
      integ = new BDBInt<GradientOperator<FeH1,2> >(curCoef, factor );
    } else {
      integ = new BDBInt<GradientOperator<FeH1,3> >(curCoef, factor );
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
    solveStep_ = new SolveStepElec(*this);
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


  // ======================================================
  // COUPLING SECTION
  // ======================================================

  void ElecPDE::InitCoupling(PDECoupling * Coupling)
  {
  REFACTOR;
  }
  


  void ElecPDE::CalcOutputCoupling()
  {
    REFACTOR;
  }


  bool ElecPDE::HasOutput(SolutionType output)
  {
  
    switch (output)
      {
      case ELEC_FORCE_VWP:
        return true;
        break;
      case ELEC_POTENTIAL:
        return true;
        break;
      case ELEC_FIELD_INTENSITY:
        return true;
        break;
      case ELEC_INTERFACE_FORCE:
        return true;
        break;
      default:
        return false;
        break;
      }
    return false;
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
      res1->unit = "V";
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
        dofName += lexical_cast<std::string>(i+1);
        res1->dofNames.Push_back( dofName );
      }
      res1->unit = "V";
      res1->definedOn = ResultInfo::ELEMENT;
      res1->entryType = ResultInfo::SCALAR;
    }
    feFunctions_[ELEC_POTENTIAL]->SetResultInfo(res1);

    res1->SetFeFunction(feFunctions_[ELEC_POTENTIAL]);
    results_.Push_back( res1 );
    availResults_.insert( res1 );
    DefineFieldResult( feFunctions_[ELEC_POTENTIAL], res1 );
    
    // Electric RHS Load
    // create new resultDof object
    shared_ptr<ResultInfo> rhs ( new ResultInfo );
    rhs->resultType = ELEC_RHS_LOAD;
    rhs->dofNames = "";
    rhs->unit = "C";
    rhs->definedOn = results_[0]->definedOn;
    rhs->entryType = ResultInfo::SCALAR;
    availResults_.insert( rhs );

    // ===================================
    // Check for non-conforming interfaces
    // ===================================
    StdVector<std::string> ncIfaceNames, ncIfaceNamesForPDE;
    StdVector<RegionIdType> ncIfaceIds;
    
    LOG_DBG2(elecpde) << "NonMatching: Checking if nonconforming "
                      << "interfaces of PDE exist in domain.";

    PtrParamNode elecPDENCIfaceListNode;
    elecPDENCIfaceListNode = param->GetByVal("sequenceStep", std::string("index"), sequenceStep_)
    ->Get("pdeList/electrostatic/ncInterfaceList", ParamNode::PASS);
    
    if(!elecPDENCIfaceListNode)
      return;

    PtrParamNode domainNCIfaceListNode;
    domainNCIfaceListNode = param->Get("domain")->Get("ncInterfaceList", ParamNode::PASS);

    if(!domainNCIfaceListNode)
    {
      EXCEPTION("No nonmatching interfaces have been specified in domain!");
    }

    ParamNodeList pdeNCIfaceNodes;
    pdeNCIfaceNodes = elecPDENCIfaceListNode->GetList("ncInterface");

    for (UInt i = 0; i < pdeNCIfaceNodes.GetSize(); i++) {
      std::string pdeIfaceName = pdeNCIfaceNodes[i]->Get("name")->As<std::string>();

      PtrParamNode domainIfaceNode = domainNCIfaceListNode
          ->GetByVal("ncInterface", "name", pdeIfaceName, ParamNode::PASS);
      if(!domainIfaceNode)
      {
        LOG_DBG2(elecpde) << "NonMatching: Nonconforming "
        << "interface '" << ncIfaceNames[i]
                                         << "' does not exist in domain.";

        EXCEPTION( "ncInterface referenced from PDE not defined in domain!");
      }

      ncIfaceNamesForPDE.Push_back(pdeIfaceName);
    }
    ptGrid_->GetRegion().Parse(ncIfaceNamesForPDE, ncIfaceIds);

    for (UInt i = 0; i < ncIfaceIds.GetSize(); i++) {
      ncIFaces_.Push_back(ncIfaceIds[i]);
    }

    // In the case of the presence of non-conforming interfaces,
    // a second resultdof object has to be created, which describes the 
    // Lagrange multiplier
    if( ncIFaces_.GetSize() > 0 ) {
      LOG_DBG2(elecpde) << "NonMatching: Defining new ResultDof Lagrange.";
      shared_ptr<ResultInfo> lagr ( new ResultInfo );
      lagr->resultType = LAGRANGE_MULT;
      lagr->dofNames = "l";
      lagr->definedOn = results_[0]->definedOn;
      results_.Push_back( lagr );
    } 
  }
  
  
  void ElecPDE::DefinePostProcResults() {

    shared_ptr<BaseFeFunction> feFct = feFunctions_[ELEC_POTENTIAL];
    
    // === ELECTRIC FIELD INTENSITY ===
    shared_ptr<ResultInfo> ef ( new ResultInfo );
    ef->resultType = ELEC_FIELD_INTENSITY;
    ef->SetVectorDOFs(dim_, isaxi_);
    ef->unit = "V/m";
    ef->definedOn = ResultInfo::ELEMENT;
    ef->entryType = ResultInfo::VECTOR;
    availResults_.insert( ef );
    shared_ptr<CoefFunctionFormBased> eFunc;
    if( isComplex_ ) {
      eFunc.reset(new CoefFunctionBOp<Complex>(feFct, ef, -1.0));
    } else {
      eFunc.reset(new CoefFunctionBOp<Double>(feFct, ef, -1.0));
    }
    DefineFieldResult( eFunc, ef );
    
    // === ELECTRIC FLUX DENSITY ===
    shared_ptr<ResultInfo> flux ( new ResultInfo );
    flux->resultType = ELEC_FLUX_DENSITY;
    flux->SetVectorDOFs(dim_, isaxi_);
    flux->unit = "C/m^2";
    flux->definedOn = ResultInfo::ELEMENT;
    flux->entryType = ResultInfo::VECTOR;
    availResults_.insert( flux );
    shared_ptr<CoefFunctionFormBased> fluxFunc;
    if( isComplex_ ) {
      fluxFunc.reset(new CoefFunctionFlux<Complex>(feFct, flux));
    } else {
      fluxFunc.reset(new CoefFunctionFlux<Double>(feFct, flux));
    }
    DefineFieldResult( fluxFunc, flux );

    // === ELECTRIC SURFACE CHARGE DENSITY ===
    shared_ptr<ResultInfo> chargeD(new ResultInfo);
    chargeD->resultType = ELEC_CHARGE_DENSITY;
    chargeD->dofNames = "";
    chargeD->unit = "C/m^2";
    chargeD->definedOn = ResultInfo::SURF_ELEM;
    chargeD->entryType = ResultInfo::SCALAR;
    availResults_.insert( chargeD );
    // the coefficient function is defined later
    PtrCoefFct sChargeDens;
    
    // === TOTAL ELECTRIC CHARGE ===
    shared_ptr<ResultInfo> charge(new ResultInfo);
    charge->resultType = ELEC_CHARGE;
    charge->dofNames = "";
    charge->unit = "C";
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

    // === ELECTRIC ENERGY DENSITY ===
    shared_ptr<ResultInfo> ed ( new ResultInfo );
    ed->resultType = ELEC_ENERGY_DENSITY;
    ed->dofNames = "";
    ed->unit = "Ws/m^3";
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
    
    // Electric energy
    shared_ptr<ResultInfo> energy( new ResultInfo );
    energy->resultType = ELEC_ENERGY;
    energy->definedOn = ResultInfo::REGION;
    energy->entryType = ResultInfo::SCALAR;
    energy->dofNames = "";
    energy->unit = "Ws";
    availResults_.insert ( energy );
    shared_ptr<ResultFunctor> energyFunc;
    if( isComplex_ ) {
      energyFunc.reset(new EnergyResultFunctor<Complex>(feFct, energy));
    } else {
      energyFunc.reset(new EnergyResultFunctor<Double>(feFct, energy));
    }
    resultFunctors_[ELEC_ENERGY] = energyFunc;
    
    // ============================
    // Initialize result functors:
    // ============================
    // collect the flux volume flux coefficient functionf for each region
    std::map<RegionIdType, PtrCoefFct> dCoefs;
    
    // 1) Loop over all BDB-integrators
    std::map<RegionIdType, BaseBDBInt*>::iterator it = bdbInts_.begin();
    for( ; it != bdbInts_.end(); ++it ) {
      RegionIdType region = it->first;
      BaseBDBInt* bdb = it->second;
      // 2) pass integrators to functors
      eFunc->AddIntegrator(bdb, region);
      fluxFunc->AddIntegrator(bdb, region);
      energyFunc->AddIntegrator(bdb, region);
      edFunc->AddIntegrator(bdb, region);
      
      dCoefs[region] = fluxFunc;
    }
    
    // define missing surface charge density
    if( isComplex_ ) {
      shared_ptr<CoefFunctionSurf<Complex> > dNormal (
          new CoefFunctionSurf<Complex>(true));
      dNormal->SetVolumeCoefs(dCoefs);
      sChargeDens= dNormal;
    } else {
      shared_ptr<CoefFunctionSurf<Double> > dNormal (
          new CoefFunctionSurf<Double>(true));
      dNormal->SetVolumeCoefs(dCoefs);
      sChargeDens = dNormal;
    }  
  }

  
  std::map<SolutionType, shared_ptr<FeSpace> >
  ElecPDE::CreateFeSpaces(const std::string& formulation, PtrParamNode infoNode) {
    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    if(formulation == "default" || formulation == "H1"){
      PtrParamNode potSpaceNode = infoNode->Get("elecPotential");
      crSpaces[ELEC_POTENTIAL] =
        FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1, ptGrid_);
      crSpaces[ELEC_POTENTIAL]->Init(solStrat_);
    }else{
      EXCEPTION("The formulation " << formulation << "of electric PDE is not known!");
    }
    return crSpaces;
  }
}
