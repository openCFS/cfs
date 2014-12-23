// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "PiezoCoupling.hh"

#include "PDE/SinglePDE.hh"
#include "CoupledPDE/BasePairCoupling.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Enum.hh"
#include "Materials/BaseMaterial.hh"
#include "Driver/FormsContexts.hh"
#include "Driver/Assemble.hh"

#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"

//transient simulations
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"
#include "Driver/TransientDriver.hh"

// include fespaces
#include "FeBasis/H1/H1Elems.hh"

// new integrator concept
#include "Forms/BiLinForms/ADBInt.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/StrainOperator.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"

namespace CoupledField {



  // ***************
  //   Constructor
  // ***************
  PiezoCoupling::PiezoCoupling( SinglePDE *pde1, SinglePDE *pde2,
                                PtrParamNode paramNode,
                                PtrParamNode infoNode,
                                shared_ptr<SimState> simState,
                                Domain* domain)
    : BasePairCoupling( pde1, pde2, paramNode, infoNode, simState, domain )
  {
    couplingName_ = "piezoDirect";
    materialClass_ = PIEZO;

    // determine subtype from mechanic pde
    pde1_->GetParamNode()->GetValue( "subType", subType_ );

    nonLin_ = false;
    
    // Initialize nonlinearities
    InitNonLin();

  }

  // **************
  //   Destructor
  // **************
  PiezoCoupling::~PiezoCoupling() {
  }


  // *********************
  //   DefineIntegrators
  // *********************
  void PiezoCoupling::DefineIntegrators() {

    // get hold of both feFunctions
    shared_ptr<BaseFeFunction> dispFct = pde1_->GetFeFunction(MECH_DISPLACEMENT);
    shared_ptr<BaseFeFunction> elecFct = pde2_->GetFeFunction(ELEC_POTENTIAL);
    shared_ptr<FeSpace> dispSpace = dispFct->GetFeSpace();
    shared_ptr<FeSpace> elecSpace = elecFct->GetFeSpace();
    
    
    // Global::ComplexPart matType = Global::REAL;
    RegionIdType actRegion;
    BaseMaterial * actSDMat = NULL;
    
    // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      
      // Set current region and material
      actRegion = it->first;
      actSDMat = it->second;
      // matType = Global::REAL;

      //transform the type
      SubTensorType tensorType;
      String2Enum(subType_,tensorType);

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( actRegion );
      
      
      // ==================================================================
      //  STANDARD COUPLING INTEGRATOR
      // ==================================================================
      BaseBDBInt * stiffInt = 
          GetStiffIntegrator( actSDMat, actRegion, complexMatData_[actRegion] );
      stiffInt->SetName("PiezoCouplingInt");
      BiLinFormContext * stiffIntDescr =
          new BiLinFormContext(stiffInt, STIFFNESS );

      stiffIntDescr->SetEntities( actSDList, actSDList );
      stiffIntDescr->SetFeFunctions( dispFct, elecFct );
      stiffIntDescr->SetCounterPart(true);

      assemble_->AddBiLinearForm( stiffIntDescr );

      // remember own bilinearform 
      bdbInts_[actRegion] = stiffInt;
      
    }
  }

  
  
  void PiezoCoupling::DefinePostProcResults() {
    
    shared_ptr<BaseFeFunction> dispFct = pde1_->GetFeFunction(MECH_DISPLACEMENT);
    shared_ptr<BaseFeFunction> elecFct = pde2_->GetFeFunction(ELEC_POTENTIAL);
    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
    MathParser * mp = domain_->GetMathParser();
    
    StdVector<std::string> stressComponents;
    if( subType_ == "3d" ) {
      stressComponents = "xx", "yy", "zz", "yz", "xz", "xy";
    } else if( subType_ == "planeStrain" ) {
      stressComponents = "xx", "yy", "xy";
    } else if( subType_ == "planeStress" ) {
      stressComponents = "xx", "yy", "xy";
    } else if( subType_ == "axi" ) {
      stressComponents = "rr", "zz", "rz", "phiphi";
    }
    
    // === Electric Flux Density ===
    shared_ptr<ResultInfo> flux ( new ResultInfo );
    flux->resultType = ELEC_FLUX_DENSITY;
    flux->SetVectorDOFs(dim_, isaxi_);
    flux->unit = "C/m^2";
    flux->definedOn = ResultInfo::ELEMENT;
    flux->entryType = ResultInfo::VECTOR;
    availResults_.insert( flux );
    
    // The electric flux density in the coupled case calculates as 
    // D = [e]*S + [eps]*E
    // a) electric part -> take from electrostatic PDE
    PtrCoefFct coefElecD = pde2_->GetCoefFct( ELEC_FLUX_DENSITY);
    
    // b) coupling part -> use own ADB-form
    Double cplFactor = 1.0;
    shared_ptr<CoefFunctionFormBased> cplFunc;
    if( isComplex_ ) {
      cplFunc.reset(new CoefFunctionFlux<Complex,true>(dispFct, flux, cplFactor));
    } else {
      cplFunc.reset(new CoefFunctionFlux<Double,true>(dispFct, flux, cplFactor));
    }
    // Build compound coefficient function for flux density
    PtrCoefFct coefFlux = CoefFunction::Generate(mp,part,
                         CoefXprBinOp(mp,coefElecD, cplFunc,
                                      CoefXpr::OP_ADD) );
    DefineFieldResult( coefFlux, flux );


    // ==== Mechanic Stress ===
    shared_ptr<ResultInfo> stress(new ResultInfo);
    stress->resultType = MECH_STRESS;
    stress->dofNames = stressComponents;
    stress->unit =  "N/m^2";
    stress->entryType = ResultInfo::TENSOR;
    stress->definedOn = ResultInfo::ELEMENT;
    availResults_.insert( stress );
    
    // The mechanic stress calculates as 
    // [sigma] = [c]*S - [e]*E
    // a) mechanic  -> take from mechanic PDE
    PtrCoefFct coefMechSigma = pde1_->GetCoefFct( MECH_STRESS );
        
    // b) coupling part [e]*E = - [e] grad(V_p) -> use own ADB-form
    shared_ptr<CoefFunctionFormBased> stressCplFunc;
    // Note: As the "B"-operator of the ADB operator only calculates 
    //       the gradient and we need the electric field (E = - grad(V_p)),
    //       we have to multiply the entries by -1.0. However, due to the 
    //       constitutive law, an additional -1 is needed, so we just
    //       need scaling factor 1.0.
    Double stressCplFactor = 1.0;
    if( isComplex_ ) {
      stressCplFunc.reset(new CoefFunctionFlux<Complex>(elecFct, stress, stressCplFactor));
    } else {
      stressCplFunc.reset(new CoefFunctionFlux<Double>(elecFct, stress, stressCplFactor));
    }
    PtrCoefFct coefStress = CoefFunction::Generate(mp,part,
                                                   CoefXprBinOp(mp,coefMechSigma, stressCplFunc, 
                                                                CoefXpr::OP_ADD) ); 
    DefineFieldResult(coefStress, stress);


    // === Electric Charge Density (surface) ===
    shared_ptr<ResultInfo> chargeD(new ResultInfo);
    chargeD->resultType = ELEC_CHARGE_DENSITY;
    chargeD->dofNames = "";
    chargeD->unit = "C/m^2";
    chargeD->definedOn = ResultInfo::SURF_ELEM;
    chargeD->entryType = ResultInfo::SCALAR;
    availResults_.insert( chargeD );
    // Note: The positive normal direction in this case is defined as the
        //       inward facing one. 
    shared_ptr<CoefFunctionSurf> sChargeDens(new CoefFunctionSurf(true, -1.0, chargeD));
    DefineFieldResult( sChargeDens, chargeD);
    
    // === Electric Charge (integrated) ===
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
      chargeFunc.reset(new ResultFunctorIntegrate<Complex>(sChargeDens, dispFct, charge ) );
    } else {
      chargeFunc.reset(new ResultFunctorIntegrate<Double>(sChargeDens, dispFct, charge ) );
    }
    resultFunctors_[ELEC_CHARGE] = chargeFunc;
    
    
    // ============================
    // Initialize result functors:
    // ============================
    // 1) Loop over all BDB-integrators
    std::map<RegionIdType, BaseBDBInt*>::iterator it = bdbInts_.begin();
    for( ; it != bdbInts_.end(); ++it ) {
      RegionIdType region = it->first;
      BaseBDBInt* bdb = it->second;

      // 2) pass integrators to functors
      cplFunc->AddIntegrator(bdb, region);
      stressCplFunc->AddIntegrator(bdb, region);
      sChargeDens->AddVolumeCoef(region, coefFlux);
    }
  }
  
  

  void PiezoCoupling::DefineAvailResults() {

//
//    // === ELECTRIC CHARGE ===
//    shared_ptr<ResultInfo> charge(new ResultInfo);
//    charge->resultType = ELEC_CHARGE;
//    charge->dofNames = "";
//    charge->unit = "C";
//    charge->definedOn = ResultInfo::SURF_ELEM;
//    charge->entryType = ResultInfo::SCALAR;
//    charge->fctType = shared_ptr<ConstFct>(new ConstFct() );
//    availResults_.insert( charge );
  }

  BaseBDBInt *
  PiezoCoupling::GetStiffIntegrator( BaseMaterial* actSDMat,
                                     RegionIdType regionId,
                                     bool isComplex ) {
    
    // Get region name
    std::string regionName = ptGrid_->GetRegion().ToString( regionId );

    SubTensorType tensorType = NO_TENSOR;
    if( subType_ == "planeStrain" || subType_ == "planeStress" ) {
      tensorType = PLANE_STRAIN;
    } else if( subType_ == "axi" ){
      tensorType = AXI;
    } else if (subType_ == "3d" ){
      tensorType = FULL;
    } else {
      EXCEPTION( "Unknown subtype '" << subType_ << "'" );
    }
    // ------------------------
    //  Obtain linear material
    // ------------------------
    shared_ptr<CoefFunction > curCoef;
    if( isComplex ) {
      curCoef = actSDMat->GetTensorCoefFnc(PIEZO_TENSOR, tensorType, 
                                          Global::COMPLEX, true  );
    } else {
      curCoef = actSDMat->GetTensorCoefFnc(PIEZO_TENSOR, tensorType, 
                                           Global::REAL, true );
    }
    
    // ----------------------------------------
    //  Determine correct stiffness integrator 
    // ----------------------------------------
    BaseBDBInt * integ = NULL;
    if ( isComplex ) {
      if( subType_ == "axi" ) {
        integ = new ADBInt<Complex>(new StrainOperatorAxi<FeH1,Complex>(),
                                    new GradientOperator<FeH1,2,1,Complex>(),
                                    curCoef, 1.0, true );
      } else if( subType_ == "planeStrain" ) {
        integ = new ADBInt<Complex>(new StrainOperator2D<FeH1,Complex>(),
                                    new GradientOperator<FeH1,2,1,Complex>(),
                                    curCoef, 1.0, true);
      } else if( subType_ == "planeStress" ) {
        integ = new ADBInt<Complex>(new StrainOperator2D<FeH1,Complex>(),
                                    new GradientOperator<FeH1,2,1,Complex>(),
                                    curCoef, 1.0, true);
      } else if( subType_ == "3d") {
        integ = new ADBInt<Complex>(new StrainOperator3D<FeH1,Complex>(),
                                    new GradientOperator<FeH1,3,1,Complex>(),
                                    curCoef, 1.0, true);
      } else {
        EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
      }
    }
    else {
      if( subType_ == "axi" ) {
        integ = new ADBInt<Double>(new StrainOperatorAxi<FeH1>(),
                                   new GradientOperator<FeH1,2>(),
                                   curCoef, 1.0, true );
      } else if( subType_ == "planeStrain" ) {
        integ = new ADBInt<Double>(new StrainOperator2D<FeH1>(),
                                   new GradientOperator<FeH1,2>(),
                                   curCoef, 1.0, true);
      } else if( subType_ == "planeStress" ) {
        integ = new ADBInt<Double>(new StrainOperator2D<FeH1>(),
                                   new GradientOperator<FeH1,2>(), 
                                   curCoef, 1.0, true);
      } else if( subType_ == "3d") {
        integ = new ADBInt<Double>(new StrainOperator3D<FeH1>(),
                                   new GradientOperator<FeH1,3>(), 
                                   curCoef, 1.0, true);
      } else {
        EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
      }
    }

    return integ;

  }
  
  void PiezoCoupling::InitTimeStepping(){

    if ( analysisType_ == BasePDE::TRANSIENT ) {
      Double dt;
      dt = dynamic_cast<TransientDriver*>(domain_->GetSingleDriver())
                ->GetDeltaT();

      //in this case we additionally need to define
      //a timestepping for the elecPDE
      shared_ptr<BaseFeFunction> elecFct = pde2_->GetFeFunction(ELEC_POTENTIAL);

      shared_ptr<BaseTimeScheme> elecScheme(new TimeSchemeGLM(GLMScheme::NEWMARK, 0) );

      elecFct->SetTimeScheme(elecScheme);
      elecFct->GetTimeScheme()->Init(elecFct->GetSingleVector(),dt);
    }
  }

}

