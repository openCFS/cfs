// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "MagnetoStrictCoupling.hh"

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
#include "Forms/Operators/CurlOperator.hh"
#include "Forms/Operators/StrainOperator.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"

namespace CoupledField {



  // ***************
  //   Constructor
  // ***************
  MagnetoStrictCoupling::MagnetoStrictCoupling( SinglePDE *pde1, SinglePDE *pde2,
                                PtrParamNode paramNode,
                                PtrParamNode infoNode,
                                shared_ptr<SimState> simState,
                                Domain* domain)
    : BasePairCoupling( pde1, pde2, paramNode, infoNode, simState, domain )
  {
    couplingName_ = "magnetoStrictDirect";
    materialClass_ = MAGNETOSTRICTIVE;
    // determine subtype from mechanic pde
    pde1_->GetParamNode()->GetValue( "subType", subType_ );

    nonLin_ = false;
    
    // Initialize nonlinearities
    InitNonLin();
  }

  // **************
  //   Destructor
  // **************
  MagnetoStrictCoupling::~MagnetoStrictCoupling() {
  }


  // *********************
  //   DefineIntegrators
  // *********************
  void MagnetoStrictCoupling::DefineIntegrators() {

    // get hold of both feFunctions
    shared_ptr<BaseFeFunction> dispFct = pde1_->GetFeFunction(MECH_DISPLACEMENT);
    shared_ptr<BaseFeFunction> magFct = pde2_->GetFeFunction(MAG_POTENTIAL);
    shared_ptr<FeSpace> dispSpace = dispFct->GetFeSpace();
    shared_ptr<FeSpace> magSpace = magFct->GetFeSpace();
    
    
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
      stiffInt->SetName("MagnetoStrictCouplingInt");
      BiLinFormContext * stiffIntDescr =
          new BiLinFormContext(stiffInt, STIFFNESS );

      stiffIntDescr->SetEntities( actSDList, actSDList );
      stiffIntDescr->SetFeFunctions( dispFct, magFct );
      stiffIntDescr->SetCounterPart(true);

      assemble_->AddBiLinearForm( stiffIntDescr );

      // remember own bilinearform 
      bdbInts_[actRegion] = stiffInt;
      
    }
  }

  
  
  void MagnetoStrictCoupling::DefinePostProcResults() {
    
    shared_ptr<BaseFeFunction> dispFct = pde1_->GetFeFunction(MECH_DISPLACEMENT);
    shared_ptr<BaseFeFunction> magFct = pde2_->GetFeFunction(MAG_POTENTIAL);
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
   
    // === MAGNETIC FIELD INTENSITY ===
    shared_ptr<ResultInfo> magIntens(new ResultInfo);
    magIntens->resultType = MAG_FIELD_INTENSITY;
    magIntens->SetVectorDOFs(dim_, isaxi_);
    magIntens->unit = "A/m";
    magIntens->definedOn = ResultInfo::ELEMENT;
    magIntens->entryType = ResultInfo::VECTOR;
    availResults_.insert( magIntens );
        
    // The magnetic field intensity in the coupled case calculates as 
    // H = -[h]*S + [nu]*B
    // a) magnetic part -> take from magnetic PDE
    PtrCoefFct coefMagH = pde2_->GetCoefFct( MAG_FIELD_INTENSITY);
   
    Double cplFactor = 1.0;
    shared_ptr<CoefFunctionFormBased> cplFunc;
    if( isComplex_ ) {
      cplFunc.reset(new CoefFunctionFlux<Complex,true>(dispFct, magIntens, cplFactor));
    } else {
      cplFunc.reset(new CoefFunctionFlux<Double,true>(dispFct, magIntens, cplFactor));
    }
  
    // Build compound coefficient function for field intensity
    PtrCoefFct coefIntens = CoefFunction::Generate(mp,part,
                         CoefXprBinOp(mp,coefMagH, cplFunc,
                                      CoefXpr::OP_SUB) );

    DefineFieldResult( coefIntens, magIntens );

    // ==== Mechanic Stress ===
    shared_ptr<ResultInfo> stress(new ResultInfo);
    stress->resultType = MECH_STRESS;
    stress->dofNames = stressComponents;
    stress->unit =  "N/m^2";
    stress->entryType = ResultInfo::TENSOR;
    stress->definedOn = ResultInfo::ELEMENT;
    availResults_.insert( stress );
    
    // The mechanic stress calculates as 
    // [sigma] = [c]*S - [h]*B
    // a) mechanic  -> take from mechanic PDE
    PtrCoefFct coefMechSigma = pde1_->GetCoefFct( MECH_STRESS );
        
    // b) coupling part [h]*B = [h] rot(A) -> use own ADB-form
    shared_ptr<CoefFunctionFormBased> stressCplFunc;
   Double stressCplFactor = 1.0;
    if( isComplex_ ) {
      stressCplFunc.reset(new CoefFunctionFlux<Complex>(magFct, stress, stressCplFactor));
    } else {
      stressCplFunc.reset(new CoefFunctionFlux<Double>(magFct, stress, stressCplFactor));
    }
    PtrCoefFct coefStress = CoefFunction::Generate(mp,part,
                            CoefXprBinOp(mp,coefMechSigma, stressCplFunc, 
                                         CoefXpr::OP_SUB ) ); 
    DefineFieldResult(coefStress, stress);

// keep as reference 
/*
    // === Electric Charge Density (surface) ===
    shared_ptr<ResultInfo> chargeD(new ResultInfo);
    chargeD->resultType = ELEC_CHARGE_DENSITY;
    chargeD->dofNames = "";
    chargeD->unit = "C/m^2";
    chargeD->definedOn = ResultInfo::SURF_ELEM;
    chargeD->entryType = ResultInfo::SCALAR;
    availResults_.insert( chargeD );
    shared_ptr<CoefFunctionSurf> sChargeDens(new CoefFunctionSurf(true, chargeD));
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
*/    
    
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
      //sChargeDens->AddVolumeCoef(region, coefFlux);
    }
  }
  
  

  void MagnetoStrictCoupling::DefineAvailResults() {
  }

  BaseBDBInt *
  MagnetoStrictCoupling::GetStiffIntegrator( BaseMaterial* actSDMat,
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
      curCoef = actSDMat->GetTensorCoefFnc(MAGNETOSTRICTION_TENSOR_h, tensorType, 
                                          Global::COMPLEX, true  );
    } else {
      curCoef = actSDMat->GetTensorCoefFnc(MAGNETOSTRICTION_TENSOR_h, tensorType, 
                                           Global::REAL, true );
    }

    // ----------------------------------------
    //  Determine correct stiffness integrator 
    // ----------------------------------------

    // NOTE: the used set of governing equations
    //   H = -[h]*S + [nu]*B
    //   [sigma] = [c]*S - [h]*B
    // couple with a negative sign, so we have to use a factor of -1 for the coupling matrices 
    double factor = -1.0;

    BaseBDBInt * integ = NULL;
    if ( isComplex ) {
      if( subType_ == "axi" ) {
        integ = new ADBInt<Complex>(new StrainOperatorAxi<FeH1,Complex>(),
						new CurlOperatorAxi<Complex>(),
                                    //new GradientOperator<FeH1,2,1,Complex>(),
                                    curCoef, factor, true );
      } else if( subType_ == "planeStrain" ) {
			
        integ = new ADBInt<Complex>(new StrainOperator2D<FeH1,Complex>(),
						new CurlOperator<FeH1,2, Complex>(),
                                    //new GradientOperator<FeH1,2,1,Complex>(),
                                    curCoef, factor, true);
      } else if( subType_ == "planeStress" ) {

        integ = new ADBInt<Complex>(new StrainOperator2D<FeH1,Complex>(),
						new CurlOperator<FeH1,2, Complex>(),
                                    //new GradientOperator<FeH1,2,1,Complex>(),
                                    curCoef, factor, true);
      } else if( subType_ == "3d") {
        integ = new ADBInt<Complex>(new StrainOperator3D<FeH1,Complex>(),
						new CurlOperator<FeH1,3, Complex>(),
                                    //new GradientOperator<FeH1,3,1,Complex>(),
                                    curCoef, factor, true);
      } else {
        EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
      }
    }
    else {
      if( subType_ == "axi" ) {
        integ = new ADBInt<Double>(new StrainOperatorAxi<FeH1>(),
        				     new CurlOperatorAxi<Double>(),
                                   //new GradientOperator<FeH1,2>(),
                                   curCoef, factor, true );
      } else if( subType_ == "planeStrain" ) {
		
        integ = new ADBInt<Double>(new StrainOperator2D<FeH1>(),
        				     new CurlOperator<FeH1,2, Double>(),
                                   //new GradientOperator<FeH1,2>(),
                                   curCoef, factor, true);
      } else if( subType_ == "planeStress" ) {
        integ = new ADBInt<Double>(new StrainOperator2D<FeH1>(),
        				     new CurlOperator<FeH1,2, Double>(),
                                   //new GradientOperator<FeH1,2>(),
                                   curCoef, factor, true);
      } else if( subType_ == "3d") {
        integ = new ADBInt<Double>(new StrainOperator3D<FeH1>(),
        				     new CurlOperator<FeH1,3, Double>(),
                                   //new GradientOperator<FeH1,2>(), 
                                   curCoef, factor, true);
      } else {
        EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
      }
    }

    return integ;

  }
  
  void MagnetoStrictCoupling::InitTimeStepping(){

    if ( analysisType_ == BasePDE::TRANSIENT ) {
      Double dt;
      dt = dynamic_cast<TransientDriver*>(domain_->GetSingleDriver())
                ->GetDeltaT();

      //in this case we additionally need to define
      //a timestepping for the magPDE
      shared_ptr<BaseFeFunction> magFct = pde2_->GetFeFunction(MAG_POTENTIAL);

      shared_ptr<BaseTimeScheme> magScheme(new TimeSchemeGLM(GLMScheme::NEWMARK, 0) );

      magFct->SetTimeScheme(magScheme);
      magFct->GetTimeScheme()->Init(magFct->GetSingleVector(),dt);
    }
  }

}

