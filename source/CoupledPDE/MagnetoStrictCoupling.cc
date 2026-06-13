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
#include "Forms/LinForms/BUInt.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/CurlOperator.hh"
#include "Forms/Operators/StrainOperator.hh"

//#include "Forms/Operators/SurfaceOperators.hh"
//#include "Forms/Operators/SurfaceNormalStressOperator.hh"
//#include "Forms/Operators/SurfaceNormalFluxDensityOperator.hh"

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
		
    // at this point, the single pdes are not set up yet
    // > we cannot find out if something is nonlinear or not
    // > postpone to definePostProcResults
    nonLin_ = false;
    isHyst_ = false;
    irrStrainsSet_ = false;

    //    old
    //		nonLin_ = pde1_->IsNonLin() || pde2_->IsNonLin();
    //
    //    isHyst_ = false;
    //    if(pde1_->IsHysteresis()){
    //      EXCEPTION("Currently only the elec PDE may be hysteretic");
    //    } else if(pde2_->IsHysteresis()){
    //      isHyst_ = true;
    //    }
    
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

    // at this point, the single PDEs have been set up and we can no check
    // if we are nonlinear or hysteretic
    nonLin_ = pde1_->IsNonLin() || pde2_->IsNonLin();
    // check if elec pde is hysteretic
    isHyst_ = false;
    if(pde1_->IsHysteresis()){
      EXCEPTION("Currently only the magnetic PDE may be hysteretic");
    } else if(pde2_->IsHysteresis()){
      isHyst_ = true;
    }
		
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

			// copied from PiezoCoupling.cc
			if(isHyst_){
				//				std::cout << "Hyst case -> check if region is hyst" << std::endl;
				BaseBDBInt* mechToMagInt = NULL;
				BaseBDBInt* magToMechInt = NULL;
				
				bool regionIsHyst = GetStiffIntegratorHyst( actSDMat,actRegion,complexMatData_[actRegion],&magToMechInt, &mechToMagInt);
				
				if(regionIsHyst){
					//					std::cout << "Region is hyst" << std::endl;
					
					//					if(mechToElecInt == NULL){
					//						std::cout << "mechToElecInt == NULL!!" << std::endl;
					//					}
					//					if(elecToMechInt == NULL){
					//						std::cout << "elecToMechInt == NULL!!" << std::endl;
					//					}
					
					magToMechInt->SetName("MagStrictCouplingIntMagToMech");
					BiLinFormContext * stiffIntDescrMagMech =
                  new BiLinFormContext(magToMechInt, STIFFNESS );
					
					stiffIntDescrMagMech->SetEntities( actSDList, actSDList );
					stiffIntDescrMagMech->SetFeFunctions( dispFct, magFct );
					stiffIntDescrMagMech->SetCounterPart(false);
					
					assemble_->AddBiLinearForm( stiffIntDescrMagMech );
					//					std::cout << "first part put into place" << std::endl;
					// remember own bilinearform 
					bdbInts_[actRegion] = magToMechInt;
					
					mechToMagInt->SetName("MagStrictCouplingIntMechToMag");
					BiLinFormContext * stiffIntDescrMechMag =
                  new BiLinFormContext(mechToMagInt, STIFFNESS );
					
					stiffIntDescrMechMag->SetEntities( actSDList, actSDList );
					stiffIntDescrMechMag->SetFeFunctions( magFct, dispFct );
					stiffIntDescrMechMag->SetCounterPart(false);
					
					assemble_->AddBiLinearForm( stiffIntDescrMechMag );
					
					// remember own bilinearform 
					bdbIntsCounterpart_[actRegion] = mechToMagInt;
					considerCounterpart_ = true;
					//					std::cout << "everything put into place" << std::endl;
				} else {
					// std linear case
					magToMechInt->SetName("MagStrictCouplingInt");
					BiLinFormContext * stiffIntDescr =
                  new BiLinFormContext(magToMechInt, STIFFNESS );
					
					stiffIntDescr->SetEntities( actSDList, actSDList );
					stiffIntDescr->SetFeFunctions( dispFct, magFct );
					stiffIntDescr->SetCounterPart(true);
					
					assemble_->AddBiLinearForm( stiffIntDescr );
					
					// remember own bilinearform 
					bdbInts_[actRegion] = magToMechInt;
				}
				
			} else {
				
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
		
		DefineRhsLoadIntegrators();
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
    // if counterpart was stored separately, it is already transposed, i.e. we
		// do not need have to tell cplFunc that we want to transpose
		bool transposed = !considerCounterpart_;
		if(transposed){
      if( isComplex_ ) {
        cplFunc.reset(new CoefFunctionFlux<Complex,true>(dispFct, magIntens, cplFactor));
      } else {
        cplFunc.reset(new CoefFunctionFlux<Double,true>(dispFct, magIntens, cplFactor));
      }
		} else {
      if( isComplex_ ) {
        cplFunc.reset(new CoefFunctionFlux<Complex>(dispFct, magIntens, cplFactor));
      } else {
        cplFunc.reset(new CoefFunctionFlux<Double>(dispFct, magIntens, cplFactor));
      }
		}

    // Build compound coefficient function for field intensity
    PtrCoefFct coefIntens = CoefFunction::Generate(mp,part,
            CoefXprBinOp(mp,coefMagH, cplFunc,
            CoefXpr::OP_SUB) );
    // H = nuB - M - hSr + hSi
		// in hyst case, we have to add [h]S_irr!
    if( !isHyst_){
      DefineFieldResult( coefIntens, magIntens );
    } else {
      if(irrStrainsSet_ == false){
        EXCEPTION("Hyst coef functions not initialized yet");
      }
      PtrCoefFct coefIntens_full = CoefFunction::Generate(mp,part,
              CoefXprBinOp(mp,coefIntens, coupledIrrStrains_,
              CoefXpr::OP_ADD) );
      DefineFieldResult( coefIntens_full, magIntens );
    }

		// === MAGNETIC ENERGY ===
		shared_ptr<ResultInfo> energy(new ResultInfo);
		energy->resultType = MAG_ENERGY;
		energy->dofNames = "";
		energy->unit = "Ws";
		energy->definedOn = ResultInfo::REGION;
		energy->entryType = ResultInfo::SCALAR;
		availResults_.insert( energy );
		
		shared_ptr<ResultFunctor> energyFunc;
		
		// calculate energy by integrate(1/2*H*B)
		PtrCoefFct coefMagB = pde2_->GetCoefFct( MAG_FLUX_DENSITY );
		
		// for comparison with pure magnetic energy take coefMagH time coefMagB
		PtrCoefFct coefBH = CoefFunction::Generate(mp, part, CoefXprBinOp(mp, coefMagH, coefMagB, CoefXpr::OP_MULT));
		//PtrCoefFct coefBH = CoefFunction::Generate(mp, part, CoefXprBinOp(mp, coefIntens, coefMagB, CoefXpr::OP_MULT));
		//std::cerr << "BH is " << coefBH->ToString();
		//std::cerr << std::endl;
		
		PtrCoefFct energyDens = CoefFunction::Generate( mp, part, CoefXprBinOp(mp, coefBH, "0.5", CoefXpr::OP_MULT) );
		
		//std::cerr << "energyDens is " << energyDens2->ToString();
		//std::cerr << std::endl;
		
		if( isComplex_ ) {
			energyFunc.reset(new ResultFunctorIntegrate<Complex>(energyDens, magFct, energy ) );
		} else {
			energyFunc.reset(new ResultFunctorIntegrate<Double>(energyDens, magFct, energy ) );
		}
		
		resultFunctors_[MAG_ENERGY] = energyFunc;
		
    // === irrevesible strains and stresses ===
    // check if at least one region has hysteretic behavior
    if ( isHyst_){
      // Note: in magPDE we have to setup storage for hyst output in definePostProcResults
      // as this function is called prior to DefineIntegrators
      // in magstrictCase, the postprocresults will be defined after the integrators, so
      // here we have to setup the storage earlier in GetStiffIntegratorHyst
      shared_ptr<ResultInfo> irrStrainInfo(new ResultInfo);
      irrStrainInfo->resultType = MECH_IRR_STRAIN;
      irrStrainInfo->dofNames = stressComponents;
      irrStrainInfo->unit =  "";
      irrStrainInfo->entryType = ResultInfo::TENSOR;
      // strains retruns them in voigt vector
      irrStrainInfo->definedOn = ResultInfo::ELEMENT;
      DefineFieldResult( irrStrains_, irrStrainInfo );
      availResults_.insert( irrStrainInfo );
      shared_ptr<ResultInfo> irrStressInfo(new ResultInfo);
      irrStressInfo->resultType = MECH_IRR_STRESS;
      irrStressInfo->dofNames = stressComponents;
      irrStressInfo->unit =  "N/m^2";
      irrStressInfo->entryType = ResultInfo::TENSOR;
      // strains retruns them in voigt vector
      irrStressInfo->definedOn = ResultInfo::ELEMENT;
      DefineFieldResult( irrStresses_, irrStressInfo );
      availResults_.insert( irrStressInfo );
    }

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

    // in case of hysteresis, stresses are computed via
    // [sigma] = [c]*S_revesible - [h]*B
    // > from mech PDE we get [sigma] = [c]*Bu = [c]*S_total = [c]*(S_reversible + S_irreversible)
    // > subtract [c]*S_irreversible
    PtrCoefFct coefStress;
    if ( isHyst_ ){
      PtrCoefFct coefStress_full = CoefFunction::Generate(mp,part,
              CoefXprBinOp(mp,coefMechSigma, stressCplFunc,
              CoefXpr::OP_ADD) );
      coefStress = CoefFunction::Generate(mp,part,
              CoefXprBinOp(mp,coefStress_full, irrStressesVector_,
              CoefXpr::OP_SUB) );
    } else {
      coefStress = CoefFunction::Generate(mp,part,
              CoefXprBinOp(mp,coefMechSigma, stressCplFunc,
              CoefXpr::OP_ADD) );
    }

    //    OLD
    //    DefineFieldResult(coefStress, stress);
    //
    //    PtrCoefFct coefStress = CoefFunction::Generate(mp,part,
    //			CoefXprBinOp(mp,coefMechSigma, stressCplFunc,
    //			CoefXpr::OP_SUB ) );
    //    DefineFieldResult(coefStress, stress);
    
    // ============================
    // Initialize result functors:
    // ============================
    // 1) Loop over all BDB-integrators
    std::map<RegionIdType, BaseBDBInt*>::iterator it = bdbInts_.begin();
    for( ; it != bdbInts_.end(); ++it ) {
      RegionIdType region = it->first;
			
			if(considerCounterpart_){
				// treat magMech and mechMag separately
				// bdbInts_ stores bdbIntegrator for mechanical pde (magToMechInt)
				BaseBDBInt* bdb = it->second;
				stressCplFunc->AddIntegrator(bdb, region);
				
				// bdbIntsCounterpart_ stores bdeIntegrator for electric pde (mechToElecInt)
				BaseBDBInt* bdbCounterpart = bdbIntsCounterpart_[region];
				cplFunc->AddIntegrator(bdbCounterpart, region);
        //				std::cout << "cpl success" << std::endl;
			} else {
				BaseBDBInt* bdb = it->second;

				// 2) pass integrators to functors
				cplFunc->AddIntegrator(bdb, region);
				stressCplFunc->AddIntegrator(bdb, region);
			}
    }
  }
  
  void MagnetoStrictCoupling::DefineAvailResults() {
  }

  // for hysteresis we need rhs loads that have to be defined in coupling
	// as it needs informattion about both pdes
  void MagnetoStrictCoupling::DefineRhsLoadIntegrators() {

    //		std::cout << "PiezoCoupling - DefineLOADIntegrators" << std::endl;
		
		shared_ptr<BaseFeFunction> dispFct = pde1_->GetFeFunction(MECH_DISPLACEMENT);
    shared_ptr<BaseFeFunction> magFct = pde2_->GetFeFunction(MAG_POTENTIAL);
      LinearForm * lin = NULL;
		LinearForm * linMag = NULL;
    // Flag, if coefficient function lives on updated geoemtry
    bool coefUpdateGeo = true;
		
		SubTensorType tensorType = NO_TENSOR;
    if( subType_ == "planeStrain" || subType_ == "planeStress" ) {
      tensorType = PLANE_STRAIN;
    } else if( subType_ == "axi" ){
      tensorType = AXI;
    } else if (subType_ == "3d" || subType_ == "2.5d"){
      tensorType = FULL;
    } else {
      EXCEPTION( "Unknown subtype '" << subType_ << "'" );
    }

		if(isHyst_){
			// get all regions with hysteresis information from elecPDE
			shared_ptr<CoefFunctionMulti> hysteresisCoefs = pde2_->GetHystCoefs();
			BaseMaterial * actSDMat = NULL;
			
      std::map<RegionIdType,PtrCoefFct > regionCoefs = hysteresisCoefs->GetRegionCoefs();
      std::map<RegionIdType, shared_ptr<CoefFunction> > ::iterator it;
      for( it = regionCoefs.begin(); it != regionCoefs.end(); it++) {

        // get regionIdType
        RegionIdType curReg = it->first;
				PtrCoefFct regionHystOperator = it->second;
				actSDMat = materials_[curReg];
				
        // get SDList
        shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
        actSDList->SetRegion( curReg );


				// per default we evaluate at each integration point
				// this behavior can be changed by setting a flag in coefFncHyst
        bool fullevaluation = true;
				Double factor = 1.0;
				
				// ------------------------
				//  Obtain linear material
				// ------------------------
				shared_ptr<CoefFunction > couplingCoef;
				if( isComplex_ ) {
					couplingCoef = actSDMat->GetTensorCoefFnc(MAGNETOSTRICTION_TENSOR_h, tensorType, 
                  Global::COMPLEX, false  );
				} else {
					couplingCoef = actSDMat->GetTensorCoefFnc(MAGNETOSTRICTION_TENSOR_h, tensorType, 
                  Global::REAL, false );
				}	
				
        //				std::cout << "couplCoef" << std::endl;
        //				std::cout << couplingCoef->ToString() << std::endl;
				
				// extract stiffness tensor from mech pde
				std::map<RegionIdType, BaseMaterial*>  mechMat;
				mechMat = pde1_->GetMaterialData();

				shared_ptr<CoefFunction > stiffCoef;
				stiffCoef = mechMat[curReg]->GetTensorCoefFnc(MECH_STIFFNESS_TENSOR, tensorType, 
                Global::REAL, true );
				
        // important: pass information about stiffness and coupling tensor to hystOperator
        // (if not already done)
        regionHystOperator->SetElastAndCouplTensor(stiffCoef, couplingCoef);

				// mech RHS
				shared_ptr<CoefFunction> mechRHS = regionHystOperator->GenerateRHSCoefFnc("IrrStressForMechPDE");
				
				if ( isComplex_ ) {
					if( subType_ == "axi" ) {
						lin = new BUIntegrator<Complex>( new StrainOperatorAxi<FeH1,Complex>(),
                    Complex(factor),mechRHS,  coefUpdateGeo, fullevaluation);

					} else if( subType_ == "planeStrain" ) {
						lin = new BUIntegrator<Complex>( new StrainOperator2D<FeH1,Complex>(),
                    Complex(factor),mechRHS,  coefUpdateGeo, fullevaluation);

					} else if( subType_ == "planeStress" ) {
						lin = new BUIntegrator<Complex>( new StrainOperator2D<FeH1,Complex>(),
                    Complex(factor),mechRHS,  coefUpdateGeo, fullevaluation);

					} else if( subType_ == "3d") {
						lin = new BUIntegrator<Complex>( new StrainOperator3D<FeH1,Complex>(),
                    Complex(factor),mechRHS,  coefUpdateGeo, fullevaluation);

					} else if( subType_ == "2.5d") {
						lin = new BUIntegrator<Complex>( new StrainOperator2p5D<FeH1,Complex>(),
                    Complex(factor),mechRHS,  coefUpdateGeo, fullevaluation);
					} else {
						EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
					}
				} else {
					if( subType_ == "axi" ) {
						lin = new BUIntegrator<Double>( new StrainOperatorAxi<FeH1>(),
                    (factor),mechRHS, coefUpdateGeo, fullevaluation);

					} else if( subType_ == "planeStrain" ) {
						lin = new BUIntegrator<Double>( new StrainOperator2D<FeH1>(),
                    (factor),mechRHS, coefUpdateGeo, fullevaluation);

					} else if( subType_ == "planeStress" ) {
						lin = new BUIntegrator<Double>( new StrainOperator2D<FeH1>(),
                    (factor),mechRHS, coefUpdateGeo, fullevaluation);
					} else if( subType_ == "3d") {
						lin = new BUIntegrator<Double>( new StrainOperator3D<FeH1>(),
                    (factor),mechRHS, coefUpdateGeo, fullevaluation);
					} else if( subType_ == "2.5d") {
						lin = new BUIntegrator<Double>( new StrainOperator2p5D<FeH1>(),
                    (factor),mechRHS, coefUpdateGeo, fullevaluation);
						} else {
						EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
					}
				}

        lin->SetName("IrrStressForMechPDE");
        lin->SetSolDependent();
        LinearFormContext *ctx = new LinearFormContext( lin );
        ctx->SetEntities( actSDList );
        ctx->SetFeFunction(dispFct);
        
				assemble_->AddLinearForm(ctx);
        // Add entity list will add nothing, if entities were already assigned
        dispFct->AddEntityList(actSDList);
				
        // ++++++++
				// mag RHS
        // > note: magPDE will get multiplied by -1 in case of magstrict coupling
        // > rhs term will also switch sign
        factor = -1.0;
				shared_ptr<CoefFunction> magRHS = regionHystOperator->GenerateRHSCoefFnc("CoupledIrrStrainForMagPDE");
				
				if(isComplex_) {
          if( dim_ == 2 ) {
            if(isaxi_){
              linMag = new BUIntegrator<Complex>( new CurlOperatorAxi<Double>(),
                      Complex(factor),magRHS,  coefUpdateGeo, fullevaluation);
            } else {
              linMag = new BUIntegrator<Complex>( new CurlOperator<FeH1,2,Double>(),
                      Complex(factor),magRHS,  coefUpdateGeo, fullevaluation);
            }
          } else {
            linMag = new BUIntegrator<Complex>( new CurlOperator<FeH1,3,Double>(),
                    Complex(factor),magRHS,  coefUpdateGeo, fullevaluation);
          }
        } else  {
          if( dim_ == 2 ) {
            if(isaxi_){
              linMag = new BUIntegrator<Double>( new CurlOperatorAxi<Double>(),
                      factor,magRHS,  coefUpdateGeo, fullevaluation);
            } else {
              linMag = new BUIntegrator<Double>( new CurlOperator<FeH1,2,Double>(),
                      factor,magRHS,  coefUpdateGeo, fullevaluation);
            }
          } else {
            linMag = new BUIntegrator<Double>( new CurlOperator<FeH1,3,Double>(),
                    factor,magRHS,  coefUpdateGeo, fullevaluation);
          }
        }

        linMag->SetName("CoupledIrrStrainForMagPDE");
        linMag->SetSolDependent();
        LinearFormContext *ctxMag = new LinearFormContext( linMag );
        ctxMag->SetEntities( actSDList );
        ctxMag->SetFeFunction(magFct);
        
				assemble_->AddLinearForm(ctxMag);
        // Add entity list will add nothing, if entities were already assigned
        magFct->AddEntityList(actSDList);

      }
    }
  }
	
  bool MagnetoStrictCoupling::GetStiffIntegratorHyst( BaseMaterial* actSDMat,
          RegionIdType regionId, bool isComplex,
          BaseBDBInt** magToMechInt, BaseBDBInt** mechToMagInt) {
    
    SubTensorType tensorType = NO_TENSOR;
    UInt numStrainEntries = 3;
    if( subType_ == "planeStrain" || subType_ == "planeStress" ) {
      tensorType = PLANE_STRAIN;
    } else if( subType_ == "axi" ){
      tensorType = AXI;
      numStrainEntries = 4;
    } else if (subType_ == "3d" || subType_ == "2.5d"){
      tensorType = FULL;
      numStrainEntries = 6;
    } else {
      EXCEPTION( "Unknown subtype '" << subType_ << "'" );
    }

    // setup arrays for output > has to be done here as we insert further below
    if(!irrStrainsSet_){
      irrStrains_.reset(new CoefFunctionMulti(CoefFunction::TENSOR, dim_,dim_,isComplex_));
      irrStresses_.reset(new CoefFunctionMulti(CoefFunction::TENSOR, dim_,dim_,isComplex_));
      irrStressesVector_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, numStrainEntries,1,isComplex_));
      coupledIrrStrains_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_,1,isComplex_));
      irrStrainsSet_ = true;
    }
       // ------------------------
    //  Obtain linear material
    // ------------------------		
		shared_ptr<CoefFunction > couplingCoef;
    if( isComplex ) {
      couplingCoef = actSDMat->GetTensorCoefFnc(MAGNETOSTRICTION_TENSOR_h, tensorType, 
              Global::COMPLEX, false  );
    } else {
      couplingCoef = actSDMat->GetTensorCoefFnc(MAGNETOSTRICTION_TENSOR_h, tensorType, 
              Global::REAL, false );
    }	
		
    //		std::cout << "couplingCoefTensor: " << couplingCoef->ToString() << std::endl;
    //
		// check if current region is hysteretic
		std::map<RegionIdType, StdVector<NonLinType> > nonLinMap;
		nonLinMap = pde2_->GetNonLinRegionTypes();

		const StdVector<NonLinType>& nonLinTypes = nonLinMap[regionId];
    if ( nonLinTypes.Find(HYSTERESIS) != -1 ){

			// extract stiffness tensor from mech pde
			std::map<RegionIdType, BaseMaterial*>  mechMat;
			mechMat = pde1_->GetMaterialData();

			shared_ptr<CoefFunction > stiffCoef;
			stiffCoef = mechMat[regionId]->GetTensorCoefFnc(MECH_STIFFNESS_TENSOR, tensorType, 
              Global::REAL, true );
			// get hyst operator for current region from elec pde
			PtrCoefFct hystOperator;
			shared_ptr<CoefFunctionMulti> hystCoefs;
			hystCoefs = pde2_->GetHystCoefs();
			hystOperator =hystCoefs->GetRegionCoef(regionId);
      
		  // important: pass information about stiffness and coupling tensor to hystOperator
      // (if not already done)
      hystCoefs->SetElastAndCouplTensor(stiffCoef, couplingCoef);

			// create hyst material
			// note: curCoef holds the coupling tensor
			// hyst case might be unsymmetric
			// > in case of deltaFormulation, elecToMech holds an addition deltaS/deltaE
			PtrCoefFct mechToMag = hystOperator->GenerateMatCoefFnc("CouplingMechToMag");
			PtrCoefFct magToMech = hystOperator->GenerateMatCoefFnc("CouplingMagToMech");
			
      //			std::cout << "check size of matcoeffnc" << std::endl;
			UInt numRows, numCols;
			mechToMag->GetTensorSize(numRows,numCols);
      //			std::cout << "mechToElec: " << numRows << " " << numCols << std::endl;
			
			magToMech->GetTensorSize(numRows,numCols);
      //			std::cout << "elecToMech: " << numRows << " " << numCols << std::endl;

      // NOTE: the used set of governing equations
      //   H = -[h]*S + [nu]*B
      //   [sigma] = [c]*S - [h]*B
      // couple with a negative sign, so we have to use a factor of -1 for the coupling matrices
      // > due to convention the MatCoefFnc will return +h, so we have to multiply by -1 as for
      //   the linear case
      double factor = -1.0;

			if ( isComplex ) {
				if( subType_ == "axi" ) {
					*magToMechInt = new ADBInt<Complex>(new StrainOperatorAxi<FeH1,Complex>(),
                  new CurlOperatorAxi<Complex>(),
                  magToMech, factor, true );
					*mechToMagInt = new ADBInt<Complex>(new CurlOperatorAxi<Complex>(),
                  new StrainOperatorAxi<FeH1,Complex>(),
                  mechToMag, factor, true );
				} else if( subType_ == "planeStrain" ) {
					*magToMechInt = new ADBInt<Complex>(new StrainOperator2D<FeH1,Complex>(),
                  new CurlOperator<FeH1,2, Complex>(),
                  magToMech, factor, true);
					*mechToMagInt = new ADBInt<Complex>(new CurlOperator<FeH1,2, Complex>(),
                  new StrainOperator2D<FeH1,Complex>(),
                  mechToMag, factor, true);
				} else if( subType_ == "planeStress" ) {
					*magToMechInt = new ADBInt<Complex>(new StrainOperator2D<FeH1,Complex>(),
                  new CurlOperator<FeH1,2, Complex>(),
                  magToMech, factor, true);
					*mechToMagInt = new ADBInt<Complex>(new CurlOperator<FeH1,2, Complex>(),
                  new StrainOperator2D<FeH1,Complex>(),
                  mechToMag, factor, true);
				} else if( subType_ == "3d") {
					*magToMechInt = new ADBInt<Complex>(new StrainOperator3D<FeH1,Complex>(),
                  new CurlOperator<FeH1,3, Complex>(),
                  magToMech, factor, true);
					*mechToMagInt = new ADBInt<Complex>(new CurlOperator<FeH1,3, Complex>(),
                  new StrainOperator3D<FeH1,Complex>(),
                  mechToMag, factor, true);
				} else if( subType_ == "2.5d") {
					EXCEPTION( "Subtype '" << subType_ << "' not implemented for magstrict yet" );
				} else {
					EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
				}
			}
			else {
				if( subType_ == "axi" ) {
					*magToMechInt = new ADBInt<Double>(new StrainOperatorAxi<FeH1>(),
                  new CurlOperatorAxi<Double>(),
                  magToMech, 1.0, true );
					*mechToMagInt = new ADBInt<Double>(new CurlOperatorAxi<Double>(),
                  new StrainOperatorAxi<FeH1>(),
                  mechToMag, 1.0, true );
				} else if( subType_ == "planeStrain" ) {
					*magToMechInt = new ADBInt<Double>(new StrainOperator2D<FeH1>(),
                  new CurlOperator<FeH1,2, Double>(),
                  magToMech, 1.0, true);
					*mechToMagInt = new ADBInt<Double>(new CurlOperator<FeH1,2, Double>(),
                  new StrainOperator2D<FeH1>(),
                  mechToMag, 1.0, true);
				} else if( subType_ == "planeStress" ) {
					*magToMechInt = new ADBInt<Double>(new StrainOperator2D<FeH1>(),
                  new CurlOperator<FeH1,2, Double>(),
                  magToMech, 1.0, true);
					*mechToMagInt = new ADBInt<Double>(new CurlOperator<FeH1,2, Double>(), 
                  new StrainOperator2D<FeH1>(),
                  mechToMag, 1.0, true);
				} else if( subType_ == "3d") {
					*magToMechInt = new ADBInt<Double>(new CurlOperator<FeH1,3, Double>(),
                  new GradientOperator<FeH1,3>(),
                  magToMech, 1.0, true);
          *mechToMagInt = new ADBInt<Double>(new GradientOperator<FeH1,3>(),
                  new CurlOperator<FeH1,3, Double>(),
                  mechToMag, 1.0, true);
				} else if( subType_ == "2.5d") {
					EXCEPTION( "Subtype '" << subType_ << "' not implemented for magstrict yet" );
				} else {
					EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
				}
			}

      PtrCoefFct irrStrains = hystOperator->GenerateOutputCoefFnc("IrrStrainsMagstrict_TensorForm");
      irrStrains_->AddRegion( regionId, irrStrains);
      PtrCoefFct irrStressesVec = hystOperator->GenerateOutputCoefFnc("IrrStressesMagstrict_VectorForm");
      irrStressesVector_->AddRegion( regionId, irrStressesVec);
      PtrCoefFct irrStresses = hystOperator->GenerateOutputCoefFnc("IrrStressesMagstrict_TensorForm");
      irrStresses_->AddRegion( regionId, irrStresses);
      PtrCoefFct coupledIrrStrains = hystOperator->GenerateOutputCoefFnc("CoupledIrrStrainsMagstrict");
      coupledIrrStrains_->AddRegion( regionId, coupledIrrStrains);

			return true;
		} else {
			// NON-HYST CASE 
			*magToMechInt = GetStiffIntegrator( actSDMat, regionId, isComplex );
			return false;
		}
  }

  BaseBDBInt *
  MagnetoStrictCoupling::GetStiffIntegrator( BaseMaterial* actSDMat,
          RegionIdType regionId,
          bool isComplex ) {

    // Get region name
    std::string regionName = ptGrid_->GetRegion().ToString( regionId );
		
		// makes no difference at the moment
    bool coordUpdate_ = true;
		
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
                curCoef, factor, coordUpdate_ );
      } else if( subType_ == "planeStrain" ) {
				
        integ = new ADBInt<Complex>(new StrainOperator2D<FeH1,Complex>(),
                new CurlOperator<FeH1,2, Complex>(),
                //new GradientOperator<FeH1,2,1,Complex>(),
                curCoef, factor, coordUpdate_);
      } else if( subType_ == "planeStress" ) {
				
        integ = new ADBInt<Complex>(new StrainOperator2D<FeH1,Complex>(),
                new CurlOperator<FeH1,2, Complex>(),
                //new GradientOperator<FeH1,2,1,Complex>(),
                curCoef, factor, coordUpdate_);
      } else if( subType_ == "3d") {
        integ = new ADBInt<Complex>(new StrainOperator3D<FeH1,Complex>(),
                new CurlOperator<FeH1,3, Complex>(),
                //new GradientOperator<FeH1,3,1,Complex>(),
                curCoef, factor, coordUpdate_);
      } else {
        EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
      }
    }
    else {
      if( subType_ == "axi" ) {
        integ = new ADBInt<Double>(new StrainOperatorAxi<FeH1>(),
                new CurlOperatorAxi<Double>(),
                //new GradientOperator<FeH1,2>(),
                curCoef, factor, coordUpdate_ );
      } else if( subType_ == "planeStrain" ) {
				
        integ = new ADBInt<Double>(new StrainOperator2D<FeH1>(),
                new CurlOperator<FeH1,2, Double>(),
                //new GradientOperator<FeH1,2>(),
                curCoef, factor, coordUpdate_);
      } else if( subType_ == "planeStress" ) {
        integ = new ADBInt<Double>(new StrainOperator2D<FeH1>(),
                new CurlOperator<FeH1,2, Double>(),
                //new GradientOperator<FeH1,2>(),
                curCoef, factor, coordUpdate_);
      } else if( subType_ == "3d") {
        integ = new ADBInt<Double>(new StrainOperator3D<FeH1>(),
                new CurlOperator<FeH1,3, Double>(),
                //new GradientOperator<FeH1,2>(),
                curCoef, factor, coordUpdate_);
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
			
      TimeSchemeGLM::NonLinType nlType = (nonLin_ || isHyst_)? TimeSchemeGLM::INCREMENTAL : TimeSchemeGLM::NONE;

      shared_ptr<BaseTimeScheme> magScheme(new TimeSchemeGLM(GLMScheme::NEWMARK, 0, nlType) );
			
      magFct->SetTimeScheme(magScheme);
      magFct->GetTimeScheme()->Init(magFct->GetSingleVector(),dt);
    }
  }
	
}

