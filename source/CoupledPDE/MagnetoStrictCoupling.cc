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
		
		nonLin_ = pde1_->IsNonLin() || pde2_->IsNonLin();
		
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
      
      // check if mag pde is hysteretic
			bool isHyst = false;
			if(pde1_->IsHysteresis()){
				EXCEPTION("Currently only the mag PDE may be hysteretic");
			} else if(pde2_->IsHysteresis()){
				isHyst = true;
			}
			
			// copied from PiezoCoupling.cc
			if(isHyst){
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
			
			if(considerCounterpart_){
				// treat magMech and mechMag separately
				// bdbInts_ stores bdbIntegrator for mechanical pde (magToMechInt)
				BaseBDBInt* bdb = it->second;
				stressCplFunc->AddIntegrator(bdb, region);
//				std::cout << "StressCpl success" << std::endl;
				
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

		// check for hysteresis
		bool isHyst = false;
		if(pde1_->IsHysteresis()){
			EXCEPTION("Currently only the elec PDE may be hysteretic");
		} else if(pde2_->IsHysteresis()){
			isHyst = true;
		}

		if(isHyst){
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
				
				
				// mech RHS
				shared_ptr<CoefFunction> mechRHS = regionHystOperator->GenerateRHSCoefFnc("MagStrictLoadForMechPDE",stiffCoef,couplingCoef);
				
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

        lin->SetName("MagStrictLoadForMechPDE");
        lin->SetSolDependent();
        LinearFormContext *ctx = new LinearFormContext( lin );
        ctx->SetEntities( actSDList );
        ctx->SetFeFunction(dispFct);
        
				assemble_->AddLinearForm(ctx);
        // Add entity list will add nothing, if entities were already assigned
        dispFct->AddEntityList(actSDList);
				
				// elec RHS
				shared_ptr<CoefFunction> magRHS = regionHystOperator->GenerateRHSCoefFnc("MagStrictLoadForMagPDE",stiffCoef,couplingCoef);
				
				if(isComplex_) {
          if( dim_ == 2 ) {
            if(isaxi_){
           // we need +factor as we put -rotM to the rhs
           // note: use BDUIntegrator, even though we have no material (dMat = Identity)
              //lin = new BDUIntegrator<CurlOperatorAxi<Double>, Complex>(Complex(factor),
                    //                                                         coef[i], reluc_, coefUpdateGeo);

              lin = new BUIntegrator<Complex>( new CurlOperatorAxi<Double>(),
                                           Complex(factor),magRHS,  coefUpdateGeo, fullevaluation);
            } else {
              lin = new BUIntegrator<Complex>( new CurlOperator<FeH1,2,Double>(),
                                              Complex(factor),magRHS,  coefUpdateGeo, fullevaluation);
            }
          } else {
            lin = new BUIntegrator<Complex>( new CurlOperator<FeH1,3,Double>(),
                                            Complex(factor),magRHS,  coefUpdateGeo, fullevaluation);
          }
        } else  {
          if( dim_ == 2 ) {
            if(isaxi_){
           // we need +factor as we put -rotH to the rhs
              lin = new BUIntegrator<Double>( new CurlOperatorAxi<Double>(),
                                           factor,magRHS,  coefUpdateGeo, fullevaluation);
            } else {
              lin = new BUIntegrator<Double>( new CurlOperator<FeH1,2,Double>(),
                                              factor,magRHS,  coefUpdateGeo, fullevaluation);
            }
          } else {
            lin = new BUIntegrator<Double>( new CurlOperator<FeH1,3,Double>(),
                                            factor,magRHS,  coefUpdateGeo, fullevaluation);
          }
        }

        linMag->SetName("MagStrictLoadForMagPDE");
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
    if( subType_ == "planeStrain" || subType_ == "planeStress" ) {
      tensorType = PLANE_STRAIN;
    } else if( subType_ == "axi" ){
      tensorType = AXI;
    } else if (subType_ == "3d" || subType_ == "2.5d"){
      tensorType = FULL;
    } else {
      EXCEPTION( "Unknown subtype '" << subType_ << "'" );
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
		StdVector<NonLinType> nonLinTypes = nonLinMap[regionId];
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
				
			// create hyst material
			// note: curCoef holds the coupling tensor
			// hyst case might be unsymmetric
			// > in case of deltaFormulation, elecToMech holds an addition deltaS/deltaE
			PtrCoefFct mechToMag = hystOperator->GenerateMatCoefFnc("CouplingMechToMag",stiffCoef,couplingCoef);
			PtrCoefFct magToMech = hystOperator->GenerateMatCoefFnc("CouplingMagToMech",stiffCoef,couplingCoef);
			
//			std::cout << "check size of matcoeffnc" << std::endl;
			UInt numRows, numCols;
			mechToMag->GetTensorSize(numRows,numCols);
//			std::cout << "mechToElec: " << numRows << " " << numCols << std::endl;
			
			magToMech->GetTensorSize(numRows,numCols);
//			std::cout << "elecToMech: " << numRows << " " << numCols << std::endl;
						
			if ( isComplex ) {
				if( subType_ == "axi" ) {
					*magToMechInt = new ADBInt<Complex>(new StrainOperatorAxi<FeH1,Complex>(),
														new CurlOperatorAxi<Complex>(),
														magToMech, 1.0, true );
					*mechToMagInt = new ADBInt<Complex>(new CurlOperatorAxi<Complex>(),
														new StrainOperatorAxi<FeH1,Complex>(),
														mechToMag, 1.0, true );
				} else if( subType_ == "planeStrain" ) {
					*magToMechInt = new ADBInt<Complex>(new StrainOperator2D<FeH1,Complex>(),
														new CurlOperator<FeH1,2, Complex>(),
														magToMech, 1.0, true);
					*mechToMagInt = new ADBInt<Complex>(new CurlOperator<FeH1,2, Complex>(),
														new StrainOperator2D<FeH1,Complex>(),
														mechToMag, 1.0, true);
				} else if( subType_ == "planeStress" ) {
					*magToMechInt = new ADBInt<Complex>(new StrainOperator2D<FeH1,Complex>(),
														new CurlOperator<FeH1,2, Complex>(),
														magToMech, 1.0, true);
					*mechToMagInt = new ADBInt<Complex>(new CurlOperator<FeH1,2, Complex>(),
														new StrainOperator2D<FeH1,Complex>(),
														mechToMag, 1.0, true);
				} else if( subType_ == "3d") {
					*magToMechInt = new ADBInt<Complex>(new StrainOperator3D<FeH1,Complex>(),
														new CurlOperator<FeH1,3, Complex>(),
														magToMech, 1.0, true);
					*mechToMagInt = new ADBInt<Complex>(new CurlOperator<FeH1,3, Complex>(),
														new StrainOperator3D<FeH1,Complex>(),
														mechToMag, 1.0, true);
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
			
      shared_ptr<BaseTimeScheme> magScheme(new TimeSchemeGLM(GLMScheme::NEWMARK, 0) );
			
      magFct->SetTimeScheme(magScheme);
      magFct->GetTimeScheme()->Init(magFct->GetSingleVector(),dt);
    }
  }
	
}

