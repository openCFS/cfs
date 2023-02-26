// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "AcouMechCoupling.hh"

#include "PDE/SinglePDE.hh"
#include "PDE/AcousticPDE.hh"
#include "PDE/MechPDE.hh"
#include "CoupledPDE/BasePairCoupling.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Enum.hh"
#include "Materials/BaseMaterial.hh"
#include "Driver/FormsContexts.hh"
#include "Driver/Assemble.hh"
#include "DataInOut/Logging/LogConfigurator.hh"


// include fespaces
#include "FeBasis/H1/H1Elems.hh"

// new integrator concept
#include "Domain/Mesh/NcInterfaces/MortarInterface.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/IdentityOperatorNormal.hh"
#include "Forms/Operators/SurfaceOperators.hh"
#include "Forms/Operators/DivOperator.hh"

#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/CoefFunction/CoefFunctionCompound.hh"
#include "Domain/CoefFunction/CoefFunctionConst.hh"
#include "Domain/CoefFunction/CoefFunctionMulti.hh"
#include "Domain/CoefFunction/CoefFunctionPML.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/CoefFunction/CoefFunctionImpedanceModel.hh"
#include "Domain/Mesh/NcInterfaces/BaseNcInterface.hh"


#include <boost/lexical_cast.hpp>
#include <cmath>


#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"
#include "Materials/AcousticMaterial.hh"

namespace CoupledField {



  // ***************
  //   Constructor
  // ***************
  AcouMechCoupling::AcouMechCoupling( SinglePDE *pde1, SinglePDE *pde2,
                                PtrParamNode paramNode, 
                                PtrParamNode infoNode,
                                shared_ptr<SimState> simState, Domain* domain)
    : BasePairCoupling( pde1, pde2, paramNode, infoNode, simState, domain )
  {
    couplingName_ = "acouMechDirect";
    materialClass_ = ACOUSTIC;

    // determine subtype from mechanic pde
    pde1_->GetParamNode()->GetValue( "subType", subType_ );

    nonLin_ = false;
    
    // Initialize nonlinearities
    InitNonLin();

    //check consistency in time stepping
    MechPDE* mechPDE = dynamic_cast<MechPDE*>(pde1_);
    AcousticPDE* acouPDE = dynamic_cast<AcousticPDE*>(pde2_);
    PtrParamNode mechParam = mechPDE->GetMyParam();
    PtrParamNode acouParam = acouPDE->GetMyParam();
    if ( mechParam->Get("timeStepAlpha")->As<Double>() != acouParam->Get("timeStepAlpha")->As<Double>() )
    	EXCEPTION("Alpha value of time stepping algorithm has to be the same in acousticPDE and mechPDE");
  }


  // **************
  //   Destructor
  // **************
  AcouMechCoupling::~AcouMechCoupling() {
  }


  // *********************
  //   DefineIntegrators
  // *********************
  void AcouMechCoupling::DefineIntegrators() {

    // get hold of both feFunctions
    MechPDE* mechPDE = dynamic_cast<MechPDE*>(pde1_);
    AcousticPDE* acouPDE = dynamic_cast<AcousticPDE*>(pde2_);
    MathParser * mp = domain_->GetMathParser();
    
    shared_ptr<BaseFeFunction> dispFct = mechPDE->GetFeFunction(MECH_DISPLACEMENT);

    SolutionType acouFormulation = acouPDE->GetFormulation();
    shared_ptr<BaseFeFunction> acouFct = acouPDE->GetFeFunction(acouFormulation);

    std::map<RegionIdType, BaseMaterial*> acouMaterials;
    acouMaterials = acouPDE->GetMaterialData();

    // Create coefficient functions for all acoustic densities
    std::map< RegionIdType, PtrCoefFct > coefFuncs;
    std::map< RegionIdType, PtrCoefFct > oneCoefFuncs;
    std::set< RegionIdType > acouRegions;
    std::map<RegionIdType, BaseMaterial*>::iterator it, end;
    it = acouMaterials.begin();
    end = acouMaterials.end();
    for( ; it != end; it++ ) {
      RegionIdType volRegId = it->first;
      acouRegions.insert(volRegId);

      // Get bulk density for acoustics
      BaseMaterial * acouMat = acouMaterials[volRegId];
      // For a complex fluid formulation, the coupling term is multiplied with 1/density,
      // which cancels out the coefficient of the coupling term (as this coefficient  is the density in
      // non-complex formulation
      if (acouPDE->IsMaterialComplex())
        coefFuncs[volRegId] = CoefFunction::Generate(mp, Global::REAL, lexical_cast<std::string>(1.0));
      else
        coefFuncs[volRegId] = acouMat->GetScalCoefFnc(DENSITY,Global::REAL);

      oneCoefFuncs[volRegId] = CoefFunction::Generate(mp, Global::REAL,
                                                   lexical_cast<std::string>(1.0));
    }

    shared_ptr<FeSpace> dispSpace = dispFct->GetFeSpace();
    shared_ptr<FeSpace> acouSpace = acouFct->GetFeSpace();
    
    for ( UInt actSD = 0, n = entityLists_.GetSize(); actSD < n; actSD++ ) {

      shared_ptr<SurfElemList> actSDList = 
          dynamic_pointer_cast<SurfElemList>(entityLists_[actSD]);

      acouFct->AddEntityList(actSDList);
      dispFct->AddEntityList(actSDList);

#if 0
      // Check which volume region belongs to acoustic PDE
      const SurfElem * surfElem = actSDList->GetSurfElem( 0 );
      RegionIdType volRegId = surfElem->ptVolElems[0]->regionId;
      if(acouMaterials.find(volRegId) == acouMaterials.end()) {
        volRegId = surfElem->ptVolElems[1]->regionId;
      }
#endif


      switch(acouFormulation) {
        // ==================
        //   ACOU_POTENTIAL
        // ==================
        case ACOU_POTENTIAL:
        {
          
          // In case of transient / harmonic simulation, the acoustic PDE is 
          // already pre-multiplied by -1, so we can define one single
          // coupling integrator for the damping matrix, which gets assembled
          // also transposed
          if( analysisType_ != BasePDE::EIGENFREQUENCY) {
            DefCouplInt( "AcouMechPotCouplingInt", true, -1.0, DAMPING, dispFct, 
                         acouFct, actSDList, coefFuncs, acouRegions ); 
          } else {
            // In case of an eigenfrequency simulation, we can only have 
            // positive definite matrices. Thus, the acoustic PDE is NOT
            // pre-multiplied by -1 and we have to define two distinct 
            // coupling integrators for the DAMPING matrix. This is also why
            // factor for the lower diagonal coupling integrator C_Phi_U 
            // has a switched sign.
          
        	DefCouplInt( "AcouMechPotCouplingInt", false, -1.0, DAMPING, dispFct,
                         acouFct, actSDList, coefFuncs, acouRegions );
        	DefCouplInt( "AcouMechPotCouplingInt_Transposed", false, 1.0,
                         DAMPING, acouFct, dispFct, actSDList, coefFuncs, acouRegions );
          }
        }
        break;

        // ==================
        //   ACOU_PRESSURE
        // ==================
        case ACOU_PRESSURE:
          
          DefCouplInt( "AcouMechPresStiffCouplingInt", false, -1.0, STIFFNESS, dispFct, 
                       acouFct, actSDList, oneCoefFuncs, acouRegions );
          
          DefCouplInt( "AcouMechPresMassCouplingInt", false, 1.0, MASS, acouFct, 
                       dispFct, actSDList, coefFuncs, acouRegions );

          DefBLCouplInt( "AcouMechPresStiffBLCouplingInt", false, 1.0, STIFFNESS, dispFct, 
                       acouFct, actSDList, coefFuncs, acouRegions );

          DefBLCouplInt( "AcouMechPresMassBLCouplingInt", false, 1.0, MASS, acouFct, 
                       dispFct, actSDList, coefFuncs, acouRegions );


         
          break;

        default:
          EXCEPTION("Unknown formulation for acoustics.");
          break;
      }
    }

    //now we do it for the non-conforming interfaces
    for ( UInt actNC = 0, n = ncIfaces_.GetSize(); actNC < n; actNC++ ) {

      // Get interface from grid and cast to MortarInterface class
      shared_ptr<BaseNcInterface> ncIf = ptGrid_->GetNcInterface(ncInterfaces_[actNC].interfaceId);

      switch(acouFormulation) {
        // ==================
        //   ACOU_POTENTIAL
        // ==================
        case ACOU_POTENTIAL:
        {

          // In case of transient / harmonic simulation, the acoustic PDE is
          // already pre-multiplied by -1, so we can define one single
          // coupling integrator for the damping matrix, which gets assembled
          // also transposed
          if( analysisType_ != BasePDE::EIGENFREQUENCY) {
            DefCouplIntNC( "AcouMechPotCouplingIntNC", true, -1.0, DAMPING, dispFct,
                         acouFct, ncIf, coefFuncs );
          } else {
            // In case of an eigenfrequency simulation, we can only have
            // positive definite matrices. Thus, the acoustic PDE is NOT
            // pre-multiplied by -1 and we have to define two distinct
            // coupling integrators for the DAMPING matrix. This is also why
            // factor for the lower diagonal coupling integrator C_Phi_U
            // has a switched sign.

        	  DefCouplIntNC( "AcouMechPotCouplingInNC", false, -1.0, DAMPING, dispFct,
                             acouFct, ncIf, coefFuncs );
        	  DefCouplIntNC( "AcouMechPotCouplingInt_TransposedNC", false, 1.0,
                             DAMPING, acouFct, dispFct, ncIf, coefFuncs );
          }
        }
        break;

        // ==================
        //   ACOU_PRESSURE
        // ==================
        case ACOU_PRESSURE:

          DefCouplIntNC( "AcouMechPresStiffCouplingIntNC", false, -1.0, STIFFNESS, dispFct,
                       acouFct, ncIf, oneCoefFuncs );

          DefCouplIntNC( "AcouMechPresMassCouplingIntNC", false, 1.0, MASS, acouFct, dispFct, ncIf, coefFuncs );
          break;

        default:
          EXCEPTION("Unknown formulation for acoustics.");
          break;
      }
    }
  }

  void AcouMechCoupling::DefCouplInt( const std::string& name,
                                      bool assembleTransposed,
                                      Double factor,
                                      FEMatrixType matType,
                                      shared_ptr<BaseFeFunction>& fnc1,
                                      shared_ptr<BaseFeFunction>& fnc2,
                                      shared_ptr<SurfElemList>& actSDList,
                                      const std::map< RegionIdType, PtrCoefFct >& coefFuncs,
                                      const std::set< RegionIdType >& acouRegions ) {
    
    // check for position of integrator
    SolutionType rowType = fnc1->GetResultInfo()->resultType;
    BiLinearForm * cplInt = NULL;
    if( dim_ == 2  ) {
      if(rowType == MECH_DISPLACEMENT) {
        cplInt = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
                                    new IdentityOperatorNormal<FeH1,2>(),
                                    coefFuncs, factor, acouRegions);
      } else {
        cplInt = new SurfaceABInt<>(new IdentityOperatorNormal<FeH1,2>(),
                                    new IdentityOperator<FeH1,2,2>(),
                                    coefFuncs, factor, acouRegions);
      }
    } else if( dim_ == 3) {
      if(rowType == MECH_DISPLACEMENT) {
        cplInt = new SurfaceABInt<>(new IdentityOperator<FeH1,3,3>(),
                                    new IdentityOperatorNormal<FeH1,3>(),
                                    coefFuncs, factor, acouRegions);
      } else {
        cplInt = new SurfaceABInt<>(new IdentityOperatorNormal<FeH1,3>(),
                                    new IdentityOperator<FeH1,3,3>(),
                                    coefFuncs, factor, acouRegions);
      }
    } else {
      EXCEPTION( "Coupling only for two and three dimensions defined" );
    }
    
    cplInt->SetName(name);
    BiLinFormContext * context =
        new BiLinFormContext(cplInt, matType );

    context->SetEntities( actSDList, actSDList );
    context->SetFeFunctions( fnc1, fnc2 );
    context->SetCounterPart(assembleTransposed);
    assemble_->AddBiLinearForm( context );
  }


  // Bilinear terms for Acou-Mech boundary layer coupling

  void AcouMechCoupling::DefBLCouplInt( const std::string& name,
                                      bool assembleTransposed,
                                      Double factor,
                                      FEMatrixType matType,
                                      shared_ptr<BaseFeFunction>& fnc1,
                                      shared_ptr<BaseFeFunction>& fnc2,
                                      shared_ptr<SurfElemList>& actSDList,
                                      const std::map< RegionIdType, PtrCoefFct >& coefFuncs,
                                      const std::set< RegionIdType >& acouRegions ) {


    // check for position of integrator
    SolutionType rowType = fnc1->GetResultInfo()->resultType;
    BiLinearForm * BL_cplInt = NULL;

    
      //========================================================================================
      // boundary Layers
      //========================================================================================

      PtrParamNode bcNode = myParam_->Get( "surfRegionList", ParamNode::PASS );

      ParamNodeList blNodes = bcNode->GetList( "surfRegion" );

      for( UInt i = 0; i < blNodes.GetSize(); ++i ) {
        
        MathParser * mp = domain_->GetMathParser();

        std::string regionName = blNodes[i]->Get("name")->As<std::string>();
        shared_ptr<EntityList> actSDList =  ptGrid_->GetEntityList( EntityList::SURF_ELEM_LIST, regionName );
        std::string volRegName = blNodes[i]->Get("volumeRegion")->As<std::string>();
        RegionIdType volRegion = ptGrid_->GetRegion().Parse(volRegName);

        // PtrCoefFct rho0 = materials_[volRegion]->GetScalCoefFnc( DENSITY, Global::REAL );
        PtrCoefFct rho0 = CoefFunction::Generate( mp,  Global::REAL, blNodes[i]->Get("rho0")->As<std::string>() );
        // PtrCoefFct rho0 = CoefFunction::Generate( mp,  Global::REAL,  "2e-4");
        PtrCoefFct dynamicViscosity = CoefFunction::Generate( mp,  Global::REAL, blNodes[i]->Get("dynamicViscosity")->As<std::string>() );
        PtrCoefFct omegaHalv = CoefFunction::Generate( mp,  Global::REAL, "pi*f");
        PtrCoefFct omega = CoefFunction::Generate( mp,  Global::REAL, "2*pi*f");

        // Kinematic viscosity dynamicViscosity/rho0
        PtrCoefFct nu = CoefFunction::Generate(mp, Global::REAL,CoefXprBinOp(mp, dynamicViscosity, rho0, CoefXpr::OP_DIV ));

        // deltaV = sqrt( 2*nu/omega )
        PtrCoefFct deltaV = CoefFunction::Generate(mp,Global::REAL, CoefXprUnaryOp(mp, CoefXprBinOp(mp, nu, omegaHalv, CoefXpr::OP_DIV ) , CoefXpr::OP_SQRT));
   
        // 1+i ... 
        PtrCoefFct onePlusI =  CoefFunction::Generate( mp,  Global::COMPLEX, "1", "1");

        // 1-i ... 
        PtrCoefFct iMinusOne =  CoefFunction::Generate( mp,  Global::COMPLEX, "-1", "1");


        // Coefficient for the mech-acou stiffness coupling integral
        PtrCoefFct coefK =  CoefFunction::Generate(
         mp,Global::COMPLEX,CoefXprBinOp(mp,CoefXprBinOp(mp, CoefXprBinOp( mp, CoefXprBinOp(mp, dynamicViscosity, iMinusOne, CoefXpr::OP_MULT), deltaV, CoefXpr::OP_DIV ),rho0, CoefXpr::OP_DIV ),omega, CoefXpr::OP_DIV )  );


        PtrCoefFct coefM =  CoefFunction::Generate(mp,Global::COMPLEX,CoefXprBinOp(mp, deltaV, onePlusI, CoefXpr::OP_MULT) );

        std::set<RegionIdType> volRegions;
        volRegions.insert(volRegion);

        if( dim_ == 2  ) {
        
        if(rowType == MECH_DISPLACEMENT) {

          BL_cplInt = new SurfaceABInt<Complex>(new IdentityOperator<FeH1,2,2>(),
                                      new SurfaceTangentialGradientOperator<FeH1,2>(),
                                      coefK, -1.0, volRegions);

        } 
        else {


          BL_cplInt = new SurfaceABInt<Complex>(new IdentityOperator<FeH1,2>(),
                                      new SurfaceTangentialDivergenceOfTangentialVector<FeH1,2,2>(),
                                       coefM, 0.5, volRegions);
        }

      } 
      else if( dim_ == 3) {
        if(rowType == MECH_DISPLACEMENT) {

          BL_cplInt = new SurfaceABInt<Complex>(new IdentityOperator<FeH1,3,3>(),
                                      new SurfaceTangentialGradientOperator<FeH1,3>(),
                                      coefK, -1.0, volRegions);
        } 
        else {

          BL_cplInt = new SurfaceABInt<Complex>(new IdentityOperator<FeH1,3>(),
                                       new SurfaceTangentialDivergenceOfTangentialVector<FeH1,3,3>(),
                                      coefM, 0.5, volRegions);
        }    

      } 
      else {
        EXCEPTION( "Coupling only for two and three dimensions defined" );
      }       
        
      }          


    BL_cplInt->SetName(name);
    
    BiLinFormContext * BL_context =
        new BiLinFormContext(BL_cplInt, matType );


    BL_context->SetEntities( actSDList, actSDList );
    BL_context->SetFeFunctions( fnc1, fnc2 );
    BL_context->SetCounterPart(assembleTransposed);
    assemble_->AddBiLinearForm( BL_context );
  
  }




  
  void AcouMechCoupling::DefCouplIntNC( const std::string& name,
                                      bool assembleTransposed,
                                      Double factor,
                                      FEMatrixType matType,
                                      shared_ptr<BaseFeFunction>& fnc1,
                                      shared_ptr<BaseFeFunction>& fnc2,
									  shared_ptr<BaseNcInterface> ncIf,
									  const std::map< RegionIdType, PtrCoefFct >& coefFuncs) {


	MortarInterface *mortarIf = dynamic_cast<MortarInterface*>(ncIf.get());
	assert(mortarIf);

	// create new entity list
  	shared_ptr<ElemList> actSDList = mortarIf->GetElemList();

  	// check for position of integrator
    SolutionType rowType = fnc1->GetResultInfo()->resultType;
    BiLinearForm * cplInt = NULL;

    NcBiLinFormContext *ncContext = NULL;

    if( dim_ == 2  ) {
      if(rowType == MECH_DISPLACEMENT) {
    	cplInt = new SurfaceMortarABIntMA<>( new IdentityOperator<FeH1,2,2>(),
    			 new IdentityOperatorNormal<FeH1,2>(),
				 coefFuncs, factor, mortarIf->IsPlanar(), geoUpdate_);
      } else {
      	cplInt = new SurfaceMortarABIntMA<>( new IdentityOperatorNormal<FeH1,2>(),
      			 new IdentityOperator<FeH1,2,2>(),
				 coefFuncs, factor, mortarIf->IsPlanar(), geoUpdate_);
      }
    } else if( dim_ == 3) {
      if(rowType == MECH_DISPLACEMENT) {
      	cplInt = new SurfaceMortarABIntMA<>( new IdentityOperator<FeH1,3,3>(),
      			 new IdentityOperatorNormal<FeH1,3>(),
				 coefFuncs, factor, mortarIf->IsPlanar(), geoUpdate_);
      } else {
        cplInt = new SurfaceMortarABIntMA<>( new IdentityOperatorNormal<FeH1,3>(),
        		 new IdentityOperator<FeH1,3,3>(),
				 coefFuncs, factor, mortarIf->IsPlanar(), geoUpdate_);
      }
    } else {
      EXCEPTION( "Coupling only for two and three dimensions defined" );
    }

    cplInt->SetName(name);
    ncContext = new NcBiLinFormContext(cplInt, matType);
    ncContext->SetEntities( actSDList, actSDList );
    ncContext->SetFeFunctions( fnc1, fnc2 );
    ncContext->SetCounterPart(assembleTransposed);
    assemble_->AddBiLinearForm( ncContext );
    ncIf->RegisterIntegrator(ncContext);
  }
}

