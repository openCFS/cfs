// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "WaterWaveMechCoupling.hh"

#include "PDE/SinglePDE.hh"
#include "PDE/WaterWavePDE.hh"
#include "PDE/MechPDE.hh"
#include "CoupledPDE/BasePairCoupling.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Enum.hh"
#include "Materials/BaseMaterial.hh"
#include "Driver/FormsContexts.hh"
#include "Driver/Assemble.hh"

// include fespaces
#include "FeBasis/H1/H1Elems.hh"

// new integrator concept
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/IdentityOperatorNormal.hh"
#include "Forms/Operators/IdentityOperatorInVector.hh"


namespace CoupledField {



  // ***************
  //   Constructor
  // ***************
  WaterWaveMechCoupling::WaterWaveMechCoupling( SinglePDE *pde1, SinglePDE *pde2,
                                PtrParamNode paramNode, 
                                PtrParamNode infoNode,
                                shared_ptr<SimState> simState, Domain* domain)
    : BasePairCoupling( pde1, pde2, paramNode, infoNode, simState, domain )
  {
    couplingName_ = "waterWaveMechDirect";
    materialClass_ = ACOUSTIC;

    // determine subtype from mechanic pde
    pde1_->GetParamNode()->GetValue( "subType", subType_ );

    nonLin_ = false;
    
    // Initialize nonlinearities
    InitNonLin();
  }


  // **************
  //   Destructor
  // **************
  WaterWaveMechCoupling::~WaterWaveMechCoupling() {
  }


  // *********************
  //   DefineIntegrators
  // *********************
  void WaterWaveMechCoupling::DefineIntegrators() {

    // get hold of both feFunctions
    MechPDE* mechPDE = dynamic_cast<MechPDE*>(pde1_);
    WaterWavePDE* waterPDE = dynamic_cast<WaterWavePDE*>(pde2_);
    MathParser * mp = domain_->GetMathParser();
    
    shared_ptr<BaseFeFunction> dispFct = mechPDE->GetFeFunction(MECH_DISPLACEMENT);
    shared_ptr<BaseFeFunction> waterFct = waterPDE->GetFeFunction(WATER_PRESSURE);
    std::map<RegionIdType, BaseMaterial*> waterMaterials;
    waterMaterials = waterPDE->GetMaterialData();

    // Create coefficient functions for all acoustic densities
    std::map< RegionIdType, PtrCoefFct > coefFuncs;
    std::map< RegionIdType, PtrCoefFct > oneCoefFuncs;
    std::set< RegionIdType > waterRegions;
    std::map<RegionIdType, BaseMaterial*>::iterator it, end;
    it = waterMaterials.begin();
    end = waterMaterials.end();
    for( ; it != end; it++ ) {
      RegionIdType volRegId = it->first;
      waterRegions.insert(volRegId);

      // Get bulk density for water waves
      BaseMaterial * waterMat = waterMaterials[volRegId];
      coefFuncs[volRegId] = waterMat->GetScalCoefFnc(DENSITY,Global::REAL);
      oneCoefFuncs[volRegId] = CoefFunction::Generate(mp, Global::REAL, "1.0");
    }

    shared_ptr<FeSpace> dispSpace = dispFct->GetFeSpace();
    shared_ptr<FeSpace> waterSpace = waterFct->GetFeSpace();
    
    for ( UInt actSD = 0, n = entityLists_.GetSize(); actSD < n; actSD++ ) {

      shared_ptr<SurfElemList> actSDList = 
          dynamic_pointer_cast<SurfElemList>(entityLists_[actSD]);

      waterFct->AddEntityList(actSDList);
      dispFct->AddEntityList(actSDList);

      // First ensure, that we have no coupled eigenfrequency simulation:
     DefCouplInt( "WaterWaveMechPresStiffCouplingInt", false, -1.0, STIFFNESS, dispFct,
          waterFct, actSDList, oneCoefFuncs, waterRegions );

      DefCouplInt( "WaterWaveMechPresMassCouplingInt", false, 1.0, MASS, waterFct,
          dispFct, actSDList, coefFuncs, waterRegions );

      //add gravity term;
      Vector<Double> gravity(dim_);
      gravity.Init(); 
      BiLinearForm * cplInt = NULL;
      if ( dim_ == 2 ) {
    	  gravity[1] = -9.81;
    	  cplInt = new SurfaceABInt<>(new IdentityOperatorInNormal<FeH1,2>(),
    	                              new IdentityOperatorInVector<FeH1,2>(gravity),
    	                              coefFuncs, -1.0, waterRegions);
      }
      else if (dim_ == 3) {
        gravity[2] = -9.81;
    	  cplInt = new SurfaceABInt<>(new IdentityOperatorInNormal<FeH1,3>(),
    	                              new IdentityOperatorInVector<FeH1,3>(gravity),
    	                              coefFuncs, -1.0, waterRegions);    
      }
      else {
        EXCEPTION( "Coupling only for two and three dimensions defined" );
      }
      cplInt->SetName("WaterWaveMechGravityInt");
      BiLinFormContext * context =
    		  new BiLinFormContext(cplInt, STIFFNESS );

      context->SetEntities( actSDList, actSDList );
      context->SetFeFunctions( dispFct, dispFct );
      assemble_->AddBiLinearForm( context );
    }
  }

  void WaterWaveMechCoupling::DefCouplInt( const std::string& name,
                                      bool assembleTransposed,
                                      Double factor,
                                      FEMatrixType matType,
                                      shared_ptr<BaseFeFunction>& fnc1,
                                      shared_ptr<BaseFeFunction>& fnc2,
                                      shared_ptr<SurfElemList>& actSDList,
                                      const std::map< RegionIdType, PtrCoefFct >& coefFuncs,
                                      const std::set< RegionIdType >& waterRegions ) {
    
    // check for position of integrator
    SolutionType rowType = fnc1->GetResultInfo()->resultType;
    BiLinearForm * cplInt = NULL;
    if( dim_ == 2  ) {
      if(rowType == MECH_DISPLACEMENT) {
        cplInt = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
                                    new IdentityOperatorNormal<FeH1,2>(),
                                    coefFuncs, factor, waterRegions);
      } else {
        cplInt = new SurfaceABInt<>(new IdentityOperatorNormal<FeH1,2>(),
                                    new IdentityOperator<FeH1,2,2>(),
                                    coefFuncs, factor, waterRegions);
      }
    } else if( dim_ == 3) {
      if(rowType == MECH_DISPLACEMENT) {
        cplInt = new SurfaceABInt<>(new IdentityOperator<FeH1,3,3>(),
                                    new IdentityOperatorNormal<FeH1,3>(),
                                    coefFuncs, factor, waterRegions);
      } else {
        cplInt = new SurfaceABInt<>(new IdentityOperatorNormal<FeH1,3>(),
                                    new IdentityOperator<FeH1,3,3>(),
                                    coefFuncs, factor, waterRegions);
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
  
}
