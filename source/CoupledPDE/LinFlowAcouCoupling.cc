// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "LinFlowAcouCoupling.hh"

#include "PDE/SinglePDE.hh"
#include "PDE/AcousticPDE.hh"
#include "PDE/LinFlowPDE.hh"
#include "CoupledPDE/BasePairCoupling.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Enum.hh"
#include "Materials/BaseMaterial.hh"
#include "Driver/FormsContexts.hh"
#include "Driver/Assemble.hh"
#include "Domain/CoefFunction/CoefFunction.hh"
#include "Domain/CoefFunction/CoefXpr.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/IdentityOperatorNormal.hh"
// include fespaces
#include "FeBasis/H1/H1Elems.hh"

// new integrator concept
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/Operators/IdentityOperator.hh"

#include "Domain/Mesh/NcInterfaces/MortarInterface.hh"

namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  LinFlowAcouCoupling::LinFlowAcouCoupling( SinglePDE *pde1, SinglePDE *pde2,
                                			PtrParamNode paramNode,
											PtrParamNode infoNode,
											shared_ptr<SimState> simState,
											Domain* domain)
    : BasePairCoupling( pde1, pde2, paramNode, infoNode, simState, domain )
  {
    couplingName_ = "linFlowAcouDirect";
    materialClass_ = FLOW;  //we need just material parameters from
                            //the two individual PDEs

    // determine subtype
    pde1_->GetParamNode()->GetValue( "subType", subType_ );

    nonLin_ = false;
    
    // Initialize nonlinearities
    InitNonLin();
  }


  // **************
  //   Destructor
  // **************
  LinFlowAcouCoupling::~LinFlowAcouCoupling() {
  }


  // *********************
  //   DefineIntegrators
  // *********************
  void LinFlowAcouCoupling::DefineIntegrators() {
    // get PDEs
    LinFlowPDE* linflowPDE = dynamic_cast<LinFlowPDE*>(pde1_);
    AcousticPDE* acouPDE = dynamic_cast<AcousticPDE*>(pde2_);
    // get math parser
    MathParser * mp = domain_->GetMathParser();

    // Get balance of momentum sign of LinFlowPDE
    double linFlowBalanceOfMomentumSign = linflowPDE->GetBalanceOfMomentumSign();

    // get hold of all feFunctions (depending on the acoustic formulation)
    shared_ptr<BaseFeFunction> velFct = linflowPDE->GetFeFunction(FLUIDMECH_VELOCITY);
    shared_ptr<BaseFeFunction> presFct = linflowPDE->GetFeFunction(FLUIDMECH_PRESSURE);

    SolutionType acouFormulation = acouPDE->GetFormulation();
    shared_ptr<BaseFeFunction> acouFct = acouPDE->GetFeFunction(acouFormulation);

    shared_ptr<FeSpace> velSpace = velFct->GetFeSpace();
    shared_ptr<FeSpace> presSpace = presFct->GetFeSpace();
    shared_ptr<FeSpace> acouSpace = acouFct->GetFeSpace();

    PtrCoefFct constOne = CoefFunction::Generate( mp, Global::REAL, "1.0");

    //get the volume regions needed by the Surface-Bilinear Forms
    //and create coefficient functions for all acoustic densities
    std::map<RegionIdType, BaseMaterial*> acouMaterials;
    acouMaterials = acouPDE->GetMaterialData();
    std::set< RegionIdType > acouRegions;
    std::map< RegionIdType, PtrCoefFct > coefFuncs;
    std::map<RegionIdType, BaseMaterial*>::iterator it, end;
    it = acouMaterials.begin();
    end = acouMaterials.end();
    for( ; it != end; it++ ) {
      RegionIdType volRegId = it->first;
      acouRegions.insert(volRegId);
      BaseMaterial * acouMat = acouMaterials[volRegId];
      coefFuncs[volRegId] = acouMat->GetScalCoefFnc(DENSITY,Global::REAL);
    }

    // check which formulation we are using
    switch(acouFormulation) {
    // ==================
    //   ACOU_POTENTIAL
    // ==================
    case ACOU_POTENTIAL:

      //loop over all conforming surface regions
      for ( UInt actSD = 0, n = entityLists_.GetSize(); actSD < n; actSD++ ) {
        shared_ptr<SurfElemList> actSDList =
            dynamic_pointer_cast<SurfElemList>(entityLists_[actSD]);
        //RegionIdType region = actSDList->GetRegion();

        velFct->AddEntityList(actSDList); //former flowFct
        acouFct->AddEntityList(actSDList);


        //
        // acoustic potential couples to linFlow velocity
        //
        BiLinearForm *cplInt1 = NULL;
        if( dim_ == 2  ) {
            cplInt1 = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
                                         new IdentityOperatorNormal<FeH1,2>(),
                                         coefFuncs, -linFlowBalanceOfMomentumSign, acouRegions);
        }
        else  {
            cplInt1 = new SurfaceABInt<>(new IdentityOperator<FeH1,3,3>(),
                                         new IdentityOperatorNormal<FeH1,3>(),
                                         coefFuncs, -linFlowBalanceOfMomentumSign, acouRegions);
        }
        cplInt1->SetName("LinFlowAcouCouplingInt1");
        BiLinFormContext *context1 = new BiLinFormContext(cplInt1, DAMPING);
        context1->SetEntities( actSDList, actSDList );
        context1->SetFeFunctions( velFct, acouFct ); //former flowFct
        assemble_->AddBiLinearForm( context1 );

        //
        // linFlow velocity couples to acoustic potential
        // (minus one due to surface normals!)
        //
        BiLinearForm *cplInt2 = NULL;
        if( dim_ == 2  ) {
            cplInt2 = new SurfaceABInt<>(new IdentityOperatorNormal<FeH1,2>(),
                                        new IdentityOperator<FeH1,2,2>(),
                                        constOne, 1.0, acouRegions);
        }
        else  {
            cplInt2 = new SurfaceABInt<>(new IdentityOperatorNormal<FeH1,3>(),
                                        new IdentityOperator<FeH1,3,3>(),
                                        constOne, 1.0, acouRegions);
        }
        cplInt2->SetName("AcouLinFlowCouplingInt2");
        BiLinFormContext *context2 = new BiLinFormContext(cplInt2, STIFFNESS);
        context2->SetEntities( actSDList, actSDList );
        context2->SetFeFunctions( acouFct, velFct ); //former flowFct
        assemble_->AddBiLinearForm( context2 );
      }

      //loop over all non-conforming interfaces
      for ( UInt actNC = 0, n = ncInterfaceIds_.GetSize(); actNC < n; actNC++ ) {

        // Get interface from grid and cast to MortarInterface class
        shared_ptr<BaseNcInterface> ncIf = ptGrid_->GetNcInterface(ncInterfaces_[actNC].interfaceId);

        MortarInterface *mortarIf = dynamic_cast<MortarInterface*>(ncIf.get());
        assert(mortarIf);

        // create new entity list
        shared_ptr<ElemList> actSDList = mortarIf->GetElemList();


        //
        // acoustic potential couples to linFlow velocity
        //
        BiLinearForm * ncCplInt1 = NULL;
        if( dim_ == 2  ) {
          ncCplInt1 = new SurfaceMortarABIntMA<>( new IdentityOperator<FeH1,2,2>(),
                                                new IdentityOperatorNormal<FeH1,2>(),
                                                coefFuncs, -linFlowBalanceOfMomentumSign, mortarIf->IsCoplanar(),
                                                geoUpdate_);
        }
        else {
          ncCplInt1 = new SurfaceMortarABIntMA<>( new IdentityOperator<FeH1,3,3>(),
                                                new IdentityOperatorNormal<FeH1,3>(),
                                                coefFuncs, -linFlowBalanceOfMomentumSign, mortarIf->IsCoplanar(),
                                                geoUpdate_);
        }
        ncCplInt1->SetName("LinFlowAcouCouplingNCInt1");
        NcBiLinFormContext *ncContext1 = NULL;
        ncContext1 = new NcBiLinFormContext(ncCplInt1, DAMPING);
        ncContext1->SetEntities( actSDList, actSDList );
        ncContext1->SetFeFunctions( velFct, acouFct ); //former flowFct
        assemble_->AddBiLinearForm( ncContext1 );


        //
        // linFlow velocity couples to acoustic potential
        // (minus one due to surface normals!)
        //
        BiLinearForm *ncCplInt2 = NULL;
        if( dim_ == 2  ) {
            ncCplInt2 = new SurfaceMortarABIntMA<>(new IdentityOperatorNormal<FeH1,2>(),
                                                 new IdentityOperator<FeH1,2,2>(),
                                                 constOne, 1.0, mortarIf->IsCoplanar(),
                                                 geoUpdate_);
        }
        else  {
            ncCplInt2 = new SurfaceMortarABIntMA<>(new IdentityOperatorNormal<FeH1,3>(),
                                                 new IdentityOperator<FeH1,3,3>(),
                                                 constOne, 1.0, mortarIf->IsCoplanar(),
                                                 geoUpdate_);
        }
        ncCplInt2->SetName("AcouLinFlowCouplingNCInt2");
        NcBiLinFormContext *ncContext2 = NULL;
        ncContext2 = new NcBiLinFormContext(ncCplInt2, STIFFNESS);
        ncContext2->SetEntities( actSDList, actSDList );
        ncContext2->SetFeFunctions( acouFct, velFct ); //former flowFct
        assemble_->AddBiLinearForm( ncContext2 );
      }

      break;

    // ==================
    //   ACOU_PRESSURE
    // ==================
    case ACOU_PRESSURE:
      //loop over all conforming surface regions
      for ( UInt actSD = 0, n = entityLists_.GetSize(); actSD < n; actSD++ ) {
        shared_ptr<SurfElemList> actSDList =
            dynamic_pointer_cast<SurfElemList>(entityLists_[actSD]);
        //RegionIdType region = actSDList->GetRegion();

        velFct->AddEntityList(actSDList);
        acouFct->AddEntityList(actSDList);

        //
        // acoustic pressure couples to linFlow velocity
        //
        BiLinearForm *cplInt1 = NULL;
        if( dim_ == 2  ) {
            cplInt1 = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
                                         new IdentityOperatorNormal<FeH1,2>(),
                                         constOne, -linFlowBalanceOfMomentumSign, acouRegions);
        }
        else  {
            cplInt1 = new SurfaceABInt<>(new IdentityOperator<FeH1,3,3>(),
                                         new IdentityOperatorNormal<FeH1,3>(),
                                         constOne, -linFlowBalanceOfMomentumSign, acouRegions);
        }
        cplInt1->SetName("LinFlowAcouCouplingInt");
        BiLinFormContext *context1 = new BiLinFormContext(cplInt1, STIFFNESS);
        context1->SetEntities( actSDList, actSDList );
        context1->SetFeFunctions( velFct, acouFct );
        assemble_->AddBiLinearForm( context1 );

        //
        // linFlow velocity couples to acoustic pressure
        // (minus one due to surface normals!)
        //
        BiLinearForm *cplInt2 = NULL;
        if( dim_ == 2  ) {
            cplInt2 = new SurfaceABInt<>(new IdentityOperatorNormal<FeH1,2>(),
                                        new IdentityOperator<FeH1,2,2>(),
                                        coefFuncs, 1.0, acouRegions);
        }
        else  {
            cplInt2 = new SurfaceABInt<>(new IdentityOperatorNormal<FeH1,3>(),
                                        new IdentityOperator<FeH1,3,3>(),
                                        coefFuncs, 1.0, acouRegions);
        }
        cplInt2->SetName("AcouLinFlowCouplingInt");
        BiLinFormContext *context2 = new BiLinFormContext(cplInt2, DAMPING);
        context2->SetEntities( actSDList, actSDList );
        context2->SetFeFunctions( acouFct, velFct );
        assemble_->AddBiLinearForm( context2 );
      }

      //loop over all non-conforming interfaces
      for ( UInt actNC = 0, n = ncInterfaceIds_.GetSize(); actNC < n; actNC++ ) {

        // Get interface from grid and cast to MortarInterface class
        shared_ptr<BaseNcInterface> ncIf = ptGrid_->GetNcInterface(ncInterfaces_[actNC].interfaceId);

        MortarInterface *mortarIf = dynamic_cast<MortarInterface*>(ncIf.get());
        assert(mortarIf);

        // create new entity list
        shared_ptr<ElemList> actSDList = mortarIf->GetElemList();

        //
        // acoustic pressure couples to linFlow velocity
        //
        BiLinearForm * cplInt1 = NULL;
        if( dim_ == 2  ) {
          cplInt1 = new SurfaceMortarABIntMA<>( new IdentityOperator<FeH1,2,2>(),
                                                new IdentityOperatorNormal<FeH1,2>(),
                                                constOne, -linFlowBalanceOfMomentumSign, mortarIf->IsCoplanar(),
                                                geoUpdate_);
        }
        else {
          cplInt1 = new SurfaceMortarABIntMA<>( new IdentityOperator<FeH1,3,3>(),
                                                new IdentityOperatorNormal<FeH1,3>(),
                                                constOne, -linFlowBalanceOfMomentumSign, mortarIf->IsCoplanar(),
                                                geoUpdate_);
        }
        cplInt1->SetName("LinFlowAcouCouplingNCInt");
        NcBiLinFormContext *ncContext1 = NULL;
        ncContext1 = new NcBiLinFormContext(cplInt1, STIFFNESS);
        ncContext1->SetEntities( actSDList, actSDList );
        ncContext1->SetFeFunctions( velFct, acouFct );
        assemble_->AddBiLinearForm( ncContext1 );

        //
        // linFlow velocity couples to acoustic pressure
        // (minus one due to surface normals!)
        //
        BiLinearForm *cplInt2 = NULL;
        if( dim_ == 2  ) {
            cplInt2 = new SurfaceMortarABIntMA<>(new IdentityOperatorNormal<FeH1,2>(),
                                                 new IdentityOperator<FeH1,2,2>(),
                                                 coefFuncs, 1.0, mortarIf->IsCoplanar(),
                                                 geoUpdate_);
        }
        else  {
            cplInt2 = new SurfaceMortarABIntMA<>(new IdentityOperatorNormal<FeH1,3>(),
                                                 new IdentityOperator<FeH1,3,3>(),
                                                 coefFuncs, 1.0, mortarIf->IsCoplanar(),
                                                 geoUpdate_);
        }
        cplInt2->SetName("AcouLinFlowCouplingIntNC");
        NcBiLinFormContext *ncContext2 = NULL;
        ncContext2 = new NcBiLinFormContext(cplInt2, DAMPING);
        ncContext2->SetEntities( actSDList, actSDList );
        ncContext2->SetFeFunctions( acouFct, velFct );
        assemble_->AddBiLinearForm( ncContext2 );
      }

      break;

    default:
      EXCEPTION("Unknown formulation for acoustics.");
      break;
    }
  }

}
