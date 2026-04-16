// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "WaterWaveAcousticsCoupling.hh"

#include "PDE/SinglePDE.hh"
#include "PDE/WaterWavePDE.hh"
#include "PDE/AcousticPDE.hh"
#include "CoupledPDE/BasePairCoupling.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Enum.hh"
#include "Materials/BaseMaterial.hh"
#include "Driver/FormsContexts.hh"
#include "Driver/Assemble.hh"

// include fespaces
#include "FeBasis/H1/H1Elems.hh"

// new integrator concept
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/Operators/IdentityOperator.hh"


namespace CoupledField {



  // ***************
  //   Constructor
  // ***************
  WaterWaveAcousticCoupling::WaterWaveAcousticCoupling( SinglePDE *pde1, SinglePDE *pde2,
                                PtrParamNode paramNode, 
                                PtrParamNode infoNode,
                                shared_ptr<SimState> simState,
                                Domain* domain)
    : BasePairCoupling( pde1, pde2, paramNode, infoNode, simState, domain )

  {
    couplingName_ = "WaterWaveAcousticDirect";
    materialClass_ = ACOUSTIC;

    nonLin_ = false;
    
    // Initialize nonlinearities
    InitNonLin();

  }


  // **************
  //   Destructor
  // **************
  WaterWaveAcousticCoupling::~WaterWaveAcousticCoupling() {
  }


  // *********************
  //   DefineIntegrators
  // *********************
  void WaterWaveAcousticCoupling::DefineIntegrators() {

    MathParser * mp = domain_->GetMathParser();

    // get hold of both feFunctions
    WaterWavePDE* waterWavePDE = dynamic_cast<WaterWavePDE*>(pde1_);
    shared_ptr<BaseFeFunction> waterFct = waterWavePDE->GetFeFunction(WATER_PRESSURE);

    AcousticPDE* acousticPDE = dynamic_cast<AcousticPDE*>(pde2_);
    shared_ptr<BaseFeFunction> acouFct = acousticPDE->GetFeFunction(ACOU_PRESSURE);

    shared_ptr<BaseFeFunction> lagrangeMultFct = feFunctions_[LAGRANGE_MULT];

    std::map<RegionIdType, BaseMaterial*> acouMaterials;
    acouMaterials = acousticPDE->GetMaterialData();

    std::map<RegionIdType, BaseMaterial*> waterMaterials;
    waterMaterials = waterWavePDE->GetMaterialData();

    // Create coefficient functions for all acoustic densities
    std::set< RegionIdType > acouRegions;
    std::map<RegionIdType, BaseMaterial*>::iterator itAcou, endAcou;
    itAcou = acouMaterials.begin();
    endAcou = acouMaterials.end();

    std::map< RegionIdType, PtrCoefFct > coefFuncsAcou;
    std::map< RegionIdType, PtrCoefFct > oneCoefFuncsAcou;
    for( ; itAcou != endAcou; itAcou++ ) {
      RegionIdType volRegId = itAcou->first;
      acouRegions.insert(volRegId);

      // Get bulk density for acoustics
      BaseMaterial * acouMat = acouMaterials[volRegId];
      coefFuncsAcou[volRegId] = acouMat->GetScalCoefFnc(DENSITY,Global::REAL);
      oneCoefFuncsAcou[volRegId] = CoefFunction::Generate(mp, Global::REAL, "1.0");
    }

    // Create coefficient functions for all water densities
    std::set< RegionIdType > waterRegions;
    std::map<RegionIdType, BaseMaterial*>::iterator itWater, endWater;
    itWater = waterMaterials.begin();
    endWater = waterMaterials.end();

    std::map< RegionIdType, PtrCoefFct > coefFuncsWater;
    std::map< RegionIdType, PtrCoefFct > oneCoefFuncsWater;
    for( ; itWater != endWater; itWater++ ) {
      RegionIdType volRegId = itWater->first;
      waterRegions.insert(volRegId);

      // Get bulk density for water
      BaseMaterial * waterMat = waterMaterials[volRegId];
      coefFuncsWater[volRegId] = waterMat->GetScalCoefFnc(DENSITY,Global::REAL);
      oneCoefFuncsWater[volRegId] = CoefFunction::Generate(mp, Global::REAL, "1.0");
    }

    shared_ptr<FeSpace> lagrangeMultSpace = lagrangeMultFct->GetFeSpace();
    for ( UInt actSD = 0, n = entityLists_.GetSize(); actSD < n; actSD++ ) {

      shared_ptr<SurfElemList> actSDList =
          dynamic_pointer_cast<SurfElemList>(entityLists_[actSD]);

      RegionIdType region = actSDList->GetRegion();
      lagrangeMultSpace->SetRegionApproximation(region, "default", "default");

      acouFct->AddEntityList(actSDList);
      waterFct->AddEntityList(actSDList);
      lagrangeMultFct->AddEntityList(actSDList);

      //Lagrange to water wave equation
      DefCouplInt( "WaterWaveLagrangeMassCouplingInt", false, 1.0, MASS, waterFct,
          lagrangeMultFct, actSDList, coefFuncsWater, waterRegions );

      //Lagrange to acoustic equation
      DefCouplInt( "AcousticLagrangeMassCouplingInt", false, -1.0, MASS, acouFct,
          lagrangeMultFct, actSDList, coefFuncsAcou, acouRegions );

      //Water to Lagrange equation
      DefCouplInt( "LagrangeWaterWaveStiffCouplingInt", false, 1.0, STIFFNESS, lagrangeMultFct,
          waterFct, actSDList, oneCoefFuncsWater, waterRegions );

      //acoustic to Lagrange equation
      DefCouplInt( "LagrangeAcousticStiffCouplingInt", false, -1.0, STIFFNESS, lagrangeMultFct,
          acouFct, actSDList, oneCoefFuncsAcou, acouRegions );

      //Lagrange to Lagrange equation
      DefCouplInt( "LagrangeLagrangeStiffCouplingInt", false, -9.81, STIFFNESS, lagrangeMultFct,
          lagrangeMultFct, actSDList, coefFuncsWater, waterRegions );

    }
  }


  void WaterWaveAcousticCoupling::DefCouplInt( const std::string& name,
                                        bool assembleTransposed,
                                        Double factor,
                                        FEMatrixType matType,
                                        shared_ptr<BaseFeFunction>& fnc1,
                                        shared_ptr<BaseFeFunction>& fnc2,
                                        shared_ptr<SurfElemList>& actSDList,
                                        const std::map< RegionIdType, PtrCoefFct >& coefFuncs,
                                        const std::set< RegionIdType >& acouRegions ) {

      BiLinearForm * cplInt = NULL;
      if( dim_ == 2  ) {
        cplInt = new SurfaceABInt<>(new IdentityOperator<FeH1,2,1>(),
                                    new IdentityOperator<FeH1,2,1>(),
                                    coefFuncs, factor, acouRegions);
      }
      else if( dim_ == 3  ) {
        cplInt = new SurfaceABInt<>(new IdentityOperator<FeH1,3,1>(),
                                            new IdentityOperator<FeH1,3,1>(),
                                            coefFuncs, factor, acouRegions);
      }
      else {
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
    
    


 void WaterWaveAcousticCoupling::CreateFeSpaces( const std::string&  type,
                                         PtrParamNode infoNode,
                                         std::map<SolutionType, shared_ptr<FeSpace> >& crSpaces) {

   //we need a lagrange multiplier
   formulation_ = LAGRANGE_MULT;

   std::string form = SolutionTypeEnum.ToString(formulation_);
   PtrParamNode spaceNode = infoNode->Get(form);

   crSpaces[formulation_] =
       FeSpace::CreateInstance(myParam_, spaceNode, FeSpace::H1, ptGrid_);

   crSpaces[formulation_]->SetLagrSurfSpace();
   crSpaces[formulation_]->Init(solStrat_);
  }


 void WaterWaveAcousticCoupling::DefineAvailResults() {
   //REFACTOR
 }

 void WaterWaveAcousticCoupling::DefinePrimaryResults() {

   // === LAGRANGE MULTIPLIER ===
   shared_ptr<ResultInfo> res1( new ResultInfo);
   res1->resultType = LAGRANGE_MULT;

   res1->dofNames = "";
   res1->unit = "m";
   res1->definedOn = ResultInfo::NODE;
   res1->entryType = ResultInfo::SCALAR;
   feFunctions_[LAGRANGE_MULT]->SetResultInfo(res1);
   results_.Push_back( res1 );
   availResults_.insert( res1 );

   res1->SetFeFunction(feFunctions_[LAGRANGE_MULT]);
   DefineFieldResult( feFunctions_[LAGRANGE_MULT], res1 );
 }
}

