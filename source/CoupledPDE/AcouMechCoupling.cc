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

// include fespaces
#include "FeBasis/H1/FeSpaceH1.hh"
#include "FeBasis/H1/H1Elems.hh"

// new integrator concept
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/IdentityOperatorNormal.hh"


namespace CoupledField {



  // ***************
  //   Constructor
  // ***************
  AcouMechCoupling::AcouMechCoupling( SinglePDE *pde1, SinglePDE *pde2,
                                PtrParamNode paramNode  )
    : BasePairCoupling( pde1, pde2, paramNode )
  {
    couplingName_ = "acouMechDirect";
    materialClass_ = FLUID;

    // determine subtype from mechanic pde
    pde1_->GetParamNode()->GetValue( "subType", subType_ );

    nonLin_ = false;
    
    // Initialize nonlinearities
    InitNonLin();
  }


  // **************
  //   Destructor
  // **************
  AcouMechCoupling::~AcouMechCoupling() {
  }


  // ***************
  //   CalcResults
  // ***************
  void AcouMechCoupling::CalcResults( shared_ptr<BaseResult> result ) {
    REFACTOR
  }


 

  // *********************
  //   DefineIntegrators
  // *********************
  void AcouMechCoupling::DefineIntegrators() {

    // get hold of both feFunctions
    MechPDE* mechPDE = dynamic_cast<MechPDE*>(pde1_);
    shared_ptr<BaseFeFunction> dispFct = mechPDE->GetFeFunction(MECH_DISPLACEMENT);
//    std::map<RegionIdType, BaseMaterial*> mechMaterials;
//    mechMaterials = mechPDE->getPDEMaterialData();

    AcousticPDE* acouPDE = dynamic_cast<AcousticPDE*>(pde2_);
    SolutionType acouFormulation = acouPDE->GetFormulation();
    shared_ptr<BaseFeFunction> acouFct = acouPDE->GetFeFunction(acouFormulation);
    std::map<RegionIdType, BaseMaterial*> acouMaterials;
    acouMaterials = acouPDE->getPDEMaterialData();

    // Create coefficient functions for all acoustic densities
    std::map< RegionIdType, shared_ptr<CoefFunction> > coefFuncs;
    std::set< RegionIdType > acouRegions;
    std::map<RegionIdType, BaseMaterial*>::iterator it, end;
    it = acouMaterials.begin();
    end = acouMaterials.end();
    for( ; it != end; it++ ) {
      RegionIdType volRegId = it->first;

      acouRegions.insert(volRegId);

      // Get bulk density for acoustics
      BaseMaterial * acouMat = acouMaterials[volRegId];
      Double density = 1.0;
      acouMat->GetScalar(density,DENSITY,Global::REAL);

      coefFuncs[volRegId] = CoefFunction::Generate(Global::REAL,
                                                   lexical_cast<std::string>(density));
    }

    shared_ptr<FeSpace> dispSpace = dispFct->GetFeSpace();
    shared_ptr<FeSpace> acouSpace = acouFct->GetFeSpace();
    
    for ( UInt actSD = 0, n = entityLists_.GetSize(); actSD < n; actSD++ ) {

      shared_ptr<SurfElemList> actSDList(dynamic_cast<SurfElemList*>(entityLists_[actSD].get()));

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
        case ACOU_POTENTIAL:
        {
          DefineDampingIntegratorsTempl<Double>(dispFct, acouFct, actSDList, coefFuncs, acouRegions);
        }
        break;

        case ACOU_PRESSURE:
          break;

        default:
          EXCEPTION("Unknown formulation for acoustics.");
          break;
      }
    }
  }

  template<class DATA_TYPE>
  void AcouMechCoupling::DefineDampingIntegratorsTempl(shared_ptr<BaseFeFunction>& dispFct,
                                                       shared_ptr<BaseFeFunction>& acouFct,
                                                       shared_ptr<SurfElemList>& actSDList,
                                                       const std::map< RegionIdType, shared_ptr<CoefFunction> >& coefFuncs,
                                                       const std::set< RegionIdType >& acouRegions){

    BiLinearForm * dampInt = NULL;

    if( subType_ == "axi" ) {
      dampInt = new SurfaceABInt< IdentityOperator<FeH1,2,2, DATA_TYPE>,
                                  IdentityOperatorNormal<FeH1,2, DATA_TYPE>,
                                  DATA_TYPE >(coefFuncs, DATA_TYPE(-1.0), acouRegions);
    } else if( subType_ == "planeStrain" ) {
      dampInt = new SurfaceABInt< IdentityOperator<FeH1,2,2, DATA_TYPE>,
                                 IdentityOperatorNormal<FeH1,2, DATA_TYPE>,
                                 DATA_TYPE >(coefFuncs, DATA_TYPE(-1.0), acouRegions);
    } else if( subType_ == "planeStress" ) {
      dampInt = new SurfaceABInt< IdentityOperator<FeH1,2,2, DATA_TYPE>,
                                 IdentityOperatorNormal<FeH1,2, DATA_TYPE>,
                                 DATA_TYPE >(coefFuncs, DATA_TYPE(-1.0), acouRegions);
    } else if( subType_ == "3d") {
      dampInt = new SurfaceABInt< IdentityOperator<FeH1,3,3, DATA_TYPE>,
                                 IdentityOperatorNormal<FeH1,3, DATA_TYPE>,
                                 DATA_TYPE >(coefFuncs, DATA_TYPE(-1.0), acouRegions);
    } else {
      EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
    }

    dampInt->SetName("AcouMechPotCouplingInt");
    BiLinFormContext * dampContext =
        new BiLinFormContext(dampInt, DAMPING );

    dampContext->SetEntities( actSDList, actSDList );
    dampContext->SetFeFunctions( dispFct, acouFct );
    dampContext->SetCounterPart(true);

    assemble_->AddBiLinearForm( dampContext );
  }

  void AcouMechCoupling::DefineAvailResults() {
    REFACTOR  
  }

  BiLinearForm *
  AcouMechCoupling::GetStiffIntegrator( BaseMaterial* actSDMat,
                                     RegionIdType regionId,
                                     bool isComplex ) {
#if 0
    // Get region name
    std::string regionName = ptGrid_->GetRegion().ToString( regionId );

    SubTensorType tensorType = NO_TENSOR;
    if( subType_ != "3d") {
      tensorType = PLANE_STRAIN;
    } else {
      tensorType = FULL;
    }
    // ------------------------
    //  Obtain linear material
    // ------------------------
    shared_ptr<CoefFunction > curCoef;
    if( isComplex ) {
      curCoef = actSDMat->GetCoefFunction(PIEZO_TENSOR,
                                          tensorType, Global::COMPLEX, true );
    } else {
      curCoef = actSDMat->GetCoefFunction(PIEZO_TENSOR,
                                          tensorType, Global::REAL, true );
    }
    
    // ----------------------------------------
    //  Determine correct stiffness integrator 
    // ----------------------------------------
    BiLinearForm * integ = NULL;
    if( subType_ == "axi" ) {
      integ = new ADBInt<StrainOperatorAxi<FeH1>,
          GradientOperator<FeH1,2> >(curCoef, 1.0);
    } else if( subType_ == "planeStrain" ) {
      integ = new ADBInt<StrainOperator2D<FeH1>,
          GradientOperator<FeH1,2> >(curCoef, 1.0);
    } else if( subType_ == "planeStress" ) {
      integ = new ADBInt<StrainOperator2D<FeH1>,
          GradientOperator<FeH1,2> >(curCoef, 1.0);
    } else if( subType_ == "3d") {
      integ = new ADBInt<StrainOperator3D<FeH1>,
          GradientOperator<FeH1,3> >(curCoef, 1.0);
    } else {
      EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
    }
    return integ;
#endif
    return NULL;
  }
  
}

