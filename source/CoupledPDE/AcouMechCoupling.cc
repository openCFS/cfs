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
                                PtrParamNode paramNode, 
                                PtrParamNode infoNode,
                                shared_ptr<SimState> simState, Domain* domain)
    : BasePairCoupling( pde1, pde2, paramNode, infoNode, simState, domain )
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
                       DAMPING, acouFct, dispFct,  
                       actSDList, coefFuncs, acouRegions );
          }
        }
        break;

        // ==================
        //   ACOU_PRESSURE
        // ==================
        case ACOU_PRESSURE:
          
          // First ensure, that we have no coupled eigenfrequency simulation:
          // This case is not tested yet and might not be solvalbe.
          if( analysisType_ == BasePDE::EIGENFREQUENCY ) {
            EXCEPTION("A coupled mechanic-acoustic simulation can only be"
                      "be performed in the acoustic potential formulation!");
          }
          DefCouplInt( "AcouMechPresStiffCouplingInt", false, -1.0, STIFFNESS, dispFct, 
                       acouFct, actSDList, oneCoefFuncs, acouRegions );
          
          DefCouplInt( "AcouMechPresMassCouplingInt", false, 1.0, MASS, acouFct, 
                       dispFct, actSDList, coefFuncs, acouRegions );
         
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
  
}
