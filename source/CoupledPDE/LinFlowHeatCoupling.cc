// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "LinFlowHeatCoupling.hh"

#include "PDE/SinglePDE.hh"
#include "PDE/HeatPDE.hh"
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

// include fespaces
#include "FeBasis/H1/H1Elems.hh"

// new integrator concept
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/ConvectiveOperator.hh"

namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  LinFlowHeatCoupling::LinFlowHeatCoupling( SinglePDE *pde1, SinglePDE *pde2,
                                            PtrParamNode paramNode,
                                            PtrParamNode infoNode,
                                            shared_ptr<SimState> simState,
                                            Domain* domain)
    : BasePairCoupling( pde1, pde2, paramNode, infoNode, simState, domain )
  {
    couplingName_ = "linFlowHeatDirect";
    materialClass_ = FLOW;  //we need just material parameters from
                            //the two individual PDEs

    // determine subtype
    pde1_->GetParamNode()->GetValue( "subType", subType_ );

    // determine whether symmetric or non symmetric formulation should be used
    isCouplingFormulationSymmetric_ = true;
    paramNode->GetValue("symmetric",isCouplingFormulationSymmetric_,ParamNode::PASS);

    nonLin_ = false;

    // Initialize nonlinearities
    InitNonLin();

    // inform linFlowPDE about coupling to HeatPDE and vice versa
    dynamic_cast<LinFlowPDE*> (pde1)->SetHeatPDECouplingFlags(isCouplingFormulationSymmetric_);
    dynamic_cast<HeatPDE*> (pde2)->SetLinFlowPDECouplingFlags(isCouplingFormulationSymmetric_);

  }


  // **************
  //   Destructor
  // **************
  LinFlowHeatCoupling::~LinFlowHeatCoupling() {
  }


  // *********************
  //   DefineIntegrators
  // *********************
  void LinFlowHeatCoupling::DefineIntegrators() {

    // get math parser
    MathParser * mp = domain_->GetMathParser();

    shared_ptr<BaseFeFunction> flowFct = pde1_->GetFeFunction(FLUIDMECH_PRESSURE);
    shared_ptr<BaseFeFunction> heatFct = pde2_->GetFeFunction(HEAT_TEMPERATURE);
    std::map<RegionIdType, BaseMaterial*> flowMaterial, heatMaterial;
    flowMaterial = pde1_->GetMaterialData();
    heatMaterial = pde2_->GetMaterialData();

    shared_ptr<FeSpace> pressSpace = flowFct->GetFeSpace();
    shared_ptr<FeSpace> heatSpace  = heatFct->GetFeSpace();

    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      // Set current region and material
      RegionIdType actRegion = it->first;

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( actRegion );

      // Get current region name
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region",
                                                                      "name",
                                                                      regionName.c_str());

      // ====================================================================
      // Existing thermo-acoustic coupling coefficients
      // ====================================================================
      // Regarding the equation of state, the symmetric coupling coefficient is
      // sqrt(rho * cp * (gamma - 1) / (K * T0)).
      // The signs are applied explicitly in the corresponding integrators.
      PtrCoefFct constMinusOne = CoefFunction::Generate( mp, Global::REAL, "-1.0");
      PtrCoefFct refTemp  = heatMaterial[actRegion]->GetScalCoefFnc(HEAT_REF_TEMPERATURE, Global::REAL);
      PtrCoefFct density  = flowMaterial[actRegion]->GetScalCoefFnc(DENSITY, Global::REAL);
      PtrCoefFct heatCapacity = heatMaterial[actRegion]->GetScalCoefFnc( HEAT_CAPACITY, Global::REAL );
      PtrCoefFct adiabaticExp = flowMaterial[actRegion]->GetScalCoefFnc(FLUID_ADIABATIC_EXPONENT, Global::REAL);
      PtrCoefFct compressionModulus  = flowMaterial[actRegion]->GetScalCoefFnc(FLUID_BULK_MODULUS, Global::REAL);

      PtrCoefFct hlp1  = CoefFunction::Generate( mp, Global::REAL,
          CoefXprBinOp(mp,density,
              CoefXprBinOp(mp,heatCapacity,
                  CoefXprBinOp(mp,constMinusOne,adiabaticExp, CoefXpr::OP_ADD ),
                  CoefXpr::OP_MULT),
                  CoefXpr::OP_MULT ));

      PtrCoefFct hlp2  = CoefFunction::Generate( mp, Global::REAL,
          CoefXprBinOp(mp,compressionModulus,refTemp,CoefXpr::OP_MULT));

      PtrCoefFct coefThermalExpansion  = CoefFunction::Generate( mp, Global::REAL,
          CoefXprUnaryOp(mp,
              CoefXprBinOp(mp,hlp1 , hlp2, CoefXpr::OP_DIV),
              CoefXpr::OP_SQRT));

      // Coefficient used for the explicit LinFlow -> Heat direction.
      // For the symmetric formulation the cross coefficient is identical in
      // both directions. In the non-symmetric formulation the existing Heat
      // equation uses the additional factor T0.
      PtrCoefFct coefFlowToHeat = coefThermalExpansion;
      if (!isCouplingFormulationSymmetric_) {
        coefFlowToHeat = CoefFunction::Generate( mp, Global::REAL,
            CoefXprBinOp(mp,coefThermalExpansion,refTemp,CoefXpr::OP_MULT));
      }

      // ====================================================================
      // Existing DAMPING coupling: Heat -> LinFlow
      //     - g * dT/dt
      // ====================================================================
      BiLinearForm *heatToFlow = NULL;
      if( dim_ == 2 ) {
        heatToFlow = new ABInt<>(new IdentityOperator<FeH1,2,1>(),
                                 new IdentityOperator<FeH1,2,1>(),
                                 coefThermalExpansion, -1.0 );
      } else {
        heatToFlow = new ABInt<>(new IdentityOperator<FeH1,3,1>(),
                                 new IdentityOperator<FeH1,3,1>(),
                                 coefThermalExpansion, -1.0 );
      }
      heatToFlow->SetName("HeatToLinFlowCoupling");

      BiLinFormContext * heatToFlowDescr =
          new BiLinFormContext(heatToFlow, DAMPING );

      heatToFlowDescr->SetEntities( actSDList, actSDList );
      heatToFlowDescr->SetFeFunctions( flowFct, heatFct );
      // In case the coupling is written in a symmetric form, the counterpart
      // creates the second damping block automatically.
      heatToFlowDescr->SetCounterPart(isCouplingFormulationSymmetric_);

      assemble_->AddBiLinearForm( heatToFlowDescr );

      // ====================================================================
      // Existing DAMPING coupling: LinFlow -> Heat (non-symmetric only)
      //     - g*T0 * dp/dt
      // ====================================================================
      if (!isCouplingFormulationSymmetric_) {
        BiLinearForm *flowToHeat = NULL;
        if( dim_ == 2 ) {
          flowToHeat = new ABInt<>(new IdentityOperator<FeH1,2,1>(),
                                   new IdentityOperator<FeH1,2,1>(),
                                   coefFlowToHeat, -1.0 );
        } else {
          flowToHeat = new ABInt<>(new IdentityOperator<FeH1,3,1>(),
                                   new IdentityOperator<FeH1,3,1>(),
                                   coefFlowToHeat, -1.0 );
        }

        flowToHeat->SetName("LinFlowToHeatCoupling");

        BiLinFormContext * flowToHeatDescr =
            new BiLinFormContext(flowToHeat, DAMPING );

        flowToHeatDescr->SetEntities( actSDList, actSDList );
        flowToHeatDescr->SetFeFunctions( heatFct, flowFct );
        flowToHeatDescr->SetCounterPart(false);

        assemble_->AddBiLinearForm( flowToHeatDescr );
      }

      // ====================================================================
      // ALE LinFlow-Heat coupling
      // ====================================================================
      // The full movingMeshLists of LinFlowPDE and HeatPDE do not have to be
      // identical. Only the moving-mesh definition on the CURRENT COUPLED
      // REGION must be consistent.
      PtrParamNode flowRegNode = pde1_->GetParamNode()->Get("regionList")->GetByVal(
          "region", "name", regionName.c_str());
      PtrParamNode heatRegNode = pde2_->GetParamNode()->Get("regionList")->GetByVal(
          "region", "name", regionName.c_str());

      std::string flowMovingMeshId = "";
      std::string heatMovingMeshId = "";
      flowRegNode->GetValue("movingMeshId", flowMovingMeshId, ParamNode::PASS);
      heatRegNode->GetValue("movingMeshId", heatMovingMeshId, ParamNode::PASS);

      const bool flowHasMovingMesh = (flowMovingMeshId != "");
      const bool heatHasMovingMesh = (heatMovingMeshId != "");

      // A coupled ALE region must be described in the same frame by both PDEs.
      if (flowHasMovingMesh != heatHasMovingMesh) {
        EXCEPTION("LinFlowHeatCoupling: inconsistent moving-mesh definition on coupled region '"
                  << regionName << "'. LinFlow movingMeshId='" << flowMovingMeshId
                  << "', Heat movingMeshId='" << heatMovingMeshId
                  << "'. Define movingMeshId for both PDEs or for neither PDE.")
      }

      // For this first implementation we use a strict and testable contract:
      // both PDEs must reference the same movingMeshId on a coupled region.
      // The two global movingMeshLists may still contain different additional
      // entries for uncoupled regions.
      if (flowHasMovingMesh && flowMovingMeshId != heatMovingMeshId) {
        EXCEPTION("LinFlowHeatCoupling: different movingMeshIds on coupled region '"
                  << regionName << "' (LinFlow='" << flowMovingMeshId
                  << "', Heat='" << heatMovingMeshId
                  << "'). The coupled ALE formulation requires one common grid velocity.")
      }

      if (flowHasMovingMesh) {
        if (analysisType_ != BasePDE::TRANSIENT) {
          EXCEPTION("LinFlowHeatCoupling: ALE LinFlow-Heat coupling is currently "
                    "implemented only for transient analysis.")
        }

        // Both PDEs register FLUIDMECH_MESH_VELOCITY as a coefficient field.
        // The LinFlow field is used as the common v_g in all coupling terms;
        // the Heat field is checked to ensure that the ALE data are available
        // on both sides of the coupling.
        PtrCoefFct flowGridVelocity = pde1_->GetCoefFct(FLUIDMECH_MESH_VELOCITY);
        PtrCoefFct heatGridVelocity = pde2_->GetCoefFct(FLUIDMECH_MESH_VELOCITY);

        if (!flowGridVelocity || !heatGridVelocity) {
          EXCEPTION("LinFlowHeatCoupling: movingMeshId is set on coupled region '"
                    << regionName << "', but FLUIDMECH_MESH_VELOCITY is not available "
                    "from both LinFlowPDE and HeatPDE.")
        }

        // ------------------------------------------------------------------
        // K_PP^ALE: coupled pressure grid transport
        //     -(gamma/K) * v_g . grad(p)
        // ------------------------------------------------------------------
        // LinFlowPDE uses gamma/K in C_PP when HeatPDE is coupled. Therefore
        // the ALE counterpart has to use the same coefficient with opposite
        // spatial transport sign.
        PtrCoefFct coefGridPressure = CoefFunction::Generate( mp, Global::REAL,
            CoefXprBinOp(mp,adiabaticExp,compressionModulus,CoefXpr::OP_DIV));

        BiLinearForm *gridPressureToPressure = NULL;
        if( dim_ == 2 ) {
          gridPressureToPressure = new ABInt<>(
              new IdentityOperator<FeH1,2,1>(),
              new ConvectiveOperator<FeH1,2,1>(),
              coefGridPressure, -1.0, true );
        } else {
          gridPressureToPressure = new ABInt<>(
              new IdentityOperator<FeH1,3,1>(),
              new ConvectiveOperator<FeH1,3,1>(),
              coefGridPressure, -1.0, true );
        }
        gridPressureToPressure->SetBCoefFunctionOpB(flowGridVelocity);
        gridPressureToPressure->SetSolDependent(true);
        gridPressureToPressure->SetName("LinFlowHeatALEGridPressurePP");

        BiLinFormContext *gridPressureToPressureDescr =
            new BiLinFormContext(gridPressureToPressure, STIFFNESS);
        gridPressureToPressureDescr->SetEntities(actSDList, actSDList);
        gridPressureToPressureDescr->SetFeFunctions(flowFct, flowFct);
        gridPressureToPressureDescr->SetCounterPart(false);
        assemble_->AddBiLinearForm(gridPressureToPressureDescr);

        // ------------------------------------------------------------------
        // K_PT^ALE: Heat -> LinFlow grid coupling
        //     + g * v_g . grad(T)
        // ------------------------------------------------------------------
        BiLinearForm *gridHeatToFlow = NULL;
        if( dim_ == 2 ) {
          gridHeatToFlow = new ABInt<>(
              new IdentityOperator<FeH1,2,1>(),
              new ConvectiveOperator<FeH1,2,1>(),
              coefThermalExpansion, 1.0, true );
        } else {
          gridHeatToFlow = new ABInt<>(
              new IdentityOperator<FeH1,3,1>(),
              new ConvectiveOperator<FeH1,3,1>(),
              coefThermalExpansion, 1.0, true );
        }
        gridHeatToFlow->SetBCoefFunctionOpB(flowGridVelocity);
        gridHeatToFlow->SetSolDependent(true);
        gridHeatToFlow->SetName("HeatToLinFlowALEGridCoupling");

        BiLinFormContext *gridHeatToFlowDescr =
            new BiLinFormContext(gridHeatToFlow, STIFFNESS);
        gridHeatToFlowDescr->SetEntities(actSDList, actSDList);
        gridHeatToFlowDescr->SetFeFunctions(flowFct, heatFct);
        // The convective counterpart is NOT generated automatically because
        // v_g.grad(T) and v_g.grad(p) are separate non-transposed operators.
        gridHeatToFlowDescr->SetCounterPart(false);
        assemble_->AddBiLinearForm(gridHeatToFlowDescr);

        // ------------------------------------------------------------------
        // K_TP^ALE: LinFlow -> Heat grid coupling
        //     + g * v_g . grad(p)          (symmetric formulation)
        //     + g*T0 * v_g . grad(p)       (non-symmetric formulation)
        // ------------------------------------------------------------------
        BiLinearForm *gridFlowToHeat = NULL;
        if( dim_ == 2 ) {
          gridFlowToHeat = new ABInt<>(
              new IdentityOperator<FeH1,2,1>(),
              new ConvectiveOperator<FeH1,2,1>(),
              coefFlowToHeat, 1.0, true );
        } else {
          gridFlowToHeat = new ABInt<>(
              new IdentityOperator<FeH1,3,1>(),
              new ConvectiveOperator<FeH1,3,1>(),
              coefFlowToHeat, 1.0, true );
        }
        gridFlowToHeat->SetBCoefFunctionOpB(flowGridVelocity);
        gridFlowToHeat->SetSolDependent(true);
        gridFlowToHeat->SetName("LinFlowToHeatALEGridCoupling");

        BiLinFormContext *gridFlowToHeatDescr =
            new BiLinFormContext(gridFlowToHeat, STIFFNESS);
        gridFlowToHeatDescr->SetEntities(actSDList, actSDList);
        gridFlowToHeatDescr->SetFeFunctions(heatFct, flowFct);
        gridFlowToHeatDescr->SetCounterPart(false);
        assemble_->AddBiLinearForm(gridFlowToHeatDescr);
      }
    }
  }
}
