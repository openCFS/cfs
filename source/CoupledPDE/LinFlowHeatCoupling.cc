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

		// determine whether symmetric or non ymmetric formulation should be used
		isCouplingFormulationSymmetric_ = false;
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

    	//bilinear form for coupling from heat to flow
    	//Regarding the equation of state:
    	//coefFct = -\frac{1}{c}\sqrt{\frac{c_{\rm p}(\gamma-1)}{T_0}}
    	//compute pre-factor
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

        PtrCoefFct hlp2  = CoefFunction::Generate( mp, Global::REAL,CoefXprBinOp(mp,compressionModulus,refTemp,CoefXpr::OP_MULT));

        PtrCoefFct coefThermalExpansion  = CoefFunction::Generate( mp, Global::REAL,
            CoefXprUnaryOp(mp,
                CoefXprBinOp(mp,hlp1 , hlp2, CoefXpr::OP_DIV),
                CoefXpr::OP_SQRT));

    	BiLinearForm *heatToFlow = NULL;
    	if( dim_ == 2 ) {
    		heatToFlow = new ABInt<>(new IdentityOperator<FeH1,2,1>(), new IdentityOperator<FeH1,2,1>(), coefThermalExpansion, -1.0 );
    	} else {
    		heatToFlow = new ABInt<>(new IdentityOperator<FeH1,3,1>(), new IdentityOperator<FeH1,3,1>(), coefThermalExpansion, -1.0 );
    	}
    	heatToFlow->SetName("HeatToLinFlowCoupling");

    	BiLinFormContext * heatToFlowDescr =
				new BiLinFormContext(heatToFlow, DAMPING );

    	heatToFlowDescr->SetEntities( actSDList, actSDList );
    	heatToFlowDescr->SetFeFunctions( flowFct, heatFct );
    	heatToFlowDescr->SetCounterPart(false);

			assemble_->AddBiLinearForm( heatToFlowDescr );

    	// bilinear form for coupling from flow to heat: coefThermalExpansion*refTemp \frac{\partial p_\ra}{\partial t}
    	// The coefficient "ThermalExpansion*refTemp" is necessary for a general fluid.
			// For an ideal gas: ThermalExpansion*refTemp = 1
			// In case the coupling is symmetric the coefficent function is divied by refTemp
			BiLinearForm *flowToHeat = NULL;
			if (isCouplingFormulationSymmetric_) {
				if( dim_ == 2 ) {
					flowToHeat = new ABInt<>(new IdentityOperator<FeH1,2,1>(), new IdentityOperator<FeH1,2,1>(), coefThermalExpansion, -1.0 );
				} else {
					flowToHeat = new ABInt<>(new IdentityOperator<FeH1,3,1>(), new IdentityOperator<FeH1,3,1>(), coefThermalExpansion, -1.0 );
				}
			} else {
				PtrCoefFct coefThermalExpansionT  = CoefFunction::Generate( mp, Global::REAL,CoefXprBinOp(mp,coefThermalExpansion,refTemp,CoefXpr::OP_MULT));
				if( dim_ == 2 ) {
					flowToHeat = new ABInt<>(new IdentityOperator<FeH1,2,1>(), new IdentityOperator<FeH1,2,1>(), coefThermalExpansionT, -1.0 );
				} else {
					flowToHeat = new ABInt<>(new IdentityOperator<FeH1,3,1>(), new IdentityOperator<FeH1,3,1>(), coefThermalExpansionT, -1.0 );
				}
			}

    	heatToFlow->SetName("LinFlowToHeatCoupling");

    	BiLinFormContext * flowToHeatDescr = new BiLinFormContext(flowToHeat, DAMPING );

    	flowToHeatDescr->SetEntities( actSDList, actSDList );
    	flowToHeatDescr->SetFeFunctions( heatFct, flowFct );
    	flowToHeatDescr->SetCounterPart(false);

			assemble_->AddBiLinearForm( flowToHeatDescr );
    }
  }
}
