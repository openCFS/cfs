#include "ThermoMechCoupling.hh"


#include "DataInOut/ParamHandling/ParamNode.hh"

#include "Driver/assemble.hh"
#include "Materials/baseMaterial.hh"
#include "PDE/SinglePDE.hh"

// integrator (bi-)linear forms
#include "Forms/linThermoMechInt.hh"
#include "Forms/linThermoMechDampInt.hh"

#include "PDE/heatCondPDE.hh"

namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  ThermoMechCoupling::ThermoMechCoupling( SinglePDE *pde1, SinglePDE *pde2,
                                      PtrParamNode paramNode  )
    : BasePairCoupling( pde1, pde2, paramNode, "thermoMechDirect") {

    materialClass_ = THERMOELASTIC;
    

	// non linearities have been not yet implemented.
	nonLinThermoMechCoupling_ = false;


    try {
		// Check the subtype of the mechanic problem
		// (full, normal or axi-symmetric)
		pde1_->GetParamNode()->GetValue( "subType", subType_ );
	}
	catch (Exception& e) {
		RETHROW_EXCEPTION(e, "Could not get the parameter subType_ on "
				<<"ThermoMechCoupling constructor " );
	}
	

  }


  // **************
  //   Destructor
  // **************
  ThermoMechCoupling::~ThermoMechCoupling() {
  }


  // ***************
  //   CalcResults
  // ***************
  void ThermoMechCoupling::CalcResults() {
  }


  // ***************
  //   Output of Results
  // ***************
  void ThermoMechCoupling::WriteResultsInFile(const UInt kstep,
                                            const Double asteptime,
                                            UInt stepOffset,
                                            Double timeOffset){
  }

  

  // *********************
  //   DefineIntegrators
  // *********************
  void ThermoMechCoupling::DefineIntegrators() {


	//Note:
	//pde1 => mechanic
	//pde2 => heatConduction


	// Calculate global thermo-mech integral
	Global::ComplexPart matType = Global::REAL;
	RegionIdType actRegion;
	// BaseMaterial * actSDMat= NULL; // TODO: Unused variable actSDMat

	try {

		// get material from mechanics 
		//** checar despues si debe estar dentro del ciclo de materiales o no
		std::map<RegionIdType, BaseMaterial*> mechMat = pde1_->getPDEMaterialData();

		// Determine which analysis is performed.
		// Flag to calculate the damping matrix
		bool dampFlag=false;
		if (pde2_->GetAnalysisType() == BasePDE::TRANSIENT) 
		{
			// create the damp thermo-mechanic block matrix
			dampFlag = true;
		}
		
		// Define integrators for "standard" materials
		std::map<RegionIdType, BaseMaterial*>::iterator it;
		for (it = materials_.begin(); it != materials_.end(); it++) {
			actRegion = it->first; // set current region
			// actSDMat = it->second; // set current material
			matType = Global::REAL;

			//transform the type
			SubTensorType tensorType;
			String2Enum(subType_, tensorType);

			// create new entity list
			shared_ptr<ElemList> actSDList(new ElemList(ptGrid_ ));
			actSDList->SetRegion(actRegion );
			
			// add stiffness matrix ----------------------------------------
			BaseForm *bilinearStiff= NULL;
			BiLinFormContext *actContextStiff= NULL;

			bilinearStiff =new LinThermoMechInt(materials_[actRegion], mechMat[actRegion], tensorType);

			bilinearStiff->SetMatDataType(matType );

			actContextStiff =new BiLinFormContext( bilinearStiff, STIFFNESS );

			actContextStiff->SetEntryType(matType );
			actContextStiff->SetPtPdes(pde1_, pde2_ );
			actContextStiff->SetResults(results1_[0], results2_[0], actSDList,
					actSDList );

			//		We don't need to set the transposed of the coupling
			//		matrix to the lower diagonal side since there is only one  
			//		coupling, here thermal to mech but not mech to thermal, the other 
			//		coupling is zero.
			actContextStiff->SetCounterPart( false);

			//		Negate the result stiffness matrix
			actContextStiff->SetNegate(true);

			assemble_->AddBiLinearForm(actContextStiff );

			//--------------------------------------------------------------


			if (dampFlag == true) {

				// add damping matrix --------------------------------------

				//std::cout << "\n working on TRANSIENT thermo-electrostatics ..... "
				//		<< std::endl;

				BaseForm *bilinearDamping= NULL;
				BiLinFormContext *actContextDamping= NULL;

				//Call the thermo-mech damp integrator 
				bilinearDamping =new LinThermoMechDampInt(materials_[actRegion], mechMat[actRegion], tensorType);

				// Read and Set the initial conditions for the heat conduction pde
				//dynamic_cast<HeatCondPDE*>(pde2_)->SetInitialCondition();
				
				// Get the heat conduction solution
				BaseNodeStoreSol * solPDE2 = pde2_->getPDESolution();
				NodeStoreSol<Double>
				* solhelp2 =dynamic_cast<NodeStoreSol<Double>*>(solPDE2);

				// keep the heat conduction solution in the thermo-electric damp integrator 
				bilinearDamping->SetSolution2(*solhelp2);

				bilinearDamping->SetMatDataType(matType );

				actContextDamping =new BiLinFormContext( bilinearDamping, DAMPING );

				actContextDamping->SetEntryType(matType );
				
				actContextDamping->SetPtPdes(pde2_, pde1_ );
				actContextDamping->SetResults(results2_[0], results1_[0],
						actSDList, actSDList );

				//		We don't need to set the transposed of the coupling
				//		matrix to the lower diagonal side since there is only 
				//		one coupling, here thermal to elec but not elec to 
				//		thermal, the other coupling is zero. 
				actContextDamping->SetCounterPart( false);

				assemble_->AddBiLinearForm(actContextDamping );
			}
		} // end linear material-iterator

	}
	catch (Exception& e) {
		RETHROW_EXCEPTION(e, "Could not define the integrator on"
				<<" ThermoMechCoupling " );
	}
}


  // ********************
  //   ReadStoreResults
  // ********************
  void ThermoMechCoupling::ReadStoreResults() {
  }




}
