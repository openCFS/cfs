#include "ThermoElectricCoupling.hh"


#include "DataInOut/ParamHandling/ParamNode.hh"

#include "Driver/assemble.hh"
#include "Materials/baseMaterial.hh"
#include "PDE/SinglePDE.hh"


// integrator (bi-)linear forms
#include "Forms/linThermoElectricInt.hh"
#include "Forms/linThermoElecDampInt.hh"


#include "PDE/heatCondPDE.hh"

namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  ThermoElectricCoupling::ThermoElectricCoupling( SinglePDE *pde1, SinglePDE *pde2,
                                      PtrParamNode paramNode  )
    : BasePairCoupling( pde1, pde2, paramNode, "thermoElectricDirect" ) {

    materialClass_ = PYROELECTRIC;

	// non linearities have been not yet implemented.
	nonLinThermoElecCoupling_ = false;

   }


  // **************
  //   Destructor
  // **************
  ThermoElectricCoupling::~ThermoElectricCoupling() {
  }


  // ***************
  //   CalcResults
  // ***************
  void ThermoElectricCoupling::CalcResults() {
  }


  // ***************
  //   Output of Results
  // ***************
  void ThermoElectricCoupling::WriteResultsInFile(const UInt kstep,
                                            const Double asteptime,
                                            UInt stepOffset,
                                            Double timeOffset){
  }



  // *********************
  //   DefineIntegrators
  // *********************
  void ThermoElectricCoupling::DefineIntegrators() {


	//Note:
	//pde1 => electrostatics
	//pde2 => heatConduction

	//std::cout << "\nworking on thermo-electrostatics ..... " << std::endl;

	// Calculate global thermo-elec integral
	Global::ComplexPart matType = Global::REAL;
	RegionIdType actRegion;
	// BaseMaterial * actSDMat= NULL; // TODO: Unused variable actSDMat

	try {
		// Determine, if which geometry is used
		std::string geometryType_;
		param->Get("domain")->GetValue("geometryType", geometryType_ );

		// convert to tensor type
		// COMPWARNING: initialized to FULL (lacking intelligent alternatives...)
		SubTensorType tensorType = FULL;
		if (geometryType_ == "3d") {
			tensorType = FULL;
		} else if (geometryType_ == "plane") {
			tensorType = PLANE_STRAIN;
		} else if (geometryType_ == "axi") {
			tensorType = AXI;
			isaxi_ = true;
		}

		// Determine which analysis is performed.
		// Flag to calculate the damping matrix
		bool dampFlag=false;
		if (pde2_->GetAnalysisType() == BasePDE::TRANSIENT) {
			// create the damp thermoelectric block matrix
			dampFlag = true;
		}

		// Define integrators for "standard" materials
		std::map<RegionIdType, BaseMaterial*>::iterator it;
		for (it = materials_.begin(); it != materials_.end(); it++) {
			actRegion = it->first; // set current region
			// actSDMat = it->second; // set current material
			matType = Global::REAL;

			// create new entity list
			shared_ptr<ElemList> actSDList(new ElemList(ptGrid_ ));
			actSDList->SetRegion(actRegion );
			
			// heat initial conditions
			//dynamic_cast<HeatCondPDE*>(pde2_)->SetInitCondition();

			// add stiffness matrix ----------------------------------------

			BaseForm *bilinearStiff= NULL;
			BiLinFormContext *actContextStiff= NULL;
			bilinearStiff =new LinThermoElectricInt(materials_[actRegion], tensorType);

			bilinearStiff->SetMatDataType(matType );

			actContextStiff =new BiLinFormContext( bilinearStiff, STIFFNESS );

			actContextStiff->SetEntryType(matType );
			actContextStiff->SetPtPdes(pde1_, pde2_ );
			actContextStiff->SetResults(results1_[0], results2_[0], actSDList,
					actSDList );

			//		We don't need to set the transposed of the coupling
			//		matrix to the lower diagonal side since there is only 
			//		one coupling, here thermal to elec but not elec to 
			//		thermal, the other coupling is zero. 
			actContextStiff->SetCounterPart( false);
			
			assemble_->AddBiLinearForm(actContextStiff );

			//--------------------------------------------------------------

			if (dampFlag == true) {

				// add damping matrix --------------------------------------


				//std::cout << "\n working on TRANSIENT thermo-electrostatics ..... "
				//		<< std::endl;

				LinThermoElecDampInt *bilinearDamping= NULL;
				BiLinFormContext *actContextDamping= NULL;

				//Call the thermo-electric damp integrator 
				bilinearDamping =new LinThermoElecDampInt(materials_[actRegion], tensorType);

				// Read and Set the initial conditions for the heat conduction pde
				//dynamic_cast<HeatCondPDE*>(pde2_)->SetInitialCondition();

				//get the heat conduction solution
				BaseNodeStoreSol * solPDE2 = pde2_->getPDESolution();
				NodeStoreSol<Double>
				* solhelp2 =dynamic_cast<NodeStoreSol<Double>*>(solPDE2);

				//keep the heat conduction solution in the thermo-electric damp integrator 
				bilinearDamping->SetSolution2(*solhelp2);

				bilinearDamping->SetMatDataType(matType );

				actContextDamping =new BiLinFormContext( bilinearDamping, DAMPING );

				actContextDamping->SetEntryType(matType );
				//			actContextDamping->SetPtPdes( pde1_, pde2_ );
				//			actContextDamping->SetResults( results1_[0], results2_[0],
				//											actSDList, actSDList );

				// now the parameters pde2_ and pde1_ are exchanged in order
				// to assembly in the row block of heat conduction pde
				actContextDamping->SetPtPdes(pde2_, pde1_ );
				actContextDamping->SetResults(results2_[0], results1_[0],
						actSDList, actSDList );

				//		We don't need to set the transposed of the coupling
				//		matrix to the lower diagonal side since there is only 
				//		one coupling, here thermal to elec but not elec to 
				//		thermal, the other coupling is zero. 
				actContextDamping->SetCounterPart( false);
		
				//		Negate the result stiffness matrix
				actContextDamping->SetNegate(true);

				assemble_->AddBiLinearForm(actContextDamping );
			}

		} // end material-iterator
	}
	catch (Exception& e) {
		RETHROW_EXCEPTION(e, "Could not define the integrator on"
				<<" ThermoElectricCoupling " );
	}

}


  // ********************
  //   ReadStoreResults
  // ********************
  void ThermoElectricCoupling::ReadStoreResults() {
  }




}
