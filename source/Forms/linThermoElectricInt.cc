#include <iostream>
#include <fstream>

#include "linThermoElectricInt.hh"

#include "Domain/domain.hh"
#include "Domain/grid.hh"

namespace CoupledField
{


  LinThermoElectricInt::LinThermoElectricInt(BaseMaterial* matData, SubTensorType type) :
     ADBInt(matData,type) 
  {

    name_ = "LinThermoElectricInt";


   if ( type == FULL ) {
       numDofsA_  = 1;
       numDofsB_  = 1;

       matDimRow_ = 3;
       matDimCol_ = 1;

    }
    else if ( type == PLANE_STRAIN || type == PLANE_STRESS ) {
      numDofsA_  = 1;
      numDofsB_  = 1; 

      matDimRow_ = 2;
      matDimCol_ = 1;
    }
    else if ( type == AXI ) {
        numDofsA_  = 2;
        numDofsB_  = 1;
        matDimRow_ = 2;
        matDimCol_ = 1;
        isaxi_     = true;
    }


  }

 
	LinThermoElectricInt::~LinThermoElectricInt()
	{
	}



    // Query material type for  D  tensor.
	MaterialType LinThermoElectricInt::getDMaterialType() {return PYROCOEFFICIENT_TENSOR;}


	//!  Compute matrix \f$ A \f$ at given integration point.
	void LinThermoElectricInt::calcAMat( Matrix<Double> &aMat, UInt ip,
						const Matrix<Double> &ptCoord ){


		try{
			//std::cout << "\nworking on the thermo-electrostatics integral ..... calcAMat" << std::endl;
		// Obtain info on problem sizes
		UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
		//std::cout << "numFncs" << numFncs << std::endl;
	
		// Set type of ansatz function , but do not recalculate
		// integration points
		ptelem->SetAnsatzFct( ansatzFct1_, false );
	
		// Set correct size of matrix A and initialize it with zeros
		aMat.Resize( numFncs * numDofsA_, matDimRow_ );
		aMat.Init();

		// Get derivatives of local shape functions with respect to global coords
		// (format: numFncs x spaceDim)
		ptelem->GetGlobDerivShFncAtIp( aMat, ip, ptCoord, it1_.GetElem() );

		//std::cerr << "LinThermoElectricInt::aMat = \n" << aMat << std::endl;
		}
		catch (Exception& e) {
			RETHROW_EXCEPTION(e, "Could not calculate A in LinThermoElectricInt" );
		}

	}

	//!  Compute matrix \f$ D \f$ at given integration point.
	void LinThermoElectricInt::calcDMat( Matrix<Double> &dMat)
	{
	  try
	  {
	    //std::cout << "\nworking on the thermo-electrostatics integral ..... calcDMat" << std::endl; 	

	    // get the thermal expansion coefficient tensor
	    Matrix<Double>  p_auxMat(matDimRow_, matDimRow_);
	    ptMaterial->GetTensor(p_auxMat,PYROCOEFFICIENT_TENSOR,matDataType_,subTensorType_);

	    dMat.Resize(matDimRow_, matDimCol_);
	    p_auxMat.GetDiagInMatrix(dMat);

	  }
	  catch (Exception& e) {
	    RETHROW_EXCEPTION(e, "Could not calculate D in LinThermoElectricInt" );
	  }
	}

	//!  Compute matrix \f$ B \f$ at given integration point.
	void LinThermoElectricInt::calcBMat( Matrix<Double> &bMat, UInt ip,
						const Matrix<Double> &ptCoord ){


		try{
		//std::cout << "\nworking on the thermo-electrostatics integral ..... calcBMat" << std::endl; 	

		// Set type of ansatz function , but do not recalculate
		// integration points (=false)
		ptelem->SetAnsatzFct( ansatzFct2_, false );
	
		// Obtain info on number of element's nodes
		UInt numFncs = ptelem->GetNumFncs( ansatzFct2_ );
	
		// Set correct size of matrix B and initialize with zeros
		bMat.Resize(1, numFncs);
		bMat.Init();

		// Get local shape functions of thermal field (get sure they're from thermal field)
		// (format: numFncs x 1)
		Vector< Double > lN_teta;
		ptelem->GetShFncAtIp( lN_teta, ip, it1_.GetElem() );

		//std::cout << "\nlN_teta=" << lN_teta.Serialize () << std::endl;

		for(UInt i=0; i < bMat.GetSizeCol(); i++ ){
			bMat[0][i] = lN_teta[i];
		}

		//std::cerr << "bMat = \n" << bMat << std::endl;
		}
		catch (Exception& e) {
			RETHROW_EXCEPTION(e, "Could not calculate B in LinThermoElectricInt" );
		}



	}



} // end namespace CoupledField
