#include <iostream>
#include <fstream>

#include "linThermoElecDampInt.hh"

#include "Domain/domain.hh"
#include "Domain/grid.hh"
//#include "Utils/result.hh"
#include "Utils/nodestoresol.hh"


namespace CoupledField
{

  LinThermoElecDampInt::LinThermoElecDampInt(BaseMaterial* matData,
		SubTensorType type) :
	BaseForm(matData, type) 
  {
    
    baseType_ = DAMPING;

    name_ = "LinThermoElecDampInt";

    pn_ =  param->Get("sequenceStep")->Get("pdeList")->
      Get("heatConduction");

    //Set that the integrator is solution dependen with respect to 
    //heat conduction pde
    isSolDependent_ = true;

	//std::cout << "\nworking on the damping thermo-electrostatics integral ..... " << std::endl; 	

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
        numDofsA_  = 1;
        numDofsB_  = 1;
        
        matDimRow_ = 1;
        matDimCol_ = 2;
        isaxi_     = true;
		
    }

  }

 
	LinThermoElecDampInt::~LinThermoElecDampInt()
	{
	}

    // Query material type for  D  tensor.
	MaterialType LinThermoElecDampInt::getDMaterialType() {return PYROCOEFFICIENT_TENSOR;}
	

	
	
	  // =====================
	  //   CalcElementMatrix
	  // =====================
	void LinThermoElecDampInt::CalcElementMatrix(Matrix<Double>& elemMat,
		EntityIterator& ent1, EntityIterator& ent2) {

	
		try {
			// Extract pointer to reference element and get coordinates
			ExtractElemInfo( ent1 );
	
			ptelem->SetAnsatzFct( ansatzFct1_ );
			UInt numFncs1 = ptelem->GetNumFncs( ansatzFct1_ );
			UInt numFncs2 = ptelem->GetNumFncs( ansatzFct2_ );
	
			const UInt nrIntPts = ptelem->GetNumIntPoints();
			const Vector<Double> & intWeights = ptelem->GetIntWeights();
			double jacDet;
	
			Matrix<Double> aMat;
			Matrix<Double> bMat;
			Matrix<Double> dMat;
			Matrix<Double> dbMat;
			Double aux;
	
			elemMat.Resize( numFncs1 * getNumDofsA(), numFncs2 * getNumDofsB() );
			elemMat.Init();
	
			// Get the temperature nodal values from the previous time step

            std::string analysis;
            if(param->Get("sequenceStep")->Get("analysis")->Has("transient") ){
                 sol2_->GetElemSolution( teta_, ent2);
            }
            else{
              teta_.Resize(numDofsA_);
              teta_.Init();
            }
            	
			// **************************************************
	
			// Setup material matrix once and for all
			//  Material matrix independent of integration point
			calcDMat( dMat );
	
			// Loop over all integration points
			for ( UInt actIntPt = 1; actIntPt <= nrIntPts; actIntPt++ ) {
	
				// Check if material has to be rotated for each integration point
				if( ptMaterial->GetCoordSys() != NULL ) {
					// Get global coordinates
					Vector<Double> * intPoints = ptelem->GetIntPoints();
					Vector<Double> globIntPoint;
	
					ptelem->Local2GlobalCoord(globIntPoint, intPoints[actIntPt-1],
							ptCoord_, ent1.GetElem() );
					ptMaterial->RotateTensorByPointCoord( globIntPoint, getDMaterialType() );
					calcDMat( dMat );
                    //					std::cerr << "dMat = \n" << dMat << std::endl;
				}
	
				//std::cerr << "*** Calculating A ****\n";
				// Setup the A matrix for current integration point
				calcAMat( aMat, actIntPt, ptCoord_ );
				//std::cerr << "aMat = \n" << aMat << std::endl;
	
				//std::cerr << "*** Calculating B ****\n";
				// Setup the B matrix for current integration point
				calcBMat( bMat, actIntPt, ptCoord_ );
				//std::cerr << "bMat = \n" << bMat << std::endl;
	
				// Compute Jacobian for integration point
				jacDet = ptelem->CalcJacobianDetAtIp( actIntPt, ptCoord_,
						ent1.GetElem() );
	
				// Perform a safety check
				if ( jacDet < 0.0 ) {
				  EXCEPTION( "ADBInt::CalcElementMatrix: Encountered "
					<< "negative Jacobian determinant!" );
				}
	
				// Special things must be done in the axi-symmetric case
				// We need to additionally scale with 2 pi r.
				//
				// NOTE: We assume here that computation is in they-z plane
				//       with z being the axis of symmetry and that y is
				//       represented by the first co-ordinate, thus
				//       2 pi r = "2 pi x"
				if ( isaxi_ ) {
					Vector<Double> ShpFncAtIp;
					ptelem->GetShFncAtIp( ShpFncAtIp, actIntPt, ent1.GetElem() );
					Double aux = 0.0;
	
					for ( UInt i = 0; i < numFncs1; i++ ) {
						aux += ptCoord_[0][i] * ShpFncAtIp[i];
					}
	
					jacDet *= 2.0 * PI * aux;
				}
	
				// Compute the matrix product D * B and store as intermediate matrix
				dbMat.Resize( dMat.GetNumRows(), bMat.GetNumCols() );
				dMat.Mult( bMat, dbMat );
	
				// We now compute A * D * B and scale it by the determinant
				// of the Jacobian and the weight of the current integration
				// point. The result is added to the element matrix
				for ( UInt i = 0; i < aMat.GetNumRows(); i++ ) {
					for ( UInt j = 0; j < dbMat.GetNumCols(); j++ ) {
	
						// Compute entry (i,j) of A * D * B
						aux = 0.0;
						for ( UInt k = 0; k < aMat.GetNumCols(); k++ ) {
							aux += aMat[i][k] * dbMat[k][j];
						}
	
						// Scale result and add to corresponding entry
						// of element matrix
						elemMat[i][j] += aux * jacDet * intWeights[actIntPt-1];
					}
				}
			}
	
			//std::cerr << "elemMat = \n" << elemMat << std::endl;
		}
		catch (Exception& e) {
			RETHROW_EXCEPTION(e, "Could not calculate CalcElementMatrix"
					<<" in LinThermoElectricDampInt" );
		}
	}
	
	

	//  Compute matrix  A  at given integration point.
	void LinThermoElecDampInt::calcAMat(Matrix<Double> &aMat, UInt ip,
		const Matrix<Double> &ptCoord) {

		try {

          Double refTemp = pn_->Has("referenceTemperature") ?
            pn_->Get("referenceTemperature")->As<Double>() : 0.0;

	
			// Set type of ansatz function , but do not recalculate
			// integration points (=false)
			ptelem->SetAnsatzFct( ansatzFct2_, false );
	
			// Obtain info on number of element's nodes
			UInt numFncs = ptelem->GetNumFncs( ansatzFct2_ );
	
			// Get local shape functions of thermal field (get sure they're from thermal field)
			// (format: numFncs x 1)
			Vector< Double > lN_teta;
			ptelem->GetShFncAtIp( lN_teta, ip, it1_.GetElem() );
	
			//std::cout << "lN_teta=" << lN_teta.Serialize () << std::endl;
	
			//multiply the transpose of the shape function matrix by the shape functions matrix
			Matrix<Double> lN_tetaXN_teta(numFncs, numFncs);
			lN_tetaXN_teta.Init();
			for(UInt i=0; i < lN_tetaXN_teta.GetNumRows(); i++ ) {
				for(UInt j=0; j < lN_tetaXN_teta.GetNumCols(); j++ ) {
					lN_tetaXN_teta[i][j] += lN_teta[i]*lN_teta[j];;
				}
			}
	
			//std::cerr << "N_tetaT * N_teta = \n" << lN_tetaXN_teta << std::endl;
	
			// Set correct size of matrix A and initialize with zeros
			aMat.Resize(numFncs, 1);
			aMat.Init();
	
			//multiply the last matrix by the last temperature
			for(UInt i=0; i < lN_tetaXN_teta.GetNumCols(); i++ ) {
              for(UInt j=0; j <teta_.GetSize(); j++ ) {
                //				for(UInt j=0; j < teta_.GetSize(); j++ ) {
                aMat[i][0] += lN_tetaXN_teta[i][j]*(teta_[j]+refTemp);
              }
			}
            
            //			std::cerr << "LinThermoElectricDampInt: aMat = \n" << aMat << std::endl;
		}
		catch (Exception& e) {
			RETHROW_EXCEPTION(e, "Could not calculate A"
					<<" in LinThermoElectricDampInt" );
		}

	}

	//  Compute matrix  D  at given integration point.
	void LinThermoElecDampInt::calcDMat(Matrix<Double> &dMat) {

	
		try {
			//std::cout << "\nworking on the thermo-electrostatics integral ..... calcDMat" << std::endl; 	
	
			// get the thermal expansion coefficient tensor
			dMat.Resize(matDimRow_, matDimCol_);
			ptMaterial->GetTensor(dMat,PYROCOEFFICIENT_TENSOR,matDataType_,subTensorType_);
	
			Matrix<Double> p_auxMat;
			if(isaxi_==false)
				p_auxMat.Resize(matDimRow_, matDimRow_);
			else
				p_auxMat.Resize(matDimCol_, matDimCol_);
			dMat.GetDiagInMatrix(p_auxMat);
	
			dMat.Init();
			p_auxMat.Transpose(dMat);
	
            //			std::cerr << "LinThermoElectricDampInt: dMat = \n" << dMat << std::endl;	
		}
		catch (Exception& e) {
			RETHROW_EXCEPTION(e, "Could not calculate D"
					<<" in LinThermoElectricDampInt" );
		}
	}

	
	//  Compute matrix  B  at given integration point.
	void LinThermoElecDampInt::calcBMat(Matrix<Double> &bMat, UInt ip,
		const Matrix<Double> &ptCoord) {

	
		try {
	
			// Obtain info on problem sizes
			UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
			//std::cout << "numFncs" << numFncs << std::endl;
	
			// Set type of ansatz function , but do not recalculate
			// integration points
			ptelem->SetAnsatzFct( ansatzFct1_, false );
	
			// Get derivatives of local shape functions with respect to global coords
			// (format: numFncs x spaceDim)
			Matrix<Double> auxMat(numFncs, matDimRow_);
			ptelem->GetGlobDerivShFncAtIp( auxMat, ip, ptCoord, it1_.GetElem() );
	
			// Set correct size of matrix B and initialize it with zeros
			bMat.Resize( matDimRow_,numFncs );
			bMat.Init();
	
			auxMat.Transpose(bMat);
	 
			//std::cerr << "LinThermoElectricDampInt: bMat = \n" << bMat << std::endl;
		}
		catch (Exception& e) {
			RETHROW_EXCEPTION(e, "Could not calculate B"
					<<" in LinThermoElectricDampInt" );
		}

	}

} // end namespace CoupledField
