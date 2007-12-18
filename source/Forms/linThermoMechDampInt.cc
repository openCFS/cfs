#include <iostream>
#include <fstream>

#include "linThermoMechDampInt.hh"
#include "Utils/nodestoresol.hh"


namespace CoupledField
{


  LinThermoMechDampInt::LinThermoMechDampInt(BaseMaterial* matData,
                                             BaseMaterial *matDataMech, SubTensorType type) :
    BaseForm(matData, type) {
	
    name_ = "LinThermoMechDampInt";
	
    baseType_ = DAMPING;

    pn_ =  param->Get("sequenceStep")->Get("pdeList")->
      Get("heatConduction");
	
    //Set that the integrator is solution dependen with respect to 
    //heat conduction pde
    isSolDependent_ = true;
	
    if (type == FULL ) {
      numDofsA_ = 1;
      numDofsB_ = 3;
      matDimRow_ = 1;
      matDimCol_ = 6;
	
    } else if (type == PLANE_STRAIN || type == PLANE_STRESS ) {
      numDofsA_ = 1;
      numDofsB_ = 2;
      matDimRow_ = 1;
      matDimCol_ = 3;
    } else if (type == AXI ) {
      numDofsA_ = 1;
      numDofsB_ = 2;
      matDimRow_ = 1;
      matDimCol_ = 4;
      isaxi_ = true;
    }
	
    try {
      // get the elastic tensor coefficients
      cMatrix_.Init();
	
      matDataMech->GetTensor(cMatrix_,MECH_STIFFNESS_TENSOR,REAL,subTensorType_);
      //ptMaterial->GetTensor(dMat,MECH_STIFFNESS_TENSOR,matDataType_,subTensorType_);
      
    }
    catch (Exception& e) {
      RETHROW_EXCEPTION(e, "Could not get the MECH_STIFFNESS_TENSOR"
                        <<" in LinThermoMechDampInt" );
    }
    
  }
  

 
  LinThermoMechDampInt::~LinThermoMechDampInt()
  {
  }


  // Query material type for  D  tensor.
  MaterialType LinThermoMechDampInt::getDMaterialType() {return THERMAL_EXPANSION_TENSOR;}


  // =====================
  //   CalcElementMatrix
  // =====================
  void LinThermoMechDampInt::CalcElementMatrix(Matrix<Double>& elemMat,
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
      if(param->Has("transient") )
        sol2_->GetElemSolution( teta_, ent2);
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
          //std::cerr << "dMat = \n" << dMat << std::endl;
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
          (*error) << "LinThermoMechDampInt::CalcElementMatrix: Encountered "
                   << "negative Jacobian determinant!";
          Error( __FILE__, __LINE__ );
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
        dbMat.Resize( dMat.GetSizeRow(), bMat.GetSizeCol() );
        dMat.Mult( bMat, dbMat );
	
        // We now compute A * D * B and scale it by the determinant
        // of the Jacobian and the weight of the current integration
        // point. The result is added to the element matrix
        for ( UInt i = 0; i < aMat.GetSizeRow(); i++ ) {
          for ( UInt j = 0; j < dbMat.GetSizeCol(); j++ ) {
	
            // Compute entry (i,j) of A * D * B
            aux = 0.0;
            for ( UInt k = 0; k < aMat.GetSizeCol(); k++ ) {
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
      RETHROW_EXCEPTION(e, "Could not get the CalcElementMatrix"
                        <<" in LinThermoMechDampInt" );
    }

  }
    
    
    
  //  Compute matrix  A  at given integration point.
  void LinThermoMechDampInt::calcAMat(Matrix<Double> &aMat, UInt ip,
                                      const Matrix<Double> &ptCoord) {

	
    try {


      Double refTemp = pn_->Has("referenceTemperature") ?
        pn_->Get("referenceTemperature")->AsDouble() : 0.0;

      
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
      for(UInt i=0; i < lN_tetaXN_teta.GetSizeRow(); i++ ) {
        for(UInt j=0; j < lN_tetaXN_teta.GetSizeCol(); j++ ) {
          lN_tetaXN_teta[i][j] += lN_teta[i]*lN_teta[j];;
        }
      }
	
      //std::cerr << "N_tetaT * N_teta = \n" << lN_tetaXN_teta << std::endl;
	
      // Set correct size of matrix A and initialize with zeros
      aMat.Resize(numFncs, numDofsA_);
      aMat.Init();
	
      
      //multiply the last matrix by the last temperature
      for(UInt i=0; i < lN_tetaXN_teta.GetSizeCol(); i++ ) {
        //        for(UInt j=0; j < teta_.GetSize(); j++ ) {
        for(UInt j=0; j < numDofsA_; j++ ) {
          aMat[i][0] += lN_tetaXN_teta[i][j]*(teta_[j]+refTemp);
        }
      }
      
    }
    catch (Exception& e) {
      RETHROW_EXCEPTION(e, "Could not get the A"
                        <<" in LinThermoMechDampInt" );
    }
    
  }

  //  Compute matrix  D  at given integration point.
  //  return the transposed of D matrix
  void LinThermoMechDampInt::calcDMat(Matrix<Double> &dMat) {

	
    try {
      // get the thermal expansion coefficient tensor
      Matrix<Double> aux;
      ptMaterial->GetTensor(aux,THERMAL_EXPANSION_TENSOR,matDataType_,subTensorType_);
	
      //std::cerr << "THERMAL_EXPANSION_TENSOR = \n" << aux << std::endl;
	
      dMat.Resize(matDimRow_, matDimCol_);
      dMat.Init();
	
      // transforming p_alphaMatrix to Voigt Notation
      Matrix<Double> l_alphaMatrix(matDimCol_,matDimRow_);
	
      if ( subTensorType_ == FULL ) {
	
        l_alphaMatrix(0,0) = aux(0,0);
        l_alphaMatrix(1,0) = aux(1,1);
        l_alphaMatrix(2,0) = aux(2,2);
        l_alphaMatrix(3,0) = 2 * aux(1,2);
        l_alphaMatrix(4,0) = 2 * aux(0,2);
        l_alphaMatrix(5,0) = 2 * aux(0,1);
	
      }
      // on the plane yz
      else if ( subTensorType_ == PLANE_STRAIN ) {
	
        l_alphaMatrix(0,0) = aux(0,0);
        l_alphaMatrix(1,0) = aux(1,1);
        l_alphaMatrix(2,0) = 2 * aux(0,1);
      }
      else if ( subTensorType_ == AXI ) {
        l_alphaMatrix(0,0) = aux(0,0);
        l_alphaMatrix(1,0) = aux(1,1);
        l_alphaMatrix(2,0) = 2 * aux(0,1);
        l_alphaMatrix(3,0) = aux(2,2);
      }
      else if ( subTensorType_ == PLANE_STRESS ) {
	
        EXCEPTION ("LinThermoMechInt: PLANE_STRESS Type not supported" );
      }
			
      //std::cerr << "l_alphaMatrix = \n" << l_alphaMatrix << std::endl;
	
      //reset aux
      aux.Resize(matDimCol_,matDimRow_);
      aux.Init();
	
      cMatrix_.Mult(l_alphaMatrix, aux);
	
      //std::cerr << "dMat = \n" << dMat << std::endl;
	
      //but we need the transposed
	
      aux.Transpose(dMat);
	
    }
    catch (Exception& e) {
      RETHROW_EXCEPTION(e, "Could not get the D"
                        <<" in LinThermoMechDampInt" );
    }
    
    
  }
  
  
  //  Compute matrix  B  at given integration point.
  void LinThermoMechDampInt::calcBMat(Matrix<Double> &bMat, UInt ip,
                                      const Matrix<Double> &ptCoord) {

	
    try {
	
      // Obtain info on problem sizes
      UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
      //std::cout << "numFncs" << numFncs << std::endl;
	
      // Set type of ansatz function , but do not recalculate
      // integration points
      ptelem->SetAnsatzFct( ansatzFct1_, false );
      //    const UInt nDofMech = 3;
      //    const UInt nRowsD   = 6;
	
      // Set correct size of matrix B and initialize it with zeros
      bMat.Resize( matDimCol_, numFncs * numDofsB_ );
      bMat.Init();
	
      // Get derivatives of local shape functions with respect to global coords
      // (format: numFncs x spaceDim)
      Matrix<Double> xiDx;
      ptelem->GetGlobDerivShFncAtIp( xiDx, ip, ptCoord, it1_.GetElem() );
	
      //std::cerr << "xiDx = \n" << xiDx << std::endl;
	
	
      // The matrix aMat can be seen as a numFncs x 1 block-vector.
      // The k-th entry of this block vector corresponds to the matrix
      // A of the ADB product evaluated at the k-th node of the finite
      // element. We assemble aMat in a top-down fashion.
      UInt actInd = 0;
      UInt actNode;
	
      if ( subTensorType_ == FULL ) {
        for( actNode = 0; actNode < numFncs; actNode++ ) {
	
          // 1st row of sub-matrix A(actNode)
          bMat[0][actInd] = xiDx[actNode][0];
          bMat[4][actInd] = xiDx[actNode][2];
          bMat[5][actInd] = xiDx[actNode][1];
          actInd++;
	
          // 2nd row of sub-matrix A(actNode)
          bMat[1][actInd] = xiDx[actNode][1];
          bMat[3][actInd] = xiDx[actNode][2];
          bMat[5][actInd] = xiDx[actNode][0];
          actInd++;
	
          // 3rd row of sub-matrix A(actNode)
          bMat[2][actInd] = xiDx[actNode][2];
          bMat[3][actInd] = xiDx[actNode][1];
          bMat[4][actInd] = xiDx[actNode][0];
          actInd++;
        }
      }
      else if ( subTensorType_ == AXI ) { // to test
	
        //EXCEPTION ("calcBMat: AXI Type not supported" );
        // we assume that the first coordinate equals r and 
        // the second z.
	
        UInt j;
        Double coordAtIp;
        Vector<Double> ShpFncAtIp;
        ptelem->GetShFncAtIp( ShpFncAtIp, ip, it1_.GetElem() );
	
        for( actNode = 0; actNode < numFncs; actNode++ ) {
	
          // 1st row of sub-matrix A(actNode)
          bMat[0][actInd] = xiDx[actNode][0]; // dN/dr
          bMat[2][actInd] = xiDx[actNode][1]; // dN/dz
	
          // For the entry 1/r things are more complicated
          coordAtIp = 0.0;
          for( j = 0; j < numFncs; j++ ) {
            coordAtIp += ptCoord[0][j] * ShpFncAtIp[j];
          }
          bMat[3][actInd] = ShpFncAtIp[actNode] / coordAtIp; // 1/r
          actInd++;
	
          // 2nd row of sub-matrix A(actNode)
          bMat[1][actInd] = xiDx[actNode][1]; // dN/dz
          bMat[2][actInd] = xiDx[actNode][0]; // dN/dr
	
          actInd++;
        }
      }
	
      else if ( subTensorType_ == PLANE_STRAIN ) {
        // we assume that the first coordinate equals y and the second z.
	
        for( actNode = 0; actNode < numFncs; actNode++ ) {
	
          // 1st row of sub-matrix A(actNode)
          bMat[0][actInd] = xiDx[actNode][0]; // dN/dx
          bMat[2][actInd] = xiDx[actNode][1]; // dN/dy
          actInd++;
	
          // 2nd row of sub-matrix A(actNode)
          bMat[1][actInd] = xiDx[actNode][1]; // dN/dy
          bMat[2][actInd] = xiDx[actNode][0]; // dN/dx
          actInd++;
        }
      }
	
      else if ( subTensorType_ == PLANE_STRESS) {
        EXCEPTION ("calcBMat: PLANE_STRESS Type not supported" );
      }
      
    }
    catch (Exception& e) {
      RETHROW_EXCEPTION(e, "Could not get the B"
                        <<" in LinThermoMechDampInt" );
    }
    
  }




} // end namespace CoupledField
