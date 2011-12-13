#include <fstream>
#include <string>

#include "Domain/entityList.hh"
#include "Elements/basefe.hh"
#include "Forms/adbInt.hh"
#include "General/environment.hh"
#include "General/exception.hh"
#include "MatVec/vector.hh"
#include "Materials/baseMaterial.hh"
#include "linThermoMechInt.hh"

//#include "Domain/domain.hh"
//#include "Domain/grid.hh"

namespace CoupledField
{


  LinThermoMechInt::LinThermoMechInt(BaseMaterial* matData, BaseMaterial *matDataMech, SubTensorType type) :
    ADBInt(matData,type) 
  {

    name_ = "LinThermoMechInt";

    // say to the parent adb that the stiffness k of
    // the coupling is thermo-mech
    //isThermoMech_ = true;


    if ( type == FULL ) {
      numDofsA_  = 3;
      numDofsB_  = 1;
      matDimRow_ = 6;
      matDimCol_ = 1;

    }
    else if ( type == PLANE_STRAIN || type == PLANE_STRESS ) {
      numDofsA_  = 2;
      numDofsB_  = 1;
      matDimRow_ = 3;
      matDimCol_ = 1;
    }
    else if ( type == AXI ) {
      numDofsA_  = 2;
      numDofsB_  = 1;
      matDimRow_ = 4;
      matDimCol_ = 1;
      isaxi_     = true;
      //EXCEPTION ("LinThermoMechInt: AXI Type not supported" );
    }

    try{

      // get the elastic tensor coefficients
      cMatrix_.Init();

      matDataMech->GetTensor(cMatrix_,MECH_STIFFNESS_TENSOR,Global::REAL,subTensorType_);
    
      //std::cerr << "\nMECH_STIFFNESS_TENSOR = \n" << cMatrix_ << std::endl;
    }
    catch (Exception& e) {
      RETHROW_EXCEPTION(e, "Could not get MECH_STIFFNESS_TENSOR in "
                        <<"LinThermoMechInt" );
    }

  }


 
  LinThermoMechInt::~LinThermoMechInt()
  {
  }


  // Query material type for  D  tensor.
  MaterialType LinThermoMechInt::getDMaterialType() {return THERMAL_EXPANSION_TENSOR;}


	
  //----------------------------------------------------------------
  //  Compute matrix  A  at given integration point.
  //----------------------------------------------------------------
  void LinThermoMechInt::calcAMat(Matrix<Double> &aMat, UInt ip,
                                  const Matrix<Double> &ptCoord) {

	
    try {
      //std::cout << "\nworking on the thermomechanics integral ..... calcAMat" << std::endl; 	
	
      // Obtain info on problem sizes
      UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
      //std::cout << "numFncs" << numFncs << std::endl;
	
      // Set type of ansatz function , but do not recalculate
      // integration points
      ptelem->SetAnsatzFct( ansatzFct1_, false );
      //    const UInt nDofMech = 3;
      //    const UInt nRowsD   = 6;
	
      // Set correct size of matrix A and initialize it with zeros
      aMat.Resize( numFncs * numDofsA_, matDimRow_ );
      aMat.Init();
	
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
          aMat[actInd][0] = xiDx[actNode][0];
          aMat[actInd][4] = xiDx[actNode][2];
          aMat[actInd][5] = xiDx[actNode][1];
          actInd++;
	
          // 2nd row of sub-matrix A(actNode)
          aMat[actInd][1] = xiDx[actNode][1];
          aMat[actInd][3] = xiDx[actNode][2];
          aMat[actInd][5] = xiDx[actNode][0];
          actInd++;
	
          // 3rd row of sub-matrix A(actNode)
          aMat[actInd][2] = xiDx[actNode][2];
          aMat[actInd][3] = xiDx[actNode][1];
          aMat[actInd][4] = xiDx[actNode][0];
          actInd++;
        }
      }
      else if ( subTensorType_ == AXI ) { 
	
        // we assume that the first coordinate equals r and 
        // the second z.
	
        UInt j;
        Double coordAtIp;
        Vector<Double> ShpFncAtIp;
        ptelem->GetShFncAtIp( ShpFncAtIp, ip, it1_.GetElem() );
	
        for( actNode = 0; actNode < numFncs; actNode++ ) {
	
          // 1st row of sub-matrix A(actNode)
          aMat[actInd][0] = xiDx[actNode][0]; // dN/dr
          aMat[actInd][2] = xiDx[actNode][1]; // dN/dz
	
          // For the entry 1/r things are more complicated
          coordAtIp = 0.0;
          for( j = 0; j < numFncs; j++ ) {
            coordAtIp += ptCoord[0][j] * ShpFncAtIp[j];
          }
          aMat[actInd][3] = ShpFncAtIp[actNode] / coordAtIp; // 1/r
          actInd++;
	
          // 2nd row of sub-matrix A(actNode)
          aMat[actInd][1] = xiDx[actNode][1]; // dN/dz
          aMat[actInd][2] = xiDx[actNode][0]; // dN/dr
	
          actInd++;
        }
      }
	
      else if ( subTensorType_ == PLANE_STRAIN ) {
        // we assume that the first coordinate equals y and the second z.
	
        for( actNode = 0; actNode < numFncs; actNode++ ) {
	
          // 1st row of sub-matrix A(actNode)
          aMat[actInd][0] = xiDx[actNode][0]; // dN/dz
          aMat[actInd][2] = xiDx[actNode][1]; // dN/dy
          actInd++;
	
          // 2nd row of sub-matrix A(actNode)
          aMat[actInd][1] = xiDx[actNode][1]; // dN/dz
          aMat[actInd][2] = xiDx[actNode][0]; // dN/dx
          actInd++;
        }
      }
	
      else if ( subTensorType_ == PLANE_STRESS) {
        EXCEPTION ("calcAMat: PLANE_STRESS Type not supported" );
      }
	
      //std::cerr << "aMat = \n" << aMat << std::endl;
    }
    catch (Exception& e) {
      RETHROW_EXCEPTION(e, "Could not calculate A in LinThermoMechInt" );
    }

  }


  //----------------------------------------------------------------
  //  Compute matrix  D  at given material parameters.
  //	D is the thermal-stress coefficient tensor
  //	calculated for anisotropic materials
  //	lambda(i,j) = c(i,j,k,l) * alpha(k,l)
  //----------------------------------------------------------------
  void LinThermoMechInt::calcDMat(Matrix<Double> &dMat) 
  {
    try {
      //std::cout << "\nworking on the thermomechanics integral ..... calcDMat" << std::endl; 	

      // get the thermal expansion coefficient tensor
      Matrix<Double> p_alphaMatrix;
      ptMaterial->GetTensor(p_alphaMatrix,THERMAL_EXPANSION_TENSOR,matDataType_,subTensorType_);

      //std::cerr << "THERMAL_EXPANSION_TENSOR = \n" << p_alphaMatrix << std::endl;

      //std::cerr << "matDimRow_ = " << matDimRow_ << "matDimCol_ ="<<matDimCol_<< std::endl;

      dMat.Resize(matDimRow_, matDimCol_);
      dMat.Init();

      // transforming p_alphaMatrix to Voigt Notation
      Matrix<Double> l_alphaMatrix(matDimRow_, matDimCol_);

      if ( subTensorType_ == FULL ) {

        l_alphaMatrix(0,0) = p_alphaMatrix(0,0);
        l_alphaMatrix(1,0) = p_alphaMatrix(1,1);
        l_alphaMatrix(2,0) = p_alphaMatrix(2,2);
        l_alphaMatrix(3,0) = 2 * p_alphaMatrix(1,2);
        l_alphaMatrix(4,0) = 2 * p_alphaMatrix(0,2);
        l_alphaMatrix(5,0) = 2 * p_alphaMatrix(0,1);

      }
      // on the plane yz
      else if ( subTensorType_ == PLANE_STRAIN ) {

        l_alphaMatrix(0,0) = p_alphaMatrix(0,0);
        l_alphaMatrix(1,0) = p_alphaMatrix(1,1);
        l_alphaMatrix(2,0) = 2 * p_alphaMatrix(0,1);

      }
      else if ( subTensorType_ == PLANE_STRESS ) {

        EXCEPTION ("LinThermoMechInt: PLANE_STRESS Type not supported" );
      }
      else if ( subTensorType_ == AXI ) {

        l_alphaMatrix(0,0) = p_alphaMatrix(0,0);
        l_alphaMatrix(1,0) = p_alphaMatrix(1,1);
        l_alphaMatrix(2,0) = 2 * p_alphaMatrix(0,1);
        l_alphaMatrix(3,0) = p_alphaMatrix(2,2);

      }

      //std::cerr << "l_alphaMatrix = \n" << l_alphaMatrix << std::endl;

      cMatrix_.Mult(l_alphaMatrix, dMat);

      //std::cerr << "dMat = \n" << dMat << std::endl;
    }
    catch (Exception& e) {
      RETHROW_EXCEPTION(e, "Could not calculate D in LinThermoMechInt" );
    }
  }
    
    
  //----------------------------------------------------------------
  //  Compute matrix  B  at given integration point.
  //----------------------------------------------------------------
  void LinThermoMechInt::CalcBMat(Matrix<Double> &bMat, UInt ip,
                                  const Matrix<Double> &ptCoord) {

	
    try {
      //std::cout << "\nworking on the thermomechanics integral ..... calcBMat" << std::endl; 	
	
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
	
      for(UInt i=0; i < bMat.GetNumCols(); i++ ) {
        bMat[0][i] = lN_teta[i];
      }
	
      //std::cerr << "bMat = \n" << bMat << std::endl;
    }
    catch (Exception& e) {
      RETHROW_EXCEPTION(e, "Could not calculate B in LinThermoMechInt" );
    }

  }




} // end namespace CoupledField
