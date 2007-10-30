#include <iostream>
#include <fstream>

#include "General/environment.hh"
#include "linHeatCondInt.hh"

namespace CoupledField {



	// ================
	//   Constructors
	// ================
	
	
	linHeatCondInt::linHeatCondInt(BaseMaterial* matData, SubTensorType type,
	                       bool coordUpdate ) 
	  : BDBInt(matData, type, coordUpdate ) {
	
	
	  name_ = "linHeatCondInt";
	  if ( type == FULL ) {
	    dim_ = 3;
	  }
	  else {
	    dim_ = 2;
	  }
	
	  if ( type == AXI ) {
	    isaxi_     = true;
	  }
	}

  // ============
  //   calcBMat
  // ============
  void linHeatCondInt::calcBMat( Matrix<Double> &bMat, UInt ip,
			     Matrix<Double> &ptCoord ) {

    // Obtain info on number of elements' funtions
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    
    // Set correct size of matrix B and initialise with zeros
    bMat.Resize( dim_, numFncs );
    bMat.Init();

    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;
    ptelem->SetAnsatzFct( ansatzFct1_ );
    ptelem->GetGlobDerivShFncAtIp( xiDx, ip, ptCoord, it1_.GetElem() );

    if ( subTensorType_ == FULL ) {
      // The matrix bMat can be seen as a 3 x numFncs block-vector.
      // The k-th entry of this block vector corresponds to the matrix
      // B of the BDB product evaluated at the k-th node of the finite
      // element. 
      for( UInt actNode = 0; actNode < numFncs; actNode++ ) {
	bMat[0][actNode] = xiDx[actNode][0];
	bMat[1][actNode] = xiDx[actNode][1];
	bMat[2][actNode] = xiDx[actNode][2];
      }
    }
    else  {
      // The matrix bMat can be seen as a 1 x numFncs block-vector.
      // The k-th entry of this block vector corresponds to the matrix
      // B of the ADB product evaluated at the k-th node of the finite
      // element. We assume that the first coordinate equals y and the
      // second z.
      for( UInt actNode = 0; actNode < numFncs; actNode++ ) {
        bMat[0][actNode] = xiDx[actNode][0];
        bMat[1][actNode] = xiDx[actNode][1];
      }
    }

    
 	Matrix<Double> auxbMat;
 	auxbMat.Init();
 	bMat.Transpose (auxbMat);
     //std::cerr << "linHeatCondInt::bMat transpose = \n" << auxbMat << std::endl;
  }


  // ============
  //   calcDMat
  // ============
  void linHeatCondInt::calcDMat( Matrix<Double> &dMat ) {

    ptMaterial->GetTensor(dMat,HEAT_CONDUCTIVITY_TENSOR,matDataType_,subTensorType_);
    
    dMat *= mParser_->Eval( mHandle_ );
    
    //std::cerr << "linHeatCondInt: dMat = \n" << dMat << std::endl;
  }


  void linHeatCondInt::SetFactor( const std::string& factor ) {
    
    //mParser_->SetExpr( mHandle_, factor );
  }
  




}
