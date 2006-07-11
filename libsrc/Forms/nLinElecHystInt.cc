#include <iostream>
#include <fstream>

#include "General/environment.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "Utils/nodestoresol.hh"
#include "nLinElecHystInt.hh"

namespace CoupledField {


  // ============
  //   calcBMat
  // ============
  void nlinElecHystInt::calcBMat( Matrix<Double> &bMat, UInt ip,
			     Matrix<Double> &ptCoord ) {

    ENTER_FCN( "nlinElecHystInt::calcBMat" );

    // Obtain info on number of element's nodes
    const UInt numNodes = ptelem->GetNumNodes();


    // Set correct size of matrix B and initialise with zeros
    bMat.Resize( dim_, numNodes );
    bMat.Init();

    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;
    ptelem->GetGlobDerivShFncAtIp( xiDx, ip, ptCoord );

    if ( subTensorType_ == FULL ) {
      // The matrix bMat can be seen as a 3 x numNodes block-vector.
      // The k-th entry of this block vector corresponds to the matrix
      // B of the BDB product evaluated at the k-th node of the finite
      // element. 
      for( UInt actNode = 0; actNode < numNodes; actNode++ ) {
	bMat[0][actNode] = xiDx[actNode][0];
	bMat[1][actNode] = xiDx[actNode][1];
	bMat[2][actNode] = xiDx[actNode][2];
      }
    }
    else  {
      // The matrix bMat can be seen as a 1 x numNodes block-vector.
      // The k-th entry of this block vector corresponds to the matrix
      // B of the ADB product evaluated at the k-th node of the finite
      // element. We assume that the first coordinate equals y and the
      // second z.
      for( UInt actNode = 0; actNode < numNodes; actNode++ ) {
        bMat[0][actNode] = xiDx[actNode][0];
        bMat[1][actNode] = xiDx[actNode][1];
      }
    }
    
  }


  // ============
  //   calcDMat
  // ============
  void nlinElecHystInt::calcDMat( Matrix<Double> &dMat ) {

    ENTER_FCN( "nlinElecHystInt::calcDMat" );

    ptMaterial->GetTensor(dMat,ELEC_PERMITTIVITY,matDataType_,subTensorType_);
    dMat *= factor_;
  }


  // ====================
  //   calcElementMatrix
  // ====================
  void nlinElecHystInt::CalcElementMatrix( Matrix<Double>& elemMat,
					   EntityIterator& ent1, 
					   EntityIterator& ent2 ) {
    ENTER_FCN( "nlinElecHystInt::CalcElementMatrix" );

    // get displacements of element
    sol_->GetElemSolution( elemPot_, ent1.GetElem()->connect );
    
    Vector<Double> LCoord, Efield;

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 ); 

    //compute electric field in the midpoint of element
    ptelem->GetCoordMidPoint(LCoord);
    EfieldOp_->CalcElemGradField( Efield,  ent1.GetElem(), LCoord, 1 );
 

    // call method of base class
    BDBInt::CalcElementMatrix( elemMat, ent1, ent2 );
  }

  // ====================
  //   set 
  // ====================
  void nlinElecHystInt::Set4Hyst(Grid* ptGrid, StdPDE* ptPDE,
				 shared_ptr<EqnMap> eqnMap,
				 shared_ptr<ResultDof> result) {
    ENTER_FCN( "nlinElecHystInt::Set4Hyst" );
    
    EfieldOp_ =  new GradientFieldOp<Double>(ptGrid, ptPDE, 
					     eqnMap, *sol_, 
					     ELEC_POTENTIAL, 
					     result, isaxi_, 
					     coordUpdate_);
  }


  // ================
  //   Constructors
  // ================


  nlinElecHystInt::nlinElecHystInt(BaseMaterial* matData, SubTensorType type,
				   bool coordUpdate ) 
    : BDBInt(matData, type, coordUpdate ) {

    ENTER_FCN( "nlinElecHystInt::nlinElecHystInt" );

    name_ = "nlinElecHystInt";
    factor_ = 1.0;

    if ( type == FULL ) {
      dim_ = 3;
    }
    else {
      dim_ = 2;
    }

    if ( type == AXI ) {
      isaxi_     = true;
    }

    isSolDependent_ = true;

  }


}
