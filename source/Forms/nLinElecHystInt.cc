// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "General/environment.hh"
#include "Utils/nodestoresol.hh"
#include "nLinElecHystInt.hh"

namespace CoupledField {


  // ============
  //   calcBMat
  // ============
  void nlinElecHystInt::CalcBMat( Matrix<Double> &bMat, UInt ip,
                                  const Matrix<Double> &ptCoord ) {


    // Obtain info on number of element's nodes
    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
   
    // Set correct size of matrix B and initialise with zeros
    bMat.Resize( dim_, numFncs );
    bMat.Init();

    // Get derivatives of local shape functions with respect to global
    // coords (format: numFncs x spaceDim)
    Matrix<Double> xiDx;
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
    
  }


  // ============
  //   calcDMat
  // ============
  void nlinElecHystInt::calcDMat( Matrix<Double> &dMat ) {


    ptMaterial->GetTensor(dMat,ELEC_PERMITTIVITY,matDataType_,subTensorType_);
    dMat.Init();
    for ( UInt i=0; i<dMat.GetNumRows(); i++ ) {
      dMat[i][i] = diffEpsVal_;
    }

    dMat *= mParser_->Eval( mHandle_ );

    //    std::cout << "dMat: \n" << dMat << std::endl;
  }


  // ====================
  //   calcElementMatrix
  // ====================
  void nlinElecHystInt::CalcElementMatrix( Matrix<Double>& elemMat,
					   EntityIterator& ent1, 
					   EntityIterator& ent2 ) {

    // Get pointer to reference element and its coordinates
    ExtractElemInfo( ent1 );

    // get scalar electric potential
    sol_->GetElemSolution( elemPot_, ent1 );

    Vector<Double> LCoord, Efield;

    //compute electric field in the midpoint of element
    ptelem->GetCoordMidPoint(LCoord);
    EfieldOp_->CalcElemGradField( Efield,  ent1, LCoord, 1 );

    Double Ecomp = Efield[dirP_]; //.NormL2(); 
 
    UInt nrEl = ent1.GetElem()->elemNum;
    diffEpsVal_ = ptMaterial->ComputeScalarDiffVal( nrEl, Ecomp );

    if (  diffEpsVal_ < 0.0 ) 
      EXCEPTION("Negative effective permittivity");

    // call method of base class
    BDBInt::CalcElementMatrix( elemMat, ent1, ent2 );
  }

  // ====================
  //   set 
  // ====================
  void nlinElecHystInt::Set4Hyst(Grid* ptGrid, StdPDE* ptPDE,
				 shared_ptr<EqnMap> eqnMap,
				 shared_ptr<ResultInfo> result) {
    
    EfieldOp_ =  new GradientFieldOp<Double>(ptGrid, ptPDE, 
					     eqnMap, *sol_, result->fctType, isaxi_, coordUpdate_);
    std::string str;
    ptMaterial->GetScalar(str, P_DIRECTION);
    Directions dir;
    String2Enum(str,dir);
    dirP_ = dir;
  }


  void nlinElecHystInt::SetFactor( const std::string& factor ) {
    mParser_->SetExpr( mHandle_, factor );
  }
  

  // ================
  //   Constructors
  // ================


  nlinElecHystInt::nlinElecHystInt(BaseMaterial* matData, SubTensorType type,
				   bool coordUpdate ) 
    : BDBInt(matData, type, coordUpdate ) {


    name_ = "nlinElecHystInt";

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
