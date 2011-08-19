// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "General/environment.hh"
#include "linElecInt.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "Elements/H1Elems.hh"

//DECLARE_LOG(forms)

namespace CoupledField {


  // ============
  //   calcBMat
  // ============

void linElecInt::calcBMat(Matrix<Double>& bMat, 
                             LocPointMapped& lp, BaseFE* ptFe )  {
  
  UInt numFncs = ptFe->GetNumFncs();
  // Set correct size of matrix B and initialise with zeros
  bMat.Resize( dim_, numFncs );
  bMat.Init();

  // Get derivatives of local shape functions with respect to global
  // coords (format: nrNodes x spaceDim)
  Matrix<Double> xiDx;
  FeH1 *feH1 = (dynamic_cast<FeH1*>(ptFe));
  feH1->GetGlobDerivShFnc( xiDx, lp, lp.shapeMap->GetElem() , 1 );
  bMat = Transpose (xiDx);
//  
//  if ( subTensorType_ == FULL ) {
//    // The matrix bMat can be seen as a 3 x numFncs block-vector.
//    // The k-th entry of this block vector corresponds to the matrix
//    // B of the BDB product evaluated at the k-th node of the finite
//    // element. 
//    for( UInt actNode = 0; actNode < numFncs; actNode++ ) {
//      bMat[0][actNode] = xiDx[actNode][0];
//      bMat[1][actNode] = xiDx[actNode][1];
//      bMat[2][actNode] = xiDx[actNode][2];
//    }
//  }
//  else  {
//    // The matrix bMat can be seen as a 1 x numFncs block-vector.
//    // The k-th entry of this block vector corresponds to the matrix
//    // B of the ADB product evaluated at the k-th node of the finite
//    // element. We assume that the first coordinate equals y and the
//    // second z.
//    for( UInt actNode = 0; actNode < numFncs; actNode++ ) {
//      bMat[0][actNode] = xiDx[actNode][0];
//      bMat[1][actNode] = xiDx[actNode][1];
//    }
//  }

      
  //  Matrix<Double> auxbMat;
  //  auxbMat.Init();
  //  bMat.Transpose (auxbMat);
  //     std::cerr << "linElecInt::bMat transpose = \n" << auxbMat << std::endl;
    
}

void linElecInt::ApplyBMat( Vector<Double>& retVec,  
                            LocPointMapped& lp, BaseFE* ptFe, 
                            const Vector<Double>& solVec ) {

   // Get derivatives of local shape functions with respect to global
   // coords (format: nrNodes x spaceDim)
   Matrix<Double> xiDx;
   FeH1 *feH1 = (dynamic_cast<FeH1*>(ptFe));
   feH1->GetGlobDerivShFnc( xiDx, lp, lp.shapeMap->GetElem(), 1 );
   retVec = Transpose(xiDx) * solVec;
}
void linElecInt::ApplyBMat( Vector<Complex>& retVec,  
                            LocPointMapped& lp, BaseFE* ptFe, 
                            const Vector<Complex>& solVec ) {
  EXCEPTION("Not implemented");
}

  void linElecInt::calcBMat( Matrix<Double> &bMat, UInt ip,
			     Matrix<Double> &ptCoord ) {

    EXCEPTION("Implement me");
//
//    // Obtain info on number of elements' funtions
//    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
//    
//    // Set correct size of matrix B and initialise with zeros
//    bMat.Resize( dim_, numFncs );
//    bMat.Init();
//
//    // Get derivatives of local shape functions with respect to global
//    // coords (format: nrNodes x spaceDim)
//    Matrix<Double> xiDx;
//    ptelem->SetAnsatzFct( ansatzFct1_ );
//    if (isSetIntPoint_) {
//      ptelem->GetGlobDerivShFnc( xiDx, intPoint_, ptCoord, it1_.GetElem() );
//    } else {
//      ptelem->GetGlobDerivShFncAtIp( xiDx, ip, ptCoord, it1_.GetElem() );
//    }
//
//    if ( subTensorType_ == FULL ) {
//      // The matrix bMat can be seen as a 3 x numFncs block-vector.
//      // The k-th entry of this block vector corresponds to the matrix
//      // B of the BDB product evaluated at the k-th node of the finite
//      // element. 
//      for( UInt actNode = 0; actNode < numFncs; actNode++ ) {
//	bMat[0][actNode] = xiDx[actNode][0];
//	bMat[1][actNode] = xiDx[actNode][1];
//	bMat[2][actNode] = xiDx[actNode][2];
//      }
//    }
//    else  {
//      // The matrix bMat can be seen as a 1 x numFncs block-vector.
//      // The k-th entry of this block vector corresponds to the matrix
//      // B of the ADB product evaluated at the k-th node of the finite
//      // element. We assume that the first coordinate equals y and the
//      // second z.
//      for( UInt actNode = 0; actNode < numFncs; actNode++ ) {
//        bMat[0][actNode] = xiDx[actNode][0];
//        bMat[1][actNode] = xiDx[actNode][1];
//      }
//    }

    
// 	Matrix<Double> auxbMat;
// 	auxbMat.Init();
// 	bMat.Transpose (auxbMat);
//     std::cerr << "linElecInt::bMat transpose = \n" << auxbMat << std::endl;
  }


  // ============
  //   calcDMat
  // ============
  void linElecInt::calcDMat( Matrix<Double> &dMat, const Elem* elem) 
  {
    ptMaterial->GetTensor(dMat,ELEC_PERMITTIVITY,matDataType_,subTensorType_);
    dMat *= mParser_->Eval( mHandle_ );

//    Double density = elem != NULL ? GetErsatzMaterialFactor(elem) : 1.0;
//    if(density != 1.0) dMat *= density;  
//    LOG_DBG3(forms) << GetName() << "::calcDMat("
//                   << (elem != NULL ? Integer(elem->elemNum) : -1)
//                    << ") -> density=" << density;
  }


  void linElecInt::SetFactor( const std::string& factor ) {
    
    mParser_->SetExpr( mHandle_, factor );
  }
  

  // ================
  //   Constructors
  // ================


  linElecInt::linElecInt(BaseMaterial* matData, SubTensorType type,
                         bool coordUpdate ) 
    : BDBInt(matData, type, coordUpdate ) {


    name_ = "linElecInt";
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


}
