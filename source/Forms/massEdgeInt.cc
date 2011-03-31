// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "massEdgeInt.hh"
#include "Elements/HCurlElemsHi.hh"
#include "Elements/fespace.hh"
namespace CoupledField
{

void MassEdgeInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                      EntityIterator& ent1, 
                                      EntityIterator& ent2 ) {
  

//  bool print = false;
//  if( ent1.GetElem()->elemNum == 10 ) print = true;
  // Extract physical element
  const Elem* ptElem = ent1.GetElem();
  // Obtain FE element from feSpace
  FeHCurlHi* ptFe = dynamic_cast<FeHCurlHi*>(ptFeSpace1_->GetFe( ent1 )); 
  
  // Special: Only use lower order functions
  ptFe->SetOnlyLowestOrder(true);
  
  
  UInt nrFncs = ptFe->BaseFE::GetNumFncs();

  // Get shape map from grid
  shared_ptr<ElemShapeMap> esm = domain->GetGrid()->GetElemShapeMap( ptElem );

  // Get integration points (shortcut: from basefe instead of 
  // IntegrationScheme class)
  StdVector<LocPoint> intPoints;
  StdVector<Double> weights;
  intScheme_->GetIntPoints( Elem::GetShapeType(ptElem->type), intPoints, weights );

  elemMat.Resize( nrFncs );
  elemMat.Init();

  // if scaling should be performed (regularization in static case),
  // we divide the jacobian determinant by h^2
  Double factor = conductivity_;
  Matrix<Double> partElemMat;
  if( scaleByEdgeSize_ ) {
    esm->GetMaxMinEdgeLength(maxEdgeLength_,minEdgeLength_);
    factor /= ( maxEdgeLength_ * maxEdgeLength_);
  }
//  if (print )
//    std::cerr << "factor of element #" << ent1.GetElem()->elemNum
//    << " is " << factor << std::endl;
  
  // Loop over all integration points
  LocPointMapped lp;
  Matrix<Double> shape;
  Double temp = 0.0;
  for( UInt i = 0; i < intPoints.GetSize(); i++  ) {
    // Calculate for each integration point the LocPointMapped
    lp.Set( intPoints[i], esm );
    
    ptFe->GetShFnc( shape, lp, lp.shapeMap->GetElem(), 0);
//    if (print )
//      std::cerr << "shape is \n" << shape << std::endl;
    
    partElemMat = Transpose(shape) * shape;
    temp = lp.jacDet * weights[i] * factor;
    elemMat += partElemMat * temp;
  }
  ptFe->SetOnlyLowestOrder(false);
}
//void MassEdgeInt::CalcElementMatrix( Matrix<Double>& elemMat,
//                                      EntityIterator& ent1, 
//                                      EntityIterator& ent2 ) {
//
//   // Extract pointer to reference element and get coordinates
//   ExtractElemInfo( ent1 );
//
//   const UInt nrIntPts = ptelem->GetNumIntPoints();
//   const UInt nrEdges  = ptelem->GetNumEdges();
//   const Vector<Double> & intWeights = ptelem->GetIntWeights();  
//   Double jacDet;
// 
//   // derivation of shape functions after global coordinates 
// 
//   Matrix<Double> shapeEdge;
//   Matrix<Double> shapeEdgeTransp;
//   Matrix<Double> partElemMat;
// 
//   // set matrix to desired size and set all elements to zero
//   elemMat.Resize(nrEdges);
//   elemMat.Init();
//   
//   // if scaling should be performed (regularization in static case),
//   // we divide the jacobian determinant by h^2
//   Double factor = conductivity_;
//   if( scaleByEdgeSize_ ) {
//     ptelem->GetMaxMinEdgeLength(ptCoord_,maxEdgeLength_,minEdgeLength_);
//     factor /= ( maxEdgeLength_ * maxEdgeLength_);
//   }
// 
//   for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
//     {
//       jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, 
//                                            ent1.GetElem() );
//     
//       ptelem->CalcEdgeShapeFncAtIp(shapeEdge, actIntPt, ptCoord_, 
//                                    ent1.GetElem());
//     
//       shapeEdge.Transpose(shapeEdgeTransp);
//       partElemMat = shapeEdge * shapeEdgeTransp;
//       partElemMat *= intWeights[actIntPt-1] * factor * jacDet; 
//       elemMat += partElemMat;
//     }
//
//    /*
//   (*debug) << "MassEdge Matrix:  "  << std::endl
//            << elemMat << std::endl << std::endl;
//    */
// }

 MassEdgeInt::MassEdgeInt( Double acond, 
                           bool scaleByEdgeSize,
                           bool coordUpdate )
   : BaseForm( NULL )
 {
   name_   = "MassEdgeInt";
   isaxi_  = false;
   conductivity_ = acond;
   coordUpdate_  = coordUpdate;
   baseType_     = MASS;
   scaleByEdgeSize_ = scaleByEdgeSize;
 }

 MassEdgeInt::~MassEdgeInt()
 {
 }
} // end namespace CoupledField
