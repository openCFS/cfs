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
  // Extract physical element
  const Elem* ptElem = ent1.GetElem();
  // Obtain FE element from feSpace
  shared_ptr<IntScheme> intScheme;
  FeHCurlHi* ptFe = dynamic_cast<FeHCurlHi*>(ptFeSpace1_->GetFe( ent1, intScheme )); 
  
  // Special: Only use lower order functions
  ptFe->SetOnlyLowestOrder(true);
  UInt nrFncs = ptFe->BaseFE::GetNumFncs();

  // Get shape map from grid
  shared_ptr<ElemShapeMap> esm = domain->GetGrid()->GetElemShapeMap( ptElem );

  // Get integration points (shortcut: from basefe instead of 
  // IntegrationScheme class)
  StdVector<LocPoint> intPoints;
  StdVector<Double> weights;
  intScheme->GetIntPoints( Elem::GetShapeType(ptElem->type), intPoints, weights );

  elemMat.Resize( nrFncs );
  elemMat.Init();

  // if scaling should be performed (regularization in static case),
  // we divide the jacobian determinant by h^2
  Double factor = conductivity_;
  Matrix<Double> partElemMat;
  if( scaleByEdgeSize_ ) {
    esm->GetMaxMinEdgeLength(maxEdgeLength_,minEdgeLength_);
    factor /= ( minEdgeLength_ * minEdgeLength_);
  }

  // Loop over all integration points
  LocPointMapped lp;
  Matrix<Double> shape;
  for( UInt i = 0; i < intPoints.GetSize(); i++  ) {
    // Calculate for each integration point the LocPointMapped
    lp.Set( intPoints[i], esm );
    
    ptFe->GetShFnc( shape, lp, lp.shapeMap->GetElem(), 0);
    partElemMat =   Transpose(shape) * shape;
    partElemMat *= lp.jacDet * weights[i] * factor;
    elemMat += partElemMat;
  }
  ptFe->SetOnlyLowestOrder(false);
}

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
