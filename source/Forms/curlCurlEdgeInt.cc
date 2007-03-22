// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "curlCurlEdgeInt.hh"

namespace CoupledField
{

  CurlCurlEdgeInt::CurlCurlEdgeInt(Double aVal)
    : BaseForm( NULL ), reluctivity_ (aVal)
  {
    ENTER_FCN( "CurlCurlEdgeInt::CurlCurlEdgeInt" );
    name_ = "CurlCurlEdgeInt";
  }


 
  CurlCurlEdgeInt::~CurlCurlEdgeInt()
  {
    ENTER_FCN( "CurlCurlEdgeInt::~CurlCurlEdgeInt" );
  }



  void CurlCurlEdgeInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                           EntityIterator& ent1, 
                                           EntityIterator& ent2 ) {
    ENTER_FCN( "CurlCurlEdgeInt::CalcElementMatrix" );

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );
  
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const UInt nrEdges = ptelem->GetNumEdges();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;  
  
  
    // derivation of shape functions after global coordinates 
    StdVector< Matrix<Double>* > xiDx;
    xiDx.Resize(nrEdges);
    for (UInt i=0; i<nrEdges; i++)
      xiDx[i] = new Matrix<Double>;
  
    Matrix<Double> curl;
    Matrix<Double> curlTransp;
    Matrix<Double> partElemMat;
  
  
  
  
    // set matrix to desired size and set all elements to zero
    elemMat.Resize(nrEdges,true); 
    elemMat.Init();
  
  
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
        // calc glob derivs of shape functions and jacobian determinante
        ptelem->GetEdgeGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_);
      
        CalcEdgeCurl(curl, xiDx);
      
        curl.Transpose(curlTransp);
      
        partElemMat = curlTransp * curl;
      
        jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, ent1.GetElem());
      
        partElemMat *= intWeights[actIntPt-1] * jacDet * reluctivity_;
      
        elemMat += partElemMat;
      }
  
  
#ifdef DEBUG 
    //      (*debug) << "CurlCurlEdgeInt: ElemMat " << std::endl
    //               << elemMat << std::endl
    //               << "\n reluctivity " << reluctivity_ << std::endl;
#endif

  }


  // calculates the curl, if the global derivates are already given in shapeDeriv
  void CurlCurlEdgeInt::CalcEdgeCurl(Matrix<Double>& curl, 
                                     const StdVector<Matrix<Double>*>& shapeDeriv)
  {
    UInt nrEdges = shapeDeriv.GetSize();
    UInt dim = shapeDeriv[0]->GetSizeRow();
    
    curl.Resize(dim, nrEdges);
    
    for (UInt actEdge=0; actEdge < nrEdges; actEdge++)
      for (UInt actDim=0; actDim < dim; actDim++)
        curl[actDim][actEdge] = 
          (*shapeDeriv[actEdge])[(actDim+2)%dim][(actDim+1)%dim] -
          (*shapeDeriv[actEdge])[(actDim+1)%dim][(actDim+2)%dim];
    
  }
  

} // end namespace CoupledField
