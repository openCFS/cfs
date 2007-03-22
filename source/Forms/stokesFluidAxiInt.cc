// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "stokesFluidAxiInt.hh"

namespace CoupledField
{

  StokesFluidAxiInt::StokesFluidAxiInt(Double density,
                                       Double dynamicViscosity)
    : StokesFluidInt(density, dynamicViscosity)
  {
    ENTER_FCN( "StokesFluidAxiInt::StokesFluidAxiInt" );

    name_ = "StokesFluidAxiInt";
  }


 
  StokesFluidAxiInt::~StokesFluidAxiInt()
  {
    ENTER_FCN( "StokesFluidAxiInt::~StokesFluidAxiInt" );
  }



  void StokesFluidAxiInt::CalcElementMatrix( Matrix<Double>& elemMat,
                          EntityIterator& ent1, 
                          EntityIterator& ent2 )
  {
    ENTER_FCN( "StokesFluidAxiInt::CalcElementMatrix" );
  
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );
    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;  
    UInt j, N;  // DOFs per Node

    // derivation of shape functions after global coordinates 
    Matrix<Double> partElemAMat, partElemATMat, locElemMat;
    Matrix<Double> xiDxDy;

    Vector<Double>  xi, xiDx, xiDy;
    Vector<Double> ShpFncAtIp;

    Vector<Double> CoordAtIP;

    N = 4; // 4 DOFs per Node

    // set matrix to desired size and set all elements to zero
    elemMat.Resize(numFncs*N); 
    elemMat.Init();
    locElemMat.Resize(numFncs*N); 
    locElemMat.Init();


    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
        jacDet = 0;
        
        ptelem->GetShFncAtIp(xi, actIntPt, ent1.GetElem() );
        ptelem->GetGlobDerivShFncAtIp(xiDxDy, actIntPt, 
                                      ptCoord_, jacDet, ent1.GetElem() );

//        xiDx = xiDxDy.get_col(0);
//        xiDy = xiDxDy.get_col(1);

        partElemATMat.Resize(N,numFncs*N); partElemATMat.Init();
        for (j=0; j<numFncs; j++)
          {
            partElemATMat[0][j]           =xiDx[j];
            partElemATMat[0][j+numFncs]   =xiDy[j];

            partElemATMat[1][j+2*numFncs]=xiDx[j];
            partElemATMat[1][j+3*numFncs]=density_ * xiDy[j];

            partElemATMat[2][j+2*numFncs] =xiDy[j];
            partElemATMat[2][j+3*numFncs] =-density_ * xiDx[j];

            partElemATMat[3][j]          =xiDy[j];
            partElemATMat[3][j+numFncs]  =-xiDx[j];
            partElemATMat[3][j+3*numFncs]=xi[j];
          }

        partElemATMat *= intWeights[actIntPt-1] * jacDet;

        partElemATMat.Transpose(partElemAMat);

        // assemble element matrix
        locElemMat += partElemAMat * partElemATMat;
      }
    ResortElementMatrix(elemMat, locElemMat, numFncs, N);
  }


} // end namespace CoupledField
