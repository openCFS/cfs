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
  }


 
  StokesFluidAxiInt::~StokesFluidAxiInt()
  {
    ENTER_FCN( "StokesFluidAxiInt::~StokesFluidAxiInt" );
  }



  void StokesFluidAxiInt::CalcElementMatrix(Matrix<Double> & ptCoord, 
                                           Matrix<Double> & elemMat)
  {
    ENTER_FCN( "StokesFluidAxiInt::CalcElementMatrix" );
  
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const UInt nrNodes = ptelem->GetNumNodes();
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
    elemMat.Resize(nrNodes*N); elemMat.Init();
    locElemMat.Resize(nrNodes*N); locElemMat.Init();


    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
        jacDet = 0;
        
        ptelem->GetShFncAtIp(xi, actIntPt);
        ptelem->GetGlobDerivShFncAtIp(xiDxDy, actIntPt, ptCoord, jacDet);

//        xiDx = xiDxDy.get_col(0);
//        xiDy = xiDxDy.get_col(1);

        partElemATMat.Resize(N,nrNodes*N); partElemATMat.Init();
        for (j=0; j<nrNodes; j++)
          {
            partElemATMat[0][j]           =xiDx[j];
            partElemATMat[0][j+nrNodes]   =xiDy[j];

            partElemATMat[1][j+2*nrNodes]=xiDx[j];
            partElemATMat[1][j+3*nrNodes]=density_ * xiDy[j];

            partElemATMat[2][j+2*nrNodes] =xiDy[j];
            partElemATMat[2][j+3*nrNodes] =-density_ * xiDx[j];

            partElemATMat[3][j]          =xiDy[j];
            partElemATMat[3][j+nrNodes]  =-xiDx[j];
            partElemATMat[3][j+3*nrNodes]=xi[j];
          }

        partElemATMat *= intWeights[actIntPt-1] * jacDet;

        partElemATMat.Transpose(partElemAMat);

        // assemble element matrix
        locElemMat += partElemAMat * partElemATMat;
      }
    ResortElementMatrix(elemMat, locElemMat, nrNodes, N);
  }


} // end namespace CoupledField
