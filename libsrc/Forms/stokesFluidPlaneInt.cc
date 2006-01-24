#include <iostream>
#include <fstream>

#include "stokesFluidPlaneInt.hh"

namespace CoupledField
{

  StokesFluidPlaneInt::StokesFluidPlaneInt(Double aVal, Double bVal)
    : StokesFluidInt()
  {
    ENTER_FCN( "StokesFluidPlaneInt::StokesFluidPlaneInt" );
    density_ = aVal;
    dynamicViscosity_ = bVal;
  }


 
  StokesFluidPlaneInt::~StokesFluidPlaneInt()
  {
    ENTER_FCN( "StokesFluidPlaneInt::~StokesFluidPlaneInt" );
  }



  void StokesFluidPlaneInt::CalcElementMatrix(Matrix<Double> & ptCoord, 
                                           Matrix<Double> & elemMat)
  {
    ENTER_FCN( "StokesFluidPlaneInt::CalcElementMatrix" );
  
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const UInt nrNodes = ptelem->GetNumNodes();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;  
    UInt j, N;  // DOFs per Node

    // derivation of shape functions after global coordinates 
    Matrix<Double> partElemAMat, partElemATMat, locElemMat;
    Matrix<Double> xiDxDy;

    Vector<Double>  xi;
    Vector<Double> ShpFncAtIp;

//    Vector<Double> CoordAtIP;

    N = 4; // 4 DOFs per Node

    // set matrix to desired size and set all elements to zero
    elemMat.Resize(nrNodes*N); elemMat.Init();
    locElemMat.Resize(nrNodes*N); locElemMat.Init();


    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
        jacDet = 0;
        
        ptelem->GetShFncAtIp(xi, actIntPt);
        ptelem->GetGlobDerivShFncAtIp(xiDxDy, actIntPt, ptCoord, jacDet);

//        xiDx = xiDxDy[][0]
//        xiDy = xiDxDy[][1];

// __                                 __
// |                                    |
// |  xiDx   xiDy  0      0             |
// |  0      0     xiDx   dynVisc*xiDy  |
// |  0      0     xiDy  -dynVisc*xiDx  |
// |  xiDy  -xiDx  0      xi            |
// |__                                __|

        partElemAMat.Resize(N,nrNodes*N); partElemAMat.Init();

        for (j=0; j<nrNodes; j++)
          {
            partElemAMat[0][j]           =  xiDxDy[j][0];
            partElemAMat[0][j+nrNodes]   =  xiDxDy[j][1];

            partElemAMat[1][j+2*nrNodes] =  xiDxDy[j][0];
            partElemAMat[1][j+3*nrNodes] =  dynamicViscosity_ * xiDxDy[j][1];

            partElemAMat[2][j+2*nrNodes] =  xiDxDy[j][1];
            partElemAMat[2][j+3*nrNodes] = -dynamicViscosity_ * xiDxDy[j][0];

            partElemAMat[3][j]           =  xiDxDy[j][1];
            partElemAMat[3][j+nrNodes]   = -xiDxDy[j][0];
            partElemAMat[3][j+3*nrNodes] =  xi[j];


          }

        partElemAMat *= intWeights[actIntPt-1] * jacDet;
        partElemAMat.Transpose(partElemATMat);
        // assemble element matrix
        locElemMat += partElemATMat * partElemAMat;
      }
      
    ResortElementMatrix(elemMat, locElemMat, nrNodes, N);
//    std::cerr << "locElemMat:" << std::endl
//              << locElemMat << std::endl;
//    std::cerr << "elemMat:" << std::endl
//              << elemMat << std::endl;
  }


} // end namespace CoupledField
