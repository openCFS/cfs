#include <iostream>
#include <fstream>

#include "stokesFluid3DInt.hh"

namespace CoupledField
{

  StokesFluid3DInt::StokesFluid3DInt(Double density,
                                     Double dynamicViscosity)
    : StokesFluidInt(density, dynamicViscosity)
  {
    ENTER_FCN( "StokesFluid3DInt::StokesFluid3DInt" );
  }


 
  StokesFluid3DInt::~StokesFluid3DInt()
  {
    ENTER_FCN( "StokesFluid3DInt::~StokesFluid3DInt" );
  }



  void StokesFluid3DInt::CalcElementMatrix(Matrix<Double> & ptCoord, 
                                           Matrix<Double> & elemMat)
  {
    ENTER_FCN( "StokesFluid3DInt::CalcElementMatrix" );
  
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const UInt nrNodes = ptelem->GetNumNodes();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;  
    UInt j, N;  // DOFs per Node

    // derivation of shape functions after global coordinates 
    Matrix<Double> partElemAMat, partElemATMat, locElemMat;
    Matrix<Double> xiDxDyDz;

    Vector<Double>  xi;
    Vector<Double> ShpFncAtIp;

//    Vector<Double> CoordAtIP;

    N = 8; // 8 DOFs per Node

    // set matrix to desired size and set all elements to zero
    elemMat.Resize(nrNodes*N); elemMat.Init();
    locElemMat.Resize(nrNodes*N); locElemMat.Init();


    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
        jacDet = 0;
        
        ptelem->GetShFncAtIp(xi, actIntPt);
        ptelem->GetGlobDerivShFncAtIp(xiDxDyDz, actIntPt, ptCoord, jacDet);

//  __                                                                    __
//  |  0      1     2      3     4         5        6        7                     |
// 0|  0      0     0      0     0        -mu*xiDz  mu*xiDy  xiDx          |
// 1|  0      0     0      0     mu*xiDz   0       -mu*xiDx  xiDy          |
// 2|  0      0     0      0    -mu*xiDy   mu*xiDx  0        xiDz          |
// 3|  0      0     0      0     xiDx     xiDy      xiDz     0             |
// 4|  0     -xiDz  xiDy   xiDx -xi       0         0        0             |
// 5|  xiDz   0    -xiDx   xiDy  0       -xi        0        0             |
// 6| -xiDy   xiDx  0      xiDz  0        0        -xi       0             |
// 7|  xiDx   xiDy  xiDz   0     0        0         0        0             |
//  |_ 0      1     2      3     4         5        6        7           __|


        partElemAMat.Resize(N,nrNodes*N); partElemAMat.Init();
        for (j=0; j<nrNodes; j++)
          {
            partElemAMat[0][j+5*nrNodes]   = -dynamicViscosity_*xiDxDyDz[j][2];
            partElemAMat[0][j+6*nrNodes]   =  dynamicViscosity_*xiDxDyDz[j][1];
            partElemAMat[0][j+7*nrNodes]   =  xiDxDyDz[j][0];

            partElemAMat[1][j+4*nrNodes]   =  dynamicViscosity_*xiDxDyDz[j][2];
            partElemAMat[1][j+6*nrNodes]   = -dynamicViscosity_*xiDxDyDz[j][0];
            partElemAMat[1][j+7*nrNodes]   =  xiDxDyDz[j][1];

            partElemAMat[2][j+4*nrNodes]   = -dynamicViscosity_*xiDxDyDz[j][1];
            partElemAMat[2][j+5*nrNodes]   =  dynamicViscosity_*xiDxDyDz[j][0];
            partElemAMat[2][j+7*nrNodes]   =  xiDxDyDz[j][2];

            partElemAMat[3][j+4*nrNodes]   =  xiDxDyDz[j][0];
            partElemAMat[3][j+5*nrNodes]   =  xiDxDyDz[j][1];
            partElemAMat[3][j+6*nrNodes]   =  xiDxDyDz[j][2];

            partElemAMat[4][j+1*nrNodes]   = -xiDxDyDz[j][2];
            partElemAMat[4][j+2*nrNodes]   =  xiDxDyDz[j][1];
            partElemAMat[4][j+3*nrNodes]   =  xiDxDyDz[j][0];
            partElemAMat[4][j+4*nrNodes]   = -xi[j];

            partElemAMat[5][j]             =  xiDxDyDz[j][2];
            partElemAMat[5][j+2*nrNodes]   = -xiDxDyDz[j][0];
            partElemAMat[5][j+3*nrNodes]   =  xiDxDyDz[j][1];
            partElemAMat[5][j+5*nrNodes]   = -xi[j];

            partElemAMat[6][j]             = -xiDxDyDz[j][1];
            partElemAMat[6][j+1*nrNodes]   =  xiDxDyDz[j][0];
            partElemAMat[6][j+3*nrNodes]   =  xiDxDyDz[j][2];
            partElemAMat[6][j+6*nrNodes]   = -xi[j];

            partElemAMat[7][j]             =  xiDxDyDz[j][0];
            partElemAMat[7][j+1*nrNodes]   =  xiDxDyDz[j][1];
            partElemAMat[7][j+2*nrNodes]   =  xiDxDyDz[j][2];

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
