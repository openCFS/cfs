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
    Double W1=1.0, W2=1.0, W3=1.0;
    UInt j, N;  // DOFs per Node

    // derivation of shape functions after global coordinates 
    Matrix<Double> partElemAMat, partElemATMat, locElemMat;

//     Matrix<Double> A, A_1, A_2, A_3;
//     Matrix<Double> B1, B1_1, B1_2, B1_3, B1_4;
//     Matrix<Double> B2, B2_1, B2_2, B2_3, B2_4;
//     Matrix<Double> B3, B3_1, B3_2, B3_3, B3_4;
//     Matrix<Double> C1, C1_1, C1_2;
//     Matrix<Double> C1T;
//     Matrix<Double> C2, C2_1, C2_2;
//     Matrix<Double> C2T;
//     Matrix<Double> C3, C3_1, C3_2;
//     Matrix<Double> C3T;
//     Matrix<Double> D1, mD1;
//     Matrix<Double> D2, mD2;
//     Matrix<Double> D3, mD3;
//     Matrix<Double> E1, mE1;
//     Matrix<Double> E2, mE2;
//     Matrix<Double> E3, mE3;
//     Matrix<Double> F1, F1_1, F1_2;
//     Matrix<Double> F2, F2_1, F2_2;
//     Matrix<Double> F3, F3_1, F3_2;
//     Matrix<Double> G1, G1_1, G1_2;
//     Matrix<Double> G2, G2_1, G2_2;
//     Matrix<Double> G3, G3_1, G3_2;
//     Matrix<Double> H1, H1_1, H1_2;
//     Matrix<Double> H2, H2_1, H2_2;
//     Matrix<Double> H3, H3_1, H3_2;
//     Matrix<Double> I1, I1_1, I1_2;
//     Matrix<Double> I2, I2_1, I2_2;
//     Matrix<Double> I3, I3_1, I3_2;

    Matrix<Double> xiDxDyDz;

    //Vector<Double>  xiDx, xiDy, xiDz;

    Vector<Double> xi;

    N = 7; // 7 DOFs per Node

    // set matrix to desired size and set all elements to zero
    elemMat.Resize(nrNodes*N); elemMat.Init(0.0);
    locElemMat.Resize(nrNodes*N); locElemMat.Init(0.0);

    //xiDx.Resize(nrNodes);
    //xiDy.Resize(nrNodes);
    //xiDz.Resize(nrNodes);

    //std::cerr << "nrIntPts=" << nrIntPts << std::endl;
    
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
        jacDet = 0;
        
        ptelem->GetShFncAtIp(xi, actIntPt);
        ptelem->GetGlobDerivShFncAtIp(xiDxDyDz, actIntPt, ptCoord, jacDet);

/////////////////////////////////////////////////////////////////////
//      Submatrix for 7 unknowns
/////////////////////////////////////////////////////////////////////
//        for (UInt i=0; i< nrNodes; i++) 
//          {
//            xiDx[i] = xiDxDyDz[i][0];
//            xiDy[i] = xiDxDyDz[i][1];
//            xiDz[i] = xiDxDyDz[i][2];
//          }
//
//        A_1.DyadicMult(xiDz,xiDz);
//        A_2.DyadicMult(xiDy,xiDy);
//        A_3.DyadicMult(xiDx,xiDx);
//        A = A_1 + A_2 + A_3;
//        A *= intWeights[actIntPt-1] * jacDet;
//
//        B1_1=(dynamicViscosity_*dynamicViscosity_) * A_1;
//        B1_2=(dynamicViscosity_*dynamicViscosity_) * A_2;
//        B1_3=A_3;
//        B1_4.DyadicMult(xi,xi);
//        B1 = B1_1 + B1_2 + B1_3 + B1_4;
//        B1 *= intWeights[actIntPt-1] * jacDet;
//
//        B2_1=B1_1;
//        B2_2= A_2;
//        B2_3=(dynamicViscosity_*dynamicViscosity_) * A_3;
//        B2_4=B1_4;
//        B2 = B2_1 + B2_2 + B2_3 + B2_4;
//        B2 *= intWeights[actIntPt-1] * jacDet;
//
//        B3_1=A_1;
//        B3_2=B1_2;
//        B3_3=B2_3;
//        B3_4=B1_4;
//        B3 = B3_1 + B3_2 + B3_3 + B3_4;
//        B3 *= intWeights[actIntPt-1] * jacDet;
//
//        C1_1.DyadicMult(xiDx,xiDy);
//        C1_2.DyadicMult(xiDy,xiDx);
//        C1=C1_2-C1_1;
//        C1 *= intWeights[actIntPt-1] * jacDet;
//
//        C1T=C1_1-C1_2;
//        C1T *= intWeights[actIntPt-1] * jacDet;
//
//        C2_1.DyadicMult(xiDx,xiDz);
//        C2_2.DyadicMult(xiDz,xiDx);
//        C2=C2_2-C2_1;
//        C2 *= intWeights[actIntPt-1] * jacDet;
//
//        C2T=C2_1-C2_2;
//        C2T *= intWeights[actIntPt-1] * jacDet;
//
//        C3_1.DyadicMult(xiDy,xiDz);
//        C3_2.DyadicMult(xiDz,xiDy);
//        C3=C3_2-C3_1;
//        C3 *= intWeights[actIntPt-1] * jacDet;
//
//        C3T=C3_1-C3_2;
//        C3T *= intWeights[actIntPt-1] * jacDet;
//
//        D1.DyadicMult(xi,xiDz);
//        D1 *= intWeights[actIntPt-1] * jacDet;
//        mD1 = (-1.0)*D1;
//
//        D2.DyadicMult(xi,xiDy);
//        D2 *= intWeights[actIntPt-1] * jacDet;
//        mD2 = (-1.0)*D2;
//
//        D3.DyadicMult(xi,xiDx);
//        D3 *= intWeights[actIntPt-1] * jacDet;
//        mD3 = (-1.0)*D3;
//
//        E1.DyadicMult(xiDz,xi);
//        E1 *= intWeights[actIntPt-1] * jacDet;
//        mE1 = (-1.0)*E1;
//
//        E2.DyadicMult(xiDy,xi);
//        E2 *= intWeights[actIntPt-1] * jacDet;
//        mE2 = (-1.0)*E2;
//
//        E3.DyadicMult(xiDx,xi);
//        E3 *= intWeights[actIntPt-1] * jacDet;
//        mE3 = (-1.0)*E3;
//
//        F1_1.DyadicMult(xiDx,xiDy);
//        F1_2.DyadicMult(xiDy,xiDx);
//        F1=((-1.0)*dynamicViscosity_*dynamicViscosity_* F1_1) + F1_2;
//        F1 *= intWeights[actIntPt-1] * jacDet;
//
//        F2_1.DyadicMult(xiDx,xiDz);
//        F2_2.DyadicMult(xiDz,xiDx);
//        F2=((-1.0)*dynamicViscosity_*dynamicViscosity_* F2_1) + F2_2;
//        F2 *= intWeights[actIntPt-1] * jacDet;
//
//        F3_1.DyadicMult(xiDy,xiDz);
//        F3_2.DyadicMult(xiDz,xiDy);
//        F3=((-1.0)*dynamicViscosity_*dynamicViscosity_* F3_1) + F3_2;
//        F3 *= intWeights[actIntPt-1] * jacDet;
//
//        G1_1.DyadicMult(xiDy,xiDx);
//        G1_2.DyadicMult(xiDx,xiDy);
//        G1=((-1.0)*dynamicViscosity_*dynamicViscosity_* G1_1) + G1_2;
//        G1 *= intWeights[actIntPt-1] * jacDet;
//
//        G2_1.DyadicMult(xiDz,xiDx);
//        G2_2.DyadicMult(xiDx,xiDz);
//        G2=((-1.0)*dynamicViscosity_*dynamicViscosity_* G2_1) + G2_2;
//        G2 *= intWeights[actIntPt-1] * jacDet;
//
//        G3_1.DyadicMult(xiDz,xiDy);
//        G3_2.DyadicMult(xiDy,xiDz);
//        G3=((-1.0)*dynamicViscosity_*dynamicViscosity_* G3_1) + G3_2;
//        G3 *= intWeights[actIntPt-1] * jacDet;
//
//        H1_1.DyadicMult(xiDy,xiDz);
//        H1_2.DyadicMult(xiDz,xiDy);
//        H1=(dynamicViscosity_ * H1_1) + ((-1.0)*dynamicViscosity_ * H1_2);
//        H1 *= intWeights[actIntPt-1] * jacDet;
//
//        H2_1.DyadicMult(xiDx,xiDz);
//        H2_2.DyadicMult(xiDz,xiDx);
//        H2=((-1.0)*dynamicViscosity_ * H2_1) + (dynamicViscosity_ * H2_2);
//        H2 *= intWeights[actIntPt-1] * jacDet;
//
//        H3_1.DyadicMult(xiDx,xiDy);
//        H3_2.DyadicMult(xiDy,xiDx);
//        H3=(dynamicViscosity_ * H3_1) + ((-1.0)*dynamicViscosity_ * H3_2);
//        H3 *= intWeights[actIntPt-1] * jacDet;
//
//        I1_1.DyadicMult(xiDz,xiDy);
//        I1_2.DyadicMult(xiDy,xiDz);
//        I1=(dynamicViscosity_ * I1_1) + ((-1.0)*dynamicViscosity_ * I1_2);
//        I1 *= intWeights[actIntPt-1] * jacDet;
//
//        I2_1.DyadicMult(xiDz,xiDx);
//        I2_2.DyadicMult(xiDx,xiDz);
//        I2=((-1.0)*dynamicViscosity_ * I2_1) + (dynamicViscosity_ * I2_2);
//        I2 *= intWeights[actIntPt-1] * jacDet;
//
//        I3_1.DyadicMult(xiDy,xiDx);
//        I3_2.DyadicMult(xiDx,xiDy);
//        I3=(dynamicViscosity_ * I3_1) + ((-1.0)*dynamicViscosity_ * I3_2);
//        I3 *= intWeights[actIntPt-1] * jacDet;
//
//
//        //Diagonal of Element-Matrix
//	locElemMat.AddSubMatrix(A       ,  0,          0);
//        locElemMat.AddSubMatrix(A       ,  nrNodes,    nrNodes);
//        locElemMat.AddSubMatrix(A       ,  2*nrNodes,  2*nrNodes);
//        locElemMat.AddSubMatrix(B1      ,  3*nrNodes,  3*nrNodes);
//        locElemMat.AddSubMatrix(B2      ,  4*nrNodes,  4*nrNodes);
//        locElemMat.AddSubMatrix(B3      ,  5*nrNodes,  5*nrNodes);
//        locElemMat.AddSubMatrix(A       ,  6*nrNodes,  6*nrNodes);
//
//        locElemMat.AddSubMatrix(C1      ,  nrNodes,    0);
//        locElemMat.AddSubMatrix(C1T     ,  0,          nrNodes);
//
//        locElemMat.AddSubMatrix(C2      ,  2*nrNodes,  0);
//        locElemMat.AddSubMatrix(C2T     ,  0,          2*nrNodes);
//
//        locElemMat.AddSubMatrix(C3      ,  2*nrNodes,  nrNodes);
//        locElemMat.AddSubMatrix(C3T     ,  nrNodes,    2*nrNodes);
//
//        locElemMat.AddSubMatrix(D1      ,  3*nrNodes,  nrNodes);
//        locElemMat.AddSubMatrix(mD1     ,  4*nrNodes,  0);
//
//        locElemMat.AddSubMatrix(D2      ,  5*nrNodes,  0);
//        locElemMat.AddSubMatrix(mD2     ,  3*nrNodes,  2*nrNodes);
//
//        locElemMat.AddSubMatrix(D3      ,  4*nrNodes,  2*nrNodes);
//        locElemMat.AddSubMatrix(mD3     ,  5*nrNodes,  nrNodes);
//
//        locElemMat.AddSubMatrix(E1      ,  nrNodes,    3*nrNodes);
//        locElemMat.AddSubMatrix(mE1     ,  0,          4*nrNodes);
//
//        locElemMat.AddSubMatrix(E2      ,  0,          5*nrNodes);
//        locElemMat.AddSubMatrix(mE2     ,  2*nrNodes,  3*nrNodes);
//
//        locElemMat.AddSubMatrix(E3      ,  2*nrNodes,  4*nrNodes);
//        locElemMat.AddSubMatrix(mE3     ,  nrNodes,    5*nrNodes);
//
//        locElemMat.AddSubMatrix(F1      ,  4*nrNodes,  3*nrNodes);
//        locElemMat.AddSubMatrix(F2      ,  5*nrNodes,  3*nrNodes);
//        locElemMat.AddSubMatrix(F3      ,  5*nrNodes,  4*nrNodes);
//
//        locElemMat.AddSubMatrix(G1      ,  3*nrNodes,  4*nrNodes);
//        locElemMat.AddSubMatrix(G2      ,  3*nrNodes,  5*nrNodes);
//        locElemMat.AddSubMatrix(G3      ,  4*nrNodes,  5*nrNodes);
//
//        locElemMat.AddSubMatrix(H1      ,  6*nrNodes,  3*nrNodes);
//        locElemMat.AddSubMatrix(H2      ,  6*nrNodes,  4*nrNodes);
//        locElemMat.AddSubMatrix(H3      ,  6*nrNodes,  5*nrNodes);
//
//        locElemMat.AddSubMatrix(I1      ,  3*nrNodes,  6*nrNodes);
//        locElemMat.AddSubMatrix(I2      ,  4*nrNodes,  6*nrNodes);
//        locElemMat.AddSubMatrix(I3      ,  5*nrNodes,  6*nrNodes);
//
/////////////////////////////////////////////////////////////////////
//      A-Matrix for 8 unknowns
/////////////////////////////////////////////////////////////////////
//  __                                                                    __
//  |  0      1     2      3     4         5        6        7             |   
// 0|  0      0     0      0     0        -mu*xiDz  mu*xiDy  xiDx          |   u
// 1|  0      0     0      0     mu*xiDz   0       -mu*xiDx  xiDy          |   v
// 2|  0      0     0      0    -mu*xiDy   mu*xiDx  0        xiDz          |   w
// 3|  0      0     0      0     xiDx     xiDy      xiDz     0             |   phi
// 4|  0     -xiDz  xiDy   xiDx -xi       0         0        0             |   omx
// 5|  xiDz   0    -xiDx   xiDy  0       -xi        0        0             |   omy
// 6| -xiDy   xiDx  0      xiDz  0        0        -xi       0             |   omy
// 7|  xiDx   xiDy  xiDz   0     0        0         0        0             |   p
//  |_ 0      1     2      3     4         5        6        7           __|

/////////////////////////////////////////////////////////////////////
//      A-Matrix for 7 unknowns
/////////////////////////////////////////////////////////////////////
//  __                                                              __
//  |  0         1        2         3         4         5        6             |
// 0|  0         0        0         0        -mu*xiDz   mu*xiDy  xiDx          |   u
// 1|  0         0        0         mu*xiDz   0        -mu*xiDx  xiDy          |   v
// 2|  0         0        0        -mu*xiDy   mu*xiDx   0        xiDz          |   w
// 3|  0         0        0         W1*xiDx   W1*xiDy   W1*xiDz  0             |   omx
// 4|  0        -W3*xiDz  W3*xiDy  -W3*       0         0        0             |   omy
// 5|  W3*xiDz   0       -W3*xiDx   0        -W3*xi     0        0             |   omz
// 6| -W3*xiDy   W3*xiDx  0         0         0        -W3*xi    0             |   p
// 7|  W2*xiDx   W2*xiDy  W2*xiDz   0         0         0        0             |
//  |_ 0         1        2         3         4         5        6           __|

/////////////////////////////////////////////////////////////////////
//      A-Matrix for 8 unknowns
/////////////////////////////////////////////////////////////////////
//        partElemAMat.Resize(N,nrNodes*N); partElemAMat.Init();
/////////////////////////////////////////////////////////////////////
//      A-Matrix for 7 unknowns
/////////////////////////////////////////////////////////////////////
        partElemAMat.Resize(N+1,nrNodes*N); partElemAMat.Init(0.0);
//
        for (j=0; j<nrNodes; j++)
          {
///////////////////////////////////////////////////////////////////////
////      A-Matrix for 8 unknowns
///////////////////////////////////////////////////////////////////////
//            partElemAMat[0][j+5*nrNodes]   = -dynamicViscosity_*xiDxDyDz[j][2];
//            partElemAMat[0][j+6*nrNodes]   =  dynamicViscosity_*xiDxDyDz[j][1];
//            partElemAMat[0][j+7*nrNodes]   =  xiDxDyDz[j][0];
//            
//            partElemAMat[1][j+4*nrNodes]   =  dynamicViscosity_*xiDxDyDz[j][2];
//            partElemAMat[1][j+6*nrNodes]   = -dynamicViscosity_*xiDxDyDz[j][0];
//            partElemAMat[1][j+7*nrNodes]   =  xiDxDyDz[j][1];
//
//            partElemAMat[2][j+4*nrNodes]   = -dynamicViscosity_*xiDxDyDz[j][1];
//            partElemAMat[2][j+5*nrNodes]   =  dynamicViscosity_*xiDxDyDz[j][0];
//            partElemAMat[2][j+7*nrNodes]   =  xiDxDyDz[j][2];
//
//            partElemAMat[3][j+4*nrNodes]   =  xiDxDyDz[j][0];
//            partElemAMat[3][j+5*nrNodes]   =  xiDxDyDz[j][1];
//            partElemAMat[3][j+6*nrNodes]   =  xiDxDyDz[j][2];
//
//            partElemAMat[4][j+1*nrNodes]   = -xiDxDyDz[j][2];
//            partElemAMat[4][j+2*nrNodes]   =  xiDxDyDz[j][1];
//            partElemAMat[4][j+3*nrNodes]   =  xiDxDyDz[j][0];
//            partElemAMat[4][j+4*nrNodes]   = -xi[j];
//
//            partElemAMat[5][j]             =  xiDxDyDz[j][2];
//            partElemAMat[5][j+2*nrNodes]   = -xiDxDyDz[j][0];
//            partElemAMat[5][j+3*nrNodes]   =  xiDxDyDz[j][1];
//            partElemAMat[5][j+5*nrNodes]   = -xi[j];
//
//            partElemAMat[6][j]             = -xiDxDyDz[j][1];
//            partElemAMat[6][j+1*nrNodes]   =  xiDxDyDz[j][0];
//            partElemAMat[6][j+3*nrNodes]   =  xiDxDyDz[j][2];
//            partElemAMat[6][j+6*nrNodes]   = -xi[j];
//
//            partElemAMat[7][j]             =  xiDxDyDz[j][0];
//            partElemAMat[7][j+1*nrNodes]   =  xiDxDyDz[j][1];
//            partElemAMat[7][j+2*nrNodes]   =  xiDxDyDz[j][2];
//
//
/////////////////////////////////////////////////////////////////////
//      A-Matrix for 7 unknowns
/////////////////////////////////////////////////////////////////////
             partElemAMat[0][j+4*nrNodes]   = -dynamicViscosity_*xiDxDyDz[j][2];
             partElemAMat[0][j+5*nrNodes]   =  dynamicViscosity_*xiDxDyDz[j][1];
             partElemAMat[0][j+6*nrNodes]   =  xiDxDyDz[j][0];

             partElemAMat[1][j+3*nrNodes]   =  dynamicViscosity_*xiDxDyDz[j][2];
             partElemAMat[1][j+5*nrNodes]   = -dynamicViscosity_*xiDxDyDz[j][0];
             partElemAMat[1][j+6*nrNodes]   =  xiDxDyDz[j][1];

             partElemAMat[2][j+3*nrNodes]   = -dynamicViscosity_*xiDxDyDz[j][1];
             partElemAMat[2][j+4*nrNodes]   =  dynamicViscosity_*xiDxDyDz[j][0];
             partElemAMat[2][j+6*nrNodes]   =  xiDxDyDz[j][2];

             partElemAMat[3][j+3*nrNodes]   =  W1*xiDxDyDz[j][0];
             partElemAMat[3][j+4*nrNodes]   =  W1*xiDxDyDz[j][1];
             partElemAMat[3][j+5*nrNodes]   =  W1*xiDxDyDz[j][2];

             partElemAMat[4][j+1*nrNodes]   = -W3*xiDxDyDz[j][2];
             partElemAMat[4][j+2*nrNodes]   =  W3*xiDxDyDz[j][1];
             partElemAMat[4][j+3*nrNodes]   = -W3*xi[j];

             partElemAMat[5][j]             =  W3*xiDxDyDz[j][2];
             partElemAMat[5][j+2*nrNodes]   = -W3*xiDxDyDz[j][0];
             partElemAMat[5][j+4*nrNodes]   = -W3*xi[j];

             partElemAMat[6][j]             = -W3*xiDxDyDz[j][1];
             partElemAMat[6][j+1*nrNodes]   =  W3*xiDxDyDz[j][0];
             partElemAMat[6][j+5*nrNodes]   = -W3*xi[j];

             partElemAMat[7][j]             =  W2*xiDxDyDz[j][0];
             partElemAMat[7][j+1*nrNodes]   =  W2*xiDxDyDz[j][1];
             partElemAMat[7][j+2*nrNodes]   =  W2*xiDxDyDz[j][2];
          }

        partElemAMat.Transpose(partElemATMat);

//        // assemble element matrix
        locElemMat += (partElemATMat * partElemAMat) * intWeights[actIntPt-1] * jacDet;
      }
    ResortElementMatrix(elemMat, locElemMat, nrNodes, N);
  }

} // end namespace CoupledField
