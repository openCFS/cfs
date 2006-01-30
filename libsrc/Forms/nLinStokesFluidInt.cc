#include <iostream>
#include <fstream>

#include "nLinStokesFluidInt.hh"

namespace CoupledField
{



  // =============================================================================
  // base class for nonlinear stokes fluid
  // =============================================================================

  // base class for nonlinear calculation of elasticity
  nLinStokesFluidInt::nLinStokesFluidInt(Double density, Double dynamicViscosity) 
    : StokesFluidInt(density, dynamicViscosity)
  {
    ENTER_FCN( "nLinStokesFluidInt::nLinStokesFluidInt" );
  }


  nLinStokesFluidInt::~nLinStokesFluidInt()
  {
    ENTER_FCN( "nLinStokesFluidInt::~nLinStokesFluidInt" );
  }


  // returns nonlinear B - matrix (first part) for BDB
//  void nLinStokesFluidInt::calcBMat(Matrix<Double> & bMat, UInt ip, Matrix<Double> & ptCoord)
//  {
//    ENTER_FCN( "nLinStokesFluidInt::calcBMat" );
//
//    
//    if (!elemVelocity_.GetSizeRow() || !elemVelocity_.GetSizeCol()) 
//      Error("Undefined velocities! ",__FILE__,__LINE__);
//
//    const UInt nrNodes = ptelem->GetNumNodes();
//    const UInt nrDofs  = 8;    
//    const UInt dimD    = 3;
////    const UInt nrDofs  = getNrDofs();    
////    const UInt dimD    = getDimD();
//
//    Matrix<Double> linBMat;    
//    //StokesFluidInt::calcBMat( linBMat, ip, ptCoord);
//
//
//    bMat.Resize(dimD, nrNodes * nrDofs);
//
//    
//    // local shape functions derived after global coords 
//    // (format: nrNodes x nrDofs)
//    Matrix<Double> xiDx;
//    ptelem->GetGlobDerivShFncAtIp(xiDx, ip, ptCoord);
//
//
//    // calculate derivates of global displacements at element nodes
//    // elemDisp_ holds displacement of the actual element (dimension: nrDofs x nrNodes)
//    // xiDx holds derivatives of shape functions after global coords in one IP (dimension: nrNodes x nrDofs)
//    // displDeriv dimension: nrDofs x nrDofs
//    Matrix<Double>  velocityDeriv;  
//    velocityDeriv = elemVelocity_ * xiDx;
//
//    // we need transposed of displacement derivatives
//    Matrix<Double> velocityDerivTransp;
//    velocityDeriv.Transpose(velocityDerivTransp);  
//
//
//    Matrix<Double> bMatOneNode;
//    // size_row() = nr of rows!!
//    bMatOneNode.Resize(linBMat.GetSizeRow(), nrDofs);
//
//    
//    Matrix<Double> helpMat;
//    for( UInt actNode=0; actNode < nrNodes; actNode++)
//      {
//        linBMat.GetSubMatrix(bMatOneNode, 0, actNode*nrDofs);
//
//        helpMat = bMatOneNode * velocityDerivTransp;
//
//        bMat.SetSubMatrix(helpMat, 0, actNode*nrDofs);
//      }
//
//
//    if (isaxi_)
//      {
//        const UInt spaceDim  = ptelem->GetDim();
//
//        Vector<Double> shpFncAtIp;
//        ptelem->GetShFncAtIp(shpFncAtIp, ip);
//
//        Vector<Double> coordAtIP;
//        coordAtIP = ptCoord * shpFncAtIp;
//
//        Double  l33=0;  // for l33, see Bathe, page 552
//        for (UInt actPos=0; actPos < shpFncAtIp.GetSize(); actPos++)
//          l33 += elemVelocity_[0][actPos] * shpFncAtIp[actPos] / coordAtIP[0];
//        
//        
//        for (UInt actNode = 0; actNode < nrNodes; actNode++)      
//          {         
//            //  (N_a/r) * (u_x/r)
//            bMat[dimD - 1][actNode * spaceDim]     = l33 * shpFncAtIp[actNode] / coordAtIP[0];
//            bMat[dimD - 1][actNode * spaceDim + 1] = 0;
//          }
//      }
//
//  }
  
  // =============================================================================
  // 3D nonlinear stokes fluid (part with convective term)
  // =============================================================================


  nLinStokesFluid3DInt_Convective::nLinStokesFluid3DInt_Convective(Double density,
                                                                   Double dynamicViscosity) 
    : nLinStokesFluidInt(density, dynamicViscosity)
  {
    ENTER_FCN( "nLinStokesFluid3DInt_Convective::nLinStokesFluid3DInt_Convective" );

    className = "nLinStokesFluid3DInt_Convective";
  }

  nLinStokesFluid3DInt_Convective::~nLinStokesFluid3DInt_Convective()
  {
    ENTER_FCN( "nLinStokesFluid3DInt_Convective::~nLinStokesFluid3DInt_Convective" );
  }

  void nLinStokesFluid3DInt_Convective::CalcElementMatrix(Matrix<Double> & ptCoord, 
                                                          Matrix<Double> & elemMat)
  {
    ENTER_FCN( "nLinStokesFluid3DInt_Convective::CalcElementMatrix" );

    if (!elemVelocity_.GetSizeRow() || !elemVelocity_.GetSizeCol()) 
      Error("Undefined velocities! ",__FILE__,__LINE__);
  
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

    N = 8; // 8 DOFs per Node

    // set matrix to desired size and set all elements to zero
    elemMat.Resize(nrNodes*N); elemMat.Init();
    locElemMat.Resize(nrNodes*N); locElemMat.Init();

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
        jacDet = 0;
        
        ptelem->GetShFncAtIp(xi, actIntPt);
        ptelem->GetGlobDerivShFncAtIp(xiDxDyDz, actIntPt, ptCoord, jacDet);

//  __                                                                                       __
//  |  0                       1                        2                      3  4  5  6  7  |

// 0|  u*xiDx+v*xiDy+w*xiDz    0                        0                      0  0  0  0  0  |
// 1|  0                       u*xiDx+v*xiDy+w*xiDz     0                      0  0  0  0  0  |
// 2|  0                       0                        u*xiDx+v*xiDy+w*xiDz   0  0  0  0  0  |
// 3|  0                       0                        0                      0  0  0  0  0  |
// 4|  0                       0                        0                      0  0  0  0  0  |
// 5|  0                       0                        0                      0  0  0  0  0  |
// 6|  0                       0                        0                      0  0  0  0  0  |
// 7|  0                       0                        0                      0  0  0  0  0 _|

        Matrix<Double>  velocityDeriv;  
        velocityDeriv = elemVelocity_ * xiDxDyDz;

        // we need transposed of velocity 
        Matrix<Double> velocityTransp;
        elemVelocity_.Transpose(velocityTransp);  

        // we need transposed of velocity derivatives
        Matrix<Double> velocityDerivTransp;
        velocityDeriv.Transpose(velocityDerivTransp);  

        partElemAMat.Resize(N,nrNodes*N); partElemAMat.Init();
        for (j=0; j<nrNodes; j++)
          {
            partElemAMat[0][j]   =           velocityTransp[j][0]*xiDxDyDz[j][0]
                                           + velocityTransp[j][1]*xiDxDyDz[j][1]
                                           + velocityTransp[j][2]*xiDxDyDz[j][2]
                                           + velocityDerivTransp[j][0]; //v_xDx;
            partElemAMat[0][j+1*nrNodes] =   velocityDerivTransp[j][1]; //v_xDy;
            partElemAMat[0][j+2*nrNodes] =   velocityDerivTransp[j][2]; //v_xDz;
            
            partElemAMat[1][j]          =    velocityDerivTransp[j+1*nrNodes][0]; //v_yDx;
            partElemAMat[1][j+1*nrNodes] =   velocityTransp[j][0]*xiDxDyDz[j][0]
                                           + velocityTransp[j][1]*xiDxDyDz[j][1]
                                           + velocityTransp[j][2]*xiDxDyDz[j][2]
                                           + velocityDerivTransp[j+1*nrNodes][1];
            partElemAMat[1][j+2*nrNodes] =   velocityDerivTransp[j+1*nrNodes][2]; //v_yDz;
            

            partElemAMat[2][j]          =    velocityDerivTransp[j+2*nrNodes][0]; //v_zDx;
            partElemAMat[2][j+1*nrNodes] =   velocityDerivTransp[j+2*nrNodes][1]; //v_zDy;
            partElemAMat[2][j+2*nrNodes] =   velocityTransp[j][0]*xiDxDyDz[j][0]
                                           + velocityTransp[j][1]*xiDxDyDz[j][1]
                                           + velocityTransp[j][2]*xiDxDyDz[j][2]
                                           + velocityDerivTransp[j+2*nrNodes][2]; //v_zDz;
          }

        partElemAMat *= intWeights[actIntPt-1] * jacDet;

        partElemAMat.Transpose(partElemATMat);

        // assemble element matrix
        locElemMat += partElemATMat * partElemAMat;
      }
    ResortElementMatrix(elemMat, locElemMat, nrNodes, N);
  }


  // =============================================================================
  // 2D nonlinear plane stokes fluid
  // =============================================================================

  

  void nLinStokesFluidPlaneInt_Convective::CalcElementMatrix(Matrix<Double> & ptCoord, 
                                                             Matrix<Double> & elemMat)
  {
    ENTER_FCN( "nLinStokesFluidPlaneInt_Convective::CalcElementMatrix" );
  
    if (!elemVelocity_.GetSizeRow() || !elemVelocity_.GetSizeCol()) 
      Error("Undefined velocities! ",__FILE__,__LINE__);
  
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

    N = 4; // 4 DOFs per Node

    // set matrix to desired size and set all elements to zero
    elemMat.Resize(nrNodes*N); elemMat.Init();
    locElemMat.Resize(nrNodes*N); locElemMat.Init();

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
        jacDet = 0;
        
        ptelem->GetShFncAtIp(xi, actIntPt);
        ptelem->GetGlobDerivShFncAtIp(xiDxDyDz, actIntPt, ptCoord, jacDet);

//  __                                                            __
//  |  0                          1                          2   3 |
// 0   0                          0                          0   0 |
// 1|  u0*xiDx + v0*xiDy + u0Dx   u0Dy    0                  0   0 |
// 2|  v0Dx                       u0*xiDx + v0*xiDy + v0Dy   0   0 |
// 3|  0                          0                          0   0 |

        Matrix<Double>  velocityDeriv;  
        velocityDeriv = elemVelocity_ * xiDxDyDz;

        // we need transposed of velocity 
        Matrix<Double> velocityTransp;
        elemVelocity_.Transpose(velocityTransp);  

        // we need transposed of velocity derivatives
        Matrix<Double> velocityDerivTransp;
        velocityDeriv.Transpose(velocityDerivTransp);  

        partElemAMat.Resize(N,nrNodes*N); partElemAMat.Init();
        for (j=0; j<nrNodes; j++)
          {
            partElemAMat[1][j]   =           velocityTransp[j][0]*xiDxDyDz[j][0]
                                           + velocityTransp[j][1]*xiDxDyDz[j][1]
                                           + velocityDerivTransp[j][0];
            partElemAMat[1][j+1*nrNodes] =   velocityDerivTransp[j][1];
            
            partElemAMat[2][j]          =    velocityDerivTransp[j+1*nrNodes][0];
            partElemAMat[2][j+1*nrNodes] =   velocityTransp[j][0]*xiDxDyDz[j][0]
                                           + velocityTransp[j][1]*xiDxDyDz[j][1]
                                           + velocityDerivTransp[j+1*nrNodes][1];
          }

        partElemAMat *= intWeights[actIntPt-1] * jacDet;

        partElemAMat.Transpose(partElemATMat);

        // assemble element matrix
        locElemMat += partElemATMat * partElemAMat;
      }
    ResortElementMatrix(elemMat, locElemMat, nrNodes, N);

  }
  
  // =============================================================================
  // 2D nonlinear axisymmetric stokes fluid
  // =============================================================================
 
  // calculated the D-matrix for the axi state
  void  nLinStokesFluidAxiInt_Convective::CalcElementMatrix(Matrix<Double> & ptCoord, 
                                                             Matrix<Double> & elemMat)
  {
    ENTER_FCN( "nLinStokesFluidAxiInt_Convective::CalcElementMatrix" );
  }


  // ===================================================================================
  // nonlinear calculation of elasticity in plane strain state 
  // ===================================================================================

  nLinStokesFluidPlaneInt_Convective::nLinStokesFluidPlaneInt_Convective(Double density, Double dynamicViscosity) 
    : nLinStokesFluidInt(density, dynamicViscosity)
  {
    ENTER_FCN( "nLinStokesFluidPlaneInt_Convective::nLinStokesFluidPlaneInt_Convective" );

    className = "nLinStokesFluidPlaneInt_Convective";
  }

  nLinStokesFluidPlaneInt_Convective::~nLinStokesFluidPlaneInt_Convective()
  {
    ENTER_FCN( "nLinStokesFluidPlaneInt_Convective::~nLinStokesFluidPlaneInt_Convective" );
  }


  // ===================================================================================
  // axisymmetric nonlinear calculation of elasticity
  // ===================================================================================

  nLinStokesFluidAxiInt_Convective::nLinStokesFluidAxiInt_Convective(Double density, Double dynamicViscosity) 
    : nLinStokesFluidInt(density, dynamicViscosity)
  {
    ENTER_FCN( "nLinStokesFluidAxiInt_Convective::nLinStokesFluidAxiInt_Convective" );

    isaxi_ = TRUE;
    className = "nLinStokesFluidAxiInt_Convective";    
  }

  nLinStokesFluidAxiInt_Convective::~nLinStokesFluidAxiInt_Convective()
  {
    ENTER_FCN( "nLinStokesFluidAxiInt_Convective::~nLinStokesFluidAxiInt_Convective" );

  }


  // ==================================================================
  // nLinStokesFluid
  // ==================================================================

  nLinStokesFluid_linFormInt::nLinStokesFluid_linFormInt(Double density,
                                                         Double dynamicViscosity,
                                                         Boolean isaxi) 
    : nLinStokesFluidInt(density, dynamicViscosity)
  {
    ENTER_FCN( "nLinStokesFluid_linFormInt::nLinStokesFluid_linFormInt" );
    isaxi_ = isaxi;
  }



  nLinStokesFluid_linFormInt::~nLinStokesFluid_linFormInt()
  {
    ENTER_FCN( "nLinStokesFluid_linFormInt::~nLinStokesFluid_linFormInt" );
  }



  void nLinStokesFluid_linFormInt::CalcElemVector(Matrix<Double>& ptCoord, 
                                           Vector<Double> & elemVec)
  {
    ENTER_FCN( "nLinStokesFluid_linFormInt::CalcElemVector" );

    Matrix<Double> xiDxDyDz;

    Double jacDet;  

    const UInt nrIntPts = ptelem->GetNumIntPoints();
    const UInt nrNodes  = ptelem->GetNumNodes();
    UInt N; 
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  

    Vector<Double> partElemVec;
  
    if (ptelem->GetDim() == 2 && !isaxi_)
      {
        N=4;
      }
    else if (ptelem->GetDim() == 2 && isaxi_)
      {
        N=4;
      }
    else if (ptelem->GetDim() == 3)
      {
        N=8;
      }
    else {
      Error("Wrong space dimension of elements! ",__FILE__,__LINE__);
    }  
  
    if (!elemVelocity_.GetSizeRow() || !elemVelocity_.GetSizeCol()) 
       Error("Undefined velocities! ",__FILE__,__LINE__);
  
    partElemVec.Resize(nrNodes * N);

    elemVec.Resize(nrNodes*N);
    elemVec *= 0;    // set elems to 0
  
//    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
//      {    
//        partElemVec *= jacDet * intWeights[actIntPt-1];
//        elemVec +=  partElemVec;
//      }
      
      
      
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
        jacDet = 0;
        
        //ptelem->GetShFncAtIp(xi, actIntPt);
        ptelem->GetGlobDerivShFncAtIp(xiDxDyDz, actIntPt, ptCoord, jacDet);

//   __                                            __                                                                   __
// 0|  v_x*(v_x*xiDx)+v_y*(v_x*xiDy)+v_z*(v_x*xiDz)  |
// 1|  v_x*(v_y*xiDx)+v_y*(v_y*xiDy)+v_z*(v_y*xiDz)  |
// 2|  v_x*(v_z*xiDx)+v_y*(v_z*xiDy)+v_z*(v_z*xiDz)  |
// 3|                      0                         |
// 4|                      0                         |
// 5|                      0                         |
// 6|                      0                         |
// 7|__                    0                       __|

//   __                              __                                                                   __
// 0|                 0                |
// 1|  v_x*(v_x*xiDx) + v_y*(v_x*xiDy) |
// 2|  v_x*(v_y*xiDx) + v_y*(v_y*xiDy) |
// 3|__               0              __|

        Matrix<Double> xiDxDyTransp;
        xiDxDyDz.Transpose(xiDxDyTransp); 
        
        Matrix<Double>  velocityDeriv;  
        velocityDeriv = elemVelocity_ * xiDxDyDz;

        // we need transposed of velocity 
        Matrix<Double> velocityTransp;
        elemVelocity_.Transpose(velocityTransp);  

        // we need transposed of velocity derivatives
        Matrix<Double> velocityDerivTransp;
        velocityDeriv.Transpose(velocityDerivTransp);  

        for (UInt j=0; j<nrNodes; j++)
          {
            partElemVec[j*N+1]   =  velocityTransp[j][0]*velocityDerivTransp[j][0]
                              + velocityTransp[j][1]*velocityDerivTransp[j][1];
//                              + velocityTransp[j][2]*velocityDerivTransp[j][2];
  
            partElemVec[j*N+2]   =  velocityTransp[j][0]*velocityDerivTransp[j+nrNodes][0]
                                + velocityTransp[j][1]*velocityDerivTransp[j+nrNodes][1];
//                              + velocityTransp[j][2]*velocityDerivTransp[j+1*nrNodes][2];

//            partElemVec[j+2*nrNodes]   =  velocityTransp[j][1]*velocityTransp[j][0]*xiDxDyDz[j][0]
//                              + velocityTransp[j][1]*velocityTransp[j][1]*xiDxDyDz[j][1]
//                              + velocityTransp[j][2]*velocityTransp[j][1]*xiDxDyDz[j][2];

            
          }

        partElemVec *= intWeights[actIntPt-1] * jacDet;
        elemVec +=  partElemVec;
      }
      
      
  }
} // end namespace CoupledField
