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

    name_ = "nLinStokesFluidInt";
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

    name_ = "nLinStokesFluid3DInt_Convective";
  }

  nLinStokesFluid3DInt_Convective::~nLinStokesFluid3DInt_Convective()
  {
    ENTER_FCN( "nLinStokesFluid3DInt_Convective::~nLinStokesFluid3DInt_Convective" );
  }

  void nLinStokesFluid3DInt_Convective::CalcElementMatrix( Matrix<Double>& elemMat,
                                                           EntityIterator& ent1, 
                                                           EntityIterator& ent2 ) 
    
  {
    ENTER_FCN( "nLinStokesFluid3DInt_Convective::CalcElementMatrix" );

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );

    if (!elemVelocity_.GetSizeRow() || !elemVelocity_.GetSizeCol()) 
      Error("Undefined velocities! ",__FILE__,__LINE__);

    ptelem->SetAnsatzFct( ansatzFct1_ );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();

    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;  
    UInt j, N;  // DOFs per Node

    // derivation of shape functions after global coordinates 
    Matrix<Double> partElemAMat, partElemATMat, locElemMat;
    Matrix<Double> xiDxDyDz;

    Vector<Double>  xi;
    Vector<Double> ShpFncAtIp;

    N = 7; // 7 DOFs per Node

    // set matrix to desired size and set all elements to zero
    elemMat.Resize(numFncs*N); 
    elemMat.Init();
    locElemMat.Resize(numFncs*N); 
    locElemMat.Init();

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
        jacDet = 0;
        
        ptelem->GetShFncAtIp(xi, actIntPt, ent1.GetElem() );
        ptelem->GetGlobDerivShFncAtIp(xiDxDyDz, actIntPt, ptCoord_, 
                                      jacDet, ent1.GetElem() );

//  __ Sucessive-Substitution                                                          __
//  __                                                                                       __
//  |  0                       1                        2                      3  4  5  6  7  |

// 0|  u*xiDx+v*xiDy+w*xiDz    0                        0                      0  0  0  0  0  |
// 1|  0                       u*xiDx+v*xiDy+w*xiDz     0                      0  0  0  0  0  |
// 2|  0                       0                        u*xiDx+v*xiDy+w*xiDz   0  0  0  0  0  |
// 3|  0                       0                        0                      0  0  0  0  0  |
// 4|  0                       0                        0                      0  0  0  0  0  |
// 5|  0                       0                        0                      0  0  0  0  0  |
// 6|  0                       0                        0                      0  0  0  0  0 _|

        Matrix<Double>  vxM, vyM, vzM, vxDxDyDzAtIP, vyDxDyDzAtIP, vzDxDyDzAtIP;  
        Vector<Double>  vxAtIP, vyAtIP, vzAtIP;

        vxM.Resize(numFncs,1);
        vyM.Resize(numFncs,1);
        vzM.Resize(numFncs,1);

        for (UInt i=0; i<numFncs; i++) {
            vxM[i][0]=elemVelocity_[0][i];
            vyM[i][0]=elemVelocity_[1][i];
            vzM[i][0]=elemVelocity_[2][i];
        }

        vxAtIP     = vxM * xi;
        vyAtIP     = vyM * xi;
        vzAtIP     = vzM * xi;

        //vxDxDyDzAtIP = vxM * xiDxDyDz;
        //vyDxDyDzAtIP = vyM * xiDxDyDz;
        //vzDxDyDzAtIP = vzM * xiDxDyDz;

        partElemAMat.Resize(N,numFncs*N); 
        partElemAMat.Init();
        for (j=0; j<numFncs; j++)
          {
//          Sucessive-Substitution                                                 
            partElemAMat[0][j]   =           vxAtIP[j]*xiDxDyDz[j][0]
                                           + vyAtIP[j]*xiDxDyDz[j][1]
                                           + vzAtIP[j]*xiDxDyDz[j][2];

            partElemAMat[1][j+1*numFncs] =   vxAtIP[j]*xiDxDyDz[j][0]
                                           + vyAtIP[j]*xiDxDyDz[j][1]
                                           + vzAtIP[j]*xiDxDyDz[j][2];
            
            partElemAMat[2][j+2*numFncs] =   vxAtIP[j]*xiDxDyDz[j][0]
                                           + vyAtIP[j]*xiDxDyDz[j][1]
                                           + vzAtIP[j]*xiDxDyDz[j][2];

//             Newton-Method
//             partElemAMat[0][j]   =           vxAtIP[j]*xiDxDyDz[j][0]
//                                            + vyAtIP[j]*xiDxDyDz[j][1]
//                                            + vzAtIP[j]*xiDxDyDz[j][2]

//             partElemAMat[0][j+1*numFncs] =   velocityDerivTransp[j][1]; //v_xDy;
//             partElemAMat[0][j+2*numFncs] =   velocityDerivTransp[j][2]; //v_xDz;
            
//             partElemAMat[1][j]          =    velocityDerivTransp[j+1*numFncs][0]; //v_yDx;
//             partElemAMat[1][j+1*numFncs] =   velocityTransp[j][0]*xiDxDyDz[j][0]
//                                            + velocityTransp[j][1]*xiDxDyDz[j][1]
//                                            + velocityTransp[j][2]*xiDxDyDz[j][2]
//                                            + velocityDerivTransp[j+1*numFncs][1];
//             partElemAMat[1][j+2*numFncs] =   velocityDerivTransp[j+1*numFncs][2]; //v_yDz;
            

//             partElemAMat[2][j]          =    velocityDerivTransp[j+2*numFncs][0]; //v_zDx;
//             partElemAMat[2][j+1*numFncs] =   velocityDerivTransp[j+2*numFncs][1]; //v_zDy;
//             partElemAMat[2][j+2*numFncs] =   velocityTransp[j][0]*xiDxDyDz[j][0]
//                                            + velocityTransp[j][1]*xiDxDyDz[j][1]
//                                            + velocityTransp[j][2]*xiDxDyDz[j][2]
//                                            + velocityDerivTransp[j+2*numFncs][2]; //v_zDz;
          }

        partElemAMat *= intWeights[actIntPt-1] * jacDet;

        partElemAMat.Transpose(partElemATMat);

        // assemble element matrix
        locElemMat += partElemATMat * partElemAMat;
      }
    ResortElementMatrix(elemMat, locElemMat, numFncs, N);
  }


  // =============================================================================
  // 2D nonlinear plane stokes fluid
  // =============================================================================

  

  void nLinStokesFluidPlaneInt_Convective::CalcElementMatrix( Matrix<Double>& elemMat,
                                                              EntityIterator& ent1, 
                                                              EntityIterator& ent2 )
  {
    ENTER_FCN( "nLinStokesFluidPlaneInt_Convective::CalcElementMatrix" );
  
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );

    if (!elemVelocity_.GetSizeRow() || !elemVelocity_.GetSizeCol()) 
      Error("Undefined velocities! ",__FILE__,__LINE__);

    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
   
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;  
    UInt j, N;  // DOFs per Node

    // derivation of shape functions after global coordinates 
    Matrix<Double> partElemAMat, partElemATMat, locElemMat;
    Matrix<Double> xiDxDy;

    Vector<Double>  xi;
    Vector<Double> ShpFncAtIp;

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
        ptelem->GetGlobDerivShFncAtIp(xiDxDy, actIntPt, ptCoord_, 
                                      jacDet, ent1.GetElem() );

//  __ Sucessive-Substitution                                                          __
//  |  0                   1                   2   3 |
// 0   0                   0                   0   0 |
// 1|  u0*xiDx + v0*xiDy   0                   0   0 |
// 2|  0                   u0*xiDx + v0*xiDy   0   0 |
// 3|  0                   0                   0   0 |

//  __ Newton-Verfahren                                                           __
//  |  0                          1                          2   3 |
// 0   0                          0                          0   0 |
// 1|  u0*xiDx + v0*xiDy + u0Dx   u0Dy    0                  0   0 |
// 2|  v0Dx                       u0*xiDx + v0*xiDy + v0Dy   0   0 |
// 3|  0                          0                          0   0 |

        Matrix<Double>  vxM, vyM, vxDxDyAtIP, vyDxDyAtIP;  
        Vector<Double>  vxAtIP, vyAtIP;

        vxM.Resize(numFncs,1);
        vyM.Resize(numFncs,1);

        for (UInt i=0; i<numFncs; i++) {
            vxM[i][0]=elemVelocity_[0][i];
            vyM[i][0]=elemVelocity_[1][i];
        }

        vxAtIP     = vxM * xi;
        vyAtIP     = vyM * xi;
        //vxDxDyAtIP = vxM * xiDxDy;
        //vyDxDyAtIP = vyM * xiDxDy;

        partElemAMat.Resize(N,numFncs*N); 
        partElemAMat.Init();
        for (j=0; j<numFncs; j++)
          {
//          Sucessive Substitution
           partElemAMat[1][j]   =           vxAtIP[j]*xiDxDy[j][0]
                                          + vyAtIP[j]*xiDxDy[j][1];
           
           partElemAMat[2][j+1*numFncs] =   vxAtIP[j]*xiDxDy[j][0]
                                          + vyAtIP[j]*xiDxDy[j][1];

//            Newton-Verfahren
//             partElemAMat[1][j]   =           vxAtIP[j]*xiDxDy[j][0]
//                                            + vyAtIP[j]*xiDxDy[j][1]
//                                            + vxDxDyAtIP[j][0];

//             partElemAMat[1][j+1*numFncs] =   vxDxDyAtIP[j][1];

//             partElemAMat[2][j]          =    vyDxDyAtIP[j][0];

//             partElemAMat[2][j+1*numFncs] =   vxAtIP[j]*xiDxDy[j][0]
//                                            + vyAtIP[j]*xiDxDy[j][1]
//                                            + vyDxDyAtIP[j][1];
          }

        partElemAMat *= intWeights[actIntPt-1] * jacDet;

        partElemAMat.Transpose(partElemATMat);

        // assemble element matrix
        locElemMat += partElemATMat * partElemAMat;
      }
    ResortElementMatrix(elemMat, locElemMat, numFncs, N);

  }
  
  // =============================================================================
  // 2D nonlinear axisymmetric stokes fluid
  // =============================================================================
 
  // calculated the D-matrix for the axi state
  void  nLinStokesFluidAxiInt_Convective::CalcElementMatrix( Matrix<Double>& stiffMat,
                                                             EntityIterator& ent1, 
                                                             EntityIterator& ent2 )
  {
    ENTER_FCN( "nLinStokesFluidAxiInt_Convective::CalcElementMatrix" );

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );
  }


  // ===================================================================================
  // nonlinear calculation of elasticity in plane strain state 
  // ===================================================================================

  nLinStokesFluidPlaneInt_Convective::
  nLinStokesFluidPlaneInt_Convective(Double density, 
                                     Double dynamicViscosity) 
    : nLinStokesFluidInt(density, dynamicViscosity)
  {
    ENTER_FCN( "nLinStokesFluidPlaneInt_Convective::nLinStokesFluidPlaneInt_Convective" );

    name_ = "nLinStokesFluidPlaneInt_Convective";
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

    name_ = "nLinStokesFluid3DInt_Convective";
    isaxi_ = true;
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
                                                         bool isaxi) 
    : nLinStokesFluidInt(density, dynamicViscosity)
  {
    ENTER_FCN( "nLinStokesFluid_linFormInt::nLinStokesFluid_linFormInt" );

    name_ = "nLinStokesFluid_linFormInt";
    isaxi_ = isaxi;
  }



  nLinStokesFluid_linFormInt::~nLinStokesFluid_linFormInt()
  {
    ENTER_FCN( "nLinStokesFluid_linFormInt::~nLinStokesFluid_linFormInt" );
  }



  void nLinStokesFluid_linFormInt::CalcElemVector(Vector<Double> & elemVec,
                                                  EntityIterator& ent )
  {
    ENTER_FCN( "nLinStokesFluid_linFormInt::CalcElemVector" );

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    Matrix<Double> xiDxDyDz;

    Double jacDet;  

    ptelem->SetAnsatzFct( ansatzFct1_ );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    

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
  
    partElemVec.Resize(numFncs * N);
    partElemVec.Init();

    elemVec.Resize(numFncs*N);
    elemVec.Init();
  
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
        jacDet = 0;
        
        //ptelem->GetShFncAtIp(xi, actIntPt);
        ptelem->GetGlobDerivShFncAtIp(xiDxDyDz, actIntPt, ptCoord_, 
                                      jacDet, ent.GetElem() );

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
        velocityDeriv = elemVelocity_ * xiDxDyDz; // hier gibts ein Problem:
        // die Velocity ist an jedem Knoten gegeben und aBleitungungen gibts nur
        // zu jedem Integrationspunkt.
        // Wenn Also Die Anzahl an Knoten und Integrationspunkte nicht
        //übereinstimmt gibts ein Problem!

        // we need transposed of velocity 
        Matrix<Double> velocityTransp;
        elemVelocity_.Transpose(velocityTransp);  

        // we need transposed of velocity derivatives
        Matrix<Double> velocityDerivTransp;
        velocityDeriv.Transpose(velocityDerivTransp);  

        for (UInt j=0; j<numFncs; j++)
          {
//            partElemVec[j*N+1]   =  velocityTransp[j][0]*velocityDerivTransp[j][0]
//                              + velocityTransp[j][1]*velocityDerivTransp[j][1];
////                              + velocityTransp[j][2]*velocityDerivTransp[j][2];
//  
//            partElemVec[j*N+2]   =  velocityTransp[j][0]*velocityDerivTransp[j+nrNodes][0]
//                                + velocityTransp[j][1]*velocityDerivTransp[j+nrNodes][1];
////                              + velocityTransp[j][2]*velocityDerivTransp[j+1*nrNodes][2];
//
////            partElemVec[j+2*nrNodes]   =  velocityTransp[j][1]*velocityTransp[j][0]*xiDxDyDz[j][0]
////                              + velocityTransp[j][1]*velocityTransp[j][1]*xiDxDyDz[j][1]
////                              + velocityTransp[j][2]*velocityTransp[j][1]*xiDxDyDz[j][2];

            
          }

        partElemVec *= intWeights[actIntPt-1] * jacDet;
        // Funzt noch net!!!
        //elemVec +=  partElemVec; 
      }
      
      
  }
} // end namespace CoupledField
