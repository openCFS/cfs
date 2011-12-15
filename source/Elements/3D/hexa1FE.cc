// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <stdlib.h>
#include <algorithm>
#include <bitset>
#include <iostream>

#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "Domain/elem.hh"
#include "Elements/3D/hexaFE.hh"
#include "Elements/basefe.hh"
#include "General/exception.hh"
#include "MatVec/exprt/xpr1.hh"
#include "MatVec/exprt/xpr2.hh"
#include "MatVec/matrix.hh"
#include "Utils/StdVector.hh"
#include "Utils/tools.hh"
//#include <General/environment.hh>
#include "hexa1FE.hh"

namespace CoupledField
{

  // declare class specific logging stream
  DECLARE_LOG(hexa1fe)
  DEFINE_LOG(hexa1fe,"hexa1fe")

  Hexa1FE::Hexa1FE():HexaFE()
  {

    // Evaluate legendre polynomials of order
    // 1-10 in the unit interval

 //    Double delta = 0.01;
//     UInt nums = (int) 2.0/delta;
//     Double x = -1.0;
//     Double val, deriv;
//     for( UInt i = 0; i <= nums; i++ ) {
//       std::cerr << x << "\t";
//       for( UInt o = 1; o <= 8; o++ ) {
//         EvalPolynom( val, deriv, o, lCoeff_[o], x  );
//         std::cerr << val << "\t";
//       }
//       x+=delta;
//       std::cerr << "\n";
//     }

    UseICModes();
    Init();
  }

  Hexa1FE::~Hexa1FE()
  {

  }

  void Hexa1FE::Init()
  {

    NumNodes_ = 8;

    CommonInit();
    SetEdgeIndices();
    SetFaceIndices();
  }

  void Hexa1FE::SetCornerCoords()
  {

    LCornerCoords_.Resize(Dim_,NumNodes_);

    LCornerCoords_[0][0] =  -1;
    LCornerCoords_[1][0] =  -1;
    LCornerCoords_[2][0] =  -1;

    LCornerCoords_[0][1] =  1;
    LCornerCoords_[1][1] =  -1;
    LCornerCoords_[2][1] =  -1;

    LCornerCoords_[0][2] =  1;
    LCornerCoords_[1][2] =  1;
    LCornerCoords_[2][2] =  -1;

    LCornerCoords_[0][3] =  -1;
    LCornerCoords_[1][3] =  1;
    LCornerCoords_[2][3] =  -1;

    LCornerCoords_[0][4] =  -1;
    LCornerCoords_[1][4] =  -1;
    LCornerCoords_[2][4] =  1;

    LCornerCoords_[0][5] =  1;
    LCornerCoords_[1][5] =  -1;
    LCornerCoords_[2][5] =  1;

    LCornerCoords_[0][6] =  1;
    LCornerCoords_[1][6] =  1;
    LCornerCoords_[2][6] =  1;

    LCornerCoords_[0][7] =  -1;
    LCornerCoords_[1][7] =  1;
    LCornerCoords_[2][7] =  1;

  }


  void Hexa1FE::SetEdgeIndices() {

    edgeIndices_ = new StdVector<UInt>[NumEdges_];
    for (UInt i=0; i<NumEdges_; i++) {
      edgeIndices_[i].Resize(2);
    }

    // Note: The orientation is taken from
    // A. Duester: High Order FEM, Lecture Nodes,
    // p. 31

    // edge 1
    edgeIndices_[0][0] = 1;
    edgeIndices_[0][1] = 2;

    // edge 2
    edgeIndices_[1][0] = 2;
    edgeIndices_[1][1] = 3;

    // edge 3
    edgeIndices_[2][0] = 4;
    edgeIndices_[2][1] = 3;

    // edge 4
    edgeIndices_[3][0] = 1;
    edgeIndices_[3][1] = 4;

    // edge 5
    edgeIndices_[4][0] = 1;
    edgeIndices_[4][1] = 5;

    // edge 6
    edgeIndices_[5][0] = 2;
    edgeIndices_[5][1] = 6;

    // edge 7
    edgeIndices_[6][0] = 3;
    edgeIndices_[6][1] = 7;

    // edge 8
    edgeIndices_[7][0] = 4;
    edgeIndices_[7][1] = 8;

    // edge 9
    edgeIndices_[8][0] = 5;
    edgeIndices_[8][1] = 6;

    // edge 10
    edgeIndices_[9][0] = 6;
    edgeIndices_[9][1] = 7;

    // edge 11
    edgeIndices_[10][0] = 8;
    edgeIndices_[10][1] = 7;

    // edge 12
    edgeIndices_[11][0] = 5;
    edgeIndices_[11][1] = 8;


  }

  void Hexa1FE::SetFaceIndices() {

    faceIndices_ = new StdVector<UInt>[NumFaces_];

    for (UInt i = 0; i < NumFaces_; i++) {
      faceIndices_[i].Resize(4);
    }


    // face 1
    faceIndices_[0][0] = 1;
    faceIndices_[0][1] = 2;
    faceIndices_[0][2] = 3;
    faceIndices_[0][3] = 4;

    // face 2
    faceIndices_[1][0] = 1;
    faceIndices_[1][1] = 2;
    faceIndices_[1][2] = 6;
    faceIndices_[1][3] = 5;

    // face 3
    faceIndices_[2][0] = 2;
    faceIndices_[2][1] = 3;
    faceIndices_[2][2] = 7;
    faceIndices_[2][3] = 6;

    // face 4
    faceIndices_[3][0] = 4;
    faceIndices_[3][1] = 3;
    faceIndices_[3][2] = 7;
    faceIndices_[3][3] = 8;

    // face 5
    faceIndices_[4][0] = 1;
    faceIndices_[4][1] = 4;
    faceIndices_[4][2] = 8;
    faceIndices_[4][3] = 5;

    // face 6
    faceIndices_[5][0] = 5;
    faceIndices_[5][1] = 6;
    faceIndices_[5][2] = 7;
    faceIndices_[5][3] = 8;

  }


  void Hexa1FE :: CalcShapeFnc(Vector<Double> & Shape,
                               const Vector<Double> & actCoord,
                               const Elem* elem, UInt dof,
                               AnsatzFct::FctEntityType fcnType )
  {

    // Check ansatzFctType
    if(  actFct_->GetType() == AnsatzFct::LAGRANGE ||
         fcnType == AnsatzFct::NODE ) {

      // ===============
      //  LAGRANGE PART
      // ===============

      Shape.Resize(NumNodes_);

      for( UInt i=0; i<NumNodes_; i++)
        Shape[i] = 0.125
          * (1 + LCornerCoords_[0][i] * actCoord[0])
          * (1 + LCornerCoords_[1][i] * actCoord[1])
          * (1 + LCornerCoords_[2][i] * actCoord[2]);
    }else if(  actFct_->GetType() == AnsatzFct::SPECTRAL) {
      // ===============
      //  Spectral PART
      // ===============
      CalcSpectralShFct(Shape,actCoord,elem,dof,fcnType);

    } else if ( actFct_->GetType() == AnsatzFct::LEGENDRE ) {

      // ===============
      //  LEGENDRE PART
      // ===============

      // Get number of ansatz functions
      UInt totalFcns = GetNumFncs( actFct_ );
      Shape.Resize( totalFcns );

      // Get number of ansatz functions
      shared_ptr<LegendreFct> legFct =
        dynamic_pointer_cast<LegendreFct, AnsatzFct>(actFct_);

      // --------------------
      //  a) nodal functions
      // --------------------
      // Offset for different functions
      UInt offset = 0;

      // Calculate nodal shape functions
      for( UInt iNode = 0; iNode < NumCorners_; iNode++, offset++ )
        Shape[offset] = 0.125
          * (1 + LCornerCoords_[0][iNode] * actCoord[0])
          * (1 + LCornerCoords_[1][iNode] * actCoord[1])
          * (1 + LCornerCoords_[2][iNode] * actCoord[2]);


      // --------------------
      //  b) edge functions
      // --------------------
      // Obtain order of element
      Integer order_psi;
      Double val, factor, deriv;
#define HEX_EDGE_FCN(edgeNum,  sign1, dir1, sign2, dir2, psi_dir)       \
      order_psi = legFct->GetMaxOrderLocDir(psi_dir);                   \
      factor = elem->edges[edgeNum-1] < 0 ? -1.0 : 1.0;                 \
      for( Integer iDof = 2; iDof <= order_psi; iDof++, offset++ ) {    \
        EvalPolynom( val, deriv, iDof, lCoeff_[iDof],                   \
                     factor*actCoord[psi_dir] );                        \
        Shape[offset] = 0.25 * ( 1.0 sign1 actCoord[dir1] )             \
                             * ( 1.0 sign2 actCoord[dir2] ) * val;      \
      }


      // Edge #1
      HEX_EDGE_FCN(1, -, 1, -, 2, 0);

      // Edge #2
      HEX_EDGE_FCN(2, +, 0, -, 2, 1);

      // Edge #3
      HEX_EDGE_FCN(3, +, 1, -, 2, 0);

      // Edge #4
      HEX_EDGE_FCN(4, -, 0, -, 2, 1);

      // Edge #5
      HEX_EDGE_FCN(5, -, 0, -, 1, 2);

      // Edge #6
      HEX_EDGE_FCN(6, +, 0, -, 1, 2);

      // Edge #7
      HEX_EDGE_FCN(7, +, 0, +, 1, 2);

      // Edge #8
      HEX_EDGE_FCN(8, -, 0, +, 1, 2);

      // Edge #9
      HEX_EDGE_FCN(9, -, 1, +, 2, 0);

      // Edge #10
      HEX_EDGE_FCN(10, +, 0, +, 2, 1);

      // Edge #11
      HEX_EDGE_FCN(11, +, 1, +, 2, 0);

      // Edge #12
      HEX_EDGE_FCN(12, -, 0, +, 2, 1);

#undef HEX_EDGE_FCN


      // --------------------
      //  c) face functions
      // --------------------

      Double val_1, deriv_1, val_2, deriv_2, val_3, deriv_3;
      Double sFlag2, sFlag3;
      Integer order_1, order_2, order_3;
      UInt d2, d3;

#define HEX_FACE_FCN(faceNum, sgn1, dir1, dir2, dir3)                   \
      order_2 = legFct->GetMaxOrderLocDir(dir2);                        \
      order_3 = legFct->GetMaxOrderLocDir(dir3);                        \
      if( elem->faceFlags[faceNum-1].test(0) == true) {                 \
        sFlag2 = (elem->faceFlags[faceNum-1].                           \
                  test(2) == true) ? 1.0 : -1.0;                        \
        sFlag3 = (elem->faceFlags[faceNum-1].                           \
                  test(1) == true) ? 1.0 : -1.0;                        \
        d2 = dir2;                                                      \
        d3 = dir3;                                                      \
      } else {                                                          \
        sFlag3 = (elem->faceFlags[faceNum-1].test(2) == true) ? 1.0 : -1.0; \
        sFlag2 = (elem->faceFlags[faceNum-1].test(1) == true) ? 1.0 : -1.0; \
        d3 = dir2;                                                      \
        d2 = dir3;                                                      \
      }                                                                 \
      for( Integer i = 2; i <= (order_2)-2; i++ ) {                     \
        for( Integer j = 2; j <= (order_3)-2; j++ ) {                   \
          if( (i + j) > std::max(order_2, order_3) ) { continue; }      \
          EvalPolynom( val_2, deriv_2, i, lCoeff_[i],                   \
                       sFlag2* actCoord[d2] );                          \
          EvalPolynom( val_3, deriv_3, j, lCoeff_[j],                   \
                       sFlag3*actCoord[d3] );                           \
          Shape[offset] = 0.5 * ( 1.0 sgn1 actCoord[dir1] )             \
                              * val_2  * val_3;                         \
          offset++;                                                     \
        } }

      // Face #1
      HEX_FACE_FCN(1, -, 2, 0, 1);

      // Face #2
      HEX_FACE_FCN(2, -, 1, 0, 2);

      // Face #3
      HEX_FACE_FCN(3, +, 0, 1, 2);

      // Face #4
      HEX_FACE_FCN(4, +, 1, 0 ,2);

      // Face #5
      HEX_FACE_FCN(5, -, 0, 1 ,2);

      // Face #6
      HEX_FACE_FCN(6, +, 2, 0 ,1);

#undef HEX_FACE_FCN



      // ----------------------
      //  d) bubble functions
      // ----------------------
      order_1 = legFct->GetMaxOrderLocDir(0);
      order_2 = legFct->GetMaxOrderLocDir(1);
      order_3 = legFct->GetMaxOrderLocDir(2);
      Integer temp, maxOrder;
      temp = std::max( order_1, order_2);
      maxOrder = std::max( temp, order_3);
      for( Integer i = 2; i <= (order_1) - 4; i++ ) {
        for( Integer j = 2; j <= (order_2) - 4; j++ ) {
          for( Integer k = 2; k <= (order_3) - 4;k++ ) {

            // Condition for trunk space: i+j in [2..max{order_i,order_j}]
            if ( (i+j+k) > maxOrder ) {continue;}

            EvalPolynom( val_1, deriv_1, i, lCoeff_[i], actCoord[0] );
            EvalPolynom( val_2, deriv_2, j, lCoeff_[j], actCoord[1] );
            EvalPolynom( val_3, deriv_3, k, lCoeff_[k], actCoord[2] );

            Shape[offset++] = val_1 * val_2 * val_3;
          }
        }
      }
    } else {
      EXCEPTION( "Approximation type not known" );
    }

  }

  void Hexa1FE::CalcLocalDerivShapeFnc( Matrix<Double> & LDeriv,
                                        const Vector<Double> & actCoord,
                                        const Elem* elem, UInt dof,
                                        AnsatzFct::FctEntityType fctType ) {

    if( actFct_->GetType() == AnsatzFct::LAGRANGE ||
        fctType == AnsatzFct::NODE ) {

      LDeriv.Resize(NumNodes_,Dim_);
      LDeriv.Init();

      for( UInt i=0; i<NumNodes_; i++) {
        LDeriv[i][0] = 0.125
          * LCornerCoords_[0][i]
          * (1 + LCornerCoords_[1][i] * actCoord[1] )
          * (1 + LCornerCoords_[2][i] * actCoord[2]);

        LDeriv[i][1] = 0.125
          * (1 + LCornerCoords_[0][i] * actCoord[0] )
          * LCornerCoords_[1][i]
          * (1 + LCornerCoords_[2][i] * actCoord[2]);

        LDeriv[i][2] = 0.125
          * (1 + LCornerCoords_[0][i] * actCoord[0])
          * (1 + LCornerCoords_[1][i] * actCoord[1])
          * LCornerCoords_[2][i];
      }

    }else if(  actFct_->GetType() == AnsatzFct::SPECTRAL) {
      // ===============
      //  Spectral PART
      // ===============
      CalcSpectralDerivFct(LDeriv,actCoord,elem,dof,fctType);
      return;

    } else if ( actFct_->GetType() == AnsatzFct::LEGENDRE ) {

      // ===============
      //  LEGENDRE PART
      // ===============

      // Get number of ansatz functions
      shared_ptr<LegendreFct> legFct =
        dynamic_pointer_cast<LegendreFct, AnsatzFct>(actFct_);
      UInt totalFcns = GetNumFncs( actFct_ );
      LDeriv.Resize( totalFcns, Dim_ );
      LDeriv.Init();

      // --------------------
      //  a) nodal functions
      // --------------------
      UInt offset = 0;

      // First of all, calculate all nodal function derivatives
      for( UInt iNode = 0; iNode < NumCorners_; iNode++,offset++ ) {
        LDeriv[iNode][0] = 0.125
          * LCornerCoords_[0][iNode]
          * (1 + LCornerCoords_[1][iNode] * actCoord[1] )
          * (1 + LCornerCoords_[2][iNode] * actCoord[2]);

        LDeriv[iNode][1] = 0.125
          * (1 + LCornerCoords_[0][iNode] * actCoord[0] )
          * LCornerCoords_[1][iNode]
          * (1 + LCornerCoords_[2][iNode] * actCoord[2]);

        LDeriv[iNode][2] = 0.125
          * (1 + LCornerCoords_[0][iNode] * actCoord[0])
          * (1 + LCornerCoords_[1][iNode] * actCoord[1])
          * LCornerCoords_[2][iNode];
      }

      // -------------------
      //  b) edge functions
      // -------------------
      // Obtain order of element
      Integer order_1, order_2, order_3;
      Double val, deriv, factor;

#define HEX_EDGE_DERIV( edgeNum, sgn_1, dir_1, sgn_2, dir_2, dir_3)     \
      order_3 = legFct->GetMaxOrderLocDir(dir_3);                       \
      factor = elem->edges[edgeNum-1] < 0 ? -1.0 : 1.0;                 \
      for( Integer iDof = 2; iDof <= order_3; iDof++, offset++ ) {      \
        EvalPolynom( val, deriv, iDof, lCoeff_[iDof],                   \
                     factor*actCoord[dir_3] );                          \
        LDeriv[offset][dir_1] = 0.25 * sgn_1(1.0 sgn_2 actCoord[dir_2]) \
          * val;                                                        \
        LDeriv[offset][dir_2] = 0.25 * sgn_2(1.0 sgn_1 actCoord[dir_1]) \
          * val;                                                        \
        LDeriv[offset][dir_3] = 0.25 * (1.0 sgn_1 actCoord[dir_1])      \
          * (1.0 sgn_2 actCoord[dir_2])                                 \
          * deriv * factor;                                             \
      }

      // EDGE #1
      HEX_EDGE_DERIV(1, -, 1, -, 2, 0 );

      // EDGE #2
       HEX_EDGE_DERIV(2, +, 0, -, 2, 1 );

       // EDGE #3
       HEX_EDGE_DERIV(3, +, 1, -, 2, 0 );

       // EDGE #4
       HEX_EDGE_DERIV(4, -, 0, -, 2, 1 );

       // EDGE #5
       HEX_EDGE_DERIV(5, -, 0, -, 1, 2 );

       // EDGE #6
       HEX_EDGE_DERIV(6, +, 0, -, 1, 2 );

       // EDGE #7
       HEX_EDGE_DERIV(7, +, 0, +, 1, 2 );

       // EDGE #8
       HEX_EDGE_DERIV(8, -, 0, +, 1, 2 );

       // EDGE #9
       HEX_EDGE_DERIV(9, -, 1, +, 2, 0 );

       // EDGE #10
       HEX_EDGE_DERIV(10, +, 0, +, 2, 1 );

       // EDGE #11
       HEX_EDGE_DERIV(11, +, 1, +, 2, 0 );

       // EDGE #12
       HEX_EDGE_DERIV(12, -, 0, +, 2, 1 );




      // -------------------
      //  c) face functions
      // -------------------

       Double val_2, deriv_2, val_3, deriv_3;
      Double sFlag2, sFlag3;
      UInt d2, d3;

#define HEX_FACE_DERIV(faceNum, sgn_1, dir_1, dir_2, dir_3)             \
      order_2 = legFct->GetMaxOrderLocDir(dir_2);                       \
      order_3 = legFct->GetMaxOrderLocDir(dir_3);                       \
      LOG_DBG3(hexa1fe) << "Treating face " << elem->faces[faceNum-1];  \
      LOG_DBG3(hexa1fe) << "bitset: " << elem->faceFlags[faceNum-1];    \
      if( elem->faceFlags[faceNum-1].test(0) == true) {                 \
        sFlag2 = (elem->faceFlags[faceNum-1].                           \
                  test(2) == true) ? 1.0 : -1.0;                        \
        sFlag3 = (elem->faceFlags[faceNum-1].                           \
                  test(1) == true) ? 1.0 : -1.0;                        \
        d2 = dir_2;                                                     \
        d3 = dir_3;                                                     \
      } else {                                                          \
        sFlag3 = (elem->faceFlags[faceNum-1].test(2) == true) ? 1.0 : -1.0; \
        sFlag2 = (elem->faceFlags[faceNum-1].test(1) == true) ? 1.0 : -1.0; \
        d3 = dir_2;                                                     \
        d2 = dir_3;                                                     \
                                                                        \
      }                                                                 \
      LOG_DBG3(hexa1fe) << "I-Flag: " << sFlag2                         \
                        << ", II-Flag: " << sFlag3;                     \
      LOG_DBG3(hexa1fe) << "I-direction: " << d2                        \
                        << ", II-direction: " << d3;                    \
      for( Integer i = 2; i <= (order_2) - 2; i++ ) {                   \
        for( Integer j = 2; j <= (order_3) - 2; j++ ) {                 \
          if( (i + j) > std::max(order_2,order_3) ) { continue; }       \
          EvalPolynom( val_2, deriv_2, i, lCoeff_[i],                   \
                       sFlag2* actCoord[d2] );                          \
          EvalPolynom( val_3, deriv_3, j, lCoeff_[j],                   \
                       sFlag3*actCoord[d3] );                           \
          LOG_DBG3(hexa1fe) << "deriv_2 = " << deriv_2                  \
                            << ", deriv_3 = " << deriv_3;               \
          LOG_DBG3(hexa1fe) << "val_2 = " << val_2                      \
                            << ", val_3 = " << val_3 << std::endl;      \
          LDeriv[offset][dir_1] = sgn_1(0.5 * val_2 * val_3);           \
          LDeriv[offset][d2] = 0.5 * ( 1.0 sgn_1 actCoord[dir_1])       \
            * deriv_2 * sFlag2 * val_3;                                 \
          LDeriv[offset][d3] = 0.5 * ( 1.0 sgn_1 actCoord[dir_1] )      \
            * val_2 * deriv_3 * sFlag3;                                 \
          offset++;                                                     \
        } }

      // Face #1
      HEX_FACE_DERIV(1, -, 2, 0, 1);

      // Face #2
      HEX_FACE_DERIV(2, -, 1, 0, 2);

      // Face #3
      HEX_FACE_DERIV(3, +, 0, 1, 2);

      // Face #4
      HEX_FACE_DERIV(4, +, 1, 0, 2);

      // Face #5
      HEX_FACE_DERIV(5, -, 0, 1, 2);

      // Face #6
      HEX_FACE_DERIV(6, +, 2, 0, 1);

      // ---------------------
      //  d) bubble functions
      // ---------------------
      order_1 = legFct->GetMaxOrderLocDir(0);
      order_2 = legFct->GetMaxOrderLocDir(1);
      order_3 = legFct->GetMaxOrderLocDir(2);
      Integer temp, maxOrder;
      temp = std::max( order_1, order_2);
      maxOrder = std::max( temp, order_3);
      Double val_i, deriv_i, val_j, deriv_j, val_k, deriv_k;
      for( Integer i = 2; i <= (order_1) - 4; i++ ) {
        for( Integer j = 2; j <= (order_2) - 4; j++ ) {
          for( Integer k = 2; k <= (order_3) - 4; k++ ) {

            // Condition for trunk space:
            // i+j+k in [6..max{order_i,order_j,order_k}]
            if ( (i+j+k) > maxOrder ) {continue;}
            EvalPolynom( val_i, deriv_i, i, lCoeff_[i], actCoord[0] );
            EvalPolynom( val_j, deriv_j, j, lCoeff_[j], actCoord[1] );
            EvalPolynom( val_k, deriv_k, k, lCoeff_[k], actCoord[2] );

            LDeriv[offset][0] = deriv_i * val_j   * val_k;
            LDeriv[offset][1] =   val_i * deriv_j * val_k;
            LDeriv[offset][2] =   val_i * val_j   * deriv_k;
            offset++;
          }
        }
      }
      //LOG_DBG2(hexa1fe) << "LDeriv = \n" << LDeriv << std::endl;
      //LOG_DBG3(hexa1fe) << "offset after bubbles: " << offset;
    } else {
      EXCEPTION( "Approximation type not known" );
    }
  }


  void Hexa1FE::CalcLocalICModesDerivShapeFnc( Matrix<Double> & LDeriv,
					       const Vector<Double> & actCoord,
					       const Elem* elem, UInt dof,
					       AnsatzFct::FctEntityType fctType ) {

    if( actFct_->GetType() == AnsatzFct::LAGRANGE ||
        fctType == AnsatzFct::NODE ) {

      LDeriv.Resize(3,Dim_);
      LDeriv.Init();
      LDeriv[0][0] = -2.0*actCoord[0];
      LDeriv[1][1] = -2.0*actCoord[1];
      LDeriv[2][2] = -2.0*actCoord[2];
    }
    else if ( actFct_->GetType() == AnsatzFct::LEGENDRE ) {
      EXCEPTION( "CalcLocalICModesDerivShapeFnc for Legendre type not implemented");
    }

  }


  void  Hexa1FE::GetNumFncs( StdVector<UInt>& numFcns,
                             const shared_ptr<AnsatzFct>& fcnType,
                             AnsatzFct::FctEntityType fctEntityType,
                             UInt dof ) {

    // Check ansatzFctType
    if( fcnType->GetType() == AnsatzFct::LAGRANGE ) {
      numFcns.Resize( NumNodes_ );
      numFcns.Init(1);
    } else if( fcnType->GetType() == AnsatzFct::NEDELEC ) {
      numFcns.Resize( NumEdges_ );
      numFcns.Init(1);
    } else if ( fcnType->GetType() == AnsatzFct::SPECTRAL ) {
       Integer order =  dynamic_pointer_cast<SpectralFct, AnsatzFct> (fcnType)->GetOrder();
      if( fctEntityType == AnsatzFct::NODE ) {
        numFcns.Resize( NumNodes_ );
        numFcns.Init( 1 );

      } else if( fctEntityType == AnsatzFct::EDGE ) {
        numFcns.Resize( NumEdges_ );
        numFcns.Init( order - 1 );
      }
      else if( fctEntityType == AnsatzFct::INTERIOR ) {
        numFcns.Resize(1);
        numFcns.Init((order-1)*(order-1)*(order-1));
      }
      else if( fctEntityType == AnsatzFct::FACE ) {
        numFcns.Resize( NumFaces_ );
        numFcns.Init((order-1)*(order-1));
      } else {
        EXCEPTION( "Not yet implemented!" );
      }

    } else if ( fcnType->GetType() == AnsatzFct::LEGENDRE ) {

      if( fcnType->IsIsotropic() == true ) {

        // *** isotropic case ***
        // Remember approximation order
        Integer order =  dynamic_pointer_cast<LegendreFct, AnsatzFct>
          (fcnType)->GetIsoOrder();

        // Check for subentity-type
        if( fctEntityType == AnsatzFct::NODE ) {
          numFcns.Resize( NumNodes_ );
          numFcns.Init( 1 );

        } else if( fctEntityType == AnsatzFct::EDGE ) {
          numFcns.Resize( NumEdges_ );
          numFcns.Init( order - 1 );
        } else if( fctEntityType == AnsatzFct::FACE ) {
          numFcns.Resize( NumFaces_ );

          // faces
          UInt numFaceFncs = 0;
          for( Integer i = 2; i<= order-2; i++ ) {
            for( Integer j = 2; j<= order-2; j++ ) {
              if( i+j <= order) {
                numFaceFncs++;
              }
            }
          }
          numFcns.Init( numFaceFncs );
        }
        else if( fctEntityType == AnsatzFct::INTERIOR ) {
          numFcns.Resize(1);
          UInt total = 0;
          for( Integer i = 2; i<= order-4; i++ ) {
            for( Integer j = 2; j<= order-4; j++ ) {
              for( Integer k = 2; k<= order-4; k++ ) {
                if( (i+j+k) <= order) {
                  total++;
                }
              }
            }
          }

          numFcns.Init(total);
        } else {
          EXCEPTION("Not yet implemented!");
        }

      } else {
        // *** anisotropic case ***
        // Remember approximation order
        shared_ptr<LegendreFct> const & legFct =
          dynamic_pointer_cast<LegendreFct, AnsatzFct> (fcnType);
        Matrix<UInt> const & order =  legFct->GetAnisoOrder();

        // Check for subentity-type
        if( fctEntityType == AnsatzFct::NODE ) {
          numFcns.Resize( NumNodes_ );
          numFcns.Init( 1 );

        } else if( fctEntityType == AnsatzFct::EDGE ) {
          numFcns.Resize( NumEdges_ );
          numFcns.Init(0);
          // xi-direction
          numFcns[0]  = order[0][dof]-1;
          numFcns[2]  = numFcns[0];
          numFcns[8]  = numFcns[0];
          numFcns[10] = numFcns[0];

          // eta-direction
          numFcns[1]  = order[1][dof]-1;
          numFcns[3]  = numFcns[1];
          numFcns[9]  = numFcns[1];
          numFcns[11] = numFcns[1];

          // zeta direction
          numFcns[4]  = order[2][dof]-1;
          numFcns[5]  = numFcns[4];
          numFcns[6]  = numFcns[4];
          numFcns[7]  = numFcns[4];
        } else if( fctEntityType == AnsatzFct::FACE ) {
          numFcns.Resize( NumFaces_ );

          UInt numFaceFncs = 0;
          Integer max = 0;
          // Face #1 and #6
          max = std::max( order[0][dof], order[1][dof] );
          for( Integer i = 2; i<= (Integer)order[0][dof]-2; i++ ) {
            for( Integer j = 2; j<= (Integer)order[1][dof]-2; j++ ) {
              if( i+j <= max) {
                numFaceFncs++;
              }
            }
          }
          numFcns[0] = numFaceFncs;
          numFcns[5] = numFaceFncs;

          // Face #2 and #4
          numFaceFncs = 0;
          max = std::max( order[0][dof], order[2][dof] );
          for( Integer i = 2; i<= (Integer)order[0][dof]-2; i++ ) {
            for( Integer j = 2; j<= (Integer)order[2][dof]-2; j++ ) {
              if( i+j <= max) {
                numFaceFncs++;
              }
            }
          }
          numFcns[1] = numFaceFncs;
          numFcns[3] = numFaceFncs;

          // Face #3 and #5
          numFaceFncs = 0;
          max = std::max( order[1][dof], order[2][dof] );
          for( Integer i = 2; i<= (Integer)order[1][dof]-2; i++ ) {
            for( Integer j = 2; j<= (Integer)order[2][dof]-2; j++ ) {
              if( i+j <= max) {
                numFaceFncs++;
              }
            }
          }
          numFcns[2] = numFaceFncs;
          numFcns[4] = numFaceFncs;

        }
        else if( fctEntityType == AnsatzFct::INTERIOR ) {

          numFcns.Resize(1);
          numFcns[0] = 0;
          UInt total = 0;
          UInt max = legFct->GetMaxOrderDof( dof );
          for( Integer i = 2; i<= (Integer)order[0][dof]-4; i++ ) {
            for( Integer j = 2; j<= (Integer)order[1][dof]-4; j++ ) {
              for( Integer k = 2; k<= (Integer)order[2][dof]-4; k++ ) {
                if( (i+j+k) <= (Integer) max) {
                  total++;
                }
              }
            }
          }

          numFcns.Init(total);
        } else {
          EXCEPTION("Not yet implemented!");
        }
      }
    } else {
      EXCEPTION("AnsatzFcnType '" << fcnType->GetType()
             << "' is not known!");
    }

  }

  UInt Hexa1FE::GetNumFncs( const shared_ptr<AnsatzFct>& type ) {

    // Check ansatzFctType
    if( type->GetType() == AnsatzFct::LAGRANGE ) {
      return NumNodes_;

    } else if ( type->GetType() == AnsatzFct::SPECTRAL ) {
      Integer order =  dynamic_pointer_cast<SpectralFct, AnsatzFct> (type)->GetOrder();
      return (order+1)*(order+1)*(order+1);

    } else if ( type->GetType() == AnsatzFct::LEGENDRE ) {

      // Cache check: if ansatz function type is the same as
      // previously calculated, simply return the last calculated number)
   //    if ( type == actFct_ ) {
//         return actNumFcns_;
//       }

      if( type->IsIsotropic() == true ) {

        Integer order =  dynamic_pointer_cast<LegendreFct, AnsatzFct>
          (type)->GetIsoOrder();

        // edge functions
        UInt numEdgeFncs = NumEdges_* (order-1);

        // faces
        UInt numFaceFncs = 0;
        for( Integer i = 2; i<= order-2; i++ ) {
          for( Integer j = 2; j<= order-2; j++ ) {
            if( i+j <= order) {
              numFaceFncs++;
            }
          }
        }
        numFaceFncs *= NumFaces_;

        LOG_DBG3(hexa1fe) << "numFaceFncs = " << numFaceFncs << std::endl;
        // determine number of bubble functions
        UInt numBubbles = 0;
        for( Integer i = 2; i<= order-4; i++ ) {
          for( Integer j = 2; j<= order-4; j++ ) {
            for( Integer k = 2; k<= order-4; k++ ) {
              if( i+j+k <= order) {
                numBubbles++;
              }
            }
          }
        }

        LOG_DBG3(hexa1fe) << "numBubbles = " << numBubbles << std::endl;
        LOG_DBG3(hexa1fe) << "total number of unknowns: "
                          <<  NumNodes_ + numFaceFncs
          + numEdgeFncs + numBubbles
                          << std::endl;

        actNumFcns_ = (NumNodes_ + numFaceFncs + numEdgeFncs + numBubbles);
        return actNumFcns_;

      } else {
        // *** anisotropic case ***
        Integer max_1, max_2, max_3;
        shared_ptr<LegendreFct> legFct =
          dynamic_pointer_cast<LegendreFct, AnsatzFct>(type);

        max_1 = legFct->GetMaxOrderLocDir( 0 );
        max_2 = legFct->GetMaxOrderLocDir( 1 );
        max_3 = legFct->GetMaxOrderLocDir( 2 );

        // a) nodes
        UInt numNodeModes = 8;

        // b) edges
        UInt numEdgeModes = 0;
        numEdgeModes += max_1 > 0 ? (max_1-1)*4 : 0;
        numEdgeModes += max_2 > 0 ? (max_2-1)*4 : 0;
        numEdgeModes += max_3 > 0 ? (max_3-1)*4 : 0;

        Matrix<UInt> const & order = legFct->GetAnisoOrder();

        // c) faces
        UInt numFaceModes = 0;
        Integer max = 0;
        // Face #1 and #6
        max = std::max( max_1, max_2 );
        for( Integer i = 2; i<= max_1 - 2; i++ ) {
          for( Integer j = 2; j<= max_2 - 2; j++ ) {
            if( i+j <= max) {
              numFaceModes++;
            }
          }
        }
        // Face #2 and #4
        max = std::max( max_1, max_3 );
        for( Integer i = 2; i<= max_1 - 2; i++ ) {
          for( Integer j = 2; j<= max_3 - 2; j++ ) {
            if( i+j <= max) {
              numFaceModes++;

            }
            }
        }

        // Face #3 and #5
        max = std::max( max_2, max_3  );
        for( Integer i = 2; i<= max_2 - 2 ; i++ ) {
          for( Integer j = 2; j<= max_3 - 2; j++ ) {
            if( i+j <= max) {
              numFaceModes++;
            }
          }
        }


        numFaceModes *= 2;

        // d) internals
        UInt numBubbleModes = 0;
        for( UInt iDof = 0; iDof < 3; iDof++ ) {
          max = legFct->GetMaxOrderDof( iDof );
          for( Integer i = 2; i<= (Integer)order[0][iDof]-4; i++ ) {
            for( Integer j = 2; j<= (Integer)order[1][iDof]-4; j++ ) {
              for( Integer k = 2; k<= (Integer)order[2][iDof]-4; k++ ) {
                if( (i+j+k) <= max) {
                  numBubbleModes++;
                }
              }
            }
          }
        }


        actNumFcns_ = numNodeModes + numEdgeModes
          + numFaceModes + numBubbleModes;
        return actNumFcns_;
      }
    }
    return 0;
  }

  void Hexa1FE::SetAnsatzFct( shared_ptr<AnsatzFct>& actFct,
                              bool setIntPoints ) {
  // Check if this ansatz fct was already set
    if( actFct_ == actFct ) {
      return;
    }

    // first of all: compute the number of ansatz functions!
    actNumFcns_ = GetNumFncs( actFct );

    actFct_ = actFct;
   if( actFct->GetType() == AnsatzFct::SPECTRAL ) {
      // If not, get order of functions
      shared_ptr<SpectralFct> lagFct =
        dynamic_pointer_cast<SpectralFct, AnsatzFct>(actFct);

      //check for higher order lagrange functions
      //right now: higher order lagrange => Spectral Elements
      IntegMethod = LOBATTO;
      IntegOrder = 2 * lagFct->GetOrder() - 1;

      // get the values by IntegMethod and IntegOrder
      SetIntPoints();

      // ... then calc shape function values at integration points
      // for subsequent calls to CalcJacobian() ....
      SetShapeFncAtIp();
      SetShapeFncDerivAtIp();

    } else if( actFct->GetType() == AnsatzFct::LEGENDRE
        && setIntPoints == true ) {

      // If not, get order of functions
      shared_ptr<LegendreFct> legFct =
        dynamic_pointer_cast<LegendreFct, AnsatzFct>(actFct);


      // check if isotropic order is present
     //  if( legFct->IsIsotropic() ) {
//         IntegMethod = ECONOMICAL;
//         // Prevent integration of first order, as this may
//         // cause non-reasonable results
//         if( legFct->GetIsoOrder() > 1 ) {
//           IntegOrder =  legFct->GetIsoOrder()*3;
//           if( IntegOrder < 5 )
//             IntegOrder = 5;
//           if ( IntegOrder > 18 )
//             IntegOrder = 19;

//         } else {
//           IntegOrder = 5;
//         }
//       } else {
        UInt max_1, max_2, max_3;
        IntegMethod = CARTESIAN;
        max_1 = legFct->GetMaxOrderLocDir(0);
        max_2 = legFct->GetMaxOrderLocDir(1);
        max_3 = legFct->GetMaxOrderLocDir(2);

        IntegOrder = EncodeCartesianOrder( max_1*2, max_2*2, max_3*2 );
        //std::cerr << "integration order: " << IntegOrder << std::endl;
    //   }

      // get the values by IntegMethod and IntegOrder
      SetIntPoints();

      // ... then calc shape function values at integration points
      // for subsequent calls to CalcJacobian() ....
      SetShapeFncAtIp();
      SetShapeFncDerivAtIp();
    }

  }

  void Hexa1FE::CalcSpectralShFct( Vector<Double> & Shape,
                                  const Vector<Double> & LCoord,
                                  const Elem* elem, UInt dof,
                                  AnsatzFct::FctEntityType type ){
    shared_ptr<SpectralFct> myFct = dynamic_pointer_cast<SpectralFct, AnsatzFct>(actFct_);
    UInt order = myFct->GetOrder();

    Shape.Resize( (order+1) * (order+1) * (order+1) );
    Shape.Init();

    //! 1D Lagrange functions at IPs for spectral mode
    Vector<Double> * sShFcnAtIp_ = new Vector<Double>[3];
    sShFcnAtIp_[0].Init();
    sShFcnAtIp_[1].Init();
    sShFcnAtIp_[2].Init();
    myFct->EvaluatePolynomial( sShFcnAtIp_[0], LCoord[0] );
    myFct->EvaluatePolynomial( sShFcnAtIp_[1], LCoord[1] );
    myFct->EvaluatePolynomial( sShFcnAtIp_[2], LCoord[2] );

    UInt c = 0;

    Shape[c++]= sShFcnAtIp_[0][0]      *  sShFcnAtIp_[1][0]       * sShFcnAtIp_[2][0];
    Shape[c++]= sShFcnAtIp_[0][order] *  sShFcnAtIp_[1][0]       * sShFcnAtIp_[2][0];
    Shape[c++]= sShFcnAtIp_[0][order] *  sShFcnAtIp_[1][order]  * sShFcnAtIp_[2][0];
    Shape[c++]= sShFcnAtIp_[0][0]      *  sShFcnAtIp_[1][order]  * sShFcnAtIp_[2][0];

    Shape[c++]= sShFcnAtIp_[0][0]      * sShFcnAtIp_[1][0]        * sShFcnAtIp_[2][order];
    Shape[c++]= sShFcnAtIp_[0][order] * sShFcnAtIp_[1][0]        * sShFcnAtIp_[2][order];
    Shape[c++]= sShFcnAtIp_[0][order] * sShFcnAtIp_[1][order]   * sShFcnAtIp_[2][order];
    Shape[c++]= sShFcnAtIp_[0][0]      * sShFcnAtIp_[1][order]   * sShFcnAtIp_[2][order];


    // --------------------
    //  b) edge functions
    // --------------------
#define FILL_EDGE( numEdge, shape, p1, p2)      \
    {                                          \
      if(elem->edges[numEdge-1] < 0){          \
        for ( UInt i= order-1; i>0 ;i-- )    \
          Shape[c++] = shape[i] * p1 * p2;    \
      }else                                    \
      {                                        \
        for ( UInt i= 1; i< order ;i++ )      \
          Shape[c++] = shape[i] * p1 * p2;    \
      }                                        \
    }

    FILL_EDGE( 1, sShFcnAtIp_[0], sShFcnAtIp_[1][0],      sShFcnAtIp_[2][0] );
    FILL_EDGE( 2, sShFcnAtIp_[1], sShFcnAtIp_[0][order], sShFcnAtIp_[2][0] );
    FILL_EDGE( 3, sShFcnAtIp_[0], sShFcnAtIp_[1][order], sShFcnAtIp_[2][0] );
    FILL_EDGE( 4, sShFcnAtIp_[1], sShFcnAtIp_[0][0],      sShFcnAtIp_[2][0] );
    FILL_EDGE( 5, sShFcnAtIp_[2], sShFcnAtIp_[0][0],      sShFcnAtIp_[1][0] );
    FILL_EDGE( 6, sShFcnAtIp_[2], sShFcnAtIp_[0][order], sShFcnAtIp_[1][0] );
    FILL_EDGE( 7, sShFcnAtIp_[2], sShFcnAtIp_[0][order], sShFcnAtIp_[1][order] );
    FILL_EDGE( 8, sShFcnAtIp_[2], sShFcnAtIp_[0][0],      sShFcnAtIp_[1][order] );
    FILL_EDGE( 9, sShFcnAtIp_[0], sShFcnAtIp_[1][0],      sShFcnAtIp_[2][order] );
    FILL_EDGE(10, sShFcnAtIp_[1], sShFcnAtIp_[0][order], sShFcnAtIp_[2][order] );
    FILL_EDGE(11, sShFcnAtIp_[0], sShFcnAtIp_[1][order], sShFcnAtIp_[2][order] );
    FILL_EDGE(12, sShFcnAtIp_[1], sShFcnAtIp_[0][0],      sShFcnAtIp_[2][order] );

    // --------------------
    //  c) face functions
    // --------------------

#define FILL_FACE( numFace, shape1, shape2, fac)                                  \
      {                                                                            \
        Integer s1,s2;                                                                \
        if(! elem->faceFlags[numFace-1].test(0)){                                  \
          std::swap(shape1,shape2);                                                \
          s1 = (elem->faceFlags[numFace-1].test(1)) ? 0 : order;                    \
          s2 = (elem->faceFlags[numFace-1].test(2)) ? 0 : order;                    \
        }else {                                                                    \
          s1 = (elem->faceFlags[numFace-1].test(2)) ? 0 : order;                    \
          s2 = (elem->faceFlags[numFace-1].test(1)) ? 0 : order;                    \
        }                                                                           \
                                                                                  \
        for ( Integer j= 1;j< (Integer)order ;j++ )                              \
        {                                                                          \
          for ( Integer i = 1; i< (Integer)order ;i++ )                              \
          {                                                                        \
            Shape[c++] = shape1[abs(s1-i)] * shape2[abs(s2-j)] * fac;              \
          }                                                                        \
        }                                                                          \
        if(! elem->faceFlags[numFace-1].test(0))                                  \
          std::swap(shape1,shape2);                                                \
      }

    FILL_FACE( 1, sShFcnAtIp_[0], sShFcnAtIp_[1], sShFcnAtIp_[2][0] );
    FILL_FACE( 2, sShFcnAtIp_[0], sShFcnAtIp_[2], sShFcnAtIp_[1][0] );
    FILL_FACE( 3, sShFcnAtIp_[2], sShFcnAtIp_[1], sShFcnAtIp_[0][order] );
    FILL_FACE( 4, sShFcnAtIp_[0], sShFcnAtIp_[2], sShFcnAtIp_[1][order] );
    FILL_FACE( 5, sShFcnAtIp_[2], sShFcnAtIp_[1], sShFcnAtIp_[0][0] );
    FILL_FACE( 6, sShFcnAtIp_[0], sShFcnAtIp_[1], sShFcnAtIp_[2][order] );


    // ----------------------
    //  d) bubble functions
    // ----------------------
    for( UInt i = 1; i < order ; i++ ) {
      for( UInt j = 1; j < order ; j++ ) {
        for( UInt k = 1; k < order ; k++ ) {
          Shape[c++] = sShFcnAtIp_[0][k] * sShFcnAtIp_[1][j] * sShFcnAtIp_[2][i];
        }
      }
    }
  }

  void Hexa1FE::CalcSpectralDerivFct( Matrix<Double> & LDeriv,
                                     const Vector<Double> & LCoord,
                                     const Elem* elem, UInt dof,
                                     AnsatzFct::FctEntityType type){
    shared_ptr<SpectralFct> myFct = dynamic_pointer_cast<SpectralFct, AnsatzFct>(actFct_);
    UInt order = myFct->GetOrder();

    UInt totalFcns = GetNumFncs(actFct_);
    LDeriv.Resize( totalFcns, Dim_ );
    LDeriv.Init();

    //! 1D Lagrange functions at IPs for spectral mode
    Vector<Double> * sShFcnAtIp_ = new Vector<Double>[3];
    Vector<Double> * sDerivAtIp_ = new Vector<Double>[3];
    sShFcnAtIp_[0].Init();
    sShFcnAtIp_[1].Init();
    sShFcnAtIp_[2].Init();
    sDerivAtIp_[0].Init();
    sDerivAtIp_[1].Init();
    sDerivAtIp_[2].Init();

    LDeriv.Resize( (order+1) * (order+1) * (order+1), Dim_ );

    myFct->EvaluatePolynomial( sShFcnAtIp_[0], LCoord[0] );
    myFct->EvaluatePolynomial( sShFcnAtIp_[1], LCoord[1] );
    myFct->EvaluatePolynomial( sShFcnAtIp_[2], LCoord[2] );
    myFct->EvaluateDerivPolynomial( sDerivAtIp_[0], LCoord[0] );
    myFct->EvaluateDerivPolynomial( sDerivAtIp_[1], LCoord[1] );
    myFct->EvaluateDerivPolynomial( sDerivAtIp_[2], LCoord[2] );

    UInt c = 0;
    Integer swamping1 = 0;
    Integer swamping2 = 0;
    for ( UInt k = 0; k<= order ; k+=order)
    {
      for(Integer i = 0; i<= 3 ; i++)
      {
        LDeriv[c][0]=   sDerivAtIp_[0][swamping1] * sShFcnAtIp_[1][swamping2] * sShFcnAtIp_[2][k];
        LDeriv[c][1]=   sShFcnAtIp_[0][swamping1] * sDerivAtIp_[1][swamping2] * sShFcnAtIp_[2][k];
        LDeriv[c++][2]= sShFcnAtIp_[0][swamping1] * sShFcnAtIp_[1][swamping2] * sDerivAtIp_[2][k];
        swamping2 = swamping1;
        swamping1 += Integer(order) * Integer(i-1) * -1;
      }
      swamping1 = 0;
      swamping2 = 0;
    }

    // -------------------
    //  b) edge functions
    // -------------------
#define HEX_E_DERIV(numEdge,idx1,idx2,idx3,dim)                                 \
      {                                                                            \
        UInt start = 1;                                                              \
        Integer inc = 1;                                                            \
        UInt end = order;                                                          \
        if(elem->edges[numEdge-1] < 0){                                                \
          start = order-1;                                                            \
          inc = -1;                                                                  \
          end = 0;                                                                  \
        }                                                                            \
        switch(dim){                                                                \
        case 1:                                                                      \
          while( start != end)                                                      \
          {                                                                          \
            LDeriv[c][0]=   sDerivAtIp_[0][start] * sShFcnAtIp_[1][idx2] * sShFcnAtIp_[2][idx3];          \
            LDeriv[c][1]=   sShFcnAtIp_[0][start] * sDerivAtIp_[1][idx2] * sShFcnAtIp_[2][idx3];          \
            LDeriv[c++][2]= sShFcnAtIp_[0][start] * sShFcnAtIp_[1][idx2] * sDerivAtIp_[2][idx3];        \
            start += inc;                                                            \
          }                                                                          \
          break;                                                                    \
        case 2:                                                                      \
          while( start != end)                                                      \
          {                                                                          \
            LDeriv[c][0]=   sDerivAtIp_[0][idx1] * sShFcnAtIp_[1][start] * sShFcnAtIp_[2][idx3];          \
            LDeriv[c][1]=   sShFcnAtIp_[0][idx1] * sDerivAtIp_[1][start] * sShFcnAtIp_[2][idx3];          \
            LDeriv[c++][2]= sShFcnAtIp_[0][idx1] * sShFcnAtIp_[1][start] * sDerivAtIp_[2][idx3];        \
            start += inc;                                                            \
          }                                                                          \
          break;                                                                    \
        case 3:                                                                      \
          while( start != end)                                                      \
          {                                                                          \
            LDeriv[c][0]  = sDerivAtIp_[0][idx1] * sShFcnAtIp_[1][idx2] * sShFcnAtIp_[2][start];          \
            LDeriv[c][1]  = sShFcnAtIp_[0][idx1] * sDerivAtIp_[1][idx2] * sShFcnAtIp_[2][start];          \
            LDeriv[c++][2]= sShFcnAtIp_[0][idx1] * sShFcnAtIp_[1][idx2] * sDerivAtIp_[2][start];        \
            start += inc;                                                            \
          }                                                                          \
          break;                                                                    \
        }                                                                            \
      }

    // EDGE #1
    HEX_E_DERIV(1, 0, 0, 0, 1);
    // EDGE #2
    HEX_E_DERIV(2, order, 0, 0, 2);
     // EDGE #3
    HEX_E_DERIV(3, 0, order, 0, 1 );
     // EDGE #4
    HEX_E_DERIV(4, 0, 0, 0, 2);
     // EDGE #5
    HEX_E_DERIV(5, 0, 0, 0, 3);
     // EDGE #6
    HEX_E_DERIV(6, order, 0, 0, 3);
     // EDGE #7
    HEX_E_DERIV(7, order, order, 0, 3);
     // EDGE #8
    HEX_E_DERIV(8, 0, order, 0, 3);
     // EDGE #9
    HEX_E_DERIV(9, 0, 0, order, 1);
     // EDGE #10
    HEX_E_DERIV(10, order, 0, order, 2);
     // EDGE #11
    HEX_E_DERIV(11, 0, order, order, 1);
    // EDGE #12
    HEX_E_DERIV(12, 0, 0, order, 2);

    // -------------------
    //  c) face functions
    // -------------------
#define HEX_F_DERIV( numFace, shape1, shape1Deriv, shape2, shape2Deriv,                     \
                        fac, facDeriv, dim)                                                    \
  {                                                                                            \
    Integer s1,s2;                                                                            \
                                                                                              \
    switch(dim){                                                                              \
      case 3:                                                                                  \
        if( elem->faceFlags[numFace-1].test(0))                                                  \
        {                                                                                        \
          s1 = (elem->faceFlags[numFace-1].test(2)) ? 0 : order;                                    \
          s2 = (elem->faceFlags[numFace-1].test(1)) ? 0 : order;                                    \
          for ( Integer j= 1;j< (Integer)order ;j++ )                                                      \
          {                                                                                      \
            for ( Integer i = 1; i< (Integer)order ;i++ )                                                  \
            {                                                                                    \
              LDeriv[c][0] = shape1Deriv[abs(s1-i)] * shape2[abs(s2-j)] * fac;                  \
              LDeriv[c][1] = shape1[abs(s1-i)] * shape2Deriv[abs(s2-j)] * fac;                  \
              LDeriv[c++][2] = shape1[abs(s1-i)] * shape2[abs(s2-j)] * facDeriv;                \
            }                                                                                    \
          }                                                                                      \
        }else{                                                                                  \
          s1 = (elem->faceFlags[numFace-1].test(1)) ? 0 : order;                                    \
          s2 = (elem->faceFlags[numFace-1].test(2)) ? 0 : order;                                    \
          for ( Integer i= 1;i< (Integer)order ;i++ )                                          \
          {                                                                                      \
            for ( Integer j = 1; j< (Integer)order ;j++ )                                      \
            {                                                                                    \
              LDeriv[c][0] = shape1Deriv[abs(s1-i)] * shape2[abs(s2-j)] * fac;                  \
              LDeriv[c][1] = shape1[abs(s1-i)] * shape2Deriv[abs(s2-j)] * fac;                  \
              LDeriv[c++][2] = shape1[abs(s1-i)] * shape2[abs(s2-j)] * facDeriv;                \
            }                                                                                    \
          }                                                                                      \
        }                                                                                        \
        break;                                                                                \
      case 2:                                                                                  \
        if( elem->faceFlags[numFace-1].test(0))                                                  \
        {                                                                                        \
          s1 = (elem->faceFlags[numFace-1].test(2)) ? 0 : order;                                    \
          s2 = (elem->faceFlags[numFace-1].test(1)) ? 0 : order;                                    \
          for ( Integer j= 1;j< (Integer)order ;j++ )                                                      \
          {                                                                                      \
            for ( Integer i = 1; i< (Integer)order ;i++ )                                                  \
            {                                                                                    \
              LDeriv[c][0] = shape1Deriv[abs(s1-i)] * shape2[abs(s2-j)] * fac;                  \
              LDeriv[c][1] = shape1[abs(s1-i)] * shape2[abs(s2-j)] * facDeriv;                  \
              LDeriv[c++][2] = shape1[abs(s1-i)] * shape2Deriv[abs(s2-j)] * fac;                \
            }                                                                                    \
          }                                                                                      \
        }else{                                                                                  \
          s1 = (elem->faceFlags[numFace-1].test(1)) ? 0 : order;                                    \
          s2 = (elem->faceFlags[numFace-1].test(2)) ? 0 : order;                                    \
          for ( Integer i= 1;i< (Integer)order ;i++ )                                                      \
          {                                                                                      \
            for ( Integer j = 1; j< (Integer)order ;j++ )                                                  \
            {                                                                                    \
              LDeriv[c][0] = shape1Deriv[abs(s1-i)] * shape2[abs(s2-j)] * fac;                  \
              LDeriv[c][1] = shape1[abs(s1-i)] * shape2[abs(s2-j)] * facDeriv;                  \
              LDeriv[c++][2] = shape1[abs(s1-i)] * shape2Deriv[abs(s2-j)] * fac;                \
            }                                                                                    \
          }                                                                                      \
        }                                                                                        \
        break;                                                                                \
      case 1:                                                                                  \
        if( elem->faceFlags[numFace-1].test(0))                                                  \
        {                                                                                        \
          s1 = (elem->faceFlags[numFace-1].test(2)) ? 0 : order;                                    \
          s2 = (elem->faceFlags[numFace-1].test(1)) ? 0 : order;                                    \
          for ( Integer j= 1;j< (Integer)order ;j++ )                                                      \
          {                                                                                      \
            for ( Integer i = 1; i< (Integer)order ;i++ )                                                  \
            {                                                                                    \
              LDeriv[c][0] = shape1[abs(s1-i)] * shape2[abs(s2-j)] * facDeriv;                  \
              LDeriv[c][1] = shape1[abs(s1-i)] * shape2Deriv[abs(s2-j)] * fac;                  \
              LDeriv[c++][2] = shape1Deriv[abs(s1-i)] * shape2[abs(s2-j)] * fac;                \
            }                                                                                    \
          }                                                                                      \
        }else{                                                                                  \
          s1 = (elem->faceFlags[numFace-1].test(1)) ? 0 : order;                                    \
          s2 = (elem->faceFlags[numFace-1].test(2)) ? 0 : order;                                    \
          for ( Integer i= 1;i< (Integer)order ;i++ )                                                      \
          {                                                                                      \
            for ( Integer j = 1; j< (Integer)order ;j++ )                                                  \
            {                                                                                    \
              LDeriv[c][0] = shape1[abs(s1-i)] * shape2[abs(s2-j)] * facDeriv;                  \
              LDeriv[c][1] = shape1[abs(s1-i)] * shape2Deriv[abs(s2-j)] * fac;                  \
              LDeriv[c++][2] = shape1Deriv[abs(s1-i)] * shape2[abs(s2-j)] * fac;                \
            }                                                                                    \
          }                                                                                      \
        }                                                                                        \
        break;                                                                                \
    }                                                                                          \
  }                                                                                            \

    HEX_F_DERIV( 1, sShFcnAtIp_[0], sDerivAtIp_[0], sShFcnAtIp_[1], sDerivAtIp_[1], sShFcnAtIp_[2][0],      sDerivAtIp_[2][0], 3 );
    HEX_F_DERIV( 2, sShFcnAtIp_[0], sDerivAtIp_[0], sShFcnAtIp_[2], sDerivAtIp_[2], sShFcnAtIp_[1][0],      sDerivAtIp_[1][0], 2 );
    HEX_F_DERIV( 3, sShFcnAtIp_[2], sDerivAtIp_[2], sShFcnAtIp_[1], sDerivAtIp_[1], sShFcnAtIp_[0][order], sDerivAtIp_[0][order], 1 );
    HEX_F_DERIV( 4, sShFcnAtIp_[0], sDerivAtIp_[0], sShFcnAtIp_[2], sDerivAtIp_[2], sShFcnAtIp_[1][order], sDerivAtIp_[1][order], 2 );
    HEX_F_DERIV( 5, sShFcnAtIp_[2], sDerivAtIp_[2], sShFcnAtIp_[1], sDerivAtIp_[1], sShFcnAtIp_[0][0],      sDerivAtIp_[0][0], 1 );
    HEX_F_DERIV( 6, sShFcnAtIp_[0], sDerivAtIp_[0], sShFcnAtIp_[1], sDerivAtIp_[1], sShFcnAtIp_[2][order], sDerivAtIp_[2][order], 3 );

    // ---------------------
    //  d) bubble functions
    // ---------------------
    for(UInt k = 1; k< order ; k++)
    {
      for(UInt j = 1; j< order ; j++)
      {
        for(UInt i = 1; i< order ; i++)
        {
          LDeriv[c][0]  = sDerivAtIp_[0][i] * sShFcnAtIp_[1][j] * sShFcnAtIp_[2][k];
          LDeriv[c][1]  = sShFcnAtIp_[0][i] * sDerivAtIp_[1][j] * sShFcnAtIp_[2][k];
          LDeriv[c++][2]= sShFcnAtIp_[0][i] * sShFcnAtIp_[1][j] * sDerivAtIp_[2][k];
        }
      }
    }
  }
  
  //=================================Edge Element==================================

   // see Ph.D. Sabine Zagelmayr
   void Hexa1FE :: CalcEdgeShapeFnc(Matrix<Double> & edgeShape, 
                                     const Vector<Double> & LCoord, 
                                     const Matrix<Double> & cornernodes,
                                     const Elem* elem)
   {

     edgeShape.Resize(NumEdges_, Dim_);
     edgeShape.Init();

     // nodal shape functions of a tet
     Vector<Double> nodeShape;
     CalcShapeFnc2(nodeShape, LCoord);


     // local derivates for special shape functions
     Matrix<Double> xDxi;  
     GetGlobDerivShFnc2(xDxi, LCoord, cornernodes, NULL, 1);    

     Double factor;
     for (UInt actEdge=0; actEdge<NumEdges_; actEdge++)
       {
         UInt node1 = edgeIndices_[actEdge][0] - 1;
         UInt node2 = edgeIndices_[actEdge][1] - 1;

         factor  = elem->edges[actEdge] < 0 ? -1.0 : 1.0;
         factor *= 0.5 * ( nodeShape[node1] + nodeShape[node2] );
         for (UInt actDim=0; actDim<Dim_; actDim++)
           for (UInt actDim=0; actDim<Dim_; actDim++)
             edgeShape[actEdge][actDim] = 
               ( xDxi[node2+NumNodes_][actDim] - xDxi[node1+NumNodes_][actDim] ) * factor;
       }  
   }


   // calculated the Nedelec shape function in an arbitrary point
   void Hexa1FE :: GetEdgeGlobalDerivShapeFnc(StdVector< Matrix<Double> > & shapeDeriv, 
                                               const Vector<Double> & lCoord,
                                               const Matrix<Double> & cornerCoords,
                                               const Elem* elem)
   {

     shapeDeriv.Resize(NumEdges_);

     Vector<Double> nodeShape;
     CalcShapeFnc2(nodeShape, lCoord);

     Matrix<Double> xDxi;  
     GetGlobDerivShFnc2(xDxi, lCoord, cornerCoords, NULL, 1);

     for (UInt actEdge=0; actEdge<NumEdges_; actEdge++) {
       shapeDeriv[actEdge].Resize(Dim_,Dim_);
       
       UInt node1 = edgeIndices_[actEdge][0] -1;
       UInt node2 = edgeIndices_[actEdge][1] -1;
       
       Double factor = elem->edges[actEdge] < 0 ? -1.0 : 1.0;  
       factor *= 0.5;

       for (UInt dim1=0; dim1<Dim_; dim1++)
         for (UInt dim2=0; dim2<Dim_; dim2++) {
           (shapeDeriv[actEdge])[dim2][dim1] = 
             ( ( xDxi[node1+NumNodes_][dim1] - xDxi[node2+NumNodes_][dim1] ) *
               ( xDxi[node1][dim2] + xDxi[node2][dim2] ) ) * factor;
         }
     }
   }


   void Hexa1FE :: GetGlobDerivShFnc2(Matrix<Double> & Deriv, 
                                       const Vector<Double> & LCoord,
                                       const Matrix<Double> & CornerCoords,
                                       const Elem * elem, 
                                       UInt dof )
   {

     CalcInvJacobian(JInv, LCoord, CornerCoords, elem );

     CalcLocalDerivShapeFnc2(LDeriv, LCoord, elem, dof); 

     Deriv = LDeriv * JInv;
   }


   void Hexa1FE :: SetShapeFnc2AtIp()
   {
   
     if (!ShFnc2AtIp_) {
       ShFnc2AtIp_ = new Vector<Double>[NumIntPoints_];
     } 
     else { 
       delete[] ShFnc2AtIp_ ;
       ShFnc2AtIp_ = new Vector<Double>[NumIntPoints_];
     }

     for( UInt i=0; i<NumIntPoints_; i++ ) {
       CalcShapeFnc2( ShFnc2AtIp_[i], IntPoints_[i]);
     }

   }

   void Hexa1FE :: CalcShapeFnc2(Vector<Double> & Shape, 
                                  const Vector<Double> & actCoord)
   {

     //"Hex Elements"
     // from Ph.D. thesis of Sabine Zagelmayer
     // 1-8: lambda_i, 9-16: sigma_i

     Shape.Resize(NumNodes_*2);
     Shape.Init();
     for( UInt i=0; i<NumNodes_; i++)
       Shape[i] = 0.125 
         * (1 + LCornerCoords_[0][i] * actCoord[0])
         * (1 + LCornerCoords_[1][i] * actCoord[1]) 
         * (1 + LCornerCoords_[2][i] * actCoord[2]);

     for( UInt i=0; i<NumNodes_; i++)
       Shape[i+NumNodes_] = 0.5 * (1 + LCornerCoords_[0][i] * actCoord[0])
         + 0.5 * (1 + LCornerCoords_[1][i] * actCoord[1]) 
         + 0.5 * (1 + LCornerCoords_[2][i] * actCoord[2]);
   }


   void Hexa1FE :: SetShapeFnc2DerivAtIp()
   {

     if( !ShFnc2DerivAtIp_) {
       ShFnc2DerivAtIp_ = new Matrix<Double>[NumIntPoints_]; 
     } 
     else{ 
       delete[] ShFnc2DerivAtIp_ ;
       ShFnc2DerivAtIp_ = new Matrix<Double>[NumIntPoints_];
     }

     for( UInt i=0; i<NumIntPoints_; i++ )
       CalcLocalDerivShapeFnc2( ShFnc2DerivAtIp_[i], IntPoints_[i], 
                                NULL, 1);
   }

   void Hexa1FE :: CalcLocalDerivShapeFnc2( Matrix<Double> & LDeriv, 
                                             const Vector<Double> & actCoord,
                                             const Elem*, UInt dof)
   {

       LDeriv.Resize(2*NumNodes_,Dim_);
       LDeriv.Init();
       
       for( UInt i=0; i<NumNodes_; i++) {
         LDeriv[i][0] = 0.125 
           * LCornerCoords_[0][i] 
           * (1 + LCornerCoords_[1][i] * actCoord[1] )
           * (1 + LCornerCoords_[2][i] * actCoord[2]);
         
         LDeriv[i][1] = 0.125 
           * (1 + LCornerCoords_[0][i] * actCoord[0] )
           * LCornerCoords_[1][i] 
           * (1 + LCornerCoords_[2][i] * actCoord[2]);
         
         LDeriv[i][2] = 0.125 
           * (1 + LCornerCoords_[0][i] * actCoord[0])
           * (1 + LCornerCoords_[1][i] * actCoord[1]) 
           * LCornerCoords_[2][i];
       }        

       for( UInt i=0; i<NumNodes_; i++) {
         LDeriv[i+NumNodes_][0] = 0.5 * LCornerCoords_[0][i];
         LDeriv[i+NumNodes_][1] = 0.5 * LCornerCoords_[1][i]; 
         LDeriv[i+NumNodes_][2] = 0.5 * LCornerCoords_[2][i];
       }   
   }
   void Hexa1FE::Global2LocalCoords(Matrix<Double> & localCoords,
                                    const Matrix<Double> & globalCoords,
                                    const Matrix<Double> & coordMat ){
     //BaseFE::Global2LocalCoords(localCoords,globalCoords,coordMat);
     //return;
     //Version 0 algorithm from duester script
     Vector<Double> delta_xi; // update for Newton-Raphson method
     Vector<Double> xi_start; // local start point for Newton-Raphson method
     Vector<Double> xi_k; // local point at iteration k
     Vector<Double> f; //current right hand side
     Vector<Double> f_start; //current right hand side
     Vector<Double> globalPoint; // global point coordinates
     //doubles for components of jacobian
     Double j11,j12,j13,j21,j22,j23,j31,j32,j33;
     UInt numPoints = globalCoords.GetNumCols(); // number of global points
     Double f_old; //storing the absolute value of search direction
     Double f_test; //storing the absolute value of search direction
     Double jac_det = 0;
     UInt k = 0; //iteration counter
     Integer l = 0; //stepping value
     Matrix<Double> J; // Jacobian at local point xi_k
     localCoords.Resize(3, numPoints);
     localCoords.Init();

     Double tolerance = 1e-14;

     //initialize everything
     xi_start.Resize(3);
     delta_xi.Resize(3);
     xi_k.Resize(3);
     J.Resize(3, 3); 
     globalPoint.Resize(3);

     //first initial guess is zero
     f.Resize(3);

    // Perform Newton-Raphson method on the list of global points
    // to find local coordinates within this element.
    for(UInt i = 0; i < numPoints; i++)
    {
      f.Init();
      globalPoint.Init();
      J.Init();
      xi_start.Init();
      xi_k.Init();
      k = 0;
      f_old = 222e20; 
      f_test = 0;
      for(UInt j = 0; j < 3; j++) {
        globalPoint[j] = globalCoords[j][i];
      }
      Local2GlobalCoord(f, xi_k, coordMat, NULL);
      f = f - globalPoint;
      f_old = f.NormL2();
      f_start = f;
      for(Double w=0.5; w <= 1.0; w+=0.5){
        for(UInt j=1;j<3;j++){
          for(UInt m=1;m<3;m++){
            for(UInt n=1;n<3;n++){
              xi_k[0] = pow(-1,j)*w;
              xi_k[1] = pow(-1,m)*w;
              xi_k[2] = pow(-1,n)*w;
              Local2GlobalCoord(f, xi_k, coordMat, NULL);
              f = f - globalPoint;
              f_test = f.NormL2();
              if(f_old > f_test){
                xi_start  = xi_k;
                f_start = f;
                f_old = f_test;
              }
            }
          }
        }
      }
      xi_k = xi_start;
      f = f_start;
      //std::cout << J << std::endl << std::endl;
      while(f_test > tolerance && k < 13){
         delta_xi.Init();
         xi_start.Init();
         // Calculate Jacobian at iteration point xi_k
         CalcJacobian(J, xi_k, coordMat, NULL );
         j11 = J[0][0]; j12 = J[0][1]; j13 = J[0][2];
         j21 = J[1][0]; j22 = J[1][1]; j23 = J[1][2];
         j31 = J[2][0]; j32 = J[2][1]; j33 = J[2][2];

         jac_det =   j11*j22*j33 - j11*j32*j23 - j21*j12*j33
                   + j32*j21*j13 - j22*j31*j13 + j31*j12*j23;
         delta_xi[0] = - j22*j33 * f[0] + j32*j23 * f[0] 
                       - j32*j13 * f[1] + j12*j33 * f[1]
                       - j12*j23 * f[2] + j22*j13 * f[2];
         delta_xi[1] = - j23*j31 * f[0] + j21*j33 * f[0] 
                       - j11*j33 * f[1] + j31*j13 * f[1]
                       - j21*j13 * f[2] + j23*j11 * f[2];
         delta_xi[2] = - j21*j32 * f[0] + j31*j22 * f[0] 
                       - j12*j31 * f[1] + j32*j11 * f[1]
                       - j11*j22 * f[2] + j21*j12 * f[2];
         delta_xi[0] /= jac_det;
         delta_xi[1] /= jac_det;
         delta_xi[2] /= jac_det;
         l = 0;
         //perform damping
         while(l<30 && f_test >= f_old){
            l++;
            Double dampFac = 1.0/pow(2.0,(UInt)(l-1.0)); 
            xi_start = xi_k + (delta_xi * dampFac); 
            Local2GlobalCoord(f, xi_start, coordMat, NULL);
            f = f - globalPoint;
            f_test = f.NormL2();
         }
         f_old = f_test;
         xi_k = xi_start;
         k++;
      }
      if( f_test > tolerance){
        std::cout << "performed " << k << " iterations to reach the point" << std::endl<< xi_k << std::endl;
        //ONLY for DEBUGGING
        Local2GlobalCoord(f, xi_k, coordMat, NULL);
        std::cout << "Calculated global point :" << std::endl << f << std::endl;
        std::cout << "Original global coord " << std::endl << globalPoint << std::endl;
        std::cout << "The error was: " << f_test << " and l: " << l << std::endl<< std::endl;
        std::cout << "The last search direction: " << delta_xi << std::endl<< std::endl;
      }

      // Put local coordinate of point into matrix.
      for(UInt j = 0; j < 3; j++) {
        if(xi_k.GetSize() == 0){
          std::cerr << "local2global messed up setting everything to zero" << std::endl;
          std::cerr << globalCoords << std::endl;
          xi_k.Resize(2);
          xi_k.Init();
        }
        localCoords[j][i] = xi_k[j];
      }
    }
  }
  

} // end of namespace

