// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

//#include <General/environment.hh>
#include "hexa1FE.hh"
#include "Domain/elem.hh"
#include "DataInOut/Logging/cfslog.hh"

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


  void Hexa1FE::CalcShapeFnc(Vector<Double> & Shape, 
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
      *error << "Approximation type not known";
      Error( __FILE__, __LINE__ );
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
      LOG_DBG2(hexa1fe) << "LDeriv = \n" << LDeriv << std::endl;
      LOG_DBG3(hexa1fe) << "offset after bubbles: " << offset;
    } else {
        *error << "Approximation type not known";
        Error( __FILE__, __LINE__ );
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
      Error("CalcLocalICModesDerivShapeFnc for Legendre type not implemented",
	    __FILE__,__LINE__);
    }

  }


  void  Hexa1FE::GetNumFncs( Vector<UInt>& numFcns,  
                             const shared_ptr<AnsatzFct>& fcnType, 
                             AnsatzFct::FctEntityType fctEntityType, 
                             UInt dof ) {

    // Check ansatzFctType
    if( fcnType->GetType() == AnsatzFct::LAGRANGE ) {
      numFcns.Resize( NumNodes_ );
      numFcns.Init(1);
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
          Error( "Not yet implemented!", __FILE__, __LINE__ );
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
          Error( "Not yet implemented!", __FILE__, __LINE__ );
        }
      }
    } else {
      *error << "AnsatzFcnType '" << fcnType->GetType() 
             << "' is not known!";
      Error( __FILE__, __LINE__ );
    } 
    
  }

  UInt Hexa1FE::GetNumFncs( const shared_ptr<AnsatzFct>& type ) {
    
    // Check ansatzFctType
    if( type->GetType() == AnsatzFct::LAGRANGE ) {
      return NumNodes_;
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
    
    if( actFct->GetType() == AnsatzFct::LEGENDRE 
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
        std::cerr << "integration order: " << IntegOrder << std::endl;
    //   }

      // get the values by IntegMethod and IntegOrder
      SetIntPoints();
      
      // ... then calc shape function values at integration points
      // for subsequent calls to CalcJacobian() ....
      SetShapeFncAtIp();
      SetShapeFncDerivAtIp();
    }

  }
  

} // end of namespace

