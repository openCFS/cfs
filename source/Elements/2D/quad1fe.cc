// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "quad1fe.hh"
#include "Domain/elem.hh"
#include "DataInOut/Logging/cfslog.hh"

namespace CoupledField
{

  // declare class specific logging stream
  DECLARE_LOG(quad1fe)
  DEFINE_LOG(quad1fe,"quad1fe")

  Quad1FE::Quad1FE() : RectangleFE()
  {

    UseICModes();
    Init();
  }
  
  Quad1FE::~Quad1FE()
  {
  }

  void Quad1FE::Init()
  {
    NumNodes_ = 4;
    
    CommonInit();
    SetEdgeIndices();
    SetFaceIndices();
  }

  void Quad1FE::SetCornerCoords()
  {

    LCornerCoords_.Resize(Dim_,NumNodes_);
  
    LCornerCoords_[0][0] =  -1;
    LCornerCoords_[1][0] =  -1;
    LCornerCoords_[0][1] = 1;
    LCornerCoords_[1][1] =  -1;
    LCornerCoords_[0][2] = 1;
    LCornerCoords_[1][2] = 1;
    LCornerCoords_[0][3] = -1;
    LCornerCoords_[1][3] = 1;


  }

  void Quad1FE::SetFaceIndices() {

    faceIndices_ = new StdVector<UInt>[NumFaces_];
    for (UInt i = 0; i < NumFaces_; i++) {
      faceIndices_[i].Resize(4);
    }
   
    // face 1
    faceIndices_[0][0] = 1;
    faceIndices_[0][1] = 2;
    faceIndices_[0][2] = 3;
    faceIndices_[0][3] = 4;
  }

  void Quad1FE::SetEdgeIndices() {
    
    edgeIndices_ = new StdVector<UInt>[NumEdges_];
    for (UInt i=0; i<NumEdges_; i++) {
      edgeIndices_[i].Resize(2);
    }

    // Note: The orientation is taken from
    // A. Duester: High Order FEM, Lecture Nodes,
    // p. 25

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
        
  }

  void Quad1FE::CalcShapeFnc( Vector<Double> & Shape, 
                                const Vector<Double> & actCoord,
                                const Elem* elem,
                                UInt dof, AnsatzFct::FctEntityType type )
  {


    // Check ansatzFctType
    if(  actFct_->GetType() == AnsatzFct::LAGRANGE ||
         type == AnsatzFct::NODE ) {

      // ===============
      //  LAGRANGE PART
      // ===============
      Shape.Resize(NumNodes_);
      
      for( UInt i=0; i<NumNodes_; i++)
        Shape[i] = 0.25 * (1 + LCornerCoords_[0][i] * actCoord[0])
          * (1 + LCornerCoords_[1][i] * actCoord[1]);
    } else {
      // Get number of ansatz functions
      shared_ptr<LegendreFct> legFct = 
        dynamic_pointer_cast<LegendreFct, AnsatzFct>(actFct_);

      // ===============
      //  LEGENDRE PART
      // ===============

      // Get number of ansatz functions
      UInt totalFcns = GetNumFncs( actFct_ );

      Shape.Resize( totalFcns );

      // --------------------
      //  a) nodal functions
      // --------------------
         // Offset for different functions
      UInt offset = 0;
      
      // First of all, calculate all nodal function derivatives
      for( UInt iNode = 0; iNode < NumCorners_; iNode++,offset++ ) {
        Shape[offset] = 0.25 * (1 + LCornerCoords_[0][iNode] * actCoord[0] )
                             * (1 + LCornerCoords_[1][iNode] * actCoord[1] );
      }

      // --------------------
      //  b) edge functions
      // --------------------
      // Obtain order of element
      Integer order = 
        dynamic_pointer_cast<LegendreFct, AnsatzFct>(actFct_)->GetIsoOrder();
      Double val, factor, deriv;      

#define QUAD_EDGE_FCN(edgeNum,  sign_1, dir_1, dir_2 )                  \
      factor = elem->edges[edgeNum-1] < 0 ? -1.0 : 1.0;                 \
      for( Integer iDof = 2; iDof <= order; iDof++, offset++ ) {        \
        EvalPolynom( val, deriv, iDof, lCoeff_[iDof],                   \
                     factor*actCoord[dir_2] );                          \
        Shape[offset] = 0.5 * ( 1 sign_1 actCoord[dir_1] ) * val;       \
      }

      //  EDGE #1
      QUAD_EDGE_FCN( 1, -, 1, 0 );

      //  EDGE #2
      QUAD_EDGE_FCN( 2, +, 0, 1 );

      //  EDGE #3
      QUAD_EDGE_FCN( 3, +, 1, 0 );

      //  EDGE #4
      QUAD_EDGE_FCN( 4, -, 0, 1 );
 
      // --------------------
      //  c) face functions
      // --------------------

      Double val_1, deriv_1, val_2, deriv_2, val_3, deriv_3;
      Double sFlag2, sFlag3;
      Integer order_1, order_2, order_3;
      UInt d2, d3;
  
      order_2 = legFct->GetMaxOrderLocDir(0);                        
      order_3 = legFct->GetMaxOrderLocDir(1);                        
      if( elem->faceFlags[0].test(0) == true) {              
        sFlag2 = (elem->faceFlags[0].                        
                  test(2) == true) ? 1.0 : -1.0;                     
        sFlag3 = (elem->faceFlags[0].                        
                  test(1) == true) ? 1.0 : -1.0;                     
        d2 = 0;                                                   
        d3 = 1;                                                   
      } else {                                                       
        sFlag3 = (elem->faceFlags[0].test(2) == true) ? 1.0 : -1.0; 
        sFlag2 = (elem->faceFlags[0].test(1) == true) ? 1.0 : -1.0; 
        d3 = 0;                                                      
        d2 = 1;                                                      
      }                                                                 
      for( Integer i = 2; i <= (order_2)-2; i++ ) {                     
        for( Integer j = 2; j <= (order_3)-2; j++ ) {                   
          if( (i + j) > std::max(order_2, order_3) ) { 
            continue; 
          }      
          EvalPolynom( val_2, deriv_2, i, lCoeff_[i],                   
                       sFlag2* actCoord[d2] );                          
          EvalPolynom( val_3, deriv_3, j, lCoeff_[j],                   
                       sFlag3*actCoord[d3] );                           
          Shape[offset] = val_2  * val_3;                         
          offset++;                                                     
        } 
      }
    }
  }


  void Quad1FE::CalcLocalDerivShapeFnc( Matrix<Double> & LDeriv, 
                                         const Vector<Double> & actCoord,
                                         const Elem* elem,
                                         UInt dof,
                                         AnsatzFct::FctEntityType type )
  {



    if( actFct_->GetType() == AnsatzFct::LAGRANGE ||
        type == AnsatzFct::NODE ) {

      // ===============
      //  LAGRANGE PART
      // ===============
      LDeriv.Resize(NumNodes_,Dim_);
      for( UInt i=0; i<NumNodes_; i++)
	{
	  LDeriv[i][0] = 0.25 * LCornerCoords_[0][i] 
	    * (1 + LCornerCoords_[1][i] * actCoord[1] );
	  LDeriv[i][1] = 0.25 * (1 + LCornerCoords_[0][i] * actCoord[0] )
	    * LCornerCoords_[1][i];
	}
    
    } else {

      // ===============
      //  LEGENDRE PART
      // ===============
      shared_ptr<LegendreFct> legFct = 
        dynamic_pointer_cast<LegendreFct, AnsatzFct>(actFct_);

      // Get number of ansatz functions
      UInt totalFcns = GetNumFncs( actFct_ );
      LDeriv.Resize( totalFcns, Dim_ );

      // --------------------
      //  a) nodal functions
      // --------------------
      UInt offset = 0;
      
      // First of all, calculate all nodal function derivatives
      for( UInt iNode = 0; iNode < NumCorners_; iNode++,offset++ ) {
        
        LDeriv[offset][0] = 0.25 * LCornerCoords_[0][iNode] 
          * (1 + LCornerCoords_[1][iNode] * actCoord[1] );
        
        LDeriv[offset][1] = 0.25 * (1 + LCornerCoords_[0][iNode] * actCoord[0] )
          * LCornerCoords_[1][iNode];
      }

      // -------------------
      //  b) edge functions
      // -------------------
      
      // Obtain order of element
      Integer order_1, order_2, order_3;
      Double val, deriv, factor;
      Integer order = 
        dynamic_pointer_cast<LegendreFct, AnsatzFct>(actFct_)->GetIsoOrder();
      
#define QUAD_EDGE_DERIV( edgeNum, sign_1, dir_1, dir_2 )                \
      factor = elem->edges[edgeNum-1] < 0 ? -1.0 : 1.0;                 \
      for( Integer iDof = 2; iDof <= order; iDof++, offset++ ) {        \
        EvalPolynom( val, deriv, iDof, lCoeff_[iDof],                   \
                     factor*actCoord[dir_2] );                          \
        LDeriv[offset][dir_1] = sign_1(0.5 * val);                      \
        LDeriv[offset][dir_2] =  0.5 * ( 1.0 sign_1 actCoord[dir_1])    \
                                     * deriv * factor;                  \
      }

      // EDGE #1
      QUAD_EDGE_DERIV( 1, -, 1, 0 );
      
      // EDGE #2
      QUAD_EDGE_DERIV( 2, +, 0, 1 );

      // EDGE #3
      QUAD_EDGE_DERIV( 3, +, 1, 0 );

      // EDGE #4
      QUAD_EDGE_DERIV( 4, -, 0, 1 );

      // -------------------
      //  c) face functions
      // -------------------
       
       Double val_2, deriv_2, val_3, deriv_3;
      Double sFlag2, sFlag3;
      UInt d2, d3;

      order_2 = legFct->GetMaxOrderLocDir(0);                       
      order_3 = legFct->GetMaxOrderLocDir(1);                       
      LOG_DBG3(quad1fe) << "Treating face " << elem->faces[0];  
      LOG_DBG3(quad1fe) << "bitset: " << elem->faceFlags[0];    
      if( elem->faceFlags[0].test(0) == true) {                 
        sFlag2 = (elem->faceFlags[0].                           
                  test(2) == true) ? 1.0 : -1.0;                        
        sFlag3 = (elem->faceFlags[0].                           
                  test(1) == true) ? 1.0 : -1.0;                        
        d2 = 0;                                                     
        d3 = 1;                                                     
      } else {                                                          
        sFlag3 = (elem->faceFlags[0].test(2) == true) ? 1.0 : -1.0; 
        sFlag2 = (elem->faceFlags[0].test(1) == true) ? 1.0 : -1.0; 
        d3 = 0;                                                     
        d2 = 1;                                                     
                                                                        
      }                                                                 
      LOG_DBG3(quad1fe) << "I-Flag: " << sFlag2                         
                        << ", II-Flag: " << sFlag3;                     
      LOG_DBG3(quad1fe) << "I-direction: " << d2                        
                        << ", II-direction: " << d3;                    
      for( Integer i = 2; i <= (order_2) - 2; i++ ) {                   
        for( Integer j = 2; j <= (order_3) - 2; j++ ) {                 
          if( (i + j) > std::max(order_2,order_3) ) { continue; }       
          EvalPolynom( val_2, deriv_2, i, lCoeff_[i],                   
                       sFlag2* actCoord[d2] );                          
          EvalPolynom( val_3, deriv_3, j, lCoeff_[j],                   
                       sFlag3*actCoord[d3] );                           
          LOG_DBG3(quad1fe) << "deriv_2 = " << deriv_2                  
                            << ", deriv_3 = " << deriv_3;               
          LOG_DBG3(quad1fe) << "val_2 = " << val_2                      
                            << ", val_3 = " << val_3 << std::endl;      
          LDeriv[offset][d2] = deriv_2 * sFlag2 * val_3;
          LDeriv[offset][d3] = val_2 * deriv_3 * sFlag3;
          offset++;                                                     
        } 
      }
    }
  }


  void Quad1FE::CalcLocalICModesDerivShapeFnc( Matrix<Double> & LDeriv, 
						const Vector<Double> & actCoord,
						const Elem* elem,
						UInt dof, AnsatzFct::FctEntityType type )
  {



    if( actFct_->GetType() == AnsatzFct::LAGRANGE ||
        type == AnsatzFct::NODE ) {

      // ===============
      //  LAGRANGE PART
      // ===============
      LDeriv.Resize(2,Dim_);
      LDeriv.Init();
      LDeriv[0][0] = -2.0*actCoord[0];
      LDeriv[1][1] = -2.0*actCoord[1];
    } 
    else {
      Error("CalcLocalICModesDerivShapeFnc for Legendre type not implemented",
	    __FILE__,__LINE__);
    }
  }


  Double Quad1FE::CalcMeanStrain(Matrix<Double> &cornerCoords, 
                                 Matrix<Double> &displacements)
  {

    Double factor;
    Double eps1, eps2, eps4, eps5,
           eps11, eps12, eps21, eps22, eps41, eps42, eps51, eps52;
    Double length1, length2, length11, length12, length21, length22;


    // calculate inital size of element
    length11 = abs(cornerCoords[0][0] - cornerCoords[0][1]);
    length12 = abs(cornerCoords[0][3] - cornerCoords[0][2]);
    length1 = (length11+length12) * 0.5;
  
    length21 = abs(cornerCoords[1][0] - cornerCoords[1][3]);
    length22 = abs(cornerCoords[1][1] - cornerCoords[1][2]);
    length2 = (length21+length22) * 0.5;

    // calculate absolute change of size
    eps11 = displacements[0][0] - displacements[0][1];
    eps12 = displacements[0][3] - displacements[0][2];
    eps1 = (eps11+eps12) * 0.5;

    eps21 = displacements[1][0] - displacements[1][3];
    eps22 = displacements[1][1] - displacements[1][2];
    eps2 = (eps21+eps22) * 0.5;
  
    eps41 = displacements[0][2] - displacements[0][1];
    eps42 = displacements[0][3] - displacements[0][0]; 
    eps4 = (eps41+eps42)*0.5;
    
    eps51 = displacements[1][1] - displacements[1][0];
    eps52 = displacements[1][3] - displacements[1][2];
    eps5= (eps51+eps52)*0.5;  

    factor =  0.2 * ((eps1*eps1/(length1*length1))
                     + (eps2*eps2/(length2*length2)) 
                     + (eps5*eps5/(length1*length2))
                     + (eps4*eps4/(length1*length1)));

    return factor;
  }

  void  Quad1FE::GetNumFncs( Vector<UInt>& numFcns,  
                             const shared_ptr<AnsatzFct>& fcnType, 
                             AnsatzFct::FctEntityType fctEntityType, 
                             UInt dof ) {

    // Check ansatzFctType
    if( fcnType->GetType() == AnsatzFct::LAGRANGE ) {
      numFcns.Resize( NumNodes_ );
      numFcns.Init(1);
    } else if ( fcnType->GetType() == AnsatzFct::LEGENDRE ) {
      
      // Remember approximation order
      Integer order =  dynamic_pointer_cast<LegendreFct, AnsatzFct>
        (fcnType)->GetIsoOrder();
      if( fcnType->IsIsotropic() == true ) {
        // Check for subentity-type
        if( fctEntityType == AnsatzFct::NODE ) {
          numFcns.Resize( NumNodes_ );
          numFcns.Init( 1 );
          
        } else if( fctEntityType == AnsatzFct::EDGE ) {
          numFcns.Resize( NumEdges_ );
          numFcns.Init( order - 1 );
        } 
        else if( fctEntityType == AnsatzFct::INTERIOR ) {
          numFcns.Resize(1);
          numFcns.Init(0);
        } else if( fctEntityType == AnsatzFct::FACE ) {
          numFcns.Resize(1);
          
          UInt inc = 1;
          UInt total = 0;
          for (Integer i = 4; i<=order; i++ ) {
            total += inc++;
          }
          numFcns.Init(total);
        } else {
          Error( "Not yet implemented!", __FILE__, __LINE__ );
        }
      }else {
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
         
          // eta-direction
          numFcns[1]  = order[1][dof]-1;
          numFcns[3]  = numFcns[1];
                   
        } else if( fctEntityType == AnsatzFct::FACE ) {
          numFcns.Resize( NumFaces_ );
          
          UInt numFaceFncs = 0;
          Integer max = 0;
           // Face #1
           // max = std::max( order[0][dof], order[1][dof] );
          for( Integer i = 2; i<= (Integer)order[0][dof]-2; i++ ) {
            for( Integer j = 2; j<= (Integer)order[1][dof]-2; j++ ) {
              if( i+j <= max) {
                numFaceFncs++;
              }
            }
          }
          numFcns[0] = numFaceFncs;
        }else if( fctEntityType == AnsatzFct::INTERIOR ) {
          numFcns.Resize(1);
          numFcns.Init(0);
        }

      } 
    } else {
      *error << "AnsatzFcnType '" << fcnType->GetType() 
             << "' is not known!";
      Error( __FILE__, __LINE__ );
    }
  }
  
  UInt Quad1FE::GetNumFncs( const shared_ptr<AnsatzFct>& type ) {
    
    // TODO: FOr anisotropic functions we have
    // to incorporate the dof in the determination
    // Check ansatzFctType
    if( type->GetType() == AnsatzFct::LAGRANGE ) {
      return NumNodes_;
    } else if ( type->GetType() == AnsatzFct::LEGENDRE ) {
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
        
        LOG_DBG3(quad1fe) << "numFaceFncs = " << numFaceFncs << std::endl;
                
        LOG_DBG3(quad1fe) << "total number of unknowns: "
                          <<  NumNodes_ + numFaceFncs + numEdgeFncs 
                          << std::endl;
        
        actNumFcns_ = (NumNodes_ + numFaceFncs + numEdgeFncs);
        return actNumFcns_;
        
      } else {
        // *** anisotropic case ***
        Integer max_1, max_2;
        shared_ptr<LegendreFct> legFct = 
          dynamic_pointer_cast<LegendreFct, AnsatzFct>(type);

        max_1 = legFct->GetMaxOrderLocDir( 0 );
        max_2 = legFct->GetMaxOrderLocDir( 1 );
        
        // a) nodes
        UInt numNodeModes = 4;
        
        // b) edges
        UInt numEdgeModes = 0;
        numEdgeModes += max_1 > 0 ? (max_1-1)*2 : 0;
        numEdgeModes += max_2 > 0 ? (max_2-1)*2 : 0;
               
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
           
        UInt actNumFcns = numNodeModes + numEdgeModes 
          + numFaceModes;
        return actNumFcns;

      }
    }
    return 0;
  }


  void Quad1FE::SetAnsatzFct( shared_ptr<AnsatzFct>& actFct, 
                              bool setIntPoints ) {

    // Check if this ansatz fct was already set
    if( actFct_ == actFct ) {
      return;
    }
    actFct_ = actFct;
    
    if( actFct->GetType() == AnsatzFct::LEGENDRE
        && setIntPoints == true ) {
      
      // If not, get order of functions
      shared_ptr<LegendreFct> legFct = 
        dynamic_pointer_cast<LegendreFct, AnsatzFct>(actFct);
      
      // check if isotropic order is present
      if( legFct->IsIsotropic() ) {
        
        // Prevent integration of first order, as this may
        // cause non-reasonable results
        //if( legFct->GetIsoOrder() > 1 ) {
        IntegMethod = CARTESIAN;
        //IntegOrder = legFct->GetIsoOrder() *2;
        IntegOrder =  EncodeCartesianOrder( legFct->GetIsoOrder() *2,
                                            legFct->GetIsoOrder() *2,
                                            legFct->GetIsoOrder() *2);
          // } else {
          //IntegOrder = 2;
          //} 
      } else {
        IntegMethod = CARTESIAN;
        IntegOrder =  EncodeCartesianOrder( legFct->GetMaxOrder() * 2,
                                            legFct->GetMaxOrder() * 2,
                                            legFct->GetMaxOrder() * 2);
      }
      
      
      // get the values by IntegMethod and IntegOrder
      SetIntPoints();
      
      // ... then calc shape function values at integration points
      // for subsequent calls to CalcJacobian() ....
      SetShapeFncAtIp();
      SetShapeFncDerivAtIp();
      
    }
  }

} // end of namespace


