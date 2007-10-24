// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "quad1fe.hh"
#include "Domain/elem.hh"

namespace CoupledField
{


  Quad1FE :: Quad1FE() : RectangleFE()
  {

    UseICModes();
    Init();
  }
  
  Quad1FE :: ~Quad1FE()
  {
  }

  void Quad1FE :: Init()
  {
    NumNodes_ = 4;
    
    CommonInit();
    SetEdgeIndices();
    SetFaceIndices();
  }

  void Quad1FE :: SetCornerCoords()
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

  void Quad1FE :: SetFaceIndices() {

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

  void Quad1FE :: SetEdgeIndices() {
    
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

  void Quad1FE :: CalcShapeFnc( Vector<Double> & Shape, 
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
 

      // ----------------------
      //  c) bubble functions
      // ----------------------
      Double val_i, val_j;
      for( Integer i = 2; i <= order-2; i++ ) {
        for( Integer j = 2; j <= order-2; j++ ) {
          
          // Condition for trunk space: i+j in [2..max{order_i,order_j}]
          if ( (i+j) > order ) {continue;}
          EvalPolynom( val_i, deriv, i, lCoeff_[i], actCoord[0] );
          EvalPolynom( val_j, deriv, j, lCoeff_[j], actCoord[1] );

          Shape[offset++] = val_i * val_j;
        }
      }
    }
  }


  void Quad1FE ::CalcLocalDerivShapeFnc( Matrix<Double> & LDeriv, 
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

      // ---------------------
      //  c) bubble functions
      // ---------------------
      Double val_i, deriv_i, val_j, deriv_j;
      for( Integer i = 2; i <= order-2; i++ ) {
        for( Integer j = 2; j <= order-2; j++ ) {
    
          // Condition for trunk space: i+j in [2..max{order_i,order_j}]
          if ( (i+j) > order ) {continue;}
          EvalPolynom( val_i, deriv_i, i, lCoeff_[i], actCoord[0] );
          EvalPolynom( val_j, deriv_j, j, lCoeff_[j], actCoord[1] );

          LDeriv[offset][0] = deriv_i * val_j;
          LDeriv[offset][1] =   val_i * deriv_j;
          offset++;
        }
      }
    }
  }


  void Quad1FE ::CalcLocalICModesDerivShapeFnc( Matrix<Double> & LDeriv, 
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
        
        UInt inc = 1;
        UInt total = 0;
        for (Integer i = 4; i<=order; i++ ) {
          total += inc++;
        }
        numFcns.Init(total);
      } else if( fctEntityType == AnsatzFct::FACE ) {
        numFcns.Clear();
      } else {
        Error( "Not yet implemented!", __FILE__, __LINE__ );
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
      Integer order =  dynamic_pointer_cast<LegendreFct, AnsatzFct>
        (type)->GetIsoOrder();
      UInt numBubbles = 0;
      UInt inc = 1;
      for (Integer i = 4; i<=order; i++ ) {
        numBubbles += inc++;
      }
      return (NumNodes_ + (NumEdges_*(order-1)) + numBubbles);
      
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
        IntegOrder =  EncodeCartesianOrder( legFct->GetMaxOrder() + 2,
                                            legFct->GetMaxOrder() + 2,
                                            legFct->GetMaxOrder() + 2);
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


