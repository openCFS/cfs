/* =====================================================================
 *
 *    Filename: H1LagrangeVar.cc 
 *
 *    Description: This file contains the implementation of the Lagrange
 *                 polinomial based tensor product elements of arbitrary order.
 *                 These are the spectral elements
 *                 Available right now:
 *                  - Line
 *                  - Quad
 *                  - Hex
 *
 *    Version:  0.1
 *    Created:  20.12.2009 11:42:58 CEST
 *    Revision:  none
 * 
 *    Author:  Andreas Hueppe
 *    Company: Universitaet Klagenfurt
 * ==========================================================================
 */

#include "H1ElemsLagVar.hh"
// header for logging
#include "DataInOut/Logging/cfslog.hh"

#include <algorithm>

namespace CoupledField {

  // declare class specific logging stream
  DECLARE_LOG(H1LagrangeVar)
  DEFINE_LOG(H1LagrangeVar,"H1LagrangeVar")

  //========================================================================
  //Base Class for Spectral elements
  //========================================================================
  
    FeH1LagrangeVar::FeH1LagrangeVar(){
    order_ = 0;
    // Precalculate the supporting Points up do order 5
      CalcAllSupportingPoints(5);
    }

    FeH1LagrangeVar::~FeH1LagrangeVar(){
    }

    void FeH1LagrangeVar::SetIntPoints( StdVector<LocPoint>& intPoints ){
    }

    void FeH1LagrangeVar::GetNumFncs( StdVector<UInt>& numFcns,
                     EntityType fctEntityType,
                     UInt dof){
      if( fctEntityType == VERTEX ) {
        numFcns.Resize( shape_.numVertices );
        numFcns.Init( dof );
      } else if ( fctEntityType == EDGE ) {
        numFcns.Resize( shape_.numEdges*(order_ - 1) );
        numFcns.Init( dof );
      } else if ( fctEntityType == FACE ) {
        numFcns.Resize( shape_.numFaces* ( (order_ - 1)*(order_ - 1)) );
        numFcns.Init( dof );
      } else if ( fctEntityType == INTERIOR ) {
        numFcns.Resize(  (order_ - 1)*(order_ - 1)*(order_ - 1)  );
        numFcns.Init( dof );
      } else {
        EXCEPTION( "Entitytype '" << fctEntityType << "' not known!");
      }
    }

    //void FeH1LagrangeVar::GetNumFncs( StdVector<UInt>& numFcns,
    //                 const shared_ptr<AnsatzFct>& fcnType,
    //                 AnsatzFct::FctEntityType fctEntityType,
    //                 UInt dof){
    //  GetNumFncs( numFcns, fctEntityType, dof);
    //}

    void FeH1LagrangeVar::GetNodalPermutation( StdVector<UInt>& fncPermutation,
                                               const Elem* ptElem,
                                               EntityType fctEntityType,
                                               UInt entNumber){
      if( fctEntityType == VERTEX ) {
        fncPermutation.Resize(1);
        fncPermutation.Init(0);
      }else if( fctEntityType == EDGE ) {
        fncPermutation.Resize(order_-1);

        //TODO> safety check if the requested entNumber is valid
        Integer factor = ptElem->edges[entNumber]; 
        if(factor < 0 ){
          for ( UInt i = 0; i < order_-1 ; i++ ) {
            fncPermutation[i] = order_-i-2;
          }
        }else{
          for ( UInt i = 0; i < order_-1 ; i++ ) {
            fncPermutation[i] = i;
          }
        }

        LOG_DBG3(H1LagrangeVar) << "Edge # " << ptElem->edges[entNumber] << " of Element #" 
                                   << ptElem->elemNum << "Has sign " << factor << " and the permutation \n"
                                   << fncPermutation << std::endl;

      }else if( fctEntityType == FACE && ptElem->faces.GetSize() > 0) {
        fncPermutation.Resize((order_-1) * (order_-1));
        Integer dI,dII;
        if(ptElem->faceFlags[entNumber].test(0)){
          //richtungI = flag(2);
          //richtungII = flag(1);
          dI = (ptElem->faceFlags[entNumber].test(2))? 0:order_-2;
          dII = (ptElem->faceFlags[entNumber].test(1))? 0:order_-2;
          for(Integer i = 0; i< order_-1 ; i++){
            for(Integer j = 0; j< order_-1 ; j++){
              fncPermutation[(i*(order_-1)) + j] = (abs(dI-i)*(order_-1)) + abs(dII-j);
            }
          }
        }else{
          //richtungI = flag(1);
          //richtungII = flag(2);
          dI = (ptElem->faceFlags[entNumber].test(1))? 0:order_-2;
          dII = (ptElem->faceFlags[entNumber].test(2))? 0:order_-2;
          for(Integer i = 0; i< order_-1 ; i++){
            for(Integer j = 0; j< order_-1 ; j++){
              fncPermutation[(i*(order_-1)) + j] = (abs(dI - j)*(order_-1)) + abs(dII - i);
            }
          }
        }
        LOG_DBG3(H1LagrangeVar) << "Face # " << ptElem->faces[entNumber] << " of Element #" 
                                << ptElem->elemNum << " Has Bitset " << ptElem->faceFlags[entNumber] << " and the permutation \n"
                                << fncPermutation << std::endl;

        //the following check for interior nodes is also limited to hexas where there are 4 faces
      }else if( fctEntityType == INTERIOR && ptElem->faces.GetSize() == 4) {
        
        //no need to check for an orientation
        UInt numIFncs = (order_-1)*(order_-1)*(order_-1);
        fncPermutation.Resize(numIFncs);

        for(UInt k = 0; k< numIFncs ; k++){
          fncPermutation[k] = k;
        }
      } else{
        fncPermutation.Resize(0);
      }
    }

    //void FeH1LagrangeVar::GetShFnc( Vector<Double> & S, const LocPoint& lp,
    //               const Elem* ptElem,  UInt comp ){
    //  CalcShFnc( S, lp.coord);
    //}

    //void FeH1LagrangeVar::GetDerivShFnc( Matrix<Double> & deriv, 
    //                    const LocPoint& lp,
    //                    const Elem * elem, 
    //                    UInt comp ){
    //  CalcDerivShFnc( deriv, lp.coord);
    //}


  //=========================================================================
  //Line Elements of arbitrary order
  //=========================================================================
  
    FeH1LagrangeLineVar::FeH1LagrangeLineVar(){
      feType_ = Elem::ET_LINE2;
      shape_ = Elem::shapes[feType_];
      //this is always element order +1
      actNumFncs_ = 2;
      order_ = 1;

    }

    FeH1LagrangeLineVar::~FeH1LagrangeLineVar(){
    }

    void FeH1LagrangeLineVar::CalcShFnc( Vector<Double>& Shape,
                                         const Vector<Double>& point,
                                         const Elem* ptElem,
                                         UInt comp  ) {
      Shape.Resize( actNumFncs_ );
      Shape.Init();
      //now get the shape functions and the derivatives for the given coordinates
      if(supportingPoints_.find(actNumFncs_-1) == supportingPoints_.end() )
       supportingPoints_[actNumFncs_-1] = CalcGaussLobattoPoints(actNumFncs_-1);

      Vector<Double> lagPoly;
      EvaluateLagrangePolynomial(lagPoly, point[0],order_ );

      UInt counter = 0;
      Shape[counter++] = lagPoly[0];
      Shape[counter++] = lagPoly[order_];

      for ( UInt i= 1 ; i< order_  ;i++ ){
        Shape[counter++] = lagPoly[i];
      }
      return;
    }

    void FeH1LagrangeLineVar::CalcLocDerivShFnc( Matrix<Double> & deriv, 
                                                 const Vector<Double>& point,
                                                 const Elem* ptElem,
                                                 UInt comp ) {
      if(supportingPoints_.find(actNumFncs_-1) == supportingPoints_.end() )
       supportingPoints_[actNumFncs_-1] = CalcGaussLobattoPoints(actNumFncs_-1);

      Vector<Double> lagPoly;

      deriv.Resize( actNumFncs_, shape_.dim );

      //now get the shape functions and the derivatives for the given coordinates
      EvaluateDerivLagrangePolynomial( lagPoly,point[0],order_ );

      deriv[0][0] = lagPoly[0];
      deriv[1][0] = lagPoly[order_];

      UInt counter = 2;
      for ( UInt i=1;i<order_;i++ )
        deriv[counter++][0] = lagPoly[order_ - i];
      return;
    }
 
  void  FeH1LagrangeLineVar::SetIsoOrder(UInt order){
    order_ = order;
    actNumFncs_ = order+1;
  }
  //=========================================================================
  //Quadrilateral elements of arbitrary order
  //=========================================================================

    FeH1LagrangeQuadVar::FeH1LagrangeQuadVar(){
      feType_ = Elem::ET_QUAD4;
      shape_ = Elem::shapes[feType_];
      actNumFncs_ = 4;
      order_ = 1;
    }

    FeH1LagrangeQuadVar::~FeH1LagrangeQuadVar(){
    }

    void FeH1LagrangeQuadVar::CalcShFnc( Vector<Double>& Shape,
                                         const Vector<Double>& point,
                                         const Elem* ptElem,
                                         UInt comp ) {
      Shape.Resize( actNumFncs_ * actNumFncs_ );
      Shape.Init();
      //now get the shape functions and the derivatives for the given coordinates
      Vector<Double> shapeX;
      Vector<Double> shapeY;

      shapeX.Resize(actNumFncs_);
      shapeY.Resize(actNumFncs_);
      shapeX.Init();
      shapeY.Init();
      EvaluateLagrangePolynomial(shapeX, point[0], order_ );
      EvaluateLagrangePolynomial( shapeY, point[1], order_ );

      //fill vertices
      Shape[0]= shapeX[0]      * shapeY[0];
      Shape[1]= shapeX[order_] * shapeY[0];
      Shape[2]= shapeX[order_] * shapeY[order_];
      Shape[3]= shapeX[0]      * shapeY[order_];

      UInt offset = shape_.numVertices;

      //fill edge 1
      for ( UInt i= 1 ; i< order_ ;i++ ){
        Shape[offset++] = shapeX[i] * shapeY[0];
      }
      //fill edge 2
      for ( UInt i= 1 ; i< order_  ;i++ ){
        Shape[offset++] = shapeX[order_] * shapeY[i];
      }
      //fill edge 3
      for ( UInt i= 1 ; i< order_  ;i++ ){
        Shape[offset++] = shapeX[i] * shapeY[order_];
      }
      //fill edge 4
      for ( UInt i= 1 ; i< order_  ;i++ ){
        Shape[offset++] = shapeX[0] * shapeY[i];
      }

      //fill the face
      for ( UInt j= 1;j< order_ ;j++ ){                                            
        for ( UInt i = 1; i< order_ ;i++ ){                                          
          Shape[offset++] = shapeX[i] * shapeY[j]; 
        }                                          
      }                                            
      //std::cout << Shape << std::endl;
      //std::cout << shapeX << std::endl;
      //std::cout << shapeY << std::endl;
      return;
    }

    void FeH1LagrangeQuadVar::CalcLocDerivShFnc( Matrix<Double> & deriv, 
                                                 const Vector<Double>& point,
                                                 const Elem* ptElem,
                                                 UInt comp ) {
      Vector<Double> shapeX;
      Vector<Double> shapeY;
      Vector<Double> shapeDerivX;
      Vector<Double> shapeDerivY;

      deriv.Resize( (order_+1) * (order_+1), shape_.dim );
      //now get the shape functions and the derivatives for the given coordinates
      EvaluateLagrangePolynomial( shapeX, point[0], order_ );
      EvaluateLagrangePolynomial( shapeY, point[1], order_ );
      EvaluateDerivLagrangePolynomial( shapeDerivX, point[0], order_ );
      EvaluateDerivLagrangePolynomial( shapeDerivY, point[1], order_ );

      deriv[0][0]= shapeDerivX[0]      * shapeY[0];
      deriv[0][1]= shapeX[0]           * shapeDerivY[0];
      deriv[1][0]= shapeDerivX[order_] * shapeY[0];
      deriv[1][1]= shapeX[order_]      * shapeDerivY[0];
      deriv[2][0]= shapeDerivX[order_] * shapeY[order_];
      deriv[2][1]= shapeX[order_]      * shapeDerivY[order_];
      deriv[3][0]= shapeDerivX[0]      * shapeY[order_];
      deriv[3][1]= shapeX[0]           * shapeDerivY[order_];

      UInt offset = shape_.numVertices;

      //now for the edge derivatives
      for ( UInt run= 1; run< order_ ;run++ ){              
        deriv[offset][0] = shapeDerivX[run] * shapeY[0];           
        deriv[offset++][1] = shapeX[run] * shapeDerivY[0];         
      }                                                    
      for ( UInt run= 1; run< order_ ;run++ ){              
        deriv[offset][0] = shapeDerivX[order_] * shapeY[run];           
        deriv[offset++][1] = shapeX[order_] * shapeDerivY[run];         
      }                                                    
      for ( UInt run= 1; run< order_ ;run++ ){              
        deriv[offset][0] = shapeDerivX[run] * shapeY[order_];           
        deriv[offset++][1] = shapeX[run] * shapeDerivY[order_];         
      }                                                    
      for ( UInt run= 1; run< order_ ;run++ ){              
        deriv[offset][0] = shapeDerivX[0] * shapeY[run];           
        deriv[offset++][1] = shapeX[0] * shapeDerivY[run];         
      }                                                    

      for ( UInt i= 1;i< order_;i++ ){
        for ( UInt j = 1; j< order_;j++ ){
          deriv[offset][0] = shapeDerivX[i] * shapeY[j];
          deriv[offset++][1] = shapeX[i] * shapeDerivY[j];
        }
      }

    }

  void  FeH1LagrangeQuadVar::SetIsoOrder(UInt order){
    order_ = order;
    actNumFncs_ = (order+1)*(order+1);
  }

  //=========================================================================
  //Hexahedral elements of arbitrary order
  //=========================================================================
  
    FeH1LagrangeHexVar::FeH1LagrangeHexVar(){
      feType_ = Elem::ET_HEXA8;
      shape_ = Elem::shapes[feType_];
    }

    FeH1LagrangeHexVar::~FeH1LagrangeHexVar(){
    }

    UInt FeH1LagrangeVar::GetNumFncsPerEntType( EntityType fctEntityType,
                                       UInt dof){
      UInt numFncs = 0;
      if( fctEntityType == VERTEX ) {
        numFncs = shape_.numVertices;
      } else if ( fctEntityType == EDGE) {
        numFncs =  order_ - 1; 
      } else if ( fctEntityType == FACE) {
        numFncs = (order_ - 1)*(order_ - 1);
      } else if ( fctEntityType == INTERIOR ) {
        numFncs = (order_ - 1)*(order_ - 1)*(order_ - 1);
      }
      return numFncs;
    }
    void FeH1LagrangeHexVar::CalcShFnc( Vector<Double>& Shape,
                                        const Vector<Double>& point,
                                        const Elem* ptElem,
                                        UInt comp ) {
    }

    void FeH1LagrangeHexVar::CalcLocDerivShFnc( Matrix<Double> & deriv, 
                                                const Vector<Double>& point,
                                                const Elem* ptElem,
                                                UInt comp ) {
    }

    void FeH1LagrangeHexVar::SetIsoOrder(UInt order){
      order_ = order;
      actNumFncs_ = (order+1)*(order+1)*(order+1);
    }
}
