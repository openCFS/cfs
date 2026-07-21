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
#include "DataInOut/Logging/LogConfigurator.hh"

#include <algorithm>

namespace CoupledField {

  // declare class specific logging stream
  DEFINE_LOG(H1LagrangeVar,"H1LagrangeVar")

  //========================================================================
  //Base Class for Spectral elements
  //========================================================================
  
    FeH1LagrangeVar::FeH1LagrangeVar() 
  : FeH1(), FeNodal() {
    order_ = 0;
    preComputShFnc_ = true;
    
    // All elements here consist of TENSOR_TYPE polynomials
    completeType_ = TENSOR_TYPE;
    
    // Precalculate the supporting Points up do order 10
    CalcAllSupportingPoints(10);
    }

    FeH1LagrangeVar::FeH1LagrangeVar(const FeH1LagrangeVar& other)
                    : FeH1(other), FeNodal(other){
      this->order_ = other.order_;
    }

    FeH1LagrangeVar::~FeH1LagrangeVar(){
    }

    //THIS is a VERSION FOR TENSOR PRODUCT ELEMENTS
    //FOR EVERY OTHER TYPE WE NEED A REIMPLEMENTATION
    void FeH1LagrangeVar::GetNumFncs( StdVector<UInt>& numFcns,
                     EntityType fctEntityType,
                     UInt dof){
      if( fctEntityType == VERTEX ) {
        numFcns.Resize( shape_.numVertices );
        numFcns.Init( 1 );
      } else if ( fctEntityType == EDGE ) {
        numFcns.Resize( shape_.numEdges);
        numFcns.Init( order_-1 );
      } else if ( fctEntityType == FACE && shape_.numFaces > 0) {
        numFcns.Resize( shape_.numFaces);
        numFcns.Init( (order_-1) * (order_-1) );
      } else if ( fctEntityType == INTERIOR && shape_.numFaces == 6 ) {
        numFcns.Resize( 1 );
        numFcns.Init( (order_-1)*(order_-1)*(order_-1) );
      } else if (fctEntityType == ALL ){
        numFcns.Resize(actNumFncs_);
      } else {
        numFcns.Resize(0);
        //EXCEPTION( "Entitytype '" << fctEntityType << "' not known!");
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
        Integer factor = ptElem->extended->edges[entNumber];
        if(factor < 0 ){
          for ( UInt i = 0; i < order_-1 ; i++ ) {
            fncPermutation[i] = order_-i-2;
          }
        }else{
          for ( UInt i = 0; i < order_-1 ; i++ ) {
            fncPermutation[i] = i;
          }
        }

        LOG_DBG3(H1LagrangeVar) << "Edge # " << ptElem->extended->edges[entNumber] << " of Element #"
                                   << ptElem->elemNum << "Has sign " << factor << " and the permutation \n"
                                   << fncPermutation << std::endl;

      }else if( fctEntityType == FACE && ptElem->extended->faces.GetSize() > 0) {
        fncPermutation.Resize((order_-1) * (order_-1));
        Integer xi,eta;
        UInt k=0;

        //check to count xi or eta backward or forward
        Integer offsetXi  =  (ptElem->extended->faceFlags[entNumber].test(0))? 0 : order_-2;
        Integer offsetEta =  (ptElem->extended->faceFlags[entNumber].test(1))? 0 : order_-2;
        if(ptElem->extended->faceFlags[entNumber].test(2)){
          //which means xi and eta are unchanged
          for(eta = 0; eta< (Integer)order_-1 ; eta++){
            for(xi = 0; xi< (Integer)order_-1 ; xi++){
              fncPermutation[k++] = abs(offsetEta-eta)*((order_-1))+abs(offsetXi-xi);
            }
          }
        }else{
          //which means xi and eta are interchanged
          for(xi = 0; xi< (Integer)order_-1 ; xi++){
            for(eta = 0; eta< (Integer)order_-1 ; eta++){
              fncPermutation[k++] = abs(offsetXi-eta)*((order_-1))+abs(offsetEta-xi);
            }
          }
        }
        LOG_DBG3(H1LagrangeVar) << "Face # " << ptElem->extended->faces[entNumber] << " of Element #"
                                << ptElem->elemNum << " Has Bitset " << ptElem->extended->faceFlags[entNumber] << " and the permutation \n"
                                << fncPermutation << std::endl;

        //the following check for interior nodes is also limited to hexas where there are 6 faces
      }else if( fctEntityType == INTERIOR && ptElem->extended->faces.GetSize() == 6) {
        
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
    
    bool FeH1LagrangeVar::operator==( const FeH1LagrangeVar& comp) const {
      bool ret = true;
      ret &= this->feType_ == comp.feType_;
      ret &= this->order_ == comp.order_;
      return ret;
    }

    
    void FeH1LagrangeVar::SetFunctionsAtIp(const StdVector<LocPoint>& iPoints){

      //! Precompute shape functions only, if we are allowed to, i.e.
      //! if we have no higher order shape functions, depending on 
      //! the global orientation of the element
      if( !preComputShFnc_ ) 
        return;
      //shapeFncsAtIp_.resize(iPoints.GetSize());
      //shapeFncDerivsAtIp_.resize(iPoints.GetSize());
      for(UInt aPoint = 0; aPoint < iPoints.GetSize();aPoint++){
        const LocPoint& lp = iPoints[aPoint];
        CalcShFnc( shapeFncsAtIp_[lp.number], lp.coord, NULL, 1);
        CalcLocDerivShFnc( shapeFncDerivsAtIp_[lp.number], lp.coord,
                           NULL, 1);
      }
    }

    void FeH1LagrangeVar::SetFunctionsAtIp(const std::map<Integer,LocPoint>& iPoints){


      // can only be performed for non-hierarchical elements
      //shapeFncsAtIp_.resize(iPoints.GetSize());
      //shapeFncDerivsAtIp_.resize(iPoints.GetSize());
      std::map<Integer,LocPoint>::const_iterator pIt = iPoints.begin();
      while(pIt != iPoints.end()){
        const LocPoint& lp = pIt->second;
      CalcShFnc( shapeFncsAtIp_[lp.number], lp.coord, NULL, 1);
      CalcLocDerivShFnc( shapeFncDerivsAtIp_[lp.number], lp.coord,
                         NULL, 1);
        pIt++;
      }
    }
    
    void FeH1LagrangeVar::ComputeMonomialCoefficients(Matrix<Integer>& P, Matrix<Double>& C){
      //this algorithm is works for the lagrange tensorial polynomials
      //CAREFUL: Other shape functions may require different algorithms

      //start by  permutating the different orders to obtain the p matrix
      P.Resize(actNumFncs_,3);
      P.Init();

      UInt dof1 = (shape_.dim>0)? order_+1 : 1;
      UInt dof2 = (shape_.dim>1)? order_+1 : 1;
      UInt dof3 = (shape_.dim>2)? order_+1 : 1;
      UInt counter = 0;
      for(UInt i=0;i<dof1;i++){
        for(UInt j=0;j<dof2;j++){
          for(UInt k=0;k<dof3;k++){
            P[counter][0] = (Integer)i;
            P[counter][1] = (Integer)j;
            P[counter++][2] = (Integer)k;
          }
        }
      }

      //now we compute a temporary matrix by evaluating the sponentials
      // at our supporting point. each line one supporting point
      Matrix<Double> C_tmp(actNumFncs_,actNumFncs_);
      C.Resize(actNumFncs_,actNumFncs_);
      C_tmp.Init();
      C.Init();

      //now get the local coordinates in the correct ordering
      //we expect a actNumFncs_x3 matrix
      Matrix<Double> coords;
      GetLocalDOFCoordinates(coords);

      //now we compute the coefficients by an application of each coordinate
      for(UInt row = 0; row<actNumFncs_ ;row++){
        for(UInt col = 0; col<actNumFncs_;col++){
          C_tmp[row][col] = std::pow(coords[row][0],P[col][0]) * std::pow(coords[row][1],P[col][1]) * std::pow(coords[row][2],P[col][2]);
        }
      }

      C_tmp.Invert_Lapack();
      C_tmp.Transpose(C);

      // none-lapack variant
      // Matrix<Double> t;
      // C_tmp.Invert(t);
      // t.Transpose(C);
    }

  //=========================================================================
  //Line Elements of arbitrary order
  //=========================================================================
  
    FeH1LagrangeLineVar::FeH1LagrangeLineVar() : FeH1LagrangeVar() {
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

  void FeH1LagrangeLineVar::GetLocalDOFCoordinates(Matrix<Double> & coordMat){
    if(supportingPoints_.find(order_) == supportingPoints_.end() )
           supportingPoints_[order_] = CalcGaussLobattoPoints(order_);
    const StdVector<Double>& supPoints1D = supportingPoints_[order_];

    //This has to be a nx3 matrix even for 1D
    coordMat.Resize(actNumFncs_,3);
    coordMat.Init();

    //vertices:
    coordMat[0][0] = supPoints1D[0];
    coordMat[1][0] = supPoints1D[order_];
    
    //edges
    for(UInt i=1;i<order_;i++){
      coordMat[i+1][0] = supPoints1D[i];
    }
  }
  //=========================================================================
  //Quadrilateral elements of arbitrary order
  //=========================================================================

    FeH1LagrangeQuadVar::FeH1LagrangeQuadVar() : FeH1LagrangeVar(){
      feType_ = Elem::ET_QUAD4;
      shape_ = Elem::shapes[feType_];
      actNumFncs_ = 4;
      order_ = 1;
      hasICModes_ = true;
    }

    FeH1LagrangeQuadVar::~FeH1LagrangeQuadVar(){
    }

    void FeH1LagrangeQuadVar::CalcShFnc( Vector<Double>& Shape,
                                         const Vector<Double>& point,
                                         const Elem* ptElem,
                                         UInt comp ) {
      Shape.Resize( actNumFncs_ );
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
          deriv[offset][0] = shapeDerivX[j] * shapeY[i];
          deriv[offset++][1] = shapeX[j] * shapeDerivY[i];
        }
      }

    }

  void  FeH1LagrangeQuadVar::SetIsoOrder(UInt order){
    order_ = order;
    actNumFncs_ = (order+1)*(order+1);
  }

  void FeH1LagrangeQuadVar::GetLocalDOFCoordinates(Matrix<Double> & coordMat){
    if(supportingPoints_.find(order_) == supportingPoints_.end() )
           supportingPoints_[order_] = CalcGaussLobattoPoints(order_);
    const StdVector<Double>& supPoints1D = supportingPoints_[order_];

    //This has to be a nx3 matrix even for 1D
    coordMat.Resize(actNumFncs_,3);
    coordMat.Init();

    //vertices:
    coordMat[0][0] = supPoints1D[0];
    coordMat[0][1] = supPoints1D[0];

    coordMat[1][0] = supPoints1D[order_];
    coordMat[1][1] = supPoints1D[0];

    coordMat[2][0] = supPoints1D[order_];
    coordMat[2][1] = supPoints1D[order_];

    coordMat[3][0] = supPoints1D[0];
    coordMat[3][1] = supPoints1D[order_];

    UInt c = 4;
    //edge1
    for(UInt i=1;i<order_;i++){
      coordMat[c][0] = supPoints1D[i];
      coordMat[c++][1] = supPoints1D[0];
    }
    //edge2
    for(UInt i=1;i<order_;i++){
      coordMat[c][0] = supPoints1D[order_];
      coordMat[c++][1] = supPoints1D[i];
    }
    //edge3
    for(UInt i=1;i<order_;i++){
      coordMat[c][0] = supPoints1D[i];
      coordMat[c++][1] = supPoints1D[order_];
    }
    //edge4
    for(UInt i=1;i<order_;i++){
      coordMat[c][0] = supPoints1D[0];
      coordMat[c++][1] = supPoints1D[i];
    }

    //inner
    for(UInt i=1;i<order_;i++){
      for(UInt j=1;j<order_;j++){
        coordMat[c][0] = supPoints1D[j];
        coordMat[c++][1] = supPoints1D[i];
      }
    }
  }

  void FeH1LagrangeQuadVar::CalcShFncICModes( Vector<Double>& shape,
		  	                                  const Vector<Double>& point,
                                              const Elem* ptElem,
                                              UInt comp  ) {
	  shape.Resize( 2 );
	  shape[0] = 1.0 - point[0] * point[0];
	  shape[1] = 1.0 - point[1] * point[1];
  }


  void FeH1LagrangeQuadVar::CalcLocDerivShFncICModes( Matrix<Double> & deriv,
                                                   const Vector<Double>& point,
                                                   const Elem* ptElem,
                                                   UInt comp ) {
	   deriv.Resize( 2, 2);
	   deriv.Init();
	   deriv[0][0] = -2.0 * point[0];
	   deriv[1][1] = -2.0 * point[1];
  }

  //=========================================================================
  //Hexahedral elements of arbitrary order
  //=========================================================================
  
    FeH1LagrangeHexVar::FeH1LagrangeHexVar(): FeH1LagrangeVar() {
      feType_ = Elem::ET_HEXA8;
      shape_ = Elem::shapes[feType_];
    }

    FeH1LagrangeHexVar::~FeH1LagrangeHexVar(){
    }

//    UInt FeH1LagrangeVar::GetNumFncsPerEntType( EntityType fctEntityType,
//                                       UInt dof){
//      UInt numFncs = 0;
//      if( fctEntityType == VERTEX ) {
//        numFncs = shape_.numVertices;
//      } else if ( fctEntityType == EDGE) {
//        numFncs =  order_ - 1; 
//      } else if ( fctEntityType == FACE) {
//        numFncs = (order_ - 1)*(order_ - 1);
//      } else if ( fctEntityType == INTERIOR ) {
//        numFncs = (order_ - 1)*(order_ - 1)*(order_ - 1);
//      }
//      return numFncs;
//    }
    void FeH1LagrangeHexVar::CalcShFnc( Vector<Double>& Shape,
                                        const Vector<Double>& point,
                                        const Elem* ptElem,
                                        UInt comp ) {

      Shape.Resize( actNumFncs_ );
      Shape.Init();
      //now get the shape functions and the derivatives for the given coordinates
      Vector<Double> shapeX;
      Vector<Double> shapeY;
      Vector<Double> shapeZ;

      EvaluateLagrangePolynomial( shapeX, point[0], order_ );
      EvaluateLagrangePolynomial( shapeY, point[1], order_ );
      EvaluateLagrangePolynomial( shapeZ, point[2], order_ );

      UInt c = 0;
      Shape[c++] = shapeX[0]     * shapeY[0]     * shapeZ[0];
      Shape[c++] = shapeX[order_] * shapeY[0]     * shapeZ[0];
      Shape[c++] = shapeX[order_] * shapeY[order_] * shapeZ[0];
      Shape[c++] = shapeX[0]     * shapeY[order_] * shapeZ[0];
      Shape[c++] = shapeX[0]     * shapeY[0]     * shapeZ[order_];
      Shape[c++] = shapeX[order_] * shapeY[0]     * shapeZ[order_];
      Shape[c++] = shapeX[order_] * shapeY[order_] * shapeZ[order_];
      Shape[c++] = shapeX[0]     * shapeY[order_] * shapeZ[order_];


         // --------------------
         //  b) edge functions
         // --------------------
     #define FILL_EDGE( numEdge, shape, p1, p2)      \
         {                                          \
           for ( UInt i= 1; i<order_ ;i++ )    \
             Shape[c++] = shape[i] * p1 * p2;    \
         }

       FILL_EDGE( 1, shapeX, shapeY[0],     shapeZ[0] );
       FILL_EDGE( 2, shapeY, shapeX[order_], shapeZ[0] );
       FILL_EDGE( 3, shapeX, shapeY[order_], shapeZ[0] );
       FILL_EDGE( 4, shapeY, shapeX[0],     shapeZ[0] );

       FILL_EDGE( 5, shapeX, shapeY[0],     shapeZ[order_] );
       FILL_EDGE( 6, shapeY, shapeX[order_], shapeZ[order_] );
       FILL_EDGE( 7, shapeX, shapeY[order_], shapeZ[order_] );
       FILL_EDGE( 8, shapeY, shapeX[0],     shapeZ[order_] );
       
       FILL_EDGE( 9, shapeZ, shapeX[0],     shapeY[0] );
       FILL_EDGE(10, shapeZ, shapeX[order_], shapeY[0] );
       FILL_EDGE(11, shapeZ, shapeX[order_], shapeY[order_] );
       FILL_EDGE(12, shapeZ, shapeX[0],     shapeY[order_] );

       // --------------------
       //  c) face functions
       // --------------------

   #define FILL_FACE( numFace, shape1, shape2, fac)                   \
         {                                                            \
           for ( UInt i= 1;i< order_ ;i++ ){                           \
             for ( UInt j = 1; j< order_ ;j++ ){                       \
               Shape[c++] = shape1[j] * shape2[i] * fac;              \
             }                                                        \
           }                                                          \
         }

       FILL_FACE( 1, shapeX, shapeZ, shapeY[0] );
       FILL_FACE( 2, shapeY, shapeZ, shapeX[order_] );
       FILL_FACE( 3, shapeX, shapeZ, shapeY[order_] );
       FILL_FACE( 4, shapeY, shapeZ, shapeX[0] );
       FILL_FACE( 5, shapeX, shapeY, shapeZ[0] );
       FILL_FACE( 6, shapeX, shapeY, shapeZ[order_] );


       // ----------------------
       //  d) bubble functions
       // ----------------------
       for( UInt i = 1; i < order_ ; i++ ) {
         for( UInt j = 1; j < order_ ; j++ ) {
           for( UInt k = 1; k < order_ ; k++ ) {
             Shape[c++] = shapeX[k] * shapeY[j] * shapeZ[i];
           }
         }
       }
    }

    void FeH1LagrangeHexVar::CalcLocDerivShFnc( Matrix<Double> & deriv, 
                                                const Vector<Double>& point,
                                                const Elem* ptElem,
                                                UInt comp ) {
      Vector<Double> shapeX;
      Vector<Double> shapeY;
      Vector<Double> shapeZ;
      Vector<Double> shapeDerivX;
      Vector<Double> shapeDerivY;
      Vector<Double> shapeDerivZ;
      deriv.Resize( actNumFncs_ , shape_.dim );
      //now get the shape functions and the derivatives for the given coordinates
      EvaluateLagrangePolynomial( shapeX, point[0], order_ );
      EvaluateLagrangePolynomial( shapeY, point[1], order_ );
      EvaluateLagrangePolynomial( shapeZ, point[2], order_ );
      EvaluateDerivLagrangePolynomial( shapeDerivX, point[0], order_ );
      EvaluateDerivLagrangePolynomial( shapeDerivY, point[1], order_ );
      EvaluateDerivLagrangePolynomial( shapeDerivZ, point[2], order_ );
      UInt c = 0;

      deriv[c][0]   = shapeDerivX[0] *  shapeY[0]       * shapeZ[0];
      deriv[c][1]   = shapeX[0]      *  shapeDerivY[0]  * shapeZ[0];
      deriv[c++][2] = shapeX[0]      *  shapeY[0]       * shapeDerivZ[0];

      deriv[c][0]   = shapeDerivX[order_] *  shapeY[0]       * shapeZ[0];
      deriv[c][1]   = shapeX[order_]      *  shapeDerivY[0]  * shapeZ[0];
      deriv[c++][2] = shapeX[order_]      *  shapeY[0]       * shapeDerivZ[0];

      deriv[c][0] = shapeDerivX[order_]      *  shapeY[order_]       * shapeZ[0];
      deriv[c][1] = shapeX[order_]      *  shapeDerivY[order_]       * shapeZ[0];
      deriv[c++][2] = shapeX[order_]      *  shapeY[order_]       * shapeDerivZ[0];

      deriv[c][0] = shapeDerivX[0]      *  shapeY[order_]       * shapeZ[0];
      deriv[c][1] = shapeX[0]      *  shapeDerivY[order_]       * shapeZ[0];
      deriv[c++][2] = shapeX[0]      *  shapeY[order_]       * shapeDerivZ[0];

      deriv[c][0] = shapeDerivX[0]      *  shapeY[0]       * shapeZ[order_];
      deriv[c][1] = shapeX[0]      *  shapeDerivY[0]       * shapeZ[order_];
      deriv[c++][2] = shapeX[0]      *  shapeY[0]       * shapeDerivZ[order_];

      deriv[c][0] = shapeDerivX[order_]      *  shapeY[0]       * shapeZ[order_];
      deriv[c][1] = shapeX[order_]      *  shapeDerivY[0]       * shapeZ[order_];
      deriv[c++][2] = shapeX[order_]      *  shapeY[0]       * shapeDerivZ[order_];

      deriv[c][0] = shapeDerivX[order_]      *  shapeY[order_]       * shapeZ[order_];
      deriv[c][1] = shapeX[order_]      *  shapeDerivY[order_]       * shapeZ[order_];
      deriv[c++][2] = shapeX[order_]      *  shapeY[order_]       * shapeDerivZ[order_];

      deriv[c][0] = shapeDerivX[0]      *  shapeY[order_]       * shapeZ[order_];
      deriv[c][1] = shapeX[0]      *  shapeDerivY[order_]       * shapeZ[order_];
      deriv[c++][2] = shapeX[0]      *  shapeY[order_]       * shapeDerivZ[order_];


      // -------------------
      //  b) edge functions
      // -------------------
  #define HEX_E_DERIV(numEdge,idx1,idx2,idx3,dim)                                 \
        {                                                                            \
          UInt start = 1;                                                              \
          Integer inc = 1;                                                            \
          UInt end = order_;                                                          \
          switch(dim){                                                                \
          case 1:                                                                      \
            while( start != end)                                                      \
            {                                                                          \
              deriv[c][0]=   shapeDerivX[start] * shapeY[idx2] * shapeZ[idx3];          \
              deriv[c][1]=   shapeX[start] * shapeDerivY[idx2] * shapeZ[idx3];          \
              deriv[c++][2]= shapeX[start] * shapeY[idx2] * shapeDerivZ[idx3];        \
              start += inc;                                                            \
            }                                                                          \
            break;                                                                    \
          case 2:                                                                      \
            while( start != end)                                                      \
            {                                                                          \
              deriv[c][0]=   shapeDerivX[idx1] * shapeY[start] * shapeZ[idx3];          \
              deriv[c][1]=   shapeX[idx1] * shapeDerivY[start] * shapeZ[idx3];          \
              deriv[c++][2]= shapeX[idx1] * shapeY[start] * shapeDerivZ[idx3];        \
              start += inc;                                                            \
            }                                                                          \
            break;                                                                    \
          case 3:                                                                      \
            while( start != end)                                                      \
            {                                                                          \
              deriv[c][0]  = shapeDerivX[idx1] * shapeY[idx2] * shapeZ[start];          \
              deriv[c][1]  = shapeX[idx1] * shapeDerivY[idx2] * shapeZ[start];          \
              deriv[c++][2]= shapeX[idx1] * shapeY[idx2] * shapeDerivZ[start];        \
              start += inc;                                                            \
            }                                                                          \
            break;                                                                    \
          }                                                                            \
        }


      // EDGE #1
      HEX_E_DERIV(1, 0, 0, 0, 1);
      // EDGE #2
      HEX_E_DERIV(2, order_, 0, 0, 2);
       // EDGE #3
      HEX_E_DERIV(3, 0, order_, 0, 1 );
       // EDGE #4
      HEX_E_DERIV(4, 0, 0, 0, 2);
      
      // EDGE #9
      HEX_E_DERIV(9, 0, 0, order_, 1);
      // EDGE #10
      HEX_E_DERIV(10, order_, 0, order_, 2);
      // EDGE #11
      HEX_E_DERIV(11, 0, order_, order_, 1);
      // EDGE #12
      HEX_E_DERIV(12, 0, 0, order_, 2);
      
       // EDGE #5
      HEX_E_DERIV(5, 0, 0, 0, 3);
       // EDGE #6
      HEX_E_DERIV(6, order_, 0, 0, 3);
       // EDGE #7
      HEX_E_DERIV(7, order_, order_, 0, 3);
       // EDGE #8
      HEX_E_DERIV(8, 0, order_, 0, 3);
    

      // -------------------
      //  c) face functions
      // -------------------
  #define HEX_F_DERIV( numFace, shape1, shape1Deriv, shape2, shape2Deriv,              \
                          fac, facDeriv, dim)                                          \
    {                                                                                  \
      switch(dim){                                                                     \
        case 3:                                                                        \
          for ( UInt i= 1;i< order_ ;i++ ){                                             \
            for ( UInt j = 1; j< order_ ;j++ ){                                         \
              deriv[c][0]   = shape1Deriv[j] * shape2[i]      * fac;                  \
              deriv[c][1]   = shape1[j]      * shape2Deriv[i] * fac;                  \
              deriv[c++][2] = shape1[j]      * shape2[i]      * facDeriv;             \
            }                                                                          \
          }                                                                            \
          break;                                                                       \
        case 2:                                                                        \
          for ( UInt i= 1;i< order_ ;i++ ){                                             \
            for ( UInt j = 1; j< order_ ;j++ ){                                         \
              deriv[c][0]   = shape1Deriv[j] * shape2[i]      * fac;                  \
              deriv[c][1]   = shape1[j]      * shape2[i]      * facDeriv;             \
              deriv[c++][2] = shape1[j]      * shape2Deriv[i] * fac;                  \
            }                                                                          \
          }                                                                            \
          break;                                                                       \
        case 1:                                                                        \
          for ( UInt i= 1;i< order_ ;i++ ){                                             \
            for ( UInt j = 1; j< order_ ;j++ ){                                         \
              deriv[c][0]   = shape1[j]      * shape2[i]      * facDeriv;             \
              deriv[c][1]   = shape1Deriv[j] * shape2[i]      * fac;                  \
              deriv[c++][2] = shape1[j]      * shape2Deriv[i] * fac;                  \
            }                                                                          \
          }                                                                            \
          break;                                                                       \
      }                                                                                \
    }                                                                                  \

      
      HEX_F_DERIV( 1, shapeX, shapeDerivX, shapeZ, shapeDerivZ, shapeY[0],      shapeDerivY[0], 2 );
      HEX_F_DERIV( 2, shapeY, shapeDerivY, shapeZ, shapeDerivZ, shapeX[order_], shapeDerivX[order_], 1 );
      HEX_F_DERIV( 3, shapeX, shapeDerivX, shapeZ, shapeDerivZ, shapeY[order_], shapeDerivY[order_], 2 );
      HEX_F_DERIV( 4, shapeY, shapeDerivY, shapeZ, shapeDerivZ, shapeX[0],      shapeDerivX[0], 1 );
      HEX_F_DERIV( 5, shapeX, shapeDerivX, shapeY, shapeDerivY, shapeZ[0],      shapeDerivZ[0], 3 );
      HEX_F_DERIV( 6, shapeX, shapeDerivX, shapeY, shapeDerivY, shapeZ[order_], shapeDerivZ[order_], 3 );

      // ---------------------
      //  d) bubble functions
      // ---------------------
      for(UInt k = 1; k< order_ ; k++) {
        for(UInt j = 1; j< order_ ; j++) {
          for(UInt i = 1; i< order_ ; i++) {
            deriv[c][0]  = shapeDerivX[i] * shapeY[j]      * shapeZ[k];
            deriv[c][1]  = shapeX[i]      * shapeDerivY[j] * shapeZ[k];
            deriv[c++][2]= shapeX[i]      * shapeY[j]      * shapeDerivZ[k];
          }
        }
      }
    }

    void FeH1LagrangeHexVar::SetIsoOrder(UInt order){
      order_ = order;
      actNumFncs_ = (order+1)*(order+1)*(order+1);
    }

    void FeH1LagrangeHexVar::GetLocalDOFCoordinates(Matrix<Double> & coordMat){
        if(supportingPoints_.find(order_) == supportingPoints_.end() )
               supportingPoints_[order_] = CalcGaussLobattoPoints(order_);
        const StdVector<Double>& supPoints1D = supportingPoints_[order_];

        //This has to be a nx3 matrix even for 1D
        coordMat.Resize(actNumFncs_,3);
        coordMat.Init();

        //vertices:
        coordMat[0][0] = supPoints1D[0];
        coordMat[0][1] = supPoints1D[0];
        coordMat[0][2] = supPoints1D[0];

        coordMat[1][0] = supPoints1D[order_];
        coordMat[1][1] = supPoints1D[0];
        coordMat[1][2] = supPoints1D[0];

        coordMat[2][0] = supPoints1D[order_];
        coordMat[2][1] = supPoints1D[order_];
        coordMat[2][2] = supPoints1D[0];

        coordMat[3][0] = supPoints1D[0];
        coordMat[3][1] = supPoints1D[order_];
        coordMat[3][2] = supPoints1D[0];

        coordMat[4][0] = supPoints1D[0];
        coordMat[4][1] = supPoints1D[0];
        coordMat[4][2] = supPoints1D[order_];

        coordMat[5][0] = supPoints1D[order_];
        coordMat[5][1] = supPoints1D[0];
        coordMat[5][2] = supPoints1D[order_];

        coordMat[6][0] = supPoints1D[order_];
        coordMat[6][1] = supPoints1D[order_];
        coordMat[6][2] = supPoints1D[order_];

        coordMat[7][0] = supPoints1D[0];
        coordMat[7][1] = supPoints1D[order_];
        coordMat[7][2] = supPoints1D[order_];

        UInt c = 8;
        //edge1
        for(UInt i=1;i<order_;i++){
          coordMat[c][0] = supPoints1D[i];
          coordMat[c][1] = supPoints1D[0];
          coordMat[c++][2] = supPoints1D[0];
        }
        //edge2
        for(UInt i=1;i<order_;i++){
          coordMat[c][0] = supPoints1D[order_];
          coordMat[c][1] = supPoints1D[i];
          coordMat[c++][2] = supPoints1D[0];
        }
        //edge3
        for(UInt i=1;i<order_;i++){
          coordMat[c][0] = supPoints1D[i];
          coordMat[c][1] = supPoints1D[order_];
          coordMat[c++][2] = supPoints1D[0];
        }
        //edge4
        for(UInt i=1;i<order_;i++){
          coordMat[c][0] = supPoints1D[0];
          coordMat[c][1] = supPoints1D[i];
          coordMat[c++][2] = supPoints1D[0];
        }

        //edge9
        for(UInt i=1;i<order_;i++){
          coordMat[c][0] = supPoints1D[i];
          coordMat[c][1] = supPoints1D[0];
          coordMat[c++][2] = supPoints1D[order_];
        }
        //edge10
        for(UInt i=1;i<order_;i++){
          coordMat[c][0] = supPoints1D[order_];
          coordMat[c][1] = supPoints1D[i];
          coordMat[c++][2] = supPoints1D[order_];
        }
        //edge11
        for(UInt i=1;i<order_;i++){
          coordMat[c][0] = supPoints1D[i];
          coordMat[c][1] = supPoints1D[order_];
          coordMat[c++][2] = supPoints1D[order_];
        }
        //edge12
        for(UInt i=1;i<order_;i++){
          coordMat[c][0] = supPoints1D[0];
          coordMat[c][1] = supPoints1D[i];
          coordMat[c++][2] = supPoints1D[order_];
        }
        
        //edge5
        for(UInt i=1;i<order_;i++){
          coordMat[c][0] = supPoints1D[0];
          coordMat[c][1] = supPoints1D[0];
          coordMat[c++][2] = supPoints1D[i];
        }
        //edge6
        for(UInt i=1;i<order_;i++){
          coordMat[c][0] = supPoints1D[order_];
          coordMat[c][1] = supPoints1D[0];
          coordMat[c++][2] = supPoints1D[i];
        }
        //edge7
        for(UInt i=1;i<order_;i++){
          coordMat[c][0] = supPoints1D[order_];
          coordMat[c][1] = supPoints1D[order_];
          coordMat[c++][2] = supPoints1D[i];
        }
        //edge8
        for(UInt i=1;i<order_;i++){
          coordMat[c][0] = supPoints1D[0];
          coordMat[c][1] = supPoints1D[order_];
          coordMat[c++][2] = supPoints1D[i];
        }

        //face1
        for(UInt i=1;i<order_;i++){
          for(UInt j=1;j<order_;j++){
            coordMat[c][0] = supPoints1D[j];
            coordMat[c][1] = supPoints1D[0];
            coordMat[c++][2] = supPoints1D[i];
          }
        }

        //face2
        for(UInt i=1;i<order_;i++){
          for(UInt j=1;j<order_;j++){
            coordMat[c][0] = supPoints1D[order_];
            coordMat[c][1] = supPoints1D[j];
            coordMat[c++][2] = supPoints1D[i];
          }
        }

        //face3
        for(UInt i=1;i<order_;i++){
          for(UInt j=1;j<order_;j++){
            coordMat[c][0] = supPoints1D[j];
            coordMat[c][1] = supPoints1D[order_];
            coordMat[c++][2] = supPoints1D[i];
          }
        }

        //face4
        for(UInt i=1;i<order_;i++){
          for(UInt j=1;j<order_;j++){
            coordMat[c][0] = supPoints1D[0];
            coordMat[c][1] = supPoints1D[j];
            coordMat[c++][2] = supPoints1D[i];
          }
        }
        
        //face5
        for(UInt i=1;i<order_;i++){
          for(UInt j=1;j<order_;j++){
            coordMat[c][0] = supPoints1D[j];
            coordMat[c][1] = supPoints1D[i];
            coordMat[c++][2] = supPoints1D[0];
          }
        }


        //face6
        for(UInt i=1;i<order_;i++){
          for(UInt j=1;j<order_;j++){
            coordMat[c][0] = supPoints1D[j];
            coordMat[c][1] = supPoints1D[i];
            coordMat[c++][2] = supPoints1D[order_];
          }
        }

        //inner
        for(UInt i=1;i<order_;i++){
          for(UInt j=1;j<order_;j++){
            for(UInt k=1;k<order_;k++){
              coordMat[c][0] = supPoints1D[j];
              coordMat[c][1] = supPoints1D[i];
              coordMat[c++][2] = supPoints1D[k];
            }
          }
        }

      }
}
