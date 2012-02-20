// =====================================================================================
// 
//       Filename:  identityOperatorNormalTrans.hh
// 
//    Description:  This class implements the identity operator just
//                  consisting of the elements shape functions evaulated at the 
//                  given local coordinate
//                  b = ( N_1*n_x, N_1*n_y, N1_*n_z, N_2*n_x, N_2*n_y, N2_*n_z, .. )
//
//                  for the basic shape functions N_1x = N_1y = N_1z
//                  for scalar unknowns the operator is just a vector
// 
//        Version:  1.0
//        Created:  10/04/2011 09:10:00 AM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================

#ifndef IDENTITYOP_NORMAL_TRANS_HH
#define IDENTITYOP_NORMAL_TRANS_HH
#include "BaseBOperator.hh"

namespace CoupledField{

//! Calculates the transposed normal-projected identity operator for vector functions

//! This class implements the transposed identity operator for vectorial shape 
//! functions,  which gets projected in normal direction. The element matrix
//! is computed as:  
//! normal direction. The element matrix is computed as 

//! b = ( N_1*n_x, N_1*n_y, N1_*n_z, N_2*n_x, N_2*n_y, N2_*n_z, .. )
//!   = n^T * N
//!
//! and is of size (1 x ( DIM_SPACE * number of functions).

  template<class FE, UInt D ,  UInt D_DOF = 1, class TYPE = Double>
  class IdentityOperatorNormalTrans : public BaseBOperator<TYPE>{

  public:
    
    // ------------------
    //  STATIC CONSTANTS 
    // ------------------
    //@{ 
    //! \name Static constants

    //! Order of differentiation
    static const UInt ORDER_DIFF = 0;

    //! Number of components of the problem (scalar, vector)
    static const UInt DIM_DOF = D_DOF;

    //! Dimension of the underlying domain / space
    static const UInt DIM_SPACE = D;

    //! Dimension of the finite element
    static const UInt DIM_ELEM = D-1;

    //! Dimension of the related material 
    static const UInt DIM_D_MAT = 1; 
    //@}


    IdentityOperatorNormalTrans(){
      return;
    }

    virtual ~IdentityOperatorNormalTrans(){
      return;
    }

    virtual void CalcOpMat(Matrix<Double> & bMat,
                           const LocPointMapped& lp, BaseFE* ptFe );

    virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                     const LocPointMapped& lp, BaseFE* ptFe );

    //avoid reimplementation of complex operator by making the bas class function
    //available
    using BaseBOperator<TYPE>::CalcOpMat;

    using BaseBOperator<TYPE>::CalcOpMatTransposed;

  protected:

};
  
  template<class FE,  UInt D, UInt D_DOF, class TYPE>
  void IdentityOperatorNormalTrans<FE,D,D_DOF,TYPE>::
  CalcOpMat(Matrix<Double> & bMat,
            const LocPointMapped& lp, BaseFE* ptFe){
    
    // ensure, that the surface information (i.e. normal direction)
    // is set at the mapped local point
    assert(lp.isSurface);

    UInt numFncs = ptFe->GetNumFncs();

    // Set correct size of matrix B and initialize with zeros
    bMat.Resize( 1, DIM_SPACE * numFncs );
    bMat.Init();

    Vector<Double> s;
    FE *fe = (static_cast<FE*>(ptFe));
    for(UInt d = 0; d < DIM_SPACE; d ++){
      fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem() , d );
      for(UInt sh = 0; sh < numFncs; sh ++){
        bMat[0][DIM_SPACE*sh+d] = s[sh] * lp.normal[d];
      }
    }

  }

  template<class FE,  UInt D, UInt D_DOF, class TYPE>
  void IdentityOperatorNormalTrans<FE,D,D_DOF,TYPE>::
  CalcOpMatTransposed(Matrix<Double> & bMat,
                      const LocPointMapped& lp, BaseFE* ptFe){
    
    // ensure, that the surface information (i.e. normal direction)
    // is set at the mapped local point
    assert(lp.isSurface);
    
    UInt numFncs = ptFe->GetNumFncs();
    // Set correct size of matrix B and initialize with zeros
    bMat.Resize( DIM_SPACE * numFncs, 1 );
    bMat.Init();

    Vector<Double> s;
    FE *fe = (static_cast<FE*>(ptFe));
    for(UInt d = 0; d < DIM_SPACE; d ++){
      fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem() , d );
      for(UInt sh = 0; sh < numFncs; sh ++){
        bMat[DIM_SPACE*sh+d][0] = s[sh] * lp.normal[d];
      }
    }

  }
} // end of namespace
#endif
