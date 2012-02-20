#ifndef IDENTITYOP_NORMAL_HH
#define IDENTITYOP_NORMAL_HH

#include "BaseBOperator.hh"

namespace CoupledField{
  
  //! Calculates the normal-projected identity operator.

  //! This class implements the identity operator, which gets projected in 
  //! normal direction. The element matrix is computed as 
  //!     / N_1*n_x  N_2*n_x  .. \
  //! b = | N_1*n_y  N_2*n_y  ..  | = \vec{n} * N
  //!     \ N_1*n_z  N_2*n_z  .. /
  //!
  //! and is of size (DIM_SPACE x Number of functions).
  template<class FE, UInt D ,  class TYPE = Double>
  class IdentityOperatorNormal : public BaseBOperator<TYPE>{

  public:
    
    // ------------------
    //  STATIC CONSTANTS 
    // ------------------
    //@{ 
    //! \name Static constants

    //! Order of differentiation
    static const UInt ORDER_DIFF = 0;

    //! Number of components of the problem (scalar, vector)
    static const UInt DIM_DOF = 1;

    //! Dimension of the underlying domain / space
    static const UInt DIM_SPACE = D;

    //! Dimension of the finite element
    static const UInt DIM_ELEM = D-1;

    //! Dimension of the related material 
    static const UInt DIM_D_MAT = 1; 
    //@}


    IdentityOperatorNormal(){
      return;
    }

    virtual ~IdentityOperatorNormal(){
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
  
  template<class FE,  UInt D, class TYPE>
  void IdentityOperatorNormal<FE,D,TYPE>::
  CalcOpMat(Matrix<Double> & bMat,
            const LocPointMapped& lp, BaseFE* ptFe){

    // ensure, that the surface information (i.e. normal direction)
    // is set at the mapped local point
    assert(lp.isSurface);
    
    UInt numFncs = ptFe->GetNumFncs();

    // Set correct size of matrix B and initialize with zeros
    bMat.Resize( DIM_SPACE, numFncs);
    bMat.Init();

    Vector<Double> s;
    FE *fe = (static_cast<FE*>(ptFe));
    for(UInt d = 0; d < DIM_SPACE; d ++){
      fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem() , d );
      for(UInt sh = 0; sh < numFncs; sh ++){
        bMat[d][sh] = s[sh] * lp.normal[d];
      }
    }

  }

  template<class FE,  UInt D, class TYPE>
  void IdentityOperatorNormal<FE,D,TYPE>::
  CalcOpMatTransposed(Matrix<Double> & bMat,
                      const LocPointMapped& lp, BaseFE* ptFe){
    // ensure, that the surface information (i.e. normal direction)
    // is set at the mapped local point
    assert(lp.isSurface);
    
    UInt numFncs = ptFe->GetNumFncs();
    // Set correct size of matrix B and initialize with zeros
    bMat.Resize( numFncs , DIM_SPACE );
    bMat.Init();

    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Vector<Double> s;
    FE *fe = (static_cast<FE*>(ptFe));
    for(UInt d = 0; d < DIM_SPACE ; d ++){
      fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem() , d );
      for(UInt sh = 0; sh < numFncs; sh ++){
        bMat[sh][d] = s[sh] * lp.normal[d];
      }
    }

  }
} // end of namespace
#endif
