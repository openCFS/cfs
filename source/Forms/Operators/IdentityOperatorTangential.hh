#ifndef IDENTITYOP_TANGENTIAL_HH
#define IDENTITYOP_TANGENTIAL_HH

#include "BaseBOperator.hh"

namespace CoupledField{
  
  //! Calculates the tangential-projected identity operator.

  //! This class implements the identity operator, which gets projected in 
  //! tangential direction (90° rotation of normal direction). The element matrix is computed as 
  //!     / -N_1*n_y  -N_2*n_y  .. \
  //! b = |                     ..  | = \vec{t} * N
  //!     \  N_1*n_x   N_2*n_x  .. /
  //!
  //! and is of size (DIM_SPACE x Number of functions).
  template<class FE, UInt D = 2 , UInt D_DOF = 1, class TYPE = Double>
  class IdentityOperatorTangential2D : public BaseBOperator{

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


    IdentityOperatorTangential2D(){
      return;
    }

    //! Copy constructor
    IdentityOperatorTangential2D(const IdentityOperatorTangential2D & other)
       : BaseBOperator(other){
    }

    //! \copydoc BaseBOperator::Clone()
    virtual IdentityOperatorTangential2D * Clone(){
      return new IdentityOperatorTangential2D(*this);
    }

    virtual ~IdentityOperatorTangential2D(){
      return;
    }

    virtual void CalcOpMat(Matrix<Double> & bMat,
                           const LocPointMapped& lp, BaseFE* ptFe );

    virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                     const LocPointMapped& lp, BaseFE* ptFe );

    //avoid reimplementation of complex operator by making the bas class function
    //available
    using BaseBOperator::CalcOpMat;

    using BaseBOperator::CalcOpMatTransposed;

    // ===============
    //  QUERY METHODS
    // ===============
    //@{ \name Query Methods
    //! \copydoc BaseBOperator::GetDiffOrder
    virtual UInt GetDiffOrder() const {
      return ORDER_DIFF;
    }

    //! \copydoc BaseBOperator::GetDimDof()
    virtual UInt GetDimDof() const {
      return DIM_DOF;
    }

    //! \copydoc BaseBOperator::GetDimSpace()
    virtual UInt GetDimSpace() const {
      return DIM_SPACE;
    }

    //! \copydoc BaseBOperator::GetDimElem()
    virtual UInt GetDimElem() const {
      return DIM_ELEM;
    }

    //! \copydoc BaseBOperator::GetDimDMat()
    virtual UInt GetDimDMat() const {
      return DIM_D_MAT;
    }
    //@}
    
  protected:

};
  
  template<class FE, UInt D, UInt D_DOF, class TYPE>
  void IdentityOperatorTangential2D<FE,D,D_DOF,TYPE>::
  CalcOpMat(Matrix<Double> & bMat,
            const LocPointMapped& lp, BaseFE* ptFe){

    EXCEPTION("IdentityOperatorTangential is not tested, please verify and remove this exception!");
    
    // ensure, that the surface information (i.e. normal direction)
    // is set at the mapped local point
    assert(lp.isSurface);
    // check if we are in 2D
    assert(D==2);
    
    const UInt numFncs = ptFe->GetNumFncs();

    // Set correct size of matrix B and initialize with zeros
    bMat.Resize( DIM_SPACE, numFncs);

    Vector<Double> s;
    FE *fe = (static_cast<FE*>(ptFe));
    fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem() , 0 );
    for(UInt sh = 0; sh < numFncs; sh ++){
      bMat[0][sh] = -s[sh] * lp.normal[1];
    }
    fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem() , 1 );
    for(UInt sh = 0; sh < numFncs; sh ++){
      bMat[1][sh] = s[sh] * lp.normal[0];
    }
  }

  template<class FE, UInt D, UInt D_DOF, class TYPE>
  void IdentityOperatorTangential2D<FE,D,D_DOF,TYPE>::
  CalcOpMatTransposed(Matrix<Double> & bMat,
                      const LocPointMapped& lp, BaseFE* ptFe){
    
    EXCEPTION("IdentityOperatorTangential is not tested, please verify and remove this exception!");
    
    // ensure, that the surface information (i.e. normal direction)
    // is set at the mapped local point
    assert(lp.isSurface);
    // check if we are in 2D
    assert(D==2);
    
    const UInt numFncs = ptFe->GetNumFncs();
    // Set correct size of matrix B and initialize with zeros
    bMat.Resize( numFncs, DIM_SPACE );

    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Vector<Double> s;
    FE *fe = (static_cast<FE*>(ptFe));
    fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem() , 0 );
    for(UInt sh = 0; sh < numFncs; sh++){
      bMat[sh][0] = -s[sh] * lp.normal[1];
    }
    fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem() , 1 );
    for(UInt sh = 0; sh < numFncs; sh++){
      bMat[sh][1] = s[sh] * lp.normal[0];
    }
  }
} // end of namespace
#endif
