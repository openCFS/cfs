#include "BaseBOperator.hh"

namespace CoupledField{
  
  //! calculates the projection of a vector valued operator onto a vector, i.e. the inner product
  
  //! The element matrix is computed as 
  //! b = ( N1_x*v_x  N1_y*n_y  .. N2_x*v_x N2_y*v_y ...) 
  //! and is of size (1 x Number of functions*D_DOF).
  
  //Inner Product with Vector:
  template<class FE, UInt D , class TYPE = Double>
  class IdentityOperatorInVector : public BaseBOperator{

  public:

    // ------------------
    //  STATIC CONSTANTS
    // ------------------
    //@{
    //! \name Static constants

    //! Order of differentiation
    static const UInt ORDER_DIFF = 0;

    //! Number of components of the problem (scalar, vector)
    static const UInt DIM_DOF = D;

    //! Dimension of the underlying domain / space
    static const UInt DIM_SPACE = D;

    //! Dimension of the finite element
    static const UInt DIM_ELEM = D-1;

    //! Dimension of the related material
    static const UInt DIM_D_MAT = 1;
    //@}


    IdentityOperatorInVector(Vector<Double> vector){
      vector_ = vector;
      return;
    }

    virtual ~IdentityOperatorInVector(){
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

    Vector<Double> vector_;

};

  template<class FE, UInt D, class TYPE>
  void IdentityOperatorInVector<FE,D,TYPE>::
  CalcOpMat(Matrix<Double> & bMat,
            const LocPointMapped& lp, BaseFE* ptFe){

    const UInt numFncs = ptFe->GetNumFncs();

    // Set correct size of matrix B and initialize with zeros
    bMat.Resize( 1, numFncs*DIM_DOF);
    
    Vector<Double> s;
    FE *fe = (static_cast<FE*>(ptFe));
    for(UInt sh = 0; sh < numFncs; sh ++){
      for (UInt idof=0; idof< DIM_DOF; idof++){
        fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem() , idof );
        bMat[0][sh*DIM_DOF+idof] = s[sh] * vector_[idof];
      }
    }
  }

  template<class FE, UInt D, class TYPE>
  void IdentityOperatorInVector<FE,D,TYPE>::
  CalcOpMatTransposed(Matrix<Double> & bMat,
                      const LocPointMapped& lp, BaseFE* ptFe){

    const UInt numFncs = ptFe->GetNumFncs();
    // Set correct size of matrix B and initialize with zeros
    bMat.Resize( numFncs*DIM_DOF, 1 );
    
    Vector<Double> s;
    FE *fe = (static_cast<FE*>(ptFe));
    for(UInt sh = 0; sh < numFncs; sh ++){
      for (UInt idof=0; idof< DIM_DOF; idof++){
        fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem() , idof );
        bMat[sh*DIM_DOF+idof][0] = s[sh] * vector_[idof];
      }
    }
  }
  
  //Inner Product with Normal Vector:
  template<class FE, UInt D , class TYPE = Double>
  class IdentityOperatorInNormal : public BaseBOperator{

  public:

    // ------------------
    //  STATIC CONSTANTS
    // ------------------
    //@{
    //! \name Static constants

    //! Order of differentiation
    static const UInt ORDER_DIFF = 0;

    //! Number of components of the problem (scalar, vector)
    static const UInt DIM_DOF = D;

    //! Dimension of the underlying domain / space
    static const UInt DIM_SPACE = D;

    //! Dimension of the finite element
    static const UInt DIM_ELEM = D-1;

    //! Dimension of the related material
    static const UInt DIM_D_MAT = 1;
    //@}


    IdentityOperatorInNormal(){
      return;
    }

    virtual ~IdentityOperatorInNormal(){
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

    Vector<Double> vector_;

};

  template<class FE, UInt D, class TYPE>
  void IdentityOperatorInNormal<FE,D,TYPE>::
  CalcOpMat(Matrix<Double> & bMat,
            const LocPointMapped& lp, BaseFE* ptFe){

    const UInt numFncs = ptFe->GetNumFncs();

    // Set correct size of matrix B and initialise with zeros
    bMat.Resize( 1, numFncs*DIM_DOF);
    
    Vector<Double> s;
    FE *fe = (static_cast<FE*>(ptFe));
    for(UInt sh = 0; sh < numFncs; sh ++){
      for (UInt idof=0; idof< DIM_DOF; idof++){
        fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem() , idof );
        bMat[0][sh*DIM_DOF+idof] = s[sh] * lp.normal[idof];
      }
    }
  }

  template<class FE, UInt D, class TYPE>
  void IdentityOperatorInNormal<FE,D,TYPE>::
  CalcOpMatTransposed(Matrix<Double> & bMat,
                      const LocPointMapped& lp, BaseFE* ptFe){

    const UInt numFncs = ptFe->GetNumFncs();
    // Set correct size of matrix B and initialise with zeros
    bMat.Resize( numFncs*DIM_DOF, 1 );
    
    Vector<Double> s;
    FE *fe = (static_cast<FE*>(ptFe));
    for(UInt sh = 0; sh < numFncs; sh ++){
      for (UInt idof=0; idof< DIM_DOF; idof++){
        fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem() , idof );
        bMat[sh*DIM_DOF+idof][0] = s[sh] * lp.normal[idof];
      }
    }
  }
  
} // end of name space
