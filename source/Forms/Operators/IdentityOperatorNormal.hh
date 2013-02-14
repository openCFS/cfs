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
  template<class FE, UInt D , class TYPE = Double>
  class IdentityOperatorNormal : public BaseBOperator{

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
  
  template<class FE, UInt D, class TYPE>
  void IdentityOperatorNormal<FE,D,TYPE>::
  CalcOpMat(Matrix<Double> & bMat,
            const LocPointMapped& lp, BaseFE* ptFe){

    // ensure, that the surface information (i.e. normal direction)
    // is set at the mapped local point
    assert(lp.isSurface);
    
    const UInt numFncs = ptFe->GetNumFncs();

    // Set correct size of matrix B and initialize with zeros
    bMat.Resize( DIM_SPACE, numFncs);

    Vector<Double> s;
    FE *fe = (static_cast<FE*>(ptFe));
    for(UInt d = 0; d < DIM_SPACE; d++){
      fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem() , d );
      for(UInt sh = 0; sh < numFncs; sh ++){
        bMat[d][sh] = s[sh] * lp.normal[d];
      }
    }

  }

  template<class FE, UInt D, class TYPE>
  void IdentityOperatorNormal<FE,D,TYPE>::
  CalcOpMatTransposed(Matrix<Double> & bMat,
                      const LocPointMapped& lp, BaseFE* ptFe){
    // ensure, that the surface information (i.e. normal direction)
    // is set at the mapped local point
    assert(lp.isSurface);
    
    const UInt numFncs = ptFe->GetNumFncs();
    // Set correct size of matrix B and initialize with zeros
    bMat.Resize( numFncs , DIM_SPACE );

    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Vector<Double> s;
    FE *fe = (static_cast<FE*>(ptFe));
    for(UInt d = 0; d < DIM_SPACE ; d ++){
      fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem() , d );
      for(UInt sh = 0; sh < numFncs; sh++){
        bMat[sh][d] = s[sh] * lp.normal[d];
      }
    }

  }
  template<class FE, UInt D ,  UInt DIM_D,  class TYPE = Double>
  class IdentityOperatorPiolaNormal : public IdentityOperatorNormal<FE,D,TYPE>{

  public:

    // ------------------
    //  STATIC CONSTANTS
    // ------------------
    //@{
    //! \name Static constants

    //! Order of differentiation
    static const UInt ORDER_DIFF = 0;

    //! Number of components of the problem (scalar, vector)
    static const UInt DIM_DOF = DIM_D;

    //! Dimension of the underlying domain / space
    static const UInt DIM_SPACE = D;

    //! Dimension of the finite element
    static const UInt DIM_ELEM = D-1;

    //! Dimension of the related material
    static const UInt DIM_D_MAT = 1;
    //@}


    IdentityOperatorPiolaNormal(){
      return;
    }

    virtual ~IdentityOperatorPiolaNormal(){
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

  protected:

};
  template<class FE,  UInt D,  UInt DIM_D, class TYPE>
  void IdentityOperatorPiolaNormal<FE,D,DIM_D,TYPE>::
  CalcOpMat(Matrix<Double> & bMat,
            const LocPointMapped& lp, BaseFE* ptFe){

    //TODO: RIGHT NOW VERY CIRCUMFENCIAL WELL SHORTEN THIS LATER
    // ensure, that the surface information (i.e. normal direction)
    // is set at the mapped local point
    assert(lp.isSurface);

    const UInt numFncs = ptFe->GetNumFncs();

    //first get usual identity OP
    Matrix<Double> bId;
    Matrix<Double> bIdPiola;
    // Set correct size of matrix B and initialize with zeros
    bId.Resize( DIM_DOF, numFncs * DIM_DOF );
    bId.Init();
    bIdPiola.Resize( DIM_DOF, numFncs * DIM_DOF );
    bIdPiola.Init();

    Vector<Double> s;
    FE *fe = (static_cast<FE*>(ptFe));
    for(UInt d = 0; d < DIM_DOF ; d ++){
      fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem() , d );
      for(UInt sh = 0; sh < numFncs; sh ++){
        bId[d][sh*DIM_DOF + d] = s[sh];
      }
    }
    //now apply piola tranform
#ifdef USE_BLAS_VERSION
    Double jacDetInv = (1.0/lp.lpmVol->jacDet);
    lp.lpmVol->jac.Mult_Blas(bId,bIdPiola,false,false,jacDetInv,0.0);
#else
    bIdPiola = lp.lpmVol->jac * bId;
    bIdPiola *= (1.0/lp.lpmVol->jacDet);
#endif

    // Set correct size of matrix B and initialize with zeros
    bMat.Resize( 1, numFncs*DIM_SPACE);
    bMat.Init();

    for(UInt d = 0; d < DIM_SPACE; d ++){
      for(UInt sh = 0; sh < numFncs; sh ++){
        for(UInt locDim = 0; locDim < DIM_SPACE; locDim++){
          bMat[0][sh*DIM_SPACE+d] += bIdPiola[d][sh*DIM_DOF + locDim] * lp.normal[locDim];
        }
      }
    }
  }

  template<class FE,  UInt D,  UInt DIM_D, class TYPE>
  void IdentityOperatorPiolaNormal<FE,D,DIM_D,TYPE>::
  CalcOpMatTransposed(Matrix<Double> & bMat,
                      const LocPointMapped& lp, BaseFE* ptFe){
    //TODO: RIGHT NOW VERY CIRCUMFENCIAL WELL SHORTEN THIS LATER
    // ensure, that the surface information (i.e. normal direction)
    // is set at the mapped local point
    assert(lp.isSurface);
    EXCEPTION("DO NOT USE THIS!");

    const UInt numFncs = ptFe->GetNumFncs();

    //first get usual identity OP
    Matrix<Double> bId;
    Matrix<Double> bIdPiola;
    // Set correct size of matrix B and initialize with zeros
    bId.Resize( numFncs * DIM_DOF, DIM_DOF );
    bId.Init();
    bIdPiola.Resize(  numFncs * DIM_DOF ,DIM_DOF );
    bIdPiola.Init();

    Vector<Double> s;
    FE *fe = (static_cast<FE*>(ptFe));
    for(UInt d = 0; d < DIM_DOF ; d ++){
      fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem() , d );
      for(UInt sh = 0; sh < numFncs; sh ++){
        bId[sh*DIM_DOF + d][d] = s[sh];
      }
    }

    //now apply piola tranform
#ifdef USE_BLAS_VERSION
          Double jacDetInv = (1.0/lp.lpmVol->jacDet);
          bId.Mult_Blas(lp.lpmVol->jac,bIdPiola,false,true,jacDetInv,0.0);
#else
          Matrix<Double> jacTmp;
          lp.lpmVol->jac.Transpose(jacTmp);
          bIdPiola = bId * jacTmp;
          bIdPiola *= (1.0/lp.lpmVol->jacDet);
#endif

    // Set correct size of matrix B and initialize with zeros
    bMat.Resize( numFncs , DIM_SPACE );
    bMat.Init();

    for(UInt d = 0; d < DIM_SPACE; d ++){
      for(UInt sh = 0; sh < numFncs; sh ++){
        for(UInt locDim = 0; locDim < DIM_SPACE; locDim++){
          bMat[sh][d] += bIdPiola[sh*DIM_DOF + d][locDim] * lp.normal[d];
        }
      }
    }
  }


} // end of namespace
#endif
