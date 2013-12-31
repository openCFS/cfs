#ifndef BaseBOperator_HH
#define BaseBOperator_HH

#include <string>
#include "MatVec/Matrix.hh"

namespace CoupledField{

//! Basic Class for differential operators. These classes are passed as a 
//! template parameter to the forms object. Thereby a PDE rather specifies
//! an operator than a specific form.
class BaseBOperator{
public:
  
  //! Constructor
  BaseBOperator(){
    isSurfOpt_ = false;
  }

  //! Destructor
  virtual ~BaseBOperator(){
  }
  
  // =====================
  //  CALCULATION METHODS
  // =====================
  //@{ \name Calculation Methods
  
  //! Calc operator matrix
  virtual void CalcOpMat(Matrix<Double> & bMat,
                         const LocPointMapped& lp, 
                         BaseFE* ptFe ) = 0;

  //! Calculate transposed operator matrix
  virtual void CalcOpMatTransposed(Matrix<Double> & bMat, 
                                   const LocPointMapped& lp, 
                                   BaseFE* ptFe ) = 0;

  //! Calc complex valued operator matrix
  virtual void CalcOpMat(Matrix<Complex> & bMat,
                         const LocPointMapped& lp, BaseFE* ptFe ){
    Matrix<Double> realMat;
    this->CalcOpMat(realMat,lp,ptFe);
    const UInt nrow = realMat.GetNumRows();
    const UInt ncol = realMat.GetNumCols();
    bMat.Resize(nrow,ncol);
    bMat.SetPart(Global::REAL, realMat, true );
  }

  //! Calculate transposed complex valued operator matrix
  virtual void CalcOpMatTransposed(Matrix<Complex> & bMat,
                                   const LocPointMapped& lp, 
                                   BaseFE* ptFe ){
    Matrix<Double> realMat;
    this->CalcOpMatTransposed(realMat,lp,ptFe);
    const UInt nrow = realMat.GetNumRows();
    const UInt ncol = realMat.GetNumCols();
    bMat.Resize(nrow,ncol);
    bMat.SetPart(Global::REAL, realMat, true);
  }

  //! Apply the operator matrix on a vector
  virtual void ApplyOp(Vector<Double>& retVec,
                       const LocPointMapped& lp,
                       BaseFE* ptFe,
                       const Vector<Double>& solVec ){
    Matrix<Double> bOp;
    CalcOpMat(bOp,lp,ptFe);
    retVec = bOp * solVec;

  }
  
  virtual void ApplyOp(Vector<Complex>& retVec,
                       const LocPointMapped& lp,
                       BaseFE* ptFe,
                       const Vector<Complex>& solVec ){
    Matrix<Double> bOp;
    CalcOpMat(bOp,lp,ptFe);
    retVec = bOp * solVec;
    }

  //! Apply the transposed operator matrix on a vector
  virtual void ApplyOpTranspose(Vector<Double>& retVec,
                                const LocPointMapped& lp,
                                BaseFE* ptFe,
                                const Vector<Double>& solVec ){
    Matrix<Double> bOp;
    CalcOpMat(bOp,lp,ptFe);
    retVec = Transpose(bOp) * solVec;
  }
  
  //! Apply the transposed operator matrix on a vector
  virtual void ApplyOpTranspose(Vector<Complex>& retVec,
                                const LocPointMapped& lp,
                                BaseFE* ptFe,
                                const Vector<Complex>& solVec ){
    Matrix<Double> bOp;
    CalcOpMat(bOp,lp,ptFe);
    retVec = Transpose(bOp) * solVec;
  }
  
  //@}
  
  // ===============
  //  QUERY METHODS
  // ===============
  //@{ \name Query Methods
  
  //! Return differential order of the operator
  virtual UInt GetDiffOrder() const = 0;
  
  //! Return number of components of the problem (scalar, vector)
  virtual UInt GetDimDof() const = 0;
  
  //! Return number of spatial dimensions of the underlying space
  virtual UInt GetDimSpace() const = 0;
  
  //! Return dimension of the finite element (different for surface elements)
  virtual UInt GetDimElem() const = 0;
  
  //! Return dimension of the related material tensor
  virtual UInt GetDimDMat() const = 0;

  //! Return name of the integrator
  std::string GetName(){
    return name_;
  }
  
  //! Set the coefficient function of the operator
  virtual void SetCoefFunction(PtrCoefFct coef){
    coef_ = coef;
  }

  //! Set the coefficient function of the operator
  virtual PtrCoefFct GetCoefFunction(){
    return coef_;
  }

  //! set operator to rurface operator
  virtual void SetOperator2SurfOperator() {
    isSurfOpt_ = true;
  }
  
protected:

  //! Name of the integrator
  std::string name_;

  //!pointer to coefficient function as used e.g. in Convective operators
  PtrCoefFct coef_;

  //! if true, then it is a surface operator
  bool isSurfOpt_;
};
}
#endif
