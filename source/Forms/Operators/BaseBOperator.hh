#ifndef BaseBOperator_HH
#define BaseBOperator_HH

#include <string>
#include "MatVec/Matrix.hh"

namespace CoupledField{

//! Basic Class for differential operators. These classes are passed as a 
//! template parameter to the forms object. Thereby a PDE rather specifies
//! an operator than a specific form.

class BaseBOperator : public CfsCopyable{
public:
  
  //! Constructor
  BaseBOperator(){
    isSurfOpt_ = false;
  }

  BaseBOperator(const BaseBOperator & other){
    this->name_ = other.name_;
    this->coef_ = other.coef_;
    this->isSurfOpt_ = other.isSurfOpt_;
  }

  virtual BaseBOperator * Clone() = 0;

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
    // bMat.Init();
    // bMat = realMat * Complex(1.0, 0.0);
    bMat.SetPart(Global::REAL, realMat, true);
  }

  //! Apply the operator matrix on a vector.
  //! The vector consists of the weights for the operator matrix
  //! and is stacked for non-scalar values
  //! For example: using identityOperator at numFcns = 4 and N_DOF = 3:
  //! [x1 y1 z1 x2 y2 z2 x3 y3 z3 x4 y4 z4]^T
  //! \param retVec (out) vector containing the operation result at the local point mapped
  //! \param lp (in) the local point mapped
  //! \param solVec (in) the stacked input vector
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
  const std::string& GetName() const {
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

  //! set operator to rurface operator
  virtual void OverrideIsSurfOperator(bool overrideIsSurfOpt) {
    overrideIsSurfOpt_ = overrideIsSurfOpt;
  }

  //! set operator to rurface operator
  virtual void ResetOperator2SurfOperator() {
    isSurfOpt_ = false;
  }

  
protected:

  //! Name of the operator
  std::string name_;

  //!pointer to coefficient function as used e.g. in Convective operators
  PtrCoefFct coef_;

  //! if true, then it is a surface operator
  bool isSurfOpt_;

  //! if true, evaluate as surface operator regardless of actual definition
  bool overrideIsSurfOpt_ = false;

};
}
#endif
