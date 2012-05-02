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
    UInt nrow = realMat.GetNumRows();
    UInt ncol = realMat.GetNumCols();
    bMat.Resize(nrow,ncol);
    bMat.Init();
    for(UInt i=0;i<nrow;++i){
      for(UInt j=0;j<ncol;++j){
        bMat[i][j] = Complex(1.0,0.0) * realMat[i][j];
      }
    }
  }

  //! Calculate transposed complex valued operator matrix
  virtual void CalcOpMatTransposed(Matrix<Complex> & bMat,
                                   const LocPointMapped& lp, 
                                   BaseFE* ptFe ){
    Matrix<Double> realMat;
    this->CalcOpMatTransposed(realMat,lp,ptFe);
    UInt nrow = realMat.GetNumRows();
    UInt ncol = realMat.GetNumCols();
    bMat.Resize(nrow,ncol);
    
    bMat.Init();
    bMat = realMat * Complex(1.0, 0.0);
//    for(UInt i=0;i<nrow;++i){
//      for(UInt j=0;j<ncol;++j){
//        bMat[i][j] = Complex(1.0,0.0) * realMat[i][j];
//      }
//    }
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
  
  //OBSOLETE!!!!
  //! Additional transformation of the Jacobian determinant (e.g. for PML)
  virtual void TransformJacDet(Double & jacLoc,
                               const LocPointMapped & lp, BaseFE* ptFe){
    return;
  }
  
  //! Additional transformation of the Jacobian determinant (e.g. for PML)
  virtual void TransformJacDet(Complex & jacLoc,
                               const LocPointMapped & lp, BaseFE* ptFe){
    return;
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

  
protected:

  //! Name of the integrator
  std::string name_;

  //!pointer to coefficient function as used e.g. in Convective operators
  PtrCoefFct coef_;
};
}
#endif
