#ifndef BaseBOperator_HH
#define BaseBOperator_HH

#include <string>
#include "MatVec/Matrix.hh"

namespace CoupledField{

//! Basic Class for differential operators. These classes are passed as a 
//! template parameter to the forms object. Thereby a PDE rather specifies
//! an operator than a specific form.
template<class DATA_TYPE>
class BaseBOperator{
public:
  
  //! Constructor
  BaseBOperator(){

  }

  //! Destructor
  virtual ~BaseBOperator(){
  }

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
    for(UInt i=0;i<nrow;++i){
      for(UInt j=0;j<ncol;++j){
        bMat[i][j] = Complex(1.0,0.0) * realMat[i][j];
      }
    }
  }

  //! Apply the operator matrix on a vector
  virtual void ApplyOp(Vector<DATA_TYPE>& retVec,
                       const LocPointMapped& lp,
                       BaseFE* ptFe,
                       const Vector<DATA_TYPE>& solVec ){
    Matrix<DATA_TYPE> bOp;
    CalcOpMat(bOp,lp,ptFe);
    retVec = bOp * solVec;

  }

  //! Apply the transposed operator matrix on a vector
  virtual void ApplyOpTranspose(Vector<DATA_TYPE>& retVec,
                                const LocPointMapped& lp,
                                BaseFE* ptFe,
                                const Vector<DATA_TYPE>& solVec ){
    
    Matrix<DATA_TYPE> bOp;
    CalcOpMat(bOp,lp,ptFe);
    retVec = Transpose(bOp) * solVec;
  }


  //! Return name of the integrator
  std::string GetName(){
    return name_;
  }

  //! Set the coefficient function of the operator
  virtual void SetCoefFunction(shared_ptr<CoefFunction> coef){
    coef_ = coef;
  }

  //! Additional transformation of the Jacobian determinant (e.g. for PML)
  virtual void TransformJacDet(DATA_TYPE & jacLoc,
                               const LocPointMapped & lp, BaseFE* ptFe){
    return;
  }
protected:

  //! Name of the integrator
  std::string name_;

  //!pointer to coefficient function as used e.g. in Convective operators
  shared_ptr<CoefFunction> coef_;
};
}
#endif
