// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_SMOOTHINT
#define FILE_SMOOTHINT

#include "Forms/bdbInt.hh"

#include "General/defs.hh"
#include "General/environment.hh"

namespace CoupledField {
class BaseMaterial;
template <class TYPE> class Matrix;
}  // namespace CoupledField

namespace CoupledField
{
  
  
  /// base class for calculation of linear elasticity
class SmoothInt : public BDBInt
{
public:

  //! Constructor
  SmoothInt(BaseMaterial* matData, SubTensorType type = FULL,
            bool coordUpdate = false );
  
  //! Destructor
  virtual ~SmoothInt();
  
protected:    
  
  /** @see BaseForm::CalcBMat() */
  virtual void CalcBMat(Matrix<Double> & bMat, UInt ip, const Matrix<Double> & ptCoord);

  //! calculate the data-matrix for 2D plain-strain
  virtual void calcDMat(Matrix<Double> & dMat, UInt ip, Matrix<Double> & ptCoord);

  //! returns dimension of D matrix
  virtual UInt getDimD() {
    return dimD_; 
  };

  //! Query material type for \f$D\f$ tensor
  MaterialType getDMaterialType() { return MECH_STIFFNESS_TENSOR; }

  
  //! returns nr. of degrees of freedom
  virtual UInt getNrDofs() {
    return nrDofs_;
  };

  private:

  //dimension of Dmatrix
  UInt dimD_;

  //! number of degrees 
  UInt nrDofs_;
};
  
} //end namespace

#endif // FILE_SMOOTHINT
