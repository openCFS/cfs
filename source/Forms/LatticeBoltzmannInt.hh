#ifndef LBMFORM_HH_
#define LBMFORM_HH_

#include "Forms/baseForm.hh"

namespace CoupledField
{

class BaseMaterial;
class EntityIterator;
struct Elem;
template <class TYPE> class Matrix;


/** This form sets up the adjoint system for sensitivity analysis by LatticeBoltzmannPDE */
class LatticeBoltzmannInt : public BaseForm
{
public:

  //! Constructor with pointer to BaseElem
  LatticeBoltzmannInt(BaseMaterial* matData, SubTensorType type = FULL,  bool coordUpdate = false);

  //! Destructor
  virtual ~LatticeBoltzmannInt() { }

  //! Compute element matrix associated to ADB form
  void CalcElementMatrix(Matrix<Double>& elemMat, EntityIterator& ent1, EntityIterator& ent2);

  void CalcElementMatrix(Matrix<Double>& elemMat, EntityIterator& ent1, EntityIterator& ent2, const DesignElement::Type direction)
  {
    CalcElementMatrix( elemMat, ent1, ent2);
  }

  //! returns dimension of D matrix
  unsigned int getDimD() { return 1; }

  //! returns nr. of degrees of freedom
  unsigned int getNrDofs();

  //! Query material type for \f$D\f$ tensor
  MaterialType getDMaterialType() { return LBM_MATERIAL; }

protected:

  //! bool for signaling that D matrix is non-constant

  //! In some cases, e.g. in non-linear computations, it may be
  //! necessary to compute the D matrix for each integration point
  //! individually. This attribute is used to signal when the latter is
  //! required.
  bool updateDMatInEveryIP_;

};

}  // end of namespace


#endif /* LBMFORM_HH_ */
