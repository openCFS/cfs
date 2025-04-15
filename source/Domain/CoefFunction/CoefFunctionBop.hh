#ifndef COEFFUNCTIONBOP_HH_
#define COEFFUNCTIONBOP_HH_

#include <boost/shared_ptr.hpp>

#include "CoefFunction.hh"
#include "Forms/BiLinForms/BiLinearForm.hh"

namespace CoupledField {

//class BiLinearForm;

/* Can be used to pass a bilinear form as a coefFunction for evalaution besides the assembly class */
class CoefFunctionBop : public CoefFunction
{
public:

  CoefFunctionBop();

  virtual ~CoefFunctionBop() { }

  virtual string GetName() const { return "CoefFunctionBop"; }

  /** set the form such that the proper transfer function can be found */
  void SetForm(BiLinearForm* form) {
    this->form_ = form;
  }

  /** only to query the name to find the proper transfer function */
  const BiLinearForm* GetForm() const  {
    return form_;
  }

  template <class T>
  void GetTensor(Matrix<T>& coefMat, const LocPointMapped& lpm, BaseFE* ptFe);

  virtual void CalcElementMatrixLpm(Matrix<Double>& elemMat, const LocPointMapped& lpm, BaseFE* ptFe, bool overrideIsSurfOpt);

  template <class T>
  void GetScalar(T& scal, const LocPointMapped& lpm);

  template <class T>
  void GetVector(Vector<T>& scal, const LocPointMapped& lpm);

protected:

  /** we store the form such that we can identify the proper transfer function */
  BiLinearForm* form_;
};

}

#endif /* COEFFUNCTIONBOP_HH_ */
