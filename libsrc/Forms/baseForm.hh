#ifndef FILE_BASEFORM
#define FILE_BASEFORM

#include <Elements/basefe.hh>

namespace CoupledField
{

//! Base class for calculation of bdb element matrices
class BaseForm 
{
public:
  //! Constructor
  BaseForm(BaseFE * aptelem);

  //! Deconstructor
  virtual ~BaseForm();

  //! Virtual function, implemented in derived classes
  virtual void CalcElementMatrix(Matrix<Double>& ptCoord, Matrix<Double> & StiffMat);

  //! Prints the bilinear form
  virtual void Print(std::ostream * out, const Matrix<Double> Result) const;

protected:

  //! Ptr to base element
  BaseFE  * ptelem;
};

} //end namespace

#endif // FILE_BASEFORM
