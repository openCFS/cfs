#ifndef FILE_NLINACOUSTIC1
#define FILE_NLINACOUSTIC1

#include "Utils/ApproxData.hh"
#include "baseForm.hh"

namespace CoupledField
{

  //! Class for calculation nonlinear term in Kuznetsov's and Westervelt's equation
class nLinAcoustic1 : public BaseForm
{
public:
  //! Constructor
  nLinAcoustic1(Double factor, Boolean axi);

  //! Destructor
  virtual ~nLinAcoustic1();

  //! Calculation of element matrix
  void CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMat);

  //! the first time derivative is needed for setting up the matrix
  virtual void SetActElemSolDeriv1(Vector<Double>& solderiv1) 
  { solderiv1_ = solderiv1;};

  //! set multiplicative factor for matrix
  virtual void SetFactor(Double factor) 
  {factor_ = factor;};

  //! i.e. fixedpoint or linesearch
  virtual void SetNonLinMethod(std::string atype);
  
  //! print out the element matrix
  virtual void Print(std::ostream * out, const Matrix<Double> Result) const;

protected: 
private:

  NonLinMethod nonLinType_;   //!< nonlinear iteration method

  Vector<Double> solderiv1_;  //!< first time derivative of solution

  Double factor_;             //!< multiplicative value for integrator
};


} // end of namespace

#endif // FILE_NLINACOUSTIC1
