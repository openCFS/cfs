#ifndef FILE_NLINACOUSTIC2
#define FILE_NLINACOUSTIC2

#include "Utils/ApproxData.hh"
#include "baseForm.hh"

namespace CoupledField
{

  //! Class for calculation nonlinear term in Kuznetsov's equation
class nLinAcoustic2 : public BaseForm
{
public:
  //! Constructor
  nLinAcoustic2(Double factor, Boolean axi);

  //! Destructor 
  virtual ~nLinAcoustic2();

  //! Calculation of element matrix
  void CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMat);

  //! the solution is needed for setting up the matrix
  virtual void SetActElemSol(Vector<Double>& asol) 
  { sol_ = asol;};

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

  Vector<Double> sol_;        //!< solution

  Double factor_;             //!< multiplicative value for integrator


};


} // end of namespace

#endif // FILE_NLINACOUSTIC2
