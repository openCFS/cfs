#ifndef FILE_NLINACOUSTIC2
#define FILE_NLINACOUSTIC2

#include "Utils/ApproxData.hh"
#include "baseForm.hh"

namespace CoupledField
{

  /// Class for calculation nonlinear term in Kuznetsov's and Westervelt's equation
class nLinAcoustic2 : public BaseForm
{
public:
  /// Constructor
  nLinAcoustic2(Boolean axi);

  /// 
  virtual ~nLinAcoustic2();

  /// Calculation of stiffmess matrix
  void CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMat);

  //
  virtual void SetNonLinMethod(std::string atype);
  
  virtual void Print(std::ostream * out, const Matrix<Double> Result) const;

protected: 
private:
  NonLinMethod nonLinType_;

};


} // end of namespace

#endif // FILE_NLINACOUSTIC2
