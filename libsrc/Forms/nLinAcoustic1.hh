#ifndef FILE_NLINACOUSTIC1
#define FILE_NLINACOUSTIC1

#include "Utils/ApproxData.hh"
#include "baseForm.hh"

namespace CoupledField
{

  /// Class for calculation nonlinear term in Kuznetsov's and Westervelt's equation
class nLinAcoustic1 : public BaseForm
{
public:
  /// Constructor
  nLinAcoustic1(Boolean axi);

  /// 
  virtual ~nLinAcoustic1();

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

#endif // FILE_NLINACOUSTIC1
