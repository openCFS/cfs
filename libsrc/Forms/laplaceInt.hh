#ifndef FILE_LAPLACEINT
#define FILE_LAPLACEINT

#include "baseForm.hh"

namespace CoupledField
{

  /// Class for calculation  element stiffness matrix
class LaplaceInt : public BaseForm
{
public:
  /// Constructor
  LaplaceInt(BaseFE * aptelem, Double laplVal);

  /// 
  virtual ~LaplaceInt();

  /// Calculation of stiffmess matrix
  void CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMat);

  virtual void Print(std::ostream * out, const Matrix<Double> Result) const;

protected: 
private:
  /// multiplicative value for laplace integration 
  Double laplVal_;
  
};

}

#endif // FILE_LAPLACEINT
