#ifndef FILE_LAPLACEINT
#define FILE_LAPLACEINT

#include "baseForm.hh"

namespace CoupledField
{

  /// Class for calculation  element stiffness matrix of a laplace operator
class LaplaceInt : public BaseForm
{
public:
  /// Constructor
  LaplaceInt(BaseFE * aptelem, Double laplVal, Boolean axi=FALSE);

  /// Constructor
  LaplaceInt(Double laplVal, Boolean axi=FALSE);

  /// 
  virtual ~LaplaceInt();

  /// Calculation of stiffmess matrix
  void CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMat);

  virtual void Print(std::ostream * out, const Matrix<Double> Result) const;

protected: 
private:
  /// multiplicative value for laplace integration 
  Double laplVal_;
  Boolean isaxi_;
  
};


}

#endif // FILE_LAPLACEINT
