#ifndef FILE_MASSINT_1
#define FILE_MASSINT_1

#include "baseForm.hh"

namespace CoupledField
{

  /// Class for calculation  element stiffness matrix
class MassInt : public BaseForm
{
public:
  /// Constructor
  MassInt(BaseFE * aptelem);

  /// 
  virtual ~MassInt();

  /// Calculation of stiffmess matrix
  void CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMat);

  virtual void Print(std::ostream * out, const Matrix<Double> Result) const;
};

}

#endif // FILE_MASSINT
