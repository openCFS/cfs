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
  MassInt(BaseFE * aptelemt, MaterialData & aMat);

  /// Destructor
  virtual ~MassInt();

  /// Calculation of stiffmess matrix
  void CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMa);


  virtual void Print(std::ostream * out, const Matrix<Double> Result) const;
};

}

#endif // FILE_MASSINT
