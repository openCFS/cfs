#ifndef FILE_LAPLACEINT_2001
#define FILE_LAPLACEINT_2001

#include "baseform.hh"

namespace CoupledField
{

/// Class for calculation  element stiffness matrix
template<Integer dim>
class LaplaceInt : public BaseForm<dim>
{
public:
  /// Constructor
  LaplaceInt(BaseFE * aptelem, const ShortInt ndofs);

  /// 
  virtual ~LaplaceInt();

  /// Calculation of stiffmess matrix
  void CalcElemMatrix(Point<dim> * ptCoord, Matrix<Double> & StiffMat);

  virtual void Print(std::ostream * out, const Matrix<Double> Result) const;

protected:

 
private:
  ShortInt DofsPerNode;
};

#ifdef __GNUC__
template class LaplaceInt<2>;
template class LaplaceInt<3>;
#endif
}

#endif // FILE_LAPLACEINT
