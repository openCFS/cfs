#ifndef FILE_LAPLACEINT
#define FILE_LAPLACEINT

#include "baseform.hh"

namespace CoupledField
{

/// Class for calculation  element stiffness matrix
template<class Dim>
class LaplaceInt : public BaseForm<Dim>
{
public:
  /// Constructor
  LaplaceInt(BaseElem * aptelem, const ShortInt ndofs);

  /// 
  virtual ~LaplaceInt();

  /// Calculation of stiffmess matrix
  void CalcElemMatrix(Dim * ptCoord, Matrix<Double> & StiffMat);

  virtual void Print(std::ostream * out, const Matrix<Double> Result) const;

protected:

 
private:
  ShortInt DofsPerNode;
};

#ifdef __GNUC__
template class LaplaceInt<Point2D>;
template class LaplaceInt<Point3D>;
#endif
}

#endif // FILE_LAPLACEINT
