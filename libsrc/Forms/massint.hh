#ifndef FILE_MASSINT
#define FILE_MASSINT

#include "baseform.hh" 

namespace CoupledField
{
/// class for calculation of Mass matrix
template<Integer dim>
class MassInt : public BaseForm<dim>
{
public:
  /// Constructor with pointer to BaseElem
  MassInt(BaseElem * aptelem, const ShortInt ndofs);

  /// Deconstructor
  virtual ~MassInt();

  /// Function for calculation mass matrix 
  virtual void CalcElemMatrix(Point<dim> * ptCoord, Matrix<Double> &);

  ///
  //  virtual void SetMaterial(Material material);

protected:

private:
  ShortInt DofsPerNode;
};

#ifdef __GNUC__
template class MassInt<2>;
template class MassInt<3>;
#endif


}

#endif // FILE_MASSINT
