#ifndef FILE_MASSINT
#define FILE_MASSINT

namespace CoupledField
{
/// class for calculation of Mass matrix
template<class Dim>
class MassInt : public BaseForm<Dim>
{
public:
  /// Constructor with pointer to BaseElem
  MassInt(BaseElem * aptelem, const ShortInt ndofs);

  /// Deconstructor
  virtual ~MassInt();

  /// Function for calculation mass matrix 
  virtual void CalcElemMatrix(Dim * ptCoord, Matrix<Double> &);

  ///
  //  virtual void SetMaterial(Material material);

protected:

private:
  ShortInt DofsPerNode;
};

#ifdef __GNUC__
template class MassInt<Point2D>;
template class MassInt<Point3D>;
#endif


}

#endif // FILE_MASSINT
