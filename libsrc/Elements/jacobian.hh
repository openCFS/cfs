#ifndef FILE_JACOBIAN_2001
#define FILE_JACOBIAN_2001

#include "matrix.hh"

namespace CoupledField
{

  //! Class for calculation Jacobian
  /*! In this class we store Jacobian of transformation element in standart, det of Jacobian and calculate inverse of Jacobian matrix */

template<class Dim>
class Jacobian
{
  public:
    Double          detJ;          //!< determinant of Jacobian
    Matrix<Double>  J;            //!< Jacobian
    Matrix<Double>  Jinv;         //!< inverse Jacobian

    //! Constructor
    Jacobian();

    //! Return colomn with derivatives respect to x in inverse Jacobian
    void GetJinvX(Vector<Double> & JinvX);

    //! Return colomn with derivatives respect to y in inverse Jacobian
    void GetJinvY(Vector<Double> & JinvY);
};

#ifdef __GNUC__
template class Jacobian<Point2D>;
template class Jacobian<Point3D>;
#endif

template <>
Jacobian<Point2D>::Jacobian()
{
  J.Resize(2,2);
  Jinv.Resize(2,2);
}

template <>
Jacobian<Point3D>::Jacobian()
{
  J.Resize(3,3);
  Jinv.Resize(3,3);
}

template<>
void Jacobian<Point3D>::GetJinvX(Vector<Double> & JinvX)
{
  JinvX.Resize(3);

  JinvX[0]=Jinv[0][0];
  JinvX[1]=Jinv[1][0];
  JinvX[2]=Jinv[2][0];
}

template<>
void Jacobian<Point2D>::GetJinvX(Vector<Double> & JinvX)
{
  JinvX.Resize(2);
  JinvX[0]=Jinv[0][0];
  JinvX[1]=Jinv[1][0];
}

template<>
void Jacobian<Point2D>::GetJinvY(Vector<Double> & JinvY)
{
  JinvY.Resize(2);

  JinvY[0]=Jinv[0][1];
  JinvY[1]=Jinv[1][1];
}

template<>
void Jacobian<Point3D>::GetJinvY(Vector<Double> & JinvY)
{
  JinvY.Resize(3);

  JinvY[0]=Jinv[0][1];
  JinvY[1]=Jinv[1][1];
  JinvY[2]=Jinv[2][1];
}
} // end of namespace
#endif
