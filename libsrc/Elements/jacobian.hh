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

} // end of namespace
#endif
