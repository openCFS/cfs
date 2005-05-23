#ifndef FILE_JACOBIAN_2001
#define FILE_JACOBIAN_2001

#include <Matrix/matrix.hh>
#include <Utils/vector.hh>

namespace CoupledField
{

  //! Class for calculation Jacobian
  /*! In this class we store Jacobian \f$ J\f$, det of Jacobian and 
    calculate inverse of Jacobian matrix \f$ J^{-1}\f$ 
  
    \f$ J = \left( \begin{array}{ccc} x_{\xi} & x_{\eta} & x_{\zeta} \\
    y_{\xi} & y_{\eta} & y_{\zeta}\\
    z_{\xi} & z_{\eta} & z_{\zeta} \end{array}\right)\f$ 
    \f$ J^{-1} = \left( \begin{array}{ccc} \xi_{x} & \xi_{y} & \xi_{z} \\
    \eta_{x} & \eta_{y} & \eta_{z}\\
    \zeta_{x} & \zeta_{y} & \zeta_{z} \end{array}\right)\f$ */

  template<Integer dim>
  class Jacobian
  {
  public:
    Double          detJ;          //!< determinant of Jacobian
    Matrix<Double>  J;            //!< Jacobian
    Matrix<Double>  Jinv;         //!< inverse Jacobian

    //! Constructor
    Jacobian();

    //! Return column for derivatives with respect to x in inverse Jacobian
    /*! 
      \param JinvX first column of inverse Jacobian (used for computation of global
      derivative with respect to x-coordiante
    */
    void GetJinvX(Vector<Double> & JinvX);

    //! Return column with derivatives respect to y in inverse Jacobian
    /*! 
      \param JinvY second column of inverse Jacobian (used for computation of global
      derivative with respect to y-coordiante
    */
    void GetJinvY(Vector<Double> & JinvY);

    //! Return column with derivatives respect to z in inverse Jacobian
    /*! 
      \param JinvZ third column of inverse Jacobian (used for computation of global
      derivative with respect to z-coordiante
    */
    void GetJinvZ(Vector<Double> & JinvZ);


    //! Returns the inverse Jacobian
    /*! 
      \param Jinv inverse Jacobian (used for computation of global derivative)
    */
    void GetJinv(Matrix<Double> & jInv);
  };

#if defined(__GNUC__)
  template class Jacobian<2>;
  template class Jacobian<3>;
#endif

} // end of namespace
#endif
